#include "WhamDriver.h"

WhamDriver::WhamDriver(const std::string& options_file):
	options_file_(options_file)
{
	setup_timer_.start();

	if ( not be_quiet_ ) {
		std::cout << "WHAM\n";
		if ( OpenMP::is_enabled() ) {
			std::cout << "  Using " << OpenMP::get_max_threads() << " threads (max.)\n";
		}
	}

	// Read input file into a ParameterPack
	InputParser input_parser;
	input_parser.parseFile(options_file_, input_parameter_pack_);

	using KeyType = ParameterPack::KeyType;
	bool found = false;

	input_parameter_pack_.readFlag("be_verbose", KeyType::Optional, be_verbose_);

	input_parameter_pack_.readNumber("Temperature", KeyType::Required, wham_options_.T);
	wham_options_.kBT = Constants::k_B * wham_options_.T;

	input_parameter_pack_.readNumber("SolverTolerance", KeyType::Required, wham_options_.tol);

	// Read the data summary
	// - This sets the number of simulations expected in all files subsequently parsed
	input_parameter_pack_.readString("DataSummaryFile", KeyType::Required, data_summary_file_);
	data_summary_ = DataSummary(data_summary_file_, input_parameter_pack_);
	const auto& data_set_labels = data_summary_.get_data_set_labels();
	const auto& t_min           = data_summary_.get_t_min();
	const auto& t_max           = data_summary_.get_t_max();
	const int num_simulations   = data_set_labels.size();

	// Register all order parameters and associated time series files
	op_registry_ = OrderParameterRegistry(input_parameter_pack_, data_summary_);
	const auto& op_names = op_registry_.get_names();
	const int num_ops    = op_registry_.getNumberOfOrderParameters();

	// Load simulation data
	if ( not be_quiet_ ) {
		std::cout << "  Loading data ...\n";
	}
	simulations_.clear();
	simulations_.reserve(num_simulations);
	for ( int i=0; i<num_simulations; ++i ) {
		simulations_.emplace_back(
			data_set_labels[i], t_min[i], t_max[i], wham_options_.T, wham_options_.floor_t, op_registry_
		);
	}

	// Bootstrap error estimation (optional)
	found = input_parameter_pack_.readNumber("num_bootstrap_samples", KeyType::Optional, num_bootstrap_samples_);
	if ( found ) {
		error_method_ = ErrorMethod::Bootstrap;
	}


	//----- Biasing Potentials -----//

	// Read the biasing parameters file 
	input_parameter_pack_.readString("BiasesLogFile", KeyType::Required, 
	                                 biases_log_file_);

	// Find the bias input packs
	ParameterPack bias_file_pack;
	input_parser.parseFile(biases_log_file_, bias_file_pack);
	std::vector<const ParameterPack*> bias_input_pack_ptrs = 
			bias_file_pack.findParameterPacks("Bias", KeyType::Required);
	int num_biases = bias_input_pack_ptrs.size();
	if ( num_biases != num_simulations ) {
		std::stringstream err_ss;
		err_ss << "Mismatch between number of biases (" << num_biases << ") and number of simulations ("
		       << num_simulations << ")\n";
		throw std::runtime_error( err_ss.str() );
	}

	for ( int i=0; i<num_biases; ++i ) {
		// TODO allow bias log and data summary to report simulations in different orders,
		//      or allow one to be a subset of the others
		// - Use data set labels to match everything
		// TODO: variable T between simulations
		biases_.push_back( Bias(*(bias_input_pack_ptrs[i]), wham_options_.kBT) );

		if ( data_set_labels[i] != biases_.back().get_data_set_label() ) {
			throw std::runtime_error("Mismatch between data set labels in biasing parameters file and data summary file");
		}
	}

#ifdef DEBUG
	std::cout << "DEBUG: " << num_biases << " biases registered\n";
	for ( int j=0; j<num_biases; ++j ) {
		int num_potentials = biases_[j].get_order_parameter_names().size();

		std::cout << "  data_set = " << biases_[j].get_data_set_label() << "\n"
		          << "  num_potentials = " << num_potentials << "\n";
		for ( int k=0; k<num_potentials; ++k ) {
			std::cout << "    op = " << biases_[j].get_order_parameter_names()[k] << "\n";
		}
	}
#endif // DEBUG


	//----- Order Parameters -----//

	// Set up order parameters
	// - Time series data is shared with Simulations
	std::vector<const ParameterPack*> op_pack_ptrs = 
			input_parameter_pack_.findParameterPacks("OrderParameter", KeyType::Required);
	order_parameters_.clear();
	for ( int p=0; p<num_ops; ++p ) {
		order_parameters_.push_back( OrderParameter(op_names[p], *(op_pack_ptrs[p]), simulations_) );
	}


	//---- Initial guess -----//

	// TODO: Option to read from input
	f_bias_guess_.assign(num_simulations, 0.0);
	std::string f_bias_guess_file;
	found = input_parameter_pack_.readString("InitialGuess", KeyType::Optional, f_bias_guess_file);
	if ( found ) {
	}


	//----- Output Options -----//

	// Default: Print all permutations of F(x) and F(x,y)
	for ( int p=0; p<num_ops; ++p ) {
		output_f_x_.push_back( p );
		for ( int q=p+1; q<num_ops; ++q ) {
			output_f_x_y_.push_back( {{p, q}} );
		}
	}

	// Check whether the user requested only certain distributions
	// TODO move to function, generalize to n-dimensional outputs, and (ideally)
	//      make this cleaner overall
	const ParameterPack* output_pack_ptr = 
			input_parameter_pack_.findParameterPack("Output", KeyType::Optional);	
	
	bool print_everything = true;
	if ( output_pack_ptr != nullptr ) {

		output_pack_ptr->readFlag("PrintAll", KeyType::Optional, print_everything);
		std::vector<const std::vector<std::string>*> vecs = 
				output_pack_ptr->findVectors("F_WHAM", KeyType::Optional);

		if ( (not print_everything) or (vecs.size() > 0) ) {
			output_f_x_.clear();
			output_f_x_y_.clear();

			int num_outputs = vecs.size();
			for ( int i=0; i<num_outputs; ++i ) {
				const std::vector<std::string>& vec = *(vecs[i]);  // convenient alias for this vector
				int num_ops_for_output = vec.size();

				// Map names to indices
				std::vector<int> op_indices(num_ops_for_output);
				for ( int j=0; j<num_ops_for_output; ++j ) {
					op_indices[j] = op_registry_.get_index(vec[j]);
				}

				// Store indices
				if ( num_ops_for_output == 0 ) {
					throw std::runtime_error("F_WHAM with no order parameters is undefined");
				}
				else if ( num_ops_for_output == 1 ) {
					output_f_x_.push_back( op_indices[0] );
				}
				else if ( num_ops_for_output == 2 ) {
					output_f_x_y_.push_back( {{ op_indices[0], op_indices[1] }} );
				}
				else {
					throw std::runtime_error("F_WHAM calcualtions only support up to 2 order parameters");
				}
			}
		}
	}

#ifdef DEBUG
	std::cout << "DEBUG: " << num_ops << " order parameters registered\n";
	for ( int i=0; i<num_ops; ++i ) {
		int num_time_series = order_parameters_[i].time_series_.size();
		std::cout << "  op = " << order_parameters_[i].name_ << ": " << num_time_series << " time series\n";
		for ( int j=0; j<num_time_series; ++j ) {
			std::cout << "    file = " << order_parameters_[i].get_time_series(j).get_file() << "\n"
			          << "    size = " << order_parameters_[i].get_time_series(j).size() << " data points\n";
		}
	}
#endif /* DEBUG */

	setup_timer_.stop();
}


void WhamDriver::run_driver()
{
	driver_timer_.start();

	if ( not be_quiet_ ) {
		std::cout << "  Solving ...\n";
	}

	int num_simulations = simulations_.size();

	//----- Solve WHAM equations -----//

	// Initial guess
	// - Shift so that f[0] = 0 (free to choose free energy offset)
	std::vector<double> f_bias_init(num_simulations, 0.0);
	double f_shift = f_bias_guess_[0];
	for ( int j=0; j<num_simulations; ++j ) {
		f_bias_init[j] = f_bias_guess_[j] - f_shift;
	}

	// Solve for optimal free energies of biasing
	solve_wham_timer_.start();
	Wham wham(data_summary_, op_registry_, simulations_, order_parameters_, biases_, f_bias_init, wham_options_.tol);
	solve_wham_timer_.stop();
	f_bias_opt_ = wham.get_f_bias_opt();

	if ( not be_quiet_ ) {
		std::cout << "  Done\n";
		std::cout << "  Optimal biasing free energies (in kBT):\n";
		for ( int i=0; i<num_simulations; ++i ) {
			std::cout << "  " << i+1 << ":  " << f_bias_opt_[i] << "  (initial: " << f_bias_guess_[i] << ")\n"; 
		}
	}


	//----- Estimate Errors -----//

	// TODO: Move this section to a different function?

	bool calc_error = ( error_method_ != ErrorMethod::None );

	std::vector<PointEstimator<double>> bootstrap_samples_f_bias;
	std::vector< std::vector<PointEstimator<double>> > bootstrap_samples_f_x;

	if ( error_method_ == ErrorMethod::Bootstrap ) {
		bootstrap_timer_.start();

		if ( not be_quiet_ ) {
			std::cout << "  Estimate errors using bootstrap subsampling ...\n";
		}

		// Allocate memory for bootstrap samples
		// - Dimensions: #outputs x #bins/output
		bootstrap_samples_f_bias.resize(num_simulations);
		int num_output_f_x = output_f_x_.size();
		bootstrap_samples_f_x.resize(num_output_f_x);
		for ( int i=0; i<num_output_f_x; ++i ) {
			OrderParameter& x = order_parameters_[ output_f_x_[i] ];

			int num_bins_x = x.get_bins().get_num_bins();
			bootstrap_samples_f_x[i].resize(num_bins_x);

			for ( int b=0; b<num_bins_x; ++b ) {
				bootstrap_samples_f_x[i][b].clear();
				bootstrap_samples_f_x[i][b].reserve(num_bootstrap_samples_);
			}
		}

		// Prepare subsamplers
		std::vector<BootstrapSubsampler> subsamplers;
		std::vector<std::vector<int>> resample_indices(num_simulations);
		for ( int j=0; j<num_simulations; ++j ) {
			int num_samples_j = simulations_[j].get_num_samples();
			// TODO: different seeds depending on debug mode flag
			subsamplers.emplace_back( num_samples_j, Random::getDebugSequence() );
			resample_indices[j].reserve( num_samples_j );
		}

		std::vector<Simulation> bootstrap_simulations(num_simulations);
		std::vector<OrderParameter> bootstrap_ops = order_parameters_;
		for ( int s=0; s<num_bootstrap_samples_; ++s ) {
			// User feedback
			if ( (not be_quiet_) and ((s == 0) or (((s+1) % 25) == 0)) ) {
				std::cout << "    sample " << s+1 << " of " << num_bootstrap_samples_ << "\n";
			}

			// Subsample each simulation
			// - TODO: OpenMP?
			subsample_timer_.start();
			for ( int j=0; j<num_simulations; ++j ) {
				subsamplers[j].generate( resample_indices[j] );
				bootstrap_simulations[j].setShuffledFromOther( simulations_[j], resample_indices[j] );
			}
			subsample_timer_.stop();

			int num_ops = order_parameters_.size();
			for ( int p=0; p<num_ops; ++p ) {
				bootstrap_ops[p].set_simulations(bootstrap_simulations);
			}

			// Re-solve WHAM equations
			// - TODO: option to re-solve with new data rather than reallocating for each loop?
			solve_wham_timer_.start();
			Wham bootstrap_wham( data_summary_, op_registry_, bootstrap_simulations, bootstrap_ops, 
													 biases_, f_bias_opt_, wham_options_.tol );
			solve_wham_timer_.stop();

			// Save bootstrap estimates for f_bias_opt
			const auto& bootstrap_f_bias_opt = bootstrap_wham.get_f_bias_opt();
			for ( int j=0; j<num_simulations; ++j ) {
				bootstrap_samples_f_bias[j].addSample( bootstrap_f_bias_opt[j] );
			}

			// Compute boostrap estimates of each F(x)
			for ( int i=0; i<num_output_f_x; ++i ) {
				OrderParameter& x = order_parameters_[ output_f_x_[i] ];
				auto bootstrap_f_x = bootstrap_wham.compute_consensus_f_x_unbiased( x.get_name() );

				int num_bins_x = x.get_bins().get_num_bins();
				for ( int b=0; b<num_bins_x; ++b ) {
					// Only save this sample for if F(x) if its value is finite (i.e. bin has samples in it)
					// - Otherwise, statistics over the samples will be corrupted
					if ( bootstrap_f_x.sample_counts[b] > 0 ) {
						bootstrap_samples_f_x[i][b].addSample( bootstrap_f_x.f_x[b] );
					}
				}
			}
		};

		// Finalize errors
		error_f_bias_opt_.resize(num_simulations);
		for ( int j=0; j<num_simulations; ++j ) {
			error_f_bias_opt_[j] = bootstrap_samples_f_bias[j].std_dev();
		}

		// TODO: Compute errors for F(x) here instead of below

		bootstrap_timer_.stop();
	} // end bootstrap resampling


	//----- Output -----//

	print_output_timer_.start();

	// Print optimal biasing free energies to file
	std::string file_name("f_bias_WHAM.out");
	std::ofstream ofs(file_name);
	ofs << "# F_bias [k_B*T]";
	if ( calc_error ) { ofs << " +/- error"; }
	ofs << " for each window after minimization\n";
	for ( int j=0; j<num_simulations; ++j ) {
		ofs << std::setprecision(7) << f_bias_opt_[j];
		if ( calc_error ) {
			ofs << "  " << std::setprecision(7) << error_f_bias_opt_[j];
		}
		ofs << "\n";
	}
	ofs.close(); ofs.clear();

	if ( error_method_ == ErrorMethod::Bootstrap ) {
		ofs.open("bootstrap_convergence.out");
		ofs << "# Convergence of bootstrap subsampling\n"
		    << "# n_B[samples]  F_bias(last_window)[kBT]\n";

		const auto& all_samples = bootstrap_samples_f_bias.back().get_samples();
		std::vector<double> samples(1, all_samples[0]);
		for ( int s=1; s<num_bootstrap_samples_; ++s ) {
			samples.push_back( all_samples[s] );
			ofs << s+1 << "  " << std::setprecision(7) << Statistics::std_dev(samples) << "\n";
		}
		ofs.close(); ofs.clear();
	}

	// 1-variable outputs
	for ( unsigned i=0; i<output_f_x_.size(); ++i ) {
		OrderParameter& x = order_parameters_[ output_f_x_[i] ];
		if ( not be_quiet_ ) {
			std::cout << "Computing F_WHAM(" << x.get_name() << ")\n";
		}

		// "Manually" unbiased distributions (non-consensus, unshifted)
		x.set_unbiased_distributions( wham.manuallyUnbiasDistributions( x.get_name() ) );

		// "Raw" distributions (i.e. using only data from each individual simulation)
		x.printRawDistributions();

		// TODO: "shifted" distributions

		// F_WHAM(x)
		auto wham_distribution = wham.compute_consensus_f_x_unbiased( x.get_name() );
		if ( calc_error and error_method_ == ErrorMethod::Bootstrap ) {
			int num_bins_x = x.get_bins().get_num_bins();
			std::vector<double> err_f_x(num_bins_x);
			for ( int b=0; b<num_bins_x; ++b ) {
				if ( bootstrap_samples_f_x[i][b].get_num_samples() >= 2 ) {
					err_f_x[b] = bootstrap_samples_f_x[i][b].std_dev();
				}
				else {
					err_f_x[b] = -1.0;  // junk value (TODO: print as nan instead?)
				}
			}
			wham_distribution.error_f_x = err_f_x;  // TODO set fxn
		}
		x.set_wham_distribution( wham_distribution );
		x.printWhamResults();

		// "Rebias" consensus histograms (for validation)
		std::vector<Distribution> f_x_rebiased;
		for ( int j=0; j<num_simulations; ++j ) {
			f_x_rebiased.push_back(
				wham.compute_consensus_f_x_rebiased( x.get_name(), data_summary_.get_data_set_label(j) )
			);
		}
		x.set_rebiased_distributions( f_x_rebiased );
		x.printRebiasedDistributions();

		x.printStats();
	}

	// 2-variable outputs
	for ( unsigned i=0; i<output_f_x_y_.size(); ++i ) {
		// F_WHAM(x,y)
		OrderParameter& x = order_parameters_[ output_f_x_y_[i][0] ];
		OrderParameter& y = order_parameters_[ output_f_x_y_[i][1] ];
		if ( not be_quiet_ ) {
			std::cout << "Computing F_WHAM(" << x.get_name() << ", " << y.get_name() << ")\n";
		}

		std::vector<std::vector<double>> p_x_y_wham, f_x_y_wham; 
		std::vector<std::vector<int>> sample_counts_x_y;

		wham.compute_consensus_f_x_y_unbiased(
			x.get_name(), y.get_name(),
			p_x_y_wham, f_x_y_wham, sample_counts_x_y
		);

		// Print results
		print_f_x_y( x, y, p_x_y_wham, f_x_y_wham, sample_counts_x_y );
	}

	print_output_timer_.stop();

	driver_timer_.stop();
}


void WhamDriver::print_f_x_y(
	const OrderParameter& x, const OrderParameter& y,
	// Consensus distributions for F_k(x,y)
	const std::vector<std::vector<double>>& p_x_y_wham,
	const std::vector<std::vector<double>>& f_x_y_wham,
	const std::vector<std::vector<int>>&    sample_counts_x_y
) const
{
	// Working variables
	const Bins& bins_x = x.get_bins();
	const Bins& bins_y = y.get_bins();
	int num_bins_x = bins_x.get_num_bins();
	int num_bins_y = bins_y.get_num_bins();
	std::string file_name, sep;
	std::ofstream ofs;

	// Print bins
	std::vector<const OrderParameter*> op_ptrs = {{ &x, &y }};
	for ( unsigned i=0; i<op_ptrs.size(); ++i ) {
		file_name = "bins_" + op_ptrs[i]->get_name() + ".out";
		ofs.open(file_name);
		ofs << "# Bins for " << op_ptrs[i]->get_name() << " for F_WHAM(" << x.get_name() << "," << y.get_name() << ")\n"
				<< "# " << op_ptrs[i]->get_name() << "\n";
		const auto& bins = op_ptrs[i]->get_bins();
		int num_bins = bins.get_num_bins();
		for ( int j=0; j<num_bins; ++j ) {
			ofs << bins[j]  << "\n";
		}
		ofs.close(); ofs.clear();
	}

	// Print F(x,y)  (TODO: shift so that F=0 at the minimum)
	file_name = "F_" + x.get_name() + "_" + y.get_name() + "_WHAM.out";
	ofs.open(file_name);
	sep = " ";
	for ( int i=0; i<num_bins_x; ++i ) {
		for ( int j=0; j<num_bins_y; ++j ) {
			Distribution::print_free_energy(ofs, f_x_y_wham[i][j], sample_counts_x_y[i][j]);
			ofs << sep;
		}
		ofs << "\n";
	}
	ofs.close(); ofs.clear();

	// Print sample distribution
	file_name = "samples_" + x.get_name() + "_" + y.get_name() + ".out";
	ofs.open(file_name);
	sep = " ";
	for ( int i=0; i<num_bins_x; ++i ) {
		for ( int j=0; j<num_bins_y; ++j ) {
			ofs << sample_counts_x_y[i][j] << sep;
		}
		ofs << "\n";
	}
	ofs.close(); ofs.clear();
}
