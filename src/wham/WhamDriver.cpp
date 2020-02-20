#include "WhamDriver.h"

WhamDriver::WhamDriver(const std::string& options_file):
	options_file_(options_file)
{
	// Read input file into a ParameterPack
	InputParser input_parser;
	input_parser.parseFile(options_file_, input_parameter_pack_);

	using KeyType = ParameterPack::KeyType;

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
	simulations_.clear();
	simulations_.reserve(num_simulations);
	for ( int i=0; i<num_simulations; ++i ) {
		simulations_.emplace_back(
			data_set_labels[i], t_min[i], t_max[i], wham_options_.T, wham_options_.floor_t, op_registry_
		);
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
		//order_parameters_.push_back( OrderParameter(*(op_pack_ptrs[p]), simulations_, wham_options_.floor_t) );
	}


	//---- Initial guess -----//

	// TODO: Option to read from input
	f_bias_guess_.assign(num_simulations, 0.0);
	std::string f_bias_guess_file;
	bool found = input_parameter_pack_.readString("InitialGuess", KeyType::Optional, f_bias_guess_file);
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


}


void WhamDriver::run_driver()
{
	int num_simulations = simulations_.size();
	bool be_verbose = true;


	//----- Solve WHAM equations -----//

	// Initial guess
	// - Shift so that f[0] = 0
	std::vector<double> f_bias_init(num_simulations, 0.0);
	double f_shift = f_bias_guess_[0];
	for ( int j=0; j<num_simulations; ++j ) {
		f_bias_init[j] = f_bias_guess_[j] - f_shift;
	}

	// Solve for optimal free energies of biasing
	Wham wham(data_summary_, op_registry_, simulations_, order_parameters_, biases_, f_bias_init, wham_options_.tol);
	f_bias_opt_ = wham.get_f_bias_opt();

	// Print optimal biasing free energies to file
	std::string file_name = "f_bias_WHAM.out";
	std::ofstream ofs(file_name);
	ofs << "# F_bias [k_B*T] for each window after minimization\n";
	for ( int i=0; i<num_simulations; ++i ) {
		ofs << std::setprecision(7) << f_bias_opt_[i] << "\n";
	}
	ofs.close(); ofs.clear();


	//----- Output -----//

	// 1-variable outputs
	for ( unsigned i=0; i<output_f_x_.size(); ++i ) {
		OrderParameter& x = order_parameters_[ output_f_x_[i] ];
		if ( be_verbose ) {
			std::cout << "Computing F_WHAM(" << x.get_name() << ")\n";
		}

		// "Manually" unbiased distributions (non-consensus, unshifted)
		x.set_unbiased_distributions( wham.manuallyUnbiasDistributions( x.get_name() ) );

		// "Raw" distributions (i.e. using only data from each individual simulation)
		x.printRawDistributions();

		// TODO: shifted distributions

		// F_WHAM(x)
		// - TODO errors
		x.set_wham_distribution( wham.compute_consensus_f_x_unbiased( x.get_name() ) );
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
		if ( be_verbose ) {
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
