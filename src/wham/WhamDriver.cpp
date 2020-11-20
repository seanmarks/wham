// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "WhamDriver.h"

#include "Bootstrap_F_x.hpp"
#include "BiasedDistributions.hpp"
#include "ManuallyUnbiasedDistributions.hpp"
#include "RebiasedDistributions.hpp"
#include "WhamStatistics.hpp"


WhamDriver::WhamDriver(const std::string& options_file):
	options_file_(options_file)
{
	setup_timer_.start();

	if ( ! be_quiet_ ) {
		std::cout << "WHAM\n";
		if ( OpenMP::is_enabled() ) {
			std::cout << "  Using " << OpenMP::get_max_threads() << " threads (max.)\n";
		}
		std::cout << std::flush;
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
	const auto& data_set_labels = data_summary_.getDataSetLabels();
	const auto& t_min           = data_summary_.get_t_min();
	const auto& t_max           = data_summary_.get_t_max();
	const int num_simulations   = data_set_labels.size();

	// Register all order parameters and associated time series files
	op_registry_ = OrderParameterRegistry(input_parameter_pack_, data_summary_);
	const auto& op_names = op_registry_.getNames();
	const int num_ops    = op_registry_.getNumberOfOrderParameters();

	// Load simulation data
	if ( ! be_quiet_ ) {
		std::cout << "  Loading data ...\n" << std::flush;
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
	input_parameter_pack_.readString("BiasesLogFile", KeyType::Required, biases_log_file_);

	// Find the bias input packs
	ParameterPack bias_file_pack;
	input_parser.parseFile(biases_log_file_, bias_file_pack);
	auto bias_input_pack_ptrs = bias_file_pack.findParameterPacks("Bias", KeyType::Required);
	int num_biases = bias_input_pack_ptrs.size();
	FANCY_ASSERT( num_biases == num_simulations, "got " << num_biases << " biases, expected " << num_simulations );

	for ( int i=0; i<num_biases; ++i ) {
		// TODO: allow bias log and data summary to report simulations in different orders,
		//      or allow one to be a subset of the others
		// - Use data set labels to match everything
		// TODO: variable T between simulations
		biases_.push_back( Bias(*(bias_input_pack_ptrs[i]), wham_options_.kBT) );
		const auto& label = biases_.back().getDataSetLabel();
		FANCY_ASSERT(label == data_set_labels[i],
			"encountered bias with data label " << label << ", expected " << data_set_labels[i]);
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
		order_parameters_.emplace_back( op_names[p], *(op_pack_ptrs[p]) );
	}


	//---- Initial guess -----//

	f_bias_guess_.assign(num_simulations, 0.0);
	std::string f_bias_guess_file;
	found = input_parameter_pack_.readString("InitialGuess", KeyType::Optional, f_bias_guess_file);
	if ( found ) {
		// TODO: option to read initial guess
	}


	//----- Output Options -----//

	// Check whether the user requested only certain distributions
	// TODO: move to function, generalize to n-dimensional outputs, and (ideally)
	//       make this cleaner overall
	const ParameterPack* output_pack_ptr = 
			input_parameter_pack_.findParameterPack("Output", KeyType::Optional);	
	if ( output_pack_ptr != nullptr ) {
		parseOutputs(*output_pack_ptr);
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

	if ( ! be_quiet_ ) {
		std::cout << "  Solving ...\n" << std::flush;
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
	Wham wham(op_registry_, simulations_, order_parameters_, biases_, f_bias_init, wham_options_.tol);
	solve_wham_timer_.stop();
	f_bias_opt_ = wham.get_f_bias_opt();

	if ( ! be_quiet_ ) {
		std::cout << "  Done\n";
		std::cout << "  Optimal biasing free energies (in kBT):\n";
		for ( int i=0; i<num_simulations; ++i ) {
			std::cout << "  " << i+1 << ":  " << f_bias_opt_[i] << "  (initial: " << f_bias_guess_[i] << ")\n"; 
		}
		std::cout << std::flush;
	}

	// Calculate single-variable outputs
	output_f_x_.reserve( est_f_x_.size() );
	for ( auto& estimator : est_f_x_ ) {
		estimator.calculate(wham);
		output_f_x_.emplace_back( estimator.get_f_x() );
	}


	//----- Estimate Errors -----//

	// TODO: Move this section to a different function?

	const bool calc_error = ( error_method_ != ErrorMethod::None );

	// TODO: MATRIX?
	std::vector<PointEstimator<double>> bootstrap_samples_f_bias;

	if ( error_method_ == ErrorMethod::Bootstrap ) {
		bootstrap_timer_.start();

		if ( ! be_quiet_ ) {
			std::cout << "  Estimate errors using bootstrap subsampling ...\n" << std::flush;
		}

		// Set up bootstrap handlers
		bootstrap_samples_f_bias.resize(num_simulations);
		const int num_output_f_x = est_f_x_.size();
		std::vector<Bootstrap_F_x> bootstrap_f_x;
		bootstrap_f_x.reserve(num_output_f_x);
		for ( int k=0; k<num_output_f_x; ++k ) {
			// Bind the output objects to where errors should ultimately be stored
			bootstrap_f_x.emplace_back( est_f_x_[k].getOrderParameter(),
				output_f_x_[k],	num_bootstrap_samples_);
		}

		// Prepare subsamplers
		std::vector<BootstrapSubsampler> subsamplers;
		std::vector<std::vector<int>> resample_indices(num_simulations);
		for ( int j=0; j<num_simulations; ++j ) {
			int num_samples_j = simulations_[j].getNumSamples();
			// TODO: different seeds depending on debug mode flag
			subsamplers.emplace_back( num_samples_j, Random::getDebugSequence() );
			resample_indices[j].reserve( num_samples_j );
		}

		std::vector<Simulation> bootstrap_simulations(num_simulations);
		std::vector<OrderParameter> bootstrap_ops = order_parameters_;
		for ( int s=0; s<num_bootstrap_samples_; ++s ) {
			// User feedback
			const int stride = 25;
			if ( (! be_quiet_) && ((s == 0) || (((s+1) % stride) == 0)) ) {
				std::cout << "    sample " << s+1 << " of " << num_bootstrap_samples_ << "\n" << std::flush;
			}

			// Subsample each simulation
			// - TODO: OpenMP worth it?
			subsample_timer_.start();
			#pragma omp parallel for
			for ( int j=0; j<num_simulations; ++j ) {
				subsamplers[j].generate( resample_indices[j] );
				bootstrap_simulations[j].setShuffledFromOther( simulations_[j], resample_indices[j] );
			}
			subsample_timer_.stop();

			// Re-solve WHAM equations
			solve_wham_timer_.start();
			Wham bootstrap_wham( op_registry_, bootstrap_simulations, bootstrap_ops, 
													 biases_, f_bias_opt_, wham_options_.tol );
			solve_wham_timer_.stop();

			// Save bootstrap estimates for f_bias_opt
			const auto& bootstrap_f_bias_opt = bootstrap_wham.get_f_bias_opt();
			for ( int j=0; j<num_simulations; ++j ) {
				bootstrap_samples_f_bias[j].addSample( bootstrap_f_bias_opt[j] );
			}

			// Compute and store boostrap estimates of each output
			// - TODO: loop over different estimators
			for ( auto& estimator : bootstrap_f_x ) {
				estimator.addSample(bootstrap_wham);
			}
		};

		// Finalize errors
		error_f_bias_opt_.resize(num_simulations);
		for ( int j=0; j<num_simulations; ++j ) {
			error_f_bias_opt_[j] = bootstrap_samples_f_bias[j].std_dev();
		}
		for ( auto& estimator : bootstrap_f_x ) {
			estimator.finalize();
		}

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
		// TODO: move to function
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
	for ( unsigned i=0; i<est_f_x_.size(); ++i ) {
		const auto& x = est_f_x_[i].getOrderParameter();
		if ( ! be_quiet_ ) {
			std::cout << "Computing F_WHAM(" << x.getName() << ")\n" << std::flush;
		}

		// F_WHAM(x)
		printWhamDistribution(x, output_f_x_[i]);

		// "Manually" unbiased distributions (non-consensus, unshifted)
		ManuallyUnbiasedDistributions f_x_unb(x, wham);
		f_x_unb.print("F_" + x.getName() + "_unbiased.out");

		// "Raw" distributions (i.e. using only data from each individual simulation)
		BiasedDistributions f_x_biased(x, simulations_);
		f_x_biased.print("F_" + x.getName() + "_biased.out");

		// "Rebiased" consensus histograms (for validation)
		RebiasedDistributions f_x_rebiased(x, wham);
		f_x_rebiased.print("F_" + x.getName() + "_rebiased.out");

		// TODO: "shifted" distributions?

		WhamStatistics wham_stats(wham, f_x_biased, f_x_rebiased);
		wham_stats.print();
	}

	// 2-variable outputs
	for ( auto& est_f_x_y : est_f_x_y_ ) {
		// F_WHAM(x,y)
		const OrderParameter& x = est_f_x_y.get_x();
		const OrderParameter& y = est_f_x_y.get_y();
		if ( ! be_quiet_ ) {
			std::cout << "Computing F_WHAM(" << x.getName() << ", " << y.getName() << ")\n" << std::flush;
		}
		est_f_x_y.calculate(wham);

		// Print output files
		est_f_x_y.saveResults();
	}

	print_output_timer_.stop();

	driver_timer_.stop();
}


void WhamDriver::parseOutputs(const ParameterPack& input_pack)
{
	using KeyType = ParameterPack::KeyType;

	bool print_everything = true;
	input_pack.readFlag("PrintAll", KeyType::Optional, print_everything);
	auto vecs = input_pack.findVectors("F_WHAM", KeyType::Optional);

	if ( (! print_everything) || (vecs.size() > 0) ) {
		est_f_x_.clear();
		est_f_x_y_.clear();

		int num_outputs = vecs.size();
		for ( int i=0; i<num_outputs; ++i ) {
			const auto& vec = *(vecs[i]);
			int num_ops_for_output = vec.size();

			// Map names to indices
			// - TODO: direct map to OPs themselves?
			std::vector<int> op_indices(num_ops_for_output);
			for ( int j=0; j<num_ops_for_output; ++j ) {
				op_indices[j] = op_registry_.nameToIndex(vec[j]);
			}

			// Store indices
			if ( num_ops_for_output == 0 ) {
				throw std::runtime_error("F_WHAM with no order parameters is undefined");
			}
			else if ( num_ops_for_output == 1 ) {
				const auto& x = order_parameters_[op_indices[0]];
				est_f_x_.emplace_back(x);
			}
			else if ( num_ops_for_output == 2 ) {
				const auto& x = order_parameters_[op_indices.front()];
				const auto& y = order_parameters_[op_indices.back()];
				est_f_x_y_.emplace_back(x,y);
			}
			else {
				throw std::runtime_error("F_WHAM calcualtions only support up to 2 order parameters");
			}
		}
	}
}


void WhamDriver::printWhamDistribution(
	const OrderParameter& x, const FreeEnergyDistribution& f, std::string file_name) const
{
	// TODO: consistency check for presence of output
	const auto& name = x.getName();
	if ( file_name.empty() ) {
		file_name = "F_" + name + "_WHAM.out";
	}
	std::ofstream ofs(file_name);

	const double temp = simulations_.front().getTemperature();

	// Header
	ofs << "# Consensus free energy distributions from WHAM: \n"
      << "#   F(" << name << ") [in k_B*T] with T = " << temp << " K\n";

	// Table header
	// FIXME: this should really be a part of the FreeEnergyDistribution class, since
	//        it controls how the distribution is printed
	ofs << "# " << name << "\tF[kBT]  NumSamples";
	if ( f.hasErrors() ) {
		ofs << "\tstderr(F)[kBT]";
	}
	ofs << "\n";

	ofs << f;
	ofs.close();
}