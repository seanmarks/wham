#include "WhamDriver.h"

WhamDriver::WhamDriver(const std::string& options_file):
	options_file_(options_file),
	// Columns in data summary file
	// TODO add input options to change these
	col_data_label_(0),
	col_t_min_(3), col_t_max_(col_t_min_+1), col_T_(-1)
{

	// Read input file into a ParameterPack
	InputParser input_parser;
	input_parser.parseFile(options_file_, input_parameter_pack_);

	using KeyType = ParameterPack::KeyType;

	input_parameter_pack_.readNumber("Temperature", KeyType::Required, wham_options_.T);
	wham_options_.kBT = K_B_ * wham_options_.T;

	input_parameter_pack_.readNumber("SolverTolerance", KeyType::Required, wham_options_.tol);

	// Read the data summary, which contains:
	//   (1) data set labels
	//   (2) production phase
	// This sets the number of simulations expected in all files subsequently parsed
	input_parameter_pack_.readString("DataSummaryFile", KeyType::Required, data_summary_file_);
	readDataSummary(data_summary_file_);
	int num_simulations = simulations_.size();


	//----- Order Parameters -----//

	// Register order parameters
	// - When they are constructed, they read in their own time series data
	std::vector<const ParameterPack*> op_pack_ptrs = 
			input_parameter_pack_.findParameterPacks("OrderParameter", KeyType::Required);
	order_parameters_.clear();
	int num_ops = op_pack_ptrs.size();
	if ( num_ops < 1 ) {
		throw std::runtime_error("no order parameters were registered");
	}
	for ( int p=0; p<num_ops; ++p ) {
		order_parameters_.push_back( OrderParameter(*(op_pack_ptrs[p]), simulations_, wham_options_.floor_t) );

		// Record the mapping from name to index
		auto ret = map_op_names_to_indices_.insert( std::make_pair( order_parameters_.back().name_, p ) );
		if ( ret.second == false ) {
			throw std::runtime_error("order parameter \"" + order_parameters_.back().name_ + "\" was defined more than once");
		}
	}

	// Check time series
	OrderParameter::checkForConsistency( order_parameters_ );

	// Use the time series of the first OP registered to count the number of samples
	// per simulation (and the total)
	num_samples_per_simulation_.resize(num_simulations);
	num_samples_total_ = 0;
	const std::vector<TimeSeries>& sample_time_series = order_parameters_[0].get_time_series();
	for ( int j=0; j<num_simulations; ++j ) {
		num_samples_per_simulation_[j] = sample_time_series[j].size();
		num_samples_total_             += num_samples_per_simulation_[j];
	}
	inv_num_samples_total_ = 1.0/static_cast<double>(num_samples_total_);

	// Fraction of total number of samples contributed by each simulation
	c_.resize(num_simulations);
	log_c_.resize(num_simulations);
	for ( int r=0; r<num_simulations; ++r ) {
		c_[r] = static_cast<double>(num_samples_per_simulation_[r]) * inv_num_samples_total_;
		log_c_[r] = log( c_[r] );
	}

	// For convenience, store the ranges corresponding to each simulation's data 
	// in the linearized arrays (e.g. 'u_bias_as_other')
	simulation_data_ranges_.resize(num_simulations);
	int first, end;
	for ( int j=0; j<num_simulations; ++j ) {
		if ( j == 0 ) {
			first = 0;
		}
		else {
			first = simulation_data_ranges_[j-1].second;
		}
		end = first + num_samples_per_simulation_[j];

		simulation_data_ranges_[j] = std::make_pair(first, end);
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
		// or allow one to be a subset of the others
		// - Use data set labels to match everything

		biases_.push_back( Bias(*(bias_input_pack_ptrs[i]), simulations_[i].kBT) );

		if ( simulations_[i].data_set_label != biases_.back().get_data_set_label() ) {
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
#endif /* DEBUG */

	/*
	// For each sample, evaluate the resulting potential under each bias
	// - Populates 'u_bias_as_other_'
	evaluateBiases();
	*/

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
				std::vector<int> op_indices;
				for ( int j=0; j<num_ops_for_output; ++j ) {
					auto pair_it = map_op_names_to_indices_.find( vec[j] );
					if ( pair_it != map_op_names_to_indices_.end() ) {
						op_indices.push_back( pair_it->second );
					}
					else {
						throw std::runtime_error( "Error setting up F_WHAM outputs: order parameter \'" + 
																			vec[j] + "\' is not registered" );
					}
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
			std::cout << "    file = " << order_parameters_[i].time_series_[j].get_file() << "\n"
			          << "    size = " << order_parameters_[i].time_series_[j].size() << " data points\n";
		}
	}
#endif /* DEBUG */


	/* FIXME
	// Manually unbias each OP
	for ( int p=0; p<num_ops; ++p ) {
		manuallyUnbiasDistributions( order_parameters_[p] );
	}
	*/
}


void WhamDriver::readDataSummary(const std::string& data_summary_file)
{
	simulations_.clear();

	// Largest column index of interest
	int max_col = std::max( { col_data_label_, col_t_min_, col_T_ }, std::less<int>() );

	// Get the location of the data summary file
	std::string data_summary_path = FileSystem::get_basename(data_summary_file);

	std::ifstream ifs(data_summary_file);
	if ( not ifs.is_open() ) {
		throw std::runtime_error("Unable to open data summary file: " + data_summary_file + "\n");
	}

	std::string line, token, rel_path;
	std::vector<std::string> tokens;

	while ( getline(ifs, line) ) {
		std::stringstream ss(line);
		ss >> token;

		if ( line.empty() or token[0] == '#' ) {
			continue;  // skip empty lines and comments
		}
		
		// Tokenize the line
		tokens = {{ token }};
		while ( ss >> token ) {
			tokens.push_back( token );
		}

		// Check that the appropriate number of tokens was parsed
		int num_tokens = tokens.size();
		if ( max_col >= num_tokens ) {
			throw std::runtime_error("a line in the data summary has too few columns");
		}

		// Create a new simulation
		simulations_.push_back( Simulation() );
		Simulation& new_sim = simulations_.back();  // convenient alias

		// Data set label
		new_sim.data_set_label = tokens[col_data_label_];

		// Production phase
		new_sim.t_min = std::stod( tokens[col_t_min_] );
		new_sim.t_max = std::stod( tokens[col_t_max_] );

		// Temperature (optional; defaults to value from WHAM options file)
		double kBT = wham_options_.kBT;
		if ( col_T_ >= 0 ) {
			simulations_.back().kBT = K_B_ * std::stod(tokens[col_T_]);
		}
		new_sim.kBT  = kBT;
		new_sim.beta = 1.0/kBT;
	}
}


/*
void WhamDriver::evaluateBiases()
{
	// First, figure out which order parameters are needed for each bias
	int num_simulations = simulations_.size();
	int num_biases      = biases_.size();
	std::vector<std::vector<int>> op_indices(num_simulations);
	for ( int r=0; r<num_biases; ++r ) {
		// Names of order parameters needed
		const std::vector<std::string>& op_names = biases_[r].get_order_parameter_names();
		int num_ops = op_names.size();

		// Search order parameter registry for matches
		op_indices[r].resize(num_ops);
		for ( int p=0; p<num_ops; ++p ) {
			auto pair_it = map_op_names_to_indices_.find( op_names[p] );
			if ( pair_it != map_op_names_to_indices_.end() ) {
				op_indices[r][p] = pair_it->second;
			}
			else {
				throw std::runtime_error( "Error evaluating biases: order parameter \'" + 
				                          op_names[p] + "\' is not registered" );
			}
		}
	}

	// Allocate memory
	u_bias_as_other_.resize( num_biases );
	for ( int r=0; r<num_biases; ++r ) {
		u_bias_as_other_[r].resize( num_samples_total_ );
	}
	u_bias_as_other_unbiased_.assign( num_samples_total_, 0.0 );

	// Now, for each biasing potential, evaluate the value of the bias that *would*
	// result if that sample were obtained from that biased simulation
	std::vector<double> args;
	for ( int r=0; r<num_biases; ++r ) {
		// Number of OPs involved in this bias
		int num_ops_for_bias = op_indices[r].size();
		args.resize( num_ops_for_bias );

		// Evaluate bias for each sample
		int sample_index = 0;
		for ( int j=0; j<num_simulations; ++j ) {
			int num_samples = num_samples_per_simulation_[j];
			for ( int i=0; i<num_samples; ++i ) {
				// Pull together the order parameter values needed for this bias
				for ( int s=0; s<num_ops_for_bias; ++s ) {
					args[s] = order_parameters_[ op_indices[r][s] ].time_series_[j][i];
				}

				u_bias_as_other_[r][sample_index] = biases_[r].evaluate( args );
				++sample_index;
			}
		}
	}
}
*/


/*
// After reading input files, use to generate histograms of raw data from biased simulations
// - TODO Move as much as possible to OrderParameter
void WhamDriver::manuallyUnbiasDistributions(OrderParameter& x)
{
	std::vector<double> u_bias_tmp;
	int num_simulations = static_cast<int>( simulations_.size() );

	for ( int i=0; i<num_simulations; ++i ) {
		const TimeSeries& samples = x.time_series_[i];
		int num_samples = samples.size();

		// Assemble the biasing potential values for this simulation's data
		// under its own bias
		u_bias_tmp.resize(num_samples);
		int index = simulation_data_ranges_[i].first;  //
		for ( int j=0; j<num_samples; ++j ) {
			u_bias_tmp[j] = u_bias_as_other_[i][index];
			++index;
		}

		// Unbiased results, using only data from this simulation
		manually_unbias_f_x( 
				samples, u_bias_tmp, simulations_[i].f_bias_guess, x.get_bins(),
				x.unbiased_distributions_[i] );
	}
}
*/


void WhamDriver::run_driver()
{
	int num_simulations = simulations_.size();
	bool be_verbose = true;

	//----- Solve WHAM equations -----//

	// Initial guess
	// - Shift so that f[0] = 0
	std::vector<double> f_init(num_simulations, 0.0);
	double f_shift = f_bias_guess_[0];
	for ( int j=0; j<num_simulations; ++j ) {
		f_init[j] = f_bias_guess_[j] - f_shift;
	}

	// Run solver
	// TODO errors

	// FIXME
	Wham wham(simulations_, order_parameters_, biases_, wham_options_.tol);
	wham.solveWhamEquations(f_init, f_bias_opt_);
	/*
	solveWhamEquations(f_init, f_bias_opt_);
	*/

	// Print optimal biasing free energies to file
	std::string file_name = "f_bias_WHAM.out";
	std::ofstream ofs(file_name);
	ofs << "# F_bias [k_B*T] for each window after minimization\n";
	for ( int i=0; i<num_simulations; ++i ) {
		ofs << std::setprecision(7) << f_bias_opt_[i] << "\n";
	}
	ofs.close(); ofs.clear();


	//----- Output -----//

	// FIXME

	/*
	// TODO Move elsewhere?
	// 1-variable outputs
	for ( unsigned i=0; i<output_f_x_.size(); ++i ) {
		OrderParameter& x = order_parameters_[ output_f_x_[i] ];
		if ( be_verbose ) {
			std::cout << "Computing F_WHAM(" << x.name_ << ")\n";
		}

		// "Raw" distributions (i.e. using only data from each individual simulation)
		x.printRawDistributions();

		// F_WHAM(x)
		compute_consensus_f_x( 
			x.time_series_, u_bias_as_other_, f_bias_opt_, u_bias_as_other_unbiased_, f_unbiased_, x.bins_,
			// Output
			x.wham_distribution_
		);

		// TODO errors
		int num_bins_x = x.bins_.get_num_bins();
		x.wham_distribution_.error_f_x.assign(num_bins_x, 0.0);

		// "Rebias" consensus histograms (for validation)
		for ( int k=0; k<num_simulations; ++k ) {
			compute_consensus_f_x( 
				x.time_series_, u_bias_as_other_, f_bias_opt_, u_bias_as_other_[k], f_bias_opt_[k], x.bins_,
				// Output
				x.rebiased_distributions_[k] );
		}

		// Relative entropy between distributions (aka Kullback-Leibler divergence)
		// - Closely related to Shannon entropy
		// - See Hummer and Zhu (J Comp Chem 2012), Eqn. 2
		// TODO move to function which takes two arbitrary distributions
		x.info_entropy_.assign(num_simulations, 0.0);
		double bin_size_x = x.bins_.get_bin_size();
		for ( int i=0; i<num_simulations; ++i ) {
			for ( int b=0; b<num_bins_x; ++b ) {
				const double& p_biased   = x.biased_distributions_[i].p_x[b];
				const double& f_biased   = x.biased_distributions_[i].f_x[b];
				const double& f_rebiased = x.rebiased_distributions_[i].f_x[b];

				if ( p_biased > 0.0 and std::isfinite(f_rebiased) ) {
					x.info_entropy_[i] += p_biased*(f_rebiased - f_biased)*bin_size_x;
				}
			}
		}

		// Print the results for this OP
		printWhamResults(x);
	}

	// 2-variable outputs
	for ( unsigned i=0; i<output_f_x_y_.size(); ++i ) {
		// F_WHAM(x,y)
		OrderParameter& x = order_parameters_[ output_f_x_y_[i][0] ];
		OrderParameter& y = order_parameters_[ output_f_x_y_[i][1] ];
		if ( be_verbose ) {
			std::cout << "Computing F_WHAM(" << x.name_ << ", " << y.name_ << ")\n";
		}

		std::vector<double> p_y_wham, f_y_wham;
		std::vector<int> sample_counts_y;
		std::vector<std::vector<double>> p_x_y_wham, f_x_y_wham; 
		std::vector<std::vector<int>> sample_counts_x_y;

		compute_consensus_f_x_y( 
			x.time_series_, y.time_series_, u_bias_as_other_, f_bias_opt_, u_bias_as_other_unbiased_, f_unbiased_,
			x.bins_, y.bins_,
			// Output
			p_x_y_wham, f_x_y_wham, sample_counts_x_y
		);

		// Print results
		print_f_x_y( x, y, p_x_y_wham, f_x_y_wham, sample_counts_x_y );
	}
	*/
}


/*
// TODO way to merge with compute_consensus_f_x?
void WhamDriver::manually_unbias_f_x(
	const TimeSeries& x, const std::vector<double>& u_bias, const double f,
	const Bins& bins_x,
	// Output
	Distribution& unbiased_distribution_x
) const
{
	// Unpack for readability below
	std::vector<double>& p_x           = unbiased_distribution_x.p_x;
	std::vector<double>& f_x           = unbiased_distribution_x.f_x;
	std::vector<int>&    sample_counts = unbiased_distribution_x.sample_counts;

	unbiased_distribution_x.bins_x = bins_x;

	// Reserve memory
	int num_samples = x.size();
	int num_bins_x = bins_x.get_num_bins();
	minus_log_sigma_k_binned_.resize(num_bins_x);
	for ( int b=0; b<num_bins_x; ++b ) {
		minus_log_sigma_k_binned_[b].resize( 0 );
		minus_log_sigma_k_binned_[b].reserve( num_samples );
	}

	// Compute and sort log(sigma_k)-values by bin
	int bin;
	for ( int i=0; i<num_samples; ++i ) {
		bin = bins_x.find_bin( x[i] );
		if ( bin >= 0 ) {
			minus_log_sigma_k_binned_[bin].push_back( u_bias[i] - f );
		}
	}

	// Compute
	f_x.resize(num_bins_x);
	p_x.resize(num_bins_x);
	sample_counts.resize(num_bins_x);
	double bin_size_x = bins_x.get_bin_size();
	double normalization = log(num_samples*bin_size_x);
	 for ( int b=0; b<num_bins_x; ++b ) {
		sample_counts[b] = minus_log_sigma_k_binned_[b].size();
		if ( sample_counts[b] > 0 ) {
			f_x[b] = normalization - log_sum_exp( minus_log_sigma_k_binned_[b] );
			p_x[b] = exp( -f_x[b] );
		}
		else {
			f_x[b] = 0.0;
			p_x[b] = 0.0;
		}
	}
}
*/


// TODO Move to OrderParameter?
void WhamDriver::printWhamResults(const OrderParameter& x) const
{
	// Working variables
	std::string file_name;
	std::ofstream ofs;
	const int num_simulations = simulations_.size();
	const int num_bins_x      = x.bins_.get_num_bins();

	// Unpack for readability below
	//const std::vector<double>& p_x_wham      = x.wham_distribution_.p_x;
	const std::vector<double>& f_x_wham      = x.wham_distribution_.f_x;
	const std::vector<int>&    sample_counts = x.wham_distribution_.sample_counts;

	// For convenience of visualizing output, shift F(x) so that F=0 at the minimum
	auto f_x_wham_shifted = Distribution::shift_f_x_to_zero( f_x_wham, sample_counts );

	// Print F_0(x)
	file_name = "F_" + x.name_ + "_WHAM.out";
	ofs.open(file_name);
	ofs << "# Consensus free energy distributions from WHAM: \n"
      << "#   F(" << x.name_ << ") [in k_B*T] with T = " << wham_options_.T << " K\n";
	ofs << "# " << x.name_ << "\tF[kBT]  NumSamples\n";  //"\t" << "\tstderr(F)\n"; TODO error estimate
	for ( int b=0; b<num_bins_x; ++b ) {
		ofs << std::setw(8) << std::setprecision(5) << x.bins_[b] << "\t";
		ofs << std::setw(8) << std::setprecision(5);
			Distribution::print_free_energy(ofs, f_x_wham_shifted[b], sample_counts[b]);
		ofs << std::setw(8) << std::setprecision(5) << sample_counts[b];
		//<< std::setw(8) << std::setprecision(5) << wham_results_.error_f[b]
		ofs << "\n";
	}
	ofs.close(); ofs.clear();

	// "Rebiased" free energy distributions
	file_name = "F_" + x.name_ + "_rebiased.out";

	std::stringstream header_stream;
	header_stream << "# \"Rebiased\" free energy distributions: "
                << " F_{rebias,i}(" << x.name_ << ") [k_B*T]\n";
	header_stream << "# Data sets (by column)\n";
	for ( int i=0; i<num_simulations; ++i ) {
		header_stream << "# " << i+2 << ": " << simulations_[i].data_set_label << "\n";
	}
	header_stream << "#\n"
	              << "# " << x.name_ << " | F(" << x.name_ << ") [kBT]\n";

	x.printDistributions( x.rebiased_distributions_, file_name, header_stream.str() );

	// Misc. stats
	file_name = "stats_" + x.name_ + ".out";
	ofs.open(file_name);
	ofs << "# data_set   avg(x)   var(x)   info_entropy(biased/rebiased)\n";
	for ( int i=0; i<num_simulations; ++i ) {
		ofs << simulations_[i].data_set_label << "\t"
		    << x.time_series_[i].average() << "\t"
		    << x.time_series_[i].variance() << "\t"
		    << x.info_entropy_[i] << "\n";
	}
	ofs.close(); ofs.clear();
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
	const Bins& bins_x = x.bins_;
	const Bins& bins_y = y.bins_;
	int num_bins_x = bins_x.get_num_bins();
	int num_bins_y = bins_y.get_num_bins();
	std::string file_name, sep;
	std::ofstream ofs;

	// Print bins
	std::vector<const OrderParameter*> op_ptrs = {{ &x, &y }};
	for ( unsigned i=0; i<op_ptrs.size(); ++i ) {
		file_name = "bins_" + op_ptrs[i]->name_ + ".out";
		ofs.open(file_name);
		ofs << "# Bins for " << op_ptrs[i]->name_ << " for F_WHAM(" << x.name_ << "," << y.name_ << ")\n"
				<< "# " << op_ptrs[i]->name_ << "\n";
		int num_bins = op_ptrs[i]->bins_.get_num_bins();
		for ( int j=0; j<num_bins; ++j ) {
			ofs << op_ptrs[i]->bins_[j]  << "\n";
		}
		ofs.close(); ofs.clear();
	}

	// Print F(x,y)  (TODO: shift so that F=0 at the minimum)
	file_name = "F_" + x.name_ + "_" + y.name_ + "_WHAM.out";
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
	file_name = "samples_" + x.name_ + "_" + y.name_ + ".out";
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
