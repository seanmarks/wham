#include "Wham.h"

Wham::Wham(const std::string& options_file):
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

	// TODO: avoid this cludgy reformatting
	using Range = std::array<double,2>;
	std::vector<Range> production_phases(num_simulations);
	for ( int i=0; i<num_simulations; ++i ) {
		production_phases[i] = {{ simulations_[i].t_min, simulations_[i].t_max }};
	}

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
	for ( int i=0; i<num_ops; ++i ) {
		order_parameters_.push_back( OrderParameter(*(op_pack_ptrs[i]), production_phases, wham_options_.floor_t) );

		// Record the mapping from name to index
		auto ret = map_op_names_to_indices_.insert( std::make_pair( order_parameters_.back().name_, i ) );
		if ( ret.second == false ) {
			throw std::runtime_error("order parameter \"" + order_parameters_.back().name_ + "\" was defined more than once");
		}
	}

	// Check time series
	OrderParameter::checkForConsistency( order_parameters_ );

	// Register output distributions requested by the input
	// TODO move to function, generalize to n-dimensional outputs, and (ideally)
	// make this cleaner
	bool print_everything = true;
	input_parameter_pack_.readFlag("PrintAll", KeyType::Optional, print_everything);
	if ( print_everything ) {
		// Print all permutations of F(x) and F(x,y)
		for ( int i=0; i<num_ops; ++i ) {
			output_f_x_.push_back( i );
			for ( int j=i+1; j<num_ops; ++j ) {
				output_f_x_y_.push_back( {{i, j}} );
			}
		}
	}
	else {
		std::vector<const std::vector<std::string>*> vecs = 
				input_parameter_pack_.findVectors("F_WHAM", KeyType::Optional);
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

	// For each sample, evaluate the resulting potential under each bias
	evaluateBiases();


	//----- Precompute some constants (TODO: move to function?) -----//

	// Total number of samples
	num_samples_total_ = 0;
	for ( int j=0; j<num_simulations; ++j ) {
		num_samples_total_ += u_bias_as_other_[j].size();
	}
	inv_num_samples_total_ = 1.0/static_cast<double>( num_samples_total_ );

	// Fraction of samples from each simulation
	c_.resize(num_simulations);
	log_c_.resize(num_simulations);
	for ( int r=0; r<num_simulations; ++r ) {
		c_[r] = static_cast<double>(u_bias_as_other_[r].size()) * inv_num_samples_total_;
		log_c_[r] = log( c_[r] );
	}

	// Analyze the raw data for each OP
	for ( int i=0; i<num_ops; ++i ) {
		analyzeRawData( order_parameters_[i] );
	}
}


void Wham::readDataSummary(const std::string& data_summary_file)
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


void Wham::evaluateBiases()
{
	// TODO linearize

	int num_simulations = simulations_.size();
	u_bias_as_other_.resize( num_simulations );

	// First, figure out which order parameters are needed for each bias
	std::vector<std::vector<int>> op_indices(num_simulations);
	for ( int j=0; j<num_simulations; ++j ) {
		// Names of order parameters needed
		const std::vector<std::string>& op_names = biases_[j].get_order_parameter_names();
		int num_ops = op_names.size();

		// Search order parameter registry for matches
		op_indices[j].resize(num_ops);
		for ( int a=0; a<num_ops; ++a ) {
			auto pair_it = map_op_names_to_indices_.find( op_names[a] );
			if ( pair_it != map_op_names_to_indices_.end() ) {
				op_indices[j][a] = pair_it->second;
			}
			else {
				throw std::runtime_error( "Error evaluating biases: order parameter \'" + 
				                          op_names[a] + "\' is not registered" );
			}
		}
	}

	// Now, for each sample, evaluate the bias experienced in each of the simulations
	std::vector<double> args;
	for ( int j=0; j<num_simulations; ++j ) {
		int num_samples = order_parameters_[0].time_series_[j].size();
		u_bias_as_other_[j].resize(num_samples);

		for ( int i=0; i<num_samples; ++i ) {
			u_bias_as_other_[j][i].resize( num_simulations );

			// Evaluate sample i from simulation j using each of the biases r
			for ( int r=0; r<num_simulations; ++r ) {
				// Pull together the order parameter values needed for this bias
				int num_ops_for_bias = op_indices[r].size();
				args.resize( num_ops_for_bias );
				for ( int s=0; s<num_ops_for_bias; ++s ) {
					args[s] = order_parameters_[ op_indices[r][s] ].time_series_[j][i];
				}

				u_bias_as_other_[j][i][r] = biases_[r].evaluate( args );
			}
		}
	}
}



// After reading input files, use to generate histograms of raw data from biased simulations
// - TODO Move to OrderParameter?
void Wham::analyzeRawData(OrderParameter& x)
{
	const Bins& bins = x.bins_;
	int num_bins = bins.get_num_bins();

	// Allocate memory
	int num_simulations = static_cast<int>( simulations_.size() );
	x.p_biased_.resize(num_simulations);
	x.f_biased_.resize(num_simulations);
	x.p_unbiased_.resize(num_simulations);
	x.f_unbiased_.resize(num_simulations);
	x.sample_counts_.resize(num_simulations);

	std::vector<double> u_bias_tmp;
	for ( int i=0; i<num_simulations; ++i ) {
		const TimeSeries& samples = x.time_series_[i];
		int num_samples = samples.size();

		// Assemble the biasing potential values for this simulation
		u_bias_tmp.resize(num_samples);
		for ( int j=0; j<num_samples; ++j ) {
			u_bias_tmp[j] = u_bias_as_other_[i][j][i];
		}

		// Unbiased results, using only data from this simulation
		manually_unbias_f_x( 
				samples, u_bias_tmp, simulations_[i].f_bias_guess, bins,
		    x.p_unbiased_[i], x.f_unbiased_[i], x.sample_counts_[i] );

		// Derive remaining histograms
		double bin_size = bins.get_bin_size();
		double num_samples_d = static_cast<double>( num_samples );
		x.p_biased_[i].resize(num_bins);
		x.f_biased_[i].resize(num_bins);
		for ( int b=0; b<num_bins; ++b ) {
			x.p_biased_[i][b] = (x.sample_counts_[i][b])/(bin_size*num_samples_d);
			x.f_biased_[i][b] = -log( x.p_biased_[i][b] );
		}
	}

	// Construct a histogram of the total number of data points in each bin, across all simulations
	x.global_sample_counts_.assign(num_bins, 0);
	for ( int i=0; i<num_simulations; ++i ) {
		for ( int b=0; b<num_bins; ++b ) {
			x.global_sample_counts_[b] += x.sample_counts_[i][b];
		}
	}
}


void Wham::run_driver()
{
	int num_simulations = simulations_.size();
	bool be_verbose = true;

	//----- Solve WHAM equations -----//

	// Initial guess
	// - Shift so that f[0] = 0
	std::vector<double> f_init(num_simulations, 0.0);
	double f_shift = simulations_[0].f_bias_guess;
	for ( int j=0; j<num_simulations; ++j ) {
		f_init[j] = simulations_[j].f_bias_guess - f_shift;
	}

	// Run solver
	// TODO errors
	std::vector<double> f_opt;
	solveWhamEquations(f_init, f_opt);
	wham_results_.f_opt = f_opt;

	// Print optimal biasing free energies to file
	std::string file_name = "f_bias_WHAM.out";
	std::ofstream ofs(file_name);
	ofs << "# F_bias [k_B*T] for each window after minimization\n";
	for ( int i=0; i<num_simulations; ++i ) {
		ofs << std::setprecision(7) << f_opt[i] << "\n";
	}
	ofs.close(); ofs.clear();


	//----- Output -----//

	// TODO Move elsewhere?
	// F_WHAM(x)
	for ( unsigned i=0; i<output_f_x_.size(); ++i ) {
		OrderParameter& x = order_parameters_[ output_f_x_[i] ];
		if ( be_verbose ) {
			std::cout << "Computing F_WHAM(" << x.name_ << ")\n";
		}

		// "Raw" distributions (i.e. using only data from each individual simulation)
		printRawDistributions(x);

		// Compute F_WHAM(x)
		compute_consensus_f_x( 
			x.time_series_, u_bias_as_other_, f_opt, unbiased_ensemble_index_, x.bins_,
			// Output
			x.p_x_wham_, x.f_x_wham_, x.global_sample_counts_
		);

		// TODO errors
		int num_bins_x = x.bins_.get_num_bins();
		x.error_f_x_wham_.assign(num_bins_x, 0.0);

		// "Rebias" consensus histograms (for validation)
		std::vector<int> sample_counts_tmp;
		x.p_rebiased_.resize(num_simulations);
		x.f_rebiased_.resize(num_simulations);
		for ( int k=0; k<num_simulations; ++k ) {
			compute_consensus_f_x( 
				x.time_series_, u_bias_as_other_, f_opt, k, x.bins_,
				// Output
			  x.p_rebiased_[k], x.f_rebiased_[k], sample_counts_tmp 
			);
		}

		// Relative entropy between distributions (aka Kullback-Leibler divergence)
		// - Closely related to Shannon entropy
		// - See Hummer and Zhu (J Comp Chem 2012), Eqn. 2
		// TODO move to function which takes two arbitrary distributions
		x.info_entropy_.assign(num_simulations, 0.0);
		double bin_size_x = x.bins_.get_bin_size();
		for ( int i=0; i<num_simulations; ++i ) {
			for ( int b=0; b<num_bins_x; ++b ) {
				const double& p_biased   = x.p_biased_[i][b];
				const double& f_biased   = x.f_biased_[i][b];
				const double& f_rebiased = x.f_rebiased_[i][b];

				if ( p_biased > 0.0 and std::isfinite(f_rebiased) ) {
					x.info_entropy_[i] += p_biased*(f_rebiased - f_biased)*bin_size_x;
				}
			}
		}

		// Print the results for this OP
		printWhamResults(x);
	}


	// F_WHAM(x,y)
	for ( unsigned i=0; i<output_f_x_y_.size(); ++i ) {
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
			x.time_series_, y.time_series_, u_bias_as_other_, f_opt, unbiased_ensemble_index_, 
			x.bins_, y.bins_,
			// Output
			p_x_y_wham, f_x_y_wham, sample_counts_x_y
		);

		// Print results
		print_f_x_y( x, y, p_x_y_wham, f_x_y_wham, sample_counts_x_y );


		// FIXME DEBUG
		int num_samples_x = std::accumulate( x.global_sample_counts_.begin(), x.global_sample_counts_.end(), 0);
		int num_samples_y = std::accumulate( y.global_sample_counts_.begin(), y.global_sample_counts_.end(), 0);

		int num_samples_x_y = 0;
		for ( unsigned i=0; i<sample_counts_x_y.size(); ++i ) {
			num_samples_x_y += std::accumulate( sample_counts_x_y[i].begin(), sample_counts_x_y[i].end(), 0 );
		}

		std::cout << "  num_samples_x = " << num_samples_x << "\n"
		          << "  num_samples_y = " << num_samples_y << "\n"
		          << "  num_samples_x_y = " << num_samples_x_y << "\n";
		// END FIXME DEBUG
	}
}


// TODO add toggle for be_verbose
void Wham::solveWhamEquations(const std::vector<double>& f_init, std::vector<double>& f_opt)
{
	bool be_verbose = true;

	// Compute differences btw. successive windows for initial guess
	int num_simulations = simulations_.size();
	int num_vars = num_simulations - 1;
	Wham::ColumnVector df(num_vars), df_init(num_vars);
	convert_f_to_df(f_init, df);
	df_init = df;  // save initial guess

	// Invoke functor wrappers for dlib calls
	WhamDlibEvalWrapper  evaluate_wrapper(*this);
	WhamDlibDerivWrapper derivatives_wrapper(*this);

	// Call dlib to perform the optimization
	double min_A = dlib::find_min( 
			dlib::bfgs_search_strategy(), 
			dlib::objective_delta_stop_strategy(wham_options_.tol).be_verbose(),
			evaluate_wrapper,
			derivatives_wrapper,
			df,
			std::numeric_limits<double>::lowest() 
	);

	// Compute optimal free energies of biasing from the differences between windows
	f_opt.resize(num_simulations);
	convert_df_to_f(df, f_opt);
	wham_results_.f_opt = f_opt;

	// Report results 
	if ( be_verbose ) {
		std::cout << "min[A] = " << min_A << "\n";
		std::cout << "Biasing free energies:\n";
		for ( int i=0; i<num_simulations; ++i ) {
			std::cout << i << ":  " << f_opt[i] << "  (initial: " << f_init[i] << ")\n"; 
		}
	}
}


void Wham::convert_f_to_df(const std::vector<double>& f, Wham::ColumnVector& df) const
{
	int num_differences = static_cast<int>( f.size() ) - 1;
	df.set_size(num_differences);
	for ( int i=0; i<num_differences; ++i ) {
		df(i) = f[i+1] - f[i];
	}	
}


void Wham::convert_df_to_f(const Wham::ColumnVector& df, std::vector<double>& f) const
{
	int num_simulations = static_cast<int>( df.size() + 1 );
	f.resize(num_simulations);
	f[0] = 0.0; // freely specified
	for ( int k=1; k<num_simulations; ++k ) {
		f[k] = f[k-1] + df(k-1);
	}	
}


// Evaluate objective function
double Wham::evalObjectiveFunction(const Wham::ColumnVector& df) const
{
	// Compute free energies of biasing from the differences between windows
	int num_simulations = static_cast<int>( simulations_.size() );
	std::vector<double> f(num_simulations);
	convert_df_to_f(df, f);

	// Compute log(sigma)
	// - These values are used here, and are stored for use when computing
	//   the derivatives later
	compute_log_sigma( f, u_bias_as_other_, unbiased_ensemble_index_, log_sigma_0_ );
	f_bias_last_ = f;

	// Compute the objective function's value
	double val = 0.0;
	for ( int j=0; j<num_simulations; ++j ) {
		int num_samples = log_sigma_0_[j].size();
		for ( int i=0; i<num_samples; ++i ) {
			val += log_sigma_0_[j][i];
		}
	}
	val *= inv_num_samples_total_;
	for ( int r=0; r<num_simulations; ++r ) {
		val -= c_[r]*f[r];
	}

#ifdef DEBUG_MODE
	std::cout << "EVALUATED: A = " << val << "\n";
#endif

	return val;
}


void Wham::compute_log_sigma(
		const std::vector<double>& f,
		const std::vector<std::vector<std::vector<double>>>& u_bias_as_other,
		const int k,
		std::vector<std::vector<double>>& log_sigma ) const
{
	int num_simulations = u_bias_as_other.size();
	log_sigma.resize(num_simulations);

	// Working variable
	std::vector<double> args(num_simulations);

	// TODO linearize
	for ( int j=0; j<num_simulations; ++j ) {
		int num_samples = u_bias_as_other[j].size();
		log_sigma[j].resize(num_samples);

		for ( int i=0; i<num_samples; ++i ) {
			// Compute log(sigma_{j,i}) for this sample
			for ( int r=0; r<num_simulations; ++r ) {
				args[r] = log_c_[r] - u_bias_as_other[j][i][r] + f[r];

				if ( k >= 0 ) {
					// Reweighting to a different ensemble
					args[r] += u_bias_as_other[j][i][k] - f[k];
				}
			}
			log_sigma[j][i] = log_sum_exp( args );
		}
	}
}


// Compute derivatives of objective function
// - Note:
//   - grad  = gradient wrt. free energy differences between neighboring windows
//   - dA_df = gradient wrt. biasing free energies themselves
const Wham::ColumnVector Wham::evalObjectiveDerivatives(const Wham::ColumnVector& df) const
{
	// Compute free energies of biasing from differences between windows
	int num_simulations = static_cast<int>( simulations_.size() );
	std::vector<double> f(num_simulations);
	convert_df_to_f(df, f);

	// Recompute log(sigma) as necessary
	if ( f != f_bias_last_ ) {
		compute_log_sigma( f, u_bias_as_other_, unbiased_ensemble_index_, log_sigma_0_ );
		f_bias_last_ = f;
	}

	// First, compute derivatives of the objective fxn (A) wrt. the biasing free
	// energies themselves (f)
	// TODO: expensive nested loops; better way? perhaps linearize
	Wham::ColumnVector dA_df(num_simulations);
	dA_df(0) = 0.0;
	args_buffer_.resize(num_samples_total_); 
	double log_sum_exp_args;
	for ( int k=1; k<num_simulations; ++k ) {
		// Big log_sum_exp (TODO: better way?)
		int index = 0;  // linear index
		for ( int j=0; j<num_simulations; ++j ) {
			int num_samples = log_sigma_0_[j].size();
			for ( int i=0; i<num_samples; ++i ) {
				args_buffer_[index] = log_c_[k] - u_bias_as_other_[j][i][k] + f[k] - log_sigma_0_[j][i];
				++index;
			}
		}
		log_sum_exp_args = log_sum_exp(args_buffer_);

		dA_df(k) = inv_num_samples_total_*exp(log_sum_exp_args) - c_[k];
	}

	// Now compute the gradient of the objective function with respect to the 
	// actual free variables: the biasing free energy differences btw. 
	// successive windows
	int num_vars = static_cast<int>( df.size() );
	Wham::ColumnVector grad(num_vars);
	for ( int i=0; i<num_vars; ++i ) {
		grad(i) = 0.0;
		for ( int j=i+1; j<num_simulations; ++j ) {
			grad(i) += dA_df(j);
		}
	}

#ifdef DEBUG_MODE
	std::cout << "dA_dj =";
	for ( int j=0; j<num_simulations; ++j ) {
		std::cout << " " << dA_df(j);
	}
	std::cout << "\n";

	std::cout << "grad =";
	for ( int j=0; j<num_vars; ++j ) {
		std::cout << " " << grad(j);
	}
	std::cout << "\n\n";
#endif

	return grad;
}


// Returns the logarithm of a sum of exponentials (input: arguments of exponentials)
double Wham::log_sum_exp(const std::vector<double>& args) const
{
	// Check input
	int num_args = static_cast<int>( args.size() );
	if ( num_args < 1 ) {
		throw std::runtime_error("Wham::log_sum_exp: No arguments supplied.\n");
	}

	// Find the largest argument
	double max_arg = args[0];
	for ( int i=1; i<num_args; ++i ) {
		if ( args[i] > max_arg ) { max_arg = args[i]; }
	}

	// Compute sum of exp(args[i] - max_arg)
	double sum = 0.0;
	for ( int i=0; i<num_args; ++i ) {
		sum += exp(args[i] - max_arg);
	}

	// Return log of desired summation
	return max_arg + log(sum);
}


// TODO way to merge with compute_consensus_f_x?
void Wham::manually_unbias_f_x(
	const TimeSeries& x, const std::vector<double>& u_bias, const double f,
	const Bins& bins_x,
	// Output
	std::vector<double>& p_x, std::vector<double>& f_x, std::vector<int>& sample_counts
) const
{
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


void Wham::compute_consensus_f_x(
	const std::vector<TimeSeries>& x, 
	const std::vector<std::vector<std::vector<double>>>& u_bias_as_other,
	const std::vector<double>& f_opt, const int k, const Bins& bins_x,
	// Output
	std::vector<double>& p_x_wham, std::vector<double>& f_x_wham, 
	std::vector<int>& sample_counts
) const
{
	// Compute log(sigma_k) using the buffer, since the array can be quite large
	compute_log_sigma(f_opt, u_bias_as_other, k, log_sigma_k_);

	// Reserve some memory for binning samples
	int num_bins_x = bins_x.get_num_bins();
	minus_log_sigma_k_binned_.resize(num_bins_x);
	for ( int b=0; b<num_bins_x; ++b ) {
		minus_log_sigma_k_binned_[b].resize( 0 );
		minus_log_sigma_k_binned_[b].reserve( x[0].size() );
	}

	// Sort log(sigma_k)-values by bin
	int bin, num_simulations = f_opt.size();
	double bin_size_x = bins_x.get_bin_size();
	for ( int j=0; j<num_simulations; ++j ) {
		int num_samples = x[j].size();
		for ( int i=0; i<num_samples; ++i ) {
			bin = bins_x.find_bin( x[j][i] );
			if ( bin >= 0 ) {
				minus_log_sigma_k_binned_[bin].push_back( -log_sigma_k_[j][i] );
			}
		}
	}

	// Compute
	f_x_wham.resize(num_bins_x);
	p_x_wham.resize(num_bins_x);
	sample_counts.resize(num_bins_x);
	double log_sum_exp_bin;
	for ( int b=0; b<num_bins_x; ++b ) {
		sample_counts[b] = minus_log_sigma_k_binned_[b].size();
		if ( sample_counts[b] > 0 ) {
			log_sum_exp_bin = log_sum_exp( minus_log_sigma_k_binned_[b] );
			f_x_wham[b] = -log(inv_num_samples_total_/bin_size_x) - log_sum_exp_bin;
			p_x_wham[b] = exp( -f_x_wham[b] );
		}
		else {
			f_x_wham[b] = 0.0;
			p_x_wham[b] = 0.0;
		}
	}
}


void Wham::compute_consensus_f_x_y(
	const std::vector<TimeSeries>& x, const std::vector<TimeSeries>& y,
	const std::vector<std::vector<std::vector<double>>>& u_bias_as_other,
	const std::vector<double>& f_opt, const int k, const Bins& bins_x, const Bins& bins_y,
	// Consensus distributions for F_k(x,y)
	std::vector<std::vector<double>>& p_x_y_wham, std::vector<std::vector<double>>& f_x_y_wham,
	std::vector<std::vector<int>>& sample_counts_x_y
) const
{
	// Input checks
	if ( y.size() != x.size() or f_opt.size() != x.size() or u_bias_as_other.size() != x.size() ) {
		throw std::runtime_error("Error in Wham::compute_consensus_f_x_y - Array size inconsistency");
	}

	// Allocate memory and set up grids
	int num_bins_x = bins_x.get_num_bins();
	p_x_y_wham.resize( num_bins_x );
	f_x_y_wham.resize( num_bins_x );
	sample_counts_x_y.resize( num_bins_x );
	int num_bins_y = bins_y.get_num_bins();
	for ( int a=0; a<num_bins_x; ++a ) {
		// Second dimension
		p_x_y_wham[a].resize( num_bins_y );
		f_x_y_wham[a].resize( num_bins_y );
		sample_counts_x_y[a].assign( num_bins_y, 0 );
	}

	// Normalization factors
	double bin_size_x     = bins_x.get_bin_size();
	double log_bin_size_x = log(bin_size_x);

	// Compute F_k(x = x_a, y) for each bin x_a
	int num_simulations = f_opt.size();
	std::vector<std::vector<double>> y_given_x_a(num_simulations);
	std::vector<std::vector<std::vector<double>>> u_bias_as_other_given_x_a(num_simulations);
	for ( int a=0; a<num_bins_x; ++a ) {
		// Reset working arrays
		for ( int j=0; j<num_simulations; ++j ) {
			y_given_x_a[j].resize(0);
			u_bias_as_other_given_x_a[j].resize(0);
		}

		// Gather all the y-samples that fall into this x-bin, and the corresponding
		// values of the biasing potential in each ensemble
		for ( int j=0; j<num_simulations; ++j ) {
			int num_samples = x[j].size();
			for ( int i=0; i<num_samples; ++i ) {
				// Check whether this samples lies in x-bin 'a'
				if ( bins_x.find_bin(x[j][i]) == a ) {
					y_given_x_a[j].push_back( y[j][i] );
					u_bias_as_other_given_x_a[j].push_back( u_bias_as_other[j][i] );
				}
			}
		}

		// TODO: avoid this cludgy reformatting
		std::vector<TimeSeries> time_series_tmp;
		for ( int j=0; j<num_simulations; ++j ) {
			time_series_tmp.push_back( TimeSeries(y_given_x_a[j]) );
		}

		// Compute consensus distribution for this slice (note: here, 'x' is actually y)
		compute_consensus_f_x( time_series_tmp, u_bias_as_other_given_x_a, f_opt, k, bins_y,
		                       p_x_y_wham[a], f_x_y_wham[a], sample_counts_x_y[a] );
		for ( int b=0; b<num_bins_y; ++b ) {
			// Include normalization for x (the above function only normalizes for bin_size_y)
			p_x_y_wham[a][b] /= bin_size_x;
			f_x_y_wham[a][b] += log_bin_size_x;
		}
	}
}


void Wham::printRawDistributions(const OrderParameter& x) const
{
	int num_simulations = static_cast<int>( simulations_.size() );
	if ( num_simulations < 1 ) {
		throw std::runtime_error("Wham::printRawDistributions: No simulations recorded.\n");
	}

	// Common header
	std::stringstream header_stream;

	std::stringstream table_header_stream;
	table_header_stream << "# Data sets (by column)\n";
	for ( int i=0; i<num_simulations; ++i ) {
		table_header_stream << "# " << i+2 << ": " << simulations_[i].data_set_label << "\n";
	}
	table_header_stream << "#\n"
	                    << "# " << x.name_ << " | F(" << x.name_ << ") [kBT]\n";

	// Working variables
	std::string file_name;
	std::ofstream ofs;
	const Bins& bins = x.bins_;
	int num_bins = bins.get_num_bins();  // All simulations are binned the same way


	//----- Print biased free energy distributions -----//

	file_name = "F_" + x.name_ + "_biased.out";
	ofs.open( file_name );
	ofs << "# Biased free energy distributions: "
      << " F_i(" << x.name_ << ") [k_B*T]\n";
	ofs << header_stream.str() << table_header_stream.str();
	for ( int b=0; b<num_bins; ++b ) {
		ofs << bins[b];
		for ( int i=0; i<num_simulations; ++i ) {
			ofs << "\t";
			ofs << std::setw(8) << std::setprecision(5) << x.f_biased_[i][b];
		}
		ofs << "\n";
	}
	ofs.close();


	//----- Print unbiased free energy distributions (non-consensus) -----//

	file_name = "F_" + x.name_ + "_unbiased.out";
	ofs.open( file_name );
	ofs << "# Unbiased free energy distributions: "
      << " F_{0,i}(" << x.name_ << ") [k_B*T]\n";
	ofs << header_stream.str() << table_header_stream.str();
	for ( int b=0; b<num_bins; ++b ) {
		ofs << bins[b];
		for ( int k=0; k<num_simulations; ++k ) {
			ofs << "\t" << std::setw(8) << std::setprecision(5);
			print_free_energy(ofs, x.f_unbiased_[k][b], x.p_unbiased_[k][b]);
		}
		ofs << "\n";
	}
	ofs.close();
}


void Wham::printWhamResults(const OrderParameter& x) const
{
	// Working variables
	std::stringstream header_stream;
	std::string file_name;
	std::ofstream ofs;
	const int num_simulations = simulations_.size();
	const int num_bins_x      = x.bins_.get_num_bins();

	// Print F_0(x)
	file_name = "F_" + x.name_ + "_WHAM.out";
	ofs.open(file_name);
	ofs << "# Consensus free energy distributions from WHAM: \n"
      << "#   F(" << x.name_ << ") [in k_B*T] with T = " << wham_options_.T << " K\n";
	ofs << "# " << x.name_ << "\tF[kBT]  NumSamples\n";  //"\t" << "\tstderr(F)\n"; TODO error estimate
	for ( int b=0; b<num_bins_x; ++b ) {
		ofs << std::setw(8) << std::setprecision(5) << x.bins_[b] << "\t";
		ofs << std::setw(8) << std::setprecision(5);
			print_free_energy(ofs, x.f_x_wham_[b], x.p_x_wham_[b]);
		ofs << std::setw(8) << std::setprecision(5) << x.global_sample_counts_[b];
		//<< std::setw(8) << std::setprecision(5) << wham_results_.error_f[b]
		ofs << "\n";
	}
	ofs.close(); ofs.clear();

	// "Rebiased" free energy distributions
	std::stringstream table_header_stream;
	table_header_stream << "# Data sets (by column)\n";
	for ( int i=0; i<num_simulations; ++i ) {
		table_header_stream << "# " << i+2 << ": " << simulations_[i].data_set_label << "\n";
	}
	table_header_stream << "#\n"
	                    << "# " << x.name_ << " | F(" << x.name_ << ") [kBT]\n";

	file_name = "F_" + x.name_ + "_rebiased.out";
	
	ofs.open( file_name );
	ofs << "# \"Rebiased\" free energy distributions: "
      << " F_{rebias,i}(" << x.name_ << ") [k_B*T]\n";
	ofs << header_stream.str() << table_header_stream.str();
	for ( int b=0; b<num_bins_x; ++b ) {
		ofs << x.bins_[b];
		for ( int k=0; k<num_simulations; ++k ) {
			ofs << "\t" << std::setw(8) << std::setprecision(5);
			print_free_energy(ofs, x.f_rebiased_[k][b], x.p_rebiased_[k][b]);
		}
		ofs << "\n";
	}
	ofs.close();

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


void Wham::print_f_x_y(
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

	// Print x-grid
	file_name = x.name_ + "_grid.out";
	ofs.open(file_name);
	sep = " ";
	for ( int i=0; i<num_bins_x; ++i ) {
		for ( int j=0; j<num_bins_y; ++j ) {
			ofs << bins_x[i] << sep;
		}
		ofs << "\n";
	}
	ofs.close(); ofs.clear();

	// Print y-grid
	file_name = y.name_ + "_grid.out";
	ofs.open(file_name);
	sep = " ";
	for ( int i=0; i<num_bins_x; ++i ) {
		for ( int j=0; j<num_bins_y; ++j ) {
			ofs << bins_y[j] << sep;
		}
		ofs << "\n";
	}
	ofs.close(); ofs.clear();

	// Print F(x,y)
	file_name = "F_" + x.name_ + "_" + y.name_ + "_WHAM.out";
	ofs.open(file_name);
	sep = " ";
	for ( int i=0; i<num_bins_x; ++i ) {
		for ( int j=0; j<num_bins_y; ++j ) {
			print_free_energy(ofs, f_x_y_wham[i][j], p_x_y_wham[i][j]);
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
