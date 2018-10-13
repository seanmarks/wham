#include "Wham.h"

Wham::Wham(
	const std::string& options_file, const std::string& data_summary_file,
	const std::string& biasing_parameters_file, const std::string& time_series_y_files_list
):
	options_file_(options_file),
	data_summary_file_(data_summary_file),
	biasing_parameters_file_(biasing_parameters_file),
	time_series_y_files_list_(time_series_y_files_list)
{
	// Read settings
	parseOptionsFile(options_file);
	readDataSummary(data_summary_file);

	// Read time series
	readRawDataFiles();
	readAuxiliaryTimeSeries(time_series_y_files_list);

	// Parse the biases present, and evaluate each sample using each
	// of the biases
	createBiases(biasing_parameters_file);
	evaluateBiases();

	//----- Precompute some constants (TODO: move to function -----//

	// Total number of samples
	num_samples_total_ = 0;
	int num_simulations = u_bias_as_other_.size();
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


	// Make histograms based on the raw data
	// - TODO Move to solve() or output functions ?
	constructRawDataHistograms();

	// FIXME DEBUG
	std::cout << "DONE SETUP\n";
	std::cout << "u_bias_as_other_.size() = " << u_bias_as_other_.size() 
	          << " (should be " << num_simulations << ")\n";
	for ( int j=0; j<num_simulations; ++j ) {
		int num_samples = u_bias_as_other_[j].size();
		std::cout << "  u_bias_as_other_[" << j << "].size() = " << num_samples
		          << " (should be " << simulations_[j].x_samples.size() << ")\n";
		for ( int i=0; i<num_samples; ++i ) {
			if ( static_cast<int>( u_bias_as_other_[j][i].size() ) != num_simulations ) {
				throw std::runtime_error("SIZE MISMATCH");
			}
		}
	}
}


void Wham::parseOptionsFile(const std::string& options_file)
{
	// Defaults
	wham_options_.round_t = true;
	wham_options_.x_name  = "x";
	wham_options_.y_name  = "y";
	wham_options_.col_x   = 1;
	wham_options_.col_y   = 1;

	// Working variables
	using BinStyle = Bins::BinStyle;
	BinStyle bin_style_x = BinStyle::Left;
	BinStyle bin_style_y = BinStyle::Left;
	double x_min, x_max, y_min, y_max;
	int    num_bins_x, num_bins_y;

	std::ifstream ifs;
	try {
		ifs.open(options_file);
	}
	catch (std::ios_base::failure failed_to_open) {
		std::stringstream err_ss;
		err_ss << "Wham::parseOptionsFile: Failed to open options file (exception: "
		       << failed_to_open.what() << ")." << "\n";
		throw std::runtime_error( err_ss.str() );
	}

	std::string line, token;
	while ( getline(ifs, line) ) {
		std::stringstream lineStream(line);

		lineStream >> token;

		if ( token[0] == '#' || line.empty() ) {
			// Comment or empty line
			continue;
		}
		else if ( token == "Temperature" ) { 
			lineStream >> wham_options_.T;
			wham_options_.kBT = K_B_ * wham_options_.T;
		}
		else if ( token == "OrderParameter_x" ) { lineStream >> wham_options_.x_name;  }
		else if ( token == "OrderParameter_y" ) { lineStream >> wham_options_.y_name;  }
		else if ( token == "DataColumn_x") { 
			// Input is indexed from 1
			lineStream >> wham_options_.col_x;
			--(wham_options_.col_x);
		}
		else if ( token == "DataColumn_y") { 
			// Input is indexed from 1
			lineStream >> wham_options_.col_y;
			--(wham_options_.col_y);
		}
		// Histogram for x
		else if ( token == "NumBins_x" )      { lineStream >> num_bins_x; }
		else if ( token == "Min_x" )          { lineStream >> x_min;   }
		else if ( token == "Max_x" )          { lineStream >> x_max;   }
		else if ( token == "BinStyle_x" ) { 
			lineStream >> token;
			bin_style_x = parseBinStyle(token);
		}
		// Histogram for y
		else if ( token == "NumBins_y" )      { lineStream >> num_bins_y; }
		else if ( token == "Min_y" )          { lineStream >> y_min;   }
		else if ( token == "Max_y" )          { lineStream >> y_max;   }
		else if ( token == "BinStyle_y" ) { 
			lineStream >> token;
			bin_style_y = parseBinStyle(token);
		}
		else {
			std::stringstream err_ss;
			err_ss << "Wham::parseOptionsFile: Error: unrecognized token \"" << token << "\".\n";
			throw std::runtime_error( err_ss.str() );
		}
	} // while getline

	ifs.close();

	// Set up bins
	bins_x_.set_bins(x_min, x_max, num_bins_x, bin_style_x);
	bins_y_.set_bins(y_min, y_max, num_bins_y, bin_style_y);
}


void Wham::readDataSummary(const std::string& data_summary_file)
{
	simulations_.clear();

	// Get the location of the data summary file
	std::string data_summary_path = get_dir(data_summary_file);

	std::ifstream ifs(data_summary_file);
	if ( not ifs.is_open() ) {
		throw std::runtime_error("Unable to open data summary file: " + data_summary_file + "\n");
	}

	std::string line, token, rel_path;

	while ( getline(ifs, line) ) {
		std::stringstream ss(line);
		ss >> token;

		if ( line.empty() or token[0] == '#' ) {
			continue;
		}

		// Create a new simulation
		simulations_.push_back( Simulation() );
		simulations_.back().data_set_label = token;

		// Samples file (path relative to data summary location)
		ss >> rel_path;  // samples file
		simulations_.back().data_file = get_realpath(data_summary_path + "/" + rel_path);

		// xtc file (ignore)
		ss >> token;

		// Production phase
		ss >> simulations_.back().t_min;
		ss >> simulations_.back().t_max;

		// TODO read from input for simulations at different temperatures?
		simulations_.back().kBT  = wham_options_.kBT;
		simulations_.back().beta = 1.0/( simulations_.back().kBT );
	}
}


void Wham::readRawDataFiles()
{
	// Working variables
	std::string line, token;
	double x, time, t_min, t_max;

	int num_simulations = simulations_.size();
	for ( int i=0; i<num_simulations; ++i ) {
		std::ifstream ifs( simulations_[i].data_file );
		if ( not ifs.is_open() ) {
			throw std::runtime_error("Failed to open data file \'" + simulations_[i].data_file + "\'");
		}

		t_min  = simulations_[i].t_min;
		t_max  = simulations_[i].t_max;

		while ( getline(ifs, line) ) {
			std::stringstream ss(line);
			ss >> token;

			// Ignore comments and blank lines
			if ( line.empty() or token[0] == '#' ) {
				continue;
			}

			// Get the time
			time = std::stod(token);
			if ( wham_options_.round_t ) { time = round(time); }

			// Check for production phase
			if ( time >= t_min and time <= t_max) {
				// Parse line
				std::vector<double> data(1, time);
				while ( ss >> token ) {
					data.push_back( std::stod(token) );
				}
				if ( static_cast<int>(data.size()) < wham_options_.col_x + 1 ) {
					throw std::runtime_error("col_x exceeds the number of columns in a data file");
				}

				// Store data
				x = data[ wham_options_.col_x ];
				simulations_[i].times.push_back( time );
				simulations_[i].x_samples.push_back( x );
				simulations_[i].nv_samples.push_back( -1 ); // FIXME
			}
		}
		ifs.close();

		// User feedback
		std::cout << " Data Set: " << simulations_[i].data_set_label << "\n";
		std::cout << "   File: " << simulations_[i].data_file << "\n";
		std::cout << "   Using " << simulations_[i].x_samples.size() << " data points.\n";
	}
}


void Wham::readAuxiliaryTimeSeries(const std::string& time_series_y_files_list)
{
	//----- Get names of files -----//

	std::ifstream list_ifs(time_series_y_files_list);
	if ( not list_ifs.is_open() ) {
		throw std::runtime_error("Failed to open file with y series \'" + time_series_y_files_list + "\'");
	}

	int i = 0;
	std::string line, token;
	while ( getline(list_ifs, line) ) {
		// Ignore possible comments and blank lines
		std::stringstream ss(line);
		ss >> token;
		if ( line.empty() or token[0] == '#' ) {
			// Comment or empty line
			continue;
		}

		ss >> simulations_[i].data_file_y;
		if ( time_series_y_files_list == data_summary_file_ ) {
			// Path is relative
			std::string data_summary_path = get_dir(data_summary_file_);
			simulations_[i].data_file_y = get_realpath(data_summary_path + "/" + simulations_[i].data_file_y);
		}

		++i;
	}
	list_ifs.close();
	// TODO check # simulations and data label


	//----- Read time series -----//

	int num_simulations = simulations_.size();
	for ( int i=0; i<num_simulations; ++i ) {
		std::ifstream ifs( simulations_[i].data_file_y );
		if ( not ifs.is_open() ) {
			throw std::runtime_error("Failed to open data file \'" + simulations_[i].data_file_y + "\'");
		}

		int index = 0;
		double time, y;
		double t_min = simulations_[i].t_min;
		double t_max = simulations_[i].t_max;

		while ( getline(ifs, line) ) {
			// Ignore possible comments and blank lines
			std::stringstream ss(line);
			ss >> token;
			if ( line.empty() or token[0] == '#' ) {
				// Comment or empty line
				continue;
			}

			// Get the time
			time = std::stod(token);
			if ( wham_options_.round_t ) { time = round(time); }

			if ( time >= t_min and time <= t_max) {
				// Production phase: check times vs. primary OP file
				if ( time != simulations_[i].times[index] ) {
					std::stringstream err_ss;
					err_ss << "Mismatch between primary OP file (" << simulations_[i].data_file << ")\n"
					       << "and secondary file (" << simulations_[i].data_file_y << ")\n"
					       << "Expected t=" << simulations_[i].times[index] << " ps, read " 
					          << time << " ps\n";
				}

				// Parse line
				std::vector<double> data(1, time);
				while ( ss >> token ) {
					data.push_back( std::stod(token) );
				}
				if ( static_cast<int>(data.size()) < wham_options_.col_y + 1 ) {
					throw std::runtime_error("Not enough samples for y in a row");
				}

				// Store data
				y = data[ wham_options_.col_y ];
				simulations_[i].y_samples.push_back( y );

				++index;
			}
		}
		ifs.close();
	}
}


void Wham::createBiases(const std::string& biasing_parameters_file)
{
	// Ensure that simulation data is available
	int num_simulations = simulations_.size();
	if ( num_simulations == 0 ) {
		throw std::runtime_error("No simulations are registered");
	}

	std::ifstream ifs(biasing_parameters_file);
	if ( not ifs.is_open() ) {
		throw std::runtime_error("Unable to open biasing parameters file: " + biasing_parameters_file);
	}

	std::string line, token;
	int index = 0;

	biases_.clear();
	std::vector<std::string> bias_tokens;

	while ( getline(ifs, line) ) {
		std::stringstream ss(line);
		ss >> token;

		if ( line.empty() or token[0] == '#' ) {
			continue;  // skip comments
		}

		// Check for internal consistency
		if ( index >= num_simulations ) {
			throw std::runtime_error("Too many simulations recorded in the biasing parameters file");
		}
		else if ( token != simulations_[index].data_set_label ) {
			std::stringstream err_ss;
			err_ss << "Expected to read data set with label \"" << simulations_[index].data_set_label 
			          << "\" in biasing parameters log,\n"
			       << "but encountered \"" << token << "\" instead\n";
			throw std::runtime_error( err_ss.str() ); 
		}

		// Create a new bias using the remaining tokens
		bias_tokens.resize(0);
		while ( ss >> token ) {
			bias_tokens.push_back( token );
		}
		biases_.push_back( Bias(bias_tokens, simulations_[index].kBT) );

		++index;
	}
}


void Wham::evaluateBiases()
{
	// TODO linearize
	// TODO need a better way to organize data for when multiple OPs are biased

	int num_simulations = simulations_.size();
	u_bias_as_other_.resize( num_simulations );

	std::vector<double> args(1);  // working variable

	for ( int j=0; j<num_simulations; ++j ) {
		int num_samples = simulations_[j].x_samples.size();
		u_bias_as_other_[j].resize(num_samples);

		for ( int i=0; i<num_samples; ++i ) {
			u_bias_as_other_[j][i].resize( num_simulations );

			// Evaluate this sample using each of the biases
			args = {{ simulations_[j].x_samples[i] }};
			for ( int r=0; r<num_simulations; ++r ) {
				u_bias_as_other_[j][i][r] = biases_[r].evaluate( args );
			}
		}
	}
}



// After reading input files, use to generate histograms of raw data from biased simulations
void Wham::constructRawDataHistograms()
{
	int num_simulations = static_cast<int>( simulations_.size() );
	int num_bins_x = bins_x_.get_num_bins();

	// Allocate memory
	for ( int i=0; i<num_simulations; ++i ) {
		simulations_[i].x_bins = bins_x_.get_bins();
		simulations_[i].sample_counts.resize(num_bins_x, 0);
		simulations_[i].p_biased.resize(num_bins_x);
		simulations_[i].f_biased.resize(num_bins_x);
		simulations_[i].f_unbiased.resize(num_bins_x);
	}

	std::vector<double> u_bias_tmp;
	for ( int i=0; i<num_simulations; ++i ) {
		// Assemble the biasing potential values for this simulation
		int num_samples = static_cast<int>( simulations_[i].x_samples.size() );
		u_bias_tmp.resize(num_samples);
		for ( int j=0; j<num_samples; ++j ) {
			u_bias_tmp[j] = u_bias_as_other_[i][j][i];
		}

		// Unbiased results, using only data from this simulation
		manually_unbias_f_x( 
				simulations_[i].x_samples, u_bias_tmp, simulations_[i].f_bias_guess, bins_x_, 
		    simulations_[i].p_unbiased, simulations_[i].f_unbiased, simulations_[i].sample_counts );

		// Derive remaining histograms
		double bin_size_x = bins_x_.get_bin_size();
		double num_samples_d = static_cast<double>( simulations_[i].x_samples.size() );
		for ( int b=0; b<num_bins_x; ++b ) {
			simulations_[i].p_biased[b] = (simulations_[i].sample_counts[b])/(bin_size_x*num_samples_d);
			simulations_[i].f_biased[b] = -log( simulations_[i].p_biased[b] );
		}

		// Statistics for this simulation
		simulations_[i].avg_x = average( simulations_[i].x_samples );
		simulations_[i].var_x = variance( simulations_[i].x_samples );
		simulations_[i].avg_nv = average( simulations_[i].nv_samples );
		simulations_[i].var_nv = variance( simulations_[i].nv_samples );
	}

	// Construct a histogram of the total number of data points in each bin, across all simulations
	sample_counts_x_.assign(num_bins_x, 0);
	for ( int i=0; i<num_simulations; ++i ) {
		for ( int b=0; b<num_bins_x; ++b ) {
			sample_counts_x_[b] += simulations_[i].sample_counts[b];
		}
	}
}


void Wham::solve()
{
	// One last sanity check
	int num_simulations = static_cast<int>( simulations_.size() );
	if ( num_simulations < 1 ) {
		throw std::runtime_error("Wham::solve(): Can't run WHAM without any simulations!\n");
	}

	// Initial guess
	// - Shift so that f[0] = 0
	std::vector<double> f_init(num_simulations, 0.0);
	double f_shift = simulations_[0].f_bias_guess;
	for ( int j=0; j<num_simulations; ++j ) {
		f_init[j] = simulations_[j].f_bias_guess - f_shift;
	}
	double tol = 1.0e-7; // TODO make input parameter

	// Compute differences btw. successive windows for initial guess
	int num_vars = num_simulations - 1;
	Wham::ColumnVector df(num_vars), df_init(num_vars);
	convert_f_to_df(f_init, df);
	df_init = df;  // save initial guess

	// Invoke functor wrappers for dlib calls
	WhamDlibEvalWrapper  evaluate_wrapper(*this);
	WhamDlibDerivWrapper derivatives_wrapper(*this);

	// Call dlib to perform the optimization
	double min_A = dlib::find_min( dlib::bfgs_search_strategy(), 
	                               dlib::objective_delta_stop_strategy(tol).be_verbose(),
	                               evaluate_wrapper,
	                               derivatives_wrapper,
	                               df,
	                               std::numeric_limits<double>::lowest() );

	// Compute optimal free energies of biasing from the differences between windows
	std::vector<double> f_opt(num_simulations);
	convert_df_to_f(df, f_opt);

	// Report results 
	std::cout << "min[A] = " << min_A << "\n";
	std::cout << "Biasing free energies:\n";
	for ( int i=0; i<num_simulations; ++i ) {
		std::cout << i << ":  " << f_opt[i] << "  (initial: " << f_init[i] << ")\n"; 
	}


	//----- Consensus histogram -----//

	int num_bins_x = bins_x_.get_num_bins();
	wham_results_.x_bins = bins_x_.get_bins();
	wham_results_.p_x_wham.assign(num_bins_x, 0.0);
	wham_results_.f_x_wham.assign(num_bins_x, 0.0);
	wham_results_.error_f.assign(num_bins_x, 0.0);
	wham_results_.sample_counts = sample_counts_x_;
	wham_results_.f_opt = f_opt;

	// FIXME wasteful to reformat data like this
	std::vector<std::vector<double>> x_samples(num_simulations);
	std::vector<std::vector<double>> y_samples(num_simulations);
	for ( int j=0; j<num_simulations; ++j ) {
		x_samples[j] = simulations_[j].x_samples;
		y_samples[j] = simulations_[j].y_samples;
	}

	compute_consensus_f_x( x_samples, u_bias_as_other_, f_opt, unbiased_ensemble_index_, bins_x_,
	                       wham_results_.p_x_wham, wham_results_.f_x_wham, wham_results_.sample_counts );


	//----- "Rebias" consensus histogram for validation -----//

	// TODO move to function
	std::vector<int> sample_counts_tmp;
	for ( int k=0; k<num_simulations; ++k ) {
		simulations_[k].f_rebiased.resize(num_bins_x);
		simulations_[k].p_rebiased.resize(num_bins_x);

		compute_consensus_f_x( x_samples, u_bias_as_other_, f_opt, k, bins_x_,
													 simulations_[k].p_rebiased, simulations_[k].f_rebiased,
		                       sample_counts_tmp );
	}


	//----- Diagnostics -----//

	// Shannon entropy (aka Kullback-Leibler divergence, relative entropy)
	// - See Hummer and Zhu (J Comp Chem 2012), Eqn. 2
	wham_results_.info_entropy.assign(num_simulations, 0.0);
	for ( int i=0; i<num_simulations; ++i ) {
		for ( int b=0; b<num_bins_x; ++b ) {
			const double& p_biased   = simulations_[i].p_biased[b];
			const double& f_rebiased = simulations_[i].f_rebiased[b];

			if ( p_biased > 0.0 and std::isfinite(f_rebiased) ) {
				wham_results_.info_entropy[i] -= p_biased*( simulations_[i].f_biased[b] - f_rebiased);
			}
		}
	}

	//----- Print WHAM results -----//

	printWhamResults();


	//----- Reweight to get F_0(y) -----//

	std::vector<double> p_y_wham, f_y_wham;
	std::vector<int> sample_counts_y;
	std::vector<std::vector<double>> p_x_y_wham, f_x_y_wham; 
	std::vector<std::vector<int>> sample_counts_x_y;

	compute_consensus_f_x_y( 
		x_samples, y_samples, u_bias_as_other_, f_opt, unbiased_ensemble_index_, bins_x_, bins_y_,
		// Output
		p_x_y_wham, f_x_y_wham, sample_counts_x_y,
		p_y_wham, f_y_wham, sample_counts_y );

	// Print results
	print_f_x_y_and_f_y(
		bins_x_, bins_y_, p_x_y_wham, f_x_y_wham, sample_counts_x_y,
		p_y_wham, f_y_wham, sample_counts_y );
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
const Wham::ColumnVector Wham::evalObjectiveDerivatives(const Wham::ColumnVector& df) const
{
	// Compute free energies of biasing from differences between windows
	int num_simulations = static_cast<int>( simulations_.size() );
	std::vector<double> f(num_simulations);
	convert_df_to_f(df, f);

	// Note:
	// - grad  = gradient wrt. free energy differences between neighboring windows
	// - dA_df = gradient wrt. biasing free energies themselves

	// Recompute log(sigma) as necessary
	if ( f != f_bias_last_ ) {
		compute_log_sigma( f, u_bias_as_other_, unbiased_ensemble_index_, log_sigma_0_ );
	}

	// First, compute derivatives of the objective fxn (A) wrt. the biasing free
	// energies themselves (f)
	// FIXME expensive nested loops
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
	const std::vector<double>& x, const std::vector<double>& u_bias, const double f,
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
	const std::vector<std::vector<double>>& x, 
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
	const std::vector<std::vector<double>>& x, const std::vector<std::vector<double>>& y,
	const std::vector<std::vector<std::vector<double>>>& u_bias_as_other,
	const std::vector<double>& f_opt, const int k, const Bins& bins_x, const Bins& bins_y,
	// Consensus distributions for F_k(x,y)
	std::vector<std::vector<double>>& p_x_y_wham, std::vector<std::vector<double>>& f_x_y_wham,
	std::vector<std::vector<int>>& sample_counts_x_y,
	// Reweighted results
	std::vector<double>& p_y_wham, std::vector<double>& f_y_wham, std::vector<int>& sample_counts_y
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
		sample_counts_x_y[a].resize( num_bins_y );
	}
	p_y_wham.resize( num_bins_y );
	f_y_wham.resize( num_bins_y );
	sample_counts_y.assign( num_bins_y, 0 );

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

		// Compute consensus distribution for this slice (note: here, 'x' is actually y)
		compute_consensus_f_x( y_given_x_a, u_bias_as_other_given_x_a, f_opt, k, bins_y,
		                       p_x_y_wham[a], f_x_y_wham[a], sample_counts_x_y[a] );
		for ( int b=0; b<num_bins_y; ++b ) {
			// Include normalization for x (the above function only normalizes for bin_size_y)
			p_x_y_wham[a][b] /= bin_size_x;
			f_x_y_wham[a][b] += log_bin_size_x;
		}

		// Accumulate sample count for y
		for ( int b=0; b<num_bins_y; ++b ) {
			sample_counts_y[b] += sample_counts_x_y[a][b];
		}
	}

	// Integrate out x to get F_k(y)
	std::vector<double> args(num_bins_x);
	for ( int b=0; b<num_bins_y; ++b ) {
		args.resize(0);
		for ( int a=0; a<num_bins_x; ++a ) {
			if ( p_x_y_wham[a][b] > 0.0 ) {
				args.push_back( -f_x_y_wham[a][b] + log_bin_size_x );
			}
		}

		if ( args.size() > 0 ) {
			f_y_wham[b] = -log_sum_exp( args );
			p_y_wham[b] = exp( -f_y_wham[b] );
		}
		else {
			p_y_wham[b] = 0.0;
			f_y_wham[b] = -log( p_y_wham[b] );
		}
	}
}


Wham::BinStyle Wham::parseBinStyle(const std::string& bin_style_token) const
{
	std::string lowercase_token = bin_style_token;
	std::transform( lowercase_token.begin(), lowercase_token.end(), 
	                lowercase_token.begin(), ::tolower );

	if ( lowercase_token == "left" ) {
		return BinStyle::Left;
	}
	else if ( lowercase_token == "center" ) {
		return BinStyle::Center;
	}
	else if ( lowercase_token == "right" ) {
		return BinStyle::Right;
	}
	else {
		throw std::runtime_error("bin style \"" + bin_style_token + "\" was not recognized");
	}
}


void Wham::printRawDistributions() const
{
	int num_simulations = static_cast<int>( simulations_.size() );
	if ( num_simulations < 1 ) {
		throw std::runtime_error("Wham::printRawDistributions: No simulations recorded.\n");
	}

	// Common header
	std::stringstream header_stream;

	std::stringstream table_header_stream;
	table_header_stream << "# " << wham_options_.x_name << " | " << wham_options_.x_name << "_star=";
	for ( int i=0; i<num_simulations; ++i ) {
		table_header_stream << "\t";
		table_header_stream << std::setw(8) << std::setprecision(5) << simulations_[i].x_star;  // FIXME xstar
	}
	table_header_stream << "\n";

	// Working variables
	std::string file_name;
	std::ofstream ofs;
	int num_bins_x = bins_x_.get_num_bins();  // All simulations are binned the same way

	//----- Print biased free energy distributions -----//

	file_name = "F_" + wham_options_.x_name + "_biased.out";
	ofs.open( file_name );
	ofs << "# Biased free energy distributions: "
      << " F_i(" << wham_options_.x_name << ") [k_B*T]\n";
	ofs << header_stream.str() << table_header_stream.str();
	for ( int b=0; b<num_bins_x; ++b ) {
		ofs << bins_x_[b];
		for ( int i=0; i<num_simulations; ++i ) {
			ofs << "\t";
			ofs << std::setw(8) << std::setprecision(5) << simulations_[i].f_biased[b];
		}
		ofs << "\n";
	}
	ofs.close();

	//----- Print unbiased free energy distributions (non-consensus) -----//

	file_name = "F_" + wham_options_.x_name + "_unbiased.out";
	ofs.open( file_name );
	ofs << "# Unbiased free energy distributions: "
      << " F_{0,i}(" << wham_options_.x_name << ") [k_B*T]\n";
	ofs << header_stream.str() << table_header_stream.str();
	for ( int b=0; b<num_bins_x; ++b ) {
		ofs << bins_x_[b];
		for ( int k=0; k<num_simulations; ++k ) {
			ofs << "\t" << std::setw(8) << std::setprecision(5);
			print_free_energy(ofs, simulations_[k].f_unbiased[b], simulations_[k].p_unbiased[b]);
		}
		ofs << "\n";
	}
	ofs.close();

	//----- Misc. Results -----//

	file_name = "info_for_plots.out";
	ofs.open( file_name );
	ofs << "# Useful information about each window for plotting\n";
	ofs << header_stream.str();
	ofs << "# " << wham_options_.x_name << "_star" << "\tavg\tvar\tavg(Nv)\tvar(Nv)\n";
	for ( int i=0; i<num_simulations; ++i ) {
		ofs << std::setw(8) << std::setprecision(5) << simulations_[i].x_star << "\t" // FIXME xstar
		    << std::setw(8) << std::setprecision(5) << simulations_[i].avg_x  << "\t"
		    << std::setw(8) << std::setprecision(5) << simulations_[i].var_x  << "\t"
		    << std::setw(8) << std::setprecision(5) << simulations_[i].avg_nv << "\t"
		    << std::setw(8) << std::setprecision(5) << simulations_[i].var_nv << "\n";
	}
	ofs.close();
}


void Wham::printWhamResults() const
{
	// Working variables
	std::stringstream header_stream;
	std::string file_name;
	std::ofstream ofs;
	const int num_simulations = simulations_.size();
	const int num_bins_x      = bins_x_.get_num_bins();

	// Print F_0(x)
	file_name = "F_" + wham_options_.x_name + "_WHAM.out";
	ofs.open(file_name);
	ofs << "# Consensus free energy distributions from WHAM: \n"
      << "#   F(" << wham_options_.x_name << ") [in k_B*T] with T = " << wham_options_.T << " K\n";
	ofs << "# " << wham_options_.x_name << "\tF[kBT]  NumSamples\n";  //"\t" << "\tstderr(F)\n"; TODO error estimate
	for ( int b=0; b<num_bins_x; ++b ) {
		ofs << std::setw(8) << std::setprecision(5) << wham_results_.x_bins[b] << "\t";
		ofs << std::setw(8) << std::setprecision(5);
			print_free_energy(ofs, wham_results_.f_x_wham[b], wham_results_.p_x_wham[b]);
		ofs << std::setw(8) << std::setprecision(5) << wham_results_.sample_counts[b];
		//<< std::setw(8) << std::setprecision(5) << wham_results_.error_f[b]
		ofs << "\n";
	}
	ofs.close(); ofs.clear();

	// Print optimal biasing free energies to file
	file_name = "f_bias_WHAM.out";
	ofs.open(file_name);
	ofs << "# F_bias [k_B*T] for each window after minimization\n";
	for ( int i=0; i<num_simulations; ++i ) {
		ofs << wham_results_.f_opt[i] << "\n";
	}
	ofs.close(); ofs.clear();

	// "Rebiased" free energy distributions
	std::stringstream table_header_stream;
	table_header_stream << "# " << wham_options_.x_name << " | " << wham_options_.x_name << "_star=";
	for ( int i=0; i<num_simulations; ++i ) {
		table_header_stream << "\t";
		table_header_stream << std::setw(8) << std::setprecision(5) << simulations_[i].x_star;  // FIXME xstar
	}
	table_header_stream << "\n";

	file_name = "F_" + wham_options_.x_name + "_rebiased.out";
	
	ofs.open( file_name );
	ofs << "# \"Rebiased\" free energy distributions: "
      << " F_{rebias,i}(" << wham_options_.x_name << ") [k_B*T]\n";
	ofs << header_stream.str() << table_header_stream.str();
	for ( int b=0; b<num_bins_x; ++b ) {
		ofs << bins_x_[b];
		for ( int k=0; k<num_simulations; ++k ) {
			ofs << "\t" << std::setw(8) << std::setprecision(5);
			print_free_energy(ofs, simulations_[k].f_rebiased[b], simulations_[k].p_rebiased[b]);
		}
		ofs << "\n";
	}
	ofs.close();

	// Misc. stats
	// FIXME xstar
	file_name = "stats_" + wham_options_.x_name + ".out";
	ofs.open(file_name);
	ofs << "# data_set   avg(x)   var(x)   avg(N)   var(N)   info_entropy\n";
	for ( int i=0; i<num_simulations; ++i ) {
		ofs << simulations_[i].data_set_label << "\t"
		    << simulations_[i].avg_x << "\t"
		    << simulations_[i].var_x << "\t"
		    << simulations_[i].avg_nv << "\t"
		    << simulations_[i].var_nv << "\t"
		    << wham_results_.info_entropy[i] << "\n";
	}
	ofs.close(); ofs.clear();

}


void Wham::print_f_x_y_and_f_y(
	// Consensus distributions for F_k(x,y)
	const Bins& bins_x, const Bins& bins_y,
	const std::vector<std::vector<double>>& p_x_y_wham,
	const std::vector<std::vector<double>>& f_x_y_wham,
	const std::vector<std::vector<int>>&    sample_counts_x_y,
	// Reweighted results for F_k(y)
	const std::vector<double>& p_y_wham, const std::vector<double>& f_y_wham,
	const std::vector<int>& sample_counts_y
) const
{
	// Working variables
	int num_bins_x = bins_x.get_num_bins();
	int num_bins_y = bins_y.get_num_bins();
	std::string file_name, sep;
	std::ofstream ofs;

	// Print x-grid
	file_name = wham_options_.x_name + "_grid.out";
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
	file_name = wham_options_.y_name + "_grid.out";
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
	file_name = "F_" + wham_options_.x_name + "_" + wham_options_.y_name + "_WHAM.out";
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
	file_name = "samples_" + wham_options_.x_name + "_" + wham_options_.y_name + ".out";
	ofs.open(file_name);
	sep = " ";
	for ( int i=0; i<num_bins_x; ++i ) {
		for ( int j=0; j<num_bins_y; ++j ) {
			ofs << sample_counts_x_y[i][j] << sep;
		}
		ofs << "\n";
	}
	ofs.close(); ofs.clear();

	// Print F(y)
	file_name = "F_" + wham_options_.y_name + "_reweighted.out";
	ofs.open(file_name);
	sep = "   ";
	ofs << "# " << wham_options_.y_name << "   F_reweighted[kBT]  NumSamples\n";
	for ( int j=0; j<num_bins_y; ++j ) {
		ofs << bins_y[j] << sep;
		print_free_energy(ofs, f_y_wham[j], p_y_wham[j]);
		ofs << sep << sample_counts_y[j] << "\n";
	}
	ofs.close(); ofs.clear();
}
