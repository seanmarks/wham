#include "Wham.h"


Wham::Wham(
	const DataSummary& data_summary, const OrderParameterRegistry& op_registry,
	const std::vector<Simulation>& simulations, const std::vector<OrderParameter>& order_parameters,
	const std::vector<Bias>& biases, const double tol
	// TODO initial guess
):
	data_summary_(data_summary),
	op_registry_(op_registry),
	simulations_(simulations),
	order_parameters_(order_parameters),
	biases_(biases),
	tol_(tol)
{
	setup();
}


void Wham::setup()
{
	int num_simulations = simulations_.size();

	// Use the time series of the first OP registered to count the number of samples
	// per simulation (and the total)
	num_samples_per_simulation_.resize(num_simulations);
	num_samples_total_ = 0;
	for ( int j=0; j<num_simulations; ++j ) {
		num_samples_per_simulation_[j] =  simulations_[j].get_num_samples(); //ref_op.get_time_series(j).size();
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

	// For each sample, evaluate the resulting potential under each bias
	// - Populates 'u_bias_as_other_'
	evaluateBiases();
}


void Wham::evaluateBiases()
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
			op_indices[r][p] = op_registry_.get_index(op_names[p]);
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
				int p;
				for ( int s=0; s<num_ops_for_bias; ++s ) {
					p = op_indices[r][s];
					args[s] = order_parameters_[p].get_time_series(j)[i];
					// FIXME cludgy
				}

				u_bias_as_other_[r][sample_index] = biases_[r].evaluate( args );
				++sample_index;
			}
		}
	}
}


// TODO add toggle for be_verbose
void Wham::solveWhamEquations(const std::vector<double>& f_init, std::vector<double>& f_bias_opt)
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
			dlib::objective_delta_stop_strategy(tol_).be_verbose(),
			evaluate_wrapper,
			derivatives_wrapper,
			df,
			std::numeric_limits<double>::lowest() 
	);

	// Compute optimal free energies of biasing from the differences between windows
	f_bias_opt.resize(num_simulations);
	convert_df_to_f(df, f_bias_opt);

	// Report results 
	if ( be_verbose ) {
		std::cout << "min[A] = " << min_A << "\n";
		std::cout << "Biasing free energies:\n";
		for ( int i=0; i<num_simulations; ++i ) {
			std::cout << i << ":  " << f_bias_opt[i] << "  (initial: " << f_init[i] << ")\n"; 
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


double Wham::evalObjectiveFunction(const Wham::ColumnVector& df) const
{
	// Compute free energies of biasing from the differences between windows
	int num_simulations = simulations_.size();
	std::vector<double> f(num_simulations);
	convert_df_to_f(df, f);

	// Compute log(sigma)
	// - These values are used here, and are stored for use when computing
	//   the derivatives later
	//   - TODO: Move buffer check elsewhere?
	compute_log_sigma( u_bias_as_other_, f, u_bias_as_other_unbiased_, f_unbiased_,
	                   log_sigma_unbiased_ );
	f_bias_last_ = f;

	// Compute the objective function's value
	double val = 0.0;
	for ( int n=0; n<num_samples_total_; ++n ) {
		val += log_sigma_unbiased_[n];
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


const Wham::ColumnVector Wham::evalObjectiveDerivatives(const Wham::ColumnVector& df) const
{
	// Compute free energies of biasing from differences between windows
	int num_simulations = simulations_.size();
	std::vector<double> f(num_simulations);
	convert_df_to_f(df, f);

	// Recompute log(sigma) as necessary
	// TODO: move buffer check to 'compute_log_sigma' instead?
	//   - Probably not safe, since the function is used on different
	//     subsets of data in compute_consensus_f_* functions
	//   - Maybe use a wrapper: "compute_log_sigma_all_data" ?
	if ( f != f_bias_last_ ) {  
		compute_log_sigma( u_bias_as_other_, f, u_bias_as_other_unbiased_, f_unbiased_,
		                   log_sigma_unbiased_ );
		f_bias_last_ = f;
	}

	// First, compute derivatives of the objective fxn (A) wrt. the biasing free
	// energies themselves (f)
	Wham::ColumnVector dA_df(num_simulations);
	dA_df(0) = 0.0;
	args_buffer_.resize(num_samples_total_); 
	double log_sum_exp_args;
	for ( int k=1; k<num_simulations; ++k ) {
		double fac = log_c_[k] + f[k];
		for ( int n=0; n<num_samples_total_; ++n ) {
			args_buffer_[n] = fac - u_bias_as_other_[k][n] - log_sigma_unbiased_[n];
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


void Wham::compute_log_sigma(
	const std::vector<std::vector<double>>& u_bias_as_other, const std::vector<double>& f,
	const std::vector<double>& u_bias_as_other_k, const double f_k,
	std::vector<double>& log_sigma
) const
{
	// Working variables
	int num_biases = u_bias_as_other_.size(); 
	std::vector<double> args(num_biases), facs(num_biases);
	for ( int r=0; r<num_biases; ++r ) {
		facs[r] = log_c_[r] - (f_k - f[r]);
	}

	// Compute log(sigma_k(x_{j,i})) for each sample
	int num_samples_total = u_bias_as_other_k.size();  // TODO check vs. each 'u_bias_as_other[r]'?
	log_sigma.resize(num_samples_total);
	for ( int n=0; n<num_samples_total; ++n ) {
		for ( int r=0; r<num_biases; ++r ) {
			args[r] = facs[r] + (u_bias_as_other_k[n] - u_bias_as_other[r][n]);
		}
		log_sigma[n] = log_sum_exp( args );
	}
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


std::vector<Distribution> Wham::manuallyUnbiasDistributions(const OrderParameter& x) const
{
	int num_simulations = static_cast<int>( simulations_.size() );
	std::vector<Distribution> unbiased_distributions(num_simulations);

	std::vector<double> u_bias_tmp;
	for ( int i=0; i<num_simulations; ++i ) {
		const TimeSeries& samples = x.get_time_series(i);
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
		manually_unbias_f_x( samples, u_bias_tmp, 0.0, x.get_bins(),
				                 unbiased_distributions[i] );
	}

	return unbiased_distributions;
}


// TODO way to merge with compute_consensus_f_x?
void Wham::manually_unbias_f_x(
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


// TODO Pass as OrderParameter 'x' instead? Then 'bins_x' is accessed from 'x'
void Wham::compute_consensus_f_x(
	const std::vector<TimeSeries>& x, 
	const std::vector<std::vector<double>>& u_bias_as_other,	const std::vector<double>& f_bias_opt,
	const std::vector<double>& u_bias_as_other_k, const double f_k,
	const Bins& bins_x,
	// Output
	Distribution& wham_distribution_x
) const
{
	// Compute log(sigma_k) using the buffer, since the array can be quite large
	compute_log_sigma(u_bias_as_other, f_bias_opt, u_bias_as_other_k, f_k, log_sigma_k_);

	// Reserve some memory for binning samples
	int num_bins_x = bins_x.get_num_bins();
	minus_log_sigma_k_binned_.resize(num_bins_x);
	for ( int b=0; b<num_bins_x; ++b ) {
		minus_log_sigma_k_binned_[b].resize( 0 );
		minus_log_sigma_k_binned_[b].reserve( x[0].size() );
	}

	// Sort log(sigma_k)-values by bin
	int bin, sample_index, num_simulations = simulations_.size();
	double bin_size_x = bins_x.get_bin_size();
	for ( int j=0; j<num_simulations; ++j ) {
		int num_samples = x[j].size();
		for ( int i=0; i<num_samples; ++i ) {
			bin = bins_x.find_bin( x[j][i] );
			if ( bin >= 0 ) {
				sample_index = simulation_data_ranges_[j].first + i;  
				minus_log_sigma_k_binned_[bin].push_back( -log_sigma_k_[sample_index] );
			}
		}
	}

	// Unpack for readability
	std::vector<double>& p_x_wham      = wham_distribution_x.p_x;
	std::vector<double>& f_x_wham      = wham_distribution_x.f_x;
	std::vector<int>&    sample_counts = wham_distribution_x.sample_counts;

	wham_distribution_x.bins_x = bins_x;

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


// TODO Pass 'x' and 'y' as OrderParameters instead (?)
void Wham::compute_consensus_f_x_y(
	const std::vector<TimeSeries>& x, const std::vector<TimeSeries>& y,
	const std::vector<std::vector<double>>& u_bias_as_other, const std::vector<double>& f_bias_opt,
	const std::vector<double>& u_bias_as_other_k, const double f_k,
	const Bins& bins_x, const Bins& bins_y,
	// Consensus distributions for F_k(x,y)
	std::vector<std::vector<double>>& p_x_y_wham, std::vector<std::vector<double>>& f_x_y_wham,
	std::vector<std::vector<int>>& sample_counts_x_y
) const
{
	// Input checks
	if ( y.size() != x.size() or f_bias_opt.size() != x.size() or u_bias_as_other.size() != x.size() ) {
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

	// TODO check vs. 'u_bias_as_other[r]'?
	int num_samples_total = u_bias_as_other_k.size();

	// Find each sample's (linearly indexed) bin
	sample_bins_.resize(num_samples_total);
	int sample_index = 0, bin_x, bin_y;
	int num_simulations = simulations_.size();
	for ( int j=0; j<num_simulations; ++j ) {
		int num_samples = x[j].size();
		for ( int i=0; i<num_samples; ++i ) {
			bin_x = bins_x.find_bin(x[j][i]);
			bin_y = bins_y.find_bin(y[j][i]);
			if ( bin_x >= 0 and bin_y >= 0 ) {
				sample_bins_[sample_index] = bin_x*num_bins_y + bin_y;
			}
			else {
				sample_bins_[sample_index] = -1;
			}
			sample_index++;
		}
	}

	// Allocate memory for working buffers
	int num_biases = u_bias_as_other_.size();
	std::vector<std::vector<double>> u_bias_as_other_binned(num_biases);
	std::vector<double> u_bias_as_other_k_binned;

	// Normalization
	double bin_area = bins_x.get_bin_size() * bins_y.get_bin_size();
	double normalization_factor = log( num_samples_total_ * bin_area );

	// Compute F_k(x_a, y_b) for each bin (a,b)
	int bin_index = 0;
	for ( int a=0; a<num_bins_x; ++a ) {
		for ( int b=0; b<num_bins_y; ++b ) {
			// Reset buffers
			for ( int r=0; r<num_biases; ++r ) {
				u_bias_as_other_binned[r].resize(0);
			}
			u_bias_as_other_k_binned.resize(0);

			// Find samples for this bin
			for ( int n=0; n<num_samples_total; ++n ) {
				if ( sample_bins_[n] == bin_index ) {
					for ( int r=0; r<num_biases; ++r ) {
						u_bias_as_other_binned[r].push_back( u_bias_as_other[r][n] );
					}
					u_bias_as_other_k_binned.push_back( u_bias_as_other_k[n] );
				}
			}

			// Compute weights
			compute_log_sigma( u_bias_as_other_binned, f_bias_opt, u_bias_as_other_k_binned, f_k,
			                   log_sigma_k_ );

			// Compute free energy and probability
			int num_samples_in_bin = log_sigma_k_.size();
			if ( num_samples_in_bin > 0 ) {
				// Perform the required log-sum-exp
				minus_log_sigma_k_.resize(num_samples_in_bin);
				for ( int s=0; s<num_samples_in_bin; ++s ) {
					minus_log_sigma_k_[s] = -log_sigma_k_[s];
				}
				double log_sum_exp_bin = log_sum_exp( minus_log_sigma_k_ );

				f_x_y_wham[a][b] = normalization_factor - log_sum_exp_bin;
				p_x_y_wham[a][b] = exp( -f_x_y_wham[a][b] );
			}
			else {
				f_x_y_wham[a][b] = 0.0;
				p_x_y_wham[a][b] = 0.0;
			}
			sample_counts_x_y[a][b] = num_samples_in_bin;

			// Go to next bin
			bin_index++;
		}
	}
}
