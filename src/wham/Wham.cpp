// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "Wham.h"
#include "Estimator_F_x.hpp"


Wham::Wham(
	const DataSummary& data_summary, const OrderParameterRegistry& op_registry,
	const std::vector<Simulation>& simulations, const std::vector<OrderParameter>& order_parameters,
	const std::vector<Bias>& biases, const std::vector<double>& f_bias_guess, const double tol
	// TODO initial guess
):
	data_summary_(data_summary),
	op_registry_(op_registry),
	simulations_(simulations),
	order_parameters_(order_parameters),
	biases_(biases),
	tol_(tol),
	f_bias_guess_(f_bias_guess),
	f_bias_opt_(f_bias_guess_)
{
	setup();

	// Solve for the unknown parameters: the free energies of turning on each bias
	f_bias_opt_ = solveWhamEquations(f_bias_guess_);
}


void Wham::setup()
{
	setup_timer_.start();
	int num_simulations = simulations_.size();

	// Use the time series of the first OP registered to count the number of samples
	// per simulation (and the total)
	num_samples_per_simulation_.resize(num_simulations);
	num_samples_total_ = 0;
	for ( int j=0; j<num_simulations; ++j ) {
		num_samples_per_simulation_[j] =  simulations_[j].getNumSamples();
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
	setup_timer_.stop();
}


void Wham::evaluateBiases()
{
	biases_timer_.start();

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
			op_indices[r][p] = op_registry_.nameToIndex(op_names[p]);
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
	#pragma omp parallel for
	for ( int r=0; r<num_biases; ++r ) {
		// Number of OPs involved in this bias
		int num_ops_for_bias = op_indices[r].size();
		std::vector<double> args(num_ops_for_bias);

		// Evaluate bias for each sample
		int sample_index = 0;
		for ( int j=0; j<num_simulations; ++j ) {
			int num_samples = num_samples_per_simulation_[j];
			for ( int i=0; i<num_samples; ++i ) {
				// Pull together the order parameter values needed for this bias
				int p;
				for ( int s=0; s<num_ops_for_bias; ++s ) {
					p = op_indices[r][s];
					args[s] = order_parameters_[p].getTimeSeries(j)[i];
					// FIXME cludgy
				}

				u_bias_as_other_[r][sample_index] = biases_[r].evaluate( args );
				++sample_index;
			}
		}
	}
	biases_timer_.stop();
}


// TODO: add toggle for be_verbose
std::vector<double> Wham::solveWhamEquations(const std::vector<double>& f_bias_guess)
{
	bool be_verbose = false;

	// Compute differences btw. successive windows for initial guess
	int num_simulations = simulations_.size();
	int num_vars = num_simulations - 1;
	Wham::ColumnVector df(num_vars), df_bias_guess(num_vars);
	convert_f_to_df(f_bias_guess, df);
	df_bias_guess = df;  // save initial guess

	// Invoke functor wrappers for dlib calls
	WhamDlibEvalWrapper  evaluate_wrapper(*this);
	WhamDlibDerivWrapper derivatives_wrapper(*this);

	// Call dlib to perform the optimization
	solve_timer_.start();
	double min_A = dlib::find_min( 
			dlib::bfgs_search_strategy(), 
			dlib::objective_delta_stop_strategy(tol_),
			//dlib::objective_delta_stop_strategy(tol_).be_verbose(),
			evaluate_wrapper,
			derivatives_wrapper,
			df,
			std::numeric_limits<double>::lowest() 
	);
	solve_timer_.stop();

	// Compute optimal free energies of biasing from the differences between windows
	convert_df_to_f(df, f_bias_opt_);

	// Current, the free energy of the first biased ensemble is zero.
	// Set the free energy for the *unbiased* ensemble as the zero point instead.
	compute_log_dhat( u_bias_as_other_, f_bias_opt_, log_dhat_ );
	/*
	FIXME: worth it?
	double delta_f = compute_consensus_f_k( u_bias_as_other_unbiased_ ); 
	for ( int i=0; i<num_simulations; ++i ) {
		f_bias_opt_[i] -= delta_f;
	}
	compute_log_dhat( u_bias_as_other_, f_bias_opt_, log_dhat_ );  // update
	*/

	if ( be_verbose ) {
		// Report results 
		std::cout << "min[A] = " << min_A << "\n";
		std::cout << "Biasing free energies:\n";
		for ( int i=0; i<num_simulations; ++i ) {
			std::cout << i << ":  " << f_bias_opt_[i] << "  (initial: " << f_bias_guess[i] << ")\n"; 
		}
	}
	
	// Precompute weights corresponding to simulation data
	w_.set_size( num_samples_total_, num_simulations );
	#pragma omp parallel for
	for ( int n=0; n<num_samples_total_; ++n ) {
		for ( int i=0; i<num_simulations; ++i ) {
			double arg = f_bias_opt_[i] - u_bias_as_other_[i][n] - log_dhat_[n];
			if ( arg > MIN_DBL_FOR_EXP ) {
				w_(n,i) = exp( arg );
			}
			else {
				w_(n,i) = 0.0;  // vanishingly small
			}
		}
	}

	/*
	// FIXME: Attempted to implement covariance matrix estimation, but the results
	// always seem to come out garbage (even checked vs. NumPy)
	wT_w_ = dlib::trans(w_) * w_;  // commonly used submatrix

	// Move to separate function(s)/only compute if error is desired?
	Matrix theta;
	compute_cov_matrix( w_, wT_w_, num_samples_per_simulation_, theta );
	error_f_bias_opt_.resize(num_simulations);
	for ( int i=0; i<num_simulations; ++i ) {
		error_f_bias_opt_[i] = sqrt( theta(i,i) );
	}
	// DEBUG
	for ( int i=0; i<num_simulations; ++i ) {
		std::cout << i << ":  " << f_bias_opt_[i] << " +/- " << error_f_bias_opt_[i] 
		          << "  (theta = cov = " << theta(i,i) << ")\n";
	}
	*/

	return f_bias_opt_;
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
	objective_timer_.start();

	// Compute free energies of biasing from the differences between windows
	int num_simulations = simulations_.size();
	std::vector<double> f(num_simulations);
	convert_df_to_f(df, f);

	// Compute log(dhat)
	// - These values are used here, and are stored for use when computing
	//   the derivatives later
	compute_log_dhat( u_bias_as_other_, f, log_dhat_tmp_ );
	f_bias_last_ = f;

	// Compute the objective function's value
	double val = 0.0;
	#pragma omp parallel for reduction(+:val)
	for ( int n=0; n<num_samples_total_; ++n ) {
		val += log_dhat_tmp_[n];
	}
	val *= inv_num_samples_total_;
	for ( int r=0; r<num_simulations; ++r ) {
		val -= c_[r]*f[r];
	}

#ifdef DEBUG_MODE
	std::cout << "EVALUATED: A = " << val << "\n";
#endif

	objective_timer_.stop();

	return val;
}


const Wham::ColumnVector Wham::evalObjectiveDerivatives(const Wham::ColumnVector& df) const
{
	gradient_timer_.start();

	// Compute free energies of biasing from differences between windows
	int num_simulations = simulations_.size();
	std::vector<double> f(num_simulations);
	convert_df_to_f(df, f);

	// Recompute log(dhat) as necessary
	// TODO: move buffer check to 'compute_log_dhat' instead?
	//   - Probably not safe, since the function is used on different
	//     subsets of data in compute_consensus_f_* functions
	//   - Maybe use a wrapper: "compute_log_dhat_all_data" ?
	if ( f != f_bias_last_ ) {  
		compute_log_dhat( u_bias_as_other_, f, log_dhat_tmp_ );
		f_bias_last_ = f;
	}

	int num_threads = OpenMP::get_max_threads();
	args_buffers_.resize(num_threads); 

	Wham::ColumnVector dA_df(num_simulations);
	dA_df(0) = 0.0;  // fix \hat{f}_1 = 0 (first ensemble's free energy)

	#pragma omp parallel
	{
		// Prepare buffer
		const int thread_id = OpenMP::get_thread_num();
		auto& args_buffer   = args_buffers_[thread_id];

		// First, compute derivatives of the objective fxn (A) wrt. the biased free
		// energies themselves (f)
		args_buffer.resize(num_samples_total_); 
		#pragma omp for
		for ( int k=1; k<num_simulations; ++k ) {
			gradient_omp_timer_.start();

			double fac = log(num_samples_per_simulation_[k]) + f[k];
			#pragma omp simd
			for ( int n=0; n<num_samples_total_; ++n ) {
				args_buffer[n] = fac - log_dhat_tmp_[n] - u_bias_as_other_[k][n];
			}
			double log_sum = log_sum_exp(args_buffer);

			dA_df(k) = inv_num_samples_total_*exp(log_sum) - c_[k];

			gradient_omp_timer_.stop();
		}
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

	gradient_timer_.stop();

	return grad;
}


void Wham::compute_log_sigma(
	const std::vector<std::vector<double>>& u_bias_as_other, const std::vector<double>& f,
	const std::vector<double>& u_bias_as_k, const double f_bias_k,
	std::vector<double>& log_sigma
) const
{
	log_sigma_timer_.start();

	int num_samples_total = u_bias_as_k.size();
	if ( num_samples_total == 0 ) {
		// Nothing to do
		log_sigma.resize(0);
		log_sigma_timer_.stop();
		return;
	}
//#ifdef WHAM_DEBUG
// TODO check consistency of # samples btw. u_bias_as_k and each 'u_bias_as_other[r]'
//#endif

	// Precompute some terms
	int num_biases = u_bias_as_other.size(); 
	std::vector<double> common_terms(num_biases);
	for ( int r=0; r<num_biases; ++r ) {
		common_terms[r] = log_c_[r] - (f_bias_k - f[r]);
	}

	// Compute log(sigma_k(x_{j,i})) for each sample
	log_sigma.resize(num_samples_total);
	#pragma omp parallel
	{
		std::vector<double> args(num_biases);

		#pragma omp for schedule(static,8)
		for ( int n=0; n<num_samples_total; ++n ) {
			log_sigma_omp_timer_.start();

			for ( int r=0; r<num_biases; ++r ) {
				args[r] = common_terms[r] + (u_bias_as_k[n] - u_bias_as_other[r][n]);
			}
			log_sigma[n] = log_sum_exp( args );

			log_sigma_omp_timer_.stop();
		}
	}

	log_sigma_timer_.stop();
}


void Wham::compute_log_dhat(
	const std::vector<std::vector<double>>& u_bias_as_other, const std::vector<double>& fhat,
	std::vector<double>& log_dhat
) const
{
	// TODO: input checks

	log_dhat_timer_.start();

	// Precompute some terms
	int num_biases = u_bias_as_other.size(); 
	std::vector<double> common_terms(num_biases);
	for ( int r=0; r<num_biases; ++r ) {
		common_terms[r] = log(num_samples_per_simulation_[r]) + fhat[r];
	}

	// Compute log(dhat(x_n)) for each sample 'n'
	int num_samples_total = u_bias_as_other[0].size();
	log_dhat.resize(num_samples_total);
	#pragma omp parallel
	{
		std::vector<double> args(num_biases);

		log_dhat_omp_timer_.start();
		#pragma omp for schedule(static,8)
		for ( int n=0; n<num_samples_total; ++n ) {
			for ( int r=0; r<num_biases; ++r ) {
				args[r] = common_terms[r] - u_bias_as_other[r][n];
			}
			log_dhat[n] = log_sum_exp( args );
		}
		log_dhat_omp_timer_.stop();
	}

	log_dhat_timer_.stop();
}


double Wham::log_sum_exp(const std::vector<double>& args) const
{
	log_sum_exp_timer_.start();

	int num_args = static_cast<int>( args.size() );
	FANCY_ASSERT( num_args > 0, "no arguments supplied" );

	// Compute sum of exp(args[i] - max_arg)
	double max_arg = *std::max_element( args.begin(), args.end() );
	double sum = 0.0;
	#pragma omp simd reduction(+: sum)
	for ( int i=0; i<num_args; ++i ) {
		sum += std::exp(args[i] - max_arg);
	}

	// Since exp(max_arg) might be huge, return log(sum_exp) instead of sum_exp
	double log_sum_exp_out = max_arg + std::log(sum);
	log_sum_exp_timer_.stop();

	return log_sum_exp_out;
}


double Wham::weighted_sum_exp(
	const std::vector<double>& args, const std::vector<double>& weights
) const
{
	weighted_sum_exp_timer_.start();

	// Check input
	int num_args = args.size();
	FANCY_ASSERT( num_args > 0, "no arguments supplied" );

	// Compute sum of exp(args[i] - max_arg)
	double max_arg = *std::max_element( args.begin(), args.end() );
	double sum = 0.0;
	#pragma omp simd reduction(+: sum)
	for ( int i=0; i<num_args; ++i ) {
		sum += weights[i]*exp(args[i] - max_arg);
	}

	// FIXME: How to handle huge values of x_max when sum < 0.0 and returning log(sum)
	// is invalid? Maybe return another variable indicating the sign?
	sum = exp(max_arg)*sum;

	weighted_sum_exp_timer_.stop();

	return sum;
}


std::vector<Distribution> Wham::manuallyUnbiasDistributions(const std::string& op_name) const
{
	int p = op_registry_.nameToIndex(op_name);
	const auto& x = order_parameters_[p];

	int num_simulations = static_cast<int>( simulations_.size() );
	std::vector<Distribution> unbiased_distributions(num_simulations);

	std::vector<double> u_bias_tmp;
	for ( int i=0; i<num_simulations; ++i ) {
		const TimeSeries& samples = x.getTimeSeries(i);
		int num_samples = samples.size();

		// Assemble the biasing potential values for this simulation's data under its own bias
		u_bias_tmp.resize(num_samples);
		int index = simulation_data_ranges_[i].first;
		for ( int j=0; j<num_samples; ++j ) {
			u_bias_tmp[j] = u_bias_as_other_[i][index];
			++index;
		}

		// Unbiased results, using only data from this simulation
		manually_unbias_f_x( samples, u_bias_tmp, 0.0, x.getBins(),
				                 unbiased_distributions[i] );
	}

	return unbiased_distributions;
}


// TODO: way to merge with compute_consensus_f_x?
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



double Wham::compute_consensus_f_k(const std::vector<double>& u_bias_as_k) const
{
	FANCY_ASSERT( static_cast<int>(u_bias_as_k.size()) == num_samples_total_, "size mismatch" );

	// TODO: optional OpenMP switch?
	std::vector<double> args(num_samples_total_);
	#pragma omp parallel for schedule(static, 8)
	for ( int n=0; n<num_samples_total_; ++n ) {
		args[n] = -u_bias_as_k[n] - log_dhat_[n];
	}

	double f_k = -log_sum_exp( args );
	return f_k;
}



Distribution Wham::compute_consensus_f_x_unbiased(const std::string& op_name) const
{
	int p = op_registry_.nameToIndex(op_name);

	Estimator_F_x est_f_x(order_parameters_[p]);
	est_f_x.calculate(*this, u_bias_as_other_unbiased_, f_unbiased_);
	return est_f_x.get_f_x();
}



Distribution Wham::compute_consensus_f_x_rebiased(const std::string& op_name, const std::string& data_set_label) const
{
	int p = op_registry_.nameToIndex(op_name);
	int j = data_summary_.dataSetLabelToIndex(data_set_label);

	Estimator_F_x est_f_x(order_parameters_[p]);
	est_f_x.calculate(*this, u_bias_as_other_[j], f_bias_opt_[j]);
	return est_f_x.get_f_x();
}



void Wham::compute_consensus_f_x_y_unbiased(
	const std::string& x_name, const std::string& y_name,
	std::vector<std::vector<double>>& p_x_y_wham, std::vector<std::vector<double>>& f_x_y_wham,
	std::vector<std::vector<int>>& sample_counts_x_y
) const
{
	int p = op_registry_.nameToIndex(x_name);
	int q = op_registry_.nameToIndex(y_name);

	compute_consensus_f_x_y( order_parameters_[p], order_parameters_[q], u_bias_as_other_, f_bias_opt_,
	                         u_bias_as_other_unbiased_, f_unbiased_, 
	                         p_x_y_wham, f_x_y_wham, sample_counts_x_y );
}


void Wham::compute_consensus_f_x_y(
	const OrderParameter& x, const OrderParameter& y,
	const std::vector<std::vector<double>>& u_bias_as_other, const std::vector<double>& f_bias_opt,
	const std::vector<double>& u_bias_as_k, const double f_bias_k,
	// Consensus distributions for F_k(x,y)
	std::vector<std::vector<double>>& p_x_y_wham, std::vector<std::vector<double>>& f_x_y_wham,
	std::vector<std::vector<int>>& sample_counts_x_y
) const
{
	f_x_y_timer_.start();
	/*
	// Input checks TODO
	if ( y.size() != x.size() or f_bias_opt.size() != x.size() or u_bias_as_other.size() != x.size() ) {
		throw std::runtime_error("Error in Wham::compute_consensus_f_x_y - Array size inconsistency");
	}
	*/

	// Compute weight factors
	compute_log_sigma(u_bias_as_other, f_bias_opt, u_bias_as_k, f_bias_k, log_sigma_k_);

	const auto& bins_x = x.getBins();
	const auto& bins_y = y.getBins();
	int num_bins_x = bins_x.get_num_bins();
	int num_bins_y = bins_y.get_num_bins();
	int num_bins_total = num_bins_x*num_bins_y;

	// Allocate memory and set up grids
	// - TODO: Wrap in a simple struct/class
	p_x_y_wham.resize( num_bins_x );
	f_x_y_wham.resize( num_bins_x );
	sample_counts_x_y.resize( num_bins_x );
	for ( int a=0; a<num_bins_x; ++a ) {
		// Second dimension
		p_x_y_wham[a].resize( num_bins_y );
		f_x_y_wham[a].resize( num_bins_y );
		sample_counts_x_y[a].resize( num_bins_y );
	}

	minus_log_sigma_k_binned_.resize(num_bins_total);
	for ( int b=0; b<num_bins_total; ++b ) {
		minus_log_sigma_k_binned_[b].resize( 0 );
		minus_log_sigma_k_binned_[b].reserve( x.getTimeSeries(0).size()/10 );
	}

	// Sort log(sigma_k)-values by bin
	//sample_bins_.resize(num_samples_total);
	// TODO: parallelize sorting?
	int sample_index, bin_x, bin_y;
	int num_simulations = simulations_.size();
	int bin_index;
	f_x_y_sort_timer_.start();
	for ( int j=0; j<num_simulations; ++j ) {
		const auto& x_j = x.getTimeSeries(j);
		const auto& y_j = y.getTimeSeries(j);
		int num_samples = x_j.size();
		for ( int i=0; i<num_samples; ++i ) {
			bin_x = bins_x.find_bin(x_j[i]);
			bin_y = bins_y.find_bin(y_j[i]);
			if ( bin_x >= 0 and bin_y >= 0 ) {
				// Sample is in the binned ranges
				bin_index    = bin_x*num_bins_y + bin_y;
				sample_index = simulation_data_ranges_[j].first + i;  
				minus_log_sigma_k_binned_[bin_index].push_back( -log_sigma_k_[sample_index] );
			}
		}
	}
	f_x_y_sort_timer_.stop();

	// Normalization
	const double bin_area = bins_x.get_bin_size() * bins_y.get_bin_size();
	const double normalization_factor = log( num_samples_total_ * bin_area );

	// Compute F_k(x_a, y_b) for each bin (a,b)
	#pragma omp parallel for collapse(2) //schedule(static,8)
	for ( int a=0; a<num_bins_x; ++a ) {
		for ( int b=0; b<num_bins_y; ++b ) {
			int bin_index = a*num_bins_y + b;

			int num_samples_in_bin = minus_log_sigma_k_binned_[bin_index].size();
			if ( num_samples_in_bin > 0 ) {
				double log_sum = log_sum_exp( minus_log_sigma_k_binned_[bin_index] );
				f_x_y_wham[a][b] = normalization_factor - log_sum;
				p_x_y_wham[a][b] = exp( -f_x_y_wham[a][b] );
			}
			
			else {
				f_x_y_wham[a][b] = 0.0;
				p_x_y_wham[a][b] = 0.0;
			}
			sample_counts_x_y[a][b] = num_samples_in_bin;
		}
	}

	f_x_y_timer_.stop();
}


void Wham::compute_error_avg_x_k(
	const std::vector<double>& x,
	const std::vector<double>& u_bias_as_k,
	// Output
	Matrix& w,                      // W_aug, augmented matrix of weights
	Matrix& wT_w,                   // W_aug^T * W_aug
	std::vector<int>& num_samples,  // per state (augmented)
	Matrix& theta
) const
{
	int num_simulations = simulations_.size();
	int num_states      = num_simulations + 2;

	int num_samples_total = x.size();
	FANCY_ASSERT( num_samples_total == num_samples_total_, "array length mismatch" );

	double f_k = compute_consensus_f_k(u_bias_as_k);

	// Compute weights for related to <x>
	ColumnVector w_bot(num_samples_total), w_top(num_samples_total);
	double sum = 0.0;
	#pragma omp parallel for schedule(static,8) reduction(+:sum)
	for ( int n=0; n<num_samples_total; ++n ) {
		w_bot(n) = exp(f_k - u_bias_as_k[n] - log_dhat_[n]);  // denominator
		w_top(n) = x[n]*w_bot(n);  // numerator
		sum += w_top(n);
	}
	w_top /= sum;

	// Note: better to copy unaugmented parts of matrices from member variables than to
	// recompute them each time

	// Augmented W
	// - Unaugmented part
	w.set_size(num_samples_total, num_states);
	dlib::set_subm( w, dlib::range(0, num_samples_total), dlib::range(0, num_simulations) ) = w_;
	// - Last two columns
	int i_top = num_simulations;
	int i_bot = i_top + 1;
	dlib::set_subm( w, dlib::range(0, num_samples_total), dlib::range(i_top, i_bot     ) ) = w_top;
	dlib::set_subm( w, dlib::range(0, num_samples_total), dlib::range(i_bot, num_states) ) = w_bot;

	// Augmented W^T*W
	//int n_top = num_samples_total;
	//int n_bot = n_top + 1;
	// - Unaugmented part
	wT_w.set_size(num_states, num_states);
	dlib::set_subm( wT_w, dlib::range(0, num_samples_total), dlib::range(0, num_simulations) ) = wT_w_;
	// - Outer edges
	RowVector r_top = dlib::trans(w_top) * w;
	RowVector r_bot = dlib::trans(w_bot) * w;
	dlib::set_subm( wT_w, dlib::range(i_top, i_bot),       dlib::range(0, num_simulations) ) = r_top;
	dlib::set_subm( wT_w, dlib::range(i_bot, num_states),  dlib::range(0, num_simulations) ) = r_bot;
	dlib::set_subm( wT_w, dlib::range(0, num_simulations), dlib::range(i_top, i_bot)           ) = dlib::trans(r_top);
	dlib::set_subm( wT_w, dlib::range(0, num_simulations), dlib::range(i_bot, num_states)    ) = dlib::trans(r_bot);

	// The fictitious states "top" and "bot" do not (formally) contribute any samples
	num_samples.assign(num_states, 0);
	for ( int i=0; i<num_simulations; ++i ) {
		num_samples[i] = num_samples_per_simulation_[i];
	}

	// Augmented covariance matrix
	compute_cov_matrix(w, wT_w, num_samples, theta);

	// TODO Propagate errors from Theta to errors in <x> using the Delta Method

	return;
}


void Wham::compute_cov_matrix(
	const Matrix& w, const Matrix& wT_w, const std::vector<int>& num_samples,
	Matrix& theta
) const
{
	int num_states = w.nc();
	// TODO size checks

	// Eigenvalue decomposition:
	//     (W^T*W) * V = V * D
	// - D = diag(lambdas)
	// - V = matrix of eigenvectors
	// - Note that W^T*W is symmetric
	//   - All eigenvalues are real and positive
	//   - Matrix of eigenvectors, V, is orthogonal
	//   - dlib::make_symmetric() returns a matrix_exp flagged as symmetric, but
	//     eigenvalue_decomposition also checks for this when given a general matrix
	using EigenvalueDecomposition = dlib::eigenvalue_decomposition<Matrix>;
	EigenvalueDecomposition eigen_decomp( wT_w );
	const Matrix&       V       = eigen_decomp.get_pseudo_v();
	const ColumnVector& lambdas = eigen_decomp.get_real_eigenvalues();

	// With the eigenvalue decomposition of W^T*W known, the necessary parts of the
	// singular value decomposition (SVD) of W are trivial
	// - W = U*Sigma*V^T, but U is not needed (a good thing, because it's large: N x N)
	// - From above: W^T * W = V*D*V^T
	Matrix sigma = dlib::zeros_matrix<double>(num_states, num_states);
	int num_lambdas = lambdas.size();
	for ( int i=0; i<num_lambdas; ++i ) {
		if ( lambdas(i) >= 0.0 ) {
			sigma(i,i) = sqrt( lambdas(i) );
		}
	}
	Matrix V_sigma = V*sigma;

	// Use the SVD of W to estimate Theta
	Matrix N = dlib::zeros_matrix<double>(num_states, num_states);
	for ( int i=0; i<num_states; ++i ) {
		N(i,i) = static_cast<double>( num_samples[i] );
	}
	Matrix I      = dlib::identity_matrix<double>(num_states);
	Matrix M      = I - dlib::trans(V_sigma)*N*V_sigma;
	Matrix M_pinv = dlib::pinv(M);  // M is square, but usually singular
	theta = V_sigma * M_pinv * dlib::trans(V_sigma);

	/*
	// DEBUG
	RowVector sum_row_w = dlib::sum_rows(w);  // sum of each column should be 1
	std::cout << "num_samples = " << w.nr() << ", num_states = " << w.nc() << std::endl;
	std::cout << "w(sum_row)=\n"
	          << "  shape = (" << sum_row_w.nr() << ", " << sum_row_w.nc() << ")\n"
	          << sum_row_w    << std::endl;

	std::cout << "V=\n"        << V            << std::endl;
	std::cout << "lambdas=\n"  << lambdas      << std::endl;
	std::cout << "sigma=\n"    << sigma        << std::endl;
	std::cout << "N=\n"        << N            << std::endl;
	std::cout << "I=\n"        << I            << std::endl;
	std::cout << "M=\n"        << M            << std::endl;
	std::cout << "det(M)="     << dlib::det(M) << std::endl;  // likely (nearly) singular
	std::cout << "pinv(M)=\n"  << M_pinv       << std::endl;
	std::cout << "theta=\n"    << theta        << std::endl;

	// DEBUG: Check major results
	Matrix V_trans = dlib::trans(V);
	Matrix D       = dlib::diagm(lambdas);
	double error_eig  = norm_frob( wT_w - V*D*V_trans );
	std::cout << "error(eig(wT_w)) = " << error_eig << std::endl;
	double error_pinv = norm_frob( M*M_pinv - I );
	std::cout << "error(pinv(M)) = " << error_pinv << std::endl;
	std::cout << std::endl;

	// DEBUG: Save key arrays
	std::ofstream ofs;
	ofs.open("DEBUG_w.out");
	ofs << "# matrix of weights, W\n" << w;
	ofs.close(); ofs.clear();
	ofs.open("DEBUG_N.out");
	ofs << "# Number of samples from each simulation, N_i\n" << dlib::diag(N);
	ofs.close(); ofs.clear();

	// DEBUG: full, nasty, expensive calculation using w directly
	Matrix i_N       = dlib::identity_matrix<double>(w.nr());
	Matrix m_big     = i_N - w*N*dlib::trans(w);
	Matrix cov_f_alt = dlib::trans(w) * dlib::pinv(m_big) * w;
	ofs.open("DEBUG_cov_f_alt.out");
	ofs << "# Covariance matrix, Theta, using W directly\n" << w;
	ofs.close(); ofs.clear();
	*/

	return;
}
