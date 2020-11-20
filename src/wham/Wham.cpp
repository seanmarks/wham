// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "Wham.h"

#include "Estimator_F_x.hpp"
#include "Estimator_F_x_y.hpp"
#include "LogSumExp.hpp"


Wham::Wham(
	const OrderParameterRegistry& op_registry,
	const std::vector<Simulation>& simulations, const std::vector<OrderParameter>& order_parameters,
	const std::vector<Bias>& biases, const std::vector<double>& f_bias_guess, const double tol
	// TODO: initial guess
):
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
	for ( int r=0; r<num_simulations; ++r ) {
		c_[r] = static_cast<double>(num_samples_per_simulation_[r]) * inv_num_samples_total_;
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
	compute_log_dhat( u_bias_as_other_, f_bias_opt_, log_dhat_opt_ );
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
	
	// TODO: Precompute weights for stored simulations in 'w_opt_'?

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
			double log_sum = logSumExp(args_buffer);

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
			log_dhat[n] = logSumExp( args );
		}
		log_dhat_omp_timer_.stop();
	}

	log_dhat_timer_.stop();
}


double Wham::logSumExp(const std::vector<double>& args) const
{
	log_sum_exp_timer_.start();
	double sum = numeric::logSumExp(args);
	log_sum_exp_timer_.stop();

	return sum;
}


double Wham::weightedSumExp(
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


Wham::Vector<Wham::Vector<double>> Wham::computeUnbiasedNonConsensusLogWeights() const
{
	const int num_simulations = simulations_.size();
	Vector<Vector<double>> log_w(num_simulations);
	for ( int i=0; i<num_simulations; ++i ) {
		const int num_samples = simulations_[i].getNumSamples();
		log_w[i].resize(num_samples);
	}

	#pragma omp parallel
	{
		for ( int i=0; i<num_simulations; ++i ) {			
			const int num_samples = simulations_[i].getNumSamples();
			const double log_N_i = std::log(num_samples);

			const auto& u_bias_ens_i = u_bias_as_other_[i];
			const int offset = simulation_data_ranges_[i].first;
			#pragma omp for simd schedule(simd:static)
			for ( int n=0; n<num_samples; ++n ) {
				log_w[i][n] = u_bias_ens_i[offset + n] - log_N_i;
			}
		}
	}

	return log_w;
}


std::vector<double> Wham::calculateWeightsForUnbiasedEnsemble() const
{
	return calculateWeights(u_bias_as_other_unbiased_, f_unbiased_);
}


std::vector<double> Wham::calculateWeightsForSimulation(const std::string& data_set_label) const
{
	// Find the index of the simulation (i.e. ensemble)
	auto it = std::find_if( simulations_.begin(), simulations_.end(),
		[=](const Simulation& s) -> bool {
			return s.getDataSetLabel() == data_set_label;
		}
	);
	const int ens = std::distance(simulations_.begin(), it);

	return calculateWeights( u_bias_as_other_[ens], f_bias_opt_[ens] );
}



std::vector<double> Wham::calculateWeights(
	const std::vector<double>& u_bias, const double f_bias) const
{
	int num_in = u_bias.size();
	FANCY_ASSERT(num_in == num_samples_total_, "unexpected number of samples in input");

	std::vector<double> weights(num_samples_total_);
	for ( int n=0; n<num_samples_total_; ++n ) {
		weights[n] = std::exp( f_bias - u_bias[n] - log_dhat_opt_[n] );
	}

	return weights;
}


double Wham::computeBiasingFreeEnergy(const std::vector<double>& u_bias_as_k) const
{
	FANCY_ASSERT( static_cast<int>(u_bias_as_k.size()) == num_samples_total_, "size mismatch" );

	// TODO: optional OpenMP switch?
	std::vector<double> args(num_samples_total_);
	#pragma omp parallel for schedule(static, 8)
	for ( int n=0; n<num_samples_total_; ++n ) {
		args[n] = -u_bias_as_k[n] - log_dhat_opt_[n];
	}

	double f_k = -logSumExp( args );
	return f_k;
}
