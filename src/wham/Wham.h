// Wham
// - Solves binless WHAM equations via log-likelihood maximization
//   - Minimization uses BFGS implementation from dlib library
// - Computes consensus estimates using optimal biasing free energies
// - References
//   - Chodera, ..., Dill (JCTC 2007)
//   - Kong, ..., Tan (J. R. Statist. Soc. B 2003)
//   - Shirts & Chodera (JCP 2008)
//   - Souaille & Roux (Comp. Phys. Comm. 2001)
//   - Tan, Gallicchio, Lapelosa, & Levy (J. Chem. Phys. 2012)
//   - Zhu & Hummer (J. Comp. Chem. 2011)

#pragma once
#ifndef WHAM_H
#define WHAM_H

// Standard headers
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>    // unique_ptr
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

// dlib
#include "dlib/matrix.h"
#include "dlib/optimization.h"

// Project headers
#include "Bias.h"
#include "Bins.h"
#include "DataSummary.h"
#include "Distribution.h"
#include "FileSystem.h"
#include "GptlWrappers.h"
#include "InputParser.h"
#include "OpenMP.h"
#include "OrderParameter.h"
#include "OrderParameterRegistry.h"
#include "Simulation.h"
#include "WhamDlibWrappers.h"

class Wham
{
 public:
	// Runtime-sized arrays
	using Matrix       = dlib::matrix<double,0,0>;
	using ColumnVector = dlib::matrix<double,0,1>;
	using RowVector    = dlib::matrix<double,1,0>;

	Wham(
		const DataSummary& data_summary,
		const OrderParameterRegistry& op_registry,
		const std::vector<Simulation>& simulations,
		const std::vector<OrderParameter>& order_parameters,
		const std::vector<Bias>& biases,
		//const std::vector<std::vector<double>>& u_bias_as_other,
		const std::vector<double>& f_bias_guess,
		const double tol
	);

	// TODO Make private and run as part of constructor?


	//----- Objective Function -----//

	double evalObjectiveFunction(const ColumnVector& df) const;

	// Compute gradient of objective function
	// - Note:
	//   - grad  = gradient wrt. free energy differences between neighboring windows
	//   - dA_df = gradient wrt. biasing free energies themselves
	const ColumnVector evalObjectiveDerivatives(const ColumnVector& df) const;


	//----- WHAM Estimators -----//

	const std::vector<double>& get_f_bias_opt()       const { return f_bias_opt_;       }
	const std::vector<double>& get_error_f_bias_opt() const { return error_f_bias_opt_; }

	// "Manually" unbias the distributions for the given OrderParameter (i.e. using only
	// each individual simulation's data, not the consensus estimates)
	std::vector<Distribution> manuallyUnbiasDistributions(const std::string& op_name) const;

	// TODO descriptions
	Distribution compute_consensus_f_x_unbiased(
		const std::string& op_name
	) const;
	Distribution compute_consensus_f_x_rebiased(
		const std::string& op_name,
		const std::string& data_set_label
	) const;

	// TODO descriptions
	void compute_consensus_f_x_y_unbiased(
		const std::string& x_name,
		const std::string& y_name,
		// Output
		std::vector<std::vector<double>>& p_x_y_wham, 
		std::vector<std::vector<double>>& f_x_y_wham,
		std::vector<std::vector<int>>& sample_counts_x_y
	) const;

 private:
	// Objects owned by the driver
	const DataSummary&                 data_summary_;
	const OrderParameterRegistry&      op_registry_;
	const std::vector<Simulation>&     simulations_;
	const std::vector<OrderParameter>& order_parameters_;
	const std::vector<Bias>&           biases_;

	// The value of the bias for each sample, evaluating using each potential (i.e. in each ensemble):
	//    u_{bias,r}( x_{j,i} )
	//      r = 0, ..., m-1    (m = number of biasing potentials/ensembles)
	//      j = 0, ..., m-1    (m = num_simulations)
	//      i = 0, ..., n_j-1  (n_j = # samples from simulation j)
	// - For each bias r = 0, ...,
	//     u_bias_as_other_[r] = [(data_sim_0), ..., (data_sim_j), ... , (data_sim_m)]
	//                                                     |
	//                                                N_j samples
	//   - Overall length of u_bias_as_other_[r]:  N = sum_{j=0}^{m} N_j
	std::vector<std::vector<double>> u_bias_as_other_;
	// - For each simulation j = 0, ..., m-1:
	//     simulation_data_ranges_[j] = indices (first, end) for the 'j'th simulation's data
	//                                  in u_bias_as_other_ (for each 'r')
	std::vector<std::pair<int,int>> simulation_data_ranges_;

	// Solver tolerance
	const double tol_ = 1.0e-7;

	// Track number of samples
	std::vector<int> num_samples_per_simulation_;   // number of samples from each simulation
	int num_samples_total_;  // Total number of samples (across all simulations)
	double inv_num_samples_total_;  // precompute for speed

	// c[r] = fraction of samples from simulation 'r'
	//      = (num. samples from simulation r)/(num. total)
	std::vector<double> c_;
	std::vector<double> log_c_;

	// u_bias_as_other for the unbiased ensemble (all zeros)
	std::vector<double> u_bias_as_other_unbiased_;
	const double f_unbiased_ = 0.0;  // Free energy of going to the *unbiased* ensemble is zero

	// Free energy of turning on each bias
	std::vector<double> f_bias_guess_;  // initial guess
	std::vector<double> f_bias_opt_;    // optimal
	std::vector<double> error_f_bias_opt_;

	//----- Working Variables -----//


	// TODO document
	std::vector<double> log_dhat_;

	Matrix w_;
	Matrix wT_w_;

	mutable std::vector<double> log_dhat_tmp_;


	// Factors related to the weight given to each data sample in the *unbiased* ensemble
	// - Computed as a log-sum-exp
	// - Formula:
	//     sigma_0(x_{j,i}) = sum_{r=1}^m exp{ log(c_r) + f_r - u_{bias,r}(x_{j,i}) }
	// - These are computed by Wham.evalObjectiveFunction and saved for re-use 
	//   in Wham.evalObjectiveDerivatives
	//   - The dlib optimization algorithm used guarantees that these functions are called in pairs,
	//     with the same order each time 
	mutable std::vector<double> f_bias_last_;
	mutable std::vector<double> log_sigma_unbiased_;

	// Buffers
	mutable std::vector<double> log_sigma_k_, minus_log_sigma_k_;
	mutable std::vector<std::vector<double>> minus_log_sigma_k_binned_;  // for binning samples
	mutable std::vector<std::vector<double>> args_buffers_;  // for log_sum_exp during optimization


	//----- Solve WHAM Equations -----//

	void setup();

	// For each sample (across all simulations), evaluate the bias that would be
	// felt under each simulation's potential
	void evaluateBiases();

	std::vector<double> solveWhamEquations(const std::vector<double>& f_bias_guess);


	//----- Helper Functions -----//

	// Compute log( sigma_k(x_{j,i}) ) for the given set of biasing free energies
	// - These correspond to the weights given to each sample x_{j,i} in the kth ensemble
	//   - f_k: free energy of turning on kth bias
	// *** NOT thread safe: encloses an OpenMP region
	void compute_log_sigma(
		const std::vector<std::vector<double>>& u_bias_as_other,
		const std::vector<double>&              f,
		const std::vector<double>&              u_bias_as_k,
		const double                            f_k,
		// Output
		std::vector<double>& log_sigma
	) const;

	// Compute log( \hat{D}(x_n) ) for the given set of biasing free energies
	// - These are the common denominator of the consensus weights given to each sample
	// *** NOT thread safe: encloses an OpenMP region
	void compute_log_dhat(
		const std::vector<std::vector<double>>& u_bias_as_other,
		const std::vector<double>&              fhat,
		// Output
		std::vector<double>& log_dhat
	) const;

	// Returns the logarithm of a sum of exponentials:
	//      log(S) = sum_{i=1}^N exp(args[i])
	// - input: arguments of exponentials
	// - THREAD_SAFE
	double log_sum_exp(const std::vector<double>& args) const;

	// Returns the weighted sum of a series of exponentials
	//      log(S) = sum_{i=1}^N weights[i]*exp(args[i])
	// - THREAD_SAFE
	double weighted_sum_exp(
		const std::vector<double>& args,    // arguments of exponentials
		const std::vector<double>& weights  // not necessarily positive
	) const;

	// Convert between the free energies of turning on the bias (f) and the free energy 
	// *differences* between windows (df), assuming f[0] = f0 = 0.0
	void convert_f_to_df(const std::vector<double>& f, ColumnVector& df) const;
	void convert_df_to_f(const ColumnVector& df, std::vector<double>& f) const;


	//----- Individual-Simulation Estimates -----//

	// TODO way to merge with compute_consensus_f_x?
	void manually_unbias_f_x(
		const TimeSeries&          x,
		const std::vector<double>& u_bias,
		const double               f,
		const Bins&                bins_x,
		// Output
		Distribution&              unbiased_distribution_x
	) const;


	//----- Consensus Estimates -----//

	// Free energy of turning on arbitrary bias 'k', given the biasing potential
	// of each sample evaluated in that ensemble
	double compute_consensus_f_k(
		const std::vector<double>& u_bias_as_k
	) const;

	// Compute the consensus distribution F_k^{WHAM}(x), where k is the index
	// of any simulation under consideraton (or use -1 to get unbiased ensemble results)
	void compute_consensus_f_x(
		const OrderParameter& x,
		const std::vector<std::vector<double>>& u_bias_as_other,
		const std::vector<double>&              f_opt,  // consensus free energies to use
		const std::vector<double>&              u_bias_as_k,
		const double                            f_k,
		// Consensus distribution for x in ensemble k
		Distribution& wham_distribution_x
	) const;

	// Compute the consensus distribution F_k^{WHAM}(x,y), where k is the index
	// of any simulation under consideraton (or use -1 to get unbiased ensemble results)
	// - Grids use 'ij' organization, i.e. grid[i][j] corresponds to (x_bins[i], y_bins[j])
	// TODO generalize to n dimensions and combine with compute_consensus_f_x
	void compute_consensus_f_x_y(
		const OrderParameter& x,
		const OrderParameter& y,
		const std::vector<std::vector<double>>& u_bias_as_other,
		const std::vector<double>&              f_bias_opt,  // consensus free energies to use
		const std::vector<double>&              u_bias_as_k,
		const double                            f_bias_k,
		// Consensus distributions for F_k(x,y)
		std::vector<std::vector<double>>& p_x_y_wham,
		std::vector<std::vector<double>>& f_x_y_wham,
		std::vector<std::vector<int>>&    sample_counts_x_y
	) const;

	// Compute the asymptotic covariance matrix, theta, for the consensus biasing free energies
	// (plus any extra observables: 'w' may be an augmented set of weights)
	void compute_cov_matrix(
		const Matrix& w,                     // W, matrix of weights
		const Matrix& wT_w,                  // W^T * W
		const std::vector<int>& num_samples,  // per state
		Matrix& theta
	) const;

	// Compute the asymptotic covariance matrix, theta, for the consensus biasing free energies
	// - Note that 'w' may be an augmented set of weights
	void compute_error_avg_x_k(
		const std::vector<double>& x,            // samples
		const std::vector<double>& u_bias_as_k,  // biasing free energies for target ensemble
		// Output
		Matrix& w,                      // W_aug, augmented matrix of weights
		Matrix& wT_w,                   // W_aug^T * W_aug
		std::vector<int>& num_samples,  // per state (augmented)
		Matrix& theta
	) const;

	// Smallest value to for which exp(x) will not underflow (about -36)
	static constexpr double MIN_DBL_FOR_EXP = log( std::numeric_limits<double>::epsilon() );


	//----- GPTL -----//

	using Timer = GPTL::Timer;

	// Core functions
	mutable Timer setup_timer_     = Timer("Wham::setup");
	mutable Timer biases_timer_    = Timer("Wham::evaluate_biases");
	mutable Timer solve_timer_     = Timer("Wham::solve");
	mutable Timer objective_timer_ = Timer("Wham::objective");
	mutable Timer gradient_timer_  = Timer("Wham::gradient");
	mutable Timer gradient_omp_timer_ = Timer("Wham::gradient_omp");

	// Output functions
	mutable Timer f_x_timer_        = Timer("Wham::consensus_f_x");
	mutable Timer f_x_sort_timer_   = Timer("Wham::consensus_f_x_sort");
	mutable Timer f_x_y_timer_      = Timer("Wham::consensus_f_x_y");
	mutable Timer f_x_y_sort_timer_ = Timer("Wham::consensus_f_x_y_sort");

	// Low-level, expensive functions
	mutable Timer log_sigma_timer_        = Timer("Wham::compute_log_sigma");
	mutable Timer log_sigma_omp_timer_    = Timer("Wham::compute_log_sigma_omp");
	mutable Timer log_sum_exp_timer_      = Timer("Wham::compute_log_sum_exp");
	mutable Timer weighted_sum_exp_timer_ = Timer("Wham::weighted_sum_exp");

	mutable Timer log_dhat_timer_     = Timer("Wham::compute_log_dhat");
	mutable Timer log_dhat_omp_timer_ = Timer("Wham::compute_log_dhat_omp");


	//----- Helper Functions -----//

	// FIXME DEBUG
	// Frobenius norm of a real-valued matrix: sqrt[ sum_{i,j}^N m(i,j)^2 ]
	static double norm_frob(const Matrix& m) {
		return sqrt( dlib::sum(dlib::pointwise_multiply(m, m)) );
	}
};

#endif // ifndef WHAM_H
