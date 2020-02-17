
#pragma once
#ifndef WHAM_H
#define WHAM_H

// Standard headers
#include <array>
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

// Library headers
#include "dlib/optimization.h"

// Project headers
#include "Bias.h"
#include "Bins.h"
#include "Distribution.h"
#include "FileSystem.h"
#include "InputParser.h"
#include "OrderParameter.h"
#include "Simulation.h"
#include "WhamDlibWrappers.h"

class Wham
{
 public:
	Wham(
		const std::vector<Simulation>& simulations,
		const std::vector<OrderParameter>& order_parameters,
		const std::vector<Bias>& biases,
		//const std::vector<std::vector<double>>& u_bias_as_other,
		const double tol
		// TODO
	);

	// TODO Make private and run as part of constructor?
	void solveWhamEquations(
		const std::vector<double>& f_init,
		std::vector<double>& f_bias_opt
	);


	//----- Objective Function -----//

	using ColumnVector = dlib::matrix<double,0,1>;

	double evalObjectiveFunction(const ColumnVector& df) const;

	// Compute gradient of objective function
	// - Note:
	//   - grad  = gradient wrt. free energy differences between neighboring windows
	//   - dA_df = gradient wrt. biasing free energies themselves
	const ColumnVector evalObjectiveDerivatives(const ColumnVector& df) const;

 private:
	const std::vector<Simulation>&     simulations_;
	const std::vector<OrderParameter>& order_parameters_;
	const std::vector<Bias>&           biases_;

	// Maps names of OrderParameters top their indices in the order_parameters_ vector
	std::map<std::string, int> map_op_names_to_indices_;

	// The value of the bias for each sample, evaluating using each potential (i.e. in each ensemble):
	//    u_{bias,r}( x_{j,i} )
	//      r = 0, ..., m-1    (m = number of biasing potentials/ensembles)
	//      j = 0, ..., m-1    (m = num_simulations)
	//      i = 0, ..., n_j-1  (n_j = # samples from simulation j)
	// - For each bias r = 0, ..., m-1:
	//     u_bias_as_other_[r] = [(data_sim_0), ..., (data_sim_j), ... , (data_sim_m)]
	//                                                     |
	//                                              [(n_j samples)]
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


	// These are computed by Wham.evalObjectiveFunction and saved for re-use 
	// in Wham.evalObjectiveDerivatives
	// - FIXME: dangerous?

	// Factors related to the weight given to each data sample in the *unbiased* ensemble
	// - Computed as a log-sum-exp
	// - Formula:
	//     sigma_0(x_{j,i}) = sum_{r=1}^m exp{ log(c_r) + f_r - u_{bias,r}(x_{j,i}) }
	mutable std::vector<double> f_bias_last_;
	mutable std::vector<double> log_sigma_unbiased_;

	// Buffers
	mutable std::vector<double> args_buffer_;               // for log_sum_exp
	mutable std::vector<double> log_sigma_k_, minus_log_sigma_k_;
	mutable std::vector<std::vector<double>> minus_log_sigma_k_binned_;  // for binning samples
	mutable std::vector<int> sample_bins_;


	//----- Setup -----//

	void setup();

	// For each sample (across all simulations), evaluate the bias that would be
	// felt under each simulation's potential
	void evaluateBiases();


	//----- Helper Functions -----//

	// Compute log( sigma_k(x_{j,i}) ) for the given set of biasing free energies
	// - k: index of ensemble to which the weights correspond
	//   - f_k: free energy of turning on kth bias
	void compute_log_sigma(
		const std::vector<std::vector<double>>& u_bias_as_other,
		const std::vector<double>&              f,
		const std::vector<double>&              u_bias_as_other_k,
		const double                            f_k,
		// Output
		std::vector<double>& log_sigma
	) const;

	// Returns the logarithm of a sum of exponentials
	// - input: arguments of exponentials
	double log_sum_exp(const std::vector<double>& args) const;

	// Convert between the free energies of turning on the bias (f) and the free energy 
	// *differences* between windows (df), assuming f[0] = f0 = 0.0
	void convert_f_to_df(const std::vector<double>& f, ColumnVector& df) const;
	void convert_df_to_f(const ColumnVector& df, std::vector<double>& f) const;


	//----- Consensus Distributions -----//

	// Compute the consensus distribution F_k^{WHAM}(x), where k is the index
	// of any simulation under consideraton (or use -1 to get unbiased ensemble results)
	void compute_consensus_f_x(
		const std::vector<TimeSeries>& x,
		const std::vector<std::vector<double>>& u_bias_as_other,
		const std::vector<double>&              f_opt,  // consensus free energies to use
		const std::vector<double>&              u_bias_as_other_k,
		const double                            f_k,
		const Bins& bins_x,
		// Consensus distribution for x in ensemble k
		Distribution& wham_distribution_x
	) const;

	// Compute the consensus distribution F_k^{WHAM}(x,y), where k is the index
	// of any simulation under consideraton (or use -1 to get unbiased ensemble results)
	// - Grids use 'ij' organization, i.e. grid[i][j] corresponds to (x_bins[i], y_bins[j])
	// TODO generalize to n dimensions and combine with compute_consensus_f_x
	void compute_consensus_f_x_y(
		const std::vector<TimeSeries>& x,
		const std::vector<TimeSeries>& y,
		const std::vector<std::vector<double>>& u_bias_as_other,
		const std::vector<double>&              f_opt,  // consensus free energies to use
		const std::vector<double>&              u_bias_as_other_k,
		const double                            f_k,
		const Bins& bins_x,
		const Bins& bins_y,
		// Consensus distributions for F_k(x,y)
		std::vector<std::vector<double>>& p_x_y_wham,
		std::vector<std::vector<double>>& f_x_y_wham,
		std::vector<std::vector<int>>&    sample_counts_x_y
	) const;
};

#endif // ifndef WHAM_H