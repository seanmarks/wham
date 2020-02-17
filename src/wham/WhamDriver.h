/* WhamDriver.h
 *
 * ABOUT: Implements the Unbinned Weighted Histogram Analysis Method (UWHAM) in N dimensions
 *   - Equations are solved using log-likelihood maximation approach
 *     - See Tan, Gallicchio, Lapelosa, & Levy (J. Chem. Phys. 2012)
 *   - See also:
 *     - Zhu & Hummer (J. Comp. Chem. 2011)
 *     - Souaille & Roux (Comp. Phys. Comm. 2001)
 * NOTES:
 *   - "x" and "y" are generic names for the order parameters in question
 *   - All energies are in units of k_B*T unless otherwise noted
 * TODO
 *   - Allow u_bias-values as input
 *   - Check for internal consistency: 
 *     - all time series for all order parameters
 *     - biasing parameters
 *
 * INPUT: (TODO update)
 *   1. wham_options.input
 *      - Key-value pairs
 *
 *   2. data_summary.input
 *      - Each line corresponds to a data set, and usually takes the following form:
 *          <data_set_label>  <time_series_file(relpath)>  <xtc_file(relpath)>  <t0>  <tf>
 *      - Each time series file contains the following columns:
 */

#pragma once
#ifndef WHAM_DRIVER_H
#define WHAM_DRIVER_H

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
#include "Wham.h"

class WhamDriver
{
 public:
	// Types
	using ColumnVector = dlib::matrix<double,0,1>;

	WhamDriver(
		const std::string& options_file
	);

	// TODO Move to standalone member variables?
	struct WhamOptions
	{
		WhamOptions(): floor_t(true), tol(1.0e-7) {};

		double T;        // default temperature (K) assumed for simulations
		double kBT;      // k_B*T
		bool   floor_t;  // whether to round times down to the nearest ps
		double tol;      // tolerance for solver
	};


	// Driver: Manages solving the WHAM equations and printing output
	void run_driver();


 private:
	// Input file
	std::string options_file_;
	ParameterPack input_parameter_pack_;

	std::string data_summary_file_;
	int col_data_label_;         // column with data label
	int col_t_min_, col_t_max_;  // columns with production phase bounds
	int col_T_;                  // column with temperature

	std::string biases_log_file_;

	std::vector<Simulation> simulations_;
	WhamDriver::WhamOptions wham_options_;

	std::vector<double> f_bias_guess_;
	std::vector<double> f_bias_opt_;

	// Organizes time series data and distributions for each OP
	std::vector<OrderParameter> order_parameters_;

	// Maps names of OrderParameters top their indices in the order_parameters_ vector
	std::map<std::string, int> map_op_names_to_indices_;

	// Objects that read in and evaluate the bias used in each simulation
	std::vector<Bias> biases_;  // [num_simulations x 1]


	//----- Precompute/save useful quantities for speed -----//

	std::vector<int> num_samples_per_simulation_;   // number of samples from each simulation
	int num_samples_total_;  // Total number of samples (across all simulations)
	double inv_num_samples_total_;  // precompute for speed

	// c[r] = fraction of samples from simulation 'r'
	//      = (num. samples from simulation r)/(num. total)
	std::vector<double> c_;
	std::vector<double> log_c_;

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

	// u_bias_as_other for the unbiased ensemble (all zeros)
	std::vector<double> u_bias_as_other_unbiased_;
	const double f_unbiased_ = 0.0;  // Free energy of going to the *unbiased* ensemble is zero


	//----- Working Variables -----//

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


	//----- Output -----//

	// OP indices for F(x) to print
	std::vector<int> output_f_x_;

	// OP indices for F(x,y) to print
	std::vector<std::array<int,2>> output_f_x_y_;


	//----- Setup -----//

	// TODO Move to DataSummary class
	// Reads the data summary
	// - Used to determine the number of simulations
	// - Contains production phase bounds [t0, tf] [ps]
	void readDataSummary(const std::string& data_summary_file);

	// After reading input files, use this to analyze the raw data
	// and populate the OrderParameter object
	//void manuallyUnbiasDistributions(OrderParameter& x);


	/*
	//----- Helper Functions -----//

	// For each sample (across all simulations), evaluate the bias that would be
	// felt under each simulation's potential
	void evaluateBiases();
	*/


	//----- Output Files -----//

	/*
	// TODO delete
	// Print sets of "raw" (i.e. non-consensus) histograms for each simulation 
	// time series,including:
	//  - Biased free energy distributions
	//  - Unbiased free energy distributions (using only same-simulation data)
	void printRawDistributions(const OrderParameter& x) const;
	*/

	// TODO move to OrderParameter?
	void printWhamResults(const OrderParameter& x) const;

	void print_f_x_y(
		const OrderParameter& x, const OrderParameter& y,
		// Consensus distributions for F(x,y)
		const std::vector<std::vector<double>>& p_x_y_wham,
		const std::vector<std::vector<double>>& f_x_y_wham,
		const std::vector<std::vector<int>>&    sample_counts_x_y
	) const;


	//----- Constants -----//

	// Boltzmann's constant [kJ/mol]
	static constexpr double K_B_ = 8.314e-3;

	// Always good to have some double-precision PI lying around, just in case
	static constexpr double PI_ = 3.14159265358979323846;
};

#endif // ifndef WHAM_DRIVER_H
