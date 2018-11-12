/* Wham.h
 *
 * ABOUT: Implements 1D, log-likelihood Weighted Histogram Analysis Method (WHAM) 
 *   - See Hummer & Zhu, J. Comp. Chem. (2011)  (TODO update)
 * NOTES:
 *   - "x" is the generic name for the order parameter in question
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
#include "FileSystem.h"
#include "InputParser.h"
#include "OrderParameter.h"
#include "WhamDlibWrappers.h"

class Wham
{
 public:
	// Types
	using ColumnVector = dlib::matrix<double,0,1>;

	Wham(
		const std::string& options_file
	);

	// Struct with simulation-specific information
	struct Simulation 
	{
		std::string data_set_label;
		double x_star;  // for labeling simulations

		double kBT, beta;    // k_B*T and beta = 1/(k_B*T)
		double t_min, t_max; // Sampling range [ps]
		double f_bias;       // Free energy of adding bias (in kBT): f_bias = -ln(Q_i/Q_0)
		double f_bias_guess; // First estimate of biasing free energy
	};

	struct WhamOptions
	{
		WhamOptions(): floor_t(true), tol(1.0e-7) {};

		double T;        // default temperature (K) assumed for simulations
		double kBT;      // k_B*T
		bool   floor_t;  // whether to round times down to the nearest ps
		double tol;      // tolerance for solver
	};

	// Stores wham.solve() results
	struct WhamResults
	{
		std::vector<double> x_bins;
		std::vector<double> p_x_wham;
		std::vector<double> f_x_wham;
		std::vector<double> error_f;        // standard errors
		std::vector<int>    sample_counts;

		std::vector<double> f_opt;         // consensus free energies
		std::vector<double> info_entropy;  // entropy between f_x_biased and f_x_rebiased
	};


	//----- Log-likelihood solution using dlib -----//

	// Solves the WHAM equations and constructs the consensus histogram
	void solve();	

	// Evaluate objective function
	// - Requires member variables log_N_ and log_M_ be set in order to work
	double evalObjectiveFunction(const ColumnVector& df) const;

	// Compute derivatives of objective function
	// - Must call evalObjectiveFunction before this one, to set log_sigma_0_
	const ColumnVector evalObjectiveDerivatives(const ColumnVector& df) const;

	//----- I/O -----//

	// TODO Make private? These rely on internal state variables
	// TODO multiple OPs
	void printRawDistributions() const {
		int num_ops = order_parameters_.size();
		for ( int i=0; i<num_ops; ++i ) {
			printRawDistributions( order_parameters_[i] );
		}
	}

	// Print sets of "raw" (i.e. non-consensus) histograms for all simulations, including:
	//  - Biased free energy distributions
	//  - Unbiased free energy distributions (using only same-simulation data)
	void printRawDistributions(const OrderParameter& x) const;

	void printWhamResults(const OrderParameter& x) const;

	void print_f_x_y(
		const OrderParameter& x, const OrderParameter& y,
		// Consensus distributions for F(x,y)
		const std::vector<std::vector<double>>& p_x_y_wham,
		const std::vector<std::vector<double>>& f_x_y_wham,
		const std::vector<std::vector<int>>&    sample_counts_x_y
	) const;

	// If the probability p > 0, print free energy f; else print "nan"
	void print_free_energy(std::ofstream& ofs, const double f, const double p) const {
		if ( p > 0.0 ) { ofs << f; }
		else           { ofs << "nan";  }
	};


 private:
	// Input file
	std::string options_file_;
	ParameterPack input_parameter_pack_;

	std::string data_summary_file_;
	int col_data_label_;         // column with data label
	int col_t_min_, col_t_max_;  // columns with production phase bounds
	int col_T_;                  // column with temperature

	std::string biasing_parameters_file_;

	std::vector<Wham::Simulation> simulations_;
	Wham::WhamOptions wham_options_;
	Wham::WhamResults wham_results_;

	// Organizes time series data and distributions for each OP
	std::vector<OrderParameter> order_parameters_;

	// Objects that read in and evaluate the bias used in each simulation
	std::vector<Bias> biases_;  // [num_simulations x 1]

	// Indices of the order parameters that were biased
	// - Used to evaluate the bias from each sample under each ensemble
	std::vector<int> biased_order_parameters_;


	//----- Precompute/save useful quantities for speed -----//

	// TODO linearize
	/*
	std::vector<int> num_samples_;         // number of samples from each simulation
	std::vector<int> simulation_offsets_;  // beginning of each sim's data in time series arrays
	*/

	// Total number of samples (across all simulations)
	int num_samples_total_;
	double inv_num_samples_total_;

	// c[r] = fraction of samples from simulation 'r'
	//      = (num. samples from simulation r)/(num. total)
	std::vector<double> c_;
	std::vector<double> log_c_;

	// The value of the bias for each sample, evaluating using each potential
	//   u_bias_as_other[j][i][r] = u_{bias,r}( x_{j,i} )
	//      j = 1, ..., m    (m = num_simulations)
	//      i = 1, ..., n_j  (n_j = # samples from simulation j)
	//      r = 1, ..., m    (m = number of biasing potentials)
	// TODO linearize for speed
	std::vector<std::vector<std::vector<double>>> u_bias_as_other_;


	//----- Working Variables -----//

	// These are computed by Wham.evalObjectiveFunction and saved for re-use 
	// in Wham.evalObjectiveDerivatives

	// Factors related to the weight given to each data sample in the
	// *unbiased* ensemble
	// - Computed as a log-sum-exp
	// - Formula:
	//     sigma_0(x_{j,i}) = sum_{r=1}^m exp{ log(c_r) + f_r - u_{bias,r}(x_{j,i}) }
	// TODO linearize for speed
	mutable std::vector<double> f_bias_last_;
	mutable std::vector<std::vector<double>> log_sigma_0_;

	// Buffers
	mutable std::vector<double> args_buffer_;
	mutable std::vector<std::vector<double>> log_sigma_k_; 
	mutable std::vector<std::vector<double>> minus_log_sigma_k_binned_;  // for binning samples


	//----- Setup -----//

	// Set WHAM options from file
	void parseOptionsFile(const std::string& options_file);

	// Reads the data summary
	// - Used to determine the number of simulations
	// - **WARNING** clears all stored data
	void readDataSummary(const std::string& data_summary_file);

	// Create Bias objects according to the given file
	void createBiases(const std::string& biasing_parameters_file);

	// After reading input files, use this to analyze the raw data
	// and populate the OrderParameter object
	void analyzeRawData(OrderParameter& x);


	//----- Helper Functions -----//

	// For each sample (across all simulations), evaluate the bias that would be
	// felt under each simulation's potential
	void evaluateBiases();

	// Compute log( sigma_k(x_{j,i}) ) for the given set of biasing free energies
	// - k: index of ensemble to which the weights correspond
	//   - k < 0 --> compute for unbiased ensemble
	// - f_k: free energy of turning on kth bias
	void compute_log_sigma(
		const std::vector<double>& f,
		const std::vector<std::vector<std::vector<double>>>& u_bias_as_other,
		const int k,
		std::vector<std::vector<double>>& log_sigma
	) const;

	// Returns the logarithm of a sum of exponentials
	// - input: arguments of exponentials
	double log_sum_exp(const std::vector<double>& args) const;

	template<typename T>
	double average(const std::vector<T>& x) const;

	template<typename T>
	double variance(const std::vector<T>& x) const;

	// Convert between the free energies of turning on the bias (f) and the free energy 
	// *differences* between windows (df), assuming f[0] = f0 = 0.0
	void convert_f_to_df(const std::vector<double>& f, Wham::ColumnVector& df) const;
	void convert_df_to_f(const Wham::ColumnVector& df, std::vector<double>& f) const;

	// Compute F_0(x) using only the data provided
	// TODO way to merge with compute_consensus_f_x?
	void manually_unbias_f_x(
		const TimeSeries& x,   // samples from a single simulation
		const std::vector<double>& u_bias,  // bias corresponding to x-samples
		const double f,                     // free energy of biasing (usually a guess)
		const Bins& bins_x,
		// Output
		std::vector<double>& p_x, 
		std::vector<double>& f_x,
		std::vector<int>& sample_counts
	) const;

	// Compute the consensus distribution F_k^{WHAM}(x), where k is the index
	// of any simulation under consideraton (or use -1 to get unbiased ensemble results)
	void compute_consensus_f_x(
		const std::vector<TimeSeries>& x,
		const std::vector<std::vector<std::vector<double>>>& u_bias_as_other,
		const std::vector<double>& f_opt,  // consensus free energies to use
		const int k,                       // index of simulation ensemble (-1 --> unbiased)
		const Bins& bins_x,
		// Consensus distributions for x in ensemble k
		std::vector<double>& p_x_wham,
		std::vector<double>& f_x_wham,
		std::vector<int>& sample_counts
	) const;

	// Compute the consensus distribution F_k^{WHAM}(x,y), where k is the index
	// of any simulation under consideraton (or use -1 to get unbiased ensemble results)
	// - Grids use 'ij' organization, i.e. grid[i][j] corresponds to (x_bins[i], y_bins[j])
	// TODO generalize to n dimensions and combine with compute_consensus_f_x
	void compute_consensus_f_x_y(
		const std::vector<TimeSeries>& x,
		const std::vector<TimeSeries>& y,
		const std::vector<std::vector<std::vector<double>>>& u_bias_as_other,
		const std::vector<double>& f_opt,  // consensus free energies to use
		const int k,  // index of simulation ensemble (-1 --> unbiased)
		const Bins& bins_x,
		const Bins& bins_y,
		// Consensus distributions for F_k(x,y)
		std::vector<std::vector<double>>& p_x_y_wham,
		std::vector<std::vector<double>>& f_x_y_wham,
		std::vector<std::vector<int>>&    sample_counts_x_y
	) const;

	//----- Constants -----//

	// Boltzmann's constant [kJ/mol]
	static constexpr double K_B_ = 8.314e-3;

	// Always good to have some double-precision PI lying around, just in case
	static constexpr double PI_ = 3.14159265358979323846;

	// For functions that take the index of an ensemble, this value is
	// used to indicate the *unbiased* ensemble
	static constexpr int unbiased_ensemble_index_ = -1;
};


//----- Templated Helper Functions -----//


template<typename T>
double Wham::average(const std::vector<T>& x) const
{
	double avg = 0.0;
	int num_values = x.size();
	if ( num_values == 0 ) {
		return avg;
	}

	for ( int i=0; i<num_values; ++i ) {
		avg += x[i];
	}
	avg /= static_cast<double>(num_values);

	return avg;
}


template<typename T>
double Wham::variance(const std::vector<T>& x) const
{
	double var = 0.0;
	int num_values = x.size();
	if ( num_values == 0.0 ) {
		return var; // FIXME appropriate to return 0 in this case?
	}

	double dx, avg_x = average(x);
	for ( int i=0; i<num_values; ++i ) {
		dx = x[i] - avg_x;
		var += dx*dx;
	}
	var /= static_cast<double>(num_values);

	return var;
}



#endif // WHAM_H
