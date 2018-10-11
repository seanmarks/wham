/* Wham.h
 *
 * ABOUT: Implements 1D, log-likelihood Weighted Histogram Analysis Method (WHAM) 
 *   - See Hummer & Zhu, J. Comp. Chem. (2011)  (TODO update)
 * NOTES:
 *   - "x" is the generic name for the order parameter in question
 *   - All energies are in units of k_B*T unless otherwise noted
 *   - Currently **assumes** a harmonic bias on the OP of choice
 * TODO
 *   - Allow u_bias-values as input
 *
 * INPUT:
 *   1. wham_options.input
 *      - Key-value pairs
 *          col_y: column of the rewighting variable "y" in the time series files
 *                       used for rewighting
 *
 *   2. data_summary.input
 *      - Each line corresponds to a data set, and takes the following form:
 *          <data_set_label>  <time_series_file(relpath)>  <xtc_file(relpath)>  <t0>  <tf>
 *      - Each time series file contains the following columns:
 *          // TODO update
 * 
 *   3. time_series_files_y.input
 *      - Each line corresponds to a data set, and takes the following form:
 *          <data_set_label>  <path_to_file>
 *      - Each time series file contains the following columns:
 *          <t(ps)>  ... <y> ...
 *      - The variable used as "y" for reweighting is indicated by "col_y" 
 *        (indexed from 1) taken from the WHAM options file
 *      - It is **assumed** that y is not biased
 */

#pragma once
#ifndef WHAM_H
#define WHAM_H

// Standard headers
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>    // unique_ptr
#include <sstream>
#include <string>
#include <vector>

// Library headers
#include "dlib/optimization.h"

// Project headers
#include "Bias.h"
#include "Bins.h"
#include "WhamDlibWrappers.h"

class Wham
{
 public:
	// Types
	using BinStyle = Bins::BinStyle;
	using ColumnVector = dlib::matrix<double,0,1>;

	Wham(
		const std::string& options_file,
		const std::string& data_summary_file,
		const std::string& biasing_parameters_file,
		const std::string& time_series_y_files_list
	);

	// Struct with everything the WHAM algorithm needs to know about a simulation
	struct Simulation 
	{
		std::string data_set_label;
		double x_star;  // for labeling simulations

		double kBT, beta;    // k_B*T and beta = 1/(k_B*T)
		double t_min, t_max; // Sampling range [ps]
		double f_bias;       // Free energy of adding bias (in kBT): f_bias = -ln(Q_i/Q_0)
		double f_bias_guess; // First estimate of biasing free energy

		// Raw output
		std::string data_file;        // File with raw data
		std::vector<double> times;     // t [ps]:            [numSamples x 1]
		std::vector<double> x_samples; // x(t):              [numSamples x 1]
		std::vector<double> u_bias;    // U_bias(t) [k_B*T]: [numSamples x 1]
		std::vector<double> nv_samples;

		// Auxiliary order parameter
		std::string         data_file_y;
		std::vector<double> y_samples;   // y(t):  [numSamples x 1]

		// Order parameter statistics (for e.g. plotting Gaussians)
		double avg_x, var_x;
		double avg_nv, var_nv;

		//----- Histograms -----//

		// - All have size [num_bins_x x 1]
		std::vector<double> x_bins;           // bin centers, x_l
		std::vector<int>    sample_counts;    // Number of samples in each bin, n_{i,l}
		// Biased free energy distribution: beta*F_{biased,i} = -ln[P_i(x_l)]
		std::vector<double> f_biased;
		std::vector<double> p_biased;
		// beta*F_{0,i} = unbiased using data only from simulation i
		std::vector<double> f_unbiased;
		std::vector<double> p_unbiased;
		// Re-apply i-th bias to consensus distribution, F_{0,WHAM}
		std::vector<double> f_rebiased;
		std::vector<double> p_rebiased;
	};

	struct WhamOptions
	{
		// Columns in input files with data (indexed from 0)
		int col_x, col_y;

		double T;  // temperature (K)
		double kBT;
		bool   round_t;  // whether to round times to the nearest ps
		std::string x_name, y_name; // Name(s) of order parameter(s)
	};

	// Stores wham.solve() results
	struct WhamHistogram
	{
		std::vector<double> x_bins;
		std::vector<double> p_wham;
		std::vector<double> f_wham;
		std::vector<double> error_f;        // standard errors
		std::vector<int>    sample_counts; // standard errors
	};


	//----- Log-likelihood solution using dlib -----//

	// Solves the WHAM equations and constructs the consensus histogram
	void solve();	

	// Evaluate objective function
	// - Requires member variables log_N_ and log_M_ be set in order to work
	double evalObjectiveFunction(const ColumnVector& df);

	// Compute derivatives of objective function
	// - Must call evalObjectiveFunction before this one, to set log_sigma_0_
	const ColumnVector evalObjectiveDerivatives(const ColumnVector& df);

	//----- I/O -----//

	// Print sets of "raw" (i.e. non-consensus) histograms for all simulations, including:
	//  - Biased free energy distributions
	//  - Unbiased free energy distributions (using only same-simulation data)
	void printRawDistributions();

 private:
	// Input files
	std::string options_file_;
	std::string data_summary_file_;
	std::string biasing_parameters_file_;
	std::string time_series_y_files_list_;

	std::vector<Wham::Simulation> simulations_;
	Wham::WhamOptions wham_options_;
	Wham::WhamHistogram wham_results_;

	// 
	std::vector<Bias> biases_;

	// Objects for managing histogram bins
	Bins bins_x_;
	Bins bins_y_;

	// Histogram of the total number of samples in each bin, across all simulations
	std::vector<int>    global_sample_counts_; // aka "M_l"


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
	//      = (num. samples from simulation r)/num_total
	std::vector<double> c_;
	std::vector<double> log_c_;

	// The value of the bias for each sample, evaluating using each potential
	//   u_bias_as_other[j][i][r] = u_{bias,r}( x_{j,i} )
	//      j = 1, ..., m    (m = num_simulations)
	//      i = 1, ..., n_j  (n_j = # samples from simulation j)
	//      r = 1, ..., m
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
	std::vector<std::vector<double>> log_sigma_0_;

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

	// Reads the biasing parameters (**assumes** harmonic bias)
	void createBiases(const std::string& biasing_parameters_file);

	// Reads the time series data for the primary (biased) variable, x
	void readRawDataFiles();

	// Reads the time series data for the auxiliary variable, y
	void readAuxiliaryTimeSeries(const std::string& time_series_y_files_list);

	// After reading input files, use to generate histograms of raw data from biased simulations
	void constructRawDataHistograms();


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

	BinStyle parseBinStyle(const std::string& bin_style_token) const;

	// Compute F_0(x) using only the data provided
	// TODO way to merge with compute_consensus_f_x?
	void manually_unbias_f_x(
		const std::vector<double>& x,       // samples from a single simulation
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
		const std::vector<std::vector<double>>& x,
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
	// TODO generalize to n dimensions and combine with compute_consensus_f_x
	void compute_consensus_f_x_y(
		const std::vector<std::vector<double>>& x,
		const std::vector<std::vector<double>>& y,
		const std::vector<double>& f_opt,  // consensus free energies to use
		const int k,  // index of simulation ensemble (-1 --> unbiased)
		const Bins& bins_x,
		const Bins& bins_y,
		// Consensus distributions for F_k(x,y)
		std::vector<std::vector<double>>& x_grid,
		std::vector<std::vector<double>>& y_grid,
		std::vector<std::vector<double>>& p_x_y_wham,
		std::vector<std::vector<double>>& f_x_y_wham,
		std::vector<std::vector<int>>&    sample_counts,
		//
		std::vector<double>& p_y_wham,
		std::vector<double>& f_y_wham
	) const;

	// If the probability p > 0, print free energy f; else print "nan"
	std::ofstream& print_free_energy(std::ofstream& ofs, const double f, const double p) const {
		if ( p > 0.0 ) { ofs << f; }
		else           { ofs << "nan";  }
		return ofs;
	};

	//----- Constants -----//

	// Boltzmann's constant [kJ/mol]
	static constexpr double K_B_ = 8.314e-3;

	// Always good to have some double-precision PI lying around, just in case
	static constexpr double PI_ = 3.14159265358979323846;

	// For functions that take the index of an ensemble, this value is
	// used to indicate the *unbiased* ensemble
	static constexpr int unbiased_ensemble_index_ = -1;

	//----- File Path Manipulations  -----//

	std::string get_realpath(const std::string& path) const 
	{
		if ( path.length() < 1 ) {
			throw std::runtime_error("get_realpath() was given an empty path");
		}

		// Use POSIX realpath()
		char* buffer = realpath(&path[0], nullptr);
		if ( buffer == nullptr ) {
			throw std::runtime_error("Error resolving path \"" + path + "\"");
		}

		// Move the path to a std string and clean up
		std::string resolved_path(buffer);
		free(buffer); buffer = nullptr;

		return resolved_path;
	}

	// Taken from:
	//   C++ Cookbook by Jeff Cogswell, Jonathan Turkanis, Christopher Diggins, D. Ryan Stephens
	std::string get_dir(const std::string& full_path) const 
	{
		// Separator
		char sep = '/';
//#ifdef _WIN32
//		sep = '\\';
//#endif

		size_t i = full_path.rfind(sep, full_path.length());
		if ( i != std::string::npos ) {
			return full_path.substr(0, i);
		}
		else {
			return ".";
		}
	}
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
