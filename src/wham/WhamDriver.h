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
#include "Constants.h"
#include "DataSummary.h"
#include "Distribution.h"
#include "FileSystem.h"
#include "InputParser.h"
#include "OrderParameter.h"
#include "OrderParameterRegistry.h"
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


	// FIXME Delete/move/use?
	struct WhamResults1D
	{
		const OrderParameter* op_ptr = nullptr;  // ptr to OP in question

		Distribution f_x_wham;  // consensus distribution

		// Free energy of turning on the bias (in kBT): Delta F_bias = -ln(Q_i/Q_0)
		std::vector<double> f_bias;

		std::vector<Distribution> f_x_unbiased;  // manually unbiased
		std::vector<Distribution> f_x_shifted;   // manually unbiased and shifted

		// Rebiased distributions and associated entropy
		std::vector<Distribution> f_x_rebiased;
		std::vector<double>       info_entropy;
	};

	// Driver: Manages solving the WHAM equations and printing output
	void run_driver();


 private:
	// Input file
	std::string options_file_;
	ParameterPack input_parameter_pack_;

	std::string data_summary_file_;
	DataSummary data_summary_;

	OrderParameterRegistry op_registry_;

	std::string biases_log_file_;

	std::vector<Simulation> simulations_;
	WhamDriver::WhamOptions wham_options_;

	std::vector<double> f_bias_guess_;
	std::vector<double> f_bias_opt_;

	// Organizes time series data and distributions for each OP
	std::vector<OrderParameter> order_parameters_;

	// Objects that read in and evaluate the bias used in each simulation
	std::vector<Bias> biases_;  // [num_simulations x 1]


	//----- Output -----//

	// OP indices for F(x) to print
	std::vector<int> output_f_x_;

	// OP indices for F(x,y) to print
	std::vector<std::array<int,2>> output_f_x_y_;


	//----- Output Files -----//

	// TODO move to OrderParameter?
	void printWhamResults(const OrderParameter& x) const;

	void print_f_x_y(
		const OrderParameter& x, const OrderParameter& y,
		// Consensus distributions for F(x,y)
		const std::vector<std::vector<double>>& p_x_y_wham,
		const std::vector<std::vector<double>>& f_x_y_wham,
		const std::vector<std::vector<int>>&    sample_counts_x_y
	) const;
};

#endif // ifndef WHAM_DRIVER_H
