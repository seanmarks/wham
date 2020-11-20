// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#pragma once
#ifndef WHAM_DRIVER_H
#define WHAM_DRIVER_H

#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

// Library headers
#include "dlib/optimization.h"

// Project headers
#include "Bias.h"
#include "Bins.h"
#include "BootstrapSubsampler.h"
#include "Constants.h"
#include "DataSummary.h"
#include "FreeEnergyDistribution.hpp"
#include "FileSystem.h"
#include "GptlWrappers.h"
#include "InputParser.h"
#include "OpenMP.h"
#include "OrderParameter.h"
#include "OrderParameterRegistry.h"
#include "PointEstimator.h"
#include "Random.h"
#include "Simulation.h"
#include "Wham.h"

#include "Estimator_F_x.hpp"
#include "Estimator_F_x_y.hpp"


// Implements the Unbinned Weighted Histogram Analysis Method (UWHAM)
// - Equations are solved using log-likelihood maximation approach
//   - See Tan, Gallicchio, Lapelosa, & Levy (J. Chem. Phys. 2012)
// - See also:
//   - Zhu & Hummer (J. Comp. Chem. 2011)
//   - Souaille & Roux (Comp. Phys. Comm. 2001)
//
// NOTES:
// - "x" and "y" are generic names for the order parameters in question
// - All energies are in units of k_B*T unless otherwise noted
//
// TODO:
// - Allow u_bias-values as input
// - Check for internal consistency: 
//   - all time series for all order parameters
//   - biasing parameters
class WhamDriver
{
 public:
	// Types
	using ColumnVector = dlib::matrix<double,0,1>;

	WhamDriver(
		const std::string& options_file
	);

	// TODO: Move to standalone member variables?
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
	DataSummary data_summary_;

	// TODO:
	OrderParameterRegistry op_registry_;

	std::string biases_log_file_;

	// Organizes data by simulaton (i.e. by ensemble) 
	std::vector<Simulation> simulations_;

	WhamDriver::WhamOptions wham_options_;

	std::vector<double> f_bias_guess_;
	std::vector<double> f_bias_opt_, error_f_bias_opt_;

	// Organizes time series data and distributions for each OP
	// - TODO: UPDATE
	std::vector<OrderParameter> order_parameters_;

	// Objects that read in and evaluate the bias used in each simulation
	std::vector<Bias> biases_;  // [num_simulations x 1]

	// Error estimation
	enum class ErrorMethod { None, Bootstrap };
	ErrorMethod error_method_ = ErrorMethod::None;
	int num_bootstrap_samples_ = 100;

	void calculateBootstrapErrors(const Wham& wham);


	//----- Output -----//

	bool be_verbose_ = false;  // extra feedback
	bool be_quiet_   = false;  // minimal/no feedback

	// Output F(x)
	std::vector<FreeEnergyDistribution> output_f_x_;

	// F(x) to compute
	std::vector<Estimator_F_x> est_f_x_;

	// F(x,y) to compute
	std::vector<Estimator_F_x_y> est_f_x_y_;

	// Determine which outputs to compute and print
	// - ex. F(x) and F(x,y) for different OPs x and/or y
	void parseOutputs(const ParameterPack& input_pack);

	void printDistributions(const Wham& wham) const;

	void printWhamDistribution(
		const OrderParameter& x,
		const FreeEnergyDistribution& f,
		std::string file_name = ""
	) const;


	//----- GPTL -----//

	using Timer = GPTL::Timer;

	mutable Timer setup_timer_      = Timer("setup");
	mutable Timer driver_timer_     = Timer("driver");
	mutable Timer solve_wham_timer_ = Timer("solve_wham");

	mutable Timer bootstrap_timer_ = Timer("bootstrap");
	mutable Timer subsample_timer_ = Timer("subsample");

	mutable Timer print_output_timer_ = Timer("print_output");
};

#endif // ifndef WHAM_DRIVER_H
