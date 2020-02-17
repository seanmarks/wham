// OrderParameter.h
// 
// ABOUT: Organizes data and variables for a single order parameter 
//        across multiple simulations
// - Wham (friend class) sets much of its internal state
//
// TODO better to convert to struct to make it clear it's really POD?
// - That's how friend class Wham sees it

#ifndef ORDER_PARAMETER_H
#define ORDER_PARAMETER_H

// Standard headers
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Project headers
#include "Bins.h"
#include "Distribution.h"
#include "FileSystem.h"
#include "InputParser.h"
#include "Simulation.h"
#include "TimeSeries.h"

class OrderParameter
{
 public:
	// FIXME
	friend class Wham;
	friend class WhamDriver;
	
	OrderParameter(
		const ParameterPack& input_pack,
		const std::vector<Simulation>& simulations,
		const bool use_floored_times  // whether to floor input time values
	);

	// Get functions
	const std::string& get_name() const { return name_; }
	const std::vector<TimeSeries>& get_time_series() const { return time_series_; }
	const Bins& get_bins() const { return bins_; }

	void printRawDistributions() const;

	// Prints a series of distributions, F_i(x), side-by-side
	void printDistributions(
		const std::vector<Distribution>& distributions,
		const std::string& file_name, 
		const std::string& header,
		const bool shift_to_zero = true
	) const;

 private:
	std::string name_;

	// File containing the list of time series data files at column file_col_
	std::string time_series_list_;
	int file_col_;

	const std::vector<Simulation>& simulations_;
	Bins bins_;

	// Time series data from each simulation
	std::vector<TimeSeries> time_series_;
	int data_col_;

	std::vector<Distribution> biased_distributions_;
	std::vector<Distribution> unbiased_distributions_;
	std::vector<Distribution> rebiased_distributions_;

	// WHAM results
	Distribution wham_distribution_;
	std::vector<double> info_entropy_;  // entropy between f_biased and f_rebiased

	// Number of samples in each bin, across all simulations
	std::vector<int> global_sample_counts_;

	// Checks the list of order parameters for consistency:
	// - Same number of time series
	// - Lengths of time series match
	// - Stored times match
	static void checkForConsistency(const std::vector<OrderParameter>& ops);
};

#endif /* ORDER_PARAMETER_H */
