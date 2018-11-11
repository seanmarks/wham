// OrderParameter.h
// 
// ABOUT: Organizes data and variables for a single order parameter 
//        across multiple simulations
// - Wham (friend class) sets much of its internal state
//
// TODO better to convert to struct to make it clear it's really POD?

#ifndef ORDER_PARAMETER_H
#define ORDER_PARAMETER_H

// Standard headers
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Project headers
#include "Bins.h"
#include "FileSystem.h"
#include "InputParser.h"
#include "TimeSeries.h"
//#include "Wham.h"

class OrderParameter
{
 public:
	friend class Wham;
	
	using Range = std::array<double,2>;
	OrderParameter(
		const ParameterPack& input_pack,
		const std::vector<Range>& production_phases,
		const bool use_floored_times  // whether to floor input time values
	);

	// Get functions
	const std::string& get_name() const { return name_; }
	const std::vector<TimeSeries>& get_time_series() const { return time_series_; }
	const Bins& get_bins() const { return bins_; }

 private:
	std::string name_;

	// File containing the list of time series data files at column file_col_
	std::string time_series_list_;
	int file_col_;

	// Time series data from each simulation
	std::vector<TimeSeries> time_series_;
	int data_col_;

	Bins bins_;

	// - Size: [ num_simulations x num_bins ]
	std::vector<std::vector<double>> p_biased_,   f_biased_;    // F_biased(x)
	std::vector<std::vector<double>> p_unbiased_, f_unbiased_;  // using only data from 1 time series
	std::vector<std::vector<double>> p_rebiased_, f_rebiased_;  // Rebias F_WHAM(x)
	std::vector<std::vector<int>>    sample_counts_;

	// WHAM results
	std::vector<double> f_x_wham_, p_x_wham_, error_f_x_wham_;  // TODO errors
	std::vector<double> info_entropy_;  // entropy between f_biased and f_rebiased

	// Number of samples in each bin, across all simulations
	std::vector<int> global_sample_counts_;
};

#endif /* ORDER_PARAMETER_H */
