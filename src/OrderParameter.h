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
#include "TimeSeries.h"
//#include "Wham.h"

class OrderParameter
{
 public:
	friend class Wham;

	// Get functions
	const std::string& get_name() const { return name_; }
	const std::vector<TimeSeries>& get_time_series() const { return time_series_; }
	const Bins& get_bins() const { return bins_; }

 private:
	std::string name_;

	// Time series data from each simulation
	std::vector<TimeSeries> time_series_;

	Bins bins_;

	// - Size: [ num_simulations x num_bins ]
	std::vector<std::vector<double>> p_biased_,   f_biased_;   // F_biased(x)
	std::vector<std::vector<double>> p_rebiased_, f_rebiased_;
	std::vector<std::vector<double>> p_unbiased_, f_unbiased_;
	std::vector<std::vector<int>>    sample_counts_;

	// Number of samples in each bin, across all simulations
	std::vector<int> global_sample_counts_;
};

#endif /* ORDER_PARAMETER_H */
