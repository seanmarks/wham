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
#include <cassert>
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
	using TimeSeriesPtr = Simulation::TimeSeriesPtr;

	// FIXME
	friend class Wham;
	friend class WhamDriver;
	
	OrderParameter(
		const std::string& name,
		const ParameterPack& input_pack,
		std::vector<Simulation>& simulations
	);

	// Get functions
	const std::string& get_name() const { return name_; }
	const TimeSeries& get_time_series(const int j) const {
#ifdef DEBUG
		assert( time_series_ptrs_[j] != nullptr );
#endif // ifdef DEBUG
		return *(time_series_ptrs_[j]);
	}
	const Bins& get_bins() const { return bins_; }

	// Set functions
	void set_unbiased_distributions(const std::vector<Distribution>& unbiased_distributions) {
		unbiased_distributions_ = unbiased_distributions;
	}
	void set_shifted_distributions(const std::vector<Distribution>& shifted_distributions) {
		shifted_distributions_ = shifted_distributions;
	}
	void set_rebiased_distributions(const std::vector<Distribution>& rebiased_distributions) {
		rebiased_distributions_ = rebiased_distributions;

		// Entropy between biased and rebiased distributions
		int num_simulations = simulations_.size();
		info_entropy_.resize(num_simulations);
		for ( int j=0; j<num_simulations; ++j ) {
			info_entropy_[j] = Distribution::computeInformationEntropy(
				biased_distributions_[j], rebiased_distributions_[j]
			);
		}
	}
	void set_wham_distribution(const Distribution& wham_distribution) {
		wham_distribution_ = wham_distribution;
	}


	//----- Printing Output -----//

	// TODO rename
	void printRawDistributions() const;

	void printRebiasedDistributions(std::string file_name = "") const;

	void printWhamResults(std::string file_name = "") const;

	void printStats(std::string file_name = "") const;

	// Prints a series of distributions, F_i(x), side-by-side
	void printDistributions(
		const std::vector<Distribution>& distributions,
		const std::string& file_name, 
		const std::string& header,
		const bool shift_to_zero = true
	) const;

 private:
	std::string name_;

	std::vector<Simulation>& simulations_;
	Bins bins_;

	// Time series data from each simulation
	std::vector<TimeSeriesPtr> time_series_ptrs_;

	// TODO: Move to driver?
	std::vector<Distribution> biased_distributions_;
	std::vector<Distribution> unbiased_distributions_;
	std::vector<Distribution> rebiased_distributions_;
	std::vector<Distribution> shifted_distributions_;

	// WHAM results
	Distribution wham_distribution_;
	std::vector<double> info_entropy_;  // entropy between f_biased and f_rebiased

	// Number of samples in each bin, across all simulations
	std::vector<int> global_sample_counts_;

	/*
	// Checks a set of OrderParameters for consistency:
	static void checkForConsistency(const std::vector<OrderParameter>& ops);
	*/
};

#endif /* ORDER_PARAMETER_H */
