// AUTHOR: Sean M. Marks (https://github.com/seanmarks) 

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
#include "FileSystem.h"
#include "FreeEnergyDistribution.hpp"
#include "InputParser.h"
#include "Simulation.h"
#include "TimeSeries.h"
#include "Assert.hpp"

// Organizes data and variables for a single order parameter 
// across multiple simulations
// - TODO: remove storage of time series data?
class OrderParameter
{
 public:
	OrderParameter(
		const std::string& name,
		const ParameterPack& input_pack,
		std::vector<Simulation>& simulations
	);

	// Returns a string representing the name of the OP
	const std::string& getName() const {
		return name_;
	}

	// FIXME: REMOVE
	const TimeSeries& getTimeSeries(const int j) const {
#ifdef DEBUG
		assert( time_series_ptrs_[j] != nullptr );
#endif // ifdef DEBUG
		return *(time_series_ptrs_[j]);
	}

	const Bins& getBins() const {
		return bins_;
	}

	// Sets the handle to the simulations to use
	void setSimulations(std::vector<Simulation>& simulations);

	// FIXME: cludgy to set this way

	void setUnbiasedDistributions(const std::vector<FreeEnergyDistribution>& unbiased_distributions) {
		unbiased_distributions_ = unbiased_distributions;
	}

	void setShiftedDistributions(const std::vector<FreeEnergyDistribution>& shifted_distributions) {
		shifted_distributions_ = shifted_distributions;
	}

	void setRebiasedDistributions(const std::vector<FreeEnergyDistribution>& rebiased_distributions) {
		rebiased_distributions_ = rebiased_distributions;

		// Entropy between biased and rebiased distributions
		int num_simulations = simulation_ptrs_.size();
		info_entropy_.resize(num_simulations);
		for ( int j=0; j<num_simulations; ++j ) {
			info_entropy_[j] = FreeEnergyDistribution::computeInformationEntropy(
				rebiased_distributions_[j], biased_distributions_[j]
			);
		}
	}

	void setWhamDistribution(const FreeEnergyDistribution& wham_distribution) {
		wham_distribution_ = wham_distribution;
	}


	//----- Printing Output -----//

	// TODO: rename
	void printRawDistributions() const;

	void printRebiasedDistributions(std::string file_name = "") const;

	void printWhamResults(std::string file_name = "") const;

	void printStats(std::string file_name = "") const;

	// Prints a series of distributions, F_i(x), side-by-side
	void printDistributions(
		const std::vector<FreeEnergyDistribution>& distributions,
		const std::string& file_name, 
		const std::string& header,
		const bool shift_to_zero = true
	) const;

 private:
	std::string name_;

	std::vector<Simulation*> simulation_ptrs_;

	Bins bins_;

	// Time series data from each simulation
	std::vector<const TimeSeries*> time_series_ptrs_;

	std::vector<FreeEnergyDistribution> biased_distributions_;
	std::vector<FreeEnergyDistribution> unbiased_distributions_;
	std::vector<FreeEnergyDistribution> shifted_distributions_;   // unbiased and shifted
	std::vector<FreeEnergyDistribution> rebiased_distributions_;  // rebias consensus distribution

	// WHAM results
	FreeEnergyDistribution wham_distribution_;
	std::vector<double> info_entropy_;  // entropy between f_biased and f_rebiased

	// Number of samples in each bin, across all simulations
	std::vector<int> global_sample_counts_;
};

#endif /* ORDER_PARAMETER_H */
