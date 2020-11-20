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

	// FIXME:
	void setInfoEntropy(const std::vector<double>& eta) {
		info_entropy_ = eta;
	}

	void setWhamDistribution(const FreeEnergyDistribution& wham_distribution) {
		wham_distribution_ = wham_distribution;
	}


	//----- Printing Output -----//

	void printWhamResults(std::string file_name = "") const;

	void printStats(std::string file_name = "") const;


 private:
	std::string name_;

	std::vector<Simulation*> simulation_ptrs_;

	Bins bins_;

	// Time series data from each simulation
	std::vector<const TimeSeries*> time_series_ptrs_;

	std::vector<FreeEnergyDistribution> rebiased_distributions_;  // rebias consensus distribution

	// WHAM results
	FreeEnergyDistribution wham_distribution_;
	std::vector<double> info_entropy_;  // entropy between f_biased and f_rebiased
};

#endif /* ORDER_PARAMETER_H */
