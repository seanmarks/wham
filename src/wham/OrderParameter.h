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

// Organizes settings for a single order parameter, especially
// its name and how its values should be binned
class OrderParameter
{
 public:
	OrderParameter(
		const std::string& name,
		const ParameterPack& input_pack
	);

	// Returns a string representing the name of the OP
	const std::string& getName() const {
		return name_;
	}

	const Bins& getBins() const {
		return bins_;
	}

	void setWhamDistribution(const FreeEnergyDistribution& wham_distribution) {
		wham_distribution_ = wham_distribution;
	}


	//----- Printing Output -----//

	void printWhamResults(std::string file_name, const double temp) const;


 private:
	std::string name_;

	Bins bins_;

	// WHAM results
	FreeEnergyDistribution wham_distribution_;
};

#endif /* ORDER_PARAMETER_H */
