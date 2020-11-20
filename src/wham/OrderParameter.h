// AUTHOR: Sean M. Marks (https://github.com/seanmarks) 

#ifndef ORDER_PARAMETER_H
#define ORDER_PARAMETER_H

#include <string>

#include "Bins.h"
#include "InputParser.h"

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


 private:
	std::string name_;

	Bins bins_;
};

#endif /* ORDER_PARAMETER_H */
