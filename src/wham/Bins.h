
#ifndef BINS_H
#define BINS_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "InputParser.h"


// A set of histogram bins
// - Currently only supports bins of equal size
class Bins
{
 public:
	// For bins that span the range [x_i, x_i + delta_x], this determines
	// the nominal x-value of the bin, x_b
	// - Left:   x_b = x_i
	// - Center: x_b = x_i + 0.5*delta_x
	// - Right:  x_b = x_i + delta_x
	enum class BinStyle { Left, Center, Right };

	Bins();

	Bins(
		const double min, 
		const double max, 
		const int num_bins,
		const BinStyle& bin_style = BinStyle::Left
	);

	Bins(const ParameterPack& input_pack);

	int getNumBins() const {
		return num_bins_;
	}

	double getBinSize() const {
		return bin_size_;
	}

	const std::vector<double>& getBins() const {
		return bins_;
	}

	void setBins(
		const double min, 
		const double max, 
		const int num_bins,
		const BinStyle& bin_style
	);

	void setBins(const ParameterPack& input_pack);

	// Returns the index of the appropriate bin for x (>=0),
	// else returns -1
	int findBin(const double x) const;


	//----- Operators ----//

	const double& operator[](const int b) const {
		return bins_[b];
	}


 private:
	double min_, max_, bin_size_;
	int num_bins_;  // store explicitly for convenience
	BinStyle bin_style_;
	std::vector<double> bins_;

	BinStyle parseBinStyle(const std::string& bin_style_token) const;
};

#endif /* BINS_H */
