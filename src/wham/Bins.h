// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

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
	
	// Typedefs
	using VectorReal     = std::vector<double>;
	using const_iterator = typename VectorReal::const_iterator;


	//----- Setup -----//

	Bins();

	Bins(
		const double min, 
		const double max, 
		const int num_bins,
		const BinStyle& bin_style = BinStyle::Left
	);

	Bins(const ParameterPack& input_pack);


	//----- Settings -----//

	int getNumBins() const noexcept {
		return num_bins_;
	}

	// Returns the number of bins
	std::size_t size() const noexcept {
		return num_bins_;
	}

	// Returns the size of each bin
	double getBinSize() const noexcept {
		return bin_size_;
	}

	// Returns the range of values spanned by the bins
	double getSpan() const noexcept {
		return max_ - min_;
	}


	//----- Bin values -----//

	// Returns the index of the appropriate bin for x (>=0),
	// else returns -1
	int findBin(const double x) const;


	//----- Set -----//

	void setBins(
		const double min, 
		const double max, 
		const int num_bins,
		const BinStyle& bin_style
	);

	// Sets the bins using the ParameterPack contents
	void setBins(const ParameterPack& input_pack);


	//----- Access values ----//

	// Returns an iterator to the first bin value
	const_iterator begin() const noexcept {
		return bins_.begin();
	}

	// Returns an iterator to the end of the bins
	const_iterator end() const noexcept {
		return bins_.end();
	}

	// Returns the location of the 'b'th bin 
	const double& operator[](const int b) const {
		return bins_[b];
	}

	// TODO: remove?
	const VectorReal& getBins() const noexcept {
		return bins_;
	}


 private:
	double min_, max_, bin_size_;
	int num_bins_;  // store explicitly for convenience
	BinStyle bin_style_;

	VectorReal bins_;

	BinStyle parseBinStyle(const std::string& bin_style_token) const;
};

#endif /* BINS_H */
