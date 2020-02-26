// TimeSeries
//
// ABOUT: Reads in and stores a single time series

#ifndef TIME_SERIES_H
#define TIME_SERIES_H

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

#include "utils.h"

class TimeSeries
{
 public:
	// Construct using data stored in a file
	TimeSeries(
		const std::string& file, 
		const int col,
		// Time range to keep
		const double t0, 
		const double tf,
		// Whether to round the input time down to the nearest integer
		const bool use_floored_input_times
	);

	// Construct using a set of plain values
	TimeSeries(const std::vector<double>& data):
		times_(data.size(), -1.0), data_(data)
	{}

	TimeSeries():
		times_(0), data_(0)
	{}

	// Extract data from another TimeSeries according to 'indices'
	void setShuffledFromOther(const TimeSeries& other, const std::vector<int>& indices);

	// Returns the number of samples in the time series
	unsigned size() const { return data_.size(); }

	const std::string& get_file() const { return file_; }

	// Access underlying data
	const double& operator[](const int i) const;
	const std::vector<double>& get_data() const { return data_; }
	std::vector<double>& access_data() { return data_; }

	const std::vector<double>& get_times() const { return times_; }

	// Statistics of an arbitrary data set
	double average(const std::vector<double>& x) const;
	double variance(const std::vector<double>& x) const;

	// Statistics of stored time series
	double average() const { return average(data_); }
	double variance() const { return variance(data_); }

 private:
	// Where data is stored, and what/how to read it in
	std::string file_ = "";
	int col_ = -1;
	double t0_ = -1.0, tf_ = -1.0;  // time range to keep

	bool use_floored_input_times_ = false;
	bool is_shuffled_ = false;

	// Underlying time series
	std::vector<double> times_;
	std::vector<double> data_;
};

inline
const double& TimeSeries::operator[](const int i) const {
	return data_[i];
}



#endif /* TIME_SERIES_H */
