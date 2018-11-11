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

class TimeSeries
{
 public:
	TimeSeries(
		const std::string& file, 
		const int col,
		// Time range to keep
		const double t0, 
		const double tf,
		// Whether to round the input time down to the nearest integer
		const bool use_floored_input_times
	);

	TimeSeries(
		const std::vector<double>& data
	): col_(-1), t0_(-1.0), tf_(-1.0), data_(data)
	{
		times_.assign( data_.size(), -1.0 );
	}

	// Returns the number of samples in the time series
	unsigned size() const { return data_.size(); }

	const std::string& get_file() const { return file_; }

	// Access underlying data
	const double& operator[](const int i) const;
	const std::vector<double>& get_data() const { return data_; }
	std::vector<double>& access_data() { return data_; }

	// Statistics of an arbitrary data set
	double average(const std::vector<double>& x) const;
	double variance(const std::vector<double>& x) const;

	// Statistics of stored time series
	double average() const { return average(data_); }
	double variance() const { return variance(data_); }

 private:
	// Where data is stored, and what/how to read it in
	std::string file_;
	int col_;
	double t0_, tf_;  // time range to keep
	bool use_floored_input_times_;

	// Underlying time series
	std::vector<double> times_;
	std::vector<double> data_;
};

inline
const double& TimeSeries::operator[](const int i) const {
	return data_[i];
}



#endif /* TIME_SERIES_H */
