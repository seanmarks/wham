#include "TimeSeries.h"

TimeSeries::TimeSeries(
		const std::string& file, const int col, 
		const double t0, const double tf, const bool use_floored_input_times
):
	file_(file), col_(col), t0_(t0), tf_(tf),
	use_floored_input_times_(use_floored_input_times)
{
	// Working variables
	std::string line, token;
	double time;
	std::vector<double> values(2);

	std::ifstream ifs( file_ );
	if ( not ifs.is_open() ) {
		throw std::runtime_error("Failed to open time series file \'" + file + "\'");
	}

	while ( getline(ifs, line) ) {
		std::stringstream ss(line);
		ss >> token;

		// Ignore comments and blank lines
		if ( line.empty() or token[0] == '#' ) {
			continue;
		}

		// Get the time
		time = std::stod(token);
		if ( use_floored_input_times_ ) { time = round(time); }

		// Check for production phase
		if ( time >= t0_ and time <= tf_ ) {
			// Parse line
			values.assign(1, time);
			while ( ss >> token ) {
				values.push_back( std::stod(token) );
			}

			if ( static_cast<int>(values.size()) < col_ + 1 ) {
				throw std::runtime_error("col exceeds the number of columns in data file " + file);
			}

			// Store data
			times_.push_back( time );
			data_.push_back( values[col] );
		}
	}
	ifs.close();
}



double TimeSeries::average(const std::vector<double>& x) const
{
	int num_values = x.size();
	if ( num_values == 0 ) {
		throw std::runtime_error("error in average() - no data provided");
	}

	double avg = 0.0;
	for ( int i=0; i<num_values; ++i ) {
		avg += x[i];
	}
	avg /= static_cast<double>(num_values);

	return avg;
}


double TimeSeries::variance(const std::vector<double>& x) const
{
	int num_values = x.size();
	if ( num_values == 0.0 ) {
		throw std::runtime_error("error in variance() - no data provided");
	}

	double var = 0.0;
	double dx, avg_x = average(x);
	for ( int i=0; i<num_values; ++i ) {
		dx = x[i] - avg_x;
		var += dx*dx;
	}
	var /= static_cast<double>(num_values);

	return var;
}
