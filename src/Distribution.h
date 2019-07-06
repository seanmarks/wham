
#pragma once
#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "Bins.h"
#include "TimeSeries.h"

class Distribution {
 public:
	Bins                bins_x;
	std::vector<double> f_x, error_f_x;
	std::vector<double> p_x;
	std::vector<int>    sample_counts;

	Distribution();

	Distribution(const Bins& bins_x_in);

	Distribution(const Bins& bins_x_in, const TimeSeries& time_series_x);

	void print(
		const std::string& file, 
		const std::string& header
	) const;
};

#endif /* ifndef DISTRIBUTION_H */
