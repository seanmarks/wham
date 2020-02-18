// Simulation

#pragma once
#ifndef SIMULATION_H
#define SIMULATION_H

#include <memory>
#include <string>

#include "Constants.h"
#include "DataSummary.h"
#include "OrderParameterRegistry.h"
#include "TimeSeries.h"

class Simulation 
{
 public:
	using TimeSeriesPtr = std::shared_ptr<TimeSeries>;

	Simulation(
		const std::string&            data_set_label,
		const double                  t_min,
		const double                  t_max,
		const double                  temperature,
		const bool                    use_floored_times,
		const OrderParameterRegistry& op_registry
	);

	const std::string& get_data_set_label() const { return data_set_label_; }
	double get_t_min() const { return t_min_; }
	double get_t_max() const { return t_max_; }

	int get_num_samples() const {
		if ( time_series_ptrs_.size() > 0 ) {
			return time_series_ptrs_[0]->size();
		}
		else {
			return 0;
		}
	}

	TimeSeriesPtr copy_time_series_ptr(const std::string& op_name) {
		int op_index = op_registry_.get_index(op_name);
		return time_series_ptrs_[op_index];
	}

 private:
	std::string data_set_label_;
	double      t_min_, t_max_;  // Sampling range [ps]
	double      temperature_;
	double      kBT_, beta_;     // k_B*T [kJ/mol] and beta = 1/(k_B*T) [mol/kJ]

	bool use_floored_times_ = false;  // default: leave input unmodified

	const OrderParameterRegistry& op_registry_;

	std::vector<TimeSeriesPtr> time_series_ptrs_;

	// Check time series for consistency
	void checkTimeSeries() const;
};

#endif /* ifndef SIMULATION_H */
