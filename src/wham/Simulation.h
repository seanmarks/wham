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
#include "utils.h"

class Simulation 
{
 public:
	using TimeSeriesPtr = std::shared_ptr<TimeSeries>;

	Simulation():
		time_series_ptrs_(0)
	{}

	Simulation(
		const std::string&            data_set_label,
		const double                  t_min,
		const double                  t_max,
		const double                  temperature,
		const bool                    use_floored_times,
		const OrderParameterRegistry& op_registry
	);

	const std::string& get_data_set_label() const { return data_set_label_; }
	double get_t_min()       const { return t_min_; }
	double get_t_max()       const { return t_max_; }
	double get_temperature() const { return temperature_; }

	int get_num_samples() const {
		if ( time_series_ptrs_.size() > 0 ) {
			return time_series_ptrs_[0]->size();
		}
		else {
			return 0;
		}
	}

	TimeSeriesPtr copy_time_series_ptr(const std::string& op_name) {
		FANCY_ASSERT(op_registry_ptr_ != nullptr, "order parameter registry is missing");
		int op_index = op_registry_ptr_->get_index(op_name);
		return time_series_ptrs_[op_index];
	}

	void setShuffledFromOther(const Simulation& other, const std::vector<int>& indices);

 private:
	std::string data_set_label_ = "";
	double      t_min_ = -1.0, t_max_ = -1.0;  // Sampling range [ps]
	double      temperature_ = -1.0;
	double      kBT_, beta_ = -1.0;            // k_B*T [kJ/mol] and beta = 1/(k_B*T) [mol/kJ]

	bool use_floored_times_ = false;  // default: leave input unmodified

	const OrderParameterRegistry* op_registry_ptr_ = nullptr;

	std::vector<TimeSeriesPtr> time_series_ptrs_ = {};

	// Check time series for consistency
	void checkTimeSeries() const;
};

#endif /* ifndef SIMULATION_H */
