// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#pragma once
#ifndef SIMULATION_H
#define SIMULATION_H

#include <memory>
#include <string>

#include "Constants.h"
#include "DataSummary.h"
#include "OrderParameterRegistry.h"
#include "TimeSeries.h"
#include "Assert.hpp"

// Stores time series data and settings corresponding to a single simulation (ensemble)
class Simulation 
{
 public:
	Simulation() = default;

	Simulation(
		const std::string&            data_set_label,
		const double                  t_min,
		const double                  t_max,
		const double                  temperature,
		const bool                    use_floored_times,
		const OrderParameterRegistry& op_registry
	);

	const std::string& getDataSetLabel() const {
		return data_set_label_; 
	}

	double get_t_min() const {
		return t_min_;
	}

	double get_t_max() const {
		return t_max_;
	}

	double getTemperature() const {
		return temperature_;
	}

	int getNumSamples() const {
		if ( time_series_.size() > 0 ) {
			return time_series_.front().size();
		}
		else {
			return 0;
		}
	}

	// Access the time series for the given OP
	// - Throws if the OP does is not registered
	const TimeSeries& getTimeSeriesForOrderParameter(const std::string& op_name) const {
		FANCY_ASSERT(op_registry_ptr_ != nullptr, "order parameter registry is missing");
		const int op_index = op_registry_ptr_->nameToIndex(op_name);
		return time_series_[op_index];
	}

	// Returns a raw (non-owning) handle to the time series of the given OP
	const TimeSeries* copyTimeSeriesPtr(const std::string& op_name) const {
		const auto& time_series = getTimeSeriesForOrderParameter(op_name);
		return &time_series;
	}

	void setShuffledFromOther(const Simulation& other, const std::vector<int>& indices);

 private:
	std::string data_set_label_ = "";
	double      t_min_ = -1.0, t_max_ = -1.0;  // Sampling range [ps]
	double      temperature_ = -1.0;
	double      kBT_, beta_ = -1.0;            // k_B*T [kJ/mol] and beta = 1/(k_B*T) [mol/kJ]

	bool use_floored_times_ = false;  // default: leave input unmodified

	const OrderParameterRegistry* op_registry_ptr_ = nullptr;

	std::vector<TimeSeries> time_series_;

	// Check time series for consistency
	void checkTimeSeries() const;
};

#endif /* ifndef SIMULATION_H */
