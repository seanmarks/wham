// Simulation

#pragma once
#ifndef SIMULATION_H
#define SIMULATION_H

#include <string>

#include "Constants.h"
#include "DataSummary.h"
#include "OrderParameterRegistry.h"

class Simulation 
{
 public:
	// TODO transition to private variables
	double f_bias;              // Free energy of adding bias (in kBT): f_bias = -ln(Q_i/Q_0)
	double f_bias_guess = 0.0;  // First estimate of biasing free energy

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

 private:
	std::string data_set_label_;
	double      t_min_, t_max_;  // Sampling range [ps]
	double      temperature_;
	double      kBT_, beta_;     // k_B*T [kJ/mol] and beta = 1/(k_B*T) [mol/kJ]

	bool use_floored_times_ = false;  // default: leave input unmodified
};

#endif /* ifndef SIMULATION_H */
