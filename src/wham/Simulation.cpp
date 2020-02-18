#include "Simulation.h"

Simulation::Simulation(
	const std::string& data_set_label, const double t_min, const double t_max,
	const double temperature, const bool use_floored_times,
	const OrderParameterRegistry& op_registry
):
	data_set_label_(data_set_label),
	t_min_(t_min), t_max_(t_max),
	temperature_(temperature),
	kBT_(Constants::k_B * temperature_), beta_(1.0/kBT_),
	use_floored_times_(use_floored_times)
{
	// TODO Read time series
}	
