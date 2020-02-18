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
	use_floored_times_(use_floored_times),
	op_registry_(op_registry)
{
	// Read time series
	const auto& time_series_files = op_registry_.get_simulation_files(data_set_label);
	const int num_ops = time_series_files.size();
	time_series_ptrs_.reserve(num_ops);
	for ( int p=0; p<num_ops; ++p ) {
		int data_col = op_registry_.get_time_series_data_col(p);

		time_series_ptrs_.push_back( 
			std::make_shared<TimeSeries>(time_series_files[p], data_col, t_min_, t_max_, use_floored_times)
		);
	}
	checkTimeSeries();
}


void Simulation::checkTimeSeries() const
{
	// Check number present
	int num_time_series = time_series_ptrs_.size();
	int num_ops         = op_registry_.getNumberOfOrderParameters();
	if ( num_time_series != num_ops ) {
		std::stringstream err_ss;
		err_ss << "Error setting up Simulation with data set label " << data_set_label_ << "\n"
		       << "  Mismatch between number of time series files parsed (" << num_time_series
		         << ") and number of order parameters (" << num_ops << "\n";
		throw std::runtime_error( err_ss.str() );
	}

	if ( num_time_series == 0 ) {
		return;  // nothing to do
	}

	const auto& time_series_ref = *(time_series_ptrs_[0]);
	for ( int p=1; p<num_time_series; ++p ) {
		// Check sampled times
		const auto& time_series_p = *(time_series_ptrs_[p]);
		if ( time_series_p.get_times() != time_series_ref.get_times() ) {
			std::stringstream err_ss;
			err_ss << "Error setting up Simulation with data set label " << data_set_label_ << "\n"
			       << "  Mismatch in times sampled for the following order parameters\n";
			std::vector<int> op_indices = {{ 0, p }};
			for ( auto j : op_indices ) {
				const auto& time_series_j = *(time_series_ptrs_[j]);

				err_ss << "    " << op_registry_.get_name(j) << ": " << time_series_j.size() << " points from file "
				                 << time_series_j.get_file() << "\n";
			}
			throw std::runtime_error( err_ss.str() );
		}
	}
}
