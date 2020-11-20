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
	op_registry_ptr_(&op_registry)
{
	// Read time series
	const auto& time_series_files = op_registry_ptr_->getSimulationFiles(data_set_label);
	const int num_ops = time_series_files.size();
	time_series_.reserve(num_ops);
	for ( int p=0; p<num_ops; ++p ) {
		int data_col = op_registry_ptr_->getTimeSeriesDataCol(p);

		time_series_.emplace_back( time_series_files[p], data_col, t_min_, t_max_, use_floored_times );
	}
	checkTimeSeries();
}


void Simulation::setShuffledFromOther(const Simulation& other, const std::vector<int>& indices)
{
	FANCY_ASSERT(this != &other, "unsupported usage");

	*this = other;
	
	const int num_time_series = other.time_series_.size();
	for ( int p=0; p<num_time_series; ++p ) {
		time_series_[p].setShuffledFromOther( other.time_series_[p], indices );
	}
}


void Simulation::checkTimeSeries() const
{
	FANCY_ASSERT(op_registry_ptr_ != nullptr, "order parameter registry is missing");

	// Check number present
	const int num_time_series = time_series_.size();
	const int num_ops         = op_registry_ptr_->getNumRegistered();
	FANCY_ASSERT( num_time_series == num_ops,
		"Error setting up Simulation with data set label " << data_set_label_ << "\n"
		<< "  Mismatch between number of time series files parsed (" << num_time_series
		<< ") and number of order parameters (" << num_ops << "\n" );

	if ( num_time_series == 0 ) {
		return;  // nothing to do
	}

	const auto& time_series_ref = time_series_.front();
	for ( int p=1; p<num_time_series; ++p ) {
		// Check sampled times
		const auto& time_series_p = time_series_[p];
		if ( time_series_p.get_times() != time_series_ref.get_times() ) {
			std::stringstream err_ss;
			err_ss << "Error setting up Simulation with data set label " << data_set_label_ << "\n"
			       << "  Mismatch in times sampled for the following order parameters\n";
			std::vector<int> op_indices = {{ 0, p }};
			for ( auto j : op_indices ) {
				const auto& time_series_j = time_series_[j];

				err_ss << "    " << op_registry_ptr_->indexToName(j) << ": " << time_series_j.size() << " points from file "
				                 << time_series_j.getFile() << "\n";
			}
			throw std::runtime_error( err_ss.str() );
		}
	}
}
