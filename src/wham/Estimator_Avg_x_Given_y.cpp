// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "Estimator_Avg_x_Given_y.hpp"

#include <numeric>


Estimator_Avg_x_Given_y::Estimator_Avg_x_Given_y(
  const OrderParameter& x, const OrderParameter& y
):
  avg_x_given_y_(x,y), est_f_y_(y)
{}



void Estimator_Avg_x_Given_y::calculateUsingStoredWeights(
  const Wham& wham
)
{  
  const auto& simulations = wham.getSimulationData();
  const auto& data_ranges = wham.getSimulationDataRanges();
  const int num_simulations = simulations.size();

	// Sort samples by bin
  const auto& bins_y = get_y().getBins();
  const int   num_bins_y = bins_y.getNumBins();
  const int   num_samples_total = wham.getNumSamplesTotal(); 
  const int   num_to_reserve = 2*(num_samples_total/num_bins_y);
  binned_weights_.resize(num_bins_y);
  for ( int b=0; b<num_bins_y; ++b ) {
    binned_weights_[b].resize(0);
    binned_weights_[b].reserve(num_to_reserve);
  }
  const auto& weights = this->getWeights();
	for ( int j=0; j<num_simulations; ++j ) {
		const auto& x_j = simulations[j].getTimeSeriesForOrderParameter(get_x().getName());
		const auto& y_j = simulations[j].getTimeSeriesForOrderParameter(get_y().getName());

		const int num_samples = y_j.size();
		for ( int i=0; i<num_samples; ++i ) {
			int b = bins_y.findBin( y_j[i] );
			if ( b >= 0 ) {
				int n = data_ranges[j].first + i;
        binned_weights_[b].push_back( x_j[i]*weights[n] );
			}
		}
	}

  // Need P(y) to condition the average
  est_f_y_.calculate(wham);
  const auto& p_y = est_f_y_.get_f_x().get_P_x();

  // Finish calculations
  std::vector<double> avg_x_given_y(num_bins_y);
  std::vector<int>    sample_counts(num_bins_y);
  for ( int b=0; b<num_bins_y; ++b ) {
    sample_counts[b] = binned_weights_[b].size();
    if ( p_y[b] > 0.0 && sample_counts[b] >= 2 ) {
      avg_x_given_y[b] = std::accumulate( binned_weights_[b].begin(), binned_weights_[b].end(), 0.0 );
      avg_x_given_y[b] /= p_y[b] * bins_y.getBinSize();
    }
    else {
      avg_x_given_y[b] = -1.0;
    }
  }

  avg_x_given_y_.set(avg_x_given_y, sample_counts);
}



void Estimator_Avg_x_Given_y::saveResults(std::string file_name) const
{
  if ( file_name.empty() ) {
    file_name = "avg_" + get_x().getName() + "_given_" + get_y().getName() + ".out";
  }
  avg_x_given_y_.save(file_name);
}