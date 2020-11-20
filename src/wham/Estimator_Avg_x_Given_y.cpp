// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "Estimator_Avg_x_Given_y.hpp"

#include <numeric>


Estimator_Avg_x_Given_y::Estimator_Avg_x_Given_y(
  const OrderParameter& x, const OrderParameter& y
):
  x_(x), y_(y), est_f_y_(y)
{}


void Estimator_Avg_x_Given_y::calculateUsingStoredWeights(
  const Wham& wham
)
{  
  const auto& simulations = wham.getSimulationData();
  const auto& data_ranges = wham.getSimulationDataRanges();
  const int num_simulations = simulations.size();

	// Sort samples by bin
  const auto& bins_y = y_.getBins();
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
		const auto& x_j = simulations[j].getTimeSeriesForOrderParameter(x_.getName());
		const auto& y_j = simulations[j].getTimeSeriesForOrderParameter(y_.getName());

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
  avg_x_given_y_.resize(num_bins_y);
  sample_counts_.resize(num_bins_y);
  for ( int b=0; b<num_bins_y; ++b ) {
    sample_counts_[b] = binned_weights_[b].size();
    if ( p_y[b] > 0.0 ) {
      avg_x_given_y_[b] = std::accumulate( binned_weights_[b].begin(), binned_weights_[b].end(), 0.0 );
      avg_x_given_y_[b] /= p_y[b] * bins_y.getBinSize();
    }
    else {
      avg_x_given_y_[b] = -1.0;
    }
  }
}



void Estimator_Avg_x_Given_y::saveResults(std::string file_name) const
{
  if ( file_name.empty() ) {
    file_name = "avg_" + x_.getName() + "_given_" + y_.getName() + ".out";
  }
  std::ofstream ofs(file_name);
  FANCY_ASSERT( ! ofs.fail(), "unable to open file: " << file_name );

  // Header
  const std::string spacer("  ");
  ofs << "# " << y_.getName() << spacer << "<" << x_.getName() << "|" << y_.getName() << ">\n";

  // Data
  const auto& bins_y = y_.getBins();
  const int num_bins = bins_y.getNumBins();
  for ( int b=0; b<num_bins; ++b ) {
    ofs << bins_y[b] << spacer << avg_x_given_y_[b] << "\n";
  }

  ofs.close();
}