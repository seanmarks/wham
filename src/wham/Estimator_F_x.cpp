// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "Estimator_F_x.hpp"


Estimator_F_x::Estimator_F_x(const OrderParameter& x):
  x_(x)
{
}


void Estimator_F_x::calculateUsingStoredWeights(
  const Wham& wham
)
{
  const auto& simulations = wham.getSimulationData();
  const auto& data_ranges = wham.getSimulationDataRanges();
  const int num_simulations = simulations.size();

	// Sort samples by bin
  const auto& bins_x = x_.getBins();
  const int   num_bins_x = bins_x.getNumBins();
  const int   num_samples_total = wham.getNumSamplesTotal(); 
  const int   num_to_reserve = 2*(num_samples_total/num_bins_x);
  binned_weights_.resize(num_bins_x);
  for ( int b=0; b<num_bins_x; ++b ) {
    binned_weights_[b].resize(0);
    binned_weights_[b].reserve(num_to_reserve);
  }
  const auto& weights = this->getWeights();
	for ( int j=0; j<num_simulations; ++j ) {
		const auto& x_j = simulations[j].getTimeSeriesForOrderParameter(x_.getName());
		const int num_samples = x_j.size();
		for ( int i=0; i<num_samples; ++i ) {
			int b = bins_x.findBin( x_j[i] );
			if ( b >= 0 ) {
				int n = data_ranges[j].first + i;
        binned_weights_[b].push_back( weights[n] );
			}
		}
	}

  // Calculate distributions using binned weights
  std::vector<double> p_x_tmp(num_bins_x);
  std::vector<double> f_x_tmp(num_bins_x);
  std::vector<int>    samples(num_bins_x);
	const double bin_size_x = bins_x.getBinSize();
  const double fac = 1.0/bin_size_x;  // normalize by bin size
  for ( int b=0; b<num_bins_x; ++b ) {
    double sum = std::accumulate( binned_weights_[b].begin(), binned_weights_[b].end(), 0.0 );
    p_x_tmp[b] = fac*sum;
    f_x_tmp[b] = -std::log(p_x_tmp[b]);
    samples[b] = binned_weights_[b].size();
  }

  f_x_ = FreeEnergyDistribution(bins_x, f_x_tmp, p_x_tmp, samples);
}