// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "Estimator_F_x.hpp"


Estimator_F_x::Estimator_F_x(const OrderParameter& x):
  x_(x)
{
}


void Estimator_F_x::calculate(
  const Wham& wham,
	const std::vector<double>& u_bias_as_r, const double f_bias_r
)
{
  const int num_samples_total = u_bias_as_r.size();
  const int num_expected      = wham.getNumSamplesTotal();
  FANCY_ASSERT( num_samples_total == num_expected,
                "got " << num_samples_total << ", expected " << num_expected );

  // Compute weights
  // - TODO: Move to WHAM class?
  const auto& log_dhat = wham.get_log_dhat_opt();
  weights_.resize(num_samples_total);
  for ( int n=0; n<num_samples_total; ++n ) {
    weights_[n] = std::exp( f_bias_r - u_bias_as_r[n] - log_dhat[n] );
  }

  const auto& simulations = wham.getSimulationData();
  const auto& data_ranges = wham.getSimulationDataRanges();
  const int num_simulations = simulations.size();

	// Sort samples by bin
  const auto& bins_x = x_.getBins();
  const int   num_bins_x = bins_x.get_num_bins();
  const int   num_to_reserve = 2*(num_samples_total/num_bins_x);
  binned_weights_.resize(num_bins_x);
  for ( int b=0; b<num_bins_x; ++b ) {
    binned_weights_[b].resize(0);
    binned_weights_[b].reserve(num_to_reserve);
  }
	for ( int j=0; j<num_simulations; ++j ) {
		const auto& x_j = simulations[j].getTimeSeriesForOrderParameter(x_.getName());
		const int num_samples = x_j.size();
		for ( int i=0; i<num_samples; ++i ) {
			int b = bins_x.find_bin( x_j[i] );
			if ( b >= 0 ) {
				int n = data_ranges[j].first + i;
        binned_weights_[b].push_back( weights_[n] );
			}
		}
	}

  std::vector<double> p_x_tmp(num_bins_x);
  std::vector<double> f_x_tmp(num_bins_x);
  std::vector<int>    samples(num_bins_x);
	const double bin_size_x = bins_x.get_bin_size();
  const double fac = 1.0/bin_size_x;  // normalize by bin size
  for ( int b=0; b<num_bins_x; ++b ) {
    double sum = std::accumulate( binned_weights_[b].begin(), binned_weights_[b].end(), 0.0 );
    p_x_tmp[b] = fac*sum;
    f_x_tmp[b] = -std::log(p_x_tmp[b]);
    samples[b] = binned_weights_[b].size();
  }

  f_x_ = FreeEnergyDistribution(bins_x, f_x_tmp, p_x_tmp, samples);
}