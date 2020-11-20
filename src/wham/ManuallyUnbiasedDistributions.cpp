// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#include "ManuallyUnbiasedDistributions.hpp"
#include "LogSumExp.hpp"


ManuallyUnbiasedDistributions::ManuallyUnbiasedDistributions(
  const OrderParameter& x, const Wham& wham
):
  FreeEnergyDistributionSet(x, wham.getSimulationData())
{
  // Compute weights
  const auto log_w = wham.computeUnbiasedNonConsensusLogWeights();
  
  const auto& simulations = getData();
  const int num_simulations = simulations.size();
  std::vector<FreeEnergyDistribution> distributions;
  distributions.reserve(num_simulations);
  
  // Set up buffer
  const auto& bins_x     = x.getBins();
  const int   num_bins_x = bins_x.getNumBins();
  
  std::vector<std::vector<double>> binned_values(num_bins_x);

  for ( int j=0; j<num_simulations; ++j ) {
    // Sort samples by bin
    for ( auto& bin : binned_values ) {
      bin.resize(0);
    }
    const auto& x_j = simulations[j].getTimeSeriesForOrderParameter(x.getName());
    const int num_samples = x_j.size();
    for ( int i=0; i<num_samples; ++i ) {
      int b = bins_x.findBin( x_j[i] );
      if ( b >= 0 ) {
        binned_values[b].push_back( log_w[j][i] );
      }
    }

    // Estimate unbiased distribution
    std::vector<int>    sample_counts(num_bins_x);
    std::vector<double> p_x(num_bins_x);
    std::vector<double> f_x(num_bins_x);
    const double fac = std::log(1.0/num_samples);
    for ( int b=0; b<num_bins_x; ++b ) {
      sample_counts[b] = binned_values[b].size();
      if ( sample_counts[b] > 0 ) {
        f_x[b] = -fac - numeric::logSumExp( binned_values[b] );
        p_x[b] = std::exp( -f_x[b] );
      }
      else {
        f_x[b] = -1.0;
        p_x[b] = 0.0;
      }
    }

    // Organize results
    distributions.emplace_back(bins_x, f_x, p_x, sample_counts);
  }

  setDistributions( std::move(distributions) );
}
  
  
 std::string ManuallyUnbiasedDistributions::getHeader() const
 {
  const auto& x_name = getOrderParameter().getName();

  std::stringstream ss;
  ss << "Unbiased free energy distributions:  F_{0,i}(" << x_name  << ") [k_B*T]";

  return ss.str();
}