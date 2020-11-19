// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#include "ManuallyUnbiasedDistributions.hpp"
#include "LogSumExp.hpp"


ManuallyUnbiasedDistributions::ManuallyUnbiasedDistributions(
  const OrderParameter& x, const Wham& wham
):
  x_(x)
{
  const auto& simulations = wham.getSimulationData();
  const int num_simulations = simulations.size();

  const auto log_w = wham.computeUnbiasedNonConsensusLogWeights();

  //Vector<double> p_x(num_bins);
  distributions_.reserve(num_simulations);
  
  const auto& bins_x     = x_.getBins();
  const int   num_bins_x = bins_x.getNumBins();
  Vector<Vector<double>> binned_values(num_bins_x);

	for ( int j=0; j<num_simulations; ++j ) {
	  // Sort samples by bin
    for ( auto& bin : binned_values ) {
      bin.resize(0);
    }
		const auto& x_j = simulations[j].getTimeSeriesForOrderParameter(x_.getName());
		const int num_samples = x_j.size();
		for ( int i=0; i<num_samples; ++i ) {
			int b = bins_x.findBin( x_j[i] );
			if ( b >= 0 ) {
        binned_values[b].push_back( log_w[j][i] );
			}
		}

    Vector<int>    sample_counts(num_bins_x);
    Vector<double> p_x(num_bins_x);
    Vector<double> f_x(num_bins_x);
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

    distributions_.push_back( 
      FreeEnergyDistribution(bins_x, f_x, p_x, sample_counts)
    );
  }
}