#include "Bootstrap_F_x.hpp"


Bootstrap_F_x::Bootstrap_F_x(
  const OrderParameter& x, FreeEnergyDistribution& f_x,
  const int num_samples_to_reserve
):
  BootstrapHandler(),
  x_(x),
  f_x_out_(f_x),
  est_f_x_(x)
{
  estimates_.resize( x_.getBins().getNumBins() );
  for ( auto& estimator : estimates_ ) {
    estimator.reserve(num_samples_to_reserve);
  }
}


void Bootstrap_F_x::addSampleImpl(const Wham& wham)
{
  est_f_x_.calculate(wham);
  const auto& bootstrap_f_x = est_f_x_.get_f_x();
  const auto& f_x    = bootstrap_f_x.get_F_x();
  const auto& counts = bootstrap_f_x.getSampleCounts();

  const int num_bins_x = estimates_.size();
  for ( int b=0; b<num_bins_x; ++b ) {
    // Only save this sample for if F(x_b) if its value is finite (i.e. bin has samples in it)
    // - Otherwise, statistics over the samples will be corrupted
    if ( counts[b] > 0 ) {
      estimates_[b].addSample( f_x[b] );
    }
  }
}


void Bootstrap_F_x::finalizeImpl()
{
  const int num_bins = estimates_.size();
  std::vector<double> err_f_x( num_bins );
  for ( int b=0; b<num_bins; ++b ) {
    err_f_x[b] = estimates_[b].std_dev();
  }

  f_x_out_.setErrors_F_x(err_f_x);
}