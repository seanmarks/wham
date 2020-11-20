// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#include "Bootstrap_BiasingFreeEnergies.hpp"


Bootstrap_BiasingFreeEnergies::Bootstrap_BiasingFreeEnergies(
  std::vector<double>& err_f_out,
  const int num_ensembles,
  const int num_samples_to_reserve
):
  BootstrapHandler(),
  err_f_out_(err_f_out)
{
  estimates_.resize(num_ensembles);
  for ( auto& estimator : estimates_ ) {
    estimator.reserve(num_samples_to_reserve);
  }
}


void Bootstrap_BiasingFreeEnergies::addSampleImpl(const Wham& wham)
{
  const auto& f_bias = wham.get_f_bias_opt();
  const int num_f = f_bias.size();

  const int num_ensembles = estimates_.size();
  FANCY_ASSERT( num_ensembles == num_f, "unexpected number of simulations" );

  for ( int i=0; i<num_ensembles; ++i ) {
    estimates_[i].addSample( f_bias[i] );
  }
}


void Bootstrap_BiasingFreeEnergies::finalizeImpl()
{
  const int num_ensembles = estimates_.size();
  err_f_out_.resize(num_ensembles);
  for ( int i=0; i<num_ensembles; ++i ) {
    err_f_out_[i] = estimates_[i].std_dev();
  }
}