// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#ifndef BOOTSTRAP_BIASING_FREE_ENERGIES_HPP
#define BOOTSTRAP_BIASING_FREE_ENERGIES_HPP

#include "BootstrapHandler.hpp"
#include "PointEstimator.h"

class Bootstrap_BiasingFreeEnergies : public BootstrapHandler
{
 public:
  Bootstrap_BiasingFreeEnergies(
    std::vector<double>& err_f_out,
    const int num_ensembles,
    const int num_samples_to_reserve = 100
  );

  const std::vector<double>& getSamples(const int k) {
    const int num_ensembles = estimates_.size();
    FANCY_ASSERT( k > 0 && k < num_ensembles, "index out of bounds" );
    return estimates_[k].get_samples();
  }

 protected:
  void addSampleImpl(const Wham& wham) override;

  void finalizeImpl() override;

 private:
  // Where to place output when finalize() is called
  std::vector<double>& err_f_out_;
  std::vector<PointEstimator<double>> estimates_;
};

#endif // ifndef BOOTSTRAP_BIASING_FREE_ENERGIES_HPP