// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#ifndef BIASED_DISTRIBUTIONS_HPP
#define BIASED_DISTRIBUTIONS_HPP

#include "OrderParameter.h"
#include "Simulation.h"


// TODO: ?
class BiasedDistributions
{
 public:
  BiasedDistributions(const OrderParameter& x, const std::vector<Simulation>& data);

  
  // Returns the distributions (one for each simulation/sampled ensemble)
  const std::vector<FreeEnergyDistribution>& getDistributions() const noexcept {
    return distributions_;
  }


 private:
  const OrderParameter& x_;

  std::vector<FreeEnergyDistribution> distributions_;
  std::vector<std::string> data_set_labels_;
};

#endif // ifndef BIASED_DISTRIBUTIONS_HPP