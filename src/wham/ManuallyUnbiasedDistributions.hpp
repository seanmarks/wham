// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#ifndef MANUALLY_UNBIASED_FREE_ENERGY_HPP
#define MANUALLY_UNBIASED_FREE_ENERGY_HPP

#include "FreeEnergyDistributionSet.hpp"
#include "OrderParameter.h"
#include "Wham.h"

// "Manually unbiased" free energy distributions, F_0^{(i)}(x),
// obtained using only data from a single simulation (non-consensus)
class ManuallyUnbiasedDistributions : public FreeEnergyDistributionSet
{
 public:
  // Constructs the set using the same simulation data used for the giventhe WHAM instance
  ManuallyUnbiasedDistributions(const OrderParameter& x, const Wham& wham);

 protected:
  virtual
  std::string getHeader() const override;
};

#endif // ifndef MANUALLY_UNBIASED_FREE_ENERGY_HPP