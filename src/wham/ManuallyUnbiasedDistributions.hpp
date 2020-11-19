
#ifndef MANUALLY_UNBIASED_FREE_ENERGY_HPP
#define MANUALLY_UNBIASED_FREE_ENERGY_HPP

#include "FreeEnergyDistribution.hpp"
#include "OrderParameter.h"
#include "Wham.h"


class ManuallyUnbiasedDistributions
{
 public:
  template<typename T>
  using Vector = Wham::Vector<T>;

  ManuallyUnbiasedDistributions(const OrderParameter& x, const Wham& wham);


 private:
  const OrderParameter& x_;

  std::vector<FreeEnergyDistribution> distributions_;
};

#endif // ifndef MANUALLY_UNBIASED_FREE_ENERGY_HPP