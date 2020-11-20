// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#ifndef REBIASED_DISTRIBUTIONS_HPP
#define REBIASED_DISTRIBUTIONS_HPP

#include "FreeEnergyDistributionSet.hpp"
#include "Wham.h"

// TODO:
class RebiasedDistributions : public FreeEnergyDistributionSet
{
 public:
  RebiasedDistributions(const OrderParameter& x, const Wham& wham);

 protected:
  virtual
  std::string getHeader() const override;
};

#endif // ifndef REBIASED_DISTRIBUTIONS_HPP