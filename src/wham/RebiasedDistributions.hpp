// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#ifndef REBIASED_DISTRIBUTIONS_HPP
#define REBIASED_DISTRIBUTIONS_HPP

#include "FreeEnergyDistributionSet.hpp"
#include "Wham.h"


// Calculates consensus biased free energy distributions using WHAM results
// - These are useful for checking the internal consistency of WHAM output
class RebiasedDistributions : public FreeEnergyDistributionSet
{
 public:
  RebiasedDistributions(const OrderParameter& x, const Wham& wham);

 protected:
  virtual
  std::string getHeader() const override;
};

#endif // ifndef REBIASED_DISTRIBUTIONS_HPP