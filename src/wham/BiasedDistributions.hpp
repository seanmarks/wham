// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#ifndef BIASED_DISTRIBUTIONS_HPP
#define BIASED_DISTRIBUTIONS_HPP

#include "FreeEnergyDistributionSet.hpp"


// TODO:
class BiasedDistributions : public FreeEnergyDistributionSet
{
 public:
  BiasedDistributions(const OrderParameter& x, const std::vector<Simulation>& data);

 protected:
  virtual
  std::string getHeader() const override;
};

#endif // ifndef BIASED_DISTRIBUTIONS_HPP