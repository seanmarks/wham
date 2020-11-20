// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#ifndef ESTIMATOR_F_X_HPP
#define ESTIMATOR_F_X_HPP

#include <string>

#include "FreeEnergyDistribution.hpp"
#include "OrderParameter.h"
#include "WhamEstimator.hpp"


// Calculates a consensus free energy distribution, F^{WHAM}(x)
class Estimator_F_x : public WhamEstimator
{
 public:
  Estimator_F_x(
    const OrderParameter& x
  );

  // Returns the calculated distribution, F^{WHAM}(x) [kBT]
  const FreeEnergyDistribution& get_f_x() const {
    return f_x_;
  }

  const OrderParameter& getOrderParameter() const noexcept {
    return x_;
  }


 protected:
  virtual
  void calculateUsingStoredWeights(const Wham& wham) override;


 private:
  const OrderParameter& x_;

  // Output
  FreeEnergyDistribution f_x_;

  // Working variables
  std::vector<WhamEstimator::Buffer> binned_weights_;
};

#endif // ifndef ESTIMATOR_F_X_HPP