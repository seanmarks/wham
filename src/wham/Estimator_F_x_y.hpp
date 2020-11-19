// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#ifndef ESTIMATOR_F_X_Y_HPP
#define ESTIMATOR_F_X_Y_HPP

#include <string>

#include "FreeEnergyDistribution2D.hpp"
#include "OrderParameter.h"
#include "OrderParameterRegistry.h"
#include "WhamEstimator.hpp"

// Estimate a consensus 2D free energy distribution, F(x,y), using WHAM results
class Estimator_F_x_y : public WhamEstimator
{
 public:
  template<typename T>
  using Matrix = FreeEnergyDistribution2D::Matrix<T>;

  Estimator_F_x_y(
    const OrderParameter& x,
    const OrderParameter& y
  );

 // Returns the consensus distribution, F^{WHAM}(x,y) [kBT]
  const FreeEnergyDistribution2D& getFxy() const {
    return f_x_y_;
  }

  // TODO: modifier to decide where files go?
  void saveResults() const;

 protected:
  virtual
  void calculateUsingStoredWeights(const Wham& wham) override;


 private:
  const OrderParameter& x_;
  const OrderParameter& y_;

  // Results
  FreeEnergyDistribution2D f_x_y_;
  
  // Buffers
  Matrix<Buffer> binned_weights_;
};

#endif // ifndef ESTIMATOR_F_X_Y_HPP