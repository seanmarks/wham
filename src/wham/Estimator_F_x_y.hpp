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

  Estimator_F_x_y() = delete;

  Estimator_F_x_y(
    const OrderParameter& x,
    const OrderParameter& y
  );

 // Returns the consensus distribution, F^{WHAM}(x,y) [kBT]
  const FreeEnergyDistribution2D& getFxy() const {
    return f_x_y_;
  }

  // Saves the estimate to files with standard names
  // - TODO:
  //   - modifier to decide where files go?
  //   - separate functions for each file?
  void saveResults() const;

  // Returns a handle to the first order parameter, 'x'
  const OrderParameter& get_x() const {
    return x_;
  }

  // Returns a handle to the second order parameter, 'y'
  const OrderParameter& get_y() const {
    return y_;
  }
  

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