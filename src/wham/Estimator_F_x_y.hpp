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

  // Perform calculations and store F(x,y)
  void calculate(
    const Wham& wham,
    const std::vector<double>& u_bias_as_r,
    const double f_bias_r
  );

  const FreeEnergyDistribution2D& get_f_x_y() const {
    return f_x_y_;
  }

  // TODO: modifier to decide where files go?
  void saveResults() const;


 private:
  const OrderParameter& x_;
  const OrderParameter& y_;

  // Working variable type
  using Buffer = std::vector<double>;

  Matrix<Buffer> binned_weights_;
  Buffer weights_;

  FreeEnergyDistribution2D f_x_y_;
};

#endif // ifndef ESTIMATOR_F_X_Y_HPP