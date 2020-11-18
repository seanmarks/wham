// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#ifndef ESTIMATOR_F_X_HPP
#define ESTIMATOR_F_X_HPP

#include <string>

#include "Distribution.h"
#include "OrderParameter.h"
#include "OrderParameterRegistry.h"
#include "WhamEstimator.hpp"


class Estimator_F_x : public WhamEstimator
{
 public:
  Estimator_F_x(
    const OrderParameter& x
    //const std::string& x_name
    //const OrderParameterRegistry& registry
  );

  // TODO:
  void calculate(
    const Wham& wham,
    const std::vector<double>& u_bias_as_k,
    const double f_bias_k
    //const std::vector<Simulation>& data
  );

  const Distribution& get_f_x() const {
    return f_x_;
  }


 private:
  const OrderParameter& x_;
  //std::string x_name_;
  Bins bins_;

  using Buffer = std::vector<double>;
  std::vector<Buffer> binned_weights_;

  std::vector<double> weights_;

  Distribution f_x_;
};

#endif // ifndef ESTIMATOR_F_X_HPP