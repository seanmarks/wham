// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#pragma once
#ifndef BOOTSTRAP_F_X_HPP
#define BOOTSTRAP_F_X_HPP

#include "BootstrapHandler.hpp"
#include "Estimator_F_x.hpp"
#include "FreeEnergyDistribution.hpp"
#include "OrderParameter.h"

class Bootstrap_F_x : public BootstrapHandler
{
 public:
  Bootstrap_F_x(
    const OrderParameter& x,
    FreeEnergyDistribution& f_x,
    const int num_samples_to_reserve = 100
  );


 protected:
  void addSampleImpl(const Wham& wham) override;

  void finalizeImpl() override;


 private:
  const OrderParameter& x_;

  // Where errors should be placed when finalize() is called
  FreeEnergyDistribution& f_x_out_;
  
  // Computes estimates of F(x) given WHAM output
  Estimator_F_x est_f_x_;

  // Organizes bootstrap samples
  std::vector<PointEstimator<double>> estimates_;
};

#endif // ifndef BOOTSTRAP_F_X_HPP