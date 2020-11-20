// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#ifndef ESTIMATOR_AVG_X_GIVEN_Y_HPP
#define ESTIMATOR_AVG_X_GIVEN_Y_HPP

#include "Estimator_F_x.hpp"
#include "FreeEnergyDistribution.hpp"
#include "OrderParameter.h"
#include "WhamEstimator.hpp"

// Estimates <x>_{0,y}: the average value of 'x' in the unbiased ensemble,
// given the fixed 'y'-value
// - Computes this average for each binned y-value
class Estimator_Avg_x_Given_y : public WhamEstimator
{
 public:
  Estimator_Avg_x_Given_y(
    const OrderParameter& x,
    const OrderParameter& y
  );

  const OrderParameter& get_x() const noexcept {
    return x_;
  }
  const OrderParameter& get_y() const noexcept {
    return y_;
  }
  const std::vector<double>& getValues() const noexcept {
    return avg_x_given_y_;
  }
  const std::vector<int>& getSampleCounts() const noexcept {
    return sample_counts_;
  }

  // TODO: move to separate class?
  void saveResults(std::string file_name = "") const;


 protected:
  virtual
  void calculateUsingStoredWeights(const Wham& wham) override;

 private:
  const OrderParameter& x_;
  const OrderParameter& y_;

  std::vector<double> avg_x_given_y_;
  std::vector<int>    sample_counts_;

  // Working variables
  Estimator_F_x est_f_y_;
  std::vector<WhamEstimator::Buffer> binned_weights_;
};

#endif // ifndef ESTIMATOR_AVG_X_GIVEN_Y_HPP