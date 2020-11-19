// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#pragma once
#ifndef WHAM_ESTIMATOR_HPP
#define WHAM_ESTIMATOR_HPP

#include <string>

#include "Wham.h"

// Calculates a consensus estimate using WHAM results
class WhamEstimator
{
 public:
  virtual ~WhamEstimator() = default;

  // Calculates the estimate in the ensemble corresponding to the simulation with the
  // given data set label
  void calculate(const Wham& wham, const std::string& data_set_label);

  // Calculates the estimate in the unbiased ensemble
  void calculate(const Wham& wham);

  // TODO: calculate using generic set of bias values?
  //       an arbitrary bias

 protected:
  // Estimates the desired quantity using weights computed and stored by
  // the base class interface
  virtual void calculateUsingStoredWeights(const Wham& wham) = 0;

  // Returns the weights for each sample in the current ensemble, k:  W_{n,k}
  // - size = N = number of samples total (for WHAM)
  // - Note that these weights sum to 1
  const std::vector<double>& getWeights() const noexcept {
    return weights_;
  }

  using Buffer = std::vector<double>;

 private:
  std::vector<double> weights_;
};

#endif // ifndef WHAM_ESTIMATOR_HPP