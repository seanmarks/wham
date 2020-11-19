#include "WhamEstimator.hpp"


void WhamEstimator::calculate(
  const Wham& wham, const std::string& data_set_label
)
{
  weights_ = wham.calculateWeightsForSimulation(data_set_label);
  calculateUsingStoredWeights(wham);
}


void WhamEstimator::calculate(const Wham& wham)
{
  weights_ = wham.calculateWeightsForUnbiasedEnsemble();
  calculateUsingStoredWeights(wham);
}