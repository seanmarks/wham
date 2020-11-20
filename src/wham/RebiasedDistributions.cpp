// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "RebiasedDistributions.hpp"
#include "Estimator_F_x.hpp"

RebiasedDistributions::RebiasedDistributions(
  const OrderParameter& x, const Wham& wham
):
  FreeEnergyDistributionSet(x, wham.getSimulationData())
{
  // Allocate space
  const auto& simulations = getData();
  const int num_simulations = simulations.size();
  std::vector<FreeEnergyDistribution> distributions;
  distributions.reserve(num_simulations);

  // Construct distributions
  Estimator_F_x est_f_x(x);
  for ( const auto& s : simulations ) {
    est_f_x.calculate(wham, s.getDataSetLabel());
    distributions.emplace_back( est_f_x.get_f_x() );
  }

  // Pass output to base class
  setDistributions( std::move(distributions) );
}


std::string RebiasedDistributions::getHeader() const
{
  const auto& x_name = getOrderParameter().getName();

  std::stringstream ss;
  ss << "\"Rebiased\" free energy distributions:  F_{rebias,i}(" << x_name << ") [k_B*T]";

  return ss.str();
}