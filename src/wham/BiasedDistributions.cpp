// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "BiasedDistributions.hpp"


BiasedDistributions::BiasedDistributions(
  const OrderParameter& x, const std::vector<Simulation>& data
):
  FreeEnergyDistributionSet(x, data)
{
  // Allocate space
  const auto& simulations = getData();
  const int num_simulations = simulations.size();
  std::vector<FreeEnergyDistribution> distributions;
  distributions.reserve(num_simulations);

  // Construct distributions
  const auto& x_name = x.getName();
  for ( int i=0; i<num_simulations; ++i ) {
    const auto& x_i = simulations[i].getTimeSeriesForOrderParameter(x_name);
    distributions.emplace_back( x.getBins(), x_i );
  }

  // Pass output to base class
  setDistributions( std::move(distributions) );
}


std::string BiasedDistributions::getHeader() const
{
  const auto& x_name = getOrderParameter().getName();

  std::stringstream ss;
	ss << "Biased free energy distributions:  F_i(" << x_name << ") [k_B*T]";

  return ss.str();
}