// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#include "WhamStatistics.hpp"


WhamStatistics::WhamStatistics(
  const Wham& wham, const BiasedDistributions& f_x_biased,
  const RebiasedDistributions& f_x_rebiased
):
  wham_(wham),
  f_x_biased_(f_x_biased),
  f_x_rebiased_(f_x_rebiased),
  x_(f_x_biased.getOrderParameter())
{
  // Consistency checks
  FANCY_ASSERT( x_.getName() == f_x_rebiased_.getOrderParameter().getName(), "inconsistent input" )
  FANCY_ASSERT( f_x_biased_.getDataSetLabels() == f_x_rebiased.getDataSetLabels(),
                "mismatch in data set labels" );

  const auto& simulations = wham.getSimulationData();
  int num_simulations = simulations.size();
  const auto& x_name = x_.getName();
  avg_x_.resize(num_simulations);
  var_x_.resize(num_simulations);
  info_entropy_.resize(num_simulations);
  for ( int i=0; i<num_simulations; ++i ) {
    const auto& x_i = simulations[i].getTimeSeriesForOrderParameter(x_name);
    avg_x_[i] = x_i.average();
    var_x_[i] = x_i.variance();

    
    const auto& f_rebiased_i = f_x_rebiased.getDistributions()[i];
    const auto& f_biased_i = f_x_biased.getDistributions()[i];
    info_entropy_[i] = FreeEnergyDistribution::computeInformationEntropy(f_rebiased_i, f_biased_i);
  }
}


void WhamStatistics::print(std::string file_name) const
{
  const auto& name = x_.getName();
  if ( file_name.empty() ) {
    file_name = "stats_" + name + ".out";
  }
  
  std::ofstream ofs(file_name);

  // Header
  const int num_simulations = info_entropy_.size();
  const std::string spacer("   ");
  ofs << "# data_set"
      << spacer << "avg(" << name << ")"
      << spacer << "var(" << name << ")"
      << spacer << "info_entropy(biased/rebiased)\n";

  const auto& labels = f_x_biased_.getDataSetLabels();
  for ( int i=0; i<num_simulations; ++i ) {
    ofs << labels[i] << "\t" << avg_x_[i] << "\t" << var_x_[i]
        << "\t" << info_entropy_[i] << "\n";
  }

  ofs.close();
}