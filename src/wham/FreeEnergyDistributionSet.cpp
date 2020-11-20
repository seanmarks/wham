// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#include "FreeEnergyDistributionSet.hpp"


FreeEnergyDistributionSet::FreeEnergyDistributionSet(
  const OrderParameter& x, const std::vector<Simulation>& data
):
  x_(x),
  data_(data)
{
  for ( const auto& s : data ) {
    data_set_labels_.push_back( s.getDataSetLabel() );
  }

  const int num_simulations = data.size();
  f_x_.reserve(num_simulations);
}


void FreeEnergyDistributionSet::print(const std::string& file_name) const
{
  std::stringstream ss;
  ss << "# " << this->getHeader() << "\n";
  
  ss << "# Data sets (by column)\n";
  const int num_simulations = data_set_labels_.size();
  for ( int i=0; i<num_simulations; ++i ) {
    ss << "# " << i+2 << ": " << data_set_labels_[i] << "\n";
  }
  const auto& name = x_.getName();
  ss << "#\n"
     << "# " << name << " | F(" << name << ") [kBT]\n";

  // FIXME: MOVE FUNCTIONALITY HERE
  FreeEnergyDistribution::printSet( f_x_, file_name, ss.str() );
}