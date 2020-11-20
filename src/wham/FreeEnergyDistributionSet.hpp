// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#ifndef FREE_ENERGY_DISTRIBUTION_SET_HPP
#define FREE_ENERGY_DISTRIBUTION_SET_HPP

#include "FreeEnergyDistribution.hpp"
#include "OrderParameter.h"


// Represents a set of single-variable free energy distributions, F_j(x),
// corresponding to a set of simulations (one each)
// - Derived classes must calculate distributions in their constructors
class FreeEnergyDistributionSet
{
 public:
  FreeEnergyDistributionSet(
    const OrderParameter& x,
    const std::vector<Simulation>& data
  );

  virtual ~FreeEnergyDistributionSet() = default;


  //----- Interface -----//

  // Returns the distributions (one for each simulation/sampled ensemble)
  const std::vector<FreeEnergyDistribution>& getDistributions() const noexcept {
    return f_x_;
  }

  // Returns a handle to the OrderParameter this set is for
  const OrderParameter& getOrderParameter() const noexcept {
    return x_;
  }

  // Returns the list of data sets used to construct the distributions
  const std::vector<std::string>& getDataSetLabels() const noexcept {
    return data_set_labels_;
  }

  // Prints the set of distributions to the file
  void print(const std::string& file_name) const;



protected:
  const std::vector<Simulation>& getData() const noexcept {
    return data_;
  }

  // Returns the distributions (one for each simulation/sampled ensemble)
  void setDistributions(std::vector<FreeEnergyDistribution>&& f_x) {
    int num_input    = f_x.size();
    int num_expected = data_set_labels_.size();
    FANCY_ASSERT( num_input == num_expected,
                  "input " << num_input << " distributions, expected " << num_expected );
    f_x_ = f_x;
  }

  // Returns a header representing this set, which should be appropriate for output files
  virtual
  std::string getHeader() const = 0;


 private:
  const OrderParameter& x_;

  const std::vector<Simulation>& data_;

  std::vector<std::string> data_set_labels_;
  std::vector<FreeEnergyDistribution> f_x_;
};

#endif // ifndef FREE_ENERGY_DISTRIBUTION_SET_HPP