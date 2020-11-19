
#ifndef FREE_ENERGY_DISTRIBUTION_2D_HPP
#define FREE_ENERGY_DISTRIBUTION_2D_HPP

#include "Bins.h"
#include "Matrix.hpp"
#include "FreeEnergyDistribution.hpp"


// Two-dimensional free energy distribution, F(x,y)
class FreeEnergyDistribution2D
{
 public:
  using Real = double;

  template<typename T>
  using Matrix = numeric::Matrix<T>;

  FreeEnergyDistribution2D() = default;

  // Construct empty
  FreeEnergyDistribution2D(const Bins& bins_x, const Bins& bins_y);

  // With specified contents
  FreeEnergyDistribution2D(
    const Bins& bins_x, const Bins& bins_y,
    const Matrix<Real>& f_x_y, const Matrix<Real>& p_x_y,
    const Matrix<int>& sample_counts
  );


  //----- Interface -----//

  const Matrix<Real>& get_p_x_y() const noexcept {
    return p_x_y_;
  }
  const Matrix<Real>& get_f_x_y() const noexcept {
    return f_x_y_;
  }
  const Matrix<int>& getSampleCounts() const noexcept {
    return sample_counts_;
  }


 private:
  Bins bins_x_, bins_y_;

  Matrix<Real> f_x_y_;  // free energy (in kBT)
  Matrix<Real> p_x_y_;
  Matrix<int>  sample_counts_;
};

#endif // ifndef FREE_ENERGY_DISTRIBUTION_2D_HPP