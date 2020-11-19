
#ifndef FREE_ENERGY_DISTRIBUTION_2D_HPP
#define FREE_ENERGY_DISTRIBUTION_2D_HPP

#include "Matrix.hpp"

// Two-dimensional free energy distribution, F(x,y)
class FreeEnergyDistribution2D
{
 public:
  using Real = double;

 private:
  
  numeric::Matrix<Real> f_x_y_;
};

#endif // ifndef FREE_ENERGY_DISTRIBUTION_2D_HPP