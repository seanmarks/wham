#include "FreeEnergyDistribution2D.hpp"
  
FreeEnergyDistribution2D::FreeEnergyDistribution2D(
  const Bins& bins_x, const Bins& bins_y
):
  bins_x_(bins_x), bins_y_(bins_y)
{
  const std::array<int,2> shape = {{ bins_x_.getNumBins(), bins_y_.getNumBins() }};
  f_x_y_.setShape(shape);
  p_x_y_.setShape(shape);
  sample_counts_.setShape(shape);
}


FreeEnergyDistribution2D::FreeEnergyDistribution2D(
  const Bins& bins_x, const Bins& bins_y,
  const Matrix<Real>& f_x_y, const Matrix<Real>& p_x_y,
  const Matrix<int>& sample_counts
):
  bins_x_(bins_x),
  bins_y_(bins_y),
  f_x_y_(f_x_y),
  p_x_y_(p_x_y),
  sample_counts_(sample_counts)
{}
	