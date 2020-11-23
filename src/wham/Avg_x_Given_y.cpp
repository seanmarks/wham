// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "Avg_x_Given_y.hpp"


Avg_x_Given_y::Avg_x_Given_y(const OrderParameter& x, const OrderParameter& y):
  OutputVectorWithBins(y.getBins()),
  x_(x), y_(y)
{
}


Avg_x_Given_y::Avg_x_Given_y(
  const OrderParameter& x, const OrderParameter& y,
  const Vector<double>& values, const Vector<int>& sample_counts
): 
  OutputVectorWithBins(y.getBins(), values, sample_counts),
  x_(x), y_(y)
{
}


Avg_x_Given_y::Avg_x_Given_y(
  const OrderParameter& x, const OrderParameter& y,
  const Vector<double>& values, const Vector<double>& errors,
  const Vector<int>& sample_counts
):
  Avg_x_Given_y(x, y, values, sample_counts)
{
  setErrors(errors);
}