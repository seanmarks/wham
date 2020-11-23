// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#pragma once
#ifndef AVG_X_GIVEN_Y_HPP
#define AVG_X_GIVEN_Y_HPP

#include "OrderParameter.h"
#include "OutputVectorWithBins.hpp"


// Represents a set of values, <x>_{0,y}: the average value of 'x' in the unbiased ensemble,
// given the fixed 'y'-value
class Avg_x_Given_y : public OutputVectorWithBins
{
 public:
  // Constructs empty without errors
  Avg_x_Given_y(
    const OrderParameter& x,
    const OrderParameter& y
  );

  // Constructs without errors
  Avg_x_Given_y(
    const OrderParameter& x,
    const OrderParameter& y,
    const Vector<double>& values,
    const Vector<int>& sample_counts
  );
  
  // Constructs with errors
  Avg_x_Given_y(
    const OrderParameter& x,
    const OrderParameter& y,
    const Vector<double>& values,
    const Vector<double>& errors,
    const Vector<int>& sample_counts
  );


  //----- Interface -----//

  const OrderParameter& get_x() const noexcept {
    return x_;
  }
  const OrderParameter& get_y() const noexcept {
    return y_;
  }


 protected:
  virtual
  std::string getNameOfBins() const override {
    return y_.getName();
  }

  virtual
  std::string getNameOfVector() const override {
    std::stringstream ss;
    ss << "<" << x_.getName() << "|" << y_.getName() << ">";
    return ss.str();
  }


 private:
  const OrderParameter& x_, y_;

};

#endif // ifndef AVG_X_GIVEN_Y_HPP