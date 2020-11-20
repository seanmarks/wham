
#pragma once
#ifndef FREE_ENERGY_DISTRIBUTION_HPP
#define FREE_ENERGY_DISTRIBUTION_HPP

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "Bins.h"
#include "TimeSeries.h"

// 1-D free energy distribution, F(x)
// TODO:
// - Better interface
class FreeEnergyDistribution {
 public:
  //----- Data -----//

  using VectorReal = std::vector<double>;

  Bins                bins_x;
  std::vector<double> f_x, error_f_x;
  std::vector<double> p_x;
  std::vector<int>    sample_counts;


  //----- Setup ----//

  FreeEnergyDistribution();

  FreeEnergyDistribution(const Bins& bins_x_in);

  FreeEnergyDistribution(const Bins& bins_x_in, const TimeSeries& time_series_x);
  
  
  FreeEnergyDistribution(
    const Bins& bins_x_in,  
    const VectorReal& f_x_in, const VectorReal& p_x_in, const std::vector<int>& samples
  ):
    bins_x(bins_x_in), f_x(f_x_in), p_x(p_x_in), sample_counts(samples)
  {
    // TODO: check for consistency of lengths
  }

  void print(
    const std::string& file, 
    const std::string& header
  ) const;

  void setErrors_F_x(const std::vector<double>& err) {
    FANCY_ASSERT(err.size() == f_x.size(), "length mismatch");
    error_f_x = err;
  }


  //----- Interface -----//

  // Prints a set of distributions to file
	// FIXME: MOVE
  static void printSet(
    const std::vector<FreeEnergyDistribution>& distributions,
    const std::string& file_name,
    const std::string& header,
    const bool shift_to_zero = true
  );

  const Bins& getBins() const noexcept {
    return bins_x;
  }

  const std::vector<int>& getSampleCounts() const noexcept {
    return sample_counts;
  }

  //----- Analysis -----//

  // Kullback-Leibler divergence (information entropy difference) between a
  // distribution and a reference distribution
  static double computeInformationEntropy(
    const FreeEnergyDistribution& dist,
    const FreeEnergyDistribution& ref
  );


  //----- Output Helper Functions -----//

  // When computed from a histogram, F(x) is finite only when num_samples > 0
  template<typename T, typename I>
  static bool is_f_x_finite(const T f_x, const I num_samples) {
    return ( num_samples > 0 ? true : false );
  }

  // Get the minimum of F(x), given that not all bins have samples
  template<typename T, typename I>
  static T get_min_f_x(const std::vector<T>& f_x, const std::vector<I>& num_samples) {
    T min      = std::numeric_limits<T>::max();
    I num_bins = num_samples.size();
    for ( int b=0; b<num_bins; ++b ) {
      if ( is_f_x_finite(f_x[b], num_samples[b]) ) {
        min = std::min( min, f_x[b] );
      }
    }
    return min;
  }

  template<typename T, typename I>
  static std::vector<T> shift_f_x_to_zero(const std::vector<T>& f_x, const std::vector<I>& num_samples) {
    T min_f_x = get_min_f_x(f_x, num_samples);

    // Subtract minimum
    std::vector<T> f_x_shifted(f_x);
    std::for_each( f_x_shifted.begin(), f_x_shifted.end(), [=](T& f) { f -= min_f_x; } );
    return f_x_shifted;
  }

  // If 'f_x' is finite, print it
  // Else print 'nan'
  // - Standardizes what is printed when a number becomes non-finite
  //   across different systems (otherwise, regtests with NaNs fail)
  // TODO Make an object and overload operator<< to make usage less clunky
  template<typename T, typename I>
  static void printFreeEnergyValue(std::ofstream& ofs, const T f_x, const I num_samples) {
    if ( is_f_x_finite(f_x, num_samples) ) { ofs << f_x;   }
    else                                   { ofs << "nan"; }    
  };
};

#endif // ifndef FREE_ENERGY_DISTRIBUTION_HPP
