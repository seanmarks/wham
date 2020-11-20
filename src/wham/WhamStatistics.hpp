// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#pragma once
#ifndef WHAM_STATISTICS_HPP
#define WHAM_STATISTICS_HPP

#include "BiasedDistributions.hpp"
#include "RebiasedDistributions.hpp"

class WhamStatistics
{
 public:
  WhamStatistics(
    const Wham&                  wham,
    const BiasedDistributions&   f_x_biased,
    const RebiasedDistributions& f_x_rebiased
  );

  void print(std::string file_name = "") const;


 private:
  const Wham&                  wham_;
  const BiasedDistributions&   f_x_biased_;
  const RebiasedDistributions& f_x_rebiased_;

  const OrderParameter& x_;
  std::vector<double> avg_x_, var_x_;
  std::vector<double> info_entropy_;
};

#endif // ifndef WHAM_STATISTICS_HPP