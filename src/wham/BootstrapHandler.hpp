// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#pragma once
#ifndef BOOTSTRAP_HANDLER_HPP
#define BOOTSTRAP_HANDLER_HPP

#include "PointEstimator.h"
#include "Wham.h"

// Handles bootstrap error calculations for some output quantity that
// is obtained from WHAM
// - Note: non-virtual interface (NVI) idiom in use
class BootstrapHandler
{
 public:
  BootstrapHandler() = default;

  virtual ~BootstrapHandler() = default;

  void addSample(const Wham& wham) {
    this->addSampleImpl(wham);
  }

  void finalize() {
    this->finalizeImpl();
  }

 protected:
  virtual
  void addSampleImpl(const Wham& wham) = 0;
  
  virtual
  void finalizeImpl() = 0;
};

#endif // ifndef BOOTSTRAP_HANDLER_HPP