// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#pragma once
#ifndef ORDER_PARAMETER_REGISTRY_H
#define ORDER_PARAMETER_REGISTRY_H

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Assert.hpp"
#include "DataSummary.h"
#include "InputParser.h"


// Tracks the order parameters (OPs) present, and the location of their time series data
//
// Terminology:
// - "time series list" = file contining a list of time series files for the OP
// - "file column"      = column in the time series list with the time series locations (default: 2nd column)
// - "time series file" = contains the actual time series data
// - "data column"      = column in the time series file with the OP data (default: 2nd column)
class OrderParameterRegistry
{
 public:
  using Strings = std::vector<std::string>;
  using const_iterator = typename Strings::const_iterator;

  OrderParameterRegistry() = default;

  OrderParameterRegistry(const ParameterPack& input_pack, const DataSummary& data_summary);

  // Returns the index corresponding to the given OP name
  int nameToIndex(const std::string& name) const;

  // Returns the OP indices corresponding to the given OP names
  std::vector<int> namesToIndices(const std::vector<std::string>& names) const;

  // Returns the number of OPs registered
  std::size_t getNumRegistered() const noexcept {
    return names_.size();
  }

  // Returns the number of OPs registered
  std::size_t size() const noexcept {
    return names_.size();
  }

  // Iterate over names
  const_iterator begin() const noexcept {
    return names_.begin();
  }
  const_iterator end() const noexcept {
    return names_.end();
  }

  // Returns a list of all OP names
  const Strings& getNames() const noexcept {
    return names_;
  }

  // Access registry info by OP name
  const Strings& getTimeSeriesFiles(const std::string& name) const {
    return getTimeSeriesFiles( nameToIndex(name) );
  }
  int getTimeSeriesDataCol(const std::string& name) const {
    return getTimeSeriesDataCol( nameToIndex(name) );
  }
  int getTimeSeriesFileCol(const std::string& name) const {
    return getTimeSeriesFileCol( nameToIndex(name) );
  }

  // Get list of files by simulation index/data set label
  const Strings& getSimulationFiles(const int i) const {
    return simulation_files_[i];
  }
  const Strings& getSimulationFiles(const std::string& data_set_label) const {
    FANCY_ASSERT( data_summary_ptr_ != nullptr, "missing DataSummary" );
    return getSimulationFiles( data_summary_ptr_->dataSetLabelToIndex(data_set_label) );
  }


  // Access registry info by OP index
  const std::string& indexToName(const int p) const { 
    return names_[p]; 
  }
  const Strings& getTimeSeriesFiles(const int p) const {
    return time_series_files_[p];
  }
  int getTimeSeriesDataCol(const int p) const {
    return data_columns_[p];
  }
  int getTimeSeriesFileCol(const int p) const {
    return file_columns_[p];
  }


 private:
  const DataSummary* data_summary_ptr_ = nullptr;

  Strings names_;

  Strings time_series_lists_;
  std::vector<int> file_columns_;

  std::vector<Strings> time_series_files_;  // outer = by OP, inner = by sim
  std::vector<int> data_columns_;

  // Files arranged by simulation
  std::vector<Strings> simulation_files_;

  // Maps names of OrderParameters top their indices in the order_parameters_ vector
  std::map<std::string, int> map_op_names_to_indices_;
};

#endif // ifndef ORDER_PARAMETER_REGISTRY_H
