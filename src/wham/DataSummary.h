// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#pragma once
#ifndef DATA_SUMMARY_H
#define DATA_SUMMARY_H

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "InputParser.h"

class DataSummary
{
 public:
  DataSummary();

  DataSummary(
    const std::string& data_summary_file,
    const ParameterPack& input_pack
  );

  std::size_t getNumSimulations() const { return data_set_labels_.size(); }

  // Data for all simulations
  const std::vector<double>& get_t_min() const { return t_min_; }
  const std::vector<double>& get_t_max() const { return t_max_; }
  const std::vector<std::string>& getDataSetLabels() const { return data_set_labels_; }

  // Data for an individual simulation (by index)
  double get_t_min(const int i) const { return t_min_[i]; }
  double get_t_max(const int i) const { return t_max_[i]; }
  const std::string& indexToDataSetLabel(const int i) const { return data_set_labels_[i]; }

  // Data for an individual simulation (by data_set_label)
  double get_t_min(const std::string& data_set_label) const {
    return get_t_min( dataSetLabelToIndex(data_set_label) );
  }
  double get_t_max(const std::string& data_set_label) const {
    return get_t_max( dataSetLabelToIndex(data_set_label) );
  }
  const std::string& get_data_set_label(const std::string& data_set_label) const {
    return indexToDataSetLabel( dataSetLabelToIndex(data_set_label) );
  }

  // Given a data set label, return the corresponding index in the data summary
  int dataSetLabelToIndex(const std::string& data_set_label) const {
    auto pair_it = map_data_set_labels_to_indices_.find( data_set_label );
    if ( pair_it != map_data_set_labels_to_indices_.end() ) {
      return pair_it->second ;
    }
    else {
      throw std::runtime_error("Error: data set \'" + data_set_label + "\' is not present");
    }
  }

  const std::map<std::string, int>& get_data_set_label_map() const {
    return map_data_set_labels_to_indices_;
  }

 private:
  std::string data_summary_file_ = "", data_summary_path_ = "";

  int col_data_label_ = 0;
  int col_t_min_      = 3;
  int col_t_max_      = 4;
  int col_T_          = -1;

  std::vector<double> t_min_, t_max_;
  std::vector<std::string> data_set_labels_;

  // Maps data set labels to simulation indices
  std::map<std::string, int> map_data_set_labels_to_indices_;

  // TODO
  //std::vector<double> temperatures_;
};

#endif // ifndef DATA_SUMMARY_H
