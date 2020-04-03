// OrderParameterRegistry
// - Tracks the order parameters (OPs) present, and the origin of their time series data

#pragma once
#ifndef ORDER_PARAMETER_REGISTRY_H
#define ORDER_PARAMETER_REGISTRY_H

// Standard headers
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Project headers
#include "DataSummary.h"
#include "FileSystem.h"
#include "InputParser.h"

class OrderParameterRegistry
{
 public:
	OrderParameterRegistry();

	OrderParameterRegistry(const ParameterPack& input_pack, const DataSummary& data_summary);

	int get_index(const std::string& name) const {
		auto pair_it = map_op_names_to_indices_.find( name );
		if ( pair_it != map_op_names_to_indices_.end() ) {
			return pair_it->second ;
		}
		else {
			throw std::runtime_error("Error \'" + name + "\' is not registered");
		}
	}

	std::size_t getNumberOfOrderParameters() const {
		return names_.size();
	}
	const std::vector<std::string>& get_names() const { 
		return names_; 
	}

	// Access registry info by OP index
	const std::string& get_name(const int p) const { 
		return names_[p]; 
	}
	const std::vector<std::string>& get_time_series_files(const int p) const {
		return time_series_files_[p];
	}
	int get_time_series_data_col(const int p) const {
		return data_columns_[p];
	}
	int get_time_series_file_col(const int p) const {
		return file_columns_[p];
	}

	// Access registry info by OP name
	const std::vector<std::string>& get_time_series_files(const std::string& name) const {
		return get_time_series_files( get_index(name) );
	}
	int get_time_series_data_col(const std::string& name) const {
		return get_time_series_data_col( get_index(name) );
	}
	int get_time_series_file_col(const std::string& name) const {
		return get_time_series_file_col( get_index(name) );
	}

	// Get list of files by simulation index/data set label
	const std::vector<std::string>& get_simulation_files(const int i) const {
		return simulation_files_[i];
	}
	const std::vector<std::string>& get_simulation_files(const std::string& data_set_label) const {
		if ( data_summary_ptr_ != nullptr ) {
			return get_simulation_files( data_summary_ptr_->get_index(data_set_label) );
		}
		else {
			throw std::runtime_error("Error: OrderParameterRegistry does not have access to a DataSummary.");
		}
	}


 private:
	const DataSummary* data_summary_ptr_ = nullptr;

	std::vector<std::string> names_;

	std::vector<std::string> time_series_lists_;
	std::vector<int> file_columns_;

	std::vector<std::vector<std::string>> time_series_files_;  // outer = by OP, inner = by sim
	std::vector<int> data_columns_;

	// Files arranged by simulation
	std::vector<std::vector<std::string>> simulation_files_;

	// Maps names of OrderParameters top their indices in the order_parameters_ vector
	std::map<std::string, int> map_op_names_to_indices_;
};

#endif // ifndef ORDER_PARAMETER_REGISTRY_H
