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


// OrderParameterRegistry
// - Tracks the order parameters (OPs) present, and the origin of their time series data
class OrderParameterRegistry
{
 public:
	OrderParameterRegistry();

	OrderParameterRegistry(const ParameterPack& input_pack, const DataSummary& data_summary);

	// Returns the index
	int nameToIndex(const std::string& name) const {
		auto pair_it = map_op_names_to_indices_.find( name );
		if ( pair_it != map_op_names_to_indices_.end() ) {
			return pair_it->second ;
		}
		else {
			throw std::runtime_error("Error \'" + name + "\' is not registered");
		}
	}

	std::vector<int> namesToIndices(const std::vector<std::string>& names) const {
		const int num_ops = names.size();
		std::vector<int> indices(num_ops);

		std::transform( names.begin(), names.end(), indices.begin(),
			[=](const std::string& name) -> int {
				return nameToIndex(name);
			}
		);

		return indices;
	}

	std::size_t getNumRegistered() const {
		return names_.size();
	}
	const std::vector<std::string>& getNames() const { 
		return names_; 
	}

	// Access registry info by OP index
	const std::string& indexToName(const int p) const { 
		return names_[p]; 
	}

	const std::vector<std::string>& getTimeSeriesFiles(const int p) const {
		return time_series_files_[p];
	}

	int getTimeSeriesDataCol(const int p) const {
		return data_columns_[p];
	}

	int getTimeSeriesFileCol(const int p) const {
		return file_columns_[p];
	}

	// Access registry info by OP name
	const std::vector<std::string>& getTimeSeriesFiles(const std::string& name) const {
		return getTimeSeriesFiles( nameToIndex(name) );
	}

	int getTimeSeriesDataCol(const std::string& name) const {
		return getTimeSeriesDataCol( nameToIndex(name) );
	}

	int getTimeSeriesFileCol(const std::string& name) const {
		return getTimeSeriesFileCol( nameToIndex(name) );
	}

	// Get list of files by simulation index/data set label
	const std::vector<std::string>& getSimulationFiles(const int i) const {
		return simulation_files_[i];
	}

	const std::vector<std::string>& getSimulationFiles(const std::string& data_set_label) const {
		if ( data_summary_ptr_ != nullptr ) {
			return getSimulationFiles( data_summary_ptr_->dataSetLabelToIndex(data_set_label) );
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
