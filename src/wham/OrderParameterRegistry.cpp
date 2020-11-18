// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "OrderParameterRegistry.h"


OrderParameterRegistry::OrderParameterRegistry()
{}


OrderParameterRegistry::OrderParameterRegistry(const ParameterPack& input_pack, const DataSummary& data_summary):
	data_summary_ptr_(&data_summary)
{
	using KeyType = ParameterPack::KeyType;

	std::vector<const ParameterPack*> op_pack_ptrs = 
			input_pack.findParameterPacks("OrderParameter", KeyType::Required);
	int num_ops = op_pack_ptrs.size();
	if ( num_ops < 1 ) {
		throw std::runtime_error("no order parameters were registered");
	}

	names_.resize(num_ops);
	time_series_lists_.resize(num_ops);
	time_series_files_.resize(num_ops);
	data_columns_.resize(num_ops);
	file_columns_.resize(num_ops);

	const int num_simulations = data_summary.getNumSimulations();

	std::string name;
	int file_col, data_col;
	for ( int p=0; p<num_ops; ++p ) {
		const auto& op_input_pack = *(op_pack_ptrs[p]);

		// Name (must be unique)
		op_input_pack.readString("Name", KeyType::Required, name);
		names_[p] = name;

		// Files list
		std::vector<std::string> tokens;
		op_input_pack.readVector("FilesList", KeyType::Required, tokens);
		time_series_lists_[p] = tokens[0];

		if ( tokens.size() > 1 ) {
			// User also provided the column in the file which has the 
			// time series file locations
			file_col = std::stoi( tokens[1] );
			--file_col;  // input is indexed from 1
			if ( file_col < 0 ) {
				std::stringstream err_ss;
				err_ss << "Error setting up order parameter " << name << "\n"
							 << "  Column number for data must be a positive integer\n";
				throw std::runtime_error( err_ss.str() );
			}
		}
		else {
			file_col = 1;
		}
		file_columns_[p] = file_col;

		// Get list of time series files
		FileSystem::readFilesList(time_series_lists_[p], file_columns_[p], time_series_files_[p]);
		int num_time_series_files = time_series_files_[p].size();
		if ( num_time_series_files != num_simulations ) {
			std::stringstream err_ss;
			err_ss << "Error setting up order parameter " << name << "\n"
						 << "  Found " << num_time_series_files << " time series files (expected " << num_simulations << ")\n";
			throw std::runtime_error( err_ss.str() );
		}

		// Column *in* each time series file that contains the time series data
		op_input_pack.readNumber("DataColumn", KeyType::Required, data_col);
		if ( data_col < 1 ) {
			std::stringstream err_ss;
			err_ss << "Error setting up order parameter " << name << "\n"
						 << "  DataColumn must be a positive integer\n";
			throw std::runtime_error( err_ss.str() );
		}
		--data_col;  // input is indexed from 1
		data_columns_[p] = data_col;

#ifdef DEBUG
		std::cout << "OrderParameter: " << names_[p] << "\n"
							<< "  From files list " << time_series_lists_[p] << "\n"
							<< "  Found " << num_time_series_files << " file names\n";
		for ( int i=0; i<num_time_series_files; ++i ) {
			std::cout << "    " << time_series_files_[p][i] << "\n";
		}
#endif /* DEBUG */

		// Record the mapping from name to index
		auto ret = map_op_names_to_indices_.insert( std::make_pair( name, p ) );
		if ( ret.second == false ) {
			throw std::runtime_error("order parameter \"" + name  + "\" was defined more than once");
		}
	}

	// Reorganize by simulation
	simulation_files_.resize(num_simulations);
	for ( int i=0; i<num_simulations; ++i ) {
		simulation_files_[i].resize(num_ops);
		for ( int p=0; p<num_ops; ++p ) {
			simulation_files_[i][p] = time_series_files_[p][i];
		}
	}
}
