#include "OrderParameter.h"

OrderParameter::OrderParameter(
		const ParameterPack& input_pack, const std::vector<Range>& production_phases,
		const bool use_floored_times
):
	file_col_(0),
	bins_()
{
	using KeyType = ParameterPack::KeyType;

	input_pack.readString("Name", KeyType::Required, name_);

	// Files list
	std::vector<std::string> tokens;
	input_pack.readVector("FilesList", KeyType::Required, tokens);
	time_series_list_ = tokens[0];
	if ( tokens.size() > 1 ) {
		// User also provided the column in the file which has the 
		// time series file locations
		file_col_ = std::stoi( tokens[1] );
		--file_col_;  // input is indexed from 1
	}

	input_pack.readNumber("DataColumn", KeyType::Required, data_col_);
	--data_col_;  // input is indexed from 1

	// Histogram settings
	const ParameterPack* bins_pack_ptr = input_pack.findParameterPack("Bins", KeyType::Required);
	bins_.set_bins( *bins_pack_ptr );

	//----- Read time series -----//

	// Get list of files
	std::vector<std::string> data_files;
	FileSystem::readFilesList(time_series_list_, file_col_, data_files);
	if ( data_files.size() != production_phases.size() ) {
		std::stringstream err_ss;
		err_ss << "Error setting up order parameter " << name_ << "\n"
		       << "  Mismatch between number of time series files parsed (" << data_files.size()
		         << ") and number expected (" << production_phases.size() << "\n";
		throw std::runtime_error( err_ss.str() );
	}
	int num_simulations = data_files.size();


#ifdef DEBUG
	std::cout << "OrderParameter: " << name_ << "\n"
	          << "  From files list " << time_series_list_ << "\n"
	          << "  Found " << num_simulations << " file names\n";
	for ( int i=0; i<num_simulations; ++i ) {
		std::cout << "    " << data_files[i] << "\n";
	}
#endif /* DEBUG */


	// Read time series data for the production phase
	time_series_.clear();
	for ( int i=0; i<num_simulations; ++i ) {
		time_series_.push_back( 
			TimeSeries( data_files[i], data_col_, production_phases[i][0], production_phases[i][1], 
			            use_floored_times )
		);
	}
}
