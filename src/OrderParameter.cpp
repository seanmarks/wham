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


void OrderParameter::checkForConsistency(const std::vector<OrderParameter>& ops)
{
	int num_ops = ops.size();
	if ( num_ops < 1 ) {
		throw std::runtime_error("there are no order parameters to check for consistency");
	}

	// Use the first OP as a reference
	const OrderParameter& ref_op = ops[0];
	int num_simulations = ref_op.time_series_.size();
	if ( num_simulations < 1 ) {
		throw std::runtime_error("order parameter \"" + ref_op.name_ + "\" has no time series data");
	}
	for ( int j=0; j<num_simulations; ++j ) {
		// Ensure each time series contains data
		if ( ref_op.time_series_[j].size() < 1 ) {
			std::stringstream err_ss;
			err_ss << "Order parameter " << ref_op.name_ << ": time series " << j+1 << " contains no data\n"
			       << "  file = " << ref_op.time_series_[j].get_file() << "\n";
			throw std::runtime_error( err_ss.str() );
		}
	}

	for ( int i=1; i<num_ops; ++i ) {
		// Check number of time series
		if ( ops[i].time_series_.size() != ref_op.time_series_.size() ) {
			std::stringstream err_ss;
			err_ss << "Order parameters " << ref_op.name_ << " and " << ops[i].name_
			       << " have different numbers of time series.\n";
		}

		// Check each time series
		for ( int j=0; j<num_simulations; ++j ) {
			// Convenient aliases
			const TimeSeries& series_i   = ops[i].time_series_[j];
			const TimeSeries& series_ref = ref_op.time_series_[j];

			// Time series from the same simulation should have the same number of points
			if ( series_i.size() != series_ref.size() ) {
				std::stringstream err_ss;
				err_ss << "Order parameters: time series length mismatch for simulation " 
				       << j+1 << " of " << num_simulations << "\n";

				std::vector<int> op_indices = {{ 0, i }};
				for ( unsigned k=0; k<op_indices.size(); ++k ) {
					int l = op_indices[k];
					const TimeSeries& time_series = ops[l].time_series_[j];
					err_ss << "  OrderParameter = " << ops[l].name_ << ": " << time_series.size() << " points stored\n"
					       << "    file = " << time_series.get_file() << "\n";
				}
				throw std::runtime_error( err_ss.str() );
			}

			// These time series should also refer to the same times sampled
			if ( series_i.get_times() != series_ref.get_times() ) {
				std::stringstream err_ss;
				err_ss << "Order parameters: stored times for each sample do not match for simulation "
				       << j+1 << " of " << num_simulations << "\n";

				std::vector<int> op_indices = {{ 0, i }};
				for ( unsigned k=0; k<op_indices.size(); ++k ) {
					int l = op_indices[k];
					const TimeSeries& time_series = ops[l].time_series_[j];
					err_ss << "  OrderParameter = " << ops[l].name_ << ":\n"
					       << "    file = " << time_series.get_file() << "\n";
				}
				throw std::runtime_error( err_ss.str() );
			}
		}
	}
}
