#include "OrderParameter.h"

OrderParameter::OrderParameter(
		const ParameterPack& input_pack, const std::vector<Simulation>& simulations,
		const bool use_floored_times):
	file_col_(0),
	simulations_(simulations),
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
		if ( file_col_ < 1 ) {
			std::stringstream err_ss;
			err_ss << "Error setting up order parameter " << name_ << "\n"
			       << "  Column number for data must be a positive integer\n";
			throw std::runtime_error( err_ss.str() );
		}
		--file_col_;  // input is indexed from 1
	}

	input_pack.readNumber("DataColumn", KeyType::Required, data_col_);
	if ( data_col_ < 1 ) {
		std::stringstream err_ss;
		err_ss << "Error setting up order parameter " << name_ << "\n"
		       << "  DataColumn must be a positive integer\n";
		throw std::runtime_error( err_ss.str() );
	}
	--data_col_;  // input is indexed from 1

	// Histogram settings
	const ParameterPack* bins_pack_ptr = input_pack.findParameterPack("Bins", KeyType::Required);
	bins_.set_bins( *bins_pack_ptr );

	//----- Read time series -----//

	// Get list of files
	std::vector<std::string> data_files;
	FileSystem::readFilesList(time_series_list_, file_col_, data_files);
	if ( data_files.size() != simulations.size() ) {
		std::stringstream err_ss;
		err_ss << "Error setting up order parameter " << name_ << "\n"
		       << "  Mismatch between number of time series files parsed (" << data_files.size()
		         << ") and number expected (" << simulations.size() << "\n";
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

	time_series_.clear();
	biased_distributions_.clear();
	for ( int i=0; i<num_simulations; ++i ) {
		// Read time series data for the production phase
		time_series_.push_back( 
			TimeSeries( data_files[i], data_col_, simulations_[i].get_t_min(), simulations_[i].get_t_max(), 
			            use_floored_times )
		);

		// Make raw biased distributions
		biased_distributions_.emplace_back( Distribution(bins_, time_series_.back()) );
	}

	// Number of samples in each bin, across all data sets
	int num_bins = bins_.get_num_bins();
	global_sample_counts_.assign(num_bins, 0);
	for ( int j=0; j<num_simulations; ++j ) {
		for ( int b=0; b<num_bins; ++b ) {
			global_sample_counts_[b] += biased_distributions_[j].sample_counts[b];
		}
	}

	// Reserve memory for later
	unbiased_distributions_.resize(num_simulations);
	rebiased_distributions_.resize(num_simulations);
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


void OrderParameter::printRawDistributions() const
{
	int num_simulations = time_series_.size();
	if ( num_simulations < 1 ) {
		throw std::runtime_error("OrderParameter::printRawDistributions: No data found.\n");
	}

	// Common header
	std::stringstream table_header_stream;
	table_header_stream << "# Data sets (by column)\n";
	for ( int i=0; i<num_simulations; ++i ) {
		table_header_stream << "# " << i+2 << ": " << simulations_[i].get_data_set_label() << "\n";
	}
	table_header_stream << "#\n"
	                    << "# " << name_ << " | F(" << name_ << ") [kBT]\n";

	// Working variables
	std::stringstream header_stream;
	std::string file_name;
	std::ofstream ofs;

	// Print biased free energy distributions
	file_name = "F_" + name_ + "_biased.out";
	header_stream.str("");  header_stream.clear();
	header_stream << "# Biased free energy distributions: "
                << " F_i(" << name_ << ") [k_B*T]\n";
	              table_header_stream.str();
	printDistributions( biased_distributions_, file_name, header_stream.str() );


	// Print unbiased free energy distributions (non-consensus)
	file_name = "F_" + name_ + "_unbiased.out";
	header_stream.str("");  header_stream.clear();
	header_stream << "# Unbiased free energy distributions: "
                << " F_{0,i}(" << name_ << ") [k_B*T]\n"
	              << table_header_stream.str();
	printDistributions( unbiased_distributions_, file_name, header_stream.str() );
}


void OrderParameter::printDistributions(
	const std::vector<Distribution>& distributions,
	const std::string& file_name, const std::string& header, const bool shift_to_zero
) const
{
	std::ofstream ofs( file_name );
	ofs << header;

	int num_distributions = distributions.size();
	std::vector<std::vector<double>> f_x_to_print(num_distributions);
	for ( int k=0; k<num_distributions; ++k ) {
		if ( shift_to_zero ) {
			// Shift distributions so that F=0 at the minimum
			f_x_to_print[k] = Distribution::shift_f_x_to_zero(
			                     distributions[k].f_x, distributions[k].sample_counts );
		}
		else {
			f_x_to_print[k] = distributions[k].f_x;
		}
	}

	int num_bins = bins_.get_num_bins();  // All simulations are binned the same way
	for ( int b=0; b<num_bins; ++b ) {
		ofs << bins_[b];
		for ( int k=0; k<num_distributions; ++k ) {
			ofs << "\t";
			ofs << std::setw(8) << std::setprecision(5);
			  Distribution::print_free_energy(ofs, f_x_to_print[k][b], distributions[k].sample_counts[b]);
		}
		ofs << "\n";
	}
	ofs.close();
}
