#include "OrderParameter.h"

OrderParameter::OrderParameter(
	const std::string& name,
	const ParameterPack& input_pack,
	std::vector<Simulation>& simulations
):
	name_(name),
	simulations_(simulations),
	bins_()
{
	using KeyType = ParameterPack::KeyType;

	// Histogram settings
	const ParameterPack* bins_pack_ptr = input_pack.findParameterPack("Bins", KeyType::Required);
	bins_.set_bins( *bins_pack_ptr );

	int num_simulations = simulations_.size();
	time_series_ptrs_.resize(num_simulations);
	biased_distributions_.clear();
	for ( int i=0; i<num_simulations; ++i ) {
		// Share time series read by Simulation objects
		time_series_ptrs_[i] = simulations_[i].copy_time_series_ptr(name_);
		
		// Make raw biased distributions
		// - TODO: Move to driver?
		biased_distributions_.emplace_back( Distribution(bins_, *(time_series_ptrs_[i])) );
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


void OrderParameter::printRawDistributions() const
{
	int num_simulations = time_series_ptrs_.size();
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
