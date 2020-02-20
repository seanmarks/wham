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
		// - TODO: Make optional?
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
	//unbiased_distributions_.resize(num_simulations);
	//rebiased_distributions_.resize(num_simulations);
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


void OrderParameter::printRebiasedDistributions(std::string file_name) const
{
	// TODO: consistency check for presence of output
	const int num_simulations = simulations_.size();

	if ( file_name.empty() ) {
		file_name = "F_" + name_ + "_rebiased.out";
	}
	std::ofstream ofs(file_name);

	// Header
	std::stringstream header_stream;
	header_stream << "# \"Rebiased\" free energy distributions: "
                << " F_{rebias,i}(" << name_ << ") [k_B*T]\n";
	header_stream << "# Data sets (by column)\n";
	for ( int j=0; j<num_simulations; ++j ) {
		header_stream << "# " << j+2 << ": " << simulations_[j].get_data_set_label() << "\n";
	}
	header_stream << "#\n"
	              << "# " << name_ << " | F(" << name_ << ") [kBT]\n";

	// Print
	printDistributions( rebiased_distributions_, file_name, header_stream.str() );
}


void OrderParameter::printWhamResults(std::string file_name) const
{
	// TODO: consistency check for presence of output
	//const int num_simulations = simulations_.size();

	if ( file_name.empty() ) {
		file_name = "F_" + name_ + "_WHAM.out";
	}
	std::ofstream ofs(file_name);

	// For convenience of visualizing output, shift F(x) so that F=0 at the minimum
	const auto& f_x           = wham_distribution_.f_x;
	const auto& sample_counts = wham_distribution_.sample_counts;
	auto f_x_shifted = Distribution::shift_f_x_to_zero(f_x, sample_counts);

	// Header
	ofs << "# Consensus free energy distributions from WHAM: \n"
      << "#   F(" << name_ << ") [in k_B*T] with T = " << simulations_[0].get_temperature() << " K\n";  // FIXME temperature
	ofs << "# " << name_ << "\tF[kBT]  NumSamples\n";  //"\t" << "\tstderr(F)\n"; TODO error estimate

	// Print F_0(x)
	const int num_bins  = bins_.get_num_bins();
	for ( int b=0; b<num_bins; ++b ) {
		ofs << std::setw(8) << std::setprecision(5) << bins_[b] << "\t";
		ofs << std::setw(8) << std::setprecision(5);
			Distribution::print_free_energy(ofs, f_x_shifted[b], sample_counts[b]);
		ofs << std::setw(8) << std::setprecision(5) << sample_counts[b];
		//<< std::setw(8) << std::setprecision(5) << wham_results_.error_f[b]
		ofs << "\n";
	}
	ofs.close(); ofs.clear();
}



void OrderParameter::printStats(std::string file_name) const
{
	if ( file_name.empty() ) {
		file_name = "stats_" + name_ + ".out";
	}

	// Check for presence of info entropy (TODO: private flag?)
	int  num_simulations   = simulations_.size();
	int  num_info_entropy  = info_entropy_.size();
	bool have_info_entropy = ( num_info_entropy == num_simulations );

	std::ofstream ofs(file_name);

	// Header
	ofs << "# data_set   avg(" << name_ << ")   var(" << name_ << ")";
	if ( have_info_entropy ) {
		ofs << "   info_entropy(biased/rebiased)";
	}
	ofs << "\n";

	// Body
	for ( int j=0; j<num_simulations; ++j ) {
		ofs << simulations_[j].get_data_set_label() << "\t"
		    << time_series_ptrs_[j]->average() << "\t"
		    << time_series_ptrs_[j]->variance();
		if ( have_info_entropy ) {
 			ofs << "\t" << info_entropy_[j] << "\n";
		}
	}
	ofs.close();
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
