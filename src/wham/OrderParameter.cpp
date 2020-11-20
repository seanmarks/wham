// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "OrderParameter.h"


OrderParameter::OrderParameter(
	const std::string& name,
	const ParameterPack& input_pack,
	std::vector<Simulation>& simulations
):
	name_(name),
	simulation_ptrs_(simulations.size(), nullptr),
	bins_()
{
	using KeyType = ParameterPack::KeyType;

	// Histogram settings
	const ParameterPack* bins_pack_ptr = input_pack.findParameterPack("Bins", KeyType::Required);
	bins_.setBins( *bins_pack_ptr );

	setSimulations(simulations);
}


void OrderParameter::setSimulations(std::vector<Simulation>& simulations)
{
	int num_simulations = simulations.size();
	simulation_ptrs_.assign(num_simulations, nullptr);
	time_series_ptrs_.resize(num_simulations);
	biased_distributions_.clear();

	// TODO: Move everything besides setting of simulation_ptrs to a separate fxn

	for ( int i=0; i<num_simulations; ++i ) {
		simulation_ptrs_[i] = &simulations[i];

		// Share time series read by Simulation objects
		time_series_ptrs_[i] = simulations[i].copyTimeSeriesPtr(name_);
		
		// Make raw biased distributions
		// - TODO: Make optional?
		biased_distributions_.emplace_back( FreeEnergyDistribution(bins_, *(time_series_ptrs_[i])) );
	}
}


void OrderParameter::printRawDistributions() const
{
	int num_simulations = time_series_ptrs_.size();
	FANCY_ASSERT(num_simulations > 0, "no simulations present");

	// Common header
	std::stringstream table_header_stream;
	table_header_stream << "# Data sets (by column)\n";
	for ( int i=0; i<num_simulations; ++i ) {
		table_header_stream << "# " << i+2 << ": " << simulation_ptrs_[i]->getDataSetLabel() << "\n";
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
	FreeEnergyDistribution::printSet( biased_distributions_, file_name, header_stream.str() );


	// Print unbiased free energy distributions (non-consensus)
	file_name = "F_" + name_ + "_unbiased.out";
	header_stream.str("");  header_stream.clear();
	header_stream << "# Unbiased free energy distributions: "
                << " F_{0,i}(" << name_ << ") [k_B*T]\n"
	              << table_header_stream.str();
	FreeEnergyDistribution::printSet( unbiased_distributions_, file_name, header_stream.str() );
}


void OrderParameter::printRebiasedDistributions(std::string file_name) const
{
	FANCY_ASSERT(simulation_ptrs_.size() == rebiased_distributions_.size(), "length mismatch");

	if ( file_name.empty() ) {
		file_name = "F_" + name_ + "_rebiased.out";
	}
	std::ofstream ofs(file_name);

	// Header
	std::stringstream header_stream;
	header_stream << "# \"Rebiased\" free energy distributions: "
                << " F_{rebias,i}(" << name_ << ") [k_B*T]\n";
	header_stream << "# Data sets (by column)\n";
	const int num_simulations = simulation_ptrs_.size();
	for ( int j=0; j<num_simulations; ++j ) {
		header_stream << "# " << j+2 << ": " << simulation_ptrs_[j]->getDataSetLabel() << "\n";
	}
	header_stream << "#\n"
	              << "# " << name_ << " | F(" << name_ << ") [kBT]\n";

	// Print
	FreeEnergyDistribution::printSet( rebiased_distributions_, file_name, header_stream.str() );
}


void OrderParameter::printWhamResults(std::string file_name) const
{
	// TODO: consistency check for presence of output
	//const int num_simulations = simulation_ptrs_.size();

	if ( file_name.empty() ) {
		file_name = "F_" + name_ + "_WHAM.out";
	}
	std::ofstream ofs(file_name);

	// For convenience of visualizing output, shift F(x) so that F=0 at the minimum
	const auto& f_x           = wham_distribution_.f_x;
	const auto& sample_counts = wham_distribution_.sample_counts;
	auto f_x_shifted = FreeEnergyDistribution::shift_f_x_to_zero(f_x, sample_counts);

	// TODO: safer check?
	bool have_error = ( wham_distribution_.error_f_x.size() == f_x.size() );

	// Header
	ofs << "# Consensus free energy distributions from WHAM: \n"
      << "#   F(" << name_ << ") [in k_B*T] with T = " << simulation_ptrs_[0]->getTemperature() << " K\n";  // FIXME temperature
	ofs << "# " << name_ << "\tF[kBT]  NumSamples";
	if ( have_error ) {
		ofs << "\tstderr(F)[kBT]";
	}
	ofs << "\n";

	// Print F_0(x)
	const int num_bins  = bins_.getNumBins();
	for ( int b=0; b<num_bins; ++b ) {
		ofs << std::setw(8) << std::setprecision(5) << bins_[b];
		ofs << "  " << std::setw(8) << std::setprecision(5);
			FreeEnergyDistribution::printFreeEnergyValue(ofs, f_x_shifted[b], sample_counts[b]);
		ofs << "  " << std::setw(8) << std::setprecision(5) << sample_counts[b];
		if ( have_error ) {
			ofs << "  " << std::setw(8) << std::setprecision(5) << wham_distribution_.error_f_x[b];
		}
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
	int  num_simulations   = simulation_ptrs_.size();
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
		ofs << simulation_ptrs_[j]->getDataSetLabel() << "\t"
		    << time_series_ptrs_[j]->average() << "\t"
		    << time_series_ptrs_[j]->variance();
		if ( have_info_entropy ) {
 			ofs << "\t" << info_entropy_[j] << "\n";
		}
	}
	ofs.close();
}
