// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "OrderParameter.h"


OrderParameter::OrderParameter(
	const std::string& name,
	const ParameterPack& input_pack
):
	name_(name),
	bins_()
{
	using KeyType = ParameterPack::KeyType;

	// Histogram settings
	const ParameterPack* bins_pack_ptr = input_pack.findParameterPack("Bins", KeyType::Required);
	bins_.setBins( *bins_pack_ptr );
}


void OrderParameter::printWhamResults(std::string file_name, const double temp) const
{
	// TODO: consistency check for presence of output
	//const int num_simulations = simulation_ptrs_.size();

	if ( file_name.empty() ) {
		file_name = "F_" + name_ + "_WHAM.out";
	}
	std::ofstream ofs(file_name);

	// For convenience of visualizing output, shift F(x) so that F=0 at the minimum
	const auto& f_x           = wham_distribution_.get_F_x();
	const auto& sample_counts = wham_distribution_.getSampleCounts();
	auto f_x_shifted = FreeEnergyDistribution::shift_f_x_to_zero(f_x, sample_counts);

	const bool have_error = wham_distribution_.hasErrors();

	// Header
	ofs << "# Consensus free energy distributions from WHAM: \n"
      << "#   F(" << name_ << ") [in k_B*T] with T = " << temp << " K\n";  // FIXME: temperature
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
			ofs << "  " << std::setw(8) << std::setprecision(5) << wham_distribution_.get_err_F_x()[b];
		}
		ofs << "\n";
	}
	ofs.close(); ofs.clear();
}