#include "WhamResults1D.h"

void WhamResults1D::print() const
{
	printConsensusResults();

	printStats();
}


void WhamResults1D::printConsensusResults() const
{
	// Working variables
	std::string file_name;
	std::ofstream ofs;
	const int num_simulations = simulations.size();

	// Unpack for readability below
	const auto& bins_x        = x.get_bins();
	const int num_bins_x      = bins_x.get_num_bins();
	const auto& x_name        = x.get_name();
	const auto& sample_counts = f_x_wham.sample_counts;
	//const std::vector<double>& p_x_wham = x.wham_distribution_.p_x;
	//const std::vector<double>& f_x_wham = x.wham_distribution_.f_x;
	// TODO error

	// For convenience of visualizing output, shift F(x) so that F=0 at the minimum
	auto f_x_wham_shifted = Distribution::shift_f_x_to_zero( f_x_wham.f_x, sample_counts );

	// Print F_0(x)
	file_name = "F_" + x_name + "_WHAM.out";
	ofs.open(file_name);
	ofs << "# Consensus free energy distributions from WHAM: \n"
      << "#   F(" << x_name << ") [in k_B*T] with T = " << simulations[0].get_temperature() << " K\n";  // FIXME temperature
	ofs << "# " << x_name << "\tF[kBT]  NumSamples\n";  //"\t" << "\tstderr(F)\n"; TODO error estimate
	for ( int b=0; b<num_bins_x; ++b ) {
		ofs << std::setw(8) << std::setprecision(5) << bins_x[b] << "\t";
		ofs << std::setw(8) << std::setprecision(5);
			Distribution::print_free_energy(ofs, f_x_wham_shifted[b], sample_counts[b]);
		ofs << std::setw(8) << std::setprecision(5) << sample_counts[b];
		//<< std::setw(8) << std::setprecision(5) << wham_results_.error_f[b]
		ofs << "\n";
	}
	ofs.close(); ofs.clear();

	// "Rebiased" free energy distributions
	file_name = "F_" + x_name + "_rebiased.out";

	std::stringstream header_stream;
	header_stream << "# \"Rebiased\" free energy distributions: "
                << " F_{rebias,i}(" << x_name << ") [k_B*T]\n";
	header_stream << "# Data sets (by column)\n";
	for ( int j=0; j<num_simulations; ++j ) {
		header_stream << "# " << j+2 << ": " << simulations[j].get_data_set_label() << "\n";
	}
	header_stream << "#\n"
	              << "# " << x_name << " | F(" << x_name << ") [kBT]\n";

	printDistributions( f_x_rebiased, file_name, header_stream.str() );
}


void WhamResults1D::printStats(std::string file_name) const
{
	// TODO check for presence of info_entropy_

	const auto& x_name = x.get_name();

	if ( file_name.empty() ) {
		file_name = "stats_" + x_name + ".out";
	}

	std::ofstream ofs(file_name);
	ofs << "# data_set   avg(" << x_name << ")   var(" << x_name << ")   "
	    << "info_entropy(biased/rebiased)\n";
	int num_simulations = simulations.size();
	for ( int j=0; j<num_simulations; ++j ) {
		ofs << simulations[j].get_data_set_label() << "\t"
		    << x.get_time_series(j).average() << "\t"
		    << x.get_time_series(j).variance() << "\t"
		    << info_entropy[j] << "\n";
	}
	ofs.close(); ofs.clear();
}


void WhamResults1D::printDistributions(
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

	const auto bins = x.get_bins();
	const int num_bins = bins.get_num_bins();  // All simulations are binned the same way
	for ( int b=0; b<num_bins; ++b ) {
		ofs << bins[b];
		for ( int k=0; k<num_distributions; ++k ) {
			ofs << "\t";
			ofs << std::setw(8) << std::setprecision(5);
			Distribution::print_free_energy(ofs, f_x_to_print[k][b], distributions[k].sample_counts[b]);
		}
		ofs << "\n";
	}
	ofs.close();
}
