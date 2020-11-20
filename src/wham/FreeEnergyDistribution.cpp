// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#include "FreeEnergyDistribution.hpp"


FreeEnergyDistribution::FreeEnergyDistribution():
	bins_x()
{}


FreeEnergyDistribution::FreeEnergyDistribution(const Bins& bins_x_in):
	bins_x(bins_x_in)
{}


FreeEnergyDistribution::FreeEnergyDistribution(const Bins& bins_x_in, const TimeSeries& time_series_x):
	bins_x(bins_x_in)
{
	// Sort x-values by bin
	int num_bins_x = bins_x.getNumBins();
	sample_counts.assign(num_bins_x, 0);
	int bin;
	int num_samples = time_series_x.size();
	for ( int i=0; i<num_samples; ++i ) {
		bin = bins_x.findBin( time_series_x[i] );
		if ( bin >= 0 ) {
			++( sample_counts[bin] );
		}
	}

	// Compute
	f_x.resize(num_bins_x);
	p_x.resize(num_bins_x);
	double bin_size_x = bins_x.getBinSize();
	double normalization = 1.0/(num_samples*bin_size_x);
	for ( int b=0; b<num_bins_x; ++b ) {
		if ( sample_counts[b] > 0 ) {
			p_x[b] = sample_counts[b]*normalization;
			f_x[b] = -log( p_x[b] );
		}
		else {
			f_x[b] = 0.0;
			p_x[b] = 0.0;
		}
	}
}


void FreeEnergyDistribution::print(const std::string& file, const std::string& header) const
{
	// TODO:
	throw std::runtime_error("Distribution::print is not yet implemented");
}


double FreeEnergyDistribution::computeInformationEntropy(const FreeEnergyDistribution& dist, const FreeEnergyDistribution& ref)
{
	// TODO: consistency checks
	double info_entropy = 0.0;

	double bin_size_x = ref.bins_x.getBinSize();
	int num_bins_x = ref.bins_x.getNumBins();

	for ( int b=0; b<num_bins_x; ++b ) {
		const auto p_ref  = ref.p_x[b];
		const auto f_ref  = ref.f_x[b];
		const auto p_dist = dist.p_x[b];
		const auto f_dist = dist.f_x[b];

		if ( p_ref > 0.0 and p_dist > 0.0 ) {
			info_entropy += p_ref*(f_dist - f_ref)*bin_size_x;
		}
	}

	return info_entropy;
}



void FreeEnergyDistribution::printSet(
	const std::vector<FreeEnergyDistribution>& distributions,
	const std::string& file_name, const std::string& header, const bool shift_to_zero
)
{
	std::ofstream ofs( file_name );
	ofs << header;

	int num_distributions = distributions.size();
	std::vector<std::vector<double>> f_x_to_print(num_distributions);
	for ( int k=0; k<num_distributions; ++k ) {
		if ( shift_to_zero ) {
			// Shift distributions so that F=0 at the minimum
			f_x_to_print[k] = shift_f_x_to_zero(
				distributions[k].f_x, distributions[k].getSampleCounts() );
		}
		else {
			f_x_to_print[k] = distributions[k].f_x;
		}
	}

	const auto& bins = distributions.front().getBins();
	int num_bins = bins.getNumBins();  // All simulations are binned the same way
	for ( int b=0; b<num_bins; ++b ) {
		ofs << bins[b];
		for ( int k=0; k<num_distributions; ++k ) {
			ofs << "\t";
			ofs << std::setw(8) << std::setprecision(5);
			printFreeEnergyValue(ofs, f_x_to_print[k][b], distributions[k].sample_counts[b]);
		}
		ofs << "\n";
	}
	ofs.close();
}