// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#include "FreeEnergyDistribution.hpp"


FreeEnergyDistribution::FreeEnergyDistribution():
	bins_x_()
{}


FreeEnergyDistribution::FreeEnergyDistribution(const Bins& bins_x_in):
	bins_x_(bins_x_in)
{}


FreeEnergyDistribution::FreeEnergyDistribution(const Bins& bins_x_in, const TimeSeries& time_series_x):
	bins_x_(bins_x_in)
{
	// Sort x-values by bin
	int num_bins_x = bins_x_.getNumBins();
	sample_counts_.assign(num_bins_x, 0);
	int bin;
	int num_samples = time_series_x.size();
	for ( int i=0; i<num_samples; ++i ) {
		bin = bins_x_.findBin( time_series_x[i] );
		if ( bin >= 0 ) {
			++( sample_counts_[bin] );
		}
	}

	// Compute
	f_x_.resize(num_bins_x);
	p_x_.resize(num_bins_x);
	double bin_size_x = bins_x_.getBinSize();
	double normalization = 1.0/(num_samples*bin_size_x);
	for ( int b=0; b<num_bins_x; ++b ) {
		if ( sample_counts_[b] > 0 ) {
			p_x_[b] = sample_counts_[b]*normalization;
			f_x_[b] = -log( p_x_[b] );
		}
		else {
			f_x_[b] = 0.0;
			p_x_[b] = 0.0;
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

	double bin_size_x = ref.bins_x_.getBinSize();
	int num_bins_x = ref.bins_x_.getNumBins();

	for ( int b=0; b<num_bins_x; ++b ) {
		const auto p_ref  = ref.p_x_[b];
		const auto f_ref  = ref.f_x_[b];
		const auto p_dist = dist.p_x_[b];
		const auto f_dist = dist.f_x_[b];

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
				distributions[k].f_x_, distributions[k].getSampleCounts() );
		}
		else {
			f_x_to_print[k] = distributions[k].f_x_;
		}
	}

	const auto& bins = distributions.front().getBins();
	int num_bins = bins.getNumBins();  // All simulations are binned the same way
	for ( int b=0; b<num_bins; ++b ) {
		ofs << bins[b];
		for ( int k=0; k<num_distributions; ++k ) {
			ofs << "\t";
			ofs << std::setw(8) << std::setprecision(5);
			printFreeEnergyValue(ofs, f_x_to_print[k][b], distributions[k].sample_counts_[b]);
		}
		ofs << "\n";
	}
	ofs.close();
}



std::ostream& operator<<(std::ostream& os, const FreeEnergyDistribution& f)
{

	// For convenience of visualizing output, shift F(x) so that F=0 at the minimum
	const auto& f_x           = f.get_F_x();
	const auto& sample_counts = f.getSampleCounts();
	const auto f_x_shifted = FreeEnergyDistribution::shift_f_x_to_zero(f_x, sample_counts);

	// Print F_0(x)
	const auto& bins = f.getBins();
	const int num_bins = bins.getNumBins();
	for ( int b=0; b<num_bins; ++b ) {
		os << std::setw(8) << std::setprecision(5) << bins[b];
		os << "  " << std::setw(8) << std::setprecision(5);
		FreeEnergyDistribution::printFreeEnergyValue(os, f_x_shifted[b], sample_counts[b]);
		os << "  " << std::setw(8) << std::setprecision(5) << sample_counts[b];
		if ( f.hasErrors() ) {
			os << "  " << std::setw(8) << std::setprecision(5) << f.get_err_F_x()[b];
		}
		os << "\n";
	}

	return os;
}