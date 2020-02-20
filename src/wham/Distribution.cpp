#include "Distribution.h"

Distribution::Distribution():
	bins_x()
{
}


Distribution::Distribution(const Bins& bins_x_in):
	bins_x(bins_x_in)
{
}


Distribution::Distribution(const Bins& bins_x_in, const TimeSeries& time_series_x):
	bins_x(bins_x_in)
{
	// Sort x-values by bin
	int num_bins_x = bins_x.get_num_bins();
	sample_counts.assign(num_bins_x, 0);
	int bin;
	int num_samples = time_series_x.size();
	for ( int i=0; i<num_samples; ++i ) {
		bin = bins_x.find_bin( time_series_x[i] );
		if ( bin >= 0 ) {
			++( sample_counts[bin] );
		}
	}

	// Compute
	f_x.resize(num_bins_x);
	p_x.resize(num_bins_x);
	double bin_size_x = bins_x.get_bin_size();
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


void Distribution::print(const std::string& file, const std::string& header) const
{
	// TODO
	throw std::runtime_error("Distribution::print is not yet implemented");
}


double Distribution::computeInformationEntropy(const Distribution& first, const Distribution& second)
{
	// TODO consistency checks
	double info_entropy = 0.0;

	double bin_size_x = first.bins_x.get_bin_size();
	int num_bins_x = first.bins_x.get_num_bins();

	for ( int b=0; b<num_bins_x; ++b ) {
		const double& p_first  = first.p_x[b];
		const double& f_first  = first.f_x[b];
		const double& f_second = second.f_x[b];

		if ( p_first > 0.0 and std::isfinite(f_second) ) {
			info_entropy += p_first*(f_second - f_first)*bin_size_x;
		}
	}

	return info_entropy;
}
