// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "Estimator_F_x_y.hpp"


Estimator_F_x_y::Estimator_F_x_y(
  const OrderParameter& x, const OrderParameter& y
):
  x_(x),
  y_(y)
{
}


void Estimator_F_x_y::calculateUsingStoredWeights(
	const Wham& wham
)
{
	// Sort samples by bin
  const auto& bins_x = x_.getBins();
  const auto& bins_y = y_.getBins();
  const int   num_bins_x = bins_x.getNumBins();
  const int   num_bins_y = bins_y.getNumBins();
  const int   num_samples_total = wham.getNumSamplesTotal();
  const int   num_to_reserve = 2*(num_samples_total/(num_bins_x*num_bins_y));
  binned_weights_.resize(num_bins_x, num_bins_y);
  for ( auto& buffer : binned_weights_ ) {
    buffer.resize(0);
    buffer.reserve(num_to_reserve);
  }
  const auto& simulations = wham.getSimulationData();
  const auto& data_ranges = wham.getSimulationDataRanges();
  const int num_simulations = simulations.size();
	const auto& weights = this->getWeights();
	for ( int j=0; j<num_simulations; ++j ) {
		const auto& x_j = simulations[j].getTimeSeriesForOrderParameter(x_.getName());
		const auto& y_j = simulations[j].getTimeSeriesForOrderParameter(y_.getName());

		const int num_samples = x_j.size();
		for ( int i=0; i<num_samples; ++i ) {
			int bx = bins_x.findBin( x_j[i] );
      int by = bins_y.findBin( y_j[i] );
      if ( bx >= 0 && by >= 0 ) {
				int n = data_ranges[j].first + i;
        binned_weights_(bx,by).push_back( weights[n] );
      }
    }
	}

	// Calculate distributions
	const double bin_size_x = bins_x.getBinSize();
	const double bin_size_y = bins_y.getBinSize();
  const double fac = 1.0/(bin_size_x*bin_size_y);  // normalize by bin size
  Matrix<double> p_x_y(num_bins_x, num_bins_y);
  Matrix<int>    samples(num_bins_x, num_bins_y);
  #pragma omp parallel for collapse(2)
  for ( int i=0; i<num_bins_x; ++i ) {
    for ( int j=0; j<num_bins_y; ++j ) {
      const auto& weights = binned_weights_(i,j);
      double sum = std::accumulate( weights.begin(), weights.end(), 0.0 );
      p_x_y(i,j) = fac*sum;
      samples(i,j) = weights.size();
    }
  }

  Matrix<double> f_x_y = numeric::log(p_x_y);
  f_x_y *= -1.0;

  f_x_y_ = FreeEnergyDistribution2D(bins_x, bins_y, f_x_y, p_x_y, samples);
}



void Estimator_F_x_y::saveResults() const
{
  // Working variables
	const Bins& bins_x = x_.getBins();
	const Bins& bins_y = y_.getBins();
	const int num_bins_x = bins_x.getNumBins();
	const int num_bins_y = bins_y.getNumBins();
	std::string file_name, sep;
	std::ofstream ofs;

	// Print bins
	std::vector<const OrderParameter*> op_ptrs = {{ &x_, &y_ }};
	for ( unsigned i=0; i<op_ptrs.size(); ++i ) {
		file_name = "bins_" + op_ptrs[i]->getName() + ".out";
		ofs.open(file_name);

		ofs << "# Bins for " << op_ptrs[i]->getName() << " for F_WHAM(" << x_.getName() << "," << y_.getName() << ")\n"
				<< "# " << op_ptrs[i]->getName() << "\n";

		const auto& bins = op_ptrs[i]->getBins();
		int num_bins = bins.getNumBins();
		for ( int j=0; j<num_bins; ++j ) {
			ofs << bins[j]  << "\n";
		}
		ofs.close(); ofs.clear();
	}

	// Print F(x,y)  (TODO: shift so that F=0 at the minimum)
	file_name = "F_" + x_.getName() + "_" + y_.getName() + "_WHAM.out";
	ofs.open(file_name);
	sep = " ";
  const auto& f_x_y_wham    = f_x_y_.get_f_x_y();
  const auto& sample_counts = f_x_y_.getSampleCounts();
	for ( int i=0; i<num_bins_x; ++i ) {
		for ( int j=0; j<num_bins_y; ++j ) {
			FreeEnergyDistribution::printFreeEnergyValue(ofs, f_x_y_wham(i,j), sample_counts(i,j));
			ofs << sep;
		}
		ofs << "\n";
	}
	ofs.close(); ofs.clear();

	// Print sample distribution
	file_name = "samples_" + x_.getName() + "_" + y_.getName() + ".out";
	ofs.open(file_name);
	sep = " ";
	for ( int i=0; i<num_bins_x; ++i ) {
		for ( int j=0; j<num_bins_y; ++j ) {
			ofs << sample_counts(i,j) << sep;
		}
		ofs << "\n";
	}
	ofs.close(); ofs.clear();
}