// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "OutputVectorWithBins.hpp"


OutputVectorWithBins::OutputVectorWithBins(const Bins& bins):
  bins_(bins),
  values_(bins.size()),
  sample_counts_(bins.size(), 0)
{
}


OutputVectorWithBins::OutputVectorWithBins(
  const Bins& bins, const Vector<double>& values, const Vector<int>& sample_counts
):
  OutputVectorWithBins(bins)
{
  set(values, sample_counts);
}


OutputVectorWithBins::OutputVectorWithBins(
  const Bins& bins, const Vector<double>& values, const Vector<double>& errors,
  const Vector<int>& sample_counts
):
  OutputVectorWithBins(bins, values, sample_counts)
{
  setErrors(errors);
}


void OutputVectorWithBins::set(const Vector<double>& values, const Vector<int>& sample_counts)
{
  const unsigned num_bins = bins_.getNumBins();
  FANCY_ASSERT( values.size() == num_bins,             "inconsistent input" );
  FANCY_ASSERT( values.size() == sample_counts.size(), "inconsistent input" );
  values_        = values;
  sample_counts_ = sample_counts;
}


void OutputVectorWithBins::save(const std::string& file_name) const
{
  FANCY_ASSERT ( ! file_name.empty(), "no file name provided" );

  std::ofstream ofs(file_name);
  FANCY_ASSERT( ! ofs.fail(), "unable to open file: " << file_name );

  // Header
  const std::string spacer("  ");
  ofs << "# " << getNameOfBins() 
      << spacer << getNameOfVector();
  if ( hasErrors() ) {
    ofs << spacer << "+/-err";
  }
  ofs << spacer << "samples"
      << '\n';

  // Data
  const int num_bins = bins_.getNumBins();
  for ( int b=0; b<num_bins; ++b ) {
    ofs << bins_[b] << spacer << values_[b];
    if ( hasErrors() ) {
      ofs << spacer << errors_[b];
    }
    ofs << spacer << sample_counts_[b]
        << '\n';
  }

  ofs.close();
}