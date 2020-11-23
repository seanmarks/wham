// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#ifndef OUTPUT_VECTOR_WITH_BINS_HPP
#define OUTPUT_VECTOR_WITH_BINS_HPP

#include "Assert.hpp"
#include "Bins.h"

#include <exception>
#include <stdexcept>


// Organizes a set of output values which are computed over a set of Bins
// - Includes the number of samples per bin, and optionally also errors
class OutputVectorWithBins
{
 public:
  template<typename T>
  using Vector = std::vector<T>;

  using value_type = double;


  //----- Setup -----//

  // Constructs empty without errors
  OutputVectorWithBins(const Bins& bins);

  // Constructs without errors
  OutputVectorWithBins(
    const Bins& bins,
    const Vector<double>& values,
    const Vector<int>& sample_counts
  );
  
  // Constructs with errors
  OutputVectorWithBins(
    const Bins& bins,
    const Vector<double>& values,
    const Vector<double>& errors,
    const Vector<int>& sample_counts
  );

  virtual ~OutputVectorWithBins() = default;

  // Sets the main contents
  void set(const Vector<double>& values, const Vector<int>& sample_counts);

  // Sets the errors associated with the values
  void setErrors(const Vector<double>& errors) {
    const int num_err = errors.size();
    FANCY_ASSERT( num_err == bins_.getNumBins(), "inconsistent input" );
    errors_ = errors;
    has_errors_ = true;
  }


  //----- Get contents -----//

  const std::size_t size() const noexcept {
    return bins_.getNumBins();
  }

  const Bins& getBins() const noexcept {
    return bins_;
  }

  const Vector<double>& getValues() const noexcept {
    return values_;
  }
  
  const Vector<int>& getSampleCounts() const noexcept {
    return sample_counts_;
  }

  // Returns true if the instance has stored errors
  bool hasErrors() const noexcept {
    return has_errors_;
  }

  // Throws if errors are not present
  const Vector<double>& getErrors() const {
    FANCY_ASSERT( hasErrors(), "invalid use of method: errors are not present" );
    return errors_;
  }

  // Returns the 'i'th value
  // - TODO: instead return a ValueWithError?
  const double& operator[](const int i) const NOEXCEPT_IF_NDEBUG {
    FANCY_DEBUG_ASSERT( i >= 0 && i < static_cast<int>(size()), "index out of range: " << i );
    return values_[i];
  }


  //----- Output -----//

  // Saves the results to file
  // TODO: move to separate class? generalize?
  void save(const std::string& file_name) const;


 protected:
  // Returns the name of the variable used for bins, which is used when writing output
  virtual
  std::string getNameOfBins() const {
    return "bin";
  }

  // Returns the name of the vector, which is used when writing output
  virtual
  std::string getNameOfVector() const {
    return "value";
  }


 private:
  Bins bins_;
  Vector<double> values_;
  Vector<int>    sample_counts_;

  bool has_errors_ = false;
  Vector<double> errors_;
};

#endif // ifndef OUTPUT_VECTOR_WITH_BINS_HPP