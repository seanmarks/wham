// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#ifndef NUMERIC_MATRIX_HPP
#define NUMERIC_MATRIX_HPP

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <type_traits>
#include <vector>

#include "Assert.hpp"
#include "OpenMP.h"


namespace numeric {

// Forward-declare friend functions
template<typename T, typename V> class Matrix;
template<typename T, typename V> Matrix<T,V> log(const Matrix<T,V>& other);


// Matrix: A flexible 2-dimensional array
// - Implemented using a 1-dimensional array (default: std::vector)
template<typename T, class Vector = std::vector<T>>
class Matrix
{
 public:
  static_assert(std::is_same<T, typename Vector::value_type>::value, "type mismatch");


  //----- Typedefs and Constants -----//

  static constexpr int N_DIM = 2;
  using Int2 = std::array<int,N_DIM>;
  
  using value_type = T;
  using size_type  = std::size_t;

  using reference       = value_type&;
  using const_reference = const value_type&;

  using pointer         = value_type*;
  using const_pointer   = const value_type*;
  using difference_type = std::ptrdiff_t;

  using iterator               = typename Vector::iterator;  //value_type*;
  using const_iterator         = typename Vector::const_iterator;  //const value_type*;
  using reverse_iterator       = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;


  //----- Constructors -----//

	// Empty
  Matrix() {
    this->clear();
  }

	// With size
  Matrix(const int num_rows, const int num_cols) {
    resize(num_rows, num_cols);
  }

	// With size
  Matrix(const Int2& shape) {
    resize(shape);
  }

	// Fill with constant value
  Matrix(const Int2& shape, const T& value) {
    setShape(shape);
    assign(value);
  }

  // Read from a file
  Matrix(const std::string& file_name);

  // Read from file (named constructor)
  static Matrix FromFile(const std::string& file_name) {
    return Matrix(file_name);
  }


  //----- Size Management -----//

  // Set size
  void resize(const int num_rows, const int num_cols) {
    FANCY_DEBUG_ASSERT(num_rows >= 0, "invalid num rows: " << num_rows);
    FANCY_DEBUG_ASSERT(num_cols >= 0, "invalid num cols: " << num_cols);

    // Allocate memory
    int len = num_rows*num_cols;
    data_.resize(len);

    num_rows_ = num_rows;
    num_cols_ = num_cols;
  }

	// Sets dimensions
	// - Does not preserve stored data
  void resize(const Int2& shape) {
    resize(shape[ROW], shape[COL]);
  }

	// Sets dimensions
	// - Does not preserve stored data
  void setShape(const int num_rows, const int num_cols) {
    resize(num_rows, num_cols);
  }

	// Sets dimensions
	// - Does not preserve stored data
  void setShape(const Int2& shape) {
    resize(shape);
  }

  // Get sizes
  Int2 getShape() const {
    return {{ num_rows_, num_cols_ }}; 
  }

  int getNumRows() const {
    return num_rows_;
  }

  int getNumCols() const {
    return num_cols_;
  }


  //----- Data Management -----//

  // Assign all entries to a particular value
  void assign(const T& value) {
    data_.assign( data_.size(), value );
  }

  void assign(const Int2& shape, const T& value) {
    this->resize(shape);
    this->assign(value);
  }

  void zero() {
    this->assign(0.0);
  }

  // Clears all contents and sets the number of rows
	// and columns to zero
  void clear() {
    data_.clear();
    num_rows_ = 0;
    num_cols_ = 0;
  }


  //----- Copying -----//

  // Copy assignment (from an array of a different type)
  template<typename U, typename VU>
  Matrix& operator=(const Matrix<U,VU>& other);


  //----- Data access -----//

  // Individual elements
  T&       operator()(const int i, const int j);
  const T& operator()(const int i, const int j) const;

  T& operator()(const Int2& indices) {
    return (*this)(indices[ROW], indices[COL]);
  }
  const T& operator()(const Int2& indices) const {
    return (*this)(indices[ROW], indices[COL]);
  }

  // Iterate over all elements in the matrix
  iterator begin() noexcept {
    return data_.begin();
  }
  iterator end() noexcept {
    return data_.end();
  }
  const_iterator begin() const noexcept {
    return data_.cbegin();
  }
  const_iterator end() const noexcept {
    return data_.cend();
  }
  const_iterator cbegin() const noexcept {
    return data_.cbegin();
  }
  const_iterator cend() const noexcept {
    return data_.cend();
  }
  

  // Iterate over a single row
  using RowIterator      = iterator;
  using ConstRowIterator = const_iterator;
  RowIterator begin(const std::size_t i) {
    return &( (*this)(i, 0) );
  }
  RowIterator end(const std::size_t i) {
    return &( (*this)(i+1, 0) );
  }

  // Underlying 1D array (use with caution!)
  Vector& data() noexcept {
    return data_;
  }
  const Vector& data() const noexcept {
    return data_;
  }


  //----- Array Properties -----//

  // Sum of all elements
  T sum() const {
    int len = this->data_.size();
    double s = 0.0;

    #pragma omp parallel for \
      default(shared) schedule(static,10) reduction(+:s)
    for ( int i=0; i<len; ++i ) {
      s += data_[i];
    }

    return s;
  }


  //----- Arithmetic Operators -----//

  // Between two arrays
  Matrix& operator+=(const Matrix& other);

  // Array and scalar
  template<typename U> Matrix& operator+=(const U& value);
  template<typename U> Matrix& operator*=(const U& value);


  //----- Misc. -----//

  // Saves the matrix to given file
  void save(
    const std::string& file_name,
    const std::string& header = ""
  ) const;


  //----- Friend Functions -----//

  template<typename U, typename V>
  friend Matrix<U,V> log(const Matrix<U,V>& other);


 protected:
  static constexpr int ROW = 0;
  static constexpr int COL = 1;

  // Map from 2D indices to 1D index of underlying array
  int getLinearIndex(const int i, const int j) const noexcept;
  int getLinearIndex(const Int2& indices)      const noexcept;

  // Adds an empty row
  void addRow() {
    data_.resize( data_.size() + num_cols_ );
    ++num_rows_;
  }

  // Adds a row of values (number of values must match number of columns)
  void addRow(const Vector& values) {
    int num_values = values.size();
    FANCY_ASSERT( num_values == num_cols_, "size mismatch" );

    data_.insert( data_.end(), values.begin(), values.end() );
    ++num_rows_;
  }

  // Tokenize a line into values
  static Vector parseLine(const std::string& line) {
    Vector values;
    std::stringstream ss(line);

    T value;
    while ( ss >> value ) {
      values.push_back(value);
    }
    
    return values;
  }

 private:
  Vector data_;  // underlying 1D array
  int    num_rows_ = 0;
  int    num_cols_ = 0;
};



template<typename T, typename V>
Matrix<T,V>::Matrix(const std::string& file_name)
{
  std::ifstream ifs(file_name);
  FANCY_ASSERT( ifs.is_open(), "unable to open file: " << file_name);

  // First, determine the number of columns
  int num_cols = 0;
  std::string line;
  while ( std::getline(ifs, line) ) {
    if ( line.empty() || line[0] == '#' ) {
      continue;
    }

    auto values = parseLine(line);
    num_cols = values.size();
    break;
  }

  // Resize as appropriate
  if ( num_cols > 0 ) {
    resize(0, num_cols);
  }
  else {
    return;  // empty matrix
  }

  // Reset stream
  ifs.clear();
  ifs.seekg(0, std::ios::beg);

  // Read data
  while ( std::getline(ifs, line) ) {
    if ( line.empty() || line[0] == '#' ) {
      continue;
    }

    auto values = parseLine(line);
    addRow(values);
  }
}


// Copy assignment (from an array of a different type)
template<typename T, typename VT>
template<typename U, typename VU>
Matrix<T,VT>& Matrix<T,VT>::operator=(const Matrix<U,VU>& other)
{
  this->setShape( other.getShape() );

  int len = data_.size();
  const auto& other_data = other.data();
  for ( int i=0; i<len; ++i ) {
    this->data_[i] = other_data[i];
  }

  return *this;
}


template<typename T, typename V>
inline
T& Matrix<T,V>::operator()(const int i, const int j)
{
  FANCY_DEBUG_ASSERT( getLinearIndex(i,j) < static_cast<int>(data_.size()),
                      "indices (" << i << "," << j << ") are out of bounds "
                      << "(" << num_rows_ << "," << num_cols_ << ")" );
  return data_[ getLinearIndex(i,j) ];
  //return data_[ i*num_cols_ + j ];
}


template<typename T, typename V>
inline
const T& Matrix<T,V>::operator()(const int i, const int j) const
{
  FANCY_DEBUG_ASSERT( getLinearIndex(i,j) < static_cast<int>(data_.size()),
                      "indices (" << i << "," << j << ") are out of bounds "
                      << "(" << num_rows_ << "," << num_cols_ << ")" );
  return data_[ getLinearIndex(i,j) ];
  //return data_[ i*num_cols_ + j ];
}


template<typename T, typename V>
inline
int Matrix<T,V>::getLinearIndex(const int i, const int j) const noexcept
{
  return i*num_cols_ + j;
}


template<typename T, typename V>
inline
int Matrix<T,V>::getLinearIndex(const Int2& indices) const noexcept
{
  return getLinearIndex( indices[ROW], indices[COL] );
}


// Save to file
template<typename T, typename V>
void Matrix<T,V>::save(
  const std::string& file_name, const std::string& header
) const
{
  std::ofstream ofs(file_name);
  FANCY_ASSERT( ofs, "error attempting to write to file: " << file_name );

  if ( ! header.empty() )  {
    ofs << "# " << header << "\n";
  }

  for ( int i=0; i<num_rows_; ++i ) {
    for ( int j=0; j<num_cols_; ++j ) {
      if ( j > 0 ) {
        ofs << " ";
      }
      ofs << (*this)(i,j);
    }
    ofs << "\n";
  }

  ofs.close();
}


// Between two arrays
template<typename T, typename V>
inline
Matrix<T,V>& Matrix<T,V>::operator+=(const Matrix& other)
{
  // TODO: DEBUG MODE: check dimensions
  int len = this->data_.size();
  #pragma omp parallel for \
    default(shared) schedule(static,10)
  for ( int i=0; i<len; ++i ) {
    this->data_[i] += other.data_[i];
  }
  return *this;
}


// Array and scalar
template<typename T, typename V>
template<typename U>
inline
Matrix<T,V>& Matrix<T,V>::operator+=(const U& value) {
  int len = this->data_.size();
  #pragma omp parallel for \
    default(shared) schedule(static,10)
  for ( int i=0; i<len; ++i ) {
    this->data_[i] += value;
  }
  return *this;
}


template<typename T, typename V>
template<typename U>
inline
Matrix<T,V>& Matrix<T,V>::operator*=(const U& value) {
  int len = this->data_.size();
  #pragma omp parallel for \
    default(shared) schedule(static,10)
  for ( int i=0; i<len; ++i ) {
    this->data_[i] *= value;
  }
  return *this;
}



//----- Friend Functions -----//

template<typename T, typename V>
Matrix<T,V> log(const Matrix<T,V>& other)
{
  Matrix<T,V> arr_out( other.getShape() );
  int len = arr_out.data_.size();

  #pragma omp parallel for schedule(static)
  for ( int i=0; i<len; ++i ) {
    arr_out.data_[i] = std::log( other.data_[i] );
  }

  return arr_out;
}

} // end namespace numeric

#endif // ifndef NUMERIC_MATRIX_HPP
