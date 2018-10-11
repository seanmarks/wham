
#ifndef MATRIX_H
#define MATRIX_H

// N-dimensional array indexing:
//  See https://docs.scipy.org/doc/numpy-1.13.0/reference/arrays.ndarray.html --> Indexing


template<typename T, std::size_t N>
class DenseMatrix
{
	// Restrict number of indices to N
	template<class... Indices, std::enable_if_t<(sizeof...(Indices) == N)>* = nullptr>
	T& operator()(const Indices&... indices) {
	}

	std::vector<T> data_;
}
#endif /* MATRIX_H */
