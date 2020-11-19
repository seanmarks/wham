// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

// TODO: Move to another file with more numeric helper functions?

#ifndef LOG_SUM_EXP_HPP
#define LOG_SUM_EXP_HPP

#include <algorithm>
#include <cmath>
#include <vector>

#include "Assert.hpp"

namespace numeric
{

// Returns the logarithm of a sum of exponentials:
//      log(S) = sum_{i=1}^N exp(args[i])
// - input: arguments of exponentials
// - THREAD_SAFE
template<typename T>
T logSumExp(const std::vector<T>& args)
{
	int num_args = static_cast<int>( args.size() );
	FANCY_ASSERT( num_args > 0, "no arguments supplied" );

	// Compute sum of exp(args[i] - max_arg)
	T max_arg = *std::max_element( args.begin(), args.end() );
	T sum = 0.0;
	#pragma omp simd reduction(+: sum)
	for ( int i=0; i<num_args; ++i ) {
		sum += std::exp(args[i] - max_arg);
	}

	// Since exp(max_arg) might be huge, return log(sum_exp) instead of sum_exp
	T log_sum_exp_out = max_arg + std::log(sum);

	return log_sum_exp_out;
}

}

#endif // ifndef LOG_SUM_EXP_HPP