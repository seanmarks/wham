// Random: namespace with helper routines for dealing with random numbers

#ifndef RANDOM_H
#define RANDOM_H

#include <algorithm>
#include <random>
#include <vector>

#include "Assert.hpp"


namespace Random
{

// Since std::seed_seq is not copyable, instead return a vector of values which can be 
// used to create a std::seed_seq
std::vector<int> generateRandomSequence(const int num_values = 10);

// Returns a fixed sequence (useful for debugging/testing purposes)
std::vector<int> getDebugSequence();

} // end namespace Random

#endif // ifndef RANDOM_H
