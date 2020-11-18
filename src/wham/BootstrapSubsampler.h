// BootstrapSubsampler
// - Use to generate a random set of indices corresponding to a set of samples,
//   for use with bootstrap methods
//
// TODO 
// - Methods to reseed RNG, change num_samples, ...
// - Default constructor
//   - Use a random sequence from Random
// - Convert Engine to template parameter

#ifndef BOOTSTRAP_SUBSAMPLER_H
#define BOOTSTRAP_SUBSAMPLER_H

#include <algorithm>
#include <functional>
#include <memory>
#include <random>

#include "Random.h"
#include "Assert.hpp"

class BootstrapSubsampler
{
 public:
	using Engine       = std::mt19937;
	using Distribution = std::uniform_int_distribution<int>;
	using SeedSequence = std::seed_seq;

	BootstrapSubsampler(
		const int num_samples,          // number of values to generate
		const std::vector<int>& seeds
	);

	// Produces a random set of indices in the range [0, num_samples_-1], which can be
	// used to create a random sample from a vector of length num_samples_
	void generate(std::vector<int>& samples);

 private:
	int num_samples_;
	std::vector<int> seeds_;

	// Since std::seed_seq is not copyable, need to store in a unique_ptr to enable reseeding
	std::unique_ptr<SeedSequence> seed_sequence_ptr_ = nullptr;

	// TODO possible to allocate statically?
	std::unique_ptr<Engine> mt_engine_ptr_ = nullptr;

	// In order for the number of samples to be able to change, the Distribution
	// must be dynamically allocaed
	std::unique_ptr<Distribution> int_distribution_ptr_ = nullptr;
};

#endif // ifndef BOOTSTRAP_SUBSAMPLER_H
