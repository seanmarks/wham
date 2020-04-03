#include "BootstrapSubsampler.h"

BootstrapSubsampler::BootstrapSubsampler(const int num_samples, const std::vector<int>& seeds):
	num_samples_(num_samples),
	seeds_(seeds),
	seed_sequence_ptr_( new SeedSequence(seeds_.begin(), seeds_.end()) ),
	mt_engine_ptr_( new Engine(*seed_sequence_ptr_) ),
	int_distribution_ptr_( new Distribution(0, num_samples-1) )
{
}

void BootstrapSubsampler::generate(std::vector<int>& samples)
{
	samples.resize(num_samples_);

	// Quick function lambda for generating random numbers
	// - TODO: Private member function?
	auto number_generator = [this]() { return (*int_distribution_ptr_)(*mt_engine_ptr_); };

	std::generate( samples.begin(), samples.end(), number_generator );
}
