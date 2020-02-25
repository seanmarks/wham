#include "Random.h"

namespace Random
{


std::vector<int> generateRandomSequence(const int num_values)
{
	std::random_device rd;

	std::vector<int> sequence(num_values);
	for ( int i=0; i<num_values; ++i ) {
		sequence[i] = rd();
	}

	return sequence;
}


std::vector<int> getDebugSequence()
{
	// Random numbers from atmospheric noise (random.org)
	// TODO: More values
	std::vector<int> debug_sequence = {
		78307901, 6985620
	};

	return debug_sequence;
}


} // end namespace Random
