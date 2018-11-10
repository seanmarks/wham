#include "Bias.h"

Bias::Bias(const std::vector<std::string>& input_tokens, const double kBT):
	kBT_(kBT), beta_(1.0/kBT_)
{
	// TODO make more sophisticated
	if ( input_tokens.size() != 3 ) {
		throw std::runtime_error("Bias expects 3 tokens (for now)\n");
	}

	// TODO assumes harmonic potential
	std::string type = input_tokens[0];
	double x_star = std::stod( input_tokens[1] );
	double kappa  = std::stod( input_tokens[2] );
	kappa *= beta_;  // convert from kJ/mol to kBT

	std::unique_ptr<Potential> tmp_ptr( new HarmonicPotential(x_star, kappa) );
	potential_ptrs_.push_back( std::move(tmp_ptr) );
}


// Evaluate the bias
double Bias::evaluate(const std::vector<double>& x) const
{
	// Check input
	unsigned num_potentials = potential_ptrs_.size();
	if ( x.size() != num_potentials ) {
		throw std::runtime_error("number of arguments doesn't match number of potentials");
	}

	// Compute the total biasing potential
	double u_bias = 0.0;
	for ( unsigned i=0; i<num_potentials; ++i  ) {
		u_bias += potential_ptrs_[i]->evaluate( x[i] );
	}

	return u_bias;
}
