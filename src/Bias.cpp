#include "Bias.h"

Bias::Bias(const std::string& input_line, const double kBT):
	kBT_(kBT), beta_(1.0/kBT_)
{	
	std::vector<std::string> input_tokens;
	StringTools string_tools;
	string_tools.split(input_line, input_tokens);

	auto start = input_tokens.begin();
	auto stop  = input_tokens.end();
	auto last = start;
	StringIt next_first;
	for ( auto first = start; first != stop; first = next_first ) {
		next_first = find_next_potential( first, stop, 
		                                  first, last );
		int index      = std::distance(start, first);
		int num_tokens = std::distance(first, last);

		if ( num_tokens < 2 ) {
			throw std::runtime_error("Error: bias terms must have at least 2 tokens");
		}
		std::string op_name = input_tokens[index++];
		std::string type    = input_tokens[index++];
		int num_tokens_left = num_tokens - 2;

		// Construct this term
		std::unique_ptr<Potential> tmp_ptr;
		if ( type == "harmonic" ) {
			if ( num_tokens_left != 2 ) {
				throw std::runtime_error("invalid harmonic bias");
			}
			double x_star = std::stod( input_tokens[index++] );
			double kappa  = std::stod( input_tokens[index++] );
			kappa *= beta_;  // convert from kJ/mol to kBT
			tmp_ptr = std::unique_ptr<Potential>(new HarmonicPotential(x_star, kappa) );
		}
		else if ( type == "none" or type == "null" or type == "zero" ) {
			if ( num_tokens_left != 0 ) {
				throw std::runtime_error("invalid null/none/zero bias");
			}
			tmp_ptr = std::unique_ptr<Potential>( new ZeroPotential() );
		}
		else {
			throw std::runtime_error("unrecognized potential type \"" + type + "\"");
		}
		/*
		TODO other potentials
		*/

		// Record the term, and the name of the order parameter it affects
		potential_ptrs_.push_back( std::move(tmp_ptr) );
		order_parameter_names_.push_back( op_name );
	}

	// Sanity checks
	if ( potential_ptrs_.size() == 0 ) {
		throw std::runtime_error("no terms added to bias");
	}
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
