#include "Bias.h"

Bias::Bias(const std::string& input_line, const double kBT):
	kBT_(kBT), beta_(1.0/kBT_)
{	
	std::vector<std::string> input_tokens;
	StringTools string_tools;
	string_tools.split(input_line, input_tokens);

	auto start = input_tokens.begin();
	auto stop  = input_tokens.end();
	auto end  = start;
	StringIt next_first;
	for ( auto first = start; first != stop; first = next_first ) {
		next_first = find_next_potential( first, stop, 
		                                  first, end );
		int index      = std::distance(start, first);
		int num_tokens = std::distance(first, end);

		if ( num_tokens < 2 ) {
			throw std::runtime_error("Error: bias terms must have at least 2 tokens");
		}
		std::string op_name = input_tokens[index++];
		std::string type    = input_tokens[index++];
		auto first_parameter = first + 2;

		// Construct this term
		// TODO use a factory
		// TODO other potentials
		std::unique_ptr<Potential> tmp_ptr;
		if ( type == "harmonic" ) {
			tmp_ptr = std::unique_ptr<Potential>( new HarmonicPotential(first_parameter, end, kBT_) );
		}
		else if ( type == "linear" ) {
			tmp_ptr = std::unique_ptr<Potential>( new LinearPotential(first_parameter, end, kBT_) );
		}
		else if ( type == "left_harmonic" ) {
			tmp_ptr = std::unique_ptr<Potential>( new LeftHarmonicPotential(first_parameter, end, kBT_) );
		}
		else if ( type == "right_harmonic" ) {
			tmp_ptr = std::unique_ptr<Potential>( new RightHarmonicPotential(first_parameter, end, kBT_) );
		}
		else if ( type == "none"  ) {
			tmp_ptr = std::unique_ptr<Potential>( new ZeroPotential(first_parameter, end, kBT_) );
		}
		else {
			throw std::runtime_error("unrecognized potential type \"" + type + "\"");
		}

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


// Harmonic potential
Bias::HarmonicPotential::HarmonicPotential(
		const StringIt& first_parameter, const StringIt& end, const double kBT):
	Potential(kBT)
{
	int num_tokens = std::distance(first_parameter, end);
	if ( num_tokens != 2 ) {
		throw std::runtime_error("invalid harmonic potential");
	}

	auto it = first_parameter;
	x_star_ = std::stod( *it );        ++it;
	kappa_  = beta_*std::stod( *it );  ++it;
}

double Bias::HarmonicPotential::evaluate(const double x) const 
{
	double delta_x = x - x_star_;
	return 0.5*kappa_*delta_x*delta_x;
}


// Linear potential
Bias::LinearPotential::LinearPotential(
		const StringIt& first_parameter, const StringIt& end, const double kBT):
	Potential(kBT)
{
	int num_tokens = std::distance(first_parameter, end);
	if ( num_tokens != 2 ) {
		throw std::runtime_error("invalid linear potential");
	}

	// Convert to kBT
	auto it = first_parameter;
	phi_ = beta_*std::stod( *it );  ++it;
	c_   = beta_*std::stod( *it );  ++it;
}
double Bias::LinearPotential::evaluate(const double x) const 
{
	return phi_*x + c_;
}


// Left one-sided harmonic potential
Bias::LeftHarmonicPotential::LeftHarmonicPotential(
		const StringIt& first_parameter, const StringIt& end, const double kBT):
	Potential(kBT)
{
	int num_tokens = std::distance(first_parameter, end);
	if ( num_tokens != 2 ) {
		throw std::runtime_error("invalid left one-sided harmonic potential");
	}

	auto it = first_parameter;
	x_left_ = std::stod( *it );        ++it;
	k_left_ = beta_*std::stod( *it );  ++it;  // convert to kBT
}
double Bias::LeftHarmonicPotential::evaluate(const double x) const
{
	if ( x < x_left_ ) {
		double delta_x = x - x_left_;
		return 0.5*k_left_*delta_x*delta_x;
	}
	else {
		return 0.0;
	}
}


// Right one-sided harmonic potential
Bias::RightHarmonicPotential::RightHarmonicPotential(
		const StringIt& first_parameter, const StringIt& end, const double kBT):
	Potential(kBT)
{
	int num_tokens = std::distance(first_parameter, end);
	if ( num_tokens != 2 ) {
		throw std::runtime_error("invalid right one-sided harmonic potential");
	}

	auto it = first_parameter;
	x_right_ = std::stod( *it );        ++it;
	k_right_ = beta_*std::stod( *it );  ++it;  // convert to kBT
}
double Bias::RightHarmonicPotential::evaluate(const double x) const
{
	if ( x > x_right_ ) {
		double delta_x = x - x_right_;
		return 0.5*k_right_*delta_x*delta_x;
	}
	else {
		return 0.0;
	}
}


// Dummy potential
Bias::ZeroPotential::ZeroPotential(
		const StringIt& first_parameter, const StringIt& end, const double kBT):
	Potential(kBT)
{
	int num_tokens = std::distance(first_parameter, end);
	if ( num_tokens != 0 ) {
		throw std::runtime_error("invalid dummy ('none') potential");
	}
}
double Bias::ZeroPotential::evaluate(const double x) const 
{
	return 0.0;
}
