#include "Bias.h"

Bias::Bias(const ParameterPack& input_pack, const double kBT):
	kBT_(kBT), beta_(1.0/kBT_)
{	
	using KeyType = ParameterPack::KeyType;

	input_pack.readString("DataSet", KeyType::Required, data_set_label_);

	// Find parameter packs for the terms in this bias
	std::vector<const ParameterPack*> potential_pack_ptrs 
			= input_pack.findParameterPacks("Potential", KeyType::Required);
	int num_potentials = potential_pack_ptrs.size();
	if ( num_potentials == 0 ) {
		std::cerr << "Warning: the bias for data set " << data_set_label_ 
		          << " has no terms in its potential\n";
	}

	for ( int i=0; i<num_potentials; ++i ) {
		const ParameterPack& potential_pack = *(potential_pack_ptrs[i]);

		std::string op_name, type;
		potential_pack.readString("order_parameter", KeyType::Required, op_name);
		potential_pack.readString("type",            KeyType::Required, type);

		// Construct this term
		// TODO use a factory
		// TODO other potentials
		std::unique_ptr<Potential> tmp_ptr;
		if ( type == "harmonic" ) {
			tmp_ptr = std::unique_ptr<Potential>( new HarmonicPotential(potential_pack, kBT_) );
		}
		else if ( type == "linear" ) {
			tmp_ptr = std::unique_ptr<Potential>( new LinearPotential(potential_pack, kBT_) );
		}
		else if ( type == "left_harmonic" ) {
			tmp_ptr = std::unique_ptr<Potential>( new LeftHarmonicPotential(potential_pack, kBT_) );
		}
		else if ( type == "right_harmonic" ) {
			tmp_ptr = std::unique_ptr<Potential>( new RightHarmonicPotential(potential_pack, kBT_) );
		}
		else if ( type == "none"  ) {
			tmp_ptr = std::unique_ptr<Potential>( new ZeroPotential(potential_pack, kBT_) );
		}
		else {
			throw std::runtime_error("unrecognized potential type \"" + type + "\"");
		}

		// Record the term, and the name of the order parameter it affects
		potential_ptrs_.push_back( std::move(tmp_ptr) );
		order_parameter_names_.push_back( op_name );
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
		const ParameterPack& input_pack, const double kBT):
	Potential(kBT)
{
	using KeyType = ParameterPack::KeyType;
	input_pack.readNumber("x_star", KeyType::Required, x_star_);
	input_pack.readNumber("kappa",  KeyType::Required, kappa_);
	kappa_ *= beta_;
}

double Bias::HarmonicPotential::evaluate(const double x) const 
{
	double delta_x = x - x_star_;
	return 0.5*kappa_*delta_x*delta_x;
}


// Linear potential
Bias::LinearPotential::LinearPotential(
		const ParameterPack& input_pack, const double kBT):
	Potential(kBT)
{
	using KeyType = ParameterPack::KeyType;
	input_pack.readNumber("phi", KeyType::Required, phi_);
	input_pack.readNumber("c",   KeyType::Required, c_);
	phi_ *= beta_;
	c_   *= beta_;
}
double Bias::LinearPotential::evaluate(const double x) const 
{
	return phi_*x + c_;
}


// Left one-sided harmonic potential
Bias::LeftHarmonicPotential::LeftHarmonicPotential(
		const ParameterPack& input_pack, const double kBT):
	Potential(kBT)
{
	using KeyType = ParameterPack::KeyType;
	input_pack.readNumber("x_left", KeyType::Required, x_left_);
	input_pack.readNumber("k_left", KeyType::Required, k_left_);
	k_left_ *= beta_;  // convert to kBT
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
Bias::RightHarmonicPotential::RightHarmonicPotential(const ParameterPack& input_pack, const double kBT):
	Potential(kBT)
{
	using KeyType = ParameterPack::KeyType;
	input_pack.readNumber("x_right", KeyType::Required, x_right_);
	input_pack.readNumber("k_right", KeyType::Required, k_right_);
	k_right_ *= beta_;  // convert to kBT
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
Bias::ZeroPotential::ZeroPotential(const ParameterPack& input_pack, const double kBT):
	Potential(kBT)
{}
double Bias::ZeroPotential::evaluate(const double x) const 
{
	return 0.0;
}
