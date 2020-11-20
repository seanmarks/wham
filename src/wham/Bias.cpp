#include "Bias.h"

#include "GenericFactory.h"


namespace PotentialRegistry {

// Register Potentials using a GenericFactory
// - Typedef for brevity
template<typename P>
using Register = RegisterInFactory<
	Potential, P, std::string,            // base class, derived class, and generating key
  const ParameterPack&, const double&   // input types
>;
static const Register<HarmonicPotential>      register_HarmonicPotential("harmonic");
static const Register<LinearPotential>        register_LinearPotential("linear");
static const Register<LeftHarmonicPotential>  register_LeftHarmonicPotential("left_harmonic");
static const Register<RightHarmonicPotential> register_RightHarmonicPotential("right_harmonic");
static const Register<ZeroPotential>          register_ZeroPotential("none");

// Convenient reference to the static GenericFactory that handles the creation of Potentials
using PotentialFactory = GenericFactory< Potential, std::string, const ParameterPack&, const double& >;
static PotentialFactory& factory = PotentialFactory::factory();

} // end namespace PotentialRegistry



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
		auto tmp_ptr = std::unique_ptr<Potential>(
			PotentialRegistry::factory.create(type, potential_pack, kBT_) 
		);

		// Record the term, and the name of the order parameter it affects
		potential_ptrs_.push_back( std::move(tmp_ptr) );
		order_parameter_names_.push_back( op_name );
	}
}



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



HarmonicPotential::HarmonicPotential(
		const ParameterPack& input_pack, const double kBT):
	Potential(kBT)
{
	using KeyType = ParameterPack::KeyType;
	input_pack.readNumber("x_star", KeyType::Required, x_star_);
	input_pack.readNumber("kappa",  KeyType::Required, kappa_);
	kappa_ *= beta_;
}



LinearPotential::LinearPotential(
		const ParameterPack& input_pack, const double kBT):
	Potential(kBT)
{
	using KeyType = ParameterPack::KeyType;
	input_pack.readNumber("phi", KeyType::Required, phi_);
	input_pack.readNumber("c",   KeyType::Required, c_);
	phi_ *= beta_;
	c_   *= beta_;
}



LeftHarmonicPotential::LeftHarmonicPotential(
		const ParameterPack& input_pack, const double kBT):
	Potential(kBT)
{
	using KeyType = ParameterPack::KeyType;
	input_pack.readNumber("x_left", KeyType::Required, x_left_);
	input_pack.readNumber("k_left", KeyType::Required, k_left_);
	k_left_ *= beta_;  // convert to kBT
}



RightHarmonicPotential::RightHarmonicPotential(const ParameterPack& input_pack, const double kBT):
	Potential(kBT)
{
	using KeyType = ParameterPack::KeyType;
	input_pack.readNumber("x_right", KeyType::Required, x_right_);
	input_pack.readNumber("k_right", KeyType::Required, k_right_);
	k_right_ *= beta_;  // convert to kBT
}



ZeroPotential::ZeroPotential(const ParameterPack& input_pack, const double kBT):
	Potential(kBT)
{}