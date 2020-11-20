// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#ifndef BIAS_H
#define BIAS_H

#include <cstdlib>
#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "Assert.hpp"
#include "InputParser.h"

class Potential;

class Bias
{
 public:
	Bias(const ParameterPack& input_pack, const double kBT = 1.0);

	// Returns the total bias at x, in kBT
	double evaluate(const std::vector<double>& x) const;

	// Returns the list of names of order parameters required to
	// evaluate the bias
	const std::vector<std::string>& getOrderParameterNames() const {
		return order_parameter_names_;
	}

	const std::string& getDataSetLabel() const { 
		return data_set_label_; 
	}

	template<typename DataSetIt, typename OutputIt>
	void evaluate(DataSetIt first, DataSetIt last, OutputIt d_first) const;


 private:
	// Terms in the total biasing potential
	std::vector<std::unique_ptr<Potential>> potential_ptrs_;

	// One order parameter per term in the potential
	std::vector<std::string> order_parameter_names_;

	std::string data_set_label_;

	// Thermodynamic beta = 1.0/(k_B*T)
	double kBT_, beta_;
};



// Base class for implementing terms in the bias
class Potential
{
 public:
	Potential(const double kBT): kBT_(kBT), beta_(1.0/kBT_) {};

	virtual ~Potential() {};

	// Returns the value of the bias (in kBT)
	virtual double evaluate(const double x) const = 0;

	// Returns the value of the bias (in kBT) for a set of input values
	virtual std::vector<double> evaluate(const std::vector<double>& x) {
		const int num = x.size();
		std::vector<double> u(num);
		for ( int i=0; i<num; ++i ) {
			u[i] = evaluate(x[i]);
		}
		return u;
	}

 protected:
	double kBT_, beta_;
};



template<typename DataSetIt, typename OutputIt>
void Bias::evaluate(DataSetIt first, DataSetIt last, OutputIt d_first) const {

	const int num_potentials = potential_ptrs_.size();
	const int num_data_sets  = std::distance(first, last);
	FANCY_ASSERT( num_potentials == num_data_sets, "length mismatch" );

	if ( num_potentials == 0 ) {
		return;
	}

	//const int num_samples = (*first)->size();
	//std::vector<double> u_tmp(num_samples);
	auto data_it = first;
	for ( int p=0; p<num_potentials; ++p, ++data_it ) {
		const auto& data = *(*data_it);
		auto u_tmp = potential_ptrs_[p]->evaluate(data);

		auto out_it = d_first;
		for ( auto it = u_tmp.begin(); it != u_tmp.end(); ++it, ++out_it ) {
			*out_it += *it;
		}
	}
}



// Linear potential
class LinearPotential : public Potential
{
 public:
	LinearPotential(const ParameterPack& input_pack, const double kBT);

	virtual double evaluate(const double x) const override;

 private:
	double phi_, c_;  // both in kBT
};



// Harmonic potential
class HarmonicPotential : public Potential
{
 public:
	HarmonicPotential(const ParameterPack& input_pack, const double kBT);
	
	virtual double evaluate(const double x) const override;
	
	virtual std::vector<double> evaluate(const std::vector<double>& x) override {
		const int num = x.size();
		std::vector<double> u(num);
		const double fac = 0.5*kappa_;
		for ( int i=0; i<num; ++i ) {
			double delta_x = x[i] - x_star_;
			u[i] = fac*delta_x*delta_x;
		}
		return u;
	}

 private:
	double x_star_, kappa_;  // kappa in kBT
};



// Left (one-sided) harmonic potential
class LeftHarmonicPotential : public Potential
{
 public:
	LeftHarmonicPotential(const ParameterPack& input_pack, const double kBT);
	
	virtual double evaluate(const double x) const override;

 private:
	double x_left_, k_left_;  // k_left in kBT
};



// Right (one-sided) harmonic potential
class RightHarmonicPotential : public Potential
{
 public:
	RightHarmonicPotential(const ParameterPack& input_pack, const double kBT);
	virtual double evaluate(const double x) const override;

 private:
	double x_right_, k_right_;  // k_right in kBT
};


// Dummy potential (for unbiased variables)
class ZeroPotential : public Potential
{
 public:
	ZeroPotential(const ParameterPack& input_pack, const double kBT);

	virtual double evaluate(const double x) const override {
		return 0.0;
	}
};



//----- Inline Functions -----//


inline
double LinearPotential::evaluate(const double x) const 
{
	return phi_*x + c_;
}


inline
double HarmonicPotential::evaluate(const double x) const 
{
	double delta_x = x - x_star_;
	return 0.5*kappa_*delta_x*delta_x;
}


inline
double LeftHarmonicPotential::evaluate(const double x) const
{
	const double delta_x = x - x_left_;
	if ( delta_x < 0.0 ) {
		return 0.5*k_left_*delta_x*delta_x;
	}
	else {
		return 0.0;
	}
}


inline
double RightHarmonicPotential::evaluate(const double x) const
{
	const double delta_x = x - x_right_;
	if ( delta_x > 0.0 ) {
		return 0.5*k_right_*delta_x*delta_x;
	}
	else {
		return 0.0;
	}
}

#endif // ifndef BIAS_H