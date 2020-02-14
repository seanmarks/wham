
#ifndef BIAS_H
#define BIAS_H

#include <cstdlib>
#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "InputParser.h"
#include "StringTools.h"

class Potential;

class Bias 
{
 public:
	Bias(const ParameterPack& input_pack, const double kBT = 1.0);

	// Returns the total bias at x, in kBT
	double evaluate(const std::vector<double>& x) const;

	const std::vector<std::string>& get_order_parameter_names() const {
		return order_parameter_names_;
	}

	const std::string& get_data_set_label() const { 
		return data_set_label_; 
	}

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

 protected:
	double kBT_, beta_;
};

// Harmonic potential
class HarmonicPotential : public Potential
{
 public:
	HarmonicPotential(const ParameterPack& input_pack, const double kBT);
	virtual double evaluate(const double x) const override;

 private:
	double x_star_, kappa_;  // kappa in kBT
};

// Linear potential
class LinearPotential : public Potential
{
 public:
	LinearPotential(const ParameterPack& input_pack, const double kBT);
	virtual double evaluate(const double x) const override;

 private:
	double phi_, c_;  // both in kBT
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
	virtual double evaluate(const double x) const override;
};

#endif /* BIAS_H */
