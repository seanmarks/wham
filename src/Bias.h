
#ifndef BIAS_H
#define BIAS_H

#include <cstdlib>
#include <exception>
#include <memory>
#include <stdexcept>
#include <vector>

class Bias 
{
 public:
	Bias(const std::vector<std::string>& input_tokens, const double kBT = 1.0);

	// Returns the total bias at x, in kBT
	double evaluate(const std::vector<double>& x) const;

 private:
	//-------------------------------//
	//----- Potential Functions -----//
	//-------------------------------//

	// Base class for implementing terms in the bias
	class Potential
	{
	 public:
		Potential() {};

		// Returns the value of the bias (in kBT)
		virtual double evaluate(const double x) const = 0;
	};

	// Harmonic potential
	class HarmonicPotential : public Potential
	{
	 public:
		HarmonicPotential(const double x_star, const double kappa) :
			Potential(), x_star_(x_star), kappa_(kappa) {};

		virtual double evaluate(const double x) const override
		{
			double delta_x = x - x_star_;
			return 0.5*kappa_*delta_x*delta_x;
		}

	 private:
		double x_star_, kappa_;  // kappa in kBT
	};

	// Left (one-sided) harmonic potential
	class LeftHarmonicPotential : public Potential
	{
	 public:
		LeftHarmonicPotential(const double x_left, const double k_left) :
			Potential(), x_left_(x_left), k_left_(k_left) {};

		virtual double evaluate(const double x) const override
		{
			if ( x > x_left_ ) {
				return 0.0;
			}
			else {
				double delta_x = x - x_left_;
				return 0.5*k_left_*delta_x*delta_x;
			}
		}

	 private:
		double x_left_, k_left_;  // k_left in kBT
	};

	// Dummy potential (for unbiased variables)
	class ZeroPotential : public Potential
	{
	 public:
		ZeroPotential() {};

		virtual double evaluate(const double x) const override 
		{
			return 0.0;
		}
	};

	//-----------------------------//
	//----- Private Variables -----//
	//-----------------------------//

	// Terms in the total biasing potential
	std::vector<std::unique_ptr<Potential>> potential_ptrs_;

	// Thermodynamic beta = 1.0/(k_B*T)
	double kBT_, beta_;
};

#endif /* BIAS_H */
