
#ifndef BIAS_H
#define BIAS_H

#include <cstdlib>
#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "StringTools.h"


class Bias 
{
 public:
	Bias(const std::string& input_line, const double kBT = 1.0);

	// Returns the total bias at x, in kBT
	double evaluate(const std::vector<double>& x) const;

	const std::vector<std::string>& get_order_parameter_names() const {
		return order_parameter_names_;
	}

 private:
	// Convenient typedefs
	using StringIt = std::vector<std::string>::iterator;


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

	// Linear potential
	class LinearPotential : public Potential
	{
	 public:
		LinearPotential(const double phi, const double c) :
			Potential(), phi_(phi), c_(c) {};

		virtual double evaluate(const double x) const override
		{
			double delta_x = x - phi_;
			return 0.5*c_*delta_x*delta_x;
		}

	 private:
		double phi_, c_;  // both in kBT
	};

	// Left (one-sided) harmonic potential
	class LeftHarmonicPotential : public Potential
	{
	 public:
		LeftHarmonicPotential(const double x_left, const double k_left) :
			Potential(), x_left_(x_left), k_left_(k_left) {};

		virtual double evaluate(const double x) const override
		{
			if ( x < x_left_ ) {
				double delta_x = x - x_left_;
				return 0.5*k_left_*delta_x*delta_x;
			}
			else {
				return 0.0;
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

	// One order parameter per term in the potential
	std::vector<std::string> order_parameter_names_;

	// Thermodynamic beta = 1.0/(k_B*T)
	double kBT_, beta_;

	//--------------------------//
	//----- Helper Methods -----//
	//--------------------------//
	
	// Starting from 'start', search the given range for the next term in
	// the potential, which is enclosed in braces, { }. If there are no
	// terms left to parse, 'first' and 'last' are set equal to 'stop'.
	// - Output: iterators bounding the tokens inside the braces
	// - Returns: iterator to the next token to check (the one beyond the closing brace)
	StringIt find_next_potential(
		const StringIt& start,  // where to start searching
		const StringIt& stop,   // where to stop searching (one past the end)
		// Output:
		StringIt& first, // first token in the next potential
		StringIt& last   // one past the last token in the next potential
	) {
		bool found_left_brace = false;
		for ( auto it = start; it != stop; ++it ) {
			if ( *it == "{" ) {
				found_left_brace = true;
				first = ++it;
			}
			else if ( *it == "}" ) {
				if ( found_left_brace ) {
					last = it;
					return ++it;  // first token for next call is the one after the '}'
				}
				else {
					throw std::runtime_error("Error parsing biases: a potential is missing its left brace");
				}
			}
		}

		if ( not found_left_brace ) {
			// No potential terms left to parse
			first = stop;
			last = stop;
			return stop;
		}
		else {
			throw std::runtime_error("Error parsing biases: a potential is missing its right brace");
		}
	}
};

#endif /* BIAS_H */
