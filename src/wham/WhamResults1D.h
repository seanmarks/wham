// WhamResults1D
// - Simple class for organizing results for a 1-D distribution
// - Includes non-consensus (single-simulation) estimates as part of a full 
//   battery of diagnostics, in addition to the conventional F_WHAM(x)
// - All energies are in kBT

#include "Distribution.h"
#include "OrderParameter.h"
#include "Simulation.h"

#pragma once
#ifndef WHAM_RESULTS_1D_H
#define WHAM_RESULTS_1D_H

class WhamResults1D
{
 public:
	//----- Data -----//

	const OrderParameter& x;
	const std::vector<Simulation>& simulations;

	Distribution f_x_wham;  // consensus distribution

	// Free energy of turning on the bias (in kBT): f_bias = Delta F_bias = -ln(Q_i/Q_0)
	//std::vector<double> f_bias;

	// Unbiased, non-consensus distributions ("manually" unbiased)
	std::vector<Distribution> f_x_unbiased;  // manually unbiased
	std::vector<Distribution> f_x_shifted;   // shifted by f_bias

	// Rebiased distributions and associated entropy
	std::vector<Distribution> f_x_rebiased;
	std::vector<double>       info_entropy;


	//----- Setup -----//

	WhamResults1D(const OrderParameter& x, const std::vector<Simulation>& simulations):
		x(x), simulations(simulations)
	{
		int num_simulations = simulations.size();
		f_x_unbiased.reserve(num_simulations);
		f_x_shifted.reserve(num_simulations);
		f_x_rebiased.reserve(num_simulations);
		info_entropy.reserve(num_simulations);
	}


	//----- Print Output -----//

	void print() const;

	void printConsensusResults() const;
	void printStats(std::string file_name = "") const;

	// Prints a series of distributions, F_i(x), side-by-side
	void printDistributions(
		const std::vector<Distribution>& distributions,
		const std::string& file_name, 
		const std::string& header,
		const bool shift_to_zero = true
	) const;

 private:
};

#endif // ifndef WHAM_RESULTS_1D_H
