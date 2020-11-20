// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#pragma once
#ifndef WHAM_H
#define WHAM_H

// Standard headers
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

// dlib
#include "dlib/matrix.h"
#include "dlib/optimization.h"

#include "Bias.h"
#include "Bins.h"
#include "FreeEnergyDistribution.hpp"
#include "FreeEnergyDistribution2D.hpp"
#include "FileSystem.h"
#include "GptlWrappers.h"
#include "InputParser.h"
#include "Matrix.hpp"
#include "OpenMP.h"
#include "OrderParameter.h"
#include "OrderParameterRegistry.h"
#include "Simulation.h"
#include "WhamDlibWrappers.h"


// Solves binless WHAM equations via log-likelihood maximization
// - Minimization uses BFGS implementation from dlib library
// - Computes consensus estimates using optimal biasing free energies
//   - TODO: move to new objects
// - References
//   - Chodera, ..., Dill (JCTC 2007)
//   - Kong, ..., Tan (J. R. Statist. Soc. B 2003)
//   - Shirts & Chodera (JCP 2008)
//   - Souaille & Roux (Comp. Phys. Comm. 2001)
//   - Tan, Gallicchio, Lapelosa, & Levy (J. Chem. Phys. 2012)
//   - Zhu & Hummer (J. Comp. Chem. 2011)
class Wham
{
 public:
  // dlib types
  using ColumnVector = dlib::matrix<double,0,1>;

  template<typename T>
  using Matrix = numeric::Matrix<T>;

	template<typename T>
	using Vector = std::vector<T>;


  Wham() = delete;

  // Performs all necessary setup, then solves the WHAM equations
  // for the given set of data. Output is then available through
  // accessor methods.
  Wham(
    //const DataSummary& data_summary,
    const OrderParameterRegistry& op_registry,
    const std::vector<Simulation>& simulations,
    const std::vector<OrderParameter>& order_parameters,
    const std::vector<Bias>& biases,
    const std::vector<double>& f_bias_guess,  // initial guess
    const double tol  // solver tolerance
  );


  //----- Objective Function -----//

  // TODO: make private/protected?
  double evalObjectiveFunction(const ColumnVector& df) const;

  // Compute gradient of objective function
  // - Note:
  //   - grad  = gradient wrt. free energy differences between neighboring windows
  //   - dA_df = gradient wrt. biasing free energies themselves
  const ColumnVector evalObjectiveDerivatives(const ColumnVector& df) const;


  //----- WHAM Estimators -----//

  /// Free energies of turning on the bias, \f$ \Delta F_{k} \f$, at the solution
  const std::vector<double>& get_f_bias_opt() const noexcept {
    return f_bias_opt_;
  }

  /// Returns the total number of samples, N, across all simulations
  int getNumSamplesTotal() const noexcept {
    return num_samples_total_;
  }

  const std::vector<Simulation>& getSimulationData() const noexcept {
    return simulations_;
  }

  const std::vector<std::pair<int,int>>& getSimulationDataRanges() const noexcept {
    return simulation_data_ranges_;
  }

  // Returns the log of the simulation-independent factors, \hat{D}(x_n)
  //   \hat{D}(x_n) = \sum{k=1}^{K} N_k * exp[f_k - u_k(x_n)]
  const std::vector<double>& get_log_dhat_opt() const {
    return log_dhat_opt_;
  }

	/// Calculates the consensus weights given to each sample in the ensemble
	/// corresponding to the given values of the bias and associated free energy
	/// @param u_bias value of the bias corresponding to each sample: size = N [kBT]
  /// @param f_bias free energy of turning on the bias, Delta_F_k = F_k - F_0 [kBT]
	/// @return weights of each sample in the given ensesmble (note: they sum to 1)
  std::vector<double> calculateWeights(
    const std::vector<double>& u_bias,
    const double f_bias = 0.0
  ) const;

  /// Calculates consensus weights in the unbiased ensemble
  std::vector<double> calculateWeightsForUnbiasedEnsemble() const;

	/// Calculates consensus weights in the ensemble corresponding to the simulation
	/// with the given data set label
  std::vector<double> calculateWeightsForSimulation(
    const std::string& data_set_label
  ) const;
  

  // "Manually" unbias the distributions for the given OrderParameter (i.e. using only
  // each individual simulation's data, not the consensus estimates)
	// - TODO: move to separate class
  //std::vector<FreeEnergyDistribution> manuallyUnbiasDistributions(const std::string& op_name) const;

	// Computes the logarithm of the weights assigned to each data point in the
	// unbiased ensemble, using only data from a single simulation
	// - Size: (K, N_j) = (num_simulations x num_samples)
	Vector<Vector<double>> computeUnbiasedNonConsensusLogWeights() const;

 private:
  // Objects owned by the driver
  const OrderParameterRegistry&      op_registry_;
  const std::vector<Simulation>&     simulations_;
  const std::vector<OrderParameter>& order_parameters_;
  const std::vector<Bias>&           biases_;

  // The value of the bias for each sample, evaluating using each potential (i.e. in each ensemble):
  //    u_{bias,r}( x_{j,i} )
  //      r = 0, ..., m-1    (m = number of biasing potentials/ensembles)
  //      j = 0, ..., m-1    (m = num_simulations)
  //      i = 0, ..., n_j-1  (n_j = # samples from simulation j)
  // - For each bias r = 0, ...,
  //     u_bias_as_other_[r] = [(data_sim_0), ..., (data_sim_j), ... , (data_sim_m)]
  //                                                     |
  //                                                N_j samples
  //   - Overall length of u_bias_as_other_[r]:  N = sum_{j=0}^{m} N_j
  std::vector<std::vector<double>> u_bias_as_other_;

  // For each simulation j = 0, ..., m-1:
  //  simulation_data_ranges_[j] = indices (first, end) for the 'j'th simulation's data
  //                               in u_bias_as_other_ (for each 'r')
  std::vector<std::pair<int,int>> simulation_data_ranges_;

  // Solver tolerance
  const double tol_ = 1.0e-7;

  // Sample counts
  std::vector<int> num_samples_per_simulation_;   // N_j = number of samples from each simulation, j
  int    num_samples_total_;      // N = total number of samples (across all simulations)
  double inv_num_samples_total_;  // precompute for speed

  // c[r] = N_r/N_tot
  //      = fraction of samples from simulation 'r'
  //      = (num. samples from simulation r)/(num. total)
  std::vector<double> c_;

  // Free energy of turning on each bias
  std::vector<double> f_bias_guess_;  // initial guess
  std::vector<double> f_bias_opt_;    // optimal

  // Ensemble-independent factors that are related to the weight given to each data sample
  // - Formula:
  //     \\hat{D}(x_n) = sum_{k=K}^m exp{ log(N_k) + f_k - u_{bias,k}(x_{j,i}) }
  // - Computed as a log-sum-exp
  std::vector<double> log_dhat_opt_;


  //----- Working Variables -----//

  // u_bias_as_other for the unbiased ensemble (all zeros)
  std::vector<double> u_bias_as_other_unbiased_;
  const double f_unbiased_ = 0.0;  // Free energy of going to the *unbiased* ensemble is zero

  // Cached results
  // - These are computed by Wham.evalObjectiveFunction and saved for 
  //   re-use  in Wham.evalObjectiveDerivatives
  //   - The dlib optimization algorithm used guarantees that these functions are
  //     called one after another for each iteration
  mutable std::vector<double> log_dhat_tmp_;
  mutable std::vector<double> f_bias_last_;  // last f-values used in iteration

  // Misc. buffers
  mutable std::vector<std::vector<double>> minus_log_sigma_k_binned_;  // for binning samples
  mutable std::vector<std::vector<double>> args_buffers_;  // for log_sum_exp during optimization


  //----- Solve WHAM Equations -----//

  void setup();

  // For each sample (across all simulations), evaluate the bias that would be
  // felt under each simulation's potential
  // - TODO: move to a completely separate class that manages biases, and simply pass in
  //   the resulting values to this object
  void evaluateBiases();

  // Solve for f_bias_opt
  std::vector<double> solveWhamEquations(const std::vector<double>& f_bias_guess);


  //----- Helper Functions -----//

  // Compute log( \hat{D}(x_n) ) for the given set of biasing free energies
  // - These are the common denominator of the consensus weights given to each sample
  // *** NOT thread safe: encloses an OpenMP region
  void compute_log_dhat(
    const std::vector<std::vector<double>>& u_bias_as_other,
    const std::vector<double>&              fhat,
    // Output
    std::vector<double>& log_dhat
  ) const;

  // Returns the logarithm of a sum of exponentials:
  //      log(S) = sum_{i=1}^N exp(args[i])
  // - input: arguments of exponentials
  // - THREAD_SAFE
  double logSumExp(const std::vector<double>& args) const;

  // Returns the weighted sum of a series of exponentials
  //      log(S) = sum_{i=1}^N weights[i]*exp(args[i])
  // - THREAD_SAFE
  double weightedSumExp(
    const std::vector<double>& args,    // arguments of exponentials
    const std::vector<double>& weights  // not necessarily positive
  ) const;

  // Convert between the free energies of turning on the bias (f) and the free energy 
  // *differences* between windows (df), assuming f[0] = f0 = 0.0
  void convert_f_to_df(const std::vector<double>& f, ColumnVector& df) const;
  void convert_df_to_f(const ColumnVector& df, std::vector<double>& f) const;


  //----- Consensus Estimates -----//

  // Free energy of turning on arbitrary bias 'k', given the biasing potential
  // of each sample evaluated in that ensemble
  double computeBiasingFreeEnergy(
    const std::vector<double>& u_bias_as_k
  ) const;

  // Smallest value to for which exp(x) will not underflow (about -36)
  static constexpr double MIN_DBL_FOR_EXP = std::log( std::numeric_limits<double>::epsilon() );


  //----- GPTL -----//

  using Timer = GPTL::Timer;

  // Core functions
  mutable Timer setup_timer_     = Timer("Wham::setup");
  mutable Timer biases_timer_    = Timer("Wham::evaluate_biases");
  mutable Timer solve_timer_     = Timer("Wham::solve");
  mutable Timer objective_timer_ = Timer("Wham::objective");
  mutable Timer gradient_timer_  = Timer("Wham::gradient");
  mutable Timer gradient_omp_timer_ = Timer("Wham::gradient_omp");

  // Low-level, expensive functions
  mutable Timer log_sum_exp_timer_      = Timer("Wham::compute_log_sum_exp");
  mutable Timer weighted_sum_exp_timer_ = Timer("Wham::weighted_sum_exp");
  mutable Timer log_dhat_timer_     = Timer("Wham::compute_log_dhat");
  mutable Timer log_dhat_omp_timer_ = Timer("Wham::compute_log_dhat_omp");
};

#endif // ifndef WHAM_H
