
#ifndef SIMULATION_H
#define SIMULATION_H

// Struct with simulation-specific information
struct Simulation 
{
	std::string data_set_label;

	double kBT, beta;    // k_B*T and beta = 1/(k_B*T)
	double t_min, t_max; // Sampling range [ps]

	double f_bias;       // Free energy of adding bias (in kBT): f_bias = -ln(Q_i/Q_0)
	double f_bias_guess; // First estimate of biasing free energy
};

#endif /* ifndef SIMULATION_H */
