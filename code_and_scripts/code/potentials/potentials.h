#ifndef POTENTIALS_H
#define POTENTIALS_H

//SPECIFIC POTENTIAL FUNCTION DEFINITIONS
//THE POTENTIALS ARE STORED IN SEPARATE CPP FILES FOR CONVENIENCE

//from ideal_cluster_potential.cpp
array_pair_and_num_elements optimize_ideal_cluster_potential(int last_step, array_pair_and_num_elements gr_data,
	double *potential_parameters, double *d_potential_parameters, gromacs_settings_class gromacs_settings,
	double *md_cutoff_pointer, double *gr_cutoff_pointer, double *unscaled_gradient_pointer, double *gr_convergence_pointer);

//from ramp_salr_cluster_potential.cpp
array_pair_and_num_elements optimize_ramp_salr_cluster_potential(int last_step, array_pair_and_num_elements gr_data,
	double *potential_parameters, double *d_potential_parameters, gromacs_settings_class gromacs_settings,
	double *md_cutoff_pointer, double *gr_cutoff_pointer, double *unscaled_gradient_pointer, double *gr_convergence_pointer);

#endif