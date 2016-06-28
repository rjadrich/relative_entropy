#ifndef POTENTIALS_H
#define POTENTIALS_H

//CLASS DEFINITION FOR THE VARIOUS POTENTIALS
class potential_data
{

public:
	//INITIALIZES THE INPUT ARRAYS NEEDED FOR EACH POTENTIAL
	bool set_potential(int potential_type_input);

	//FETCHES PRIVATE VARIABLES
	int get_num_parameters();
	int get_num_d_parameters();

	//THIS IS JUST A WRAPPER FUNCTION THAT USES THE EXPLICIT POTENTIAL DEFINITIONS BELOW
	array_pair_and_num_elements optimize_potential(int last_step, array_pair_and_num_elements gr_data,
		double *potential_parameters, double *d_potential_parameters, gromacs_settings_class gromacs_settings,
		double *md_cutoff_pointer, double *unscaled_gradient_pointer, double *gr_convergence_pointer);

private:
	int potential_type;
	int num_parameters;
	int num_d_parameters;
	int num_parameters_rows;
	int num_d_parameters_rows;

	///////////////   ||   /////////////////////////   ||   ///////////////
	///////////////   ||   /////////////////////////   ||   ///////////////
	///////////////   ||   ///ADD POTENTIALS HERE///   ||   ///////////////
	///////////////  \  /  /////////////////////////  \  /  ///////////////
	///////////////   \/   /////////////////////////   \/   ///////////////

	//from ideal_cluster_potential.cpp
	array_pair_and_num_elements optimize_ideal_cluster_potential(int last_step, array_pair_and_num_elements gr_data,
		double *potential_parameters, double *d_potential_parameters, gromacs_settings_class gromacs_settings,
		double *md_cutoff_pointer, double *unscaled_gradient_pointer, double *gr_convergence_pointer);

	//from ramp_salr_cluster_potential.cpp
	array_pair_and_num_elements optimize_ramp_salr_cluster_potential(int last_step, array_pair_and_num_elements gr_data,
		double *potential_parameters, double *d_potential_parameters, gromacs_settings_class gromacs_settings,
		double *md_cutoff_pointer, double *unscaled_gradient_pointer, double *gr_convergence_pointer);

	//from crystal_potential.cpp
	array_pair_and_num_elements optimize_crystal_potential(int last_step, array_pair_and_num_elements gr_data,
		double *potential_parameters, double *d_potential_parameters, gromacs_settings_class gromacs_settings,
		double *md_cutoff_pointer, double *unscaled_gradient_pointer, double *gr_convergence_pointer);

	//from splined_potential.cpp
	array_pair_and_num_elements optimize_splined_potential(int last_step, array_pair_and_num_elements gr_data,
		double *potential_parameters, double *d_potential_parameters, gromacs_settings_class gromacs_settings,
		double *md_cutoff_pointer, double *unscaled_gradient_pointer, double *gr_convergence_pointer);

	//from splined_potential_standard.cpp
	array_pair_and_num_elements optimize_splined_potential_standard(int last_step, array_pair_and_num_elements gr_data,
		double *potential_parameters, double *d_potential_parameters, gromacs_settings_class gromacs_settings,
		double *md_cutoff_pointer, double *unscaled_gradient_pointer, double *gr_convergence_pointer);

	//NEW POTENTIALS...

};

#endif