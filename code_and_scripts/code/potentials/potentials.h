#ifndef POTENTIALS_H
#define POTENTIALS_H

///////////////////////////////////////////////////////
//////CLASS DEFINITION FOR THE VARIOUS POTENTIALS//////
///////////////////////////////////////////////////////
/*class potential
{
private:
	int num_parameters;
	int num_d_parameters;
	double *potential_parameters; 
	double *d_potential_parameters; 

public:
	void set_potential(int potential_type)
	{
		if (potential_type == 0) //ideal_cluster_potential
		{
			num_parameters = 6;
			num_d_parameters = 4;
		}
		if (potential_type == 1) //ramp_salr_cluster_potential
		{
			num_parameters = 6;
			num_d_parameters = 4;
		}
		
		//ALLOCATE MEMORY FOR READING IN THE PARAMETERS
		potential_parameters = (double*)malloc(sizeof(double) * num_parameters);
		d_potential_parameters = (double*)malloc(sizeof(double) * num_d_parameters);
	}
};*/

///////////////////////////////////////////////////////
///////SPECIFIC POTENTIAL FUNCTION DEFINITIONS/////////
///////////////////////////////////////////////////////

//THE POTENTIALS ARE STORED IN INDIVIDUAL CPP FILES FOR CONVENIENCE

//from ideal_cluster_potential.cpp
array_pair_and_num_elements optimize_ideal_cluster_potential(int last_step, array_pair_and_num_elements gr_data,
	double *potential_parameters, double *d_potential_parameters, gromacs_settings_class gromacs_settings,
	double *md_cutoff_pointer, double *gr_cutoff_pointer, double *unscaled_gradient_pointer, double *gr_convergence_pointer);

//from ramp_salr_cluster_potential.cpp
array_pair_and_num_elements optimize_ramp_salr_cluster_potential(int last_step, array_pair_and_num_elements gr_data,
	double *potential_parameters, double *d_potential_parameters, gromacs_settings_class gromacs_settings,
	double *md_cutoff_pointer, double *gr_cutoff_pointer, double *unscaled_gradient_pointer, double *gr_convergence_pointer);

#endif