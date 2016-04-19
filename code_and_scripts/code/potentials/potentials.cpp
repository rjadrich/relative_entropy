///////////////////////////////////////////////////////
#include "../dtypes.h" //need this in every potential file
#include "potentials.h"//need this in every potential file
///////////////////////////////////////////////////////

//DEFINES THE NUMBER OF PARAMETERS NEEDED FOR EACH POTENTIAL AND RETURN TRUE IF SET CORRECTLY
bool potential_data::set_potential(int potential_type_input)
{
	potential_type = potential_type_input;

	///////////////   ||   /////////////////////////   ||   ///////////////
	///////////////   ||   /////////////////////////   ||   ///////////////
	///////////////   ||   ///ADD POTENTIALS HERE///   ||   ///////////////
	///////////////  \  /  /////////////////////////  \  /  ///////////////
	///////////////   \/   /////////////////////////   \/   ///////////////

	if (potential_type == 0) //ideal_cluster_potential
	{
		num_parameters = 7;
		num_d_parameters = 4;
		return true;
	}
	if (potential_type == 1) //ramp_salr_cluster_potential
	{
		num_parameters = 6;
		num_d_parameters = 4;
		return true;
	}
	else
	{
		return false;
	}

	//NEW POTENTIALS...

}

//FUNCTIONS TO FETCH PRIVATE PARAMETERS
int potential_data::get_num_parameters()
{
	return num_parameters;
}
int potential_data::get_num_d_parameters()
{
	return num_d_parameters;
}

//WRAPPER FOR THE VARIOUS POTENTIALS
array_pair_and_num_elements potential_data::optimize_potential(int last_step, array_pair_and_num_elements gr_data,
	double *potential_parameters, double *d_potential_parameters, gromacs_settings_class gromacs_settings,
	double *md_cutoff_pointer, double *unscaled_gradient_pointer, double *gr_convergence_pointer)
{

	///////////////   ||   /////////////////////////   ||   ///////////////
	///////////////   ||   /////////////////////////   ||   ///////////////
	///////////////   ||   ///ADD POTENTIALS HERE///   ||   ///////////////
	///////////////  \  /  /////////////////////////  \  /  ///////////////
	///////////////   \/   /////////////////////////   \/   ///////////////

	if (potential_type == 0) //ideal_cluster_potential
	{
		return optimize_ideal_cluster_potential(last_step, gr_data,
			potential_parameters, d_potential_parameters, gromacs_settings,
			md_cutoff_pointer, unscaled_gradient_pointer, gr_convergence_pointer);
	}
	if (potential_type == 1) //ramp_salr_cluster_potential
	{
		return optimize_ramp_salr_cluster_potential(last_step, gr_data,
			potential_parameters, d_potential_parameters, gromacs_settings,
			md_cutoff_pointer, unscaled_gradient_pointer, gr_convergence_pointer);
	}

	//NEW POTENTIALS...

}