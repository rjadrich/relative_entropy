///////////////////////////////////////////////////////
#include "../dtypes.h" //need this in every potential file
#include "potentials.h"//need this in every potential file
///////////////////////////////////////////////////////

//USE GLOBAL LOG FILE DEFINED IN MAIN CODE
#include <fstream>
#include <stdlib.h>
using namespace std;
extern ofstream log_filestream;

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
	else if (potential_type == 1) //ramp_salr_cluster_potential
	{
		num_parameters = 6;
		num_d_parameters = 4;
		return true;
	}
	else if (potential_type == 2) //crystal_potential
	{
		num_parameters = 9;
		num_d_parameters = 9;
		return true;
	}
	else if (potential_type == -1) //akima splined potential (for this I read in choice from a file)
	{
		ifstream spline_filestream;
		spline_filestream.open("num_spline_parameters.txt");
		if (spline_filestream >> num_parameters)
		{
			num_d_parameters = num_parameters;
			if (num_parameters < 6) //chose 6 sort of arbitrarily
			{
				log_filestream << "not enough spline parameters (minimum is set to 6): " << num_parameters << " -> killing!" << endl;
				log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
				exit(EXIT_FAILURE);
			}
			else
			{
				log_filestream << "set the number of spline parameters: " << num_parameters << endl;
			}
		}
		else
		{
			log_filestream << "could not read in the number of spline parameters (file missing?) -> killing!" << endl;
			log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
			exit(EXIT_FAILURE);
		}
		spline_filestream.close();
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
	else if (potential_type == 1) //ramp_salr_cluster_potential
	{
		return optimize_ramp_salr_cluster_potential(last_step, gr_data,
			potential_parameters, d_potential_parameters, gromacs_settings,
			md_cutoff_pointer, unscaled_gradient_pointer, gr_convergence_pointer);
	}
	else if (potential_type == 2) //ramp_salr_cluster_potential
	{
		return optimize_crystal_potential(last_step, gr_data,
			potential_parameters, d_potential_parameters, gromacs_settings,
			md_cutoff_pointer, unscaled_gradient_pointer, gr_convergence_pointer);
	}
	else if (potential_type == -1) //akima splined potential
	{
		log_filestream << "akima splined potential not coded yet -> killing!" << endl;
		log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
		exit(EXIT_FAILURE);
	}

	//NEW POTENTIALS...

}