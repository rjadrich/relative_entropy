///////////////////////////////////////////////////////
#include "../dtypes.h" //need this in every potential file
#include "potentials.h"//need this in every potential file
///////////////////////////////////////////////////////

//INCLUDES FOR THIS POTENTIAL
#include <math.h> //for functions
#include <stdlib.h> //for malloc
#include <iostream>
#include <fstream>
#include "akima.h"

using namespace std;

//DEFINITION OF AUXILLARY DATA EXTRACTION AND OUTPUT FUNCTIONS
void extract_data_standard(/*input*/ double *potential_parameters, int num_parameters,
	/*output*/ double *R, double *U, int *State);
void format_data_standard(/*output*/ double *potential_parameters,
	/*input*/int num_parameters, double *R, double *U, int *State);


//ACTUAL POTENTIAL OPTIMIZER FUNCTION
array_pair_and_num_elements potential_data::optimize_splined_potential_standard(int last_step, array_pair_and_num_elements gr_data,
	double *potential_parameters, double *d_potential_parameters, gromacs_settings_class gromacs_settings,
	double *md_cutoff_pointer, double *unscaled_gradient_pointer, double *gr_convergence_pointer)
{
	//OPTIMIZATION PARAMETERS
	//for this potential we just work with the potential_parameters directly

	//PARAMETERS FOR MAKING SPLINE
	double *R, *U, *grad_U;
	int *State;
	Maths::Interpolation::Akima Spline(num_parameters); //potential spline
	Maths::Interpolation::Akima Spline_L(num_parameters); //spline for derivative calc (left side)
	Maths::Interpolation::Akima Spline_R(num_parameters); //spline for derivative calc (right side)

	//ALLOCATE MEMORY FOR X, Y AND STATE (FIXED OR NOT FIXED)
	R = (double*)malloc(sizeof(double) * num_parameters); 
	U = (double*)malloc(sizeof(double) * num_parameters); //for mental simplicity the last element exists but is just ignored
	grad_U = (double*)malloc(sizeof(double) * num_parameters); //stores the various gradients
	State = (int*)malloc(sizeof(int) * num_parameters); //for mental simplicity the last element exists but is just ignored

	//SOME RELEVANT OPTIMIZATION PARAMETERS
	double delta_U, dUdu_I, dUdu_II;

	//GENERIC DISTANCE VARIABLE
	double r;
	//GENERIC LOOP INTEGER
	int i, j, k;
	//TABLE CREATION PARAMETERS
	array_pair_and_num_elements table_data; //stores the generated table data to pass back
	double *table_u, *table_du; //table arrays
	int num_table_entries; //number of elements in table arrays
	//RDF STORAGE
	int num_lines_gr; //again, better than using gr_data directly
	//double grad_rcut;
		/////////////////////////////////////////////////////////////////
	double grad_magnitude;
	double r_weight_I, r_weight_II;
	double r_I, r_II;
	double gr_I, gr_II;
	double gr_tgt_I, gr_tgt_II;
		/////////////////////////////////////////////////////////////////
	//GR CONVERGENCE
	double gr_convergence;

	
	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////UPDATING PARAMETERS///////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	//EXTRACT (FORMAT) THE DATA FOR USE IN THE AKIMA SPLINING
	extract_data_standard(/*input*/ potential_parameters, num_parameters, /*output*/ R, U, State);

	//EXTRACT MEMORY ADDRESS FROM GR_DATA
	num_lines_gr = gr_data.num_elements;

	//IF LAST STEP IS NOT ZERO CALCULATE UPDATES (THIS IS FOR INTEGRALS)
	if (last_step != 0)
	{
		//INITIALIZE THE GRADIENT ARRAY
		for (i = 0; i < num_parameters; i++)
		{
			*(grad_U + i) = 0.0;
		}

		//INTEGRAL FOR FINDING THE GRADIENT TOWARDS LARGER PROBABILITY
		//LOOP OVER THE PARAMETERS TO CALCULATE THE DERIVATIVES AND INTEGRATE
		//loop for selecting which parameter to perturb
		gr_convergence = 0.0; 
		delta_U = 0.01; //just fixed for now to test
		for (i = 0; i < num_parameters; i++)
		{
			if (*(State + i) == 1)
			{
				//make the splines for derivative calculations
				*(U + i) = *(U + i) - delta_U; //move to the left by one unit
				Spline_L.Make_Akima(R, U);
				*(U + i) = *(U + i) + 2.0 * delta_U; //to the right by one unit (requires two steps)
				Spline_R.Make_Akima(R, U);
				*(U + i) = *(U + i) - delta_U; //move back to original location

				//now compute the derivative and update
				for (k = 0; k < num_lines_gr - 1; k++)
				{
					//using trapezoids so need  two values
					r_I = (double)k*gromacs_settings.delta_r;
					r_II = (double)(k + 1)*gromacs_settings.delta_r;
					//convert to spherical weights based on dimensionality
					r_weight_I = pow(r_I, gromacs_settings.dimensions - 1);
					r_weight_II = pow(r_II, gromacs_settings.dimensions - 1);

					//fill in rdf data
					gr_I = *(gr_data.array_1 + k);
					gr_II = *(gr_data.array_1 + k + 1);
					gr_tgt_I = *(gr_data.array_2 + k);
					gr_tgt_II = *(gr_data.array_2 + k + 1);

					//calculate the derivatives
					dUdu_I = (Spline_R.getValue(r_I) - Spline_L.getValue(r_I)) / delta_U;
					dUdu_II = (Spline_R.getValue(r_II) - Spline_L.getValue(r_II)) / delta_U;

					//actually do the damn integral
					if (isfinite(dUdu_I) && isfinite(dUdu_II))
					{
						*(grad_U + i) = *(grad_U + i) + 0.5*((r_weight_I*dUdu_I*(gr_I - gr_tgt_I)) + (r_weight_II*dUdu_II*(gr_II - gr_tgt_II)))*gromacs_settings.delta_r;
					}

					//calculate the gr_convergence (only once)
					if (i == 0)
					{
						gr_convergence = gr_convergence + 0.5*((r_weight_I*(gr_I - gr_tgt_I)*(gr_I - gr_tgt_I)) + (r_weight_II*(gr_II - gr_tgt_II)*(gr_II - gr_tgt_II)))*gromacs_settings.delta_r;
					}
				}
			}
		}

		//prepare to pass back the convergence criterion
		*(gr_convergence_pointer) = gr_convergence;

		//calculate unscaled gradient and scale the gradient back by input value
		//grad_magnitude = 0.0;
		//for (i = 0; i < num_parameters; i++)
		//{
		//	grad_magnitude = grad_magnitude + (*(grad_U_step + i))*(*(grad_U_step + i)); //add square for this parameter
		//	*(grad_U_step + i) = gromacs_settings.gradient_scale*(*(grad_U_step + i)); //scale this parameters gradient
		//}
		//grad_magnitude = sqrt(grad_magnitude); //get the square root of the total gradient
		//*unscaled_gradient_pointer = grad_magnitude;



		//calculate unscaled gradient and only include what contribution will not violate a constraint
		//the component is not removed entirely as a more holistic, continuous measure is better
		grad_magnitude = 0.0;
		for (i = 0; i < num_parameters; i++)
		{
			grad_magnitude = grad_magnitude + (*(grad_U + i))*(*(grad_U + i)); //add square for this parameter
			*(grad_U + i) = gromacs_settings.gradient_scale*(*(grad_U + i)); //scale this parameters gradient
		}
		grad_magnitude = sqrt(grad_magnitude); //get the square root of the total gradient
		*unscaled_gradient_pointer = grad_magnitude;




		//MOMENTUM DISABLED FOR THIS POTENTIAL
		//TO BE UPDATED AFTER TESTING

		
		//any constraints desired...



		//UPDATE THE PARAMETERS AND FORMAT FOR DATA FILE BY STORING IN POTENTIAL_PARAMETERS
		//note: the new spline calculation happens outside this region

		for (i = 0; i < num_parameters; i++)
		{
			*(U + i) = *(U + i) + *(grad_U + i);
		}

		format_data_standard(/*output*/ potential_parameters,
			/*input*/num_parameters, R, U, State);

	}

	

	//MAKE THE SPLINE USING PARAMETERS THAT MAY OR MAY NOT HAVE BEEN UPDATED (IF STEP 0 NO)
	Spline.Make_Akima(R, U);

	


	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////TABLE AND MD CUTOFF STUFF/////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	//ALLOCATE TABLE ARRAYS (FETCHING THE FIRST ADDRESS OF A CONTIGUOUS BLOCK OF MEMORY ON THE HEAP)
	num_table_entries = 1 + (int)(gromacs_settings.end_table / gromacs_settings.delta_r);
	table_u = (double*)malloc(sizeof(double) * num_table_entries);
	table_du = (double*)malloc(sizeof(double) * num_table_entries);

	//FILL IN POTENTIAL
	for (i = 0; i < num_table_entries; i++)
	{
		r = (double)i*gromacs_settings.delta_r;

		*(table_u + i) = Spline.getValue(r);

		if (!isfinite(*(table_u + i)) || *(table_u + i) > 1.0e5)
		{
			*(table_u + i) = 1.0e5;
		}
	}

	

	//FIND THE CUTOFF
		//for this potential we just have it and it is the last R value
	*(md_cutoff_pointer) = *(R + num_parameters - 1);


	//MANUALLY FILL IN THE FIRST AND LAST FORCE ELEMENTS AND USE FINITE DIFFERENCE
	*(table_du + 0) = 0.0;
	*(table_du + num_table_entries - 1) = 0.0;
	for (i = 1; i < num_table_entries - 1; i++)
	{
		*(table_du + i) = -1.0*(*(table_u + i + 1) - *(table_u + i - 1)) / (2.0*gromacs_settings.delta_r);
	}

	

	//ASSIGN RETURN DATA
	table_data.array_1 = table_u;
	table_data.array_2 = table_du;
	table_data.num_elements = num_table_entries;

	

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////FINAL STUFF///////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////




	//AS OF NOW ALL THE NECCESSARY FILES FOR THE NEXT STEP SHOULD HAVE BEEN GENERATED
	//IN THE LAST STEP DIRECTORY WITH THE OUT FILE NAME FLAG FOR COPYING (AND RENAMING) TO
	//THE NEXT STEP DIRECTORY

	return table_data;
}

//ACCEPTS ADDRESS OF POTENTIAL PARAMETERS AND EXTRACT THE INDIVIDUAL COMPONENTS
void extract_data_standard(/*input*/ double *potential_parameters, int num_parameters, 
				/*output*/ double *R, double *U, int *State)
{
	//load in the points until the last r-space entry according the the input format
	for (int i = 0; i < num_parameters - 1; i++)
	{
		*(R + i) = *(potential_parameters + 0 + 3 * i);
		*(U + i) = *(potential_parameters + 1 + 3 * i);
		*(State + i) = *(potential_parameters + 2 + 3 * i);
	}

	//load in the last r-space entry manually since technically this is not associated with a step parameter
	*(R + num_parameters - 1) = *(potential_parameters + 0 + 3 * (num_parameters - 1));
	*(U + num_parameters - 1) = 0.0;
	*(State + num_parameters - 1) = 0;
}

//THIS DOES THE INVERSE OF THE ABOVE FUNCTION AND PREPARES THE NEWLY UPDATE PARAMETERS FOR FILE WRITING
void format_data_standard(/*output*/ double *potential_parameters,
	/*input*/int num_parameters, double *R, double *U, int *State)
{
	//set the first points until the last r-space entry according the the input format
	for (int i = 0; i < num_parameters - 1; i++)
	{
		*(potential_parameters + 0 + 3 * i) = *(R + i);
		*(potential_parameters + 1 + 3 * i) = *(U + i);
		*(potential_parameters + 2 + 3 * i) = *(State + i);
	}

	//set the last r-space entry manually since technically this is not associated with a step parameter
	*(potential_parameters + 0 + 3 * (num_parameters - 1)) = *(R + num_parameters - 1);
}