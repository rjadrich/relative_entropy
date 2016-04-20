///////////////////////////////////////////////////////
#include "../dtypes.h" //need this in every potential file
#include "potentials.h"//need this in every potential file
///////////////////////////////////////////////////////

//INCLUDES FOR THIS POTENTIAL
#include <math.h> //for functions
#include <stdlib.h> //for malloc

double sech(double x);

//ACTUAL POTENTIAL OPTIMIZER FUNCTION
array_pair_and_num_elements potential_data::optimize_crystal_potential(int last_step, array_pair_and_num_elements gr_data,
	double *potential_parameters, double *d_potential_parameters, gromacs_settings_class gromacs_settings,
	double *md_cutoff_pointer, double *unscaled_gradient_pointer, double *gr_convergence_pointer)
{
	//OPTIMIZATION PARAMETERS
	double A = *(potential_parameters + 0); double n = *(potential_parameters + 1);
	double L1 = *(potential_parameters + 2); double k1 = *(potential_parameters + 3); double d1 = *(potential_parameters + 4);
	double L2 = *(potential_parameters + 5); double k2 = *(potential_parameters + 6); double d2 = *(potential_parameters + 7);
	double rcut = *(potential_parameters + 8);
	//NON_OPTIMIZED PARAMETERS SET BY USER
			//none right now
	//STORAGE POTENTIAL AND DERIVATIVE (ALL ARE IN UNITS OF KBT)
	double u, dudr;
	//GENERIC DISTANCE VARIABLE
	double r;
	//GENERIC LOOP INTEGER
	int i;
	//TABLE CREATION PARAMETERS
	array_pair_and_num_elements table_data; //stores the generated table data to pass back
	int md_cutoff_index; //location in array of the last finite point (cutoff)
	double u_at_cutoff; //value of potential at cutoff to shift it with
	double *table_u, *table_du; //table arrays
	int num_table_entries; //number of elements in table arrays
	//RDF STORAGE
	int num_lines_gr; //again, better than using gr_data directly
	//GRADIENT CALCULATION
	double grad_A, grad_n;
	double grad_L1, grad_k1, grad_d1;
	double grad_L2, grad_k2, grad_d2;
	double grad_rcut;
		/////////////////////////////////////////////////////////////////
	double grad_magnitude, inverse_grad_magnitude;
	double r_I, r_II;
	double gr_I, gr_II;
	double gr_tgt_I, gr_tgt_II;
		/////////////////////////////////////////////////////////////////
	double dudA_I, dudA_II, dudn_I, dudn_II;
	double dudL1_I, dudL1_II, dudk1_I, dudk1_II, dudd1_I, dudd1_II;
	double dudL2_I, dudL2_II, dudk2_I, dudk2_II, dudd2_I, dudd2_II;
	double dudrcut_I, dudrcut_II;
	//GR CONVERGENCE
	double gr_convergence;
	//MOMENTUM
	double step_A, step_n;
	double step_L1, step_k1, step_d1;
	double step_L2, step_k2, step_d2;
	double step_rcut;
	//CONSTRAINED AMPLITUDE
		//none right now
	//SOME OTHER VARIABLES USED IN POTENTIAL
	double P, Q, R;

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////UPDATING PARAMETERS///////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	//EXTRACT MEMORY ADDRESS FROM GR_DATA
	num_lines_gr = gr_data.num_elements;

	//IF LAST STEP IS NOT ZERO CALCULATE UPDATES (THIS IS FOR INTEGRALS)
	if (last_step != 0)
	{
		//INTEGRAL FOR FINDING THE GRADIENT TOWARDS LARGER PROBABILITY
		grad_A = 0.0;
		grad_n = 0.0;
		grad_L1 = 0.0;
		grad_k1 = 0.0;
		grad_d1 = 0.0;
		grad_L2 = 0.0;
		grad_k2 = 0.0;
		grad_d2 = 0.0;
		grad_rcut = 0.0;
		gr_convergence = 0.0;
		for (i = 0; i < num_lines_gr - 1; i++)
		{
			//using trapezoids so need  two values
			r_I = (double)i*gromacs_settings.delta_r;
			r_II = (double)(i + 1)*gromacs_settings.delta_r;

			//fill in rdf data
			gr_I = *(gr_data.array_1 + i);
			gr_II = *(gr_data.array_1 + i + 1);
			gr_tgt_I = *(gr_data.array_2 + i);
			gr_tgt_II = *(gr_data.array_2 + i + 1);

			//first set of potential derivatives
			dudA_I = -((1.0 + n)*(2.0 + n)) / (2.0*pow(rcut, n)) + n*(2.0 + n)*pow(rcut, -1.0 - n)*r_I - (n*(1.0 + n)*pow(rcut, -2.0 - n)*pow(r_I, 2.0)) / 2.0 + pow(r_I, -n);
			dudn_I = (pow(rcut, -2.0 - n)*(A*pow(r_I, n)*((-rcut + r_I)*((3.0 + 2.0 * n)*rcut - (1.0 + 2.0 * n)*r_I) +
				((1.0 + n)*(2.0 + n)*pow(rcut, 2.0) - 2.0 * n*(2.0 + n)*rcut*r_I + n*(1.0 + n)*pow(r_I, 2.0))*log(rcut)) - 2.0 * A*pow(rcut, 2.0 + n)*log(r_I))) / (2.0*pow(r_I, n));
			dudL1_I = -tanh(k1*(d1 - rcut)) + k1*(rcut - r_I)*pow(sech(k1*(-d1 + rcut)), 2.0)*(-1.0 + k1*(rcut - r_I)*tanh(k1*(d1 - rcut))) - tanh(k1*(-d1 + r_I));
			dudk1_I = (L1*(6.0 * pow(k1, 2.0)*(d1 - rcut)*pow(rcut - r_I, 2.0)*pow(sech(k1*(d1 - rcut)), 4.0) + 2.0 * (d1 - r_I)*pow(sech(k1*(-d1 + r_I)), 2.0) +
				2.0 * pow(sech(k1*(d1 - rcut)), 2.0)*(-(d1*(1.0 + 2.0 * pow(k1, 2.0)*pow(rcut - r_I, 2.0))) + 2.0 * pow(k1, 2.0)*rcut*pow(rcut - r_I, 2.0) + r_I +
				2.0 * k1*(d1 - r_I)*(rcut - r_I)*tanh(k1*(d1 - rcut))))) / 2.0;
			dudd1_I = k1*L1*(-pow(sech(k1*(d1 - rcut)), 2.0) + pow(sech(k1*(-d1 + r_I)), 2.0) -
				k1*(rcut - r_I)*pow(sech(k1*(-d1 + rcut)), 4.0)*(k1*(rcut - r_I)*(-2.0 + cosh(2.0 * k1*(-d1 + rcut))) + sinh(2.0 * k1*(-d1 + rcut))));
			dudL2_I = -tanh(k2*(d2 - rcut)) + k2*(rcut - r_I)*pow(sech(k2*(-d2 + rcut)), 2.0)*(-1.0 + k2*(rcut - r_I)*tanh(k2*(d2 - rcut))) - tanh(k2*(-d2 + r_I));
			dudk2_I = (L2*(6.0 * pow(k2, 2.0)*(d2 - rcut)*pow(rcut - r_I, 2.0)*pow(sech(k2*(d2 - rcut)), 4.0) + 2.0 * (d2 - r_I)*pow(sech(k2*(-d2 + r_I)), 2.0) +
				2.0 * pow(sech(k2*(d2 - rcut)), 2.0)*(-(d2*(1.0 + 2.0 * pow(k2, 2.0)*pow(rcut - r_I, 2.0))) + 2.0 * pow(k2, 2.0)*rcut*pow(rcut - r_I, 2.0) + r_I +
				2.0 * k2*(d2 - r_I)*(rcut - r_I)*tanh(k2*(d2 - rcut))))) / 2.0;
			dudd2_I = k2*L2*(-pow(sech(k2*(d2 - rcut)), 2.0) + pow(sech(k2*(-d2 + r_I)), 2.0) -
				k2*(rcut - r_I)*pow(sech(k2*(-d2 + rcut)), 4.0)*(k2*(rcut - r_I)*(-2.0 + cosh(2.0 * k2*(-d2 + rcut))) + sinh(2.0 * k2*(-d2 + rcut))));
			dudrcut_I = (2.0 * A*n*pow(rcut, -1.0 - n) + 3.0 * A*pow(n, 2.0)*pow(rcut, -1.0 - n) + A*pow(n, 3.0)*pow(rcut, -1.0 - n) - 4.0 * A*n*pow(rcut, -2.0 - n)*r_I -
				6.0 * A*pow(n, 2.0)*pow(rcut, -2.0 - n)*r_I - 2.0 * A*pow(n, 3.0)*pow(rcut, -2.0 - n)*r_I + 2.0 * A*n*pow(rcut, -3.0 - n)*pow(r_I, 2.0) +
				3.0 * A*pow(n, 2.0)*pow(rcut, -3.0 - n)*pow(r_I, 2.0) + A*pow(n, 3.0)*pow(rcut, -3.0 - n)*pow(r_I, 2.0) + 2.0 * k1*L1*pow(sech(k1*(d1 - rcut)), 2.0) +
				2.0 * k2*L2*pow(sech(k2*(d2 - rcut)), 2.0) - 2.0 * k1*L1*pow(sech(k1*(-d1 + rcut)), 2.0) - 2.0 * pow(k1, 3.0)*L1*pow(rcut, 2.0)*pow(sech(k1*(-d1 + rcut)), 4.0) +
				4.0 * pow(k1, 3.0)*L1*rcut*r_I*pow(sech(k1*(-d1 + rcut)), 4.0) - 2.0 * pow(k1, 3.0)*L1*pow(r_I, 2.0)*pow(sech(k1*(-d1 + rcut)), 4.0) -
				2.0 * k2*L2*pow(sech(k2*(-d2 + rcut)), 2.0) - 2.0 * pow(k2, 3.0)*L2*pow(rcut, 2.0)*pow(sech(k2*(-d2 + rcut)), 4.0) +
				4.0 * pow(k2, 3.0)*L2*rcut*r_I*pow(sech(k2*(-d2 + rcut)), 4.0) - 2.0 * pow(k2, 3.0)*L2*pow(r_I, 2.0)*pow(sech(k2*(-d2 + rcut)), 4.0) +
				4.0 * pow(k1, 3.0)*L1*pow(rcut, 2.0)*pow(sech(k1*(-d1 + rcut)), 2.0)*pow(tanh(k1*(-d1 + rcut)), 2.0) -
				8.0 * pow(k1, 3.0)*L1*rcut*r_I*pow(sech(k1*(-d1 + rcut)), 2.0)*pow(tanh(k1*(-d1 + rcut)), 2.0) +
				4.0 * pow(k1, 3.0)*L1*pow(r_I, 2.0)*pow(sech(k1*(-d1 + rcut)), 2.0)*pow(tanh(k1*(-d1 + rcut)), 2.0) +
				4.0 * pow(k2, 3.0)*L2*pow(rcut, 2.0)*pow(sech(k2*(-d2 + rcut)), 2.0)*pow(tanh(k2*(-d2 + rcut)), 2.0) -
				8.0 * pow(k2, 3.0)*L2*rcut*r_I*pow(sech(k2*(-d2 + rcut)), 2.0)*pow(tanh(k2*(-d2 + rcut)), 2.0) +
				4.0 * pow(k2, 3.0)*L2*pow(r_I, 2.0)*pow(sech(k2*(-d2 + rcut)), 2.0)*pow(tanh(k2*(-d2 + rcut)), 2.0)) / 2.0;

			//second set of potential derivatives
			dudA_II = -((1.0 + n)*(2.0 + n)) / (2.0*pow(rcut, n)) + n*(2.0 + n)*pow(rcut, -1.0 - n)*r_II - (n*(1.0 + n)*pow(rcut, -2.0 - n)*pow(r_II, 2.0)) / 2.0 + pow(r_II, -n);
			dudn_II = (pow(rcut, -2.0 - n)*(A*pow(r_II, n)*((-rcut + r_II)*((3.0 + 2.0 * n)*rcut - (1.0 + 2.0 * n)*r_II) +
				((1.0 + n)*(2.0 + n)*pow(rcut, 2.0) - 2.0 * n*(2.0 + n)*rcut*r_II + n*(1.0 + n)*pow(r_II, 2.0))*log(rcut)) - 2.0 * A*pow(rcut, 2.0 + n)*log(r_II))) / (2.0*pow(r_II, n));
			dudL1_II = -tanh(k1*(d1 - rcut)) + k1*(rcut - r_II)*pow(sech(k1*(-d1 + rcut)), 2.0)*(-1.0 + k1*(rcut - r_II)*tanh(k1*(d1 - rcut))) - tanh(k1*(-d1 + r_II));
			dudk1_II = (L1*(6.0 * pow(k1, 2.0)*(d1 - rcut)*pow(rcut - r_II, 2.0)*pow(sech(k1*(d1 - rcut)), 4.0) + 2.0 * (d1 - r_II)*pow(sech(k1*(-d1 + r_II)), 2.0) +
				2.0 * pow(sech(k1*(d1 - rcut)), 2.0)*(-(d1*(1.0 + 2.0 * pow(k1, 2.0)*pow(rcut - r_II, 2.0))) + 2.0 * pow(k1, 2.0)*rcut*pow(rcut - r_II, 2.0) + r_II +
				2.0 * k1*(d1 - r_II)*(rcut - r_II)*tanh(k1*(d1 - rcut))))) / 2.0;
			dudd1_II = k1*L1*(-pow(sech(k1*(d1 - rcut)), 2.0) + pow(sech(k1*(-d1 + r_II)), 2.0) -
				k1*(rcut - r_II)*pow(sech(k1*(-d1 + rcut)), 4.0)*(k1*(rcut - r_II)*(-2.0 + cosh(2.0 * k1*(-d1 + rcut))) + sinh(2.0 * k1*(-d1 + rcut))));
			dudL2_II = -tanh(k2*(d2 - rcut)) + k2*(rcut - r_II)*pow(sech(k2*(-d2 + rcut)), 2.0)*(-1.0 + k2*(rcut - r_II)*tanh(k2*(d2 - rcut))) - tanh(k2*(-d2 + r_II));
			dudk2_II = (L2*(6.0 * pow(k2, 2.0)*(d2 - rcut)*pow(rcut - r_II, 2.0)*pow(sech(k2*(d2 - rcut)), 4.0) + 2.0 * (d2 - r_II)*pow(sech(k2*(-d2 + r_II)), 2.0) +
				2.0 * pow(sech(k2*(d2 - rcut)), 2.0)*(-(d2*(1.0 + 2.0 * pow(k2, 2.0)*pow(rcut - r_II, 2.0))) + 2.0 * pow(k2, 2.0)*rcut*pow(rcut - r_II, 2.0) + r_II +
				2.0 * k2*(d2 - r_II)*(rcut - r_II)*tanh(k2*(d2 - rcut))))) / 2.0;
			dudd2_II = k2*L2*(-pow(sech(k2*(d2 - rcut)), 2.0) + pow(sech(k2*(-d2 + r_II)), 2.0) -
				k2*(rcut - r_II)*pow(sech(k2*(-d2 + rcut)), 4.0)*(k2*(rcut - r_II)*(-2.0 + cosh(2.0 * k2*(-d2 + rcut))) + sinh(2.0 * k2*(-d2 + rcut))));
			dudrcut_II = (2.0 * A*n*pow(rcut, -1.0 - n) + 3.0 * A*pow(n, 2.0)*pow(rcut, -1.0 - n) + A*pow(n, 3.0)*pow(rcut, -1.0 - n) - 4.0 * A*n*pow(rcut, -2.0 - n)*r_II -
				6.0 * A*pow(n, 2.0)*pow(rcut, -2.0 - n)*r_II - 2.0 * A*pow(n, 3.0)*pow(rcut, -2.0 - n)*r_II + 2.0 * A*n*pow(rcut, -3.0 - n)*pow(r_II, 2.0) +
				3.0 * A*pow(n, 2.0)*pow(rcut, -3.0 - n)*pow(r_II, 2.0) + A*pow(n, 3.0)*pow(rcut, -3.0 - n)*pow(r_II, 2.0) + 2.0 * k1*L1*pow(sech(k1*(d1 - rcut)), 2.0) +
				2.0 * k2*L2*pow(sech(k2*(d2 - rcut)), 2.0) - 2.0 * k1*L1*pow(sech(k1*(-d1 + rcut)), 2.0) - 2.0 * pow(k1, 3.0)*L1*pow(rcut, 2.0)*pow(sech(k1*(-d1 + rcut)), 4.0) +
				4.0 * pow(k1, 3.0)*L1*rcut*r_II*pow(sech(k1*(-d1 + rcut)), 4.0) - 2.0 * pow(k1, 3.0)*L1*pow(r_II, 2.0)*pow(sech(k1*(-d1 + rcut)), 4.0) -
				2.0 * k2*L2*pow(sech(k2*(-d2 + rcut)), 2.0) - 2.0 * pow(k2, 3.0)*L2*pow(rcut, 2.0)*pow(sech(k2*(-d2 + rcut)), 4.0) +
				4.0 * pow(k2, 3.0)*L2*rcut*r_II*pow(sech(k2*(-d2 + rcut)), 4.0) - 2.0 * pow(k2, 3.0)*L2*pow(r_II, 2.0)*pow(sech(k2*(-d2 + rcut)), 4.0) +
				4.0 * pow(k1, 3.0)*L1*pow(rcut, 2.0)*pow(sech(k1*(-d1 + rcut)), 2.0)*pow(tanh(k1*(-d1 + rcut)), 2.0) -
				8.0 * pow(k1, 3.0)*L1*rcut*r_II*pow(sech(k1*(-d1 + rcut)), 2.0)*pow(tanh(k1*(-d1 + rcut)), 2.0) +
				4.0 * pow(k1, 3.0)*L1*pow(r_II, 2.0)*pow(sech(k1*(-d1 + rcut)), 2.0)*pow(tanh(k1*(-d1 + rcut)), 2.0) +
				4.0 * pow(k2, 3.0)*L2*pow(rcut, 2.0)*pow(sech(k2*(-d2 + rcut)), 2.0)*pow(tanh(k2*(-d2 + rcut)), 2.0) -
				8.0 * pow(k2, 3.0)*L2*rcut*r_II*pow(sech(k2*(-d2 + rcut)), 2.0)*pow(tanh(k2*(-d2 + rcut)), 2.0) +
				4.0 * pow(k2, 3.0)*L2*pow(r_II, 2.0)*pow(sech(k2*(-d2 + rcut)), 2.0)*pow(tanh(k2*(-d2 + rcut)), 2.0)) / 2.0;

			//actually do the damn integral
			grad_A = grad_A + 0.5*((r_I*r_I*dudA_I*(gr_I - gr_tgt_I)) + (r_II*r_II*dudA_II*(gr_II - gr_tgt_II)))*gromacs_settings.delta_r;
			grad_n = grad_n + 0.5*((r_I*r_I*dudn_I*(gr_I - gr_tgt_I)) + (r_II*r_II*dudn_II*(gr_II - gr_tgt_II)))*gromacs_settings.delta_r;
			grad_L1 = grad_L1 + 0.5*((r_I*r_I*dudL1_I*(gr_I - gr_tgt_I)) + (r_II*r_II*dudL1_II*(gr_II - gr_tgt_II)))*gromacs_settings.delta_r;
			grad_k1 = grad_k1 + 0.5*((r_I*r_I*dudk1_I*(gr_I - gr_tgt_I)) + (r_II*r_II*dudk1_II*(gr_II - gr_tgt_II)))*gromacs_settings.delta_r;
			grad_d1 = grad_d1 + 0.5*((r_I*r_I*dudd1_I*(gr_I - gr_tgt_I)) + (r_II*r_II*dudd1_II*(gr_II - gr_tgt_II)))*gromacs_settings.delta_r;
			grad_L2 = grad_L2 + 0.5*((r_I*r_I*dudL2_I*(gr_I - gr_tgt_I)) + (r_II*r_II*dudL2_II*(gr_II - gr_tgt_II)))*gromacs_settings.delta_r;
			grad_k2 = grad_k2 + 0.5*((r_I*r_I*dudk2_I*(gr_I - gr_tgt_I)) + (r_II*r_II*dudk2_II*(gr_II - gr_tgt_II)))*gromacs_settings.delta_r;
			grad_d2 = grad_d2 + 0.5*((r_I*r_I*dudd2_I*(gr_I - gr_tgt_I)) + (r_II*r_II*dudd2_II*(gr_II - gr_tgt_II)))*gromacs_settings.delta_r;
			grad_rcut = grad_rcut + 0.5*((r_I*r_I*dudrcut_I*(gr_I - gr_tgt_I)) + (r_II*r_II*dudrcut_II*(gr_II - gr_tgt_II)))*gromacs_settings.delta_r;

			gr_convergence = gr_convergence + 0.5*((r_I*r_I*(gr_I - gr_tgt_I)*(gr_I - gr_tgt_I)) + (r_II*r_II*(gr_II - gr_tgt_II)*(gr_II - gr_tgt_II)))*gromacs_settings.delta_r;
		}
		//prepare to pass back the convergence criterion
		*(gr_convergence_pointer) = gr_convergence;

		//calculate unscaled gradient
		grad_magnitude = sqrt((grad_A*grad_A) + (grad_n*grad_n) 
			+ (grad_L1*grad_L1) + (grad_k1*grad_k1) + (grad_d1*grad_d1)
			+ (grad_L2*grad_L2) + (grad_k2*grad_k2) + (grad_d2*grad_d2)
			+ (grad_rcut*grad_rcut));
		*unscaled_gradient_pointer = grad_magnitude;

		//scale the gradient
		grad_A = gromacs_settings.gradient_scale*grad_A;
		grad_n = gromacs_settings.gradient_scale*grad_n;
		grad_L1 = gromacs_settings.gradient_scale*grad_L1;
		grad_k1 = gromacs_settings.gradient_scale*grad_k1;
		grad_d1 = gromacs_settings.gradient_scale*grad_d1;
		grad_L2 = gromacs_settings.gradient_scale*grad_L2;
		grad_k2 = gromacs_settings.gradient_scale*grad_k2;
		grad_d2 = gromacs_settings.gradient_scale*grad_d2;
		grad_rcut = gromacs_settings.gradient_scale*grad_rcut;

		//add in the momentum contribution first and then reset the step for passing back via the pointer
		step_A = grad_A + gromacs_settings.momentum_scale * (*(d_potential_parameters + 0)); *(d_potential_parameters + 0) = step_A;
		step_n = grad_n + gromacs_settings.momentum_scale * (*(d_potential_parameters + 1)); *(d_potential_parameters + 1) = step_n;
		step_L1 = grad_L1 + gromacs_settings.momentum_scale * (*(d_potential_parameters + 2)); *(d_potential_parameters + 2) = step_L1;
		step_k1 = grad_k1 + gromacs_settings.momentum_scale * (*(d_potential_parameters + 3)); *(d_potential_parameters + 3) = step_k1;
		step_d1 = grad_d1 + gromacs_settings.momentum_scale * (*(d_potential_parameters + 4)); *(d_potential_parameters + 4) = step_d1;
		step_L2 = grad_L2 + gromacs_settings.momentum_scale * (*(d_potential_parameters + 5)); *(d_potential_parameters + 5) = step_L2;
		step_k2 = grad_k2 + gromacs_settings.momentum_scale * (*(d_potential_parameters + 6)); *(d_potential_parameters + 6) = step_k2;
		step_d2 = grad_d2 + gromacs_settings.momentum_scale * (*(d_potential_parameters + 7)); *(d_potential_parameters + 7) = step_d2;
		step_rcut = grad_rcut + gromacs_settings.momentum_scale * (*(d_potential_parameters + 8)); *(d_potential_parameters + 8) = step_rcut;

		
		//any constraints desired...
		

		//update the parameters and pass back
		A = A + step_A; *(potential_parameters + 0) = A;
		n = n + step_n; *(potential_parameters + 1) = n;
		L1 = L1 + step_L1; *(potential_parameters + 2) = L1;
		k1 = k1 + step_k1; *(potential_parameters + 3) = k1;
		d1 = d1 + step_d1; *(potential_parameters + 4) = d1;
		L2 = L2 + step_L2; *(potential_parameters + 5) = L2;
		k2 = k2 + step_k2; *(potential_parameters + 6) = k2;
		d2 = d2 + step_d2; *(potential_parameters + 7) = d2;
		rcut = rcut + step_rcut; *(potential_parameters + 8) = rcut;

	}


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

		if (r <= rcut)
		{
			P = -(A*n*(1.0 + n)*pow(rcut, -2.0 - n)) / 2.0 - pow(k1, 2.0)*L1*pow(sech(k1*(-d1 + rcut)), 2.0)*tanh(k1*(-d1 + rcut)) -
				pow(k2, 2.0)*L2*pow(sech(k2*(-d2 + rcut)), 2.0)*tanh(k2*(-d2 + rcut));
			Q = pow(rcut, -1.0 - n)*(A*n*(2.0 + n) + pow(rcut, 1.0 + n)*(k1*L1*pow(sech(k1*(-d1 + rcut)), 2.0)*(1.0 + 2.0 * k1*rcut*tanh(k1*(-d1 + rcut))) +
				k2*L2*pow(sech(k2*(-d2 + rcut)), 2.0)*(1.0 + 2.0 * k2*rcut*tanh(k2*(-d2 + rcut)))));
			R = (-(A*(1.0 + n)*(2.0 + n)) - 2.0 * pow(rcut, n)*(L1 + L2 + L1*tanh(k1*(d1 - rcut)) + L2*tanh(k2*(d2 - rcut)) +
				k1*L1*rcut*pow(sech(k1*(-d1 + rcut)), 2.0)*(1.0 + k1*rcut*tanh(k1*(-d1 + rcut))) +
				k2*L2*rcut*pow(sech(k2*(-d2 + rcut)), 2.0)*(1.0 + k2*rcut*tanh(k2*(-d2 + rcut))))) / (2.0*pow(rcut, n));
			*(table_u + i) = A / pow(r, n) + L1*(1.0 - tanh(k1*(-d1 + r))) + L2*(1.0 - tanh(k2*(-d2 + r))) + P*r*r + Q*r + R;

			if (!isfinite(*(table_u + i)) || *(table_u + i) > 1.0e5)
			{
				*(table_u + i) = 1.0e5;
			}
		}
		else
		{
			*(table_u + i) = 0.0;
		}
	}

	//FIND THE CUTOFF
		//for this potential we just have it
	*(md_cutoff_pointer) = rcut;

	//FIND THE CUTOFF
	/*md_cutoff_index = 0; //location in array of the last finite point (cutoff)
	u_at_cutoff = 0.0; //magnitude at cutoff to shift potential by
	for (i = num_table_entries - 2; i > 0; i--)
	{
		r = (double)i*gromacs_settings.delta_r;
		dudr = (*(table_u + i + 1) - *(table_u + i - 1)) / (2.0*gromacs_settings.delta_r);

		if (fabs(dudr) > gromacs_settings.md_cutoff_magnitude)
		{
			md_cutoff_index = i;
			u_at_cutoff = *(table_u + md_cutoff_index);
			*(md_cutoff_pointer) = r;
			break;
		}
	}*/

	//EMPLOY THE CUTOFF
	/*for (i = 0; i <= md_cutoff_index; i++)
	{
		*(table_u + i) = *(table_u + i) - u_at_cutoff;
	}
	for (i = md_cutoff_index + 1; i < num_table_entries; i++)
	{
		*(table_u + i) = 0.0;
	}*/

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


double sech(double x)
{
	return 1.0 / cosh(x);
}