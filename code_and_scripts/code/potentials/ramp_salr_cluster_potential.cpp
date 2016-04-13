///////////////////////////////////////////////////////
#include "../dtypes.h" //need this in every potential file
///////////////////////////////////////////////////////

//INCLUDES FOR THIS POTENTIAL
#include <math.h> //for functions
#include <stdlib.h> //for malloc

//ACTUAL POTENTIAL OPTIMIZER FUNCTION
array_pair_and_num_elements optimize_ramp_salr_cluster_potential(int last_step, array_pair_and_num_elements gr_data,
	double *potential_parameters, double *d_potential_parameters, gromacs_settings_class gromacs_settings,
	double *md_cutoff_pointer, double *gr_cutoff_pointer, double *unscaled_gradient_pointer, double *gr_convergence_pointer)
{
	//OPTIMIZATION PARAMETERS
	double e1 = *(potential_parameters + 0), a1 = *(potential_parameters + 1);
	double e2 = *(potential_parameters + 2), a2 = *(potential_parameters + 3);
	//NON_OPTIMIZED PARAMETERS SET BY USER
	int n1 = (int)(*(potential_parameters + 4) + 0.5); //read in as a double so need to round and type cast
	int n2 = (int)(*(potential_parameters + 5) + 0.5); //read in as a double so need to round and type cast
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
	double grad_e1, grad_a1;
	double grad_e2, grad_a2;
	double grad_magnitude, inverse_grad_magnitude;
	double r_I, r_II;
	double gr_I, gr_II;
	double gr_tgt_I, gr_tgt_II;
	double dude1_I, dude1_II, duda1_I, duda1_II;
	double dude2_I, dude2_II, duda2_I, duda2_II;
	//GR CONVERGENCE
	double gr_convergence;
	//MOMENTUM
	double step_e1, step_a1, step_e2, step_a2;


	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////UPDATING PARAMETERS///////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	//EXTRACT MEMORY ADDRESS FROM GR_DATA
	num_lines_gr = gr_data.num_elements;

	//IF LAST STEP IS NOT ZERO CALCULATE UPDATES (THIS IS FOR INTEGRALS)
	if (last_step != 0)
	{
		//INTEGRAL FOR FINDING THE GRADIENT TOWARDS LARGER PROBABILITY
		grad_e1 = 0.0;
		grad_a1 = 0.0;
		grad_e2 = 0.0;
		grad_a2 = 0.0;
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
			dude1_I = -(1.0 / (exp(pow(-((2.0 + a1 - 2.0 * r_I) / a1), n1))*gromacs_settings.kB_T));
			duda1_I = e2 / (a2*exp(pow(2.0, n2)*pow((-1.0 - a1 - a2 / 2.0 + r_I) / a2, n2))*gromacs_settings.kB_T) +
				(pow(2.0, n2)*(double)n2*pow((-1.0 - a1 - a2 / 2.0 + r_I) / a2, -1 + n2)*(((1.0 + a1 + a2)*e2) / a2 - (e2*r_I) / a2)) /
				(a2*exp(pow(2.0, n2)*pow((-1.0 - a1 - a2 / 2.0 + r_I) / a2, n2))*gromacs_settings.kB_T) +
				(pow(2.0, n1)*e1*(double)n1*pow((-1.0 - a1 / 2.0 + r_I) / a1, -1 + n1)*(-1.0 / (2.0*a1) - (-1.0 - a1 / 2.0 + r_I) / pow(a1, 2))) /
				(exp(pow(2.0, n1)*pow((-1.0 - a1 / 2.0 + r_I) / a1, n1))*gromacs_settings.kB_T);
			dude2_I = (1.0 + a1 + a2 - r_I) / (a2*exp(pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), n2))*gromacs_settings.kB_T);
			duda2_I = (e2 / a2 - ((1.0 + a1 + a2)*e2) / pow(a2, 2) + (e2*r_I) / pow(a2, 2)) /
				(exp(pow(2.0, n2)*pow((-1.0 - a1 - a2 / 2.0 + r_I) / a2, n2))*gromacs_settings.kB_T) -
				(pow(2.0, n2)*(double)n2*pow((-1.0 - a1 - a2 / 2.0 + r_I) / a2, -1 + n2)*(((1.0 + a1 + a2)*e2) / a2 - (e2*r_I) / a2)*
				(-1.0 / (2.0*a2) - (-1.0 - a1 - a2 / 2.0 + r_I) / pow(a2, 2))) /
				(exp(pow(2.0, n2)*pow((-1.0 - a1 - a2 / 2.0 + r_I) / a2, n2))*gromacs_settings.kB_T);

			//second set of potential derivatives
			dude1_I = -(1.0 / (exp(pow(-((2.0 + a1 - 2.0 * r_II) / a1), n1))*gromacs_settings.kB_T));
			duda1_I = e2 / (a2*exp(pow(2.0, n2)*pow((-1.0 - a1 - a2 / 2.0 + r_II) / a2, n2))*gromacs_settings.kB_T) +
				(pow(2.0, n2)*(double)n2*pow((-1.0 - a1 - a2 / 2.0 + r_II) / a2, -1 + n2)*(((1.0 + a1 + a2)*e2) / a2 - (e2*r_II) / a2)) /
				(a2*exp(pow(2.0, n2)*pow((-1.0 - a1 - a2 / 2.0 + r_II) / a2, n2))*gromacs_settings.kB_T) +
				(pow(2.0, n1)*e1*(double)n1*pow((-1.0 - a1 / 2.0 + r_II) / a1, -1 + n1)*(-1.0 / (2.0*a1) - (-1.0 - a1 / 2.0 + r_II) / pow(a1, 2))) /
				(exp(pow(2.0, n1)*pow((-1.0 - a1 / 2.0 + r_II) / a1, n1))*gromacs_settings.kB_T);
			dude2_I = (1.0 + a1 + a2 - r_II) / (a2*exp(pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_II) / a2), n2))*gromacs_settings.kB_T);
			duda2_I = (e2 / a2 - ((1.0 + a1 + a2)*e2) / pow(a2, 2) + (e2*r_II) / pow(a2, 2)) /
				(exp(pow(2.0, n2)*pow((-1.0 - a1 - a2 / 2.0 + r_II) / a2, n2))*gromacs_settings.kB_T) -
				(pow(2.0, n2)*(double)n2*pow((-1.0 - a1 - a2 / 2.0 + r_II) / a2, -1 + n2)*(((1.0 + a1 + a2)*e2) / a2 - (e2*r_II) / a2)*
				(-1.0 / (2.0*a2) - (-1.0 - a1 - a2 / 2.0 + r_II) / pow(a2, 2))) /
				(exp(pow(2.0, n2)*pow((-1.0 - a1 - a2 / 2.0 + r_II) / a2, n2))*gromacs_settings.kB_T);

			//integral_data << setprecision(15) << r_I << "," << dude1_I << "," << duda1_I << "," << dude2_I << "," << duda2_I << "," << r_II << "," << dude1_II << "," << duda1_II << "," << dude2_II << "," << duda2_II << endl;

			//actually do the damn integral
			grad_e1 = grad_e1 + 0.5*((r_I*r_I*dude1_I*(gr_I - gr_tgt_I)) + (r_II*r_II*dude1_II*(gr_II - gr_tgt_II)))*gromacs_settings.delta_r;
			grad_a1 = grad_a1 + 0.5*((r_I*r_I*duda1_I*(gr_I - gr_tgt_I)) + (r_II*r_II*duda1_II*(gr_II - gr_tgt_II)))*gromacs_settings.delta_r;
			grad_e2 = grad_e2 + 0.5*((r_I*r_I*dude2_I*(gr_I - gr_tgt_I)) + (r_II*r_II*dude2_II*(gr_II - gr_tgt_II)))*gromacs_settings.delta_r;
			grad_a2 = grad_a2 + 0.5*((r_I*r_I*duda2_I*(gr_I - gr_tgt_I)) + (r_II*r_II*duda2_II*(gr_II - gr_tgt_II)))*gromacs_settings.delta_r;

			gr_convergence = gr_convergence + 0.5*((r_I*r_I*(gr_I - gr_tgt_I)*(gr_I - gr_tgt_I)) + (r_II*r_II*(gr_II - gr_tgt_II)*(gr_II - gr_tgt_II)))*gromacs_settings.delta_r;
		}
		//prepare to pass back the convergence criterion
		*(gr_convergence_pointer) = gr_convergence;

		//scale the e's by the appropriate kBT factor
		grad_e1 = grad_e1*(gromacs_settings.kB_T)*(gromacs_settings.kB_T);
		grad_e2 = grad_e2*(gromacs_settings.kB_T)*(gromacs_settings.kB_T);

		//calculate unscaled gradient - need to get rid of one power of kBT to get dimensionless gradient
		grad_magnitude = sqrt((grad_e1*grad_e1) / (gromacs_settings.kB_T*gromacs_settings.kB_T)
			+ (grad_a1*grad_a1) + (grad_e2*grad_e2) / (gromacs_settings.kB_T*gromacs_settings.kB_T) + (grad_a2*grad_a2));
		*unscaled_gradient_pointer = grad_magnitude;

		//scale the gradient
		grad_e1 = gromacs_settings.gradient_scale*grad_e1;
		grad_a1 = gromacs_settings.gradient_scale*grad_a1;
		grad_e2 = gromacs_settings.gradient_scale*grad_e2;
		grad_a2 = gromacs_settings.gradient_scale*grad_a2;

		//add in the momentum contribution first and then reset the step for passing back via the pointer
		step_e1 = grad_e1 + gromacs_settings.momentum_scale * (*(d_potential_parameters + 0)); *(d_potential_parameters + 0) = step_e1;
		step_a1 = grad_a1 + gromacs_settings.momentum_scale * (*(d_potential_parameters + 1)); *(d_potential_parameters + 1) = step_a1;
		step_e2 = grad_e2 + gromacs_settings.momentum_scale * (*(d_potential_parameters + 2)); *(d_potential_parameters + 2) = step_e2;
		step_a2 = grad_a2 + gromacs_settings.momentum_scale * (*(d_potential_parameters + 3)); *(d_potential_parameters + 3) = step_a2;

		//scale the gradient vector so it is not too large of a jump
		/*grad_magnitude = sqrt((grad_e1*grad_e1) + (grad_a1*grad_a1) + (grad_e2*grad_e2) + (grad_a2*grad_a2));
		if (grad_magitude)
		{
		inverse_grad_magnitude = 1.0 / sqrt((grad_e1*grad_e1) + (grad_a1*grad_a1) + (grad_e2*grad_e2) + (grad_a2*grad_a2));
		grad_e1 = grad_e1*inverse_grad_magnitude*gromacs_settings.gradient_scale;
		grad_a1 = grad_a1*inverse_grad_magnitude*gromacs_settings.gradient_scale;
		grad_e2 = grad_e2*inverse_grad_magnitude*gromacs_settings.gradient_scale;
		grad_a2 = grad_a2*inverse_grad_magnitude*gromacs_settings.gradient_scale;
		}

		//pass back the unscaled gradient to asess convergence
		*unscaled_gradient_pointer = 1.0 / inverse_grad_magnitude;*/

		//update the parameters and pass back
		e1 = e1 + /*grad*/ step_e1; *(potential_parameters + 0) = e1;
		a1 = a1 + /*grad*/ step_a1; *(potential_parameters + 1) = a1; //SALR SO KEEPING SHORT RANGED
		e2 = e2 + /*grad*/ step_e2; *(potential_parameters + 2) = e2;
		a2 = a2 + /*grad*/ step_a2; *(potential_parameters + 3) = a2;
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

		//"HARDCORE" WCA PORTION (WE WANT THIS TO SCALE WITH TEMPERATURE SUCH THAT RESCALING T TRIVIALLY AFFECTS OPTIMIZED PARAMETERS)
		if (r <= pow(2.0, 1.0 / 10.0))
		{
			//effective hardcore WCA portion
			*(table_u + i) = gromacs_settings.kB_T*(4.0*(pow(1.0 / r, 20) - pow(1.0 / r, 10)) + 1.0);

			//optimized portion
			*(table_u + i) = *(table_u + i) +
				-(e1 / exp(pow(2.0, n1)*pow((-1.0 - a1 / 2.0 + r) / a1, n1))) +
				(((1.0 + a1 + a2)*e2) / a2 - (e2*r) / a2) / exp(pow(2.0, n2)*pow((-1.0 - a1 - a2 / 2.0 + r) / a2, n2));

			//make sure there is no effective divergence and if so scale back to something GROMACS can handle
			if (!isfinite(*(table_u + i)) || *(table_u + i) > 1.0e5*gromacs_settings.kB_T)
			{
				*(table_u + i) = 1.0e5*gromacs_settings.kB_T;
			}
		}
		else
		{
			//non-existent WCA portion
			*(table_u + i) = 0.0;

			//optimized portion
			*(table_u + i) = *(table_u + i) +
				-(e1 / exp(pow(2.0, n1)*pow((-1.0 - a1 / 2.0 + r) / a1, n1))) +
				(((1.0 + a1 + a2)*e2) / a2 - (e2*r) / a2) / exp(pow(2.0, n2)*pow((-1.0 - a1 - a2 / 2.0 + r) / a2, n2));
		}
	}

	//FIND THE CUTOFF
	md_cutoff_index = 0; //location in array of the last finite point (cutoff)
	u_at_cutoff = 0.0; //magnitude at cutoff to shift potential by
	for (i = num_table_entries - 2; i > 0; i--)
	{
		r = (double)i*gromacs_settings.delta_r;
		dudr = (*(table_u + i + 1) - *(table_u + i - 1)) / (2.0*gromacs_settings.delta_r);

		if (fabs(dudr / gromacs_settings.kB_T) > gromacs_settings.md_cutoff_magnitude)
		{
			md_cutoff_index = i;
			u_at_cutoff = *(table_u + md_cutoff_index);
			*(md_cutoff_pointer) = r;
			break;
		}
	}

	//EMPLOY THE CUTOFF
	for (i = 0; i <= md_cutoff_index; i++)
	{
		*(table_u + i) = *(table_u + i) - u_at_cutoff;
	}
	for (i = md_cutoff_index + 1; i < num_table_entries; i++)
	{
		*(table_u + i) = 0.0;
	}

	//MANUALLY FILL IN THE FIRST AND LAST FORCE ELEMENTS AND USE FINITE DIFFERENCE ON CUT AND SHIFTED POTENTIAL
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





	/*((e2*(2.0 + 2.0 * a1 + a2 - 2.0 * (double)n2*pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), n2)*(1.0 + a1 + a2 - r_I) - 2.0 * r_I)) /
	(a2*exp(pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), n2))*(2.0 + 2.0 * a1 + a2 - 2.0 * r_I)) -
	(2.0 * e1*(double)n1*pow(-((2.0 + a1 - 2.0 * r_I) / a1), -1 + n1)*(-1.0 + r_I)) /
	(pow(a1, 2)*exp(pow(-((2.0 + a1 - 2.0 * r_I) / a1), n1)))) / gromacs_settings.kB_T;*/

	/*(e2*(1.0 + a1 - r_I)*(-2.0 - 2.0 * a1 - a2 + 2.0 * (double)n2*pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), n2)*(1.0 + a1 + a2 - r_I) + 2.0 * r_I)) /
	(pow(a2, 2)*exp(pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), n2))*gromacs_settings.kB_T*(2.0 + 2.0 * a1 + a2 - 2.0 * r_I));*/

	//AS OF NOW ALL THE NECCESSARY FILES FOR THE NEXT STEP SHOULD HAVE BEEN GENERATED
	//IN THE LAST STEP DIRECTORY WITH THE OUT FILE NAME FLAG FOR COPYING (AND RENAMING) TO
	//THE NEXT STEP DIRECTORY

	return table_data;
}

