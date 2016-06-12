#ifndef DTYPES_H
#define DTYPES_H

//CLASSES FOR COMPACT DESCRIPTIONS OF COMMONLY NEEDED DATA
class gromacs_settings_class
{
public:
	double delta_r; //spacing for g(r) and u(r) to integrate 
	int num_gr_calcs; //number of parallelized gr calcs
	double equil_time; //equilibration time for calculating g(r)
	double final_time; //last time to use for calculating g(r)
	double kB_T; //need this for gromacs to create appropriately scaled table file with u(r) and -du(r)/dr
	double end_table; //extends the table by specified amount to ensure that gromacs does not throw a fit
	double gradient_scale; //amplitude to scale back the step in the direction of the gradient by to stabilize iterative scheme
	double momentum_scale; //amplitude to scale back the momentum part of the update by (set to zero for normal gradient descent)
	double md_cutoff_magnitude; //magnitude where the potential can be truncated and set to zero
	double buffer_size; //amount to buffer beyond the cutoff to allow for infrequent neighbor table updates
	int rlist; //line number for rlist
	int rcoulomb; //line number for rcoulomb
	int rvdw; //line number for rvdw
	int dimensions; //line number for selecting if 1D, 2D, or 3D
	int init_conf; //-1 for previous initial and 1 for previous output
};

class array_pair_and_num_elements
{
public:
	int num_elements;
	double *array_1;
	double *array_2;
};

#endif