// relative_entropy.cpp : Defines the entry point for the console application.
//

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctime>

using namespace std;

const int rdf_gromacs_skip = 0; //number of lines to skip in gromacs RDF file (seems to be 13 for default RDF settings)

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
	double gr_cutoff_magnitude; //magnitude determining where potential derivatives with respect to parameters is small
	double buffer_size; //amount to buffer beyond the cutoff to allow for infrequent neighbor table updates
	int rlist; //line number for rlist
	int rcoulomb; //line number for rcoulomb
	int rvdw; //line number for rvdw
};

class array_pair_and_num_elements
{
public:
	int num_elements;
	double *array_1;
	double *array_2;
};

//LOG FILE TO WRITE OUT CURRENT STEP INFO AND ANY ERRORS ENCOUNTERED IN ANY FUNCTION
ofstream log_filestream("log.txt", ios::app);

//GENERAL FUNCTION DEFINITIONS
int fetch_max_steps(); //extracts the user defined maximum number of iterations allowed and returns it
string convert_int_to_string(int convert_in); //generic integer-to-string converter for use in file system code
int find_last_step(int max_steps); //finds the last complete step by finding the most recent step with the a "done.txt" file in it
void archive_incomplete_steps(int last_step, int max_steps); //archives any steps that may occur after the last complete step (i.e., crash and restart)
void create_next_step(int last_step); //post processes last step if no "post_process_done.txt" and then create the next step 
void copy_file(string initialFilePath, string outputFilePath); //generic copy file code for filesystem
bool last_step_post_process_status(int last_step); //code that actually checks for "post_process_done.txt" in last step
void post_process_last_step(int last_step); //actually post processes the data from the last step if not already done 
void update_auxillary_script(int last_step); //updates some commands to be executed by primary script after creating new step (i.e., grdf)
void new_step_grompp_files_script(int last_step); //creates script to grompp the files for simulation
void new_step_mdrun_gromacs_script(int last_step); //creates script to start the gromacs simulation
void new_step_rdf_gromacs_script(int last_step); //creates script to start the rdf calculation
void new_step_gromacs_scripts(int last_step); //calls the three gromacs creation scripts to be self contained and make modifying easier for other sim. packages
void copy_gromacs_files(string last_step_directory, string next_step_directory); //copies over gromacs required files to new step directory
array_pair_and_num_elements fetch_gr_data(string last_step_directory, gromacs_settings_class gromacs_settings); //extracts the RDF data and only reads in as much as is in the new gr file
void create_table_file(int last_step, string last_step_directory, array_pair_and_num_elements table_data, gromacs_settings_class gromacs_settings); //generates post processed table file
void create_new_grad_file(int last_step, double unscaled_gradient, double gr_convergence); //accumulates the unscaled gradient values from all iterations
void create_new_parameters_file(int last_step, int potential_type, double *potential_parameters, int num_parameters); //writes out the new set of updated potential parameters (zeroth step just reprints input)
void create_new_grompp_file(int last_step, double md_cutoff, gromacs_settings_class gromacs_settings); //creates new grompp file with adjusted md_cutoff
void create_new_grompp_files_script(int last_step); //creates new script to grompp files in next step
void create_new_mdrun_gromacs_script(int last_step); //creates new script to execute gromacs in next step
void create_new_rdf_gromacs_script(int last_step, gromacs_settings_class gromacs_settings); //creates new  script to calc the rdf with parallelization
void create_new_gromacs_scripts(int last_step, gromacs_settings_class gromacs_settings); //calls the above three to generate scripts
void create_post_process_done_file(int last_step); //creates post process done file
void create_new_momentum_file(int last_step, double *d_potential_parameters, int num_d_parameters); //creates the new momentum parameters file

//SPECIFIC POTENTIAL FUNCTION DEFINITIONS
array_pair_and_num_elements optimize_ideal_cluster_potential(int last_step, array_pair_and_num_elements gr_data,
	double *potential_parameters, double *d_potential_parameters, gromacs_settings_class gromacs_settings,
	double *md_cutoff_pointer, double *gr_cutoff_pointer, double *unscaled_gradient_pointer, double *gr_convergence_pointer);
array_pair_and_num_elements optimize_ramp_salr_cluster_potential(int last_step, array_pair_and_num_elements gr_data,
	double *potential_parameters, double *d_potential_parameters, gromacs_settings_class gromacs_settings,
	double *md_cutoff_pointer, double *gr_cutoff_pointer, double *unscaled_gradient_pointer, double *gr_convergence_pointer);


int main()
{
	//MAIN VARIABLES
	int last_step; //last step fully completed so start from that folder and contained data
	int max_steps; //maximum number of allowed iterations (if exceeded just clean up old stuff and re-label folders or change number to be larger)
	
	//FETCH THE MAXIMUM NUMBER OF STEPS
	max_steps=fetch_max_steps();

	//WRITE OUT DEMARCATION LINE FOR OUTPUT FROM CODE (THERE IS SCRIPT OUTPUT TOO)
	log_filestream << "///////////////////BEGIN RE CODE///////////////////" << endl;

	//ALL FUNCTIONS HAVE THEIR OWN ERROR CHECKING AND LOG FILE WRITING
	//SO NO MAIN BODY ERROR LOGGING OR KILLING OF PROGRAM IS FOUND HERE
	//IF AN ERROR IS TRIGGERED A KILL STATEMENT IS WRITTEN TO THE KILL SCRIPT CALLED BY THE MAIN SCRIPT

	//FIND THE LAST COMPLETE STEP
	last_step = find_last_step(max_steps);

	//ARCHIVE ANY FOLDERS FROM UNFINISHED STEPS
	archive_incomplete_steps(last_step, max_steps);

	//EXTRACT SETTINGS, CREATE THE NEXT STEP DIRECTORY, POST PROCESS LAST STEP DATA IF NEEDED AND SET UP SIMULATION FILES 
	create_next_step(last_step);

	//UPDATE AUXILLARY SCRIPT IN MAIN DIRECTORY TO SET GROMACS RDF CALCULATION TO WORKING STEP DIRECTORY
	update_auxillary_script(last_step);

	//CLOSE OUT LOG FILE - HOPEFULLY NOT TOO MANY ERRORS
	log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
	log_filestream.close();

	return 0;
}

//GENERAL FUNCTION DEFINITIONS
int fetch_max_steps()
{
	int max_steps;
	ifstream max_steps_filestream;

	max_steps_filestream.open("max_steps.txt");

	max_steps_filestream >> max_steps;

	//ERROR CHECK AND SEND ISSUES TO LOG FILE
	if (max_steps_filestream.fail() != true)
	{
		if (max_steps >= 0)
		{
			log_filestream << "retrieved maximum allowed steps: " << max_steps << endl;
		}
		else
		{
			log_filestream << "invalid maximum allowed steps: " << max_steps << " -> killing!" << endl;
			log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		log_filestream << "could not retrieve maximum allowed steps, check permissions -> killing!" << endl;
		log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
		exit(EXIT_FAILURE);
	}

	max_steps_filestream.close();

	return max_steps;
}

string convert_int_to_string(int convert_in)
{
	string convert_out; //string to hold the character version of the input integer
	stringstream convert_to_string; //temporary string converter

	convert_to_string.str(std::string()); 
	convert_to_string << convert_in; 
	convert_out = convert_to_string.str(); 
	
	return convert_out;
}

int find_last_step(int max_steps) //FIND LAST FULL STEP REGARDLESS IF MULTIPLE EARLIER FOLDERS WERE DELETED
{
	int last_step=-1; //last step fully completed so start from that folder and contained data
	int check_step; //step number to check
	ifstream check_step_filestream; //filestream for checking of the done file exists
	string check_step_done_file; //string holding address to done file for current check step iteration

	for (check_step = max_steps; check_step >= 0; check_step--)
	{
		check_step_done_file = "./step_" + convert_int_to_string(check_step) + "/done.txt";
		check_step_filestream.open(check_step_done_file.c_str());

		if (check_step_filestream.fail()!=true)
		{
			last_step = check_step; //if 0 this will cause a failure (as desired) and if finite it will continue on from that step
			check_step = -1; //break this loop
		}

		check_step_filestream.close();
	}

	//WRITE LAST COMPLETE STEP EXTRACTED TO LOG FILE
	log_filestream << "last complete step: " << last_step << endl;

	//CHECK LAST STEP AND IF INVALID WRITE ERROR TO LOG FILE AND KILL
	if (last_step >= 0 && last_step < max_steps)
	{
		log_filestream << "valid step -> proceeding" << endl;
	}
	else if (last_step == max_steps)
	{
		log_filestream << "reached maximum steps; please reset folder indices -> killing!" << endl;
		log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		log_filestream << "invalid step -> killing!" << endl;
		log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
		exit(EXIT_FAILURE);
	}

	return last_step;
}

void archive_incomplete_steps(int last_step, int max_steps)
{
	//FOR EXTRACTING THE SYSTEM TIME
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, 80, "%d-%m-%Y_%I:%M:%S", timeinfo);
	std::string str(buffer);

	std::cout << str;

	//ARCHIVE INCOMPLETE FOLDERS ACCORDING TO SYSTEM TIME AND PASS BACK FALSE IF FAILED TO KILL PROGRAM AND LOG ERROR
	bool archived_files = true; //test if archival was successful
	int rename_step; //start renaming folders from here
	string old_name, new_name; //old and new folder names
	struct stat info; //used to determine if a folder exists before attempting to rename it

	for (rename_step = last_step + 1; rename_step <= max_steps; rename_step++)
	{
		old_name = "./step_" + convert_int_to_string(rename_step);
		new_name = old_name + "_" + str;

		if (rename(old_name.c_str(), new_name.c_str()) != 0 && stat(old_name.c_str(), &info) == 0)
		{
			archived_files = false;
		}	
	}

	log_filestream << "archiving any unfinished steps" << endl;

	if (archived_files == true)
	{
		log_filestream << "unfinished steps successfully archived -> proceeding" << endl;
	}
	else
	{
		log_filestream << "could not archive unfinished steps, check file privileges -> killing!" << endl ;
		log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
		exit(EXIT_FAILURE);
	}
}

void create_next_step(int last_step)
{
	int next_step=last_step+1;
	string last_step_directory;
	string next_step_directory;

	//CREATE THE NEXT STEP FOLDER
	next_step_directory = "./step_" + convert_int_to_string(next_step);
	mkdir(next_step_directory.c_str(), 0700);

	//GENERATE LAST STEP DIRECTORY STRING FOR USAGE IN COPYING FILES
	last_step_directory = "./step_" + convert_int_to_string(last_step);

	//CHECK TO SEE IF THE LAST STEP WAS POST PROCESSED YET OR NOT (I.E. AFTER A KILL OR CRASH)
	if (last_step_post_process_status(last_step) == false)
	{
		post_process_last_step(last_step);
		//THIS USES THE LOCAL FOLDER SETTINGS FOR CONSISTENCY
		//NOTE: FOR THE ZEROTH STEP THIS DOES NO CALCULATION BUT IT DOES GENERATE
		//(1) PARAMETERS_OUT (JUST AN EXACT COPY)
		//(2) GROMMP_OUT (USES PARAMETERS TO INFORM CUTOFF DECISION)
		//(3) TABLE_OUT (FILL OUT TABLE WITH ZEROS TO HALF OF BOX LENGTH)
	}
	
	//COPY THE POST PROCESSED, UPDATED INFO FOR USE IN NEW STEP FOLDER SIMULATION
	copy_gromacs_files(last_step_directory, next_step_directory);

	//COPY MAIN DIRECTORY SETTINGS FILE TO DEMARCATE WHAT WAS USED
	//log_filestream << "copying current settings file to new step directory for reference" << endl;
	//copy_file("./settings.txt", next_step_directory + "/settings.txt");

	//CREATE NEW GROMACS SPECIFIC SCRIPTS FOR GROMPP, MDRUN AND RDF
	//new_step_gromacs_scripts(last_step);

}

void copy_file(string initialFilePath, string outputFilePath)  //COPY ANY FILES THAT I NEED TO MOVE FOR VARIOUS ITERATIONS
{
	ifstream initialFile(initialFilePath.c_str(), ios::in | ios::binary); //input stream for reading
	ofstream outputFile(outputFilePath.c_str(), ios::out | ios::binary); //output stream for writing
	string line; //string to pipe file data into

	while (getline(initialFile, line))
	{
		if (!(outputFile << line << endl))
		{
			log_filestream << "could not copy from " << initialFilePath << " to " << outputFilePath << " -> killing!" << endl;
			log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
			exit(EXIT_FAILURE);
		}
	}

	initialFile.close();
	outputFile.close();
}

bool last_step_post_process_status(int last_step)
{
	bool post_process_status; //returned true or false for a fully post processes data set
	ifstream check_post_process_filestream; //filestream for checking of the post processing done file exists
	string check_post_process_done_file; //string holding address to post processing done file for current check step iteration
	
	check_post_process_done_file = "./step_" + convert_int_to_string(last_step) + "/post_process_done.txt";
	check_post_process_filestream.open(check_post_process_done_file.c_str());

	log_filestream << "checking post process status of the last step" << endl;

	if (check_post_process_filestream.fail() == true)
	{
		post_process_status = false;
		log_filestream << "need to post process data from last step" << endl;
	}
	else
	{
		post_process_status = true;
		log_filestream << "data from last step already post processed" << endl;
	}

	check_post_process_filestream.close();

	return post_process_status;
}

void post_process_last_step(int last_step)
{
	//LAST STEP DIRECTORY NEEDED FOR EVERYTHING
	string last_step_directory = "./step_" + convert_int_to_string(last_step);
	string setting_file_address = last_step_directory + "/settings.txt";
	string potential_parameters_file_address = last_step_directory + "/parameters.txt";
	string d_potential_parameters_file_address = last_step_directory + "/d_parameters.txt";

	//GROMACS SPECIFIC VARIABLES
	gromacs_settings_class gromacs_settings; //will store the needed gromacs settings
	ifstream gromacs_settings_filestream; //filestream to read in the gromacs specific settings

	//POTENTIAL SPECIFIC VARIABLES
	int potential_type; //specifies the type of potential (first line of parameters file) and thus the number of parameters, N
	int num_parameters; //number of parameters used for this potential
	int num_d_parameters; //number of momentum variables
	double *potential_parameters; //stores the potential parameters [0]-[N-1] to pass to appropriate function
	double *d_potential_parameters; //stores the update applied to be used in the momentum addition step
	ifstream potential_parameters_filestream; //filestream to read in the potential type and various parameters
	ifstream d_potential_parameters_filestream; //filestream to read in the potential type and various parameters
	int i; //loop variable to read in parameters

	//TABLE ARRAY AND NUMBER OF ENTRIES (WILL NO BE NORMALIZED BY KBT)
	array_pair_and_num_elements table_data; //brings back address to two table arrays and the number of elements
	
	//GR ARRAY AND NUMBER OF ENTRIES
	array_pair_and_num_elements gr_data; //brings back the address to two gr arrays and the number of elements

	//CUTOFFS TO BE DETERMINED BY OPTIMIZATION STEP
	double md_cutoff, gr_cutoff; //md and rdf calculation cutoffs
	double unscaled_gradient; //asesses the convergence of the optimization scheme
	double gr_convergence; //asesses the convergence of the optimization scheme but with gr

	//OPEN OUR THREE FILESTREAMS
	gromacs_settings_filestream.open(setting_file_address.c_str());
	potential_parameters_filestream.open(potential_parameters_file_address.c_str());
	d_potential_parameters_filestream.open(d_potential_parameters_file_address.c_str());

	//READ IN SETTINGS
	if (!(gromacs_settings_filestream >> gromacs_settings.delta_r) ||
		!(gromacs_settings_filestream >> gromacs_settings.num_gr_calcs) ||
		!(gromacs_settings_filestream >> gromacs_settings.equil_time) ||
		!(gromacs_settings_filestream >> gromacs_settings.final_time) ||
		!(gromacs_settings_filestream >> gromacs_settings.end_table) ||
		!(gromacs_settings_filestream >> gromacs_settings.gradient_scale) ||
		!(gromacs_settings_filestream >> gromacs_settings.momentum_scale) ||
		!(gromacs_settings_filestream >> gromacs_settings.kB_T) ||
		!(gromacs_settings_filestream >> gromacs_settings.md_cutoff_magnitude) ||
		!(gromacs_settings_filestream >> gromacs_settings.gr_cutoff_magnitude) ||
		!(gromacs_settings_filestream >> gromacs_settings.buffer_size) ||
		!(gromacs_settings_filestream >> gromacs_settings.rlist) ||
		!(gromacs_settings_filestream >> gromacs_settings.rcoulomb) ||
		!(gromacs_settings_filestream >> gromacs_settings.rvdw))
	{
		log_filestream << "local last step settings could not be extracted -> killing!" << endl;
		log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
		exit(EXIT_FAILURE);
	}

	//READ IN THE NEW AND TARGET RDFS IF NOT THE LAST STEP
	if (last_step != 0)
	{
		gr_data = fetch_gr_data(last_step_directory, gromacs_settings);
	}
	
	
	//READ IN THE FIRST LINE OF PARAMETERS FILE TO ESTABLISH THE POTENTIAL TYPE (AND THUS # OF PARAMETERS)
	if (!(potential_parameters_filestream >> potential_type))
	{
		log_filestream << "could not extract potential type from last step parameters file -> killing!" << endl;
		log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
		exit(EXIT_FAILURE);
	}

	//ALLOCATE APPROPRIATE ARRAY SIZE, READ IN POTENTIAL PARAMETERS AND FINALLY OPTIMIZE USING POTENTIAL SPECIFIC FUNCTIONS
	//THESE FUNCTIONS ALSO GENERATE THE GROMACS TABLE ARRAYS
	if (potential_type == 0) //ideal cluster potential with two energy and range specific parameters
	{
		num_parameters = 7;
		num_d_parameters = 4;
		potential_parameters = (double*)malloc(sizeof(double) * num_parameters);
		d_potential_parameters = (double*)malloc(sizeof(double) * num_d_parameters);

		for (i = 0; i < num_parameters; i++)
		{
			potential_parameters_filestream >> *(potential_parameters + i);
		}
		for (i = 0; i < num_d_parameters; i++)
		{
			d_potential_parameters_filestream >> *(d_potential_parameters + i);
		}

		table_data = optimize_ideal_cluster_potential(last_step, gr_data,
			potential_parameters, d_potential_parameters, gromacs_settings, &md_cutoff, &gr_cutoff, &unscaled_gradient, &gr_convergence);
	}
	else if (potential_type == 1)
	{
		num_parameters = 6;
		num_d_parameters = 4;
		potential_parameters = (double*)malloc(sizeof(double) * num_parameters);
		d_potential_parameters = (double*)malloc(sizeof(double) * num_d_parameters);

		for (i = 0; i < num_parameters; i++)
		{
			potential_parameters_filestream >> *(potential_parameters + i);
		}
		for (i = 0; i < num_d_parameters; i++)
		{
			d_potential_parameters_filestream >> *(d_potential_parameters + i);
		}

		table_data = optimize_ramp_salr_cluster_potential(last_step, gr_data,
			potential_parameters, d_potential_parameters, gromacs_settings, &md_cutoff, &gr_cutoff, &unscaled_gradient, &gr_convergence);
	}
	else
	{
		log_filestream << "invalid potential type -> killing!" << endl;
		log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
		exit(EXIT_FAILURE);
	} 

	//APPLY THE BUFFER
	md_cutoff = md_cutoff + gromacs_settings.buffer_size;

	//CREATE GRADIENT CONVERGENCE FILE
	create_new_grad_file(last_step, unscaled_gradient, gr_convergence);

	//CREATE MOMENTUM FILE FOR UPDATE
	create_new_momentum_file(last_step, d_potential_parameters, num_d_parameters);

	//CREATE NEW TABLE FILE
	create_table_file(last_step, last_step_directory, table_data, gromacs_settings);

	//CREATE NEW GROMACS FILES USING CUTOFFS
	create_new_grompp_file(last_step, md_cutoff, gromacs_settings);

	//CREATE NEW PARAMETERS FILE
	create_new_parameters_file(last_step, potential_type, potential_parameters, num_parameters);

	//CREATE NEW SCRIPT FILES
	create_new_gromacs_scripts(last_step, gromacs_settings);

	//CLOSE FILESTREAMS
	gromacs_settings_filestream.close();
	potential_parameters_filestream.close();

	//CREATE NEW DONE_POST_PROCESS.TXT
	create_post_process_done_file(last_step);

	//DEALLOCATE ALL ARRAYS

}

void update_auxillary_script(int last_step)
{
	string new_step_directory = "./step_" + convert_int_to_string(last_step + 1); //contains string with new step directory
	ofstream auxillary_script_filestream("./auxillary.sh", ios::out | ios::trunc); //for writing commands to be run by the global script (for gromacs just g rdf)

	//CHANGE TO APPROPRIATE DIRECTORY
	auxillary_script_filestream << "cd " << new_step_directory << endl;
	auxillary_script_filestream << endl;

	//BEGINNING HEADER FOR AUXILLARY SCRIPT TO WRITE TO LOG FILE
	auxillary_script_filestream << "echo \"//////////////////BEGIN AUX SCRIPT/////////////////\" >> ../log.txt" << endl;
	auxillary_script_filestream << "echo \"preparing simulation for step " + convert_int_to_string(last_step + 1) << "\" >> ../log.txt" << endl;
	auxillary_script_filestream << endl;

	//JUST WRITING HEADER FOR CLARITY IN SCRIPT FILE
	auxillary_script_filestream << "##############################################" << endl;
	auxillary_script_filestream << "################GROMMP SCRIPT#################" << endl;
	auxillary_script_filestream << "##############################################" << endl;
	auxillary_script_filestream << endl;

	//CREATE COMMAND TO EXECUTE GROMPP
	auxillary_script_filestream << "echo \"executing grompp to prepare for gromacs run\" >> ../log.txt" << endl;
	auxillary_script_filestream << "chmod 744 ./grompp_files.sh" << endl;
	auxillary_script_filestream << "./grompp_files.sh" << endl;
	auxillary_script_filestream << endl;
	//ERROR CHECK GROMPP SCRIPT AND SEND BACK FAILURE IF NEED BE
	auxillary_script_filestream << "grompp_exit=$?" << endl;
	auxillary_script_filestream << endl;
	auxillary_script_filestream << "if [[ $grompp_exit != 0 ]]" << endl;
	auxillary_script_filestream << "then" << endl;
	auxillary_script_filestream << "   echo \"grompp was unsuccessful -> killing!\" >> ../log.txt" << endl;
	auxillary_script_filestream << "   echo \"///////////////////END AUX SCRIPT//////////////////\" >> ../log.txt" << endl;
	auxillary_script_filestream << "   echo \" \" >> ../log.txt" << endl;
	auxillary_script_filestream << "   exit $grompp_exit" << endl;
	auxillary_script_filestream << "fi" << endl;
	auxillary_script_filestream << endl;




	//AGAIN, WRITING HEADER FOR CLARITY IN SCRIPT FILE
	auxillary_script_filestream << "##############################################" << endl;
	auxillary_script_filestream << "################MDRUN SCRIPT##################" << endl;
	auxillary_script_filestream << "##############################################" << endl;
	auxillary_script_filestream << endl;

	//CREATE COMMAND TO EXECUTE MDRUN
	auxillary_script_filestream << "echo \"executing gromacs run\" >> ../log.txt" << endl;
	auxillary_script_filestream << "chmod 744 ./mdrun_gromacs.sh" << endl;
	auxillary_script_filestream << "./mdrun_gromacs.sh" << endl;
	auxillary_script_filestream << endl;
	//ERROR CHECK GROMPP SCRIPT AND SEND BACK FAILURE IF NEED BE
	auxillary_script_filestream << "mdrun_exit=$?" << endl;
	auxillary_script_filestream << endl;
	auxillary_script_filestream << "if [[ $mdrun_exit != 0 ]]" << endl;
	auxillary_script_filestream << "then" << endl;
	auxillary_script_filestream << "   echo \"mdrun was unsuccessful -> killing!\" >> ../log.txt" << endl;
	auxillary_script_filestream << "   echo \"///////////////////END AUX SCRIPT//////////////////\" >> ../log.txt" << endl;
	auxillary_script_filestream << "   echo \" \" >> ../log.txt" << endl;
	auxillary_script_filestream << "   exit $mdrun_exit" << endl;
	auxillary_script_filestream << "fi" << endl;
	auxillary_script_filestream << endl;




	//AND YET AGAIN, WRITING HEADER FOR CLARITY IN SCRIPT FILE
	auxillary_script_filestream << "##############################################" << endl;
	auxillary_script_filestream << "#################RDF SCRIPT###################" << endl;
	auxillary_script_filestream << "##############################################" << endl;
	auxillary_script_filestream << endl;

	//CREATE COMMAND TO EXECUTE MDRUN
	auxillary_script_filestream << "echo \"calculating rdf with gromacs\" >> ../log.txt" << endl;
	auxillary_script_filestream << "chmod 744 ./rdf_gromacs.sh" << endl;
	auxillary_script_filestream << "./rdf_gromacs.sh" << endl;
	auxillary_script_filestream << endl;
	//ERROR CHECK GROMPP SCRIPT AND SEND BACK FAILURE IF NEED BE
	auxillary_script_filestream << "rdf_exit=$?" << endl;
	auxillary_script_filestream << endl;
	auxillary_script_filestream << "if [[ $rdf_exit != 0 ]]" << endl;
	auxillary_script_filestream << "then" << endl;
	auxillary_script_filestream << "   echo \"rdf calculation was unsuccessful -> killing!\" >> ../log.txt" << endl;
	auxillary_script_filestream << "   echo \"///////////////////END AUX SCRIPT//////////////////\" >> ../log.txt" << endl;
	auxillary_script_filestream << "   echo \" \" >> ../log.txt" << endl;
	auxillary_script_filestream << "   exit $rdf_exit" << endl;
	auxillary_script_filestream << "fi" << endl;
	auxillary_script_filestream << endl;




	//WRITE THE DONE FILE SEGMENT OF SCRIPT
	auxillary_script_filestream << "##############################################" << endl;
	auxillary_script_filestream << "#################WRITE DONE###################" << endl;
	auxillary_script_filestream << "##############################################" << endl;
	auxillary_script_filestream << endl;
	auxillary_script_filestream << "echo \"creating done file\" >> ../log.txt" << endl;
	auxillary_script_filestream << "echo \" \" >> ./done.txt" << endl;
	auxillary_script_filestream << endl;
	//CHECK TO ENSURE IT WAS WRITTEN
	auxillary_script_filestream << "wd_exit=$?" << endl;
	auxillary_script_filestream << endl;
	auxillary_script_filestream << "if [[ $wd_exit != 0 ]]" << endl;
	auxillary_script_filestream << "then" << endl;
	auxillary_script_filestream << "   echo \"done file write unsuccessful -> killing!\" >> ../log.txt" << endl;
	auxillary_script_filestream << "   echo \"///////////////////END AUX SCRIPT//////////////////\" >> ../log.txt" << endl;
	auxillary_script_filestream << "   echo \" \" >> ../log.txt" << endl;
	auxillary_script_filestream << "   exit $wd_exit" << endl;
	auxillary_script_filestream << "fi" << endl;
	auxillary_script_filestream << endl;


	//ENDING HEADER FOR AUXILLARY SCRIPT TO WRITE TO LOG FILE
	auxillary_script_filestream << "echo \"///////////////////END AUX SCRIPT//////////////////\" >> ../log.txt" << endl;
	auxillary_script_filestream << "echo \" \" >> ../log.txt" << endl;
	auxillary_script_filestream << endl;


	//CHANGE BACK TO MAIN DIRECTORY
	auxillary_script_filestream << "cd .." << endl;
	auxillary_script_filestream << endl;


	//LOG ANY ERROR AND KILL IF NECCESSARY
	if (!auxillary_script_filestream.is_open())
	{
		log_filestream << "couldn't open or write to ./auxillary.sh -> killing!" << endl;
		log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
		exit(EXIT_FAILURE);
	}
	
	auxillary_script_filestream.close();
}

void new_step_grompp_files_script(int last_step)
{
	string file_address = "./step_" + convert_int_to_string(last_step + 1) + "/grompp_files.sh"; //contains string with new file directory
	ofstream grompp_files_script_filestream(file_address.c_str(), ios::out | ios::trunc); //for writing commands to be run by the grompp script

	grompp_files_script_filestream << "grompp -f grompp.mdp -po md_out.mdp -c conf.gro -n index.ndx -p topol.top -o topol.tpr" << endl;
	grompp_files_script_filestream << endl;
	grompp_files_script_filestream << "grompp_exit=$?" << endl;
	grompp_files_script_filestream << endl;
	grompp_files_script_filestream << "exit $grompp_exit" << endl;

	grompp_files_script_filestream.close();
}

void new_step_mdrun_gromacs_script(int last_step)
{
	string file_address = "./step_" + convert_int_to_string(last_step + 1) + "/mdrun_gromacs.sh"; //contains string with new file directory
	ofstream mdrun_gromacs_script_filestream(file_address.c_str(), ios::out | ios::trunc); //for writing commands to be run by the mdrun script

	mdrun_gromacs_script_filestream << "mdrun -s topol.tpr -c conf_out.gro -o traj.trr -x traj.xtc" << endl;
	mdrun_gromacs_script_filestream << endl;
	mdrun_gromacs_script_filestream << "mdrun_exit=$?" << endl;
	mdrun_gromacs_script_filestream << endl;
	mdrun_gromacs_script_filestream << "exit $mdrun_exit" << endl;

	mdrun_gromacs_script_filestream.close();
}

void new_step_rdf_gromacs_script(int last_step)
{
	string file_address = "./step_" + convert_int_to_string(last_step + 1) + "/rdf_gromacs.sh"; //contains string with new file directory
	ofstream rdf_gromacs_script_filestream(file_address.c_str(), ios::out | ios::trunc); //for writing commands to be run by the rdf script

	rdf_gromacs_script_filestream << "g_rdf -f traj.xtc -s topol.tpr -n index.ndx -o rdf.xvg << 'EOF'" << endl;
	rdf_gromacs_script_filestream << "2" << endl;
	rdf_gromacs_script_filestream << "2" << endl;
	rdf_gromacs_script_filestream << "'EOF'" << endl;
	rdf_gromacs_script_filestream << endl;
	rdf_gromacs_script_filestream << "rdf_exit=$?" << endl;
	rdf_gromacs_script_filestream << endl;
	rdf_gromacs_script_filestream << "exit $rdf_exit" << endl;

	rdf_gromacs_script_filestream.close();
}

void new_step_gromacs_scripts(int last_step)
{
	new_step_grompp_files_script(last_step);
	new_step_mdrun_gromacs_script(last_step);
	new_step_rdf_gromacs_script(last_step);
}

void copy_gromacs_files(string last_step_directory, string next_step_directory)
{
	/*//COPY THE POST PROCESSED, UPDATED INFO FOR USE IN NEW STEP FOLDER
	//new starting state configuration, parameters (plus table file) and grompp file with cutoff adjustment
	log_filestream << "copying last step simulation files to new step directory" << endl;
	copy_file(last_step_directory + "/conf_out.gro", next_step_directory + "/conf.gro"); //CHANGE TO CONF_OUT
	copy_file(last_step_directory + "/parameters.txt", next_step_directory + "/parameters.txt"); //CHANGE TO PARAMETERS_OUT
	copy_file(last_step_directory + "/table_out.xvg", next_step_directory + "/table.xvg"); //CHANGE TO TABLE_OUT
	copy_file(last_step_directory + "/grompp.mdp", next_step_directory + "/grompp.mdp"); //CHANGE TO GROMPP_OUT
	//generic files for gromacs that will not be modified
	copy_file(last_step_directory + "/index.ndx", next_step_directory + "/index.ndx");
	copy_file(last_step_directory + "/topol.top", next_step_directory + "/topol.top");
	//target RDF that can be changed mid run by adjusting in last complete step
	copy_file(last_step_directory + "/rdf_tgt.xvg", next_step_directory + "/rdf_tgt.xvg");*/


	//COPY THE POST PROCESSED, UPDATED INFO FOR USE IN NEW STEP FOLDER
	//new starting state configuration, parameters (plus table file) and grompp file with cutoff adjustment
	log_filestream << "copying last step simulation files to new step directory" << endl;
	copy_file(last_step_directory + "/conf_out.gro", next_step_directory + "/conf.gro");
	copy_file(last_step_directory + "/parameters_out.txt", next_step_directory + "/parameters.txt");
	copy_file(last_step_directory + "/d_parameters_out.txt", next_step_directory + "/d_parameters.txt");
	copy_file(last_step_directory + "/table_out.xvg", next_step_directory + "/table.xvg");
	copy_file(last_step_directory + "/grompp_out.mdp", next_step_directory + "/grompp.mdp"); //CHANGE TO GROMPP_OUT
	copy_file(last_step_directory + "/rdf_gromacs_out.sh", next_step_directory + "/rdf_gromacs.sh"); 
	copy_file(last_step_directory + "/mdrun_gromacs_out.sh", next_step_directory + "/mdrun_gromacs.sh");
	copy_file(last_step_directory + "/grompp_files_out.sh", next_step_directory + "/grompp_files.sh");
	//generic files for gromacs that will not be modified
	copy_file(last_step_directory + "/index.ndx", next_step_directory + "/index.ndx");
	copy_file(last_step_directory + "/topol.top", next_step_directory + "/topol.top");
	//target RDF that can be changed mid run by adjusting in last complete step
	copy_file(last_step_directory + "/rdf_tgt.xvg", next_step_directory + "/rdf_tgt.xvg");
	copy_file(last_step_directory + "/settings.txt", next_step_directory + "/settings.txt");
}

array_pair_and_num_elements fetch_gr_data(string last_step_directory, gromacs_settings_class gromacs_settings)
{
	string rdf_file_address = last_step_directory + "/rdf.xvg"; //address for new simulated rdf file
	string rdf_tgt_file_address = last_step_directory + "/rdf_tgt.xvg"; //address for target rdf file
	int num_lines_gr, num_lines_gr_tgt; //stores the number of lines--needed for allocating the RDF arrays (target gr is always longer than needed)
	string line; //dummy variable to store string in while counting number of lines
	ifstream gr_filestream, gr_tgt_filestream; //filestreams for reading in the RDF data
	int i; //generic loop variable
	double r_temp; //temporary storage of the r values for reading in RDF files
	double *r, *r_tgt; //r values read in for RDFS to check if read properly
	double check_alignment; //for checking the difference between the stream alignment
	double check_spacing_1, check_spacing_2; //for checking that the file spacing matches the settings file
	double *gr, *gr_tgt; //for allocating array and referencing before loading into gr_data
	array_pair_and_num_elements gr_data; //passes back array address and number of array elements
	double garbage; //stores garbage extra input from VOTCA rdf calculation

	//LOG CURRENT PROCESS
	log_filestream << "extracting RDF data for calculating update" << endl;

	//OPEN THE TWO RDF FILESTREAMS
	gr_filestream.open(rdf_file_address.c_str());
	gr_tgt_filestream.open(rdf_tgt_file_address.c_str());

	//DETERMINE THE TOTAL NUMBER OF POINTS TO READ IN FOR ALLOCATING ARRAYS
	num_lines_gr = 0; /*-rdf_gromacs_skip*/; //gromacs rdf defaults to having 13 lines of header text (this is set here by a constant global integer)
	while (getline(gr_filestream, line))
	{
		++num_lines_gr;
	}
	num_lines_gr_tgt = 0;
	while (getline(gr_tgt_filestream, line))
	{
		++num_lines_gr_tgt;
	}

	//ERROR CHECK TO MAKE SURE TARGET RDF HAS MORE DATA THAN THE NEW SIMULATED RDF
	if (num_lines_gr_tgt < num_lines_gr)
	{
		num_lines_gr = num_lines_gr_tgt;
		//log_filestream << num_lines_gr << "," << num_lines_gr_tgt << endl;
		//log_filestream << "new RDF has more points than target -> killing!" << endl;
		//log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
		//exit(EXIT_FAILURE);
	}

	//CLEAR AND RESET FILESTREAMS
	gr_filestream.clear(); 
	gr_filestream.seekg(0, ios::beg);
	gr_tgt_filestream.clear();
	gr_tgt_filestream.seekg(0, ios::beg);

	//ALLOCATE THE ARRAYS BASED ON AMOUNT OF DATA IN NEW RDF FILE AND SET THE ARRAY LENGTH TO BE PASSED BACK VIA THE POINTER
	gr = (double*)malloc(sizeof(double) * num_lines_gr);
	gr_tgt = (double*)malloc(sizeof(double) * num_lines_gr);
	r = (double*)malloc(sizeof(double) * num_lines_gr);
	r_tgt = (double*)malloc(sizeof(double) * num_lines_gr);

	//READ IN RDF DATA AND STORE IN ARRAYS
	for (i = 0; i < rdf_gromacs_skip; i++) //skip dummy lines in gromacs generated file
	{
		getline(gr_filestream, line);
	}
	for (i = 0; i < num_lines_gr; i++) //now fill in with actual gr data
	{
		gr_filestream >> scientific >> *(r+i) >> *(gr + i) >> garbage;
		gr_tgt_filestream >> scientific >> *(r_tgt + i) >> *(gr_tgt + i);

		//CHECK ALIGNMENT
		check_alignment = fabs(*(r + i) - *(r_tgt + i)) / gromacs_settings.delta_r;

		if (check_alignment > 0.000000001)
		{
			log_filestream << "target and current rdf inputs seem unaligned -> killing" << endl;
			log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
			exit(EXIT_FAILURE);
		}

		//CHECK IF SPACING MATCHES THAT IN SETTINGS FILE
		if (i > 0)
		{
			check_spacing_1 = fabs(1.0 - (*(r + i) - *(r + i - 1)) / gromacs_settings.delta_r);
			check_spacing_2 = fabs(1.0 - (*(r_tgt + i) - *(r_tgt + i - 1)) / gromacs_settings.delta_r);
			
			if (check_spacing_1 > 0.000000001 || check_spacing_2 > 0.000000001)
			{
				log_filestream << "spacing in target and/or new rdf seems wrong -> killing" << endl;
				log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
				exit(EXIT_FAILURE);
			}
		}
	}

	//ERROR CHECK FILESTREAM AND THAT THE GR DATA BEGINS AT 0
	if (gr_filestream.fail() == true || gr_filestream.fail() == true)
	{
		log_filestream << "failed reading in RDF data -> killing!" << endl;
		log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
		exit(EXIT_FAILURE);
	}
	else if (*(r + 0) > 0.000000001 || *(r_tgt + 0) > 0.000000001)
	{
		log_filestream << "one or both RDF files do not seem to begin at r=0 -> killing!" << endl;
		log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
		exit(EXIT_FAILURE);
	}
	
	
	//CLOSE FILESTREAMS
	gr_filestream.close();
	gr_tgt_filestream.close();

	//NOTE: THE ARRAYS ARE DEALLOCATED BACK IN THE POST_PROCESS_LAST_STEP FUNCTION


	//FOR TESTING PURPOSES ONLY (WRITE OUT NEW GR FILES TO EXAMINE)
	ofstream both_gr_filestream_out;
	string both_gr_file_address = last_step_directory + "/rdf_combined.xvg";
	both_gr_filestream_out.open(both_gr_file_address.c_str());

	for (i = 0; i < num_lines_gr; i++)
	{
		both_gr_filestream_out << *(r + i) << "," << *(gr + i) << "," << *(r_tgt+i) << "," << *(gr_tgt + i) << endl;
	}

	both_gr_filestream_out.close();

	//LOAD IN THE GR DATA
	gr_data.array_1 = gr;
	gr_data.array_2 = gr_tgt;
	gr_data.num_elements = num_lines_gr;

	return gr_data;
}

void create_table_file(int last_step, string last_step_directory, array_pair_and_num_elements table_data, gromacs_settings_class gromacs_settings)
{
	int i; //generic loop variable
	double r; //generic distance variable
	ofstream table_file_stream; //table file stream
	string table_file_address = last_step_directory + "/table_out.xvg"; //address for new simulated rdf file

	//OPEN FILE STREAM AND WRITE HEADER
	table_file_stream.open(table_file_address.c_str());
	table_file_stream << "#post processed table file for step " << last_step << endl;

	for (i = 0; i < table_data.num_elements; i++)
	{
		r = (double)i*gromacs_settings.delta_r;

		table_file_stream.precision(10);

		table_file_stream << scientific << r << "   " << 0.0 << " " << 0.0 << "   "
			<< 0.0 << " " << 0.0 << "   " << *(table_data.array_1 + i) << " " << *(table_data.array_2 + i) << endl;
	}

	table_file_stream.close();
}

void create_new_grad_file(int last_step, double unscaled_gradient, double gr_convergence)
{
	ofstream unscaled_gradient_filestream;

	//IF LAST STEP IS ZERO THEN CREATE BLANK GRAD FILE IN FOLDER 0 OTHERWISE COPY ONE FROM LAST_STEP-1 AND APPEND
	if (last_step == 0)
	{
		string last_step_file_path = "./step_0/unscaled_gradient.txt";
		unscaled_gradient_filestream.open(last_step_file_path.c_str(), ios::out | ios::trunc);
		unscaled_gradient_filestream.precision(8);
		if (!(unscaled_gradient_filestream << "last step" << "," << "unscaled gradient" << "," << "gr convergence" << endl))
		{
			log_filestream << "could not initialize first gradient file -> killing" << endl;
			log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
			exit(EXIT_FAILURE);
		}
		unscaled_gradient_filestream.close();
	}
	else
	{
		string last_last_step_file_path = "./step_" + convert_int_to_string(last_step-1) + "/unscaled_gradient.txt";
		string last_step_file_path = "./step_" + convert_int_to_string(last_step) + "/unscaled_gradient.txt";;
		copy_file(last_last_step_file_path, last_step_file_path);
		unscaled_gradient_filestream.open(last_step_file_path.c_str(), ios::out | ios::app);
		unscaled_gradient_filestream.precision(8);
		if (!(unscaled_gradient_filestream << last_step << "," << unscaled_gradient << "," << gr_convergence << endl))
		{
			log_filestream << "could not write data to gradient file -> killing" << endl;
			log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
			exit(EXIT_FAILURE);
		}
		unscaled_gradient_filestream.close();
	}
}

void create_new_parameters_file(int last_step, int potential_type, double *potential_parameters, int num_parameters)
{
	int i; //generic loop variable
	string parameters_file_directory = "./step_" + convert_int_to_string(last_step) + "/parameters_out.txt"; //parameters file address
	ofstream parameters_filestream(parameters_file_directory.c_str()); //filestream for writing parameters file

	parameters_filestream.precision(8);

	//WRITE OUT THE POTENTIAL TYPE
	if (!(parameters_filestream << potential_type << endl))
	{
		log_filestream << "could not write the new parameters file potential type-> killing" << endl;
		log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
		exit(EXIT_FAILURE);
	}

	//WRITE OUT THE UPDATED PARAMETERS (NOT REALLY UPDATED IF LAST_STEP=0)
	for (i = 0; i < num_parameters; i++)
	{
		if (!(parameters_filestream << *(potential_parameters + i) << endl))
		{
			log_filestream << "could not write the new parameters file updated parameters -> killing" << endl;
			log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
			exit(EXIT_FAILURE);
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void create_new_grompp_file(int last_step, double md_cutoff, gromacs_settings_class gromacs_settings)
{
	string grompp_address = "./step_" + convert_int_to_string(last_step) + "/grompp.mdp"; //contains string with new file directory
	string grompp_out_address = "./step_" + convert_int_to_string(last_step) + "/grompp_out.mdp"; //contains string with new file directory
	ifstream grompp_filestream(grompp_address.c_str()); //for writing commands to be run by the grompp script
	ofstream grompp_out_filestream(grompp_out_address.c_str(), ios::out | ios::trunc); //for writing commands to be run by the grompp script
	string line; //generic string to store line of file in to search for keyword
	int cur_line = 0; //current line in the file

	grompp_out_filestream.precision(4);

	while (getline(grompp_filestream, line)) 
	{
		if (cur_line == gromacs_settings.rlist-1)
		{
			grompp_out_filestream << "rlist                    = " << md_cutoff << endl;
		}
		else if (cur_line == gromacs_settings.rcoulomb-1)
		{
			grompp_out_filestream << "rcoulomb                 = " << md_cutoff << endl;
		}
		else if (cur_line == gromacs_settings.rvdw-1)
		{
			grompp_out_filestream << "rvdw                     = " << md_cutoff << endl;
		}
		else
		{
			grompp_out_filestream << line << endl;
		}
		cur_line++;
	}

	grompp_filestream.close();
	grompp_out_filestream.close();
}

void create_new_grompp_files_script(int last_step)
{
	string file_address = "./step_" + convert_int_to_string(last_step) + "/grompp_files_out.sh"; //contains string with new file directory
	ofstream grompp_files_script_filestream(file_address.c_str(), ios::out | ios::trunc); //for writing commands to be run by the grompp script

	grompp_files_script_filestream << "grompp -f grompp.mdp -po md_out.mdp -c conf.gro -n index.ndx -p topol.top -o topol.tpr" << endl;
	grompp_files_script_filestream << endl;
	grompp_files_script_filestream << "grompp_exit=$?" << endl;
	grompp_files_script_filestream << endl;
	grompp_files_script_filestream << "exit $grompp_exit" << endl;

	grompp_files_script_filestream.close();
}

void create_new_mdrun_gromacs_script(int last_step)
{
	string file_address = "./step_" + convert_int_to_string(last_step) + "/mdrun_gromacs_out.sh"; //contains string with new file directory
	ofstream mdrun_gromacs_script_filestream(file_address.c_str(), ios::out | ios::trunc); //for writing commands to be run by the mdrun script

	mdrun_gromacs_script_filestream << "mdrun -s topol.tpr -c conf_out.gro -o traj.trr -x traj.xtc" << endl;
	mdrun_gromacs_script_filestream << endl;
	mdrun_gromacs_script_filestream << "mdrun_exit=$?" << endl;
	mdrun_gromacs_script_filestream << endl;
	mdrun_gromacs_script_filestream << "exit $mdrun_exit" << endl;

	mdrun_gromacs_script_filestream.close();
}

void create_new_rdf_gromacs_script(int last_step, gromacs_settings_class gromacs_settings)
{
	string file_address = "./step_" + convert_int_to_string(last_step) + "/rdf_gromacs_out.sh"; //contains string with new file directory
	ofstream rdf_gromacs_script_filestream(file_address.c_str(), ios::out | ios::trunc); //for writing commands to be run by the rdf script

	rdf_gromacs_script_filestream << "../multi_g_rdf" <<
		/*" -N " << gromacs_settings.num_gr_calcs <<*/
		" -b " << gromacs_settings.equil_time <<
		" -e " << gromacs_settings.final_time <<
		" --" <<
		" -bin " << gromacs_settings.delta_r <<
		" << 'EOF' " << endl;
	rdf_gromacs_script_filestream << "2" << endl;
	rdf_gromacs_script_filestream << "2" << endl;
	rdf_gromacs_script_filestream << "'EOF'" << endl;
	rdf_gromacs_script_filestream << endl;
	rdf_gromacs_script_filestream << "rdf_exit=$?" << endl;
	rdf_gromacs_script_filestream << endl;
	rdf_gromacs_script_filestream << "exit $rdf_exit" << endl;

	rdf_gromacs_script_filestream.close();
}

void create_new_gromacs_scripts(int last_step, gromacs_settings_class gromacs_settings)
{
	create_new_grompp_files_script(last_step);
	create_new_mdrun_gromacs_script(last_step);
	create_new_rdf_gromacs_script(last_step, gromacs_settings);
}

void create_post_process_done_file(int last_step)
{
	string file_address = "./step_" + convert_int_to_string(last_step) + "/post_process_done.txt"; //contains string with new file directory
	ofstream post_process_done_filestream(file_address.c_str(), ios::out | ios::trunc); //for writing commands to be run by the rdf script

	log_filestream << "writing post process done indicator file" << endl;

	if (!(post_process_done_filestream << " " << endl))
	{
		log_filestream << "could not write the post process done indicator file -> killing" << endl;
		log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
		exit(EXIT_FAILURE);
	}

	post_process_done_filestream.close();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void create_new_momentum_file(int last_step, double *d_potential_parameters, int num_d_parameters)
{
	int i; //generic loop variable
	string d_parameters_file_directory = "./step_" + convert_int_to_string(last_step) + "/d_parameters_out.txt"; //parameters file address
	ofstream d_parameters_filestream(d_parameters_file_directory.c_str()); //filestream for writing parameters file

	d_parameters_filestream.precision(8);

	//WRITE OUT THE UPDATED MOMENTUM (NOT REALLY UPDATED IF LAST_STEP=0)
	for (i = 0; i < num_d_parameters; i++)
	{
		if (!(d_parameters_filestream << *(d_potential_parameters + i) << endl))
		{
			log_filestream << "could not write the new momentum file updated parameters -> killing" << endl;
			log_filestream << "///////////////////END RE CODE/////////////////////" << endl << endl;
			exit(EXIT_FAILURE);
		}
	}
}

//SPECIFIC POTENTIAL FUNCTION DEFINITIONS
array_pair_and_num_elements optimize_ideal_cluster_potential(int last_step, array_pair_and_num_elements gr_data, 
	double *potential_parameters, double *d_potential_parameters, gromacs_settings_class gromacs_settings, 
	double *md_cutoff_pointer, double *gr_cutoff_pointer, double *unscaled_gradient_pointer, double *gr_convergence_pointer)
{
	//OPTIMIZATION PARAMETERS
	double e1 = *(potential_parameters + 0), a1 = *(potential_parameters + 1);
	double e2 = *(potential_parameters + 2), a2 = *(potential_parameters + 3);
	//NON_OPTIMIZED PARAMETERS SET BY USER
	int n1 = (int)(*(potential_parameters + 4)+0.5); //read in as a double so need to round and type cast
	int n2 = (int)(*(potential_parameters + 5)+0.5); //read in as a double so need to round and type cast
	double e1_e2_max = *(potential_parameters + 6); //this the the constraint on how big e1 + e2 can be
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
	//CONSTRAINED AMPLITUDE
	double scale_contact;
	double step_e1_contact, step_a1_contact;
	double step_e2_contact, step_a2_contact;
	double step_e1_overstep, step_a1_overstep;
	double step_e2_overstep, step_a2_overstep;
	double step_e1_overstep_p, step_a1_overstep_p;
	double step_e2_overstep_p, step_a2_overstep_p;

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
			r_II = (double)(i+1)*gromacs_settings.delta_r;

			//fill in rdf data
			gr_I = *(gr_data.array_1+i);
			gr_II = *(gr_data.array_1+i+1);
			gr_tgt_I = *(gr_data.array_2+i);
			gr_tgt_II = *(gr_data.array_2+i+1);

			//first set of potential derivatives
			dude1_I = -(1.0 / (exp(pow(-((2.0 + a1 - 2.0 * r_I) / a1), n1))*gromacs_settings.kB_T));
			duda1_I = (2.0 * ((e2*(double)n2*pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), -1 + n2)) / (a2*exp(pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), n2))) -
				(e1*(double)n1*pow(-((2.0 + a1 - 2.0 * r_I) / a1), -1 + n1)*(-1.0 + r_I)) /
				(pow(a1, 2)*exp(pow(-((2.0 + a1 - 2.0 * r_I) / a1), n1))))) / gromacs_settings.kB_T;
			dude2_I = 1.0 / (exp(pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), n2))*gromacs_settings.kB_T);
			duda2_I = (-2.0 * e2*(double)n2*pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), -1 + n2)*(1.0 + a1 - r_I)) /
				(pow(a2, 2)*exp(pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), n2))*gromacs_settings.kB_T);

			//second set of potential derivatives
			dude1_II = -(1.0 / (exp(pow(-((2.0 + a1 - 2.0 * r_II) / a1), n1))*gromacs_settings.kB_T));
			duda1_II = (2.0 * ((e2*(double)n2*pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_II) / a2), -1 + n2)) / (a2*exp(pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_II) / a2), n2))) -
				(e1*(double)n1*pow(-((2.0 + a1 - 2.0 * r_II) / a1), -1 + n1)*(-1.0 + r_II)) /
				(pow(a1, 2)*exp(pow(-((2.0 + a1 - 2.0 * r_II) / a1), n1))))) / gromacs_settings.kB_T;
			dude2_II = 1.0 / (exp(pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_II) / a2), n2))*gromacs_settings.kB_T);
			duda2_II = (-2.0 * e2*(double)n2*pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_II) / a2), -1 + n2)*(1.0 + a1 - r_II)) /
				(pow(a2, 2)*exp(pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_II) / a2), n2))*gromacs_settings.kB_T);

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

		//calculate unscaled gradient
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



		//check if constraint is satisfied and if not perform adjustment (project overstep vector onto constraint plane)
		if ((e1 + step_e1 + e2 + step_e2) >= e1_e2_max)
		{
			//amount to scale back by to just hit constraint
			if ((e1_e2_max > (e1 + e2)))
				scale_contact = (e1_e2_max - e1 - e2) / (step_e1 + step_e2);
			else
				scale_contact = 0.0; //if on boundary

			//constraint contact vector
			step_e1_contact = scale_contact * step_e1;
			step_a1_contact = scale_contact * step_a1;
			step_e2_contact = scale_contact * step_e2;
			step_a2_contact = scale_contact * step_a2;

			//constraint overstep vector
			step_e1_overstep = (1.0-scale_contact) * step_e1;
			step_a1_overstep = (1.0-scale_contact) * step_a1;
			step_e2_overstep = (1.0-scale_contact) * step_e2;
			step_a2_overstep = (1.0-scale_contact) * step_a2;

			//perform the projection of the overstep vector
			step_e1_overstep_p = (1.0 / 2.0) * (step_e1_overstep - step_e2_overstep);
			step_a1_overstep_p = step_a1_overstep * 1.0;
			step_e2_overstep_p = (1.0 / 2.0) * (-(step_e1_overstep - step_e2_overstep));
			step_a2_overstep_p = step_a2_overstep * 1.0;

			//new combined contacted and then projected vector update
			step_e1 = step_e1_contact + step_e1_overstep_p;
			step_a1 = step_a1_contact + step_a1_overstep_p;
			step_e2 = step_e2_contact + step_e2_overstep_p;
			step_a2 = step_a2_contact + step_a2_overstep_p;
		}


		//update the parameters and pass back
		e1 = e1 + /*grad_e1*/ step_e1; *(potential_parameters + 0) = e1;
		a1 = a1 + /*grad_a1*/ step_a1; *(potential_parameters + 1) = a1;
		e2 = e2 + /*grad_e2*/ step_e2; *(potential_parameters + 2) = e2;
		a2 = a2 + /*grad_a2*/ step_a2; *(potential_parameters + 3) = a2;

	}



	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////FIND THE GR CUTOFF////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	//NOTE: THIS IS NOT CURRENTLY BEING USED BUT I HAVE KEPT IT AS A OPTION IN CASE OF FUTURE DEVELOPMENTS
	
	//FIND GR CUTOFF BASED ON UPDATED (OR NOT) PARAMETERS USING KBT NORMALIZATION
	//NOTE!!! THIS USES NUM_TABLE_ENTRIES TO ENSURE IT GOES OUT FAR ENOUGH
	for (i = num_table_entries - 1; i > 0; i--)
	{
		r_I = (double)i*gromacs_settings.delta_r;

		dude1_I = -(1.0 / (exp(pow(-((2.0 + a1 - 2.0 * r_I) / a1), n1))*gromacs_settings.kB_T));
		duda1_I = (2.0 * ((e2*(double)n2*pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), -1 + n2)) / (a2*exp(pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), n2))) -
			(e1*(double)n1*pow(-((2.0 + a1 - 2.0 * r_I) / a1), -1 + n1)*(-1.0 + r_I)) /
			(pow(a1, 2)*exp(pow(-((2.0 + a1 - 2.0 * r_I) / a1), n1))))) / gromacs_settings.kB_T;
		dude2_I = 1.0 / (exp(pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), n2))*gromacs_settings.kB_T);
		duda2_I = (-2.0 * e2*(double)n2*pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), -1 + n2)*(1.0 + a1 - r_I)) /
			(pow(a2, 2)*exp(pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), n2))*gromacs_settings.kB_T);
		
		if (fabs(dude1_I) > gromacs_settings.gr_cutoff_magnitude ||
			fabs(duda1_I) > gromacs_settings.gr_cutoff_magnitude ||
			fabs(dude2_I) > gromacs_settings.gr_cutoff_magnitude ||
			fabs(duda2_I) > gromacs_settings.gr_cutoff_magnitude)
		{
			*(gr_cutoff_pointer) = r_I;
			break;
		}
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
				-(e1 / (exp(pow(2.0, n1)*pow((-1.0 - a1 / 2.0 + r) / a1, n1)))) +
				e2 / (exp(pow(2.0, n2)*pow((-1.0 - a1 - a2 / 2.0 + r) / a2, n2)));

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
				-(e1 / (exp(pow(2.0, n1)*pow((-1.0 - a1 / 2.0 + r) / a1, n1)))) +
				e2 / (exp(pow(2.0, n2)*pow((-1.0 - a1 - a2 / 2.0 + r) / a2, n2)));
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




	//AS OF NOW ALL THE NECCESSARY FILES FOR THE NEXT STEP SHOULD HAVE BEEN GENERATED
	//IN THE LAST STEP DIRECTORY WITH THE OUT FILE NAME FLAG FOR COPYING (AND RENAMING) TO
	//THE NEXT STEP DIRECTORY

	return table_data;
}

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
	///////////////////////////////////FIND THE GR CUTOFF////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	//NOTE: THIS IS NOT CURRENTLY BEING USED BUT I HAVE KEPT IT AS A OPTION IN CASE OF FUTURE DEVELOPMENTS

	//FIND GR CUTOFF BASED ON UPDATED (OR NOT) PARAMETERS USING KBT NORMALIZATION
	//NOTE!!! THIS USES NUM_TABLE_ENTRIES TO ENSURE IT GOES OUT FAR ENOUGH
	for (i = num_table_entries - 1; i > 0; i--)
	{
		r_I = (double)i*gromacs_settings.delta_r;

		dude1_I = -(1.0 / (exp(pow(-((2.0 + a1 - 2.0 * r_I) / a1), n1))*gromacs_settings.kB_T));
		duda1_I = (2.0 * ((e2*(double)n2*pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), -1 + n2)) / (a2*exp(pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), n2))) -
			(e1*(double)n1*pow(-((2.0 + a1 - 2.0 * r_I) / a1), -1 + n1)*(-1.0 + r_I)) /
			(pow(a1, 2)*exp(pow(-((2.0 + a1 - 2.0 * r_I) / a1), n1))))) / gromacs_settings.kB_T;
		dude2_I = 1.0 / (exp(pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), n2))*gromacs_settings.kB_T);
		duda2_I = (-2.0 * e2*(double)n2*pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), -1 + n2)*(1.0 + a1 - r_I)) /
			(pow(a2, 2)*exp(pow(-((2.0 + 2.0 * a1 + a2 - 2.0 * r_I) / a2), n2))*gromacs_settings.kB_T);

		if (fabs(dude1_I) > gromacs_settings.gr_cutoff_magnitude ||
			fabs(duda1_I) > gromacs_settings.gr_cutoff_magnitude ||
			fabs(dude2_I) > gromacs_settings.gr_cutoff_magnitude ||
			fabs(duda2_I) > gromacs_settings.gr_cutoff_magnitude)
		{
			*(gr_cutoff_pointer) = r_I;
			break;
		}
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
