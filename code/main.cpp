#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <limits>
//#include <format> // not available for GCC !!!
#include <cstdio>

#include "global.h"

#include "components.h"
#include "units.h"
#include "simulation_setup.h"


using namespace std;



/*
Return values of the program:
0	Normal execution, no errors occured
1	Wrong parameters
2	A required file was not found or error during database connections
3	Errors duing the simulation
4	Erroneous input files
*/

int main(int argc, char* argv[]) {
	
	//
	// parsing command line arguments
	//
	int scenario_id;
	if (argc == 1) {
		scenario_id = 1;
	} else if (argc == 2) {
		scenario_id = stoi(argv[1]);
	} else {
		cout << "Error: Wrong number of arguments!" << endl;
		cout << "Usage: main[-dbg] [scenario id]" << endl;
		return 1;
	}

	cout << "Initializing the simulation for scenario ID " << scenario_id << endl;

	//
	// creating variables and initializing the
	// static attributes in the classes
	//
	float expansion_matrix[16][16];
	Global::InitializeStaticVariables();

	//
	// open and parse the simulation scenario csv file
	//
	if (!configld::parse_scenario_file(scenario_id)) {
		return 2;
	}

	//
	// Load the expansion matrix
	//
	if (!expansion::load_expansion_matrix(expansion_matrix)) {
		return 2;
	}

	//
	// Check the expansion matrix
	//
	if (!expansion::verify_expansion_matrix(expansion_matrix)) {
		return 4;
	}

	//
	// Create the control units and load the data for the connected measurement units
	//
	if (!configld::load_data_from_central_database("../data/input/Merged_Information.db")) {
		return 2;
	}

	ComponentPV cPV(22.3);
	// this is not allowed anymore, as the object are deleted twice!!!
	//MeasurementUnit mu(1, NULL, 0/*meloID,string melo, locID*/);
	// do this:
	MeasurementUnit* mu = new MeasurementUnit(1, NULL, 0/*meloID,string melo, locID*/);
	mu->load_data("../data/input/SeparatedSmartMeterData/2.csv");

	//
	// Load central solar radation and wind profiles
	//

	//
	// bevore starting the simulation:
	// check if all global variables are set -> if not, error!
	//
	if (!global::all_variables_initialized() || !Global::AllVariablesInitialized()) {
		cout << "Some global variables are not initialized!" << endl;
		return 3;
	}
	// TODO: and then make all global variables final!

	//
	// Add expansion[s] to the control units
	//

	//
	// Run the simulation
	//
	
	//
	// clean up
	//
	MeasurementUnit::VacuumInstancesAndStaticVariables();
	ControlUnit::VacuumInstancesAndStaticVariables();
	Substation::VacuumInstancesAndStaticVariables();
	Global::DeleteStaticVariables();
	global::vacuum();

	return 0;
}

