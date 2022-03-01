#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <sqlite3.h>
#include <limits>
//#include <format> // not available for GCC !!!
#include <cstdio>

#include "global.h"

#include "components.h"
#include "units.h"


using namespace std;




static int callback(void* data, int argc, char** argv, char** colName) {
	/* argc: holds the number of results, argv: holds each value in array, colName: holds each column returned in array, */
	string comparisonStr1 = "n_timesteps";
	string comparisonStr2 = "n_substations";
	if (argc != 2) {
		cout << "Number of arguments not equal to 2 for one column!";
		return 1;
	}
	
	if        (comparisonStr1.compare(argv[0]) == 0) {
		global::n_timesteps        = stoi(argv[1]);
		global::n_timesteps_init   = true;
	} else if (comparisonStr2.compare(colName[0]) == 0) {
		global::n_substations      = stoi(argv[1]);
		global::n_substations_init = true;
	}

	return 0;
}

/*
Return values of the program:
0	Normal execution, no errors occured
1	Wrong parameters
2	A required file was not found or error during database connections
3	Errors duing the simulation
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

	cout << "Initializing the simulation!\n";

	//
	// open and parse the simulation scenario csv file
	//
	const char* scenarios_input_path = "../config/simulation_scenarios.csv";
	ifstream scenarios_input;
	scenarios_input.open(scenarios_input_path);
	if (!scenarios_input.good()) {
		cout << "Error when connecting to the simulation scenario file with path " << scenarios_input_path << endl;
		return 2;
	} else {
		string currLineString;
		getline( scenarios_input, currLineString ); // jump first line, as this is the header
		for (int r = 0; r < scenario_id; r++) {
			// iterate over every row
			getline( scenarios_input, currLineString );
			stringstream currLineStream( currLineString );
			string currLineSplitted[9];
			for (int col = 0; col < 9; col++) {
				// split this row on the ","
				getline( currLineStream, currLineSplitted[col], ',' );
			}
			// convert individual strings to int / float / char
			int currentLineID = stoi( currLineSplitted[0] );
			if (currentLineID == scenario_id) {
				global::ts_start_str    = new string( currLineSplitted[1] ); // copy constructor
				global::ts_end_str      = new string( currLineSplitted[2] ); // copy constructor
				global::tsteps_per_hour = stoi( currLineSplitted[3] );
				global::expansion_scenario_id = stoi( currLineSplitted[4] );
				global::exp_pv_kWp      = stof( currLineSplitted[5] );
				global::exp_bess_kW     = stof( currLineSplitted[6] );
				global::exp_bess_kWh    = stof( currLineSplitted[7] );
				global::exp_bess_start_soc    = stof( currLineSplitted[8] );

				global::ts_start_str_init    = true;
				global::ts_end_str_init      = true;
				global::tsteps_per_hour_init = true;
				global::expansion_scenario_id_init = true;
				global::exp_pv_kWp_init      = true;
				global::exp_bess_kW_init     = true;
				global::exp_bess_kWh_init    = true;
				global::exp_bess_start_soc_init    = true;

				break;
			} else {
				continue;
			}
		}
		scenarios_input.close();
	}

	//
	// bevore starting the simulation:
	// check if all global variables are set -> if not, error!
	//
	if (!global::all_variables_initialized()) {
		cout << "Some variables are not initialized!" << endl;
		return 3;
	}
	// TODO: and then make all global variables final!

	//
	// Load the expansion matrix
	//

	//
	// Check the expansion matrix
	//

	//
	// Load central solar radation and wind profiles
	//

	//
	// Create the control units and load the data for the connected measurement units
	//

	sqlite3* dbcon;
	int rc = sqlite3_open("/home/daniel/Daten/current/31_export_for_simulation/Merged_Information.db", &dbcon);

	if (rc == 0) {
		cout << "Database opened" << endl;
		
		string sql_query = "SELECT key, value FROM general_data_information;";
		char* sqlErrorMsg;
		int ret_val = sqlite3_exec(dbcon, sql_query.c_str(), callback, NULL, &sqlErrorMsg);
		if (ret_val != 0) {
			cout << "Error when reading the SQL-Table: " << sqlErrorMsg;
			sqlite3_free(sqlErrorMsg);
		}
		
		sqlite3_close(dbcon);
	} else {
		cout << "Error when connecting to the database" << endl;
		return 2;
	}
	
	if (!global::n_timesteps_init) {
		cout << "g_n_timesteps is not initialized! Stopping execution." << endl;
		return 3;
	}

	ComponentPV cPV(22.3);
	MeasurementUnit mu(0, NULL, 0/*meloID,string melo, locID*/);
	mu.load_data("/home/daniel/Daten/current/31_export_for_simulation/SeparatedSmartMeterData/2.csv");

	//
	// Add expansion[s] to the control units
	//

	//
	// Run the simulation
	//
	
	//
	// clean up
	//
	delete global::ts_start_str;  global::ts_start_str = 0;
	delete global::ts_end_str;    global::ts_end_str   = 0;

	return 0;
}

