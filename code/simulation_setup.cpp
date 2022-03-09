#include "simulation_setup.h"

using namespace expansion;


#include <fstream>
#include <iomanip>
#include <iostream>
#include <sqlite3.h>
#include <sstream>
#include <string>
#include <vector>

using namespace std;


#include "global.h"
#include "units.h"



//
// open and parse the simulation scenario csv file
//
bool configld::parse_scenario_file(int scenario_id) {
    const char* scenarios_input_path = "../config/simulation_scenarios.csv";
	ifstream scenarios_input;
	scenarios_input.open(scenarios_input_path);
	if (!scenarios_input.good()) {
		cerr << "Error when connecting to the simulation scenario file with path " << scenarios_input_path << endl;
		return false;
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
				Global::set_ts_start_str(   new string( currLineSplitted[1] )); // copy constructor
				Global::set_ts_end_str(     new string( currLineSplitted[2] )); // copy constructor
				Global::set_tsteps_per_hour(stoi( currLineSplitted[3] ));
				Global::set_expansion_scenario_id(stoi( currLineSplitted[4] ));
				Global::set_exp_pv_kWp(        stof( currLineSplitted[5] ));
				Global::set_exp_bess_kW(       stof( currLineSplitted[6] ));
				Global::set_exp_bess_kWh(      stof( currLineSplitted[7] ));
				Global::set_exp_bess_start_soc(stof( currLineSplitted[8] ));
				break;
			} else {
				continue;
			}
		}
		scenarios_input.close();
	}
    return true;
}


int load_data_from_central_database_callbackA(void* data, int argc, char** argv, char** colName) {
	/* 
	 * This is the callback function for reading global attributes
	 * (like timesteps and substations available in the data)
	 * 
	 * Parameter information taken from documentation:
	 *   argc: holds the number of results
	 *   argv: holds each value in array
	 *   colName: holds each column returned in array
	 */
	string comparisonStr1 = "n_timesteps";
	string comparisonStr2 = "n_substations";
	string comparisonStr3 = "n_control_units";
	string comparisonStr4 = "n_measurement_units";
	if (argc != 2) {
		cerr << "Number of arguments not equal to 2 for one row!" << endl;
		return 1;
	}
	
	if        (comparisonStr1.compare(argv[0]) == 0) {
		Global::set_n_timesteps(   stoi(argv[1]) );
	} else if (comparisonStr2.compare(argv[0]) == 0) {
		Global::set_n_substations( stoi(argv[1]) );
	} else if (comparisonStr3.compare(argv[0]) == 0) {
		Global::set_n_CUs(          stoi(argv[1]) );
	} else if (comparisonStr4.compare(argv[0]) == 0) {
		Global::set_n_MUs(          stoi(argv[1]) );
	}

	return 0;
}
int load_data_from_central_database_callbackB(void* data, int argc, char** argv, char** colName) {
	/* 
	 * This is the callback function for loading the time indices
	 */
    static int callcounter = 1;
    int pos = callcounter - 1; // current position in the array is one behind the count of calls
	if (argc != 3) {
		cerr << "Number of arguments not equal to 3 for one row!" << endl;
		return 1;
	}
	int current_time_index = stoi(argv[0]);
	if (current_time_index != callcounter) {
		cerr << "Time indices are not ordered sequentially!" << endl;
		return 1;
	}
    global::time_timestep_id[pos] = current_time_index;
    global::time_localtime_str->push_back(string(argv[1]));
    global::time_localtimezone_str->push_back(string(argv[2]));
    callcounter++;
	return 0;
}
int load_data_from_central_database_callbackC(void* data, int argc, char** argv, char** colName) {
	/* 
	 * This is the callback function for geeting information about the substations
	 * out of the database.
	 * This function also creates the substations.
	 */
	if (argc != 2) {
		cerr << "Number of arguments not equal to 2 for one row!" << endl;
		return 1;
	}
	int current_station_id = stoi(argv[0]);
	string* stationName = new string(argv[1]);
	try {
    	/*Substation* newSubstation = */ new Substation(current_station_id, stationName);
	} catch (const char* msg) {
		return 1;
	}
	return 0;
}
int load_data_from_central_database_callbackD(void* data, int argc, char** argv, char** colName) {
	/* 
	 * This is the callback function for geeting information about the control units
	 * out of the database.
	 * This function also creates the control units.
	 */
	if (argc != 2) {
		cerr << "Number of arguments not equal to 2 for one row!" << endl;
		return 1;
	}
	int current_cu_id = stoi(argv[0]);
	int conn_to_subst_id = stoi(argv[1]);
	try {
    	/*ControlUnit* newCU = */ new ControlUnit(current_cu_id, conn_to_subst_id);
	} catch (const char* msg) {
		return 1;
	}
	return 0;
}
int load_data_from_central_database_callbackE(void* data, int argc, char** argv, char** colName) {
	/* 
	 * This is the callback function for geeting information about the measurement units
	 * out of the database.
	 * This function also creates the measurement units.
	 * 
	 * Columns:
	 * 0        1       2     3           4           5
	 * MELO_ID, UnitID, MELO, has_demand, has_feedin, has_pv_residential
	 * 6                  7         8       9        10
	 * has_pv_open_space, has_bess, has_hp, has_chp, LocID
	 */
	if (argc != 11) {
		cerr << "Number of arguments not equal to 10 for one row!" << endl;
		return 1;
	}
	int current_mu_id = stoi(argv[0]);
	int conn_to_unitID = stoi(argv[1]);
	string* melo_str   = new string(argv[2]);
	bool has_demand    = argv[3][0] == '1';
	bool has_feedin    = argv[4][0] == '1';
	bool has_pv_resid  = argv[5][0] == '1';
	bool has_pv_opens  = argv[6][0] == '1';
	bool has_bess      = argv[7][0] == '1';
	bool has_hp        = argv[8][0] == '1';
	bool has_wb        = false; // TODO: wallboxes not implemented jet
	bool has_chp       = argv[9][0] == '1';
	int locID          = stoi(argv[10]);

	stringstream data_input_path;
	data_input_path << "../data/input/SeparatedSmartMeterData/";
	data_input_path << current_mu_id;
	data_input_path << ".csv";

	try {
		MeasurementUnit* newMU =
			new MeasurementUnit(current_mu_id, conn_to_unitID, melo_str, locID,
								has_demand, has_feedin, has_pv_resid, has_pv_opens,
								has_bess,   has_hp,     has_wb,       has_chp);
		newMU->load_data(data_input_path.str().c_str());
	} catch (const char* msg) {
		return 1;
	}
	return 0;
}
bool configld::load_data_from_central_database(const char* filepath) {
    sqlite3* dbcon;
	int rc = sqlite3_open(filepath, &dbcon);

	if (rc == 0) {
		cout << "Retrieving information from Merged_Information.db ..." << endl;
		
        //
        // Load general data
        //
		string sql_queryA = "SELECT key, value FROM general_data_information;";
		char* sqlErrorMsgA;
		int ret_valA = sqlite3_exec(dbcon, sql_queryA.c_str(), load_data_from_central_database_callbackA, NULL, &sqlErrorMsgA);
		if (ret_valA != 0) {
			cerr << "Error when reading the SQL-Table: " << sqlErrorMsgA;
			sqlite3_free(sqlErrorMsgA);
			return false;
		}

        //
        // Initialize global time list
        //
        global::time_timestep_id = new int[Global::get_n_timesteps()];
        global::time_localtime_str = new vector<string>();
        global::time_localtimezone_str = new vector<string>();
        //
        // Load time indices
        //
        string sql_queryB = "SELECT TimestepID, local_time, local_time_zone FROM time_indices ORDER BY TimestepID;";
		char* sqlErrorMsgB;
		int ret_valB = sqlite3_exec(dbcon, sql_queryB.c_str(), load_data_from_central_database_callbackB, NULL, &sqlErrorMsgB);
		if (ret_valB != 0) {
			cerr << "Error when reading the SQL-Table: " << sqlErrorMsgB;
			sqlite3_free(sqlErrorMsgB);
			return false;
		}
        //
        global::time_info_init = true;
		
		//
        // Initialize global list of units (i.e. control unit, measurement unit and substation)
        //
		Substation::InitializeStaticVariables(Global::get_n_substations());
		ControlUnit::InitializeStaticVariables(Global::get_n_CUs());
		MeasurementUnit::InitializeStaticVariables(Global::get_n_MUs());
        //
        // Load component information and create units
        // 1. per substation
        // 2. per control unit
        // 3. per measurement unit
        //
		// 1. substations
		string sql_queryC = "SELECT substation_id, substation_name FROM substation_information ORDER BY substation_id;";
		char* sqlErrorMsgC;
		int ret_valC = sqlite3_exec(dbcon, sql_queryC.c_str(), load_data_from_central_database_callbackC, NULL, &sqlErrorMsgC);
		if (ret_valC != 0) {
			cerr << "Error when reading the SQL-Table: " << sqlErrorMsgC;
			sqlite3_free(sqlErrorMsgC);
			return false;
		}
		// 2. CUs
		string sql_queryD = "SELECT UnitID, substation_id FROM melo_information ORDER BY UnitID;";
		char* sqlErrorMsgD;
		int ret_valD = sqlite3_exec(dbcon, sql_queryD.c_str(), load_data_from_central_database_callbackD, NULL, &sqlErrorMsgD);
		if (ret_valD != 0) {
			cerr << "Error when reading the SQL-Table: " << sqlErrorMsgD;
			sqlite3_free(sqlErrorMsgD);
			return false;
		}
		// 3. MUs
		string sql_queryE = "SELECT MELO_ID, UnitID, MELO, has_demand, has_feedin, has_pv_residential, has_pv_open_space, has_bess, has_hp, has_chp, LocID FROM melo_information ORDER BY MELO_ID;";
		char* sqlErrorMsgE;
		int ret_valE = sqlite3_exec(dbcon, sql_queryE.c_str(), load_data_from_central_database_callbackE, NULL, &sqlErrorMsgE);
		if (ret_valE != 0) {
			cerr << "Error when reading the SQL-Table: " << sqlErrorMsgE;
			sqlite3_free(sqlErrorMsgE);
			return false;
		}
        
		//
        // Load address data
        //
		
		sqlite3_close(dbcon);
        return true;
	} else {
		cerr << "Error when connecting to the database" << endl;
		return false;
	}
}





bool expansion::is_expansion_combination_possible(int current_scenario_number, int b) {
	/*
	This function returns if a expansion from current_scenario_number to b (i.e. the future scenario number) is possible or not
	*/
	if (current_scenario_number == EXPMAT_POS_) {
		if (b == EXPMAT_POS_PV || b == EXPMAT_POS_BS || b == EXPMAT_POS_HP || b == EXPMAT_POS_WB || b == EXPMAT_POS_PV_BS || b == EXPMAT_POS_PV_HP || b == EXPMAT_POS_PV_WB || b == EXPMAT_POS_BS_HP || b == EXPMAT_POS_BS_WB || b == EXPMAT_POS_HP_WB || b == EXPMAT_POS_PV_BS_HP || b == EXPMAT_POS_PV_BS_WB || b == EXPMAT_POS_PV_HP_WB || b == EXPMAT_POS_BS_HP_WB || b == EXPMAT_POS_PV_BS_HP_WB)
			return true;
	} else if (current_scenario_number == EXPMAT_POS_PV) {
		if (b == EXPMAT_POS_PV_BS || b == EXPMAT_POS_PV_HP || b == EXPMAT_POS_PV_WB || b == EXPMAT_POS_PV_BS_HP || b == EXPMAT_POS_PV_BS_WB || b == EXPMAT_POS_PV_HP_WB || b == EXPMAT_POS_PV_BS_HP_WB)
			return true;
	} else if (current_scenario_number == EXPMAT_POS_BS) {
		if (b == EXPMAT_POS_PV_BS || b == EXPMAT_POS_BS_HP || b == EXPMAT_POS_BS_WB || b == EXPMAT_POS_PV_BS_HP || b == EXPMAT_POS_PV_BS_WB || b == EXPMAT_POS_BS_HP_WB || b == EXPMAT_POS_PV_BS_HP_WB)
			return true;
	} else if (current_scenario_number == EXPMAT_POS_HP) {
		if (b == EXPMAT_POS_PV_HP || b == EXPMAT_POS_BS_HP || b == EXPMAT_POS_HP_WB || b == EXPMAT_POS_PV_BS_HP || b == EXPMAT_POS_PV_HP_WB || b == EXPMAT_POS_BS_HP_WB || b == EXPMAT_POS_PV_BS_HP_WB)
			return true;
	} else if (current_scenario_number == EXPMAT_POS_WB) {
		if (b == EXPMAT_POS_PV_WB || b == EXPMAT_POS_BS_WB || b == EXPMAT_POS_HP_WB || b == EXPMAT_POS_PV_BS_WB || b == EXPMAT_POS_PV_HP_WB || b == EXPMAT_POS_BS_HP_WB || b == EXPMAT_POS_PV_BS_HP_WB)
			return true;
	} else if (current_scenario_number == EXPMAT_POS_PV_BS) {
		if (b == EXPMAT_POS_PV_BS_HP || b == EXPMAT_POS_PV_BS_WB || b == EXPMAT_POS_PV_BS_HP_WB)
			return true;
	} else if (current_scenario_number == EXPMAT_POS_PV_HP) {
		if (b == EXPMAT_POS_PV_BS_HP || b == EXPMAT_POS_PV_HP_WB || b == EXPMAT_POS_PV_BS_HP_WB)
			return true;
	} else if (current_scenario_number == EXPMAT_POS_PV_WB) {
		if (b == EXPMAT_POS_PV_BS_WB || b == EXPMAT_POS_PV_HP_WB || b == EXPMAT_POS_PV_BS_HP_WB)
			return true;
	} else if (current_scenario_number == EXPMAT_POS_BS_HP) {
		if (b == EXPMAT_POS_PV_BS_HP || b == EXPMAT_POS_BS_HP_WB || b == EXPMAT_POS_PV_BS_HP_WB)
			return true;
	} else if (current_scenario_number == EXPMAT_POS_BS_WB) {
		if (b == EXPMAT_POS_PV_BS_WB || b == EXPMAT_POS_BS_HP_WB || b == EXPMAT_POS_PV_BS_HP_WB)
			return true;
	} else if (current_scenario_number == EXPMAT_POS_HP_WB) {
		if (b == EXPMAT_POS_PV_HP_WB || b == EXPMAT_POS_BS_HP_WB || b == EXPMAT_POS_PV_BS_HP_WB)
			return true;
	} else if (current_scenario_number == EXPMAT_POS_PV_BS_HP) {
		if (b == EXPMAT_POS_PV_BS_HP_WB)
			return true;
	} else if (current_scenario_number == EXPMAT_POS_PV_BS_WB) {
		if (b == EXPMAT_POS_PV_BS_HP_WB)
			return true;
	} else if (current_scenario_number == EXPMAT_POS_PV_HP_WB) {
		if (b == EXPMAT_POS_PV_BS_HP_WB)
			return true;
	} else if (current_scenario_number == EXPMAT_POS_BS_HP_WB) {
		if (b == EXPMAT_POS_PV_BS_HP_WB)
			return true;
	} else if (current_scenario_number == EXPMAT_POS_PV_BS_HP_WB) {
		return false;
	} else {
		cerr << "Warning: A wrong value has been passed to function is_expansion_combination_possible." << endl;
	}
	return false;
}


/*
Loads the expansion matrix into the first argument.

Returns false if an error occurs, else true.
*/
bool expansion::load_expansion_matrix(float expansion_matrix[16][16]) {
	ifstream expmat_input;
	stringstream expmat_input_path;
	expmat_input_path << "../config/expansion_scenarios/";
	expmat_input_path << setw(4) << setfill('0') << Global::get_expansion_scenario_id();
	expmat_input_path << ".csv";
	expmat_input.open(expmat_input_path.str());
	if (!expmat_input.good()) {
		cerr << "Error when loading the expansion matrix with path " << expmat_input_path.str() << endl;
		return false;
	} else {
		cout << "Opening expansion matrix file " << expmat_input_path.str() << endl;
		string currLineString;
		getline( expmat_input, currLineString ); // jump first line, as this is the header
		for (int i = 0; i < 16; i++) {
			getline( expmat_input, currLineString );
			stringstream currLineStream( currLineString );
			string currLineSplittedElement;
			float current_value;
			getline( currLineStream, currLineSplittedElement, ',' ); // jump first row, as this is the index
			for (int col = 0; col < 16; col++) {
				current_value = 0.0;
				// split this row on the ","
				if ( getline( currLineStream, currLineSplittedElement, ',' ) ) {
					if (currLineSplittedElement.length() > 0) {
						// remove \r if it occurs at the end
						if (currLineSplittedElement.back() == '\r') {
							currLineSplittedElement.erase(currLineSplittedElement.length()-1, 1);
						}
						// convert the value if there still is a element remaining
						if (currLineSplittedElement.length() > 0)
							current_value = stof( currLineSplittedElement );
					}
				} else {
					cout << "Warning: End of line reached before 16 elements are parsed in the expansion scenario file!" << endl;
				}
				expansion_matrix[i][col] = current_value;
			}
		}
		expmat_input.close();
	}
	#ifdef DEBUG_EXTRA_OUTPUT
	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 16; j++) {
			cout << setw(4) << setfill(' ') << expansion_matrix[i][j] << " ";
		}
		cout << endl;
	}
	#endif
    return true;
}




bool expansion::verify_expansion_matrix(float expansion_matrix[16][16]) {
    //
	// Check the expansion matrix
	// 1. 0 for impossible combinations (and diagonal, as this is computed below)
	//    and all values between 0 and 1
	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 16; j++) {
			if (expansion::is_expansion_combination_possible(i,j)) {
				// check if between 0 and 1
				if (expansion_matrix[i][j] < 0.0 || expansion_matrix[i][j] > 1.0) {
					cerr << "Error in expansion matrix: value at position " << i << ", " << j << " is greater 1 or below 0!" << endl;
					return false;
				}
			} else {
				// check if 0
				if (expansion_matrix[i][j] != 0.0) {
					cerr << "Error in expansion matrix: value at position " << i << ", " << j << " is " << expansion_matrix[i][j] << " instead of 0!" << endl;
					return false;
				}
			}
		}
	}
	// 2. calculate the diagonal values (i.e. percentage of unchanged values)
	//    and check if row sum is smaller or equal 1
	for (int i = 0; i < 16; i++) {
		float row_sum = 0.0;
		for (int j = i+1; j < 16; j++) {
			row_sum += expansion_matrix[i][j];
		}
		if (row_sum > 1) {
			cerr << "Error in expansion matrix: row sum > 1 for line " << i << endl;
			return false;
		}
		expansion_matrix[i][i] = 1 - row_sum;
	}

    return true;
}

