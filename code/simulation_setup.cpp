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
    return true;
}


int load_data_from_central_database_callbackA(void* data, int argc, char** argv, char** colName) {
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
	} else if (comparisonStr2.compare(argv[0]) == 0) {
		global::n_substations      = stoi(argv[1]);
		global::n_substations_init = true;
	}

	return 0;
}
int load_data_from_central_database_callbackB(void* data, int argc, char** argv, char** colName) {
	/* argc: holds the number of results, argv: holds each value in array, colName: holds each column returned in array, */
    static int callcounter = 1;
    int pos = callcounter - 1; // current position in the array is one behind the count of calls
	if (argc != 3) {
		cout << "Number of arguments not equal to 3 for one column!";
		return 1;
	}
    global::time_timestep_id[pos] = stoi(argv[0]);
    global::time_localtime_str->push_back(string(argv[1]));
    global::time_localtimezone_str->push_back(string(argv[2]));
    callcounter++;
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
			cout << "Error when reading the SQL-Table: " << sqlErrorMsgA;
			sqlite3_free(sqlErrorMsgA);
		}

        //
        // Initialize global time list
        //
        global::time_timestep_id = new int[global::n_timesteps];
        global::time_localtime_str = new vector<string>();
        global::time_localtimezone_str = new vector<string>();
        //
        // Load time indices
        //
        string sql_queryB = "SELECT TimestepID, local_time, local_time_zone FROM time_indices ORDER BY TimestepID;";
		char* sqlErrorMsgB;
		int ret_valB = sqlite3_exec(dbcon, sql_queryB.c_str(), load_data_from_central_database_callbackB, NULL, &sqlErrorMsgB);
		if (ret_valB != 0) {
			cout << "Error when reading the SQL-Table: " << sqlErrorMsgB;
			sqlite3_free(sqlErrorMsgB);
		}
        //
        global::time_info_init = true;
        // TODO: Check time indices, are they continuos?
		
        //
        // Load component information and create units
        // 1. per substation
        // 2. per control unit
        // 3. per measurement unit
        //
        
		//
        // Load address data
        //
		
		sqlite3_close(dbcon);
        return true;
	} else {
		cout << "Error when connecting to the database" << endl;
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
	expmat_input_path << setw(4) << setfill('0') << global::expansion_scenario_id;
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

