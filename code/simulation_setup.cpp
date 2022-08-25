#include "simulation_setup.h"

using namespace expansion;


#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <memory>
#include <numeric>
#include <sqlite3.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

using namespace std;
namespace bpt = boost::property_tree;


#include "global.h"
#include "helper.h"
#include "units.h"


//
// loads the global config file
//
bool configld::load_config_file(int scenario_id, string& filepath) {
    //
    // parse json
    bpt::ptree tree_root;
    try {
        bpt::read_json(filepath, tree_root);
    } catch (bpt::json_parser_error& j) {
        cerr << "Error when reading json file: " << j.what() << endl;
        return false;
    }
    try {
        string str_data_ipt = "";  bool str_data_ipt_set = false;
        string str_data_opt = "";  bool str_data_opt_set = false;
        string start_str    = "";  bool start_str_set    = false;
        string end_str      = "";  bool end_str_set      = false;
        int    ts_per_hour  = 1;   bool ts_per_hour_set  = false;
        int    expansionID  = 0;   bool expansionID_set  = false;
        float  exp_bs_kW    = 0.0; bool exp_bs_kW_set    = false;
        float  exp_bs_kWh   = 0.0; bool exp_bs_kWh_set   = false;
        float  exp_bs_iSOC  = 0.0; bool exp_bs_iSOC_set  = false;
        float  os_pv_kWp    = 0.0; bool os_pv_kWp_set    = false;
        float  os_wind_kWp  = 0.0; bool os_wind_kWp_set  = false;
        // variables for PV expansion
        float  exp_pv_kWp_static  = 0.0; bool exp_pv_kWp_static_set   = false;
        float  exp_pv_kWp_m2_roof = 0.0; bool exp_pv_kWp_m2_roof_set  = false;
        float  exp_pv_min_kWp     = 0.0; bool exp_pv_min_kWp_set      = false;
        float  exp_pv_max_kWp     = 0.0; bool exp_pv_max_kWp_set      = false;
        bool   exp_pv_mode_static = false; bool exp_pv_mode_static_set= false;

        //
        // read default values
        for (auto& scenario_dict_all : tree_root.get_child("Default Scenario Values")) {
            string element_name = scenario_dict_all.first;
            if ( element_name.compare("data input path") == 0 ) {
                str_data_ipt    = scenario_dict_all.second.get_value<string>();
                str_data_ipt_set= true;
            } else if ( element_name.compare("data output path") == 0 ) {
                str_data_opt    = scenario_dict_all.second.get_value<string>();
                str_data_opt_set= true;
            } else if ( element_name.compare("start") == 0 ) {
                start_str       = scenario_dict_all.second.get_value<string>();
                start_str_set   = true;
            } else if ( element_name.compare("end") == 0 ) {
                end_str         = scenario_dict_all.second.get_value<string>();
                end_str_set     = true;
            } else if ( element_name.compare("time steps per hour") == 0 ) {
                ts_per_hour     = scenario_dict_all.second.get_value<int>();
                ts_per_hour_set = true;
            }
        }

        //
        // search the correct scenario dictionary
        // and read all variables from there, overwrite defaults if it necessary
        bool scenario_found = false;
        for (auto& scenario_dict_all : tree_root.get_child("Scenarios")) {
            auto scenario_dict = scenario_dict_all.second;
            // if we have found the correct entry ...
            if (scenario_dict.get<int>("id") == scenario_id) {
                // ... we read all variables
                for (auto& s : scenario_dict) {
                    string element_name = s.first;
                    if ( element_name.compare("data input path") == 0 ) {
                        str_data_ipt    = scenario_dict.get<string>("data input path");
                        str_data_ipt_set= true;
                    } else if ( element_name.compare("data output path") == 0 ) {
                        str_data_opt    = scenario_dict.get<string>("data output path");
                        str_data_opt_set= true;
                    } else if ( element_name.compare("start") == 0 ) {
                        start_str       = scenario_dict.get<string>("start");
                        start_str_set   = true;
                    } else if ( element_name.compare("end") == 0 ) {
                        end_str         = scenario_dict.get<string>("end");
                        end_str_set     = true;
                    } else if ( element_name.compare("time steps per hour") == 0 ) {
                        ts_per_hour     = scenario_dict.get<int>("time steps per hour");
                        ts_per_hour_set = true;
                    } else if ( element_name.compare("expansion id") == 0 ) {
                        expansionID     = scenario_dict.get<int>("expansion id");
                        expansionID_set = true;
                    } else if ( element_name.compare("expansion BS P in kW") == 0 ) {
                        exp_bs_kW       = scenario_dict.get<float>("expansion BS P in kW");
                        exp_bs_kW_set   = true;
                    } else if ( element_name.compare("expansion BS E in kWh") == 0 ) {
                        exp_bs_kWh      = scenario_dict.get<float>("expansion BS E in kWh");
                        exp_bs_kWh_set  = true;
                    } else if ( element_name.compare("expansion BS initial SOC") == 0 ) {
                        exp_bs_iSOC     = scenario_dict.get<float>("expansion BS initial SOC");
                        exp_bs_iSOC_set = true;
                    } else if ( element_name.compare("open space PV kWp") == 0 ) {
                        os_pv_kWp       = scenario_dict.get<float>("open space PV kWp");
                        os_pv_kWp_set   = true;
                    } else if ( element_name.compare("open space wind kWp") == 0 ) {
                        os_wind_kWp     = scenario_dict.get<float>("open space wind kWp");
                        os_wind_kWp_set = true;
                    } else if ( element_name.compare("expansion PV kWp static") == 0 ) {
                        exp_pv_kWp_static      = scenario_dict.get<float>("expansion PV kWp static");
                        exp_pv_kWp_static_set  = true;
                    } else if ( element_name.compare("expansion PV min kWp for section usage") == 0 ) {
                        exp_pv_min_kWp         = scenario_dict.get<float>("expansion PV min kWp for section usage");
                        exp_pv_min_kWp_set     = true;
                    } else if ( element_name.compare("expansion PV max inst kWp per section") == 0 ) {
                        exp_pv_max_kWp         = scenario_dict.get<float>("expansion PV max inst kWp per section");
                        exp_pv_max_kWp_set     = true;
                    } else if ( element_name.compare("expansion PV kWp per roof area in m2") == 0 ) {
                        exp_pv_kWp_m2_roof     = scenario_dict.get<float>("expansion PV kWp per roof area in m2");
                        exp_pv_kWp_m2_roof_set = true;
                    } else if ( element_name.compare("expansion PV kWp static mode") == 0 ) {
                        exp_pv_mode_static     = scenario_dict.get<bool>("expansion PV kWp static mode");
                        exp_pv_mode_static_set = true;
                    }
                }
                //
                scenario_found = true;
                break;
            }
        }
        #ifdef DEBUG
        cout << "str_data_ipt = " << str_data_ipt << endl;
        cout << "str_data_opt = " << str_data_opt << endl;
        cout << "start_str = " << start_str << endl;
        cout << "end_str = " << end_str << endl;
        cout << "ts_per_hour = " << ts_per_hour << endl;
        #endif

        //
        // is a parameter variation selected?
        // if yes, load the variables to vary
        bool pvar_scenario_found = false;
        if (Global::is_parameter_variation()) {
            vector<pair<string,shared_ptr<vector<float>>>> parameter_var_list; // this list holds the parsed parameters (string) with values (vector<float>) for the parameter variation
            for (auto& scenario_dict_all : tree_root.get_child("Parameter Variation")) {
                auto scenario_dict = scenario_dict_all.second;
                // if we have found the correct entry
                if (scenario_dict.get<int>("id") == Global::get_parameter_varID()) {
                    // ... we read all variables
                    for (auto& s : scenario_dict.get_child("variations")) {
                        // read variables from file
                        string varname = s.second.get<string>("variable");
                        float  step    = s.second.get<float>("step");
                        float  rStart  = s.second.get<float>("range start");
                        float  rEnd    = s.second.get<float>("range end");
                        // generate linspace
                        if (step < 0 || rEnd <= rStart) {
                            cerr << "Error when parsing parameter variation!" << endl;
                            cerr << "step < 0 || rEnd <= rStart" << endl;
                            return false;
                        }
                        shared_ptr<vector<float>> linspace ( new vector<float>() );
                        float cVal = rStart;
                        while (cVal <= rEnd) {
                            linspace->push_back(cVal);
                            #ifdef DEBUG
                            cout << "Linspace step with value " << cVal << endl;
                            #endif
                            cVal += step;
                        }
                        // add linspace (with varname) to list of parameter variations
                        parameter_var_list.emplace_back(varname,linspace);
                    }
                    //
                    pvar_scenario_found = true;
                    break;
                }
            }
            //
            // now, if the parameter variation scenario was found
            // we generate all combinations that should be simulated;
            // i.e. we compute the cartesian product
            if (pvar_scenario_found) {
                global::parameter_var_list = cartesian_product(parameter_var_list);
            }

            #ifdef DEBUG
            // test only
            for (auto elem : *(global::parameter_var_list)) {
                for (auto elem1 : elem) {
                    cout << setw(14) << setfill(' ') << elem1.first;
                }
                cout << endl;
            }
            for (auto elem : *(global::parameter_var_list)) {
                for (auto elem1 : elem) {
                    cout << setw(6) << setfill(' ') << elem1.second;
                }
                cout << endl;
            }
            // test end
            #endif
        }

        //
        // check, if path ends with an "/", add it, if not
        if (str_data_ipt.back() != '/') {
            str_data_ipt += "/";
        }
        if (str_data_opt.back() != '/') {
            str_data_opt += "/";
        }

        //
        // Finally, add to global variable collection
        if (scenario_found) {
            if (str_data_ipt_set) Global::set_input_path(str_data_ipt);
            if (str_data_opt_set) Global::set_output_path(str_data_opt);
            // for start / end time: convert value bevore setting global variable
            if (start_str_set && end_str_set) {
                struct tm* tm_start = new struct tm;
                struct tm* tm_end   = new struct tm;
                stringstream stream_val_t_start( start_str );
                stringstream stream_val_t_end(   end_str );
                stream_val_t_start >> get_time(tm_start, "%Y-%m-%d %H:%M:%S");
                stream_val_t_end   >> get_time(tm_end,   "%Y-%m-%d %H:%M:%S");
                Global::set_ts_start_tm( tm_start );
                Global::set_ts_end_tm(   tm_end   );
            }
            // other values
            if (ts_per_hour_set)  Global::set_tsteps_per_hour(ts_per_hour);
            if (expansionID_set)  Global::set_expansion_scenario_id(expansionID);
            if (exp_bs_kW_set)    Global::set_exp_bess_kW(exp_bs_kW);
            if (exp_bs_kWh_set)   Global::set_exp_bess_kWh(exp_bs_kWh);
            if (exp_bs_iSOC_set)  Global::set_exp_bess_start_soc(exp_bs_iSOC);
            if (os_pv_kWp_set)    Global::set_open_space_pv_kWp(os_pv_kWp);
            if (os_wind_kWp_set)  Global::set_wind_kWp(os_wind_kWp);
            if (exp_pv_mode_static_set) Global::set_exp_pv_mode(exp_pv_mode_static);
            if (exp_pv_kWp_static_set)  Global::set_exp_pv_kWp_static(exp_pv_kWp_static);
            if (exp_pv_kWp_m2_roof_set) Global::set_exp_pv_kWp_per_m2(exp_pv_kWp_m2_roof);
            if (exp_pv_min_kWp_set)     Global::set_exp_pv_min_kWp_roof_sec(exp_pv_min_kWp);
            if (exp_pv_max_kWp_set)     Global::set_exp_pv_max_kWp_roof_sec(exp_pv_max_kWp);
            //
            // change current working dir to the location of the config file
            filesystem::path config_fp = filepath;
            filesystem::current_path( config_fp.parent_path() );
            //
            //
            return true;
        }

    } catch (bpt::ptree_bad_path& j) {
        cerr << "Error when parsing json file: " << j.what() << endl;
        return false;
    }

    cerr << "Scenario ID " << scenario_id << " not found in the simulation configuration JSON file!" << endl;
    return false; // as we have not found what we searched
}

/*
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
            string currLineSplitted[11];
            for (int col = 0; col < 11; col++) {
                // split this row on the ","
                getline( currLineStream, currLineSplitted[col], ',' );
            }
            // convert individual strings to int / float / char
            int currentLineID = stoi( currLineSplitted[0] );
            if (currentLineID == scenario_id) {
                // read and parse time info
                //struct tm
                struct tm* tm_start = new struct tm;
                struct tm* tm_end   = new struct tm;
                stringstream stream_val_t_start( currLineSplitted[1] );
                stringstream stream_val_t_end(   currLineSplitted[2] );
                stream_val_t_start >> get_time(tm_start, "%Y-%m-%d %H:%M:%S");
                stream_val_t_end   >> get_time(tm_end,   "%Y-%m-%d %H:%M:%S");
                Global::set_ts_start_tm( tm_start );
                Global::set_ts_end_tm(   tm_end   );
                // read other values
                Global::set_tsteps_per_hour(stoi( currLineSplitted[3] ));
                Global::set_expansion_scenario_id(stoi( currLineSplitted[4] ));
                Global::set_exp_pv_kWp(        stof( currLineSplitted[5] ));
                Global::set_exp_bess_kW(       stof( currLineSplitted[6] ));
                Global::set_exp_bess_kWh(      stof( currLineSplitted[7] ));
                Global::set_exp_bess_start_soc(stof( currLineSplitted[8] ));
                Global::set_open_space_pv_kWp( stof( currLineSplitted[9] ));
                Global::set_wind_kWp(          stof( currLineSplitted[10]));
                break;
            } else {
                continue;
            }
        }
        scenarios_input.close();
    }
    return true;
}
*/


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
    string comparisonStr5 = "n_profiles_hp";
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
    } else if (comparisonStr5.compare(argv[0]) == 0) {
        Global::set_n_heatpump_profiles( stoul(argv[1]) );
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
    // convert time values
    struct tm* time_value = new struct tm;
    stringstream stream_time_value(argv[1]);
    stream_time_value >> get_time(time_value, "%Y-%m-%d %H:%M:%S");
    // add time values to global list
    global::time_timestep_id[pos] = current_time_index;
    global::time_localtime_str->push_back(time_value);
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
        Substation::InstantiateNewSubstation(current_station_id, stationName);
    } catch (runtime_error& e) {
        cerr << "Error when creating a substation:" << endl;
        cerr << e.what() << endl;
        return 1;
    }
    return 0;
}
int load_data_from_central_database_callbackD(void* data, int argc, char** argv, char** colName) {
    /* 
     * This is the callback function for geeting information about the control units
     * out of the database.
     * This function also creates the control units.
     * 
     * Columns:
	 * 0       1              2
	 * UnitID  substation_id  LocID
     * 
     */
    if (argc != 3) {
        cerr << "Number of arguments not equal to 3 for one row!" << endl;
        return 1;
    }
    int current_cu_id = stoi(argv[0]);
    int conn_to_subst_id = stoi(argv[1]);
    int location_id      = stoul(argv[2]);
    try {
        ControlUnit::InstantiateNewControlUnit(current_cu_id, conn_to_subst_id, location_id);
    } catch (runtime_error& e) {
        cerr << "Error when creating a control unit:" << endl;
        cerr << e.what() << endl;
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
	 * 0      1       2             3           4           5
	 * MeUID, UnitID, MeterPointID, has_demand, has_feedin, has_pv_residential
	 * 6                  7         8       9        10
	 * has_pv_open_space, has_bess, has_hp, has_chp, LocID
	 */
	if (argc != 11) {
		cerr << "Number of arguments not equal to 10 for one row!" << endl;
		return 1;
	}
	int current_mu_id = stoi(argv[0]);
	int conn_to_unitID = stoi(argv[1]);
	string* mPointIDStr= new string(argv[2]);
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
	data_input_path << Global::get_input_path() << "SeparatedSmartMeterData/";
	data_input_path << current_mu_id;
	data_input_path << ".csv";

	try {
        MeasurementUnit* newMU = MeasurementUnit::InstantiateNewMeasurementUnit(
                                current_mu_id, conn_to_unitID, mPointIDStr, locID,
								has_demand, has_feedin, has_pv_resid, has_pv_opens,
								has_bess,   has_hp,     has_wb,       has_chp);
		newMU->load_data(data_input_path.str().c_str());
    } catch (runtime_error& e) {
        cerr << "Error when creating a control unit:" << endl;
        cerr << e.what() << endl;
		return 1;
	}
	return 0;
}
int load_data_from_central_database_callback_PV_info(void* data, int argc, char** argv, char** colName) {
    /*
     * This is the callback function for geeting information about available
     * orientations and their counts of the global PV profiles.
     *
     * The first argument (data) holds the reference to the target array, where
     * data should be written into.
     *
     * Columns:
     * 0            1
     * orientation  number_of_ts
     */
    if (argc != 2) {
        cerr << "Number of arguments not equal to 2 for one row!" << endl;
        return 1;
    }

    string orientation  = argv[0];
    size_t number_of_ts = stoul(argv[1]);
    global::pv_profiles_information[orientation] = number_of_ts;

    return 0;
}
int load_data_from_central_database_callback_PV(void* data, int argc, char** argv, char** colName) {
    /*
     * This is the callback function for geeting the global PV profiles.
     *
     * The first argument (data) holds the reference to the target array, where
     * data should be written into.
     *
     * Columns:
     * 0           1             2            3
     * TimestepID  Value_Feedin  Orientation  SameOrientationTimeSeriesIndex
     */
    static int callcounter_timestepID = 1; // internal counter for the timestep ID
    static unsigned long callcounter_timeseries = 0; // internal counter for the timeseries
    static string last_orientation;
    static int last_same_orientation_ts_index = -1;

    // check callcounter
    if (callcounter_timeseries >= Global::get_n_pv_profiles()) {
        cerr << "Error in time series table for PV profiles.\n";
        cerr << "It is said that there are " << Global::get_n_pv_profiles() << " PV profiles, but therer are more available in the table." << endl;
        return 1;
    }

    if (argc != 4) {
        cerr << "Number of arguments not equal to 4 for one row!" << endl;
        return 1;
    }

    // parse the current values
    int current_timestepID;
    int current_same_orientation_ts_index;
    string current_orientation;
    try {
        current_timestepID = stoi(argv[0]);
        current_same_orientation_ts_index = stoi(argv[3]);
        current_orientation = argv[2];
    } catch (const exception& ex) {
        cerr << "An error occured during the reading of the PV profiles:" << endl;
        cerr << ex.what() << endl;
        return 1;
    }

    // initialize last_orientation and last_same_orientation_tsIndex
    // at the beginning of a new time series
    if (callcounter_timestepID == 1) {
        last_same_orientation_ts_index = current_same_orientation_ts_index;
        last_orientation = current_orientation;
    } else {
        // else check the values
        if (last_orientation != current_orientation || last_same_orientation_ts_index != current_same_orientation_ts_index) {
            cerr << "There is at least one global PV profile where values are missing!" << endl;
            return 1;
        }
    }

    if (current_timestepID != callcounter_timestepID) {
        cerr << "Wrong ordering of the global PV profile values!" << endl;
        return 1;
    }

    int pos = callcounter_timestepID - 1; // the current position is one behind the callcounter
    ((float**) data)[callcounter_timeseries][pos] = stof(argv[1]);

    if (callcounter_timestepID < Global::get_n_timesteps()) {
        callcounter_timestepID++;
    } else {
        // time series is loaded completly
        // 1. add this time series to the global list
            /*
            This is done automatically
        auto search_result = global::pv_profiles_per_ori.find(current_orientation);
        if (search_result == global::pv_profiles_per_ori.end())
            global::pv_profiles_per_ori[current_orientation] = vector<const float *>();
            */
        global::pv_profiles_per_ori[current_orientation].push_back( ((float**) data)[callcounter_timeseries] );
        // 2. reset counters
        callcounter_timestepID = 1;
        callcounter_timeseries++;
    }
    return 0;
}
int load_data_from_central_database_callback_Wind(void* data, int argc, char** argv, char** colName) {
    /*
     * This is the callback function for geeting the global Wind profile.
     *
     * The first argument (data) holds the reference to the target array, where
     * data should be written into.
     *
     * Columns:
     * 0           1
     * TimestepID  wind_profile_value
     */
    static int callcounter = 1;
    int pos = callcounter - 1; // the current position is one behind the callcounter
    if (argc != 2) {
        cerr << "Number of arguments not equal to 2 for one row!" << endl;
        return 1;
    }

    if (stoi(argv[0]) != callcounter) {
        cerr << "Wrong ordering of the global PV profile values!" << endl;
        return 1;
    }
    ((float*) data)[pos] = stof(argv[1]);

    callcounter++;
    return 0;
}
int load_data_from_central_database_callback_HP(void* data, int argc, char** argv, char** colName) {
    /*
     * This is the callback function for geeting the global heat pump profiles.
     *
     * The first argument (data) holds the reference to the target array, where
     * data should be written into.
     *
     * Columns:
     * 0           1             2
     * TimestepID  Value_Demand  TimeSeriesIndex
     */
    static int callcounter_timestepID = 1; // internal counter for the timestep ID
    static unsigned long callcounter_timeseries = 0; // internal counter for the timeseries, both are used to identify missing values

    // check callcounter
    if (callcounter_timeseries >= Global::get_n_heatpump_profiles()) {
        cerr << "Error in data for heat pumps.\n";
        cerr << "It is said that there are " << Global::get_n_heatpump_profiles() << " heat pump profiles, but therer are more available in the table." << endl;
        return 1;
    }

    if (argc != 3) {
        cerr << "Number of arguments not equal to 3 for one row!" << endl;
        return 1;
    }

    if (stoi(argv[0]) != callcounter_timestepID) {
        cerr << "There is one row missing for at least one timestep in the list of heat pump profiles!" << endl;
        return 1;
    }
    if (stoul(argv[2]) != callcounter_timeseries) {
        cerr << "There are missing values for at least one time series in the list of heat pump profiles!" << endl;
        return 1;
    }
    int pos = callcounter_timestepID - 1; // the current position is one behind the callcounter
    ((float**) data)[callcounter_timeseries][pos] = stof(argv[1]);

    if (callcounter_timestepID < Global::get_n_timesteps()) {
        callcounter_timestepID++;
    } else {
        callcounter_timestepID = 1;
        callcounter_timeseries++;
    }
    return 0;
}
int load_data_from_central_database_callback_address_data_A(void* data, int argc, char** argv, char** colName) {
    /*
     * This is the callback function for geeting the yearly heat pump electricity demand in kWh per Location ID.
     *
     * Columns:
     * 0      1
     * LocID  YearlyHPHeatDemand_kWh
     */
    if (argc != 2) {
        cerr << "Number of arguments not equal to 2 for one row!" << endl;
        return 1;
    }
    global::yearly_hp_energy_demand_kWh[ stoul(argv[0]) ] = stof(argv[1]);
    return 0;
}
int load_data_from_central_database_callback_address_data_B(void* data, int argc, char** argv, char** colName) {
    /*
     * This is the callback function for geeting the roof orientations per location.
     *
     * Columns:
     * 0      1           2
     * LocID  Area_in_m2  Orientation
     */
    if (argc != 3) {
        cerr << "Number of arguments not equal to 3 for one row!" << endl;
        return 1;
    }
    //
    // Assume a flat roof to be a south orientation where modules are installed with a 45-degree pitch
    //     in this case, we assume a 45 degree pitch of the installed modules, so the acual area is 
    //     acual_usabel_area = given_area / cos(45 deg) = given_area / 0.707107
    float usable_area = stof(argv[1]);
    string orientation ( argv[2] );
    if (orientation == "flat_roof" || orientation == "f") {
        // compute new usable area as stated above
        usable_area = usable_area / 0.707107;
        orientation = "S";
    }
    global::roof_section_orientations[ stoul(argv[0]) ].push_back( pair<float, std::string>(usable_area, orientation) );
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
            cerr << "Error when reading the SQL-Table: " << sqlErrorMsgA << endl;
            sqlite3_free(sqlErrorMsgA);
            return false;
        }

        //
        // Initialize global time list
        //
        global::time_timestep_id = new int[Global::get_n_timesteps()];
        global::time_localtime_str = new vector<struct tm*>();
        global::time_localtimezone_str = new vector<string>();
        //
        // Load time indices
        //
        string sql_queryB = "SELECT TimestepID, local_time, local_time_zone FROM time_indices ORDER BY TimestepID;";
        char* sqlErrorMsgB;
        int ret_valB = sqlite3_exec(dbcon, sql_queryB.c_str(), load_data_from_central_database_callbackB, NULL, &sqlErrorMsgB);
        if (ret_valB != 0) {
            cerr << "Error when reading the SQL-Table: " << sqlErrorMsgB << endl;
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
        string sql_queryC = "SELECT substation_id, substation_name FROM list_of_substations ORDER BY substation_id;";
        char* sqlErrorMsgC;
        int ret_valC = sqlite3_exec(dbcon, sql_queryC.c_str(), load_data_from_central_database_callbackC, NULL, &sqlErrorMsgC);
        if (ret_valC != 0) {
            cerr << "Error when reading the SQL-Table: " << sqlErrorMsgC << endl;
            sqlite3_free(sqlErrorMsgC);
            return false;
        }
        // 2. CUs
        string sql_queryD = "SELECT UnitID, substation_id, LocID FROM list_of_control_units ORDER BY UnitID;";
        char* sqlErrorMsgD;
        int ret_valD = sqlite3_exec(dbcon, sql_queryD.c_str(), load_data_from_central_database_callbackD, NULL, &sqlErrorMsgD);
        if (ret_valD != 0) {
            cerr << "Error when reading the SQL-Table: " << sqlErrorMsgD << endl;
            sqlite3_free(sqlErrorMsgD);
            return false;
        }
        // 3. MUs
        string sql_queryE = "SELECT MeUID, UnitID, MeterPointID, has_demand, has_feedin, has_pv_residential, has_pv_open_space, has_bess, has_hp, has_chp, LocID FROM list_of_measurement_units ORDER BY MeUID;";
        char* sqlErrorMsgE;
        int ret_valE = sqlite3_exec(dbcon, sql_queryE.c_str(), load_data_from_central_database_callbackE, NULL, &sqlErrorMsgE);
        if (ret_valE != 0) {
            cerr << "Error when reading the SQL-Table: " << sqlErrorMsgE << endl;
            sqlite3_free(sqlErrorMsgE);
            return false;
        }

        //
        // Load central solar radation profiles
        //
        // 1. Load metadata
        string sql_query = "SELECT orientation, number_of_ts FROM global_profiles_pv_info;";
        char* sqlErrorMsgF;
        int ret_valF = sqlite3_exec(dbcon, sql_query.c_str(), load_data_from_central_database_callback_PV_info, NULL, &sqlErrorMsgF);
        if (ret_valF != 0) {
            cerr << "Error when reading the SQL-Table: " << sqlErrorMsgF << endl;
            sqlite3_free(sqlErrorMsgF);
            return false;
        }
        unsigned long number_of_pv_profiles = accumulate(
                begin(global::pv_profiles_information),
                end(global::pv_profiles_information),
                0,
                [](const size_t prev, const auto& elem){ return prev + elem.second; });
        Global::set_n_pv_profiles(number_of_pv_profiles);
        // 2. Load the concrete profiles
        float** new_pv_array = new float*[Global::get_n_pv_profiles()];
        for (size_t pv_idx = 0; pv_idx < Global::get_n_pv_profiles(); pv_idx++) {
            new_pv_array[pv_idx] = new float[Global::get_n_timesteps()];
        }
        sql_query = "SELECT TimestepID,Value_Feedin,Orientation,SameOrientationTimeSeriesIndex FROM global_profiles_pv ORDER BY Orientation,SameOrientationTimeSeriesIndex,TimestepID;";
        ret_valF = sqlite3_exec(dbcon, sql_query.c_str(), load_data_from_central_database_callback_PV, new_pv_array/*Reference to the new array*/, &sqlErrorMsgF);
        if (ret_valF != 0) {
            cerr << "Error when reading the SQL-Table: " << sqlErrorMsgF << endl;
            sqlite3_free(sqlErrorMsgF);
            return false;
        }
        global::pv_profiles_data = new_pv_array;
        //
        // initialize the global open space pv unit
        global::unit_open_space_pv = new OpenSpacePVOrWind(Global::get_open_space_pv_kWp(), OpenSpacePVOrWindType::PV);

        //
        // Load central wind profile
        //
        float* new_wind_array = new float[Global::get_n_timesteps()];
        sql_query = "SELECT TimestepID,wind_profile_value FROM global_profile_wind;";
        ret_valF  = sqlite3_exec(dbcon, sql_query.c_str(), load_data_from_central_database_callback_Wind, new_wind_array/*Reference to the new array*/, &sqlErrorMsgF);
        if (ret_valF != 0) {
            cerr << "Error when reading the SQL-Table: " << sqlErrorMsgF << endl;
            sqlite3_free(sqlErrorMsgF);
            return false;
        }
        global::wind_profile = new_wind_array;
        //
        // initialize the global wind turbine unit
        global::unit_open_space_wind = new OpenSpacePVOrWind(Global::get_wind_kWp(), OpenSpacePVOrWindType::Wind);

        //
        // Load central heat pump profiles
        //
        float** new_hp_profile_array = new float*[Global::get_n_heatpump_profiles()];
        for (unsigned long hp_idx = 0; hp_idx < Global::get_n_heatpump_profiles(); hp_idx++) {
            new_hp_profile_array[hp_idx] = new float[Global::get_n_timesteps()];
        }
        sql_query = "SELECT TimestepID,Value_Demand,TimeSeriesIndex FROM global_profiles_heatpumps ORDER BY TimeSeriesIndex,TimestepID;";
        ret_valF  = sqlite3_exec(dbcon, sql_query.c_str(), load_data_from_central_database_callback_HP, new_hp_profile_array/*Reference to the new array*/, &sqlErrorMsgF);
        if (ret_valF != 0) {
            cerr << "Error when reading the SQL-Table: " << sqlErrorMsgF << endl;
            sqlite3_free(sqlErrorMsgF);
            return false;
        }
        global::hp_profiles = new_hp_profile_array;

        //
        // Load address data
        //
        // 1. yearly heat demand for heat pumps
        sql_query = "SELECT LocID, YearlyHPElectricityDemand_kWh FROM address_data ORDER BY LocID;";
        ret_valF = sqlite3_exec(dbcon, sql_query.c_str(), load_data_from_central_database_callback_address_data_A, NULL, &sqlErrorMsgF);
        if (ret_valF != 0) {
            cerr << "Error when reading the SQL-Table: " << sqlErrorMsgF << endl;
            sqlite3_free(sqlErrorMsgF);
            return false;
        }
        // 2. get information about roof sections
        sql_query = "SELECT LocID, Area_in_m2, Orientation FROM address_roof_data ORDER BY LocID;";
        ret_valF = sqlite3_exec(dbcon, sql_query.c_str(), load_data_from_central_database_callback_address_data_B, NULL, &sqlErrorMsgF);
        if (ret_valF != 0) {
            cerr << "Error when reading the SQL-Table: " << sqlErrorMsgF << endl;
            sqlite3_free(sqlErrorMsgF);
            return false;
        }

        sqlite3_close(dbcon);
        return true;
    } else {
        cerr << "Error when connecting to the database" << endl;
        return false;
    }
}





int expansion::expCombiMatrixOrderToBitRepr(int indexMatO) {
    /*
     * This function maps an expansion combination number (as orderd in the
     * expansion matrix) to its bitwise representation.
     */
    switch (indexMatO) {
        case  0: return MaskNothing;
        case  1: return MaskPV;
        case  2: return        MaskBS;
        case  3: return               MaskHP;
        case  4: return                      MaskWB;
        case  5: return MaskPV|MaskBS;
        case  6: return MaskPV|       MaskHP;
        case  7: return MaskPV|              MaskWB;
        case  8: return        MaskBS|MaskHP;
        case  9: return        MaskBS|       MaskWB;
        case 10: return               MaskHP|MaskWB;
        case 11: return MaskPV|MaskBS|MaskHP;
        case 12: return MaskPV|MaskBS|       MaskWB;
        case 13: return MaskPV|       MaskHP|MaskWB;
        case 14: return        MaskBS|MaskHP|MaskWB;
        case 15: return MaskPV|MaskBS|MaskHP|MaskWB;
    }
    throw logic_error("Impossible index passed to function!");
    return 0;
}

int expansion::expCombiBitReprToMatrixOrder(int bitRepr) {
    /*
     * This function maps an bitwise representation of the expansion 
     * combination to the number (as orderd in the expansion matrix).
     * It is the inverse to expCombiMatrixOrderToBitRepr()
     */
    if      (bitRepr ==  MaskNothing)
        return  0;
    else if (bitRepr ==  MaskPV)
        return  1;
    else if (bitRepr ==         MaskBS)
        return  2;
    else if (bitRepr ==                MaskHP)
        return  3;
    else if (bitRepr ==                       MaskWB)
        return  4;
    else if (bitRepr == (MaskPV|MaskBS)               )
        return  5;
    else if (bitRepr == (MaskPV|       MaskHP)        )
        return  6;
    else if (bitRepr == (MaskPV|              MaskWB) )
        return  7;
    else if (bitRepr == (       MaskBS|MaskHP)        )
        return  8;
    else if (bitRepr == (       MaskBS|       MaskWB) )
        return  9;
    else if (bitRepr == (              MaskHP|MaskWB) )
        return 10;
    else if (bitRepr == (MaskPV|MaskBS|MaskHP)        )
        return 11;
    else if (bitRepr == (MaskPV|MaskBS|       MaskWB) )
        return 12;
    else if (bitRepr == (MaskPV|       MaskHP|MaskWB) )
        return 13;
    else if (bitRepr == (       MaskBS|MaskHP|MaskWB) )
        return 14;
    else if (bitRepr == (MaskPV|MaskBS|MaskHP|MaskWB) )
        return 15;

    throw logic_error("Invalid bit representation!");
}

int expansion::genExpCombiAsBitRepr(bool has_pv, bool has_bs, bool has_hp, bool has_wb) {
	/*
	 * This function returns the binary representation of the expansion
	 */
	int retval = 0;
	if (has_pv)
		retval = retval | MaskPV;
	if (has_bs)
		retval = retval | MaskBS;
	if (has_hp)
		retval = retval | MaskHP;
	if (has_wb)
		retval = retval | MaskWB;
	return retval;
}

bool expansion::isExpCombiPossible(int currExpNumber, int newExpNumber) {
	/*
	This function returns if a expansion from currExpNumber to newExpNumber (i.e. the future scenario number) is possible or not
	*/
	if (currExpNumber == newExpNumber)
		return false; // this is not a expansion by definition!
	int currExpBitRepr = expCombiBitReprToMatrixOrder(currExpNumber);
	int newExpBitRepr  = expCombiBitReprToMatrixOrder(newExpNumber);
	// A combination is possible, if and only if
	//  - the bitwise and of curr and new state returns in the curr state
	// and
	//  - the bitwise or of curr and new state returns in the new state.
	// One can see, that it is enough to check one only.
	return (currExpBitRepr & newExpBitRepr) == currExpBitRepr;
}


/*
Loads the expansion matrix into the first argument.

Returns false if an error occurs, else true.
*/
bool expansion::load_expansion_matrix(float expansion_matrix[16][16]) {
	ifstream expmat_input;
	stringstream expmat_input_path;
	expmat_input_path << "./expansion_scenarios/";
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
			if (expansion::isExpCombiPossible(i,j)) {
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

void expansion::add_expansion_to_units(float expansion_matrix_rel_freq[16][16], int expansion_matrix_abs_freq[16][16], int scenario_id) {
	/*
	 * This function adds the expansion given als relative counts in expansion_matrix_rel_freq
	 * to the control units.
	 * @param expansion_matrix_abs_freq contains the absolute counts afterwards
	 */

	int currExpCountsBitIndexed[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // order as given from bitwise representation (as this represents a sequential integer as well)
	int currExpCountsMatIndexed[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // order as in expansion matrix
	int newExpCountsMatIndexed[16]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	vector<list<ControlUnit*>> cuRefLstVectBitOrder(16); // vector of 16 lists, that contains the references to the CUs, where the index of the list in the vector corresponds to the expansion in bitwise order

	//
	// 1. count current expansion status (as it is in the given data)
	//    and store the references to the CUs for a given current expansion in a list
	ControlUnit *const * unit_list = ControlUnit::GetArrayOfInstances();
	const int n_CUs = ControlUnit::GetNumberOfInstances();
	for (int i = 0; i < n_CUs; i++) {
		ControlUnit* current_unit = unit_list[i];
		int expCombi = current_unit->get_exp_combi_bit_repr();
		currExpCountsBitIndexed[ expCombi ]++;
		cuRefLstVectBitOrder[ expCombi ].push_back( current_unit ); // add unit to reference list
	}
	// change order
	for (int i = 0; i < 16; i++) {
		currExpCountsMatIndexed[i] = currExpCountsBitIndexed[ expCombiMatrixOrderToBitRepr(i) ];
	}

	//
	// 2. calculate expansion matrix with absolute values based
	//    on the abolute counts from 1 by a row-wise multiplication
	//    Note: Impossible values will be jumped, and diagonal values
	//          will be computed afterwards as a difference because
	//          we the absolute numbers have to be integer (half units
	//          are impossible, obviously)
	for (int i = 0; i < 16; i++) {
		for (int j = i+1; j < 16; j++) {
			expansion_matrix_abs_freq[i][j] = expansion_matrix_rel_freq[i][j] * currExpCountsMatIndexed[i];
		}
	}
	// calculate diagonal values
	for (int i = 0; i < 16; i++) {
		int sum_expanded_i = 0;
		for (int j = i+1; j < 16; j++){
			sum_expanded_i += expansion_matrix_abs_freq[i][j];
		}
		expansion_matrix_abs_freq[i][i] = currExpCountsMatIndexed[i] - sum_expanded_i;
	}

	//
	// 3. count new expansion status, that will be simulated now
	//    This is the column sum of the absolute expansion matrix
	for (int j = 0; j < 16; j++) { // sum over j (cols) instead of i (rows) first
		for (int i = 0; i < 16; i++) {
			newExpCountsMatIndexed[j] += expansion_matrix_abs_freq[i][j];
		}
	}

	//
	// 4. plan and execute expansion
	for (int iMatO = 0; iMatO < 16; iMatO++) {
		int iBitO = expCombiMatrixOrderToBitRepr( iMatO ); // get index in Bitwise Order (BitO)
		list<ControlUnit*>* listOfCUs = &(cuRefLstVectBitOrder[ iBitO ]);
		list<ControlUnit*>::iterator iter = listOfCUs->begin();
		// loop over all current expansion states
		for (int jExpTargetMatO = 0; jExpTargetMatO < 16; jExpTargetMatO++) {
			// get number of CUs that get the current expansion
			int numThisCombi_i_j = expansion_matrix_abs_freq[iMatO][jExpTargetMatO];
			// find out, which units we have to add for this i/j-combination
			int iBitRepr = expCombiMatrixOrderToBitRepr(iMatO);
			int jBitRepr = expCombiMatrixOrderToBitRepr(jExpTargetMatO);
			int ijXOR = iBitRepr ^ jBitRepr;
			bool expPV = false;
			bool expBS = false;
			bool expHP = false;
			bool expWB = false;
			if (ijXOR & MaskPV) expPV = true;
			if (ijXOR & MaskBS) expBS = true;
			if (ijXOR & MaskHP) expHP = true;
			if (ijXOR & MaskWB) expWB = true;
			// loop over this number
			for (int n = 0; n < numThisCombi_i_j; n++) {
				if (iter == listOfCUs->end()) {
					cerr << "Warning: end of list for expansion reached before all expansion planing were fulfilled." << endl;
					goto outer_loop_end;
				}
				// 1. add components
				if (expPV) (*iter)->add_exp_pv();
				if (expBS) (*iter)->add_exp_bs();
				if (expHP) (*iter)->add_exp_hp();
				if (expWB) (*iter)->add_exp_wb();
				// 2. remove from list (would be good, but not required)
				iter++;
			}
		}
		outer_loop_end:;
	}

	//
	// finally: write expansion information to file
	// A. output expansion matrix with absolute numbers
    filesystem::path info_path_A (*(global::current_global_output_dir)); // for the expansion matrix it is acceptable to use the static output path (that does not change over time for different parameter variations, as expansion cannot change)
	info_path_A /= "expansion-matrix-abs-values.csv";
	ofstream output_exp_mat(info_path_A, std::ofstream::out);
	output_exp_mat << ",0. Nothing,1. PV,2. BS,3. HP,4. WB,5. PV+BS,6. PV+HP,7. PV+WB,8. BS+HP,9. BS+WB,10. HP+WB,11. PV+BS+HP,12. PV+BS+WB,13. PV+HP+WB,14. BS+HP+WB,15. PV+BS+HP+WB,Sum as in data" << endl;
	const char * first_column[16] = {"0. Nothing","1. PV","2. BS","3. HP","4. WB","5. PV+BS","6. PV+HP","7. PV+WB","8. BS+HP","9. BS+WB","10. HP+WB","11. PV+BS+HP","12. PV+BS+WB","13. PV+HP+WB","14. BS+HP+WB","15. PV+BS+HP+WB"};
	for (int i = 0; i < 16; i++) {
		output_exp_mat << first_column[i];
		for (int j = 0; j < 16; j++) {
			output_exp_mat << "," << expansion_matrix_abs_freq[i][j];
		}
		output_exp_mat << "," << currExpCountsMatIndexed[i] << endl;
	}
	output_exp_mat << "Sum as simulated";
	for (int i = 0; i < 16; i++)
		output_exp_mat << "," << newExpCountsMatIndexed[i];
	output_exp_mat << "," << endl;
	output_exp_mat.close();
	//
	// B. output information about added components per MeUID
    filesystem::path info_path_B (*(global::current_global_output_dir)); // same argument as 21 lines above
    info_path_B /= "expansion-per-cu.csv";
	ofstream output_per_cu(info_path_B, std::ofstream::out);
	output_per_cu << "UnitID,n_MUs,pv_orig,pv_added,bs_orig,bs_added,hp_orig,hp_added,wb_orig,wb_added,added_pv_kWp,added_bess_E_kWh,added_bess_P_kW" << endl;
	// n_CUs and unit_list defined above, at 1.
	for (int i = 0; i < n_CUs; i++) {
		ControlUnit* current_unit = unit_list[i];
		int expCombiAsInData    = current_unit->get_exp_combi_bit_repr_from_MUs();
		int expCombiAsSimulated = current_unit->get_exp_combi_bit_repr_sim_added();
		// output information
		output_per_cu <<        current_unit->get_unitID();
		output_per_cu << "," << current_unit->get_n_MUs();
		output_per_cu << "," << (0 < (expansion::MaskPV & expCombiAsInData));
		output_per_cu << "," << (0 < (expansion::MaskPV & expCombiAsSimulated));
		output_per_cu << "," << (0 < (expansion::MaskBS & expCombiAsInData));
		output_per_cu << "," << (0 < (expansion::MaskBS & expCombiAsSimulated));
		output_per_cu << "," << (0 < (expansion::MaskHP & expCombiAsInData));
		output_per_cu << "," << (0 < (expansion::MaskHP & expCombiAsSimulated));
		output_per_cu << "," << (0 < (expansion::MaskWB & expCombiAsInData));
		output_per_cu << "," << (0 < (expansion::MaskWB & expCombiAsSimulated));
		output_per_cu << "," << endl;
	}
	output_per_cu.close();
}
