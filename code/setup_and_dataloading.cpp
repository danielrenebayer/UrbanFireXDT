#include "setup_and_dataloading.h"


#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
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
#include "sac_planning.h"
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
        float  exp_bs_E_P   = 0.0; bool exp_bs_E_P_set   = false;
        float  os_pv_kWp    = 0.0; bool os_pv_kWp_set    = false;
        float  os_wind_kWp  = 0.0; bool os_wind_kWp_set  = false;
        string exp_profile_mode    = ""; bool exp_profile_mode_set    = false;
        string sac_planning_mode   = ""; bool sac_planning_mode_set   = false;
        string exp_bs_P_comp_mode  = ""; bool exp_bs_P_comp_mode_set  = false;
        // variables for PV expansion
        float  exp_pv_kWp_static  = 0.0; bool exp_pv_kWp_static_set   = false;
        float  exp_pv_kWp_m2_roof = 0.0; bool exp_pv_kWp_m2_roof_set  = false;
        float  exp_pv_min_kWp     = 0.0; bool exp_pv_min_kWp_set      = false;
        float  exp_pv_max_kWp     = 0.0; bool exp_pv_max_kWp_set      = false;
        bool   exp_pv_mode_static = false; bool exp_pv_mode_static_set= false;
        float  exp_pv_max_total_kWp = 0.0; bool exp_pv_max_total_kWp_set = false;
        float  exp_pv_max_unit_kWp  = 0.0; bool exp_pv_max_unit_kWp_set  = false;
        // variables of NPV calculation
        float  feed_in_tariff     = 0.0; bool feed_in_tariff_set      = false;
        float  demand_tariff      = 0.0; bool demand_tariff_set       = false;
        float  npv_discount_rate  = 0.0; bool npv_discount_rate_set   = false;
        uint   npv_time_horizon   = 0.0; bool npv_time_horizon_set    = false;
        float  icost_PV_per_kWp   = 0.0; bool icost_PV_per_kWp_set    = false;
        float  icost_BS_per_kWh   = 0.0; bool icost_BS_per_kWh_set    = false;
        // variables, that always default to a value
        bool   use_BS_for_SSR_list = false;

        //
        // define internal functions (here i.e. a lambda function with complete capture-by-reference)
        auto parse_element = [&](string& element_name, boost::property_tree::ptree& scenario_dict) -> void {
            if ( element_name.compare("data input path") == 0 ) {
                str_data_ipt    = scenario_dict.get_value<string>();
                str_data_ipt_set= true;
            } else if ( element_name.compare("data output path") == 0 ) {
                str_data_opt    = scenario_dict.get_value<string>();
                str_data_opt_set= true;
            } else if ( element_name.compare("start") == 0 ) {
                start_str       = scenario_dict.get_value<string>();
                start_str_set   = true;
            } else if ( element_name.compare("end") == 0 ) {
                end_str         = scenario_dict.get_value<string>();
                end_str_set     = true;
            } else if ( element_name.compare("time steps per hour") == 0 ) {
                ts_per_hour     = scenario_dict.get_value<int>();
                ts_per_hour_set = true;
            } else if ( element_name.compare("expansion id") == 0 ) {
                expansionID     = scenario_dict.get_value<int>();
                expansionID_set = true;
            } else if ( element_name.compare("expansion BS P in kW") == 0 ) {
                exp_bs_kW       = scenario_dict.get_value<float>();
                exp_bs_kW_set   = true;
            } else if ( element_name.compare("expansion BS E in kWh") == 0 ) {
                exp_bs_kWh      = scenario_dict.get_value<float>();
                exp_bs_kWh_set  = true;
            } else if ( element_name.compare("expansion BS initial SOC") == 0 ) {
                exp_bs_iSOC     = scenario_dict.get_value<float>();
                exp_bs_iSOC_set = true;
            } else if (element_name.compare("expansion BS E:P ratio") == 0) {
                exp_bs_E_P      = scenario_dict.get_value<float>();
                exp_bs_E_P_set  = true;
            } else if ( element_name.compare("open space PV kWp") == 0 ) {
                os_pv_kWp       = scenario_dict.get_value<float>();
                os_pv_kWp_set   = true;
            } else if ( element_name.compare("open space wind kWp") == 0 ) {
                os_wind_kWp     = scenario_dict.get_value<float>();
                os_wind_kWp_set = true;
            } else if ( element_name.compare("expansion profile selection") == 0 ) {
                exp_profile_mode       = scenario_dict.get_value<string>();
                exp_profile_mode_set   = true;
            } else if ( element_name.compare("CU selection mode for comp. add.") == 0 ) {
                sac_planning_mode      = scenario_dict.get_value<string>();
                sac_planning_mode_set  = true;
            } else if ( element_name.compare("expansion PV kWp static") == 0 ) {
                exp_pv_kWp_static      = scenario_dict.get_value<float>();
                exp_pv_kWp_static_set  = true;
            } else if ( element_name.compare("expansion PV min kWp for section usage") == 0 ) {
                exp_pv_min_kWp         = scenario_dict.get_value<float>();
                exp_pv_min_kWp_set     = true;
            } else if ( element_name.compare("expansion PV max inst kWp per section") == 0 ) {
                exp_pv_max_kWp         = scenario_dict.get_value<float>();
                exp_pv_max_kWp_set     = true;
            } else if ( element_name.compare("expansion PV kWp per roof area in m2") == 0 ) {
                exp_pv_kWp_m2_roof     = scenario_dict.get_value<float>();
                exp_pv_kWp_m2_roof_set = true;
            } else if ( element_name.compare("expansion PV kWp static mode") == 0 ) {
                exp_pv_mode_static     = scenario_dict.get_value<bool>();
                exp_pv_mode_static_set = true;
            } else if ( element_name.compare("tariff feed-in per kWh") == 0 ) {
                feed_in_tariff         = scenario_dict.get_value<float>();
                feed_in_tariff_set     = true;
            } else if ( element_name.compare("tariff demand per kWh") == 0 ) {
                demand_tariff          = scenario_dict.get_value<float>();
                demand_tariff_set      = true;
            } else if ( element_name.compare("net present value discount rate") == 0 ) {
                npv_discount_rate      = scenario_dict.get_value<float>();
                npv_discount_rate_set  = true;
            } else if ( element_name.compare("net present value time horizon in years") == 0 ) {
                npv_time_horizon       = scenario_dict.get_value<unsigned int>();
                npv_time_horizon_set   = true;
            } else if ( element_name.compare("installation cost PV per kWp") == 0 ) {
                icost_PV_per_kWp       = scenario_dict.get_value<float>();
                icost_PV_per_kWp_set   = true;
            } else if ( element_name.compare("installation cost BS per kWh") == 0 ) {
                icost_BS_per_kWh       = scenario_dict.get_value<float>();
                icost_BS_per_kWh_set   = true;
            } else if ( element_name.compare("expansion BS power computation mode") == 0 ) {
                exp_bs_P_comp_mode       = scenario_dict.get_value<string>();
                exp_bs_P_comp_mode_set   = true;
            } else if ( element_name.compare("use BS for creating SSR list for SAC planning") == 0 ) {
                use_BS_for_SSR_list     = scenario_dict.get_value<bool>();
            } else if ( element_name.compare("expansion PV max inst kWp per unit") == 0 ) {
                exp_pv_max_unit_kWp     = scenario_dict.get_value<float>();
                exp_pv_max_unit_kWp_set = true;
            } else if ( element_name.compare("expansion PV max total kWp addition") == 0 ) {
                exp_pv_max_total_kWp    = scenario_dict.get_value<float>();
                exp_pv_max_total_kWp_set= true;
            } else if ( element_name.compare("id") == 0 ) {
            } else if ( element_name.compare("comment") == 0 ) {
            } else {
                cout << "Unknonw config parameter " << element_name << endl;
            }
            return;
        };

        //
        // read default values
        for (auto& scenario_dict_all : tree_root.get_child("Default Scenario Values")) {
            string element_name = scenario_dict_all.first;
            parse_element(element_name, scenario_dict_all.second);
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
                    parse_element(element_name, s.second);
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
        // transform parameters if required (e.g. for exp_profile_mode)
        //   a) expansion profile allocation mode
        global::ExpansionProfileAllocationMode exp_profile_mode_transf = global::ExpansionProfileAllocationMode::Uninitialized;
        if (exp_profile_mode_set) {
            if      (exp_profile_mode == "as in data")
                exp_profile_mode_transf = global::ExpansionProfileAllocationMode::AsInData;
            else if (exp_profile_mode == "random")
                exp_profile_mode_transf = global::ExpansionProfileAllocationMode::Random;
            else {
                cerr << "Parameter 'expansion profile selection' is defined as '" << exp_profile_mode << "' in config-json, but this value is unknown." << endl;
                throw runtime_error("Parameter 'expansion profile selection' as defined in config-json is unknown.");
            }
        }
        //   b) Control Unit Selectio Mode For Component Addition
        global::CUSModeFCA sac_planning_mode_transl = global::CUSModeFCA::Uninitialized;
        if (sac_planning_mode_set) {
            if (sac_planning_mode == "as in data")  {
                sac_planning_mode_transl = global::CUSModeFCA::OrderAsInData;
            } else if (sac_planning_mode == "random")  {
                sac_planning_mode_transl = global::CUSModeFCA::RandomSelection;
            } else if (sac_planning_mode == "best SSR")  {
                sac_planning_mode_transl = global::CUSModeFCA::BestSSR;
            } else if (sac_planning_mode == "best NPV")  {
                sac_planning_mode_transl = global::CUSModeFCA::BestNPV;
            } else {
                cerr << "Parameter 'CU selection mode for comp. add' is defined as '" << sac_planning_mode << "' in config-json, but this value is unknown." << endl;
                throw runtime_error("Parameter 'CU selection mode for comp. add' as defined in config-json is unknown.");
            }
        }
        //   c) Battery Power Computation Mode
        global::BatteryPowerComputationMode exp_bs_P_comp_mode_transl = global::BatteryPowerComputationMode::AsDefinedByConfigVar;
        if (exp_bs_P_comp_mode_set) {
            if (exp_bs_P_comp_mode == "Power as given") {
                exp_bs_P_comp_mode_transl = global::BatteryPowerComputationMode::AsDefinedByConfigVar;
            } else if (exp_bs_P_comp_mode == "Use E:P-ratio") {
                exp_bs_P_comp_mode_transl = global::BatteryPowerComputationMode::UseEOverPRatio;
            } else {
                cerr << "Parameter 'expansion BS power computation mode' is defined as '" << exp_bs_P_comp_mode << "' in config-json, but this value is unknown." << endl;
                throw runtime_error("Parameter 'expansion BS power computation mode' as defined in config-json is unknown.");
            }
        }

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
            if (exp_bs_E_P_set)   Global::set_exp_bess_E_P_ratio(exp_bs_E_P);
            if (os_pv_kWp_set)    Global::set_open_space_pv_kWp(os_pv_kWp);
            if (os_wind_kWp_set)  Global::set_wind_kWp(os_wind_kWp);
            if (exp_profile_mode_set)   Global::set_exp_profile_mode(exp_profile_mode_transf);
            if (sac_planning_mode_set)  Global::set_cu_selection_mode_fca(sac_planning_mode_transl);
            if (exp_pv_mode_static_set) Global::set_exp_pv_mode(exp_pv_mode_static);
            if (exp_pv_kWp_static_set)  Global::set_exp_pv_kWp_static(exp_pv_kWp_static);
            if (exp_pv_kWp_m2_roof_set) Global::set_exp_pv_kWp_per_m2(exp_pv_kWp_m2_roof);
            if (exp_pv_max_unit_kWp_set)  Global::set_exp_pv_max_kWp_per_unit(exp_pv_max_unit_kWp);
            if (exp_pv_max_total_kWp_set) Global::set_exp_pv_max_kWp_total(exp_pv_max_total_kWp);
            if (exp_pv_min_kWp_set)     Global::set_exp_pv_min_kWp_roof_sec(exp_pv_min_kWp);
            if (exp_pv_max_kWp_set)     Global::set_exp_pv_max_kWp_roof_sec(exp_pv_max_kWp);
            if (feed_in_tariff_set)     Global::set_feed_in_tariff(feed_in_tariff);
            if (demand_tariff_set)      Global::set_demand_tariff(demand_tariff);
            if (npv_discount_rate_set)  Global::set_npv_discount_rate(npv_discount_rate);
            if (npv_time_horizon_set)   Global::set_npv_time_horizon(npv_time_horizon);
            if (icost_PV_per_kWp_set)   Global::set_inst_cost_PV_per_kWp(icost_PV_per_kWp);
            if (icost_BS_per_kWh_set)   Global::set_inst_cost_BS_per_kWh(icost_BS_per_kWh);
            if (exp_bs_P_comp_mode_set) Global::set_battery_power_computation_mode(exp_bs_P_comp_mode_transl);
            //
            Global::set_use_BS_for_SSR_list(use_BS_for_SSR_list);
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


//
// Switch off unused parameter warning for the following block,
// as the parameter "data" is ignored most of the time
//
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"


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
    static size_t callcounter = 1;
    size_t pos = callcounter - 1; // current position in the array is one behind the count of calls
    if (argc != 3) {
        cerr << "Number of arguments not equal to 3 for one row!" << endl;
        return 1;
    }
    size_t current_time_index = stoul(argv[0]);
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
    size_t current_station_id = stoul(argv[0]);
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
    unsigned long current_cu_id    = stoul(argv[0]);
    unsigned long conn_to_subst_id = stoul(argv[1]);
    unsigned long location_id      = stoul(argv[2]);
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
    size_t current_mu_id  = stoul(argv[0]);
    size_t conn_to_unitID = stoul(argv[1]);
    string* mPointIDStr= new string(argv[2]);
    bool has_demand    = argv[3][0] == '1';
    bool has_feedin    = argv[4][0] == '1';
    bool has_pv_resid  = argv[5][0] == '1';
    bool has_pv_opens  = argv[6][0] == '1';
    bool has_bess      = argv[7][0] == '1';
    bool has_hp        = argv[8][0] == '1';
    bool has_wb        = false; // TODO: wallboxes not implemented jet
    bool has_chp       = argv[9][0] == '1';
    size_t locID       = stoul(argv[10]);

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
    static unsigned long callcounter_timestepID = 1; // internal counter for the timestep ID
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
    unsigned long current_timestepID;
    int current_same_orientation_ts_index;
    string current_orientation;
    try {
        current_timestepID = stoul(argv[0]);
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

    size_t pos = callcounter_timestepID - 1; // the current position is one behind the callcounter
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
     * This is the callback function can be used for getting the global Wind profile.
     *
     * The first argument (data) holds the reference to the target array, where
     * data should be written into.
     *
     * Columns:
     * 0           1
     * TimestepID  wind_profile_value
     */
    static size_t callcounter = 1;
    size_t pos = callcounter - 1; // the current position is one behind the callcounter
    if (argc != 2) {
        cerr << "Number of arguments not equal to 2 for one row!" << endl;
        return 1;
    }

    if (stoul(argv[0]) != callcounter) {
        cerr << "Wrong ordering of the global wind profile values!" << endl;
        return 1;
    }
    try {
        if (argv[1] == NULL) {
            ((float*) data)[pos] = 0.0;
        } else {
            ((float*) data)[pos] = stof(argv[1]);
        }
    } catch (exception& e) {
        cerr << "An error happened during the parsing of the wind profile.\n";
        cerr << " - More details: At time step " << callcounter << endl;
        return 1;
    }

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
    static unsigned long callcounter_timestepID = 1; // internal counter for the timestep ID
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

    if (stoul(argv[0]) != callcounter_timestepID) {
        cerr << "There is one row missing for at least one timestep in the list of heat pump profiles!" << endl;
        return 1;
    }
    if (stoul(argv[2]) != callcounter_timeseries) {
        cerr << "There are missing values for at least one time series in the list of heat pump profiles!" << endl;
        return 1;
    }
    size_t pos = callcounter_timestepID - 1; // the current position is one behind the callcounter
    ((float**) data)[callcounter_timeseries][pos] = stof(argv[1]);

    if (callcounter_timestepID < Global::get_n_timesteps()) {
        callcounter_timestepID++;
    } else {
        callcounter_timestepID = 1;
        callcounter_timeseries++;
    }
    return 0;
}
int load_data_from_central_database_callback_ResGridload(void* data, int argc, char** argv, char** colName) {
    /*
     * This is the callback function can be used for getting the global residual grid load (i.e. load that is not measured by simulated meters).
     *
     * The first argument (data) holds the reference to the target array, where
     * data should be written into.
     *
     * Columns:
     * 0           1
     * TimestepID  P_residual_gridload
     */
    static size_t callcounter = 1;
    size_t pos = callcounter - 1; // the current position is one behind the callcounter
    if (argc != 2) {
        cerr << "Number of arguments not equal to 2 for one row!" << endl;
        return 1;
    }

    if (stoul(argv[0]) != callcounter) {
        cerr << "Wrong ordering of the residual grid load values!" << endl;
        return 1;
    }
    try {
        if (argv[1] == NULL) {
            ((float*) data)[pos] = 0.0;
        } else {
            ((float*) data)[pos] = stof(argv[1]);
        }
    } catch (exception& e) {
        cerr << "An error happened during the parsing of the residual grid load series.\n";
        cerr << " - More details: At time step " << callcounter << endl;
        return 1;
    }

    callcounter++;
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
        usable_area = usable_area / 0.707107f;
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
        global::time_timestep_id = new unsigned long[Global::get_n_timesteps()];
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
        // Load residual netload
        //
        float* new_res_gridload_array = new float[Global::get_n_timesteps()];
        sql_query = "SELECT TimestepID,P_residual_gridload FROM residual_grid_load ORDER BY TimestepID;";
        ret_valF  = sqlite3_exec(dbcon, sql_query.c_str(), load_data_from_central_database_callback_ResGridload, new_res_gridload_array/*Reference to the new array*/, &sqlErrorMsgF);
        if (ret_valF != 0) {
            cerr << "Error when reading the SQL-Table: " << sqlErrorMsgF << endl;
            sqlite3_free(sqlErrorMsgF);
            return false;
        }
        global::residual_gridload_kW = new_res_gridload_array;

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


#pragma GCC diagnostic pop /* activate all warnings again */

