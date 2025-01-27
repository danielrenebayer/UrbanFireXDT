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


bool helper_read_ev_info_json_data(const std::string&);
bool helper_read_ev_profile_data(const std::string&);


//
// loads the global config file
//
bool configld::load_config_file(unsigned long scenario_id, string& filepath) {
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
        // variables that need to be translated to another type
        string start_str    = "";  bool start_str_set    = false;
        string end_str      = "";  bool end_str_set      = false;
        string exp_profile_mode    = ""; bool exp_profile_mode_set    = false;
        string sac_planning_mode   = ""; bool sac_planning_mode_set   = false;
        string exp_bs_P_comp_mode  = ""; bool exp_bs_P_comp_mode_set  = false;
        // variables for PV expansion
        bool   exp_pv_mode_static = false; bool exp_pv_mode_static_set= false;

        //
        // define internal functions (here i.e. a lambda function with complete capture-by-reference)
        auto parse_element = [&](string& element_name, boost::property_tree::ptree& scenario_dict) -> void {
            if      ( element_name.compare("data input path")           == 0 )
            {
                string path = scenario_dict.get_value<string>();
                Global::set_input_path( &path );
            }
            else if ( element_name.compare("data output path")          == 0 )
            {
                string path = scenario_dict.get_value<string>();
                Global::set_output_path( &path );
            }
            else if ( element_name.compare("database name")             == 0 )
            {
                string fname = scenario_dict.get_value<string>();
                Global::set_structure_database_name( &fname );
            }
            else if ( element_name.compare("start")                     == 0 )
            {
                start_str       = scenario_dict.get_value<string>();
                start_str_set   = true;
            }
            else if ( element_name.compare("end")                       == 0 )
            {
                end_str         = scenario_dict.get_value<string>();
                end_str_set     = true;
            }
            else if ( element_name.compare("time steps per hour")       == 0 )
            {
                Global::set_tsteps_per_hour( scenario_dict.get_value<int>() );
            }
            else if ( element_name.compare("expansion id")              == 0 )
            {
                Global::set_expansion_scenario_id( scenario_dict.get_value<unsigned long>() );
            }
            else if ( element_name.compare("expansion BS P in kW")      == 0 )
            {
                Global::set_exp_bess_kW( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("expansion BS E in kWh")     == 0 )
            {
                Global::set_exp_bess_kWh( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("expansion BS initial SOC")  == 0 )
            {
                Global::set_exp_bess_start_soc( scenario_dict.get_value<float>() );
            }
            else if (element_name.compare("expansion BS E:P ratio")     == 0)
            {
                Global::set_exp_bess_E_P_ratio( scenario_dict.get_value<float>() );
            }
            else if (element_name.compare("expansion BS capacity computation mode")     == 0)
            {
                string selection = scenario_dict.get_value<string>();
                if (selection == "const") {
                    Global::set_battery_capacity_computation_mode(global::BatteryCapacityComputationMode::Constant);
                } else if (selection == "use PV power") {
                    Global::set_battery_capacity_computation_mode(global::BatteryCapacityComputationMode::BasedOnNominalPVPower);
                } else if (selection == "use mean annual consumption") {
                    Global::set_battery_capacity_computation_mode(global::BatteryCapacityComputationMode::BasedOnAnnualConsumption);
                } else if (selection == "use mean annual consumption with heat pump") {
                    Global::set_battery_capacity_computation_mode(global::BatteryCapacityComputationMode::BasedOnAnnualConsumptionWithHeatPump);
                } else {
                    cerr << "Parameter 'expansion BS capacity computation mode' is defined as '" << selection << "' in config-json, but this value is unknown." << endl;
                    throw runtime_error("Parameter 'expansion BS capacity computation mode' as defined in config-json is unknown.");
                }
            }
            else if (element_name.compare("expansion BS max capacity") == 0)
            {
                Global::set_exp_bess_max_capacity( scenario_dict.get_value<float>() );
            }
            else if (element_name.compare("expansion BS capacity sizing factor for PV") == 0)
            {
                Global::set_exp_bess_sizingE_boPV( scenario_dict.get_value<float>() );
            }
            else if (element_name.compare("expansion BS efficiency in")  == 0)
            {
                Global::set_exp_bess_effi_in( scenario_dict.get_value<float>() );
            }
            else if (element_name.compare("expansion BS efficiency out") == 0)
            {
                Global::set_exp_bess_effi_out( scenario_dict.get_value<float>() );
            }
            else if (element_name.compare("expansion BS self-discharge per ts") == 0)
            {
                Global::set_exp_bess_self_ds_ts( scenario_dict.get_value<float>() );
            }
            else if (element_name.compare("expansion BS power for SOC 0")             == 0)
            {
                Global::set_exp_bess_P_for_SOC_0( scenario_dict.get_value<float>() );
            }
            else if (element_name.compare("expansion BS power for SOC 1")             == 0)
            {
                Global::set_exp_bess_P_for_SOC_1( scenario_dict.get_value<float>() );
            }
            else if (element_name.compare("expansion BS max total E addition")        == 0)
            {
                Global::set_exp_bess_max_E_total( scenario_dict.get_value<float>() );
            }
            else if (element_name.compare("expansion BS max total P addition")        == 0)
            {
                Global::set_exp_bess_max_P_total( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("open space PV kWp")         == 0 )
            {
                Global::set_open_space_pv_kWp( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("open space wind kWp")       == 0 )
            {
                Global::set_wind_kWp( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("expansion profile selection")             == 0 )
            {
                exp_profile_mode       = scenario_dict.get_value<string>();
                exp_profile_mode_set   = true;
            }
            else if ( element_name.compare("CU selection mode for comp. add.")        == 0 )
            {
                sac_planning_mode      = scenario_dict.get_value<string>();
                sac_planning_mode_set  = true;
            }
            else if ( element_name.compare("expansion PV kWp static")                 == 0 )
            {
                Global::set_exp_pv_kWp_static( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("expansion PV min kWp for section usage")  == 0 )
            {
                Global::set_exp_pv_min_kWp_roof_sec( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("expansion PV max inst kWp per section")   == 0 )
            {
                Global::set_exp_pv_max_kWp_roof_sec( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("expansion PV kWp per roof area in m2")    == 0 )
            {
                Global::set_exp_pv_kWp_per_m2( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("expansion PV kWp static mode")            == 0 )
            {
                exp_pv_mode_static     = scenario_dict.get_value<bool>();
                exp_pv_mode_static_set = true;
            }
            else if ( element_name.compare("tariff feed-in per kWh")                  == 0 )
            {
                Global::set_feed_in_tariff( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("tariff demand per kWh")                   == 0 )
            {
                Global::set_demand_tariff( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("emissions per kWh")                       == 0 )
            {
                Global::set_emissions_g_CO2eq_per_kWh( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("net present value discount rate")         == 0 )
            {
                Global::set_npv_discount_rate( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("net present value time horizon in years") == 0 )
            {
                Global::set_npv_time_horizon( scenario_dict.get_value<unsigned int>() );
            }
            else if ( element_name.compare("installation cost PV per kWp")            == 0 )
            {
                Global::set_inst_cost_PV_per_kWp( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("installation cost BS per kWh")            == 0 )
            {
                Global::set_inst_cost_BS_per_kWh( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("expansion BS power computation mode")     == 0 )
            {
                exp_bs_P_comp_mode       = scenario_dict.get_value<string>();
                exp_bs_P_comp_mode_set   = true;
            }
            else if ( element_name.compare("expansion PV max inst kWp per unit")      == 0 )
            {
                Global::set_exp_pv_max_kWp_per_unit( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("expansion PV max total kWp addition")     == 0 )
            {
                Global::set_exp_pv_max_kWp_total( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("expansion PV static mode profile orientation") == 0 )
            {
                string value = scenario_dict.get_value<string>();
                Global::set_exp_pv_static_profile_orientation( &value );
            }
            else if ( element_name.compare("expansion PV static mode profile index") == 0 )
            {
                Global::set_exp_pv_static_profile_idx( scenario_dict.get_value<int>() );
            }
            else if ( element_name.compare("HP flexibility in ts") == 0 )
            {
                Global::set_hp_flexibility_in_ts( scenario_dict.get_value<unsigned int>() );
            }
            else if ( element_name.compare("th. E to HP el. E conversion factor") == 0 )
            {
                Global::set_heat_demand_thermalE_to_hpE_conv_f( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("Heat consumption apo building V: slope") == 0 )
            {
                Global::set_heat_cons_bobv_slope( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("Heat consumption apo building V: intercept") == 0 )
            {
                Global::set_heat_cons_bobv_intercept( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("EV plugin probability") == 0)
            {
                Global::set_ev_plugin_probability( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("EV battery size kWh") == 0)
            {
                Global::set_ev_battery_size_kWh( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("EV consumption kWh per km") == 0)
            {
                Global::set_ev_consumption_kWh_km( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("EV max charging power") == 0)
            {
                Global::set_ev_max_charging_power_kW( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("EV charging efficiency") == 0)
            {
                Global::set_ev_charging_effi( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("ev data path")         == 0 )
            {
                string value = scenario_dict.get_value<string>();
                Global::set_ev_data_path( &value );
            }
            else if ( element_name.compare("use emission time series ia")             == 0 )
            {
                Global::set_use_emission_time_series_ia( scenario_dict.get_value<bool>() );
            }
            else if ( element_name.compare("use prices time series ia")               == 0 )
            {
                Global::set_use_prices_time_series_ia( scenario_dict.get_value<bool>() );
            }
            else if ( element_name.compare("controller mode")                         == 0 )
            {
                string selection = scenario_dict.get_value<string>();
                if (selection == "rule-based") {
                    Global::set_controller_mode( global::ControllerMode::RuleBased );
                } else if (selection == "opti with perfect forecast") {
                    Global::set_controller_mode( global::ControllerMode::OptimizedWithPerfectForecast );
                } else {
                    cerr << "Parameter 'controller mode' is defined as '" << selection << "' in config-json, but this value is unknown." << endl;
                    throw runtime_error("Parameter 'controller mode' as defined in config-json is unknown.");
                }
            }
            else if ( element_name.compare("control horizon in ts")                  == 0 )
            {
                Global::set_control_horizon_in_ts( scenario_dict.get_value<unsigned int>() );
            }
            else if ( element_name.compare("control update freq in ts")               == 0 )
            {
                Global::set_control_update_freq_in_ts( scenario_dict.get_value<unsigned int>() );
            }
            else if ( element_name.compare("controller allow bs charging from grid")  == 0 )
            {
                Global::set_controller_allow_bs_grid_charging( scenario_dict.get_value<bool>() );
            }
            else if ( element_name.compare("controller optimization target")          == 0 )
            {
                string selection = scenario_dict.get_value<string>();
                if (selection == "electricity costs") {
                    Global::set_controller_optimization_target( global::ControllerOptimizationTarget::ElectricityCosts );
                } else if (selection == "peak load") {
                    Global::set_controller_optimization_target( global::ControllerOptimizationTarget::PeakLoad );
                } else if (selection == "emissions") {
                    Global::set_controller_optimization_target( global::ControllerOptimizationTarget::Emissions );
                } else {
                    cerr << "Parameter 'controller optimization target' is defined as '" << selection << "' in config-json, but this value is unknown." << endl;
                    throw runtime_error("Parameter 'controller optimization target' as defined in config-json is unknown.");
                }
            }
            else if ( element_name.compare("select buildings with given heat demand only") == 0 )
            {
                Global::set_select_buildings_wg_heatd_only( scenario_dict.get_value<bool>() );
            }
            else if ( element_name.compare("annual heat demand limit for selection") == 0 )
            {
                Global::set_annual_heat_demand_limit_fsac( scenario_dict.get_value<float>() );
            }
            else if ( element_name.compare("expansion HP max total addition") == 0 )
            {
                Global::set_exp_hp_max_n_addition( scenario_dict.get_value<unsigned long>() );
            }
            else if ( element_name.compare("expansion EV max total addition") == 0 )
            {
                Global::set_exp_ev_max_n_addition( scenario_dict.get_value<unsigned long>() );
            }
            else if ( element_name.compare("break SAC loop if limit reached") == 0 )
            {
                Global::set_break_sac_loop_if_limit_reached( scenario_dict.get_value<bool>() );
            }
            else if ( element_name.compare("select only residential buildings") == 0 )
            {
                Global::set_select_only_residential_buildings( scenario_dict.get_value<bool>() );
            }
            else if ( element_name.compare("id") == 0 )
            {}
            else if ( element_name.starts_with("comment") == 0 )
            {}
            else if ( element_name.starts_with("__disabled ") == 0 )
            {}
            else if ( element_name.compare("inherits from") == 0 )
            {}
            else
            {
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
        bool scenario_found = true;
        // define lambda function for searching and parsing the scenario entries for the correct ID
        auto find_sceario_id_and_parse = [&](unsigned long scenario_id_to_find) -> bool {
            for (auto& scenario_dict_all : tree_root.get_child("Scenarios")) {
                auto scenario_dict = scenario_dict_all.second;
                // if we have found the correct entry ...
                if (scenario_dict.get<unsigned long>("id") == scenario_id_to_find) {
                    // ... we read all variables
                    for (auto& s : scenario_dict) {
                        string element_name = s.first;
                        parse_element(element_name, s.second);
                    }
                    //
                    return true;
                }
            }
            return false;
        };
        list<unsigned long> scenarios_to_load;
        scenarios_to_load.push_front(scenario_id);
        // get all scenario IDs from which the selected one inherits
        bool inheritance_ended = false;
        unsigned long current_search_scenario_id = scenario_id;
        while (!inheritance_ended) {
            inheritance_ended = true;
            for (auto& scenario_dict_all : tree_root.get_child("Scenarios")) {
                auto scenario_dict = scenario_dict_all.second;
                // if we have found the correct entry ...
                if (scenario_dict.get<unsigned long>("id") == current_search_scenario_id) {
                    auto e = scenario_dict.get_optional<int>("inherits from");
                    // check if there is a scenario from which we inherited
                    if (e.is_initialized()) {
                        int upper_scenario = e.get();
                        if (find(scenarios_to_load.begin(), scenarios_to_load.end(), upper_scenario) != scenarios_to_load.end()) {
                            cerr << "Error in config file: Ring closure in the inheritance for scenario ID " << upper_scenario << "!" << endl;
                            return false;
                        }
                        cout << "Reading settings for inherited scenario with ID " << upper_scenario << endl;
                        scenarios_to_load.push_front(upper_scenario);
                        current_search_scenario_id = upper_scenario;
                        inheritance_ended = false;
                    } else {
                        inheritance_ended = true;
                    }
                    break; // quit the inner loop
                }
            }
        }
        // load all required scenario definitions
        for (unsigned long s : scenarios_to_load) {
            if (! find_sceario_id_and_parse(s) ) {
                cerr << "Scenario " << s << " was not found in the config file!" << endl;
                scenario_found = false;
                break;
            }
        }

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
        // Finally, add to global variable collection
        if (scenario_found) {
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
            if (exp_profile_mode_set)   Global::set_exp_profile_mode(exp_profile_mode_transf);
            if (sac_planning_mode_set)  Global::set_cu_selection_mode_fca(sac_planning_mode_transl);
            if (exp_pv_mode_static_set) Global::set_exp_pv_mode(exp_pv_mode_static);
            if (exp_bs_P_comp_mode_set) Global::set_battery_power_computation_mode(exp_bs_P_comp_mode_transl);
            //
            Global::LockAllVariables();
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
     * The first parameter (data) is a reference to a function - see query_list
     * 
     * Parameter information taken from documentation:
     *   argc: holds the number of results
     *   argv: holds each value in array
     *   colName: holds each column returned in array
     */
    if (argc != 1) {
        cerr << "Number of arguments not equal to 1 for one row (in callback A)!" << endl;
        return 1;
    }

    // convert the first argument to the correct type
    void (*fp) (unsigned long) = (void (*)(unsigned long)) data;
    unsigned long val = std::stoul( argv[0] );
    // call the anonymous function / method
    (*fp)(val);

    return 0;
}
int load_data_from_central_database_callbackB(void* data, int argc, char** argv, char** colName) {
    /* 
     * This is the callback function for loading the time indices
     *
     * Columns:
     * 0           1              2              3
     * TimestepID  local_time_ra  local_time_la  local_time_zone
     * 
     */
    static size_t callcounter = 1;
    size_t pos = callcounter - 1; // current position in the array is one behind the count of calls
    if (argc != 4) {
        cerr << "Number of arguments not equal to 4 for one row!" << endl;
        return 1;
    }
    size_t current_time_index = stoul(argv[0]);
    if (current_time_index != callcounter) {
        cerr << "Time indices are not ordered sequentially!" << endl;
        return 1;
    }
    // convert time values
    struct tm* time_value_ra = new struct tm;
    struct tm* time_value_la = new struct tm;
    stringstream stream_time_value_ra(argv[1]);
    stringstream stream_time_value_la(argv[2]);
    stream_time_value_ra >> get_time(time_value_ra, "%Y-%m-%d %H:%M:%S");
    stream_time_value_la >> get_time(time_value_la, "%Y-%m-%d %H:%M:%S");
    string timezone_str { argv[3] };
    if (timezone_str == "CEST") {
        time_value_ra->tm_isdst = 1;
        time_value_la->tm_isdst = 1;
    } else {
        time_value_ra->tm_isdst = 0;
        time_value_la->tm_isdst = 0;
    }
    // add time values to global list
    global::time_timestep_id[pos] = current_time_index;
    global::time_localtime_r->push_back(time_value_ra);
    global::time_localtime_l->push_back(time_value_la);
    global::time_localtimezone_str->push_back(timezone_str);
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
        if (!Substation::InstantiateNewSubstation(current_station_id, stationName)) {
            cerr << "Error when creating substation with id " << current_station_id << endl;
            cerr << "Is the ID of the substation unique?" << endl;
            return 1;
        }
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
     * 0       1              2      3                          4
     * UnitID  substation_id  LocID  has_residential_buildings  n_flats
     * 
     */
    if (argc != 5) {
        cerr << "Number of arguments not equal to 5 for one row!" << endl;
        return 1;
    }
    static unsigned long last_valid_cu_id = 0;
    try {
        unsigned long current_cu_id    = stoul(argv[0]);
        unsigned long conn_to_subst_id = stoul(argv[1]);
        unsigned long location_id      = stoul(argv[2]);
        unsigned int  n_flats          =  stoi(argv[4]);
        bool residential = (location_id > 0) && (strlen(argv[3]) > 0) && (stoi(argv[3]) > 0);
        if (!ControlUnit::InstantiateNewControlUnit(current_cu_id, conn_to_subst_id, location_id, residential, n_flats)) {
            cerr << "Error when creating control unit with id " << current_cu_id << endl;
            cerr << "Is the ID of the control unit unique?" << endl;
            return 1;
        }
        last_valid_cu_id = current_cu_id;
    } catch (std::invalid_argument const& e) {
        cerr << "Non-parsable integer or unsigned long detected in list_of_control_units. Last valid ID = " << last_valid_cu_id << std::endl;
        return 1;
    } catch (std::out_of_range const& e) {
        cerr << "An integer is out of bounds in list_of_control_units. Last valid ID = " << last_valid_cu_id << std::endl;
        return 1;
    } catch (runtime_error& e) {
        cerr << "Error when creating control unit. Last valid ID = " << last_valid_cu_id << endl;
        cerr << "Details:" << endl;
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
     * 6                  7         8       9        10     11        12
     * has_pv_open_space, has_bess, has_hp, has_chp, LocID, has_wind, has_biomass
     * 13        14
     * has_evcs  has_public_evcs
     */
    if (argc != 15) {
        cerr << "Number of arguments not equal to 15 for one row!" << endl;
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
    bool has_chp       = argv[9][0] == '1';
    bool has_wind      = argv[11][0] == '1';
    bool has_evchst    = argv[13][0] == '1';
  //bool has_pub_evchst= argv[14][0] == '1'; // TODO -> implement public EV charging stations
    size_t locID       = stoul(argv[10]);

    stringstream data_input_path;
    data_input_path << Global::get_input_path() << "SeparatedSmartMeterData/";
    data_input_path << current_mu_id;
    data_input_path << ".csv";

    try {
        bool res = MeasurementUnit::InstantiateNewMeasurementUnit(
                                current_mu_id, conn_to_unitID, mPointIDStr, locID,
                                has_demand, has_feedin, has_pv_resid, has_pv_opens,
                                has_bess,   has_hp,     has_wind,     has_evchst,
                                has_chp, data_input_path.str());
        if (!res) {
            cerr << "Error when creating measurement unit with id " << current_mu_id << endl;
            cerr << "Is the ID of the control unit unique?" << endl;
            return 1;
        }
        /*MeasurementUnit* newMU = */ MeasurementUnit::GetInstancePublicID(current_mu_id);
    } catch (runtime_error& e) {
        cerr << "Error when creating measurement unit with id " << current_mu_id << endl;
        cerr << "Details:" << endl;
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
unsigned long callcounter_callback_HP = 0;
int load_data_from_central_database_callback_HP(void* data, int argc, char** argv, char** colName) {
    /*
     * This is the callback function for geeting the global heat pump profiles.
     *
     * The first argument (data) holds the reference to the target array, where
     * data should be written into.
     *
     * Columns:
     *   0           1                     2
     *   TimestepID  ShiftableDemand_kW    TimeSeriesIndex
     * or
     *               UnshiftableDemand_kW
     */
    static unsigned long callcounter_timestepID = 1; // internal counter for the timestep ID
    //static unsigned long callcounter_timeseries = 0; // internal counter for the timeseries, both are used to identify missing values
    // the latter callcounter is defined outside to make it resettable!

    // check callcounter
    if (callcounter_callback_HP >= Global::get_n_heatpump_profiles()) {
        cerr << "Error in data for heat pumps.\n";
        cerr << "It is said that there are " << Global::get_n_heatpump_profiles() << " heat pump profiles, but therer are more available in the table." << endl;
        return 1;
    }

    if (argc != 3) {
        cerr << "Number of arguments not equal to 3 for one row!" << endl;
        return 1;
    }

    try {
        if (stoul(argv[0]) != callcounter_timestepID) {
            cerr << "There is one row missing for at least one timestep in the list of heat pump profiles!" << endl;
            return 1;
        }
        if (stoul(argv[2]) != callcounter_callback_HP) {
            cerr << "There are missing values for at least one time series in the list of heat pump profiles!" << endl;
            return 1;
        }
        size_t pos = callcounter_timestepID - 1; // the current position is one behind the callcounter
        ((float**) data)[callcounter_callback_HP][pos] = stof(argv[1]);
    } catch (exception& e) {
        cerr << "An error happened during the parsing of the heat pump profiles.\n";
        cerr << " - More details: At time step " << callcounter_timestepID << " for time series " << callcounter_callback_HP << endl;
        return 1;
    }

    if (callcounter_timestepID < Global::get_n_timesteps()) {
        callcounter_timestepID++;
    } else {
        callcounter_timestepID = 1;
        callcounter_callback_HP++;
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
     * 0      1            2           3
     * LocID  n_buildings  max_volume  Heat Demand in kWh per year
     */
    if (argc != 4) {
        cerr << "Number of arguments not equal to 4 for one row!" << endl;
        return 1;
    }
    
    size_t locationID = 0;
    float  annual_heat_demand_kWh = -1.0;
    float  volume = 0.0;

    try {
        locationID = stoul(argv[0]);

        if (argv[2] != NULL) {
            volume = stof(argv[2]);
        }
        if (argv[3] != NULL) {
            annual_heat_demand_kWh = stof(argv[3]);
        }
    } catch (exception& e) {
        cerr << "An error happened during the parsing of the heat demand abd building volume information.\n";
        cerr << " - More details: " << e.what() << endl;
        return 1;
    }

    global::annual_heat_demand_kWh[ locationID ] = annual_heat_demand_kWh;
    global::building_volumes_m3[    locationID ] = volume;
    global::locations_with_geodata.insert( locationID );

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
int load_data_from_central_database_callback_emissions(void* data, int argc, char** argv, char** colName) {
    /*
     * This is the callback function for geeting the emission time series.
     *
     * Columns:
     * 0           1
     * TimestepID  emissions_g_kWh
     */
    static size_t callcounter = 1;
    size_t pos = callcounter - 1; // the current position is one behind the callcounter
    if (argc != 2) {
        cerr << "Number of arguments not equal to 2 for one row!" << endl;
        return 1;
    }

    if (stoul(argv[0]) != callcounter) {
        cerr << "Wrong ordering, duplicated or missing items in the emissions time series!" << endl;
        return 1;
    }
    try {
        if (argv[1] == NULL) {
            ((float*) data)[pos] = 0.0;
        } else {
            ((float*) data)[pos] = stof(argv[1]);
        }
    } catch (exception& e) {
        cerr << "An error happened during the parsing of the emissions time series.\n";
        cerr << " - More details: At time step " << callcounter << endl;
        return 1;
    }
    callcounter++;
    return 0;
}
int load_data_from_central_database_callback_prices(void* data, int argc, char** argv, char** colName) {
    /*
     * This is the callback function for geeting the emission time series.
     *
     * Columns:
     * 0           1            2
     * TimestepID  local_price  spotmarket_price
     */
    static size_t callcounter = 1;
    size_t pos = callcounter - 1; // the current position is one behind the callcounter
    if (argc != 3) {
        cerr << "Number of arguments not equal to 3 for one row!" << endl;
        return 1;
    }

    if (stoul(argv[0]) != callcounter) {
        cerr << "Wrong ordering, duplicated or missing items in the prices time series!" << endl;
        return 1;
    }
    try {
        if (argv[1] == NULL) {
            ((float*) (((float**)data)[0]) )[pos] = 0.0;
        } else {
            float local_price = stof(argv[1]);
            if (local_price < Global::get_feed_in_tariff()) {
                cerr << "Warning: local energy price is lower than the feed-in tariff at time step " << callcounter << ".\n";
                cerr << "Setting it to 1.05%% of the feed-in tariff." << endl;
                local_price = Global::get_feed_in_tariff() * 1.05f;
            }
            ((float*) (((float**)data)[0]) )[pos] = local_price;
        }
        if (argv[2] == NULL) {
            ((float*) (((float**)data)[1]) )[pos] = 0.0;
        } else {
            ((float*) (((float**)data)[1]) )[pos] = stof(argv[2]);
        }
    } catch (exception& e) {
        cerr << "An error happened during the parsing of the prices time series.\n";
        cerr << " - More details: At time step " << callcounter << endl;
        return 1;
    }
    callcounter++;
    return 0;
}
int sql_check_if_table_exists_callback(void* data, int argc, char** argv, char** colName) {
    /*
     * This is the callback function for checking, if a table exists or not
     *
     * @param data: A reference to a boolean variable that will be set to true if the table exists
     *
     * Columns:
     * 0
     * table_name
     */
    *((bool*) data) = true;
    return 0;
}
bool configld::load_data_from_central_database(const char* filepath) {
    if (! filesystem::exists( filesystem::path(filepath) ) ) {
        cerr << "Structure database file " << filepath << " not found!" << endl;
        return false;
    }

    sqlite3* dbcon;
    int rc = sqlite3_open(filepath, &dbcon);

    if (rc == 0) {
        cout << "Loading system structure from database ..." << endl;
        cout << "Using database " << filepath << endl;

        //
        // Load general data
        //
        string sql_queryA;
        char*  sqlErrorMsgA;
        int    ret_valA;
        // the tables to query and the functions of the calls where to store this information
        // pair< pair< TABLE_NAME_TO_QUERY, TRUE_IF_ONLY_SELECT_COUNT >, ANONYMOUS_FUNCTION_NAME_TO_SAFE_THE_RESULT >
        std::list<std::pair<std::pair<const char *, bool>, void (*)(unsigned long)>>
        /*auto*/ query_list = {
            std::make_pair(std::make_pair("time_indices", true), &Global::set_n_timesteps),
            std::make_pair(std::make_pair("list_of_substations", true), &Global::set_n_substations),
            std::make_pair(std::make_pair("list_of_control_units", true), &Global::set_n_CUs),
            std::make_pair(std::make_pair("list_of_measurement_units", true), &Global::set_n_MUs),
            std::make_pair(std::make_pair("global_profiles_heatpumps", false), &Global::set_n_heatpump_profiles)
        };
        //auto u = Global::set_n_timesteps;
        // execute the query for every entry in the query_list
        for (auto e : query_list) {
            // check, which query should be executed
            if (e.first.second) {
                sql_queryA  = "SELECT count(*) FROM ";
                sql_queryA += e.first.first;
                sql_queryA += ";";
            } else {
                sql_queryA  = "SELECT count(*) FROM ( SELECT DISTINCT TimeSeriesIndex FROM ";
                sql_queryA += e.first.first;
                sql_queryA += " );";
            }
            // execute the query and call the function of the second argument of the pair to save the results
            ret_valA   = sqlite3_exec(dbcon, sql_queryA.c_str(), load_data_from_central_database_callbackA, (void*) e.second, &sqlErrorMsgA);
            if (ret_valA != 0) {
                cerr << "Error when executing command '" << sql_queryA << "': " << sqlErrorMsgA << endl;
                sqlite3_free(sqlErrorMsgA);
                return false;
            }
        }

        //
        // Initialize global time list
        //
        global::time_timestep_id = new unsigned long[Global::get_n_timesteps()];
        global::time_localtime_r = new vector<struct tm*>();
        global::time_localtime_l = new vector<struct tm*>();
        global::time_localtimezone_str = new vector<string>();
        //
        // Load time indices
        //
        string sql_queryB = "SELECT TimestepID, local_time_ra, local_time_la, local_time_zone FROM time_indices ORDER BY TimestepID;";
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
        // Process time-related data
        // -- check if time-related variables have been initialized
        if (!Global::CheckTimeRelatedVariablesInitState()) {
            std::cerr << "A time-related variable has not been initialized! Stopping." << std::endl;
            return false;
        }
        // Calculation of fist and last time step as given in data
        unsigned long n_tsteps = Global::get_n_timesteps();
        struct tm* tm_start = Global::get_ts_start_tm();
        struct tm* tm_end   = Global::get_ts_end_tm();
        unsigned long first_sim_ts = 0;
        unsigned long last_sim_ts = 0;
        bool sim_started  = false; // gets true, if simulation range (as given by tm_start) has been reached
        bool sim_finished = false; // gets true, if the enf of the simulation range (as given by tm_end) has been reached
        for (unsigned long tsID = 0; tsID < n_tsteps; tsID++) {
            // get current time as struct tm
            struct tm* current_tm = global::time_localtime_r->at(tsID);
            // jump time steps if they are not inside the simulation range
            if (sim_started) {
                if (compare_struct_tm(current_tm, tm_end) >= 0) {
                    last_sim_ts = tsID + 1;
                    sim_finished = true;
                    break;
                }
            } else {
                if (compare_struct_tm(current_tm, tm_start) >= 0) {
                    sim_started = true;
                    first_sim_ts = tsID + 1;
                } else {
                    continue;
                }
            }
        }
        // set last time step to end of simulation if not reached previously
        if (!sim_finished) {
            last_sim_ts = n_tsteps;
        }
        std::cout << "Last step in simulation run(s)  = " << std::setw(6) << last_sim_ts << "\n";
        std::cout << "First step in simulation run(s) = " << std::setw(6) << first_sim_ts << "\n";
        Global::set_first_timestep(first_sim_ts);
        Global::set_last_timestep(last_sim_ts);

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
        string sql_queryD = "SELECT A.UnitID, A.substation_id, A.LocID, B.has_residential_buildings, A.n_flats FROM list_of_control_units AS A LEFT JOIN (SELECT LocID, has_residential_buildings FROM address_data) AS B ON A.LocID = B.LocID ORDER BY UnitID;";
        char* sqlErrorMsgD;
        int ret_valD = sqlite3_exec(dbcon, sql_queryD.c_str(), load_data_from_central_database_callbackD, NULL, &sqlErrorMsgD);
        if (ret_valD != 0) {
            cerr << "Error when reading the SQL-Table: " << sqlErrorMsgD << endl;
            sqlite3_free(sqlErrorMsgD);
            return false;
        }
        // 3. MUs
        string sql_queryE = "SELECT MeUID, UnitID, MeterPointID, has_demand, has_feedin, has_pv_residential, has_pv_open_space, has_bess, has_hp, has_chp, LocID, has_wind, has_biomass, has_evcs, has_public_evcs FROM list_of_measurement_units ORDER BY MeUID;";
        char* sqlErrorMsgE;
        int ret_valE = sqlite3_exec(dbcon, sql_queryE.c_str(), load_data_from_central_database_callbackE, NULL, &sqlErrorMsgE);
        if (ret_valE != 0) {
            cerr << "Error when reading the SQL-Table: " << sqlErrorMsgE << endl;
            sqlite3_free(sqlErrorMsgE);
            return false;
        }
        // 3b. Load MU data
        if (!MeasurementUnit::LoadDataForAllInstances()) {
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
            // initialize with 0 by default
            for (unsigned long l = 0; l < Global::get_n_timesteps(); l++)
                new_pv_array[pv_idx][l] = 0;
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
        // initialize with 0 by default
        for (unsigned long l = 0; l < Global::get_n_timesteps(); l++)
            new_wind_array[l] = 0;
        // run query
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
        float** new_hp_profile_s_array   = new float*[Global::get_n_heatpump_profiles()];
        float** new_hp_profile_nos_array = new float*[Global::get_n_heatpump_profiles()];
        double** new_hp_profile_s_cumsum  = new double*[Global::get_n_heatpump_profiles()];
        for (unsigned long hp_idx = 0; hp_idx < Global::get_n_heatpump_profiles(); hp_idx++) {
            new_hp_profile_s_array[hp_idx] = new float[Global::get_n_timesteps()];
            new_hp_profile_nos_array[hp_idx] = new float[Global::get_n_timesteps()];
            new_hp_profile_s_cumsum[hp_idx]  = new double[Global::get_n_timesteps()];
            // initialize with 0 by default
            for (unsigned long l = 0; l < Global::get_n_timesteps(); l++) {
                new_hp_profile_s_array[hp_idx][l]   = 0;
                new_hp_profile_nos_array[hp_idx][l] = 0;
                new_hp_profile_s_cumsum[hp_idx][l]  = 0;
            }
        }
        callcounter_callback_HP = 0;
        sql_query = "SELECT TimestepID,ShiftableDemand_kW,TimeSeriesIndex FROM global_profiles_heatpumps ORDER BY TimeSeriesIndex,TimestepID;";
        ret_valF  = sqlite3_exec(dbcon, sql_query.c_str(), load_data_from_central_database_callback_HP, new_hp_profile_s_array  /*Reference to the new array*/, &sqlErrorMsgF);
        if (ret_valF != 0) {
            cerr << "Error when reading the SQL-Table: " << sqlErrorMsgF << endl;
            sqlite3_free(sqlErrorMsgF);
            return false;
        }
        callcounter_callback_HP = 0;
        sql_query = "SELECT TimestepID,UnshiftableDemand_kW,TimeSeriesIndex FROM global_profiles_heatpumps ORDER BY TimeSeriesIndex,TimestepID;";
        ret_valF  = sqlite3_exec(dbcon, sql_query.c_str(), load_data_from_central_database_callback_HP, new_hp_profile_nos_array/*Reference to the new array*/, &sqlErrorMsgF);
        if (ret_valF != 0) {
            cerr << "Error when reading the SQL-Table: " << sqlErrorMsgF << endl;
            sqlite3_free(sqlErrorMsgF);
            return false;
        }
        // compute cumulative sum per profile
        for (unsigned long hp_idx = 0; hp_idx < Global::get_n_heatpump_profiles(); hp_idx++) {
            double cumsum = 0.0;
            for (unsigned long l = 0; l < Global::get_n_timesteps(); l++) {
                // If first simulation time step >= current time step l, set added value to 0
                if (l < Global::get_first_timestep() - 1) { // l is the ID, not the timestep number (starting at 1)
                    cumsum = 0.0;
                } else {
                    cumsum += new_hp_profile_s_array[hp_idx][l];
                }
                new_hp_profile_s_cumsum[hp_idx][l] = cumsum;
            }
        }
        // allocate to globale variables
        global::hp_profiles_shiftable = new_hp_profile_s_array;
        global::hp_profiles_not_shift = new_hp_profile_nos_array;
        global::hp_profiles_s_cumsum  = new_hp_profile_s_cumsum;

        //
        // Load residual netload
        //
        float* new_res_gridload_array = new float[Global::get_n_timesteps()];
        // initialize with 0 by default
        for (unsigned long l = 0; l < Global::get_n_timesteps(); l++)
            new_res_gridload_array[l] = 0;
        // run query
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
        // 1. annual heat demand for heat pumps AND buildings with available geo data?
        sql_query = "SELECT A.LocID, A.n_buildings, A.max_volume, B.MeanHeatEnergy_kWh FROM address_data as A LEFT JOIN heat_demand_per_location as B ON A.LocID = B.LocID WHERE n_buildings >= 1 ORDER BY A.LocID;";
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

        //
        // emissions time series (if available)
        //
        bool table_exists = false;
        sql_query = "SELECT name FROM sqlite_master WHERE type='table' AND name='electricity_emissions';";
        sqlite3_exec(dbcon, sql_query.c_str(), sql_check_if_table_exists_callback, &table_exists, &sqlErrorMsgF);
        if (table_exists && Global::get_use_emission_time_series_ia()) {
            cout << "Loading time series on emission data as this table is present.\n";
            float* new_array = new float[Global::get_n_timesteps()];
            // initialize with 0 by default
            for (unsigned long l = 0; l < Global::get_n_timesteps(); l++)
                new_array[l] = 0;
            // run query
            sql_query = "SELECT TimestepID,emissions_g_kWh FROM electricity_emissions ORDER BY TimestepID;";
            ret_valF = sqlite3_exec(dbcon, sql_query.c_str(), load_data_from_central_database_callback_emissions, new_array, &sqlErrorMsgF);
            if (ret_valF != 0) {
                cerr << "Error when reading the SQL-Table: " << sqlErrorMsgF << endl;
                sqlite3_free(sqlErrorMsgF);
                return false;
            }
            global::emission_ts = new_array;
        }

        //
        // prices time series (if available)
        //
        /*bool*/ table_exists = false;
        sql_query = "SELECT name FROM sqlite_master WHERE type='table' AND name='electricity_prices';";
        sqlite3_exec(dbcon, sql_query.c_str(), sql_check_if_table_exists_callback, &table_exists, &sqlErrorMsgF);
        if (table_exists && Global::get_use_prices_time_series_ia()) {
            cout << "Loading time series on electricity prices as this table is present.\n";
            float* new_arrayA = new float[Global::get_n_timesteps()];
            float* new_arrayB = new float[Global::get_n_timesteps()];
            // initialize with 0 by default
            for (unsigned long l = 0; l < Global::get_n_timesteps(); l++) {
                new_arrayA[l] = 0;
                new_arrayB[l] = 0;
            }
            float* new_arrays[2] = {new_arrayA, new_arrayB};
            // run query
            sql_query = "SELECT TimestepID,local_price,spotmarket_price FROM electricity_prices ORDER BY TimestepID;";
            ret_valF = sqlite3_exec(dbcon, sql_query.c_str(), load_data_from_central_database_callback_prices, new_arrays, &sqlErrorMsgF);
            if (ret_valF != 0) {
                cerr << "Error when reading the SQL-Table: " << sqlErrorMsgF << endl;
                sqlite3_free(sqlErrorMsgF);
                return false;
            }
            global::eprices_local_ts = new_arrayA;
            global::eprices_spotm_ts = new_arrayB;
        }

        //
        // Load EV profiles
        //
        if (Global::get_ev_data_path().length() > 0) {
            cout << "    Start parsing mobility data ...\n";
            // A) Load car data
            if (!helper_read_ev_info_json_data(Global::get_ev_data_path() + "car_info.json"))
                return false;
            // B) Load home-ceneterd car driving profiles
            if (!helper_read_ev_profile_data(Global::get_ev_data_path() + "car_tours.json"))
                return false;
            cout << "    ... finished loading mobility data.\n";
        }

        cout << "... finished loading data and system structure.\n";
        cout << global::output_section_delimiter << endl;

        sqlite3_close(dbcon);
        return true;
    } else {
        cerr << "Error when connecting to the database" << endl;
        return false;
    }
}


#pragma GCC diagnostic pop /* activate all warnings again */


//
// Helper function for parsing the EV data from JSON files
bool helper_read_ev_info_json_data(const std::string& filepath) {
    bpt::ptree tree_root;
    try {
        bpt::read_json(filepath, tree_root);
    } catch (bpt::json_parser_error& j) {
        cerr << "Error when reading EV json file: " << j.what() << endl;
        return false;
    }
    try {
        //
        // parse json file per element
        for (auto& scenario_dict_all : tree_root) {
            if (scenario_dict_all.first != "") {
                cerr << "Error when parsing EV json file: Wrong formatting of file 'car_info.json'!" << endl;
                return false;
            }
            unsigned long carID = scenario_dict_all.second.get_child("carID").get_value<unsigned long>();
            unsigned long cuID  = scenario_dict_all.second.get_child("ControlUnitID").get_value<unsigned long>();
            ControlUnit::GetInstancePublicIDWE(cuID)->add_ev(carID);
        }
    } catch (bpt::ptree_bad_path& j) {
        cerr << "Error when parsing EV json file: " << j.what() << endl;
        return false;
    }
    return true;
}

//
// Helper function for parsing home-ceneterd car driving profiles
bool helper_read_ev_profile_data(const std::string& filepath) {
    bpt::ptree tree_root;
    try {
        bpt::read_json(filepath, tree_root);
    } catch (bpt::json_parser_error& j) {
        cerr << "Error when reading EV-profiles json file: " << j.what() << endl;
        return false;
    }
    try {
        //
        // parse json file per element
        for (auto& scenario_dict_all : tree_root) {
            if (scenario_dict_all.first != "") {
                cerr << "Error when parsing EV-profiles json file: Wrong formatting of file 'car_tours.json'!" << endl;
                return false;
            }
            unsigned long carID          = scenario_dict_all.second.get_child("carID").get_value<unsigned long>();
                     short weekday       = scenario_dict_all.second.get_child("weekday").get_value<short>();
            unsigned long hour_departure = scenario_dict_all.second.get_child("hour_departure").get_value<unsigned long>();
            unsigned long duration_in_ts = scenario_dict_all.second.get_child("duration_in_ts").get_value<unsigned long>();
            double        tour_length_km = scenario_dict_all.second.get_child("tour_length_km").get_value<double>();
            short         with_work      = scenario_dict_all.second.get_child("with_work").get_value<short>();
            //
            EVFSM::AddWeeklyTour(carID, weekday, hour_departure, duration_in_ts, tour_length_km, with_work > 0);
        }
    } catch (bpt::ptree_bad_path& j) {
        cerr << "Error when parsing EV-profiles json file: " << j.what() << endl;
        return false;
    }
    return true;
}





#define PRINT_VAR(varname) cout << "    " << std::setw(44) << std::left << #varname << " = " << varname << "\n"
#define PRINT_TM_VAR(varname) cout << "    " << std::setw(44) << #varname << " = " << std::put_time(varname, "%F %R") << "\n"
#define PRINT_ENUM_VAR(varname, lambda_expr) cout << "    " << std::setw(44) << #varname << " = " << lambda_expr(varname) << "\n"

//
// Implementation of configld::output_variable_values()
//
void configld::output_variable_values() {
    cout << "Simulation information:\n";
    cout << "    Simulation build at " << __DATE__ << " " << __TIME__ <<  "\n";
    #ifdef __GNUC__
    cout << "    GCC was used as compiler.\n    GCC Version = " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << "\n";
    #endif
    #ifdef __VERSION__
    cout << "    Compiler version = " << __VERSION__ << "\n";
    #endif
    #ifdef __OPTIMIZE__
    cout << "    Optimization was enabled during compile time.\n";
    #endif
    cout << "    C++ standard = " << __cplusplus << "\n\n";
    cout << "List of parameter settings:\n";
    // Scenario selection
    cout << "  Scenario selection and simulation flow control:\n";
    PRINT_VAR(Global::get_expansion_scenario_id());
    PRINT_VAR(Global::is_parameter_variation());
    PRINT_VAR(Global::get_parameter_varID());
    PRINT_VAR(Global::get_repetitions_selected());
    PRINT_VAR(Global::is_seed_set());
    if (Global::is_seed_set()) { PRINT_VAR(Global::get_seed()); }
    PRINT_VAR(Global::get_n_repetitions());
    PRINT_VAR(Global::get_n_threads());
    PRINT_TM_VAR(Global::get_ts_start_tm());
    PRINT_TM_VAR(Global::get_ts_end_tm());
    // Data
    cout << "  Data:\n";
    PRINT_VAR(Global::get_output_path());
    PRINT_VAR(Global::get_input_path());
    PRINT_VAR(Global::get_structure_database_name());
    PRINT_VAR(Global::get_n_timesteps());
    PRINT_VAR(Global::get_tsteps_per_hour());
    PRINT_VAR(Global::get_n_substations());
    PRINT_VAR(Global::get_n_CUs());
    PRINT_VAR(Global::get_n_MUs());
    PRINT_VAR(Global::get_n_pv_profiles());
    PRINT_VAR(Global::get_n_heatpump_profiles());
    PRINT_VAR(Global::get_ev_data_path());
    PRINT_VAR(Global::get_use_emission_time_series_ia());
    PRINT_VAR(Global::get_use_prices_time_series_ia());
    // Control strategy settings
    cout << "  Control strategy settings:\n";
    PRINT_ENUM_VAR(Global::get_controller_mode(), [](auto var){switch(var){case global::ControllerMode::RuleBased: return "RuleBased"; case global::ControllerMode::OptimizedWithPerfectForecast: return "OptimizedWithPerfectForecast"; default: return "";}});
    PRINT_VAR(Global::get_control_horizon_in_ts());
    PRINT_VAR(Global::get_control_update_freq_in_ts());
    PRINT_VAR(Global::get_controller_allow_bs_grid_charging);
    PRINT_ENUM_VAR(Global::get_controller_optimization_target(), [](auto var){switch(var){case global::ControllerOptimizationTarget::ElectricityCosts: return "ElectricityCosts"; case global::ControllerOptimizationTarget::PeakLoad: return "PeakLoad"; case global::ControllerOptimizationTarget::Emissions: return "Emissions"; default: return "";}});
    // Selection settings
    cout << "  Selection settings:\n";
    PRINT_ENUM_VAR(Global::get_exp_profile_mode(),  [](auto var){switch(var){case global::ExpansionProfileAllocationMode::Uninitialized: return "Uninitialized"; case global::ExpansionProfileAllocationMode::AsInData: return "AsInData"; case global::ExpansionProfileAllocationMode::Random: return "Random"; default: return "";}});
    PRINT_ENUM_VAR(Global::get_cu_selection_mode_fca(),   [](auto var){switch(var){case global::CUSModeFCA::Uninitialized: return "Uninitialized"; case global::CUSModeFCA::OrderAsInData: return "OrderAsInData"; case global::CUSModeFCA::RandomSelection: return "RandomSelection"; case global::CUSModeFCA::BestSSR: return "BestSSR"; case global::CUSModeFCA::BestNPV: return "BestNPV"; default: return "";}});
    PRINT_ENUM_VAR(Global::get_battery_power_computation_mode(), [](auto var){switch(var){case global::BatteryPowerComputationMode::AsDefinedByConfigVar: return "AsDefinedByConfigVar"; case global::BatteryPowerComputationMode::UseEOverPRatio: return "UseEOverPRatio"; default: return "";}});
    PRINT_ENUM_VAR(Global::get_battery_capacity_computation_mode(), [](auto var){switch(var){case global::BatteryCapacityComputationMode::Constant: return "Constant"; case global::BatteryCapacityComputationMode::BasedOnNominalPVPower: return "BasedOnNominalPVPower"; case global::BatteryCapacityComputationMode::BasedOnAnnualConsumption: return "BasedOnAnnualConsumption"; default: return "";}});
    PRINT_VAR(Global::get_annual_heat_demand_limit_fsac());
    PRINT_VAR(Global::get_select_buildings_wg_heatd_only());
    PRINT_VAR(Global::get_break_sac_loop_if_limit_reached());
    PRINT_VAR(Global::get_select_only_residential_buildings());
    // Scenario settings
    cout << "  Scenario settings:\n";
    PRINT_VAR(Global::get_exp_pv_static_mode());
    PRINT_VAR(Global::get_exp_pv_kWp_static());
    PRINT_VAR(Global::get_exp_pv_kWp_per_m2());
    PRINT_VAR(Global::get_exp_pv_min_kWp_roof_sec());
    PRINT_VAR(Global::get_exp_pv_max_kWp_roof_sec());
    PRINT_VAR(Global::get_exp_pv_max_kWp_per_unit());
    PRINT_VAR(Global::get_exp_pv_max_kWp_total());
    PRINT_VAR(Global::get_exp_pv_static_profile_orientation());
    PRINT_VAR(Global::get_exp_pv_static_profile_idx());
    PRINT_VAR(Global::get_exp_bess_kW());
    PRINT_VAR(Global::get_exp_bess_kWh());
    PRINT_VAR(Global::get_exp_bess_E_P_ratio());
    PRINT_VAR(Global::get_exp_bess_max_capacity());
    PRINT_VAR(Global::get_exp_bess_sizingE_boPV());
    PRINT_VAR(Global::get_exp_bess_effi_in());
    PRINT_VAR(Global::get_exp_bess_effi_out());
    PRINT_VAR(Global::get_exp_bess_self_ds_ts());
    PRINT_VAR(Global::get_exp_bess_start_soc());
    PRINT_VAR(Global::get_exp_bess_max_E_total());
    PRINT_VAR(Global::get_exp_bess_max_P_total());
    PRINT_VAR(Global::is_exp_hp_max_n_addition_set());
    PRINT_VAR(Global::get_exp_hp_max_n_addition());
    PRINT_VAR(Global::is_exp_ev_max_n_addition_set());
    PRINT_VAR(Global::get_exp_ev_max_n_addition());
    PRINT_VAR(Global::get_open_space_pv_kWp());
    PRINT_VAR(Global::get_wind_kWp());
    PRINT_VAR(Global::get_feed_in_tariff());
    PRINT_VAR(Global::get_demand_tariff());
    PRINT_VAR(Global::get_emissions_g_CO2eq_per_kWh());
    PRINT_VAR(Global::get_inst_cost_PV_per_kWp());
    PRINT_VAR(Global::get_inst_cost_BS_per_kWh());
    PRINT_VAR(Global::get_npv_discount_rate());
    PRINT_VAR(Global::get_npv_time_horizon());
    PRINT_VAR(Global::get_hp_flexibility_in_ts());
    PRINT_VAR(Global::get_heat_demand_thermalE_to_hpE_conv_f());
    PRINT_VAR(Global::get_heat_cons_bobv_slope());
    PRINT_VAR(Global::get_heat_cons_bobv_intercept());
    PRINT_VAR(Global::get_ev_plugin_probability());
    PRINT_VAR(Global::get_ev_battery_size_kWh());
    PRINT_VAR(Global::get_ev_consumption_kWh_km());
    PRINT_VAR(Global::get_ev_max_charging_power_kW());
    PRINT_VAR(Global::get_ev_charging_effi());
    // Output settings
    cout << "  Output settings:\n";
    PRINT_VAR(Global::get_compute_weekly_metrics());
    PRINT_ENUM_VAR(Global::get_output_mode_per_cu(), [](auto var){switch(var){case global::OutputModePerCU::IndividualFile: return "IndividualFile"; case global::OutputModePerCU::SingleFile: return "SingleFile"; case global::OutputModePerCU::NoOutput: return "NoOutput"; default: return "";}});
    PRINT_VAR(global::n_ts_between_flushs);
    PRINT_VAR(Global::get_create_substation_output());
    cout << global::output_section_delimiter << "\n";
}

