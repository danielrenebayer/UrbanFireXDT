#include "global.h"

using namespace global;


#include <ctime>
#include <iostream>

using namespace std;



bool global::all_variables_initialized() {
    if (time_info_init &&
        pv_profiles_data != NULL &&
        hp_profiles  != NULL &&
        wind_profile != NULL &&
        unit_open_space_pv   != NULL &&
        unit_open_space_wind != NULL &&
        residual_gridload_kW != NULL)
    {
        return true;
    } else {
        return false;
    }
}

void global::vacuum() {
    // delete global arrays
    for (struct tm* t : *time_localtime_str) {
        delete t;
    }

    // delete heatpump profiles
    for (size_t hpIdx = 0; hpIdx < Global::get_n_heatpump_profiles(); hpIdx++)
        delete[] hp_profiles[hpIdx];
    delete[] hp_profiles; hp_profiles = NULL;
    // delete pv profiles
    for (size_t pvIdx = 0; pvIdx < Global::get_n_pv_profiles(); pvIdx++)
        delete[] pv_profiles_data[pvIdx];
    delete[] pv_profiles_data; pv_profiles_data = NULL;

    delete[] time_timestep_id;       time_timestep_id       = NULL;
    delete   time_localtime_str;     time_localtime_str     = NULL;
    delete   time_localtimezone_str; time_localtimezone_str = NULL;
    delete   unit_open_space_pv;     unit_open_space_pv   = NULL;
    delete   unit_open_space_wind;   unit_open_space_wind = NULL;
    delete   parameter_var_list;     parameter_var_list   = NULL;
    delete   current_output_dir;     current_output_dir   = NULL;
    delete   current_output_dir_prefix;current_output_dir_prefix   = NULL;
    delete   current_global_output_dir;current_global_output_dir   = NULL;
    delete[] wind_profile; wind_profile = NULL;
    delete[] residual_gridload_kW;   residual_gridload_kW = NULL;
}

void global::print_uninitialized_variables() {
    if (!time_info_init) {
        cout << "Variable time_info not initialized." << endl;
    }
    if (pv_profiles_data == NULL) {
        cout << "PV profiles (variable pv_profiles_data) undefined." << endl;
    }
    if (hp_profiles == NULL) {
        cout << "HP profiles (variable hp_profiles) undefined." << endl;
    }
    if (wind_profile == NULL) {
        cout << "Wind profiles (variable wind_profile) undefined." << endl;
    }
    if (unit_open_space_pv == NULL) {
        cout << "Open-space PV installations (variable unit_open_space_pv) undefined." << endl;
    }
    if (unit_open_space_wind == NULL) {
        cout << "Open-space wind farms (variable unit_open_space_wind) undefined." << endl;
    }
    if (residual_gridload_kW == NULL) {
        cout << "List of residual grid load in kW per time step is undefined." << endl;
    }
}





// ----------------------------- //
//      Implementation of        //
//            Global             //
// ----------------------------- //

unsigned long Global::n_timesteps = 0;
unsigned long Global::n_substations = 0;
unsigned long Global::n_CUs       = 0;
unsigned long Global::n_MUs       = 0;
unsigned long Global::n_pv_ts     = 0;
unsigned long Global::n_hp_ts     = 0;
//bool Global::comp_eval_metrics    = 0;
bool Global::pvar_selected        = false;
int  Global::pvar_id              = 0;
bool Global::repetitions_selected = false;
uint Global::n_repetitions        = 0;
struct tm* Global::ts_start_tm    = NULL;
struct tm* Global::ts_end_tm      = NULL;
int Global::tsteps_per_hour       = 1;
int Global::expansion_scenario_id = 0;
float Global::time_step_size_in_h   = 0.0;
bool  Global::exp_pv_kWp_static_mode  = false;
float Global::exp_pv_kWp_static       = 0.0;
float Global::exp_pv_kWp_per_m2       = 0.0;
float Global::exp_pv_min_kWp_roof_sec = 0.0;
float Global::exp_pv_max_kWp_per_sec  =-1.0;
float Global::exp_pv_max_kWp_per_unit =-1.0;
float Global::exp_pv_max_kWp_total    =-1.0;
float Global::exp_bess_kW           = 0.0;
float Global::exp_bess_kWh          = 0.0;
float Global::exp_bess_E_P_ratio    = 0.0;
float Global::exp_bess_start_soc    = 0.0;
float Global::open_space_pv_kWp     = 0.0;
float Global::wind_kWp              = 0.0;
float Global::feed_in_tariff        = 0.0;
float Global::demand_tariff         = 0.0;
float Global::inst_cost_PV_per_kWp  = 0.0;
float Global::inst_cost_BS_per_kWh  = 0.0;
float Global::npv_discount_rate     = 0.0;
unsigned int Global::npv_time_horizon  = 0;
bool  Global::use_BS_for_SSR_list = false;
string Global::input_path         = "";
string Global::output_path        = "";
OutputModePerCU Global::output_mode_per_cu = OutputModePerCU::IndividualFile;
ExpansionProfileAllocationMode Global::exp_profile_mode = ExpansionProfileAllocationMode::AsInData;
global::CUSModeFCA Global::cu_selection_mode_fca        = global::CUSModeFCA::OrderAsInData;
global::BatteryPowerComputationMode Global::bat_power_comp_mode = global::BatteryPowerComputationMode::AsDefinedByConfigVar;
//
bool Global::n_timesteps_init      = false;
bool Global::n_substations_init    = false;
bool Global::n_CUs_init            = false;
bool Global::n_MUs_init            = false;
bool Global::n_pv_ts_init          = false;
bool Global::n_hp_ts_init          = false;
//bool Global::comp_eval_metrics_init= false;
bool Global::pvar_set              = false;
bool Global::repetitions_selected_set     = false;
bool Global::n_repetitions_set     = false;
bool Global::ts_start_str_init     = false;
bool Global::ts_end_str_init       = false;
bool Global::tsteps_per_hour_init  = false;
bool Global::expansion_scenario_id_init = false;
bool Global::exp_pv_kWp_static_mode_init  = false;
bool Global::exp_pv_kWp_static_init       = false;
bool Global::exp_pv_kWp_per_m2_init       = false;
bool Global::exp_pv_min_kWp_roof_sec_init = false;
bool Global::exp_pv_max_kWp_per_sec_init  = false;
bool Global::exp_pv_max_kWp_per_unit_init = false;
bool Global::exp_pv_max_kWp_total_init    = false;
bool Global::exp_bess_kW_init      = false;
bool Global::exp_bess_kWh_init     = false;
bool Global::exp_bess_E_P_ratio_init    = false;
bool Global::exp_bess_start_soc_init    = false;
bool Global::open_space_pv_kWp_init= false;
bool Global::wind_kWp_init         = false;
bool Global::feed_in_tariff_set    = false;
bool Global::demand_tariff_set     = false;
bool Global::inst_cost_PV_per_kWp_set     = false;
bool Global::inst_cost_BS_per_kWh_set     = false;
bool Global::npv_discount_rate_set = false;
bool Global::npv_time_horizon_set  = false;
bool Global::use_BS_for_SSR_list_set      = false;
bool Global::input_path_init       = false;
bool Global::output_path_init      = false;
bool Global::output_mode_per_cu_init    = false;
bool Global::exp_profile_mode_init      = false;
bool Global::cu_selection_mode_fca_init = false;
bool Global::bat_power_comp_mode_init   = false;

void Global::InitializeStaticVariables() {
    // nothing to do anymore
}

void Global::DeleteStaticVariables() {
    delete ts_start_tm;
    delete ts_end_tm;
}

bool Global::AllVariablesInitialized() {
    if (n_timesteps_init && 
        n_substations_init &&
        n_CUs_init &&
        n_MUs_init &&
        n_pv_ts_init &&
        n_hp_ts_init &&
      //comp_eval_metrics_init &&
        pvar_set &&
        repetitions_selected_set &&
        ts_start_str_init &&
        ts_end_str_init &&
        tsteps_per_hour_init &&
        expansion_scenario_id_init &&
        exp_pv_kWp_static_mode_init &&
        exp_bess_kWh_init &&
        exp_bess_start_soc_init &&
        input_path_init &&
        output_path_init &&
        output_mode_per_cu_init &&
        exp_profile_mode_init &&
        cu_selection_mode_fca_init &&
        feed_in_tariff_set &&
        demand_tariff_set &&
        npv_discount_rate_set &&
        npv_time_horizon_set &&
        inst_cost_PV_per_kWp_set &&
        inst_cost_BS_per_kWh_set &&
        input_path_init &&
        bat_power_comp_mode_init)
    {
        if ((
             ( exp_pv_kWp_static_mode && exp_pv_kWp_static_init) ||
             (!exp_pv_kWp_static_mode && exp_pv_kWp_per_m2_init
                                      && exp_pv_min_kWp_roof_sec_init)
            ) && (
             ( bat_power_comp_mode == global::BatteryPowerComputationMode::AsDefinedByConfigVar && exp_bess_kW_init ) ||
             ( bat_power_comp_mode == global::BatteryPowerComputationMode::UseEOverPRatio       && exp_bess_E_P_ratio_init )
            ) && (
                !repetitions_selected || ( repetitions_selected && n_repetitions_set )
            )) {

                if (exp_pv_max_kWp_per_unit >= 0.0 && exp_pv_kWp_static_mode) {
                    cerr << "Warning: The combination of setting a maximal kWp per unit and using static mode does not make sens! The first will be ignored!" << endl;
                }

            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
        
}

void Global::PrintUninitializedVariables() {
    if (!n_timesteps_init) {
        cout << "Variable n_timesteps not initialized." << endl;
    }
    if (!n_substations_init) {
        cout << "Variable n_substations not initialized." << endl;
    }
    if (!n_CUs_init) {
        cout << "Variable n_CUs not initialized." << endl;
    }
    if (!n_MUs_init) {
        cout << "Variable n_MUs not initialized." << endl;
    }
    if (!n_pv_ts_init) {
        cout << "Variable n_pv_ts not initialized." << endl;
    }
    if (!n_hp_ts_init) {
        cout << "Variable n_hp_ts not initialized." << endl;
    }
    //if (!comp_eval_metrics_init) {
    //    cout << "Variable comp_eval_metrics not initialized." << endl;
    //}
    if (!pvar_set) {
        cout << "Variable pvar not initialized." << endl;
    }
    if (!repetitions_selected_set) {
        cout << "Variable repetitions_selected not initialized." << endl;
    }
    if (repetitions_selected && !n_repetitions_set) {
        cout << "Variable n_repetitions_set not initialized." << endl;
    }
    if (!ts_start_str_init) {
        cout << "Variable ts_start_str not initialized." << endl;
    }
    if (!ts_end_str_init) {
        cout << "Variable ts_end_str not initialized." << endl;
    }
    if (!tsteps_per_hour_init) {
        cout << "Variable tsteps_per_hour not initialized." << endl;
    }
    if (!expansion_scenario_id_init) {
        cout << "Variable expansion_scenario_id not initialized." << endl;
    }
    if (!exp_pv_kWp_static_mode_init) {
        cout << "Variable exp_pv_kWp_static_mode not initialized." << endl;
    }
    if (!exp_bess_kWh_init) {
        cout << "Variable exp_bess_kWh not initialized." << endl;
    }
    if (!exp_bess_start_soc_init) {
        cout << "Variable exp_bess_start_soc not initialized." << endl;
    }
    if (!input_path_init) {
        cout << "Variable input_path not initialized." << endl;
    }
    if (!output_path_init) {
        cout << "Variable output_path not initialized." << endl;
    }
    if (!output_mode_per_cu_init) {
        cout << "Variable output_mode_per_cu not initialized." << endl;
    }
    if (!exp_profile_mode_init) {
        cout << "Variable exp_profile_mode not initialized." << endl;
    }
    if (!cu_selection_mode_fca_init) {
        cout << "Variable cu_selection_mode_fca not initialized." << endl;
    }
    if (exp_pv_kWp_static_mode && !exp_pv_kWp_static_init) {
        cout << "Variable exp_pv_kWp_static not initialized." << endl;
    }
    if (!exp_pv_kWp_static_mode && !exp_pv_kWp_per_m2_init) {
        cout << "Variable exp_pv_kWp_per_m2 not initialized." << endl;
    }
    if (!exp_pv_kWp_static_mode && !exp_pv_min_kWp_roof_sec_init) {
        cout << "Variable exp_pv_min_kWp_roof_sec not initialized." << endl;
    }
    if (!exp_pv_kWp_static_mode && !exp_pv_max_kWp_per_sec_init) {
        cout << "Variable exp_pv_max_kWp_per_sec not initialized." << endl;
    }
    if (!feed_in_tariff_set) {
        cout << "Variable feed_in_tariff_set not initialized." << endl;
    }
    if (!demand_tariff_set) {
        cout << "Variable demand_tariff_set not initialized." << endl;
    }
    if (!npv_discount_rate_set) {
        cout << "Variable npv_discount_rate_set not initialized." << endl;
    }
    if (!npv_time_horizon_set) {
        cout << "Variable npv_time_horizon_set not initialized." << endl;
    }
    if (!inst_cost_PV_per_kWp_set) {
        cout << "Variable inst_cost_PV_per_kWp_set not initialized." << endl;
    }
    if (!inst_cost_BS_per_kWh_set) {
        cout << "Variable inst_cost_BS_per_kWh_set not initialized." << endl;
    }
    if (!bat_power_comp_mode_init) {
        cout << "Variable bat_power_comp_mode not initialized." << endl;
    }
    if (bat_power_comp_mode == global::BatteryPowerComputationMode::AsDefinedByConfigVar && !exp_bess_kW_init) {
        cout << "Variable exp_bess_kW not initialized." << endl;
    }
    if (bat_power_comp_mode == global::BatteryPowerComputationMode::UseEOverPRatio && !exp_bess_E_P_ratio_init) {
        cout << "Variable exp_bess_E_P_ratio not initialized." << endl;
    }
}

/*
inline int Global::get_n_timesteps() {
	return n_timesteps;
}
inline int Global::get_n_substations() {
	return n_substations;
}
inline int Global::get_n_CUs() {
	return n_CUs;
}
inline int Global::get_n_MUs() {
	return n_MUs;
}
inline std::string* Global::get_ts_start_str() {
	return ts_start_str;
}
inline std::string* Global::get_ts_end_str() {
	return ts_end_str;
}
inline int Global::get_tsteps_per_hour() {
	return tsteps_per_hour;
}
inline int Global::get_expansion_scenario_id() {
	return expansion_scenario_id;
}
inline float Global::get_exp_pv_kWp() {
	return exp_pv_kWp;
}
inline float Global::get_exp_bess_kW() {
	return exp_bess_kW;
}
inline float Global::get_exp_bess_kWh() {
	return exp_bess_kWh;
}
inline float Global::get_exp_bess_start_soc() {
	return exp_bess_start_soc;
}
*/

void Global::set_n_timesteps(unsigned long n_timesteps) {
    if (n_timesteps_init) {
        cerr << "Global variable n_timesteps is already initialized!" << endl;
    } else {
        Global::n_timesteps = n_timesteps;
        Global::n_timesteps_init = true;
    }
}
void Global::set_n_substations(unsigned long n_substations) {
    if (n_substations_init) {
        cerr << "Global variable n_substations is already initialized!" << endl;
    } else {
        Global::n_substations = n_substations;
        Global::n_substations_init = true;
    }
}
void Global::set_n_CUs(unsigned long n_CUs) {
    if (n_CUs_init) {
        cerr << "Global variable n_CUs is already initialized!" << endl;
    } else {
        Global::n_CUs = n_CUs;
        Global::n_CUs_init = true;
    }
}
void Global::set_n_pv_profiles(unsigned long value) {
    if (n_pv_ts_init) {
        cerr << "Global variable n_pv_profiles is already initialized!" << endl;
    } else {
        Global::n_pv_ts = value;
        Global::n_pv_ts_init = true;
    }
}
void Global::set_n_heatpump_profiles(unsigned long n_hp_ts) {
    if (n_hp_ts_init) {
        cerr << "Global variable n_heatpump_profiles is already initialized!" << endl;
    } else {
        Global::n_hp_ts = n_hp_ts;
        Global::n_hp_ts_init = true;
    }
}
void Global::set_n_MUs(unsigned long n_MUs) {
    if (n_MUs_init) {
        cerr << "Global variable n_MUs is already initialized!" << endl;
    } else {
        Global::n_MUs = n_MUs;
        Global::n_MUs_init = true;
    }
}
/*void Global::set_comp_eval_metrics(bool value) {
    if (comp_eval_metrics_init) {
        cerr << "Global variable comp_eval_metrics is already initialized!" << endl;
    } else {
        Global::comp_eval_metrics = value;
        Global::comp_eval_metrics_init = true;
    }
}*/
void Global::set_pvar_vals(bool pvar_val, int pvarID) {
    if (pvar_set) {
        cerr << "Values for parameter variation are already set!" << endl;
    } else {
        Global::pvar_selected = pvar_val;
        Global::pvar_id  = pvarID;
        Global::pvar_set = true;
    }
}
void Global::set_repetitions_selected(bool value) {
    if (repetitions_selected_set){
        cerr << "Global variable repetitions_selected is already initialized!" << endl;
    } else {
        Global::repetitions_selected = value;
        Global::repetitions_selected_set = true;
    }
}
void Global::set_n_repetitions(unsigned int value) {
    if (n_repetitions_set){
        cerr << "Global variable n_repetitions is already initialized!" << endl;
    } else {
        Global::n_repetitions = value;
        Global::n_repetitions_set = true;
    }
}
void Global::set_ts_start_tm(struct tm* ts_start_tm) {
    if (ts_start_str_init) {
        cerr << "Global variable ts_start_str is already initialized!" << endl;
    } else {
        Global::ts_start_tm = ts_start_tm;
        Global::ts_start_str_init = true;
    }
}
void Global::set_ts_end_tm(struct tm* ts_end_tm) {
    if (ts_end_str_init) {
        cerr << "Global variable ts_end_str is already initialized!" << endl;
    } else {
        Global::ts_end_tm = ts_end_tm;
        Global::ts_end_str_init = true;
    }
}
void Global::set_tsteps_per_hour(int tsteps_per_hour) {
    if (tsteps_per_hour_init) {
        cerr << "Global variable tsteps_per_hour is already initialized!" << endl;
    } else {
        if (tsteps_per_hour < 0) {
            cerr << "tsteps_per_hour should be set to a value < 0. This is not allowed." << endl;
            return;
        }
        Global::tsteps_per_hour = tsteps_per_hour;
        Global::tsteps_per_hour_init = true;
        // set set_tsteps_per_hour accordingly
        Global::time_step_size_in_h = 1.0f / (float)tsteps_per_hour;
        #ifdef DEBUG
        cout << "tsteps_per_hour = " << tsteps_per_hour << endl;
        cout << "Global::time_step_size_in_h = " << Global::time_step_size_in_h << endl;
        #endif
    }
}
void Global::set_expansion_scenario_id(int expansion_scenario_id) {
    if (expansion_scenario_id_init) {
        cerr << "Global variable expansion_scenario_id is already initialized!" << endl;
    } else {
        Global::expansion_scenario_id = expansion_scenario_id;
        Global::expansion_scenario_id_init = true;
    }
}
void Global::set_exp_pv_mode(bool mode) {
    if (exp_pv_kWp_static_mode_init) {
        cerr << "Global variable exp_pv_kWp_static_mode is already initialized!" << endl;
    } else {
        Global::exp_pv_kWp_static_mode      = mode;
        Global::exp_pv_kWp_static_mode_init = true;
    }
}
void Global::set_exp_pv_kWp_static(float value) {
    if (exp_pv_kWp_static_init) {
        cerr << "Global variable exp_pv_kWp_static is already initialized!" << endl;
    } else {
        Global::exp_pv_kWp_static      = value;
        Global::exp_pv_kWp_static_init = true;
    }
}
void Global::set_exp_pv_kWp_per_m2(float value) {
    if (exp_pv_kWp_per_m2_init) {
        cerr << "Global variable exp_pv_kWp_per_m2 is already initialized!" << endl;
    } else {
        Global::exp_pv_kWp_per_m2      = value;
        Global::exp_pv_kWp_per_m2_init = true;
    }
}
void Global::set_exp_pv_min_kWp_roof_sec(float value) {
    if (exp_pv_min_kWp_roof_sec_init) {
        cerr << "Global variable exp_pv_min_kWp_roof_sec is already initialized!" << endl;
    } else {
        Global::exp_pv_min_kWp_roof_sec      = value;
        Global::exp_pv_min_kWp_roof_sec_init = true;
    }
}
void Global::set_exp_pv_max_kWp_roof_sec(float value) {
    if (exp_pv_max_kWp_per_sec_init) {
        cerr << "Global variable exp_pv_max_kWp_per_sec is already initialized!" << endl;
    } else {
        Global::exp_pv_max_kWp_per_sec      = value;
        Global::exp_pv_max_kWp_per_sec_init = true;
    }
}
void Global::set_exp_pv_max_kWp_per_unit(float value) {
    if (Global::exp_pv_max_kWp_per_unit_init) {
        cerr << "Global variable exp_pv_max_kWp_per_unit is already initialized!" << endl;
    } else {
        Global::exp_pv_max_kWp_per_unit      = value;
        Global::exp_pv_max_kWp_per_unit_init = true;
    }
}
void Global::set_exp_pv_max_kWp_total(float value) {
    if (Global::exp_pv_max_kWp_total_init) {
        cerr << "Global variable exp_pv_max_kWp_total is already initialized!" << endl;
    } else {
        Global::exp_pv_max_kWp_total      = value;
        Global::exp_pv_max_kWp_total_init = true;
    }
}
void Global::set_exp_bess_kW(float exp_bess_kW) {
    if (exp_bess_kW_init) {
        cerr << "Global variable exp_bess_kW is already initialized!" << endl;
    } else {
        Global::exp_bess_kW = exp_bess_kW;
        Global::exp_bess_kW_init = true;
    }
}
void Global::set_exp_bess_kWh(float exp_bess_kWh) {
    if (exp_bess_kWh_init) {
        cerr << "Global variable exp_bess_kWh is already initialized!" << endl;
    } else {
        Global::exp_bess_kWh = exp_bess_kWh;
        Global::exp_bess_kWh_init = true;
    }
}
void Global::set_exp_bess_E_P_ratio(float value) {
    if (exp_bess_E_P_ratio_init) {
        cerr << "Global variable exp_bess_E_P_ratio is already initialized!" << endl;
    } else {
        Global::exp_bess_E_P_ratio = value;
        Global::exp_bess_E_P_ratio_init = true;
    }
}
void Global::set_exp_bess_start_soc(float exp_bess_start_soc) {
    if (exp_bess_start_soc_init) {
        cerr << "Global variable exp_bess_start_soc is already initialized!" << endl;
    } else {
        Global::exp_bess_start_soc = exp_bess_start_soc;
        Global::exp_bess_start_soc_init = true;
    }
}
void Global::set_open_space_pv_kWp(float open_space_kWp) {
    if (open_space_pv_kWp_init) {
        cerr << "Global variable open_space_pv_kWp is already initialized!" << endl;
    } else {
        Global::open_space_pv_kWp = open_space_kWp;
        open_space_pv_kWp_init = true;
    }
}
void Global::set_wind_kWp(float wind_kWp) {
    if (wind_kWp_init) {
        cerr << "Global variable wind_kWp is already initialized!" << endl;
    } else {
        Global::wind_kWp = wind_kWp;
        wind_kWp_init = true;
    }
}
void Global::set_feed_in_tariff(float value) {
    if (feed_in_tariff_set) {
        cerr << "Global variable feed_in_tariff is already initialized!" << endl;
    } else {
        Global::feed_in_tariff = value;
        feed_in_tariff_set = true;
    }
}
void Global::set_demand_tariff(float value) {
    if (demand_tariff_set) {
        cerr << "Global variable demand_tariff is already initialized!" << endl;
    } else {
        Global::demand_tariff = value;
        demand_tariff_set = true;
    }
}
void Global::set_npv_discount_rate(float value) {
    if (npv_discount_rate_set) {
        cerr << "Global variable npv_discount_rate is already initialized!" << endl;
    } else {
        Global::npv_discount_rate = value;
        npv_discount_rate_set = true;
    }
}
void Global::set_npv_time_horizon(unsigned int value) {
    if (npv_time_horizon_set) {
        cerr << "Global variable npv_time_horizon is already initialized!" << endl;
    } else {
        Global::npv_time_horizon = value;
        npv_time_horizon_set = true;
    }
}
void Global::set_inst_cost_PV_per_kWp(float value) {
    if (inst_cost_PV_per_kWp_set) {
        cerr << "Global variable inst_cost_PV_per_kWp is already initialized!" << endl;
    } else {
        Global::inst_cost_PV_per_kWp = value;
        inst_cost_PV_per_kWp_set = true;
    }
}
void Global::set_inst_cost_BS_per_kWh(float value) {
    if (inst_cost_BS_per_kWh_set) {
        cerr << "Global variable inst_cost_BS_per_kWh is already initialized!" << endl;
    } else {
        Global::inst_cost_BS_per_kWh = value;
        inst_cost_BS_per_kWh_set = true;
    }
}
void Global::set_use_BS_for_SSR_list(bool mode) {
    if (use_BS_for_SSR_list_set) {
        cerr << "Global variable use_BS_for_SSR_list is already initialized!" << endl;
    } else {
        Global::use_BS_for_SSR_list = mode;
        use_BS_for_SSR_list_set = true;
    }
}
void Global::set_input_path(string& path) {
    if (input_path_init) {
        cerr << "Input path already set!" << endl;
    } else {
        Global::input_path = path;
        Global::input_path_init = true;
    }
}
void Global::set_output_path(string& path) {
    if (output_path_init) {
        cerr << "Output path already set!" << endl;
    } else {
        Global::output_path = path;
        Global::output_path_init = true;
    }
}
void Global::set_output_mode_per_cu(OutputModePerCU mode) {
    if (output_mode_per_cu_init) {
        cerr << "Output mode already set!" << endl;
    } else {
        Global::output_mode_per_cu = mode;
        Global::output_mode_per_cu_init = true;
    }
}
void Global::set_exp_profile_mode(ExpansionProfileAllocationMode mode) {
    if (exp_profile_mode_init) {
        cerr << "Expansion profile mode already set!" << endl;
    } else {
        Global::exp_profile_mode = mode;
        Global::exp_profile_mode_init = true;
    }
}
void Global::set_cu_selection_mode_fca(global::CUSModeFCA mode) {
    if (cu_selection_mode_fca_init) {
        cerr << "CU selection mode for component addition already set!" << endl;
    } else {
        Global::cu_selection_mode_fca = mode;
        Global::cu_selection_mode_fca_init = true;
    }
}
void Global::set_battery_power_computation_mode(global::BatteryPowerComputationMode mode) {
    if (bat_power_comp_mode_init) {
        cerr << "Battery power computation mode is already set!" << endl;
    } else {
        Global::bat_power_comp_mode = mode;
        Global::bat_power_comp_mode_init = true;
    }
}

