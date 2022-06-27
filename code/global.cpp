#include "global.h"

using namespace global;


#include <ctime>
#include <iostream>

using namespace std;



bool global::all_variables_initialized() {
    if (time_info_init &&
        pv_profile   != NULL &&
        wind_profile != NULL &&
        unit_open_space_pv   != NULL &&
        unit_open_space_wind != NULL)
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

    delete[] time_timestep_id;       time_timestep_id       = NULL;
    delete   time_localtime_str;     time_localtime_str     = NULL;
    delete   time_localtimezone_str; time_localtimezone_str = NULL;
    delete   unit_open_space_pv;     unit_open_space_pv   = NULL;
    delete   unit_open_space_wind;   unit_open_space_wind = NULL;
    delete   parameter_var_list;     parameter_var_list   = NULL;
    delete   current_output_dir;     current_output_dir   = NULL;
    delete   current_output_dir_prefix;current_output_dir_prefix   = NULL;
    delete   current_global_output_dir;current_global_output_dir   = NULL;
    delete[] pv_profile;   pv_profile   = NULL;
    delete[] wind_profile; wind_profile = NULL;
}





// ----------------------------- //
//      Implementation of        //
//            Global             //
// ----------------------------- //

int Global::n_timesteps           = 0;
int Global::n_substations         = 0;
int Global::n_CUs                 = 0;
int Global::n_MUs                 = 0;
bool Global::comp_eval_metrics    = 0;
bool Global::pvar_selected        = false;
int  Global::pvar_id              = 0;
struct tm* Global::ts_start_tm    = NULL;
struct tm* Global::ts_end_tm      = NULL;
int Global::tsteps_per_hour       = 1;
int Global::expansion_scenario_id = 0;
float Global::time_step_size_in_h   = 0.0;
float Global::exp_pv_kWp            = 0.0;
float Global::exp_bess_kW           = 0.0;
float Global::exp_bess_kWh          = 0.0;
float Global::exp_bess_start_soc    = 0.0;
float Global::open_space_pv_kWp     = 0.0;
float Global::wind_kWp              = 0.0;
string Global::input_path         = "";
string Global::output_path        = "";
OutputModePerCU Global::output_mode_per_cu = OutputModePerCU::IndividualFile;
//
bool Global::n_timesteps_init      = false;
bool Global::n_substations_init    = false;
bool Global::n_CUs_init            = false;
bool Global::n_MUs_init            = false;
bool Global::comp_eval_metrics_init= false;
bool Global::pvar_set              = false;
bool Global::ts_start_str_init     = false;
bool Global::ts_end_str_init       = false;
bool Global::tsteps_per_hour_init  = false;
bool Global::expansion_scenario_id_init = false;
bool Global::exp_pv_kWp_init       = false;
bool Global::exp_bess_kW_init      = false;
bool Global::exp_bess_kWh_init     = false;
bool Global::exp_bess_start_soc_init    = false;
bool Global::open_space_pv_kWp_init= false;
bool Global::wind_kWp_init         = false;
bool Global::input_path_init       = false;
bool Global::output_path_init      = false;
bool Global::output_mode_per_cu_init    = false;

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
        comp_eval_metrics_init &&
        pvar_set &&
        ts_start_str_init &&
        ts_end_str_init &&
        tsteps_per_hour_init &&
        expansion_scenario_id_init &&
        exp_pv_kWp_init &&
        exp_bess_kW_init &&
        exp_bess_kWh_init &&
        exp_bess_start_soc_init &&
        input_path_init &&
        output_path_init &&
        output_mode_per_cu_init)
    {
        return true;
    } else {
        return false;
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

void Global::set_n_timesteps(int n_timesteps) {
    if (n_timesteps_init) {
        cerr << "Global variable n_timesteps is already initialized!" << endl;
    } else {
        Global::n_timesteps = n_timesteps;
        Global::n_timesteps_init = true;
    }
}
void Global::set_n_substations(int n_substations) {
    if (n_substations_init) {
        cerr << "Global variable n_substations is already initialized!" << endl;
    } else {
        Global::n_substations = n_substations;
        Global::n_substations_init = true;
    }
}
void Global::set_n_CUs(int n_CUs) {
    if (n_CUs_init) {
        cerr << "Global variable n_CUs is already initialized!" << endl;
    } else {
        Global::n_CUs = n_CUs;
        Global::n_CUs_init = true;
    }
}
void Global::set_n_MUs(int n_MUs) {
    if (n_MUs_init) {
        cerr << "Global variable n_MUs is already initialized!" << endl;
    } else {
        Global::n_MUs = n_MUs;
        Global::n_MUs_init = true;
    }
}
void Global::set_comp_eval_metrics(bool value) {
    if (comp_eval_metrics_init) {
        cerr << "Global variable comp_eval_metrics is already initialized!" << endl;
    } else {
        Global::comp_eval_metrics = value;
        Global::comp_eval_metrics_init = true;
    }
}
void Global::set_pvar_vals(bool pvar_val, int pvarID) {
    if (pvar_set) {
        cerr << "Values for parameter variation are already set!" << endl;
    } else {
        Global::pvar_selected = pvar_val;
        Global::pvar_id  = pvarID;
        Global::pvar_set = true;
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
        Global::time_step_size_in_h = 1.0/tsteps_per_hour;
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
void Global::set_exp_pv_kWp(float exp_pv_kWp) {
    if (exp_pv_kWp_init) {
        cerr << "Global variable exp_pv_kWp is already initialized!" << endl;
    } else {
        Global::exp_pv_kWp = exp_pv_kWp;
        Global::exp_pv_kWp_init = true;
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

