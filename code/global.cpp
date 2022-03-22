#include "global.h"

using namespace global;


#include <ctime>
#include <iostream>

using namespace std;



bool global::all_variables_initialized() {
    if (time_info_init && substation_output_init)
    {
        return true;
    } else {
        return false;
    }
}

void global::vacuum() {
    // close output file for substations
    if (substation_output_init) {
        substation_output->close();
        delete substation_output;
    }

    // delete global arrays
    for (struct tm* t : *time_localtime_str) {
        delete t;
    }

    delete[] time_timestep_id;       time_timestep_id       = NULL;
    delete   time_localtime_str;     time_localtime_str     = NULL;
    delete   time_localtimezone_str; time_localtimezone_str = NULL;
}





// ----------------------------- //
//      Implementation of        //
//            Global             //
// ----------------------------- //

int Global::n_timesteps           = 0;
int Global::n_substations         = 0;
int Global::n_CUs                 = 0;
int Global::n_MUs                 = 0;
struct tm* Global::ts_start_tm    = NULL;
struct tm* Global::ts_end_tm      = NULL;
int Global::tsteps_per_hour       = 1;
int Global::expansion_scenario_id = 0;
float Global::exp_pv_kWp            = 0.0;
float Global::exp_bess_kW           = 0.0;
float Global::exp_bess_kWh          = 0.0;
float Global::exp_bess_start_soc    = 0.0;
//
bool Global::n_timesteps_init      = false;
bool Global::n_substations_init    = false;
bool Global::n_CUs_init            = false;
bool Global::n_MUs_init            = false;
bool Global::ts_start_str_init     = false;
bool Global::ts_end_str_init       = false;
bool Global::tsteps_per_hour_init  = false;
bool Global::expansion_scenario_id_init = false;
bool Global::exp_pv_kWp_init       = false;
bool Global::exp_bess_kW_init      = false;
bool Global::exp_bess_kWh_init     = false;
bool Global::exp_bess_start_soc_init    = false;

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
		ts_start_str_init &&
	    ts_end_str_init &&
	    tsteps_per_hour_init &&
		expansion_scenario_id_init &&
        exp_pv_kWp_init &&
        exp_bess_kW_init &&
        exp_bess_kWh_init &&
        exp_bess_start_soc_init)
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
        Global::tsteps_per_hour = tsteps_per_hour;
        Global::tsteps_per_hour_init = true;
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
