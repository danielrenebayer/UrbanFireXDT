#include "global.h"

using namespace global;
    
bool global::all_variables_initialized() {
    if (n_timesteps_init && 
	    n_substations_init &&
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

void global::vacuum() {
    delete   ts_start_str;           ts_start_str           = NULL;
    delete   ts_end_str;             ts_end_str             = NULL;
    delete[] time_timestep_id;       time_timestep_id       = NULL;
    delete   time_localtime_str;     time_localtime_str     = NULL;
    delete   time_localtimezone_str; time_localtimezone_str = NULL;
}


