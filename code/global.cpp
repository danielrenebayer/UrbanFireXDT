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


