/*
 *
 * global.h
 *
 * Contains a namespace where all global variables are stored
 *
 * */

#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>

namespace global {

    inline int n_timesteps   = 0; ///< Total number of timesteps for which data is available
    inline int n_substations = 0; ///< Total number of substations for which data is available
    inline std::string* ts_start_str = NULL; ///< String of the start date
    inline std::string* ts_end_str   = NULL; ///< String of the end date
    inline int tsteps_per_hour      = 1;    ///< Time steps per hour in the simulation (and the data!)
    inline int expansion_scenario_id= 0;    ///< ID of the expansion scenario
    inline float exp_pv_kWp         = 0.0;  ///< kWp of in the simulation added PV installations
    inline float exp_bess_kW        = 0.0;  ///< P [kW] of in the simulation added BESS installations
    inline float exp_bess_kWh       = 0.0;  ///< E [kWh] of in the simulation added BESS installations
    inline float exp_bess_start_soc = 0.0;  ///< SOC at the beginning of the simulation for newly added BESS installations

    inline bool n_timesteps_init     = false;
    inline bool n_substations_init   = false;
    inline bool ts_start_str_init    = false;
    inline bool ts_end_str_init      = false;
    inline bool tsteps_per_hour_init = false;
    inline bool expansion_scenario_id_init = false;
    inline bool exp_pv_kWp_init      = false;
    inline bool exp_bess_kW_init     = false;
    inline bool exp_bess_kWh_init    = false;
    inline bool exp_bess_start_soc_init    = false;

    //
    // Checks if all variables are initialized
    //
    bool all_variables_initialized() {
        if (n_substations_init && 
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

}

#endif

