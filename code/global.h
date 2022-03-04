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
#include <vector>

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
    // Time info (read from time_indices table in central values database)
    inline int* time_timestep_id = NULL; ///< Reference to the list of time steps
    inline std::vector<std::string>* time_localtime_str     = NULL; ///< Reference to the list of the time as string - alignment fits to time_timestep_id
    inline std::vector<std::string>* time_localtimezone_str = NULL; ///< Reference to the list of the time zone as string - alignment fits to time_timestep_id

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
    inline bool time_info_init       = false;


    bool all_variables_initialized(); ///< Checks if all variables are initialized
    void vacuum();  ///< Deletes all global variables in the end

}

#endif

