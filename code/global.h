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

#include "units.h"

namespace global {
    /*
     * Namespace global
     *
     * It contains all global attributes and variables that
     * might change during simulation execution.
     * 
     * Attention: There is no access protection for these variables!
     * For access protection use class Global.
     * 
     * Attention: Do not confuse with class Global (mind the capital "G")!
     */

    // inline float current_energy_feedin_price = 0.0;
    // inline float current_energy_supply_price = 0.0;
    // Time info (read from time_indices table in central values database)
    inline int* time_timestep_id = NULL; ///< Reference to the list of time steps
    inline std::vector<std::string>* time_localtime_str     = NULL; ///< Reference to the list of the time as string - alignment fits to time_timestep_id
    inline std::vector<std::string>* time_localtimezone_str = NULL; ///< Reference to the list of the time zone as string - alignment fits to time_timestep_id

    inline bool time_info_init       = false;


    bool all_variables_initialized(); ///< Checks if all variables are initialized
    void vacuum();  ///< Deletes all global variables in the end

}

class Global {
    /*
     * class Global
     *
     * This class contains all global variables that cannot change
     * after they have been set once.
     * 
     * Attention: Not to be confused with namespace global (mind the lower case "g").
     */
    public:
        static void InitializeStaticVariables();
        static void DeleteStaticVariables();
        //
        static bool AllVariablesInitialized();
        //
        // getter methods
        /*
        static inline int get_n_timesteps();
        static inline int get_n_substations();
        static inline int get_n_CUs();
        static inline int get_n_MUs();
        static inline std::string* get_ts_start_str();
        static inline std::string* get_ts_end_str();
        static inline int get_tsteps_per_hour();
        static inline int get_expansion_scenario_id();
        static inline float get_exp_pv_kWp();
        static inline float get_exp_bess_kW();
        static inline float get_exp_bess_kWh();
        static inline float get_exp_bess_start_soc();
        */
        static int get_n_timesteps()   {  return n_timesteps; }
        static int get_n_substations() {  return n_substations; }
        static int get_n_CUs()         {  return n_CUs; }
        static int get_n_MUs()         {  return n_MUs; }
        static std::string* get_ts_start_str() { return ts_start_str;    }
        static std::string* get_ts_end_str()   { return ts_end_str;    }
        static int get_tsteps_per_hour()       { return tsteps_per_hour;    }
        static int get_expansion_scenario_id() { return expansion_scenario_id;    }
        static float get_exp_pv_kWp()          { return exp_pv_kWp;    }
        static float get_exp_bess_kW()         { return exp_bess_kW;    }
        static float get_exp_bess_kWh()        { return exp_bess_kWh;    }
        static float get_exp_bess_start_soc()  { return exp_bess_start_soc;    }
        // setter methods
        static void set_n_timesteps(int n_timesteps);
        static void set_n_substations(int n_substations);
        static void set_n_CUs(int n_CUs);
        static void set_n_MUs(int n_MUs);
        static void set_ts_start_str(std::string* ts_start_str);
        static void set_ts_end_str(std::string* ts_end_str);
        static void set_tsteps_per_hour(int tsteps_per_hour);
        static void set_expansion_scenario_id(int expansion_scenario_id);
        static void set_exp_pv_kWp(float exp_pv_kWp);
        static void set_exp_bess_kW(float exp_bess_kW);
        static void set_exp_bess_kWh(float exp_bess_kWh);
        static void set_exp_bess_start_soc(float exp_bess_start_soc);
    private:
        Global(); ///< Global cannot be initialized, it is a static only class
        // variables
        static int n_timesteps;            ///< Total number of timesteps for which data is available
        static int n_substations;          ///< Total number of substations for which data is available
        static int n_CUs;                  ///< Total number of control units for which data is available
        static int n_MUs;                  ///< Total number of meausrement units for which data is available
        static std::string* ts_start_str;  ///< String of the start date
        static std::string* ts_end_str;    ///< String of the end date
        static int tsteps_per_hour;        ///< Time steps per hour in the simulation (and the data!)
        static int expansion_scenario_id;  ///< ID of the expansion scenario
        static float exp_pv_kWp;           ///< kWp of in the simulation added PV installations
        static float exp_bess_kW;          ///< P [kW] of in the simulation added BESS installations
        static float exp_bess_kWh;         ///< E [kWh] of in the simulation added BESS installations
        static float exp_bess_start_soc;   ///< SOC at the beginning of the simulation for newly added BESS installations
        // boolean values holding information if the correspoding 
        // variable has been set or not
        static bool n_timesteps_init;
        static bool n_substations_init;
        static bool n_CUs_init;
        static bool n_MUs_init;
        static bool ts_start_str_init;
        static bool ts_end_str_init;
        static bool tsteps_per_hour_init;
        static bool expansion_scenario_id_init;
        static bool exp_pv_kWp_init;
        static bool exp_bess_kW_init;
        static bool exp_bess_kWh_init;
        static bool exp_bess_start_soc_init;
};

#endif

