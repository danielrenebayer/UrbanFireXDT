/*
 *
 * global.h
 *
 * Contains a namespace where all global variables are stored
 *
 * */

#ifndef GLOBAL_H
#define GLOBAL_H

#include <ctime>
#include <list>
#include <string>
#include <utility>
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
    inline std::vector<struct tm*>* time_localtime_str     = NULL; ///< Reference to the list of the time as struct tm - alignment fits to time_timestep_id
    inline std::vector<std::string>* time_localtimezone_str = NULL; ///< Reference to the list of the time zone as string - alignment fits to time_timestep_id
    inline const float* pv_profile   = NULL; ///< Reference to the list of the global PV profile values
    inline const float* wind_profile = NULL; ///< Reference to the list of the global wind profile values
    inline OpenSpacePVOrWind* unit_open_space_pv   = NULL; ///< Reference to the global open space pv unit
    inline OpenSpacePVOrWind* unit_open_space_wind = NULL; ///< Reference to the global open space wind unit
    inline unsigned long n_ts_between_flushs = 1000; ///< Number of timesteps between the flush of the output buffers

    inline std::list<std::list<std::pair<string,float>>>* parameter_var_list = NULL; ///< List of parameters variation settings (i.e. the list contains a list of lists, where the inner lists represent a setting of ONE parameter variation setting (variable name, variable value))

    inline unsigned int curr_param_vari_combi_index = 0; ///< The index of the current parameter variation combination, that is simulated (0, if no parameter variation is selected)
    inline filesystem::path* current_output_dir  = NULL; ///< Reference to the object holding the current output path (maybe changed due to different parameter variations)
    inline filesystem::path* current_output_dir_prefix = NULL; ///< Reference to the path of the output path where all parameter variations can be found (i.e. the top level of current_output_dir if param vari is selected; i.e.2. one level below current_global_output_dir)
    inline filesystem::path* current_global_output_dir = NULL; ///< Reference to the object holding the current output dir for global information (i.e. information that does not change during parameter variations)

    inline bool time_info_init       = false;


    bool all_variables_initialized(); ///< Checks if all variables are initialized
    void vacuum();  ///< Deletes all global variables in the end


    /*
     * This enum defines different output modes per CU.
     * It corresponds to the --cu-output cmd line parameter.
     */
    enum struct OutputModePerCU : short {
        IndividualFile,
        SingleFile,
        NoOutput
    };


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
        static bool get_comp_eval_metrics()    { return comp_eval_metrics; }
        static bool is_parameter_variation()   { return pvar_selected;  }
        static int  get_parameter_varID()      { return pvar_id;        }
        static struct tm* get_ts_start_tm()    { return ts_start_tm;    }
        static struct tm* get_ts_end_tm()      { return ts_end_tm;    }
        static int get_tsteps_per_hour()       { return tsteps_per_hour;    }
        static int get_expansion_scenario_id() { return expansion_scenario_id;    }
        static float get_exp_pv_kWp()          { return exp_pv_kWp;    }
        static float get_exp_bess_kW()         { return exp_bess_kW;    }
        static float get_exp_bess_kWh()        { return exp_bess_kWh;    }
        static float get_exp_bess_start_soc()  { return exp_bess_start_soc;    }
        static float get_open_space_pv_kWp()   { return open_space_pv_kWp; }
        static float get_wind_kWp()            { return wind_kWp; }
        static const std::string& get_input_path()  { return input_path;  }
        static const std::string& get_output_path() { return output_path; }
        static global::OutputModePerCU get_output_mode_per_cu() { return output_mode_per_cu; }
        // setter methods
        static void set_n_timesteps(int n_timesteps);
        static void set_n_substations(int n_substations);
        static void set_n_CUs(int n_CUs);
        static void set_n_MUs(int n_MUs);
        static void set_comp_eval_metrics(bool value);
        static void set_pvar_vals(bool pvar_set, int pvarID);
        static void set_ts_start_tm(struct tm* ts_start_tm);
        static void set_ts_end_tm(struct tm* ts_end_tm);
        static void set_tsteps_per_hour(int tsteps_per_hour);
        static void set_expansion_scenario_id(int expansion_scenario_id);
        static void set_exp_pv_kWp(float exp_pv_kWp);
        static void set_exp_bess_kW(float exp_bess_kW);
        static void set_exp_bess_kWh(float exp_bess_kWh);
        static void set_exp_bess_start_soc(float exp_bess_start_soc);
        static void set_open_space_pv_kWp(float open_space_kWp);
        static void set_wind_kWp(float wind_kWp);
        static void set_input_path(std::string& path);
        static void set_output_path(std::string& path);
        static void set_output_mode_per_cu(global::OutputModePerCU mode);
    private:
        Global(); ///< Global cannot be initialized, it is a static only class
        // variables
        static int n_timesteps;            ///< Total number of timesteps for which data is available
        static int n_substations;          ///< Total number of substations for which data is available
        static int n_CUs;                  ///< Total number of control units for which data is available
        static int n_MUs;                  ///< Total number of meausrement units for which data is available
        static bool comp_eval_metrics;     ///< True, if evaluation metrics (like SSC,SSR) should be computed directly
        static bool pvar_selected;         ///< True, if a parameter variation is selected
        static int  pvar_id;               ///< ID of the parameter variation
        static struct tm* ts_start_tm;     ///< struct tm of the start date
        static struct tm* ts_end_tm;       ///< struct tm of the end date
        static int tsteps_per_hour;        ///< Time steps per hour in the simulation (and the data!)
        static int expansion_scenario_id;  ///< ID of the expansion scenario
        static float exp_pv_kWp;           ///< kWp of in the simulation added PV installations
        static float exp_bess_kW;          ///< P [kW] of in the simulation added BESS installations
        static float exp_bess_kWh;         ///< E [kWh] of in the simulation added BESS installations
        static float exp_bess_start_soc;   ///< SOC at the beginning of the simulation for newly added BESS installations
        static float open_space_pv_kWp;    ///< kWp of the open space PV installations (complete)
        static float wind_kWp;             ///< kWp of the wind turbines
        static std::string input_path;     ///< reference to the string holding the input path of the data
        static std::string output_path;    ///< reference to the string holding the output path of the data
        static global::OutputModePerCU output_mode_per_cu; ///< Variable storing the selected output mode per CU
        // boolean values holding information if the correspoding 
        // variable has been set or not
        static bool n_timesteps_init;
        static bool n_substations_init;
        static bool n_CUs_init;
        static bool n_MUs_init;
        static bool comp_eval_metrics_init;
        static bool pvar_set;
        static bool ts_start_str_init;
        static bool ts_end_str_init;
        static bool tsteps_per_hour_init;
        static bool expansion_scenario_id_init;
        static bool exp_pv_kWp_init;
        static bool exp_bess_kW_init;
        static bool exp_bess_kWh_init;
        static bool exp_bess_start_soc_init;
        static bool open_space_pv_kWp_init;
        static bool wind_kWp_init;
        static bool input_path_init;
        static bool output_path_init;
        static bool output_mode_per_cu_init;
};

#endif

