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
#include <filesystem>
#include <list>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "sac_planning.h"
#include "units.h"

/*!
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
namespace global {

    // inline float current_energy_feedin_price = 0.0;
    // inline float current_energy_supply_price = 0.0;
    // Time info (read from time_indices table in central values database)
    inline unsigned long* time_timestep_id = NULL; ///< Reference to the list of time steps
    inline std::vector<struct tm*>* time_localtime_r        = NULL; ///< Reference to the list of right-aligned time stamps as struct tm. Right-aligned means the time at time step ts denotes the time at the ending of the interval - alignment fits to time_timestep_id
    inline std::vector<struct tm*>* time_localtime_l        = NULL; ///< Reference to the list of left-aligned time stamps as struct tm. Left-aligned means the time at time step ts denotes the time at the start of the interval - alignment fits to time_timestep_id
    inline std::vector<std::string>* time_localtimezone_str = NULL; ///< Reference to the list of the time zone as string - alignment fits to time_timestep_id
    inline std::map<std::string, size_t> pv_profiles_information; ///< Map (orientation, number of time series) where the number of pv-profiles per orientation are given
    inline std::map<std::string, std::vector<const float*>> pv_profiles_per_ori; ///< Map (orientation, vector of references to the time series) where the available time series per orientation are given - same ordering as in data
    inline const float* const* pv_profiles_data = NULL; ///< Reference to the list of global PV profiles (Two dimensional array with [n_pv_profiles, n_timesteps])
    inline const float* const* hp_profiles = NULL; ///< Reference to the list of global heat pump profiles (Two dimensional array with [n_heatpump_profiles, n_timesteps])
    inline const double* const* hp_profiles_cumsum = NULL; ///< Reference to the list of cumulated global heat pump profiles (Two dimensional array with [n_heatpump_profiles, n_timesteps])
    inline const float* residual_gridload_kW = NULL; ///< Residual netload, i.e. amount of load that has to be added to the final netload, this is a load that is not measured by smart meters occuring in this simulation
    inline const float* wind_profile = NULL; ///< Reference to the list of the global wind profile values
    inline OpenSpacePVOrWind* unit_open_space_pv   = NULL; ///< Reference to the global open space pv unit
    inline OpenSpacePVOrWind* unit_open_space_wind = NULL; ///< Reference to the global open space wind unit
    inline unsigned long n_ts_between_flushs = 1000; ///< Number of timesteps between the flush of the output buffers
    inline std::map<unsigned long, float> annual_heat_demand_kWh; ///< Map storing the annaul heat demand for the buildings in kWh (thermal) [per location id]
    inline std::map<unsigned long, float> building_volumes_m3;    ///< Map storing the volume of the biggest building at a given location in cubic-meters (m^3)
    inline std::map<unsigned long, std::vector<std::pair<float, std::string>>> roof_section_orientations; ///< Map storing a complet list of roof sections per location ID. Roof sections are tuples/pairs with the information (roof section area, roof section orientation)
    inline std::set<unsigned long> locations_with_geodata; ///< List / Set storing all locations for which geodata is available
    inline unsigned int current_repetition_counter = 0; ///< Repetition counter, if repetitions is set as a command line argument
    inline const float* emission_ts      = NULL; ///< Reference to the emission time series for grid demanded energy (if present) in g CO2eq per kWh
    inline const float* eprices_local_ts = NULL; ///< Reference to the time series of the energy prices (for residential customers) for grid demanded energy (if present)
    inline const float* eprices_spotm_ts = NULL; ///< Reference to the time series of the energy prices at the spot market (for non-residential customers, without vat or taxes) (if present)
    inline std::vector<unsigned long> unitIDs_selected_for_output; ///< A vector containing a list of all unitIDs selected for output. If this vector is empty (as it is by default), all units will generate an output.

    inline std::list<std::list<std::pair<std::string,float>>>* parameter_var_list = NULL; ///< List of parameters variation settings (i.e. the list contains a list of lists, where the inner lists represent a setting of ONE parameter variation setting (variable name, variable value))

    inline unsigned int curr_param_vari_combi_index = 0; ///< The index of the current parameter variation combination, that is simulated (0, if no parameter variation is selected)
    inline std::filesystem::path* current_output_dir  = NULL; ///< Reference to the object holding the current output path (maybe changed due to different parameter variations)
    inline std::filesystem::path* current_output_dir_prefix = NULL; ///< Reference to the path of the output path where all parameter variations can be found (i.e. the top level of current_output_dir if param vari is selected; i.e.2. one level below current_global_output_dir)
    inline std::filesystem::path* current_global_output_dir = NULL; ///< Reference to the object holding the current output dir for global information (i.e. information that does not change during parameter variations)

    inline bool time_info_init       = false;

    inline std::chrono::time_point<std::chrono::system_clock> time_of_simulation_start; ///< The time of the simulation start

    inline std::string structure_db_column_name_energy_prices = "local_price"; ///< The name of the column containing the local energy prices in EUR/kWh


    bool all_variables_initialized(); ///< Checks if all variables are initialized
    void vacuum();  ///< Deletes all global variables in the end
    void print_uninitialized_variables(); ///< Prints all variable names to stdout, that are not initialized


    /*!
     * This enum defines different output modes per CU.
     * It corresponds to the --cu-output cmd line parameter.
     */
    enum struct OutputModePerCU : short {
        IndividualFile,
        SingleFile,
        NoOutput
    };

    /*!
     * This enum defines the different expansion profile allocation modes.
     */
    enum struct ExpansionProfileAllocationMode : short {
        Uninitialized,
        AsInData,
        Random
    };

    /*!
     * This struct defines the
     * Control Unit Selection Mode For Component Addition
     * i.e. it defines how the control units are selected, that get simulatively added components (like PV, BS, HP, ...)
     */
    enum struct CUSModeFCA {
        Uninitialized,
        OrderAsInData,
        RandomSelection,
        BestSSR,
        BestNPV,
        UseList
    };

    /*!
     * This enum defines how the battery storage power
     * should be computet.
     * There are two options:
     * 1. by setting a given value, i.e. config-variable 'expansion BS P in kW' will be used.
     * 2. by setting the E/P-ratio of the battery, i.e. battery P = E/P-ratio * selected E
     */
    enum struct BatteryPowerComputationMode {
        AsDefinedByConfigVar,
        UseEOverPRatio
    };

    /*!
     * This enum defines how the battery capacity sizing mode.
     * There are the following options:
     * 1. by using the constant value given in config-variable 'expansion BS E in kWh'
     * 2. by using the following formula (only, if a simulated PV is present): E_BS = P_PV * config-variable 'expansion BS capacity sizing factor for PV'
     * 3. by using the annual electricity consumption - only of the connected measurement units
     * 4. by using the annual electricity consumption - only of the connected measurement units + the heat pump (but not the EVs!)
     * 5. determine the optimal size based on a optimization with a perfect forecast
     */
    enum struct BatteryCapacityComputationMode {
        Constant,
        BasedOnNominalPVPower,
        BasedOnAnnualConsumption,
        BasedOnAnnualConsumptionWithHeatPump,
        Optimized,
    };

    /**
     * This enum defines the different PV sizing modes. This is the analogon to the enum struct `BatteryCapacityComputationMode`.
     * There are the following options:
     * 1. use the maximum available roof area
     * 2. use a constant, i.e., static, PV size for all simulatively added PV installations with the same (configurable) orientation
     * 3. determine the optimal size based on a optimization with a perfect forecast
     */
    enum struct PVSizingMode : short {
        MaxAvailableRoofArea,
        StaticPVSize,
        Optimized,
    };

    /**
     * This enum defines the possible modes used inside the
     * control units for making decissions.
     */
    enum struct ControllerMode {
        RuleBased, ///< Rule-based control with the focus on maximizing PV self-consumption
        OptimizedWithPerfectForecast, ///< Controller using a linear optimization (gurobi) with perfect forecast
    };

    /**
     * This enum defines the possible battery storage (BS) charging and discharging options.
     */
    enum struct ControllerBSGridChargingMode {
        NoGridCharging, ///< Do not allow charging from the grid and discharging into the grid
        OnlyGridCharging, ///< Only allow charging from the grid, but no discharging into the grid
        GridChargingAndDischarging, ///< Allow both, charging from and discharging into the grid
    };

    /**
     * This enum defines the possible optimization targets used inside the
     * control units for making decissions if a optimized controller is used.
     */
    enum struct ControllerOptimizationTarget {
        ElectricityCosts, ///< minimize the costs for electricity consumption (minus revenue for feedin) per control unit
        PeakLoad,  ///< minimize the peak load per control unit
        Emissions, ///< minimize the CO2 emissions caused by grid demand
    };

    /*!
     * The string to delmitit output sections
     */
    const char* const output_section_delimiter = "*********************************************************************************";

}


/*
 * class Global
 *
 * This class contains all global variables that cannot change
 * after they have been set once.
 * 
 * Attention: Not to be confused with namespace global (mind the lower case "g").
 */
class Global {
    public:
        static void InitializeStaticVariables();
        static void DeleteStaticVariables();
        //
        static bool AllVariablesInitialized();
        static void PrintUninitializedVariables(); ///< Prints all variable names to stdout, that are not initialized
        static bool CheckTimeRelatedVariablesInitState(); ///< Returns true if all time-related variables have been initialized. If not, it will return false and prints out all missing variables. It does not include Global::first_timestep and Global::last_timestep (as this checking function is considered to be executed before the calculation of both values). Moreover, it will print out warnings for impossible combinations.
        //
        static void LockAllVariables();   ///< No (set) variable can be overwritten after this call, unset variables can still be set / overwritten once
        static void UnlockAllVariables(); ///< All variables can now be overwritten
        //
        // getter methods
        static unsigned long get_n_timesteps()          { return n_timesteps; }
        static unsigned long get_first_timestep()       { return first_timestep;} //!< Returns the first timestep (starting at 1) in the simulation range
        static unsigned long get_last_timestep()        { return last_timestep; } //!< Returns the last (inclusive) timestep (starting at 1) in the simulation range
        static unsigned long get_n_substations()        { return n_substations; }
        static unsigned long get_n_CUs()                { return n_CUs; }
        static unsigned long get_n_MUs()                { return n_MUs; }
        static unsigned long get_n_pv_profiles()        { return n_pv_ts; }
        static unsigned long get_n_heatpump_profiles()  { return n_hp_ts; }
        static unsigned int  get_seed()                 { return seed; }
        static bool          is_seed_set()              { return seed_set; }
        static bool get_compute_weekly_metrics()        { return compute_weekly_metrics; }
        static bool is_parameter_variation()   { return pvar_selected;  }
        static int  get_parameter_varID()      { return pvar_id;        }
        static bool get_repetitions_selected() { return repetitions_selected; }
        static uint get_n_repetitions()        { return n_repetitions;  }
        static uint get_n_threads()            { return n_threads;      }
        static bool get_work_stealing()        { return work_stealing;  } ///< True, if work stealing is selected between the running worker threads
        static bool get_stop_on_cc_err()       { return stop_on_cc_err; } ///< Stores the value of command line option 'stop_on_cc_err'. True, if the simulation execution should be stopped on if an computation error occurs inside an optimization-based controller in any control unit.
        static struct tm* get_ts_start_tm()    { return ts_start_tm;    }
        static struct tm* get_ts_end_tm()      { return ts_end_tm;    }
        static int get_tsteps_per_hour()       { return tsteps_per_hour;    }
        static unsigned long get_expansion_scenario_id() { return expansion_scenario_id; }
        static float get_time_step_size_in_h() { return time_step_size_in_h; }
        static bool  get_break_sac_loop_if_limit_reached() { return break_sac_loop_if_limit_reached; } //!< Returns true if the SAC loop should be stopped for an individual combination (like PV + HP) if one of the limits is reached (either PV or HP) (even though HP components should still be added)
        static float get_exp_pv_kWp_static()    { return exp_pv_kWp_static;      }
        static float get_exp_pv_kWp_per_m2()    { return exp_pv_kWp_per_m2;      }
        static float get_exp_pv_min_kWp_roof_sec()   { return exp_pv_min_kWp_roof_sec; }  //!< Returns the min. kWp per roof section - smaller sections will be ignored, defaults to 0.0
        static float get_exp_pv_max_kWp_roof_sec()   { return exp_pv_max_kWp_per_sec;  }  //!< Returns the max. kWp per roof section, or -1.0 if it is not set
        static float get_exp_pv_max_kWp_per_unit()   { return exp_pv_max_kWp_per_unit; }  //!< Returns the max. kWp per unit, or -1.0 if it is not set
        static float get_exp_pv_max_kWp_total()      { return exp_pv_max_kWp_total;    }  //!< Returns the max. kWp for all expanded units, or -1.0 if it is not set
        static float get_exp_bess_kW()         { return exp_bess_kW;    }
        static float get_exp_bess_kWh()        { return exp_bess_kWh;    }
        static float get_exp_bess_E_P_ratio()  { return exp_bess_E_P_ratio; }
        static float get_exp_bess_max_capacity(){return exp_bess_max_capacity; } ///< Returns the maxium capacity of a battery storage (important for sizing with other modes than constant). Returns -1 if it is not set.
        static float get_exp_bess_sizingE_boPV(){return exp_bess_sizingE_PV_ratio; } ///!< Returns the factor multiplied with PV nominal power if the battery sizing mode is set to use the PV capacity
        static float get_exp_bess_start_soc()  { return exp_bess_start_soc;    }
        static float get_exp_bess_effi_in()    { return exp_bess_effi_in;    }
        static float get_exp_bess_effi_out()   { return exp_bess_effi_out;   }
        static float get_exp_bess_self_ds_ts() { return exp_bess_self_ds_ts; }
        static float get_exp_bess_P_for_SOC_0(){ return exp_bess_P_for_SOC_0;}
        static float get_exp_bess_P_for_SOC_1(){ return exp_bess_P_for_SOC_1;}
        static float get_exp_bess_max_E_total(){ return exp_bess_max_E_total;}  //!< Returns the upper limit of (resid.) battery capacity for installation, or -1.0 if it is not set
        static float get_exp_bess_max_P_total(){ return exp_bess_max_P_total;}  //!< Returns the upper limit of (resid.) battery power for installation, or -1.0 if it is not set
        static bool          is_exp_hp_max_n_addition_set() { return exp_hp_max_n_addition_set; } //!< Returns true iif a upper limit for the number of added heat pumps is set. The concrete value would be given by Global::get_exp_hp_max_n_addition().
        static unsigned long get_exp_hp_max_n_addition() { return exp_hp_max_n_addition; } //!< Upper limit for the number of added simulated heat pumps
        static bool          is_exp_ev_max_n_addition_set() { return exp_ev_max_n_addition_set; } //!< Returns true iif a upper limit for the number of added EVs is set. The concrete value would be given by Global::get_exp_ev_max_n_addition().
        static unsigned long get_exp_ev_max_n_addition() { return exp_ev_max_n_addition; } //!< Upper limit for the (total) number of added EVs (summed over all control units)
        static unsigned long get_exp_cs_max_ev_per_cs()  { return exp_cs_max_ev_per_cs;  } ///< 0 if not selected. Returns the upper limit for the number of added EVs per control unit / EV charging station. Units with more EVs are not selected.
        static float get_open_space_pv_kWp()   { return open_space_pv_kWp; }
        static float get_wind_kWp()            { return wind_kWp; }
        static float get_feed_in_tariff()      { return feed_in_tariff; }
        static float get_demand_tariff()       { return demand_tariff; }
        static float get_emissions_g_CO2eq_per_kWh() { return emissions_per_kWh; } ///!< Returns the mean emissions in g CO2eq for one kWh of grid demanded energy. If a time series is available, use the latter one!
        static float get_inst_cost_PV_per_kWp(){ return inst_cost_PV_per_kWp; }
        static float get_inst_cost_BS_per_kWh(){ return inst_cost_BS_per_kWh; }
        static float get_inst_cost_HP_per_kW() { return inst_cost_HP_per_kW;  }
        static float get_inst_cost_CS_per_unit(){return inst_cost_CS_per_unit;}
        static float get_npv_discount_rate()   { return npv_discount_rate; }
        static unsigned int get_npv_time_horizon()    { return npv_time_horizon; }
        static double get_npv_factor_if_const()       { return npv_factor_if_const; } ///< Returns the factor for multiplication in the NPV formula for the discounted revenues if the cashflow is constant
        static double get_annuity_factor()            { return 1.0 / npv_factor_if_const; } ///< Returns the annuity factor for economic assesment
        static unsigned int get_hp_flexibility_in_ts()        { return hp_flexibility_in_ts;       } ///< Returns the flexibility of the heat pump profiles in time steps
        static float get_heat_demand_thermalE_to_hpE_conv_f() { return heat_demand_thermalE_to_hpE_conv_f; }
        static float get_heat_cons_bobv_slope()               { return heat_cons_bobv_slope;        } ///< Parameter of linear regression (slope) for estimating annual heat consumption based on the building volume (bobv)
        static float get_heat_cons_bobv_intercept()           { return heat_cons_bobv_intercept;    } ///< Parameter of linear regression (intercept) for estimating annual heat consumption based on the building volume (bobv)
        static float get_ev_plugin_probability()              { return ev_plugin_probability;   }
        static float get_ev_battery_size_kWh()                { return ev_battery_size_kWh;         } ///< The battery capacity in kWh for simulated EVs
        static float get_ev_consumption_kWh_km()              { return ev_consumption_kWh_km;       } ///< The electricity consumption of an EV for driving 1 km
        static float get_ev_max_charging_power_kW()           { return ev_max_charging_power_kW;    } ///< The maximum charging power of a simulated EV
        static float get_cs_max_charging_power_kW()           { return cs_max_charging_power_kW;    } ///< Returns the maximum charging power per simulated ComponentCS or -1, if not set
        static bool  get_use_emission_time_series_ia()        { return use_emission_time_series_ia; } ///< Return weather the emission time series (if it is available in the data) should be used or not
        static bool  get_use_prices_time_series_ia()          { return use_prices_time_series_ia;   } ///< Return weather the electricity prices time series (if it is available in the data) should be used or not
        static bool  get_select_only_residential_buildings()  { return select_only_residential_buildings; } ///< Should only be residential buildings be considered for simulated component addition (SAC planning)?
        static uint  get_control_horizon_in_ts()              { return control_horizon_in_ts;       } ///< The control horizon in time steps (only if controller_mode is set to a value including optimization)
        static uint  get_control_update_freq_in_ts()          { return control_update_freq_in_ts;   } ///< The update frequency of the optimization inside the control units (only if controller_mode is set to a value including optimization)
        static unsigned long get_max_parallel_opti_vars()     { return max_parallel_opti_vars;      } ///< Returns the number of maximum parallel optimization variables, or 0, if not set
        static const std::string& get_input_path()  { return input_path;  }
        static const std::string& get_output_path() { return output_path; }
        static const std::string& get_structure_database_name() { return system_db_name; }
        static const std::string& get_ev_data_path()            { return ev_data_path;   }
        static const std::string& get_cache_dir_path()          { return cache_dir_path_init ? cache_dir_path : input_path ; }///< Returns the path for the cache directory - it this was NOT set, it will return the default value which is the input path
        static global::OutputModePerCU get_output_mode_per_cu() { return output_mode_per_cu; }
        static global::ExpansionProfileAllocationMode get_exp_profile_mode() { return exp_profile_mode; }
        static global::CUSModeFCA get_cu_selection_mode_fca() { return cu_selection_mode_fca; }
        static global::BatteryPowerComputationMode get_battery_power_computation_mode() { return bat_power_comp_mode; }
        static global::BatteryCapacityComputationMode get_battery_capacity_computation_mode() { return bat_capacity_comp_mode; }
        static global::PVSizingMode   get_exp_pv_sizing_mode()   { return exp_pv_sizing_mode; }
        static global::ControllerMode get_controller_mode() { return controller_mode; }
        static global::ControllerBSGridChargingMode get_controller_bs_grid_charging_mode()  { return controller_bs_grid_charging_mode; } ///< Return whether the battery can be charged from the grid / discharged into the grid - only effective if an optimized charging strategy is selected
        static global::ControllerOptimizationTarget get_controller_optimization_target() { return controller_optimization_target; } ///< Returns the selected optimization target (only valid if a controller mode with optimization is selected)
        static const std::set<unsigned long>* get_cu_list_for_sac_planning() { return cu_list_for_sac_planning; } ///< Returns the std::set of control unit IDs (if existing, otherwise NULL) that contains the IDs of the units that should be selected for SAC addition
        static float get_annual_heat_demand_limit_fsac()        { return annual_heat_demand_limit_fsac;  } ///< Returns the upper limit for selection of a control unit for simulative addition based on the annual heat demand in kWh; -1 if this value is not set (thus no limit is given; default)
        static bool get_select_buildings_wg_heatd_only()        { return select_buildings_wg_heatd_only; } ///< True, if only buildings are to be selected for the simulated addition for which an exact specified heat demand is given in the input data - Defaults to false
        static bool get_create_substation_output() { return create_substation_output; } ///< Returns whether a output for the substation time series should be created or not
        static bool get_create_control_cmd_output() { return create_control_cmd_output; } ///< Returns wheather a output for detailed control commands per time step and CU should be created
        static bool get_create_ev_detailed_output() { return create_ev_detailed_output; } ///< Returns wheather a output for detailed EV states per time step and EV should be created
        static const std::string& get_exp_pv_static_profile_orientation() { return exp_pv_static_profile_orientation; }
        static int                get_exp_pv_static_profile_idx()         { return exp_pv_static_profile_idx;         }
        // special setter methods
        static void increment_seed(); ///< Thus class method increments the current seed at anytime during the simulation, no lock is applied to this function.
        // setter methods
        static void set_n_timesteps(unsigned long n_timesteps);
        static void set_first_timestep(unsigned long ts);
        static void set_last_timestep(unsigned long ts);
        static void set_n_substations(unsigned long n_substations);
        static void set_n_CUs(unsigned long n_CUs);
        static void set_n_MUs(unsigned long n_MUs);
        static void set_n_pv_profiles(unsigned long n_pv_ts);
        static void set_n_heatpump_profiles(unsigned long n_hp_ts);
        static void set_seed(unsigned int value);
        static void set_compute_weekly_metrics(bool mode);
        static void set_pvar_vals(bool pvar_set, int pvarID);
        static void set_repetitions_selected(bool value);
        static void set_n_repetitions(unsigned int value);
        static void set_n_threads(unsigned int value);
        static void set_work_stealing(bool value);
        static void set_stop_on_cc_err(bool value);
        static void set_ts_start_tm(struct tm* ts_start_tm);
        static void set_ts_end_tm(struct tm* ts_end_tm);
        static void set_tsteps_per_hour(int tsteps_per_hour);
        static void set_expansion_scenario_id(unsigned long expansion_scenario_id);
        static void set_break_sac_loop_if_limit_reached(bool value);
        static void set_exp_pv_sizing_mode(global::PVSizingMode mode);
        static void set_exp_pv_kWp_static(float value);
        static void set_exp_pv_kWp_per_m2(float value);
        static void set_exp_pv_min_kWp_roof_sec(float value);
        static void set_exp_pv_max_kWp_roof_sec(float value);
        static void set_exp_pv_max_kWp_per_unit(float value);
        static void set_exp_pv_max_kWp_total(float value);
        static void set_exp_bess_kW(float exp_bess_kW);
        static void set_exp_bess_kWh(float exp_bess_kWh);
        static void set_exp_bess_E_P_ratio(float value);
        static void set_exp_bess_max_capacity(float value);
        static void set_exp_bess_sizingE_boPV(float value);
        static void set_exp_bess_start_soc(float exp_bess_start_soc);
        static void set_exp_bess_effi_in(float value);
        static void set_exp_bess_effi_out(float value);
        static void set_exp_bess_self_ds_ts(float value);
        static void set_exp_bess_P_for_SOC_0(float value);
        static void set_exp_bess_P_for_SOC_1(float value);
        static void set_exp_bess_max_E_total(float value);
        static void set_exp_bess_max_P_total(float value);
        static void set_exp_hp_max_n_addition(unsigned long value);
        static void set_exp_ev_max_n_addition(unsigned long value);
        static void set_exp_cs_max_ev_per_cs(unsigned long value);
        static void set_open_space_pv_kWp(float open_space_kWp);
        static void set_wind_kWp(float wind_kWp);
        static void set_feed_in_tariff(float value);
        static void set_demand_tariff(float value);
        static void set_emissions_g_CO2eq_per_kWh(float value);
        static void set_inst_cost_PV_per_kWp(float value);
        static void set_inst_cost_BS_per_kWh(float value);
        static void set_inst_cost_HP_per_kW(float value);
        static void set_inst_cost_CS_per_unit(float value);
        static void set_npv_discount_rate(float value);
        static void set_npv_time_horizon(unsigned int value);
        static void set_hp_flexibility_in_ts(unsigned int value);
        static void set_heat_demand_thermalE_to_hpE_conv_f(float value);
        static void set_heat_cons_bobv_slope(float value);
        static void set_heat_cons_bobv_intercept(float value);
        static void set_ev_plugin_probability(float value);
        static void set_ev_battery_size_kWh(float value);
        static void set_ev_consumption_kWh_km(float value);
        static void set_ev_max_charging_power_kW(float value);
        static void set_cs_max_charging_power_kW(float value);
        static void set_use_emission_time_series_ia(bool use);
        static void set_use_prices_time_series_ia(bool use);
        static void set_select_only_residential_buildings(bool value);
        static void set_control_horizon_in_ts(unsigned int value);
        static void set_control_update_freq_in_ts(unsigned int value);
        static void set_max_parallel_opti_vars(unsigned long value);
        static void set_input_path(std::string* path);
        static void set_output_path(std::string* path);
        static void set_structure_database_name(std::string* fname);
        static void set_ev_data_path(std::string* path);
        static void set_cache_dir_path(const std::string* path);
        static void set_output_mode_per_cu(global::OutputModePerCU mode);
        static void set_exp_profile_mode(global::ExpansionProfileAllocationMode mode);
        static void set_cu_selection_mode_fca(global::CUSModeFCA mode);
        static void set_battery_power_computation_mode(global::BatteryPowerComputationMode mode);
        static void set_battery_capacity_computation_mode(global::BatteryCapacityComputationMode mode);
        static void set_controller_mode(global::ControllerMode mode);
        static void set_controller_bs_grid_charging_mode(global::ControllerBSGridChargingMode mode);
        static void set_controller_optimization_target(global::ControllerOptimizationTarget mode);
        static void set_cu_list_for_sac_planning(const std::set<unsigned long>* selected_cuIDs);
        static void set_annual_heat_demand_limit_fsac(float value);
        static void set_select_buildings_wg_heatd_only(bool value);
        static void set_create_substation_output(bool value);
        static void set_create_control_cmd_output(bool value);
        static void set_create_ev_detailed_output(bool value);
        static void set_exp_pv_static_profile_orientation(std::string* value);
        static void set_exp_pv_static_profile_idx(int value);
    private:
        Global(); ///< Global cannot be initialized, it is a static only class
        static bool is_locked;             ///< if set to true, values cannot be changed anymore
        // variables
        static unsigned long n_timesteps;  ///< Total number of timesteps for which data is available
        static unsigned long first_timestep;
        static unsigned long last_timestep;
        static unsigned long n_substations;///< Total number of substations for which data is available
        static unsigned long n_CUs;        ///< Total number of control units for which data is available
        static unsigned long n_MUs;        ///< Total number of meausrement units for which data is available
        static unsigned long n_pv_ts;      ///< Total number of available normalized pv feedin time series that can be used for simulating new pv installations
        static unsigned long n_hp_ts;      ///< Total number of available normalized heat pump time series that can be used for simulating new heat pumps
        static unsigned int  seed;         ///< The seed for all random number generators
        static bool compute_weekly_metrics;///< Should weekly metrics be computed?
        static bool pvar_selected;         ///< True, if a parameter variation is selected
        static int  pvar_id;               ///< ID of the parameter variation
        static bool repetitions_selected;  ///< True, if complete twin (sac planning + param. vari. + simulation) should run more than once
        static unsigned int n_repetitions; ///< Number of repetitions if repetitions_selected == True
        static unsigned int n_threads;     ///< Number of parallel threads
        static bool work_stealing;
        static bool stop_on_cc_err;
        static struct tm* ts_start_tm;     ///< struct tm of the start date
        static struct tm* ts_end_tm;       ///< struct tm of the end date
        static int tsteps_per_hour;        ///< Time steps per hour in the simulation (and the data!)
        static unsigned long expansion_scenario_id;  ///< ID of the expansion scenario
        static float time_step_size_in_h;  ///< time step size in hours, defines how long a simulation time step is in reality - attention, this global variable is set automatically by set_tsteps_per_hour, it has no own setter
        static bool  break_sac_loop_if_limit_reached; ///< true, if the SAC loop should be stopped for an individual combination (like PV+HP) if one of the limits is reached (either PV or HP)
        static global::PVSizingMode exp_pv_sizing_mode;
        static float exp_pv_kWp_static;    ///< kWp of in the simulation added PV installations, if all units should get the same kWp in the end
        static float exp_pv_kWp_per_m2;    ///< kWp per m2 roof area, if dynamic kWp calculation is selected
        static float exp_pv_min_kWp_roof_sec; ///< minimal size in kWp that a PV section on a given roof section must have so that this section is used for expansion, only applicable with dynamic kWp calculation
        static float exp_pv_max_kWp_per_sec;  ///< maximal size in kWp of a PV section on a given roof section, bigger roof sections will be cut to this value; only applicable when dynamic kWp calculation is selected
        static float exp_pv_max_kWp_per_unit; ///< maximal size of a PV unit (sum over all sections); only applicable when dynamic kWp calculation is selected
        static float exp_pv_max_kWp_total;    ///< maximal size of all added PV units; only applicable when dynamic kWp calculation is selected
        static float exp_bess_kW;          ///< P [kW] of in the simulation added BESS installations
        static float exp_bess_kWh;         ///< E [kWh] of in the simulation added BESS installations
        static float exp_bess_E_P_ratio;   ///< E:P-ratio for new battery storages, this or exp_bess_kW has to be defined!
        static float exp_bess_max_capacity;///< The maxium capacity of a battery storage (important for sizing with other modes than constant)
        static float exp_bess_sizingE_PV_ratio; ///< If the battery capacity computation mode is set to use the nominal PV power, this variable defines the ratio
        static float exp_bess_start_soc;   ///< SOC at the beginning of the simulation for newly added BESS installations
        static float exp_bess_effi_in;     ///< efficiency for charging
        static float exp_bess_effi_out;    ///< efficiency for discharging
        static float exp_bess_self_ds_ts;  ///< self-discharge per time step of the battery
        static float exp_bess_P_for_SOC_0; ///< power consumption of battery if SOC is 0 in kW
        static float exp_bess_P_for_SOC_1; ///< power consumption of battery if SOC is 1 in kW
        static float exp_bess_max_E_total; ///< Upper limit of installed battery storage capacity (only residential)
        static float exp_bess_max_P_total; ///< Upper limit of installed battery storage power (only residential)
        static unsigned long exp_hp_max_n_addition;///< Upper limit for the number of added simulated heat pumps
        static unsigned long exp_ev_max_n_addition;///< Upper limit for the number of added EVs
        static unsigned long exp_cs_max_ev_per_cs;
        static float open_space_pv_kWp;    ///< kWp of the open space PV installations (complete)
        static float wind_kWp;             ///< kWp of the wind turbines
        static float feed_in_tariff;       ///< Tariff for feed in of energy into the grid
        static float demand_tariff;        ///< Tariff for demand from the grid
        static float emissions_per_kWh;    ///< The emissions in CO2eq for one kWh of energy taken from the grid
        static float inst_cost_PV_per_kWp; ///< Installation cost of a PV installation per kWp installed
        static float inst_cost_BS_per_kWh; ///< Installation cost of a battery storage per kWh capacity
        static float inst_cost_HP_per_kW;
        static float inst_cost_CS_per_unit;
        static float npv_discount_rate;    ///< Discount rate for the net present value computation
        static unsigned int npv_time_horizon; ///< Time horizont for the net present value computation
        static double npv_factor_if_const; ///< The NPV constant multiplication factor - no setter available
        static unsigned int hp_flexibility_in_ts;
        static float heat_demand_thermalE_to_hpE_conv_f; ///< Factor for converting thermal energy to heat pump el. energy
        static float heat_cons_bobv_slope;    ///< Parameter of linear regression (coefficient) for estimating annual heat consumption of a building based on its volume
        static float heat_cons_bobv_intercept;///< Parameter of linear regression (intercept) for estimating annual heat consumption of a building based on its volume
        static float ev_plugin_probability;   ///< The probability of plugin in an EV
        static float ev_battery_size_kWh;     ///< The battery capacity in kWh for simulated EVs
        static float ev_consumption_kWh_km;   ///< The electricity consumption of an EV for driving 1 km
        static float ev_max_charging_power_kW;
        static float cs_max_charging_power_kW;
        static bool  use_emission_time_series_ia; ///< Yes, if the emission time series (if it is available in the data) should be used or not
        static bool  use_prices_time_series_ia;   ///< Yes, if the electricity prices time series (if it is available in the data) should be used or not
        static bool  select_only_residential_buildings;
        static uint  control_horizon_in_ts;      ///< The control horizont in time steps (only if controller_mode is set to a value including optimization)
        static uint  control_update_freq_in_ts;   ///< The update frequency of the optimization inside the control units (only if controller_mode is set to a value including optimization)
        static unsigned long max_parallel_opti_vars;
        static std::string input_path;     ///< reference to the string holding the input path of the data
        static std::string output_path;    ///< reference to the string holding the output path of the data
        static std::string system_db_name; ///< String holding the name of the database that contains the system structure
        static std::string ev_data_path;   ///< reference to the string holding the path with EV related data (not required)
        static std::string cache_dir_path;
        static global::OutputModePerCU output_mode_per_cu; ///< Variable storing the selected output mode per CU
        static global::ExpansionProfileAllocationMode exp_profile_mode; ///< Variable storing the selected mode for assigning profiles to PV sections or heat pumps
        static global::CUSModeFCA cu_selection_mode_fca; ///< The selected mode for selecting control units that get sim. added components
        static global::BatteryPowerComputationMode bat_power_comp_mode; ///< The selected mode for computing the (maximal) power of all batteries
        static global::BatteryCapacityComputationMode bat_capacity_comp_mode; ///< The selected mode for computing the battery capacity
        static global::ControllerMode controller_mode; ///< The mode how the control units make decisions
        static global::ControllerBSGridChargingMode controller_bs_grid_charging_mode;
        static global::ControllerOptimizationTarget controller_optimization_target;
        static const std::set<unsigned long>* cu_list_for_sac_planning;
        static float annual_heat_demand_limit_fsac; ///< Select only buildings where the heat demand is lower or equal than the given limit; set to -1 (default) if no limit should be choosen
        static bool select_buildings_wg_heatd_only; ///< Only select buildings with heat demand given in the input data
        static bool create_substation_output; ///< Should an output be created for outputting the substation time series?
        static bool create_control_cmd_output;
        static bool create_ev_detailed_output;
        static std::string exp_pv_static_profile_orientation; ///< fixed orientation for PV static selection mode (ignoring roof data)
        static int exp_pv_static_profile_idx;                 ///< fixed profile ID for PV static selection mode (-1 if not defined)
        // boolean values holding information if the correspoding 
        // variable has been set or not
        static bool n_timesteps_init;
        static bool first_timestep_init;
        static bool last_timestep_init;
        static bool n_substations_init;
        static bool n_CUs_init;
        static bool n_MUs_init;
        static bool n_pv_ts_init;
        static bool n_hp_ts_init;
        static bool seed_set;
        static bool pvar_set;
        static bool repetitions_selected_set;
        static bool n_repetitions_set;
        static bool n_threads_set;
        static bool ts_start_str_init;
        static bool ts_end_str_init;
        static bool tsteps_per_hour_init;
        static bool expansion_scenario_id_init;
        static bool exp_pv_kWp_static_init;
        static bool exp_pv_kWp_per_m2_init;
        static bool exp_pv_min_kWp_roof_sec_init;
        static bool exp_pv_max_kWp_per_sec_init;
        static bool exp_pv_max_kWp_per_unit_init;
        static bool exp_pv_max_kWp_total_init;
        static bool exp_bess_kW_init;
        static bool exp_bess_kWh_init;
        static bool exp_bess_E_P_ratio_init;
        static bool exp_bess_start_soc_init;
        static bool exp_hp_max_n_addition_set;
        static bool exp_ev_max_n_addition_set;
        static bool open_space_pv_kWp_init;
        static bool wind_kWp_init;
        static bool feed_in_tariff_set;
        static bool demand_tariff_set;
        static bool inst_cost_PV_per_kWp_set;
        static bool inst_cost_BS_per_kWh_set;
        static bool npv_discount_rate_set;
        static bool npv_time_horizon_set;
        static bool heat_demand_thermalE_to_hpE_conv_f_set;
        static bool input_path_init;
        static bool output_path_init;
        static bool cache_dir_path_init;
        static bool output_mode_per_cu_init;
        static bool exp_profile_mode_init;
        static bool cu_selection_mode_fca_init;
        static bool bat_power_comp_mode_init;
};

#endif

