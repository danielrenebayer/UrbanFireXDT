#include "global.h"

using namespace global;


#include <cmath>
#include <ctime>
#include <iostream>

using namespace std;



bool global::all_variables_initialized() {
    if (time_info_init &&
        pv_profiles_data != NULL &&
        hp_profiles != NULL &&
        hp_profiles_cumsum != NULL &&
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
    for (struct tm* t : *time_localtime_r) {
        delete t;
    }
    for (struct tm* t : *time_localtime_l) {
        delete t;
    }

    // delete heatpump profiles
    for (size_t hpIdx = 0; hpIdx < Global::get_n_heatpump_profiles(); hpIdx++) {
        delete[] hp_profiles[hpIdx];
        delete[] hp_profiles_cumsum[hpIdx];
    }
    delete[] hp_profiles; hp_profiles = NULL;
    delete[] hp_profiles_cumsum; hp_profiles_cumsum = NULL;
    // delete pv profiles
    for (size_t pvIdx = 0; pvIdx < Global::get_n_pv_profiles(); pvIdx++)
        delete[] pv_profiles_data[pvIdx];
    delete[] pv_profiles_data; pv_profiles_data = NULL;

    delete[] time_timestep_id;       time_timestep_id       = NULL;
    delete   time_localtime_r;       time_localtime_r       = NULL;
    delete   time_localtime_l;       time_localtime_l       = NULL;
    delete   time_localtimezone_str; time_localtimezone_str = NULL;
    delete   unit_open_space_pv;     unit_open_space_pv   = NULL;
    delete   unit_open_space_wind;   unit_open_space_wind = NULL;
    delete   parameter_var_list;     parameter_var_list   = NULL;
    delete   current_output_dir;     current_output_dir   = NULL;
    delete   current_output_dir_prefix;current_output_dir_prefix   = NULL;
    delete   current_global_output_dir;current_global_output_dir   = NULL;
    delete[] wind_profile; wind_profile = NULL;
    delete[] residual_gridload_kW;   residual_gridload_kW = NULL;

    if (emission_ts      != NULL) { delete[] emission_ts;      emission_ts      = NULL; }
    if (eprices_local_ts != NULL) { delete[] eprices_local_ts; eprices_local_ts = NULL; }
    if (eprices_spotm_ts != NULL) { delete[] eprices_spotm_ts; eprices_spotm_ts = NULL; }
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
    if (hp_profiles_cumsum == NULL) {
        cout << "HP profiles (variable hp_profiles_cumsum) undefined." << endl;
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

bool Global::is_locked            = false;
unsigned long Global::n_timesteps = 0;
unsigned long Global::first_timestep= 0;
unsigned long Global::last_timestep = 0;
unsigned long Global::n_substations = 0;
unsigned long Global::n_CUs       = 0;
unsigned long Global::n_MUs       = 0;
unsigned long Global::n_pv_ts     = 0;
unsigned long Global::n_hp_ts     = 0;
unsigned int  Global::seed        = 0;
bool Global::compute_weekly_metrics = false;
bool Global::pvar_selected        = false;
int  Global::pvar_id              = 0;
bool Global::repetitions_selected = false;
uint Global::n_repetitions        = 0;
uint Global::n_threads            = 3;
bool Global::work_stealing        = false;
bool Global::stop_on_cc_err       = false;
struct tm* Global::ts_start_tm    = NULL;
struct tm* Global::ts_end_tm      = NULL;
int Global::tsteps_per_hour       = 1;
unsigned long Global::expansion_scenario_id = 0;
float Global::time_step_size_in_h   = 0.0;
bool  Global::break_sac_loop_if_limit_reached = true;
float Global::exp_pv_kWp_static       = 0.0;
float Global::exp_pv_kWp_per_m2       = 0.0;
float Global::exp_pv_min_kWp_roof_sec = 0.0;
float Global::exp_pv_max_kWp_per_sec  =-1.0;
float Global::exp_pv_max_kWp_per_unit =-1.0;
float Global::exp_pv_max_kWp_total    =-1.0;
float Global::exp_bess_kW           = 0.0;
float Global::exp_bess_kWh          = 0.0;
float Global::exp_bess_E_P_ratio    = 0.0;
float Global::exp_bess_max_capacity =-1.0;
float Global::exp_bess_sizingE_PV_ratio = 1.0;
float Global::exp_bess_start_soc    = 0.0;
float Global::exp_bess_effi_in      = 1.0;
float Global::exp_bess_effi_out     = 1.0;
float Global::exp_bess_self_ds_ts   = 0.0;
float Global::exp_bess_P_for_SOC_0  = 0.0;
float Global::exp_bess_P_for_SOC_1  = 0.0;
float Global::exp_bess_max_E_total  =-1.0;
float Global::exp_bess_max_P_total  =-1.0;
ulong Global::exp_hp_max_n_addition = 0;
ulong Global::exp_ev_max_n_addition = 0;
ulong Global::exp_cs_max_ev_per_cs  = 0;
float Global::open_space_pv_kWp     = 0.0;
float Global::wind_kWp              = 0.0;
float Global::feed_in_tariff        = 0.0;
float Global::demand_tariff         = 0.0;
float Global::emissions_per_kWh     = 100.0;
float Global::inst_cost_PV_per_kWp  = 0.0;
float Global::inst_cost_BS_per_kWh  = 0.0;
float Global::npv_discount_rate     = 0.0;
unsigned int Global::npv_time_horizon  = 0;
double Global::npv_factor_if_const  = 0.0;
unsigned int Global::hp_flexibility_in_ts = 1;
float Global::heat_demand_thermalE_to_hpE_conv_f = 0.0;
float Global::heat_cons_bobv_slope     = 0.0;
float Global::heat_cons_bobv_intercept = 0.0;
float Global::ev_plugin_probability   = 0.25;
float Global::ev_battery_size_kWh    = 30.0;
float Global::ev_consumption_kWh_km   = 0.2f;
float Global::ev_max_charging_power_kW=11.0f;
float Global::ev_charging_effi        = 1.0f;
float Global::cs_max_charging_power_kW=-1.0f;
bool  Global::use_emission_time_series_ia = true;
bool  Global::use_prices_time_series_ia   = true;
bool  Global::select_only_residential_buildings = false;
uint  Global::control_horizon_in_ts       =   24;
uint  Global::control_update_freq_in_ts   =    1;
ulong Global::max_parallel_opti_vars      =    0;
string Global::input_path         = "";
string Global::output_path        = "";
string Global::system_db_name     = "SystemStructure.db";
string Global::ev_data_path       = "";
OutputModePerCU Global::output_mode_per_cu = OutputModePerCU::IndividualFile;
ExpansionProfileAllocationMode Global::exp_profile_mode = ExpansionProfileAllocationMode::AsInData;
global::CUSModeFCA Global::cu_selection_mode_fca        = global::CUSModeFCA::OrderAsInData;
global::BatteryPowerComputationMode Global::bat_power_comp_mode = global::BatteryPowerComputationMode::AsDefinedByConfigVar;
global::BatteryCapacityComputationMode Global::bat_capacity_comp_mode = global::BatteryCapacityComputationMode::Constant;
global::PVSizingMode Global::exp_pv_sizing_mode = global::PVSizingMode::MaxAvailableRoofArea;
global::ControllerMode Global::controller_mode = global::ControllerMode::RuleBased;
global::ControllerBSGridChargingMode Global::controller_bs_grid_charging_mode = global::ControllerBSGridChargingMode::NoGridCharging;
global::ControllerOptimizationTarget Global::controller_optimization_target = global::ControllerOptimizationTarget::ElectricityCosts;
const std::set<unsigned long>* Global::cu_list_for_sac_planning = NULL;
float Global::annual_heat_demand_limit_fsac  = -1;
bool Global::select_buildings_wg_heatd_only  = false;
bool Global::create_substation_output = true;
bool Global::create_control_cmd_output       = false;
bool Global::create_ev_detailed_output       = false;
string Global::exp_pv_static_profile_orientation = "";
int Global::exp_pv_static_profile_idx            = -1;
//
bool Global::n_timesteps_init      = false;
bool Global::first_timestep_init   = false;
bool Global::last_timestep_init    = false;
bool Global::n_substations_init    = false;
bool Global::n_CUs_init            = false;
bool Global::n_MUs_init            = false;
bool Global::n_pv_ts_init          = false;
bool Global::n_hp_ts_init          = false;
bool Global::seed_set              = false;
bool Global::pvar_set              = false;
bool Global::repetitions_selected_set     = false;
bool Global::n_repetitions_set     = false;
bool Global::n_threads_set         = false;
bool Global::ts_start_str_init     = false;
bool Global::ts_end_str_init       = false;
bool Global::tsteps_per_hour_init  = false;
bool Global::expansion_scenario_id_init = false;
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
bool Global::exp_hp_max_n_addition_set  = false;
bool Global::exp_ev_max_n_addition_set  = false;
bool Global::open_space_pv_kWp_init= false;
bool Global::wind_kWp_init         = false;
bool Global::feed_in_tariff_set    = false;
bool Global::demand_tariff_set     = false;
bool Global::inst_cost_PV_per_kWp_set     = false;
bool Global::inst_cost_BS_per_kWh_set     = false;
bool Global::npv_discount_rate_set = false;
bool Global::npv_time_horizon_set  = false;
bool Global::heat_demand_thermalE_to_hpE_conv_f_set = false;
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
    if (Global::cu_list_for_sac_planning != NULL) {
        delete Global::cu_list_for_sac_planning;
        Global::cu_list_for_sac_planning = NULL;
    }
}

bool Global::AllVariablesInitialized() {
    // Check for impossible combinations
    if (cs_max_charging_power_kW > 0.0 && controller_mode != ControllerMode::RuleBased ) {
        std::cerr << "Warning: Setting a upper limit for max charging power per control unit will be ignored if controller mode is other than rule-based." << std::endl;
    }
    if (control_horizon_in_ts < control_update_freq_in_ts) {
        std::cerr << "Error: Control horizon < control update frequency!" << std::endl;
        return false;
    }
    if ((
            exp_pv_sizing_mode == global::PVSizingMode::Optimized  ||
            bat_capacity_comp_mode == global::BatteryCapacityComputationMode::Optimized
        ) && controller_mode != global::ControllerMode::OptimizedWithPerfectForecast
    ) {
        std::cerr << "Error: PV sizing mode or BS capacity computation mode ist set to be optimized, but the controller mode is not set to use an optimized strategy. This is an impossible combination." << std::endl;
        return false;
    }
    // Check for initialized variables
    if (n_timesteps_init && 
        n_substations_init &&
        n_CUs_init &&
        n_MUs_init &&
        n_pv_ts_init &&
        n_hp_ts_init &&
        pvar_set &&
        repetitions_selected_set &&
        n_threads_set &&
        ts_start_str_init &&
        ts_end_str_init &&
        tsteps_per_hour_init &&
        expansion_scenario_id_init &&
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
        heat_demand_thermalE_to_hpE_conv_f_set &&
        inst_cost_PV_per_kWp_set &&
        inst_cost_BS_per_kWh_set &&
        input_path_init &&
        bat_power_comp_mode_init)
    {
        bool exp_pv_kWp_static_mode = exp_pv_sizing_mode == global::PVSizingMode::StaticPVSize; // only a helper variable in this function
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

                if ( (bat_capacity_comp_mode == global::BatteryCapacityComputationMode::Optimized ||
                      exp_pv_sizing_mode == global::PVSizingMode::Optimized )
                    &&
                     (controller_optimization_target != global::ControllerOptimizationTarget::ElectricityCosts ||
                      controller_mode != global::ControllerMode::OptimizedWithPerfectForecast)
                ) {
                    std::cerr << "Error: When the battery or the PV installation should be sized to optimize the electricity costs, the configuration parameter 'controller optimization target' must be set to 'electricity costs' and 'controller mode' must be set to 'opti with perfect forecast'!" << std::endl;
                    return false;
                }

                if (exp_pv_sizing_mode == global::PVSizingMode::Optimized && exp_pv_max_kWp_total >= 0.0) {
                    std::cerr << "Warning: Parameter 'expansion PV max total kWp addition' has no effect when using 'expansion PV sizing mode' == 'Optimized'!" << std::endl;
                }
                if (bat_capacity_comp_mode == global::BatteryCapacityComputationMode::Optimized && ( exp_bess_max_P_total >= 0.0 || exp_bess_max_E_total >= 0.0 )) {
                    std::cerr << "Warning: Parameter 'expansion PV max total kWp addition' and 'expansion BS max total E addition' have no effect when using 'expansion BS capacity computation mode' == 'Optimized'!" << std::endl;
                }
                if (bat_capacity_comp_mode == global::BatteryCapacityComputationMode::Optimized && Global::get_battery_power_computation_mode() != global::BatteryPowerComputationMode::UseEOverPRatio) {
                    std::cerr << "Error: Configuration parameter 'expansion BS capacity computation mode' == 'Optimized' can only be used in conjunction with 'expansion BS power computation mode' == 'Use E:P-ratio'!" << std::endl;
                    return false;
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
    if (!pvar_set) {
        cout << "Variable pvar not initialized." << endl;
    }
    if (!repetitions_selected_set) {
        cout << "Variable repetitions_selected not initialized." << endl;
    }
    if (repetitions_selected && !n_repetitions_set) {
        cout << "Variable n_repetitions_set not initialized." << endl;
    }
    if (!n_threads_set) {
        cout << "Variable n_threads not initialized." << endl;
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
    if (!exp_pv_kWp_static_init && Global::get_exp_pv_sizing_mode() == global::PVSizingMode::StaticPVSize) {
        cout << "Variable exp_pv_kWp_static not initialized." << endl;
    }
    if (!exp_pv_kWp_per_m2_init && Global::get_exp_pv_sizing_mode() != global::PVSizingMode::StaticPVSize) {
        cout << "Variable exp_pv_kWp_per_m2 not initialized." << endl;
    }
    if (!exp_pv_min_kWp_roof_sec_init && Global::get_exp_pv_sizing_mode() != global::PVSizingMode::StaticPVSize) {
        cout << "Variable exp_pv_min_kWp_roof_sec not initialized." << endl;
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
    if (!heat_demand_thermalE_to_hpE_conv_f_set) {
        cout << "Variable heat_demand_thermalE_to_hpE_conv_f not initialized." << endl;
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

bool Global::CheckTimeRelatedVariablesInitState() {
    if (!n_timesteps_init) {
        cout << "Variable n_timesteps not initialized." << endl;
        return false;
    }
    if (!ts_start_str_init) {
        cout << "Variable ts_start_str not initialized." << endl;
        return false;
    }
    if (!ts_end_str_init) {
        cout << "Variable ts_end_str not initialized." << endl;
        return false;
    }
    return true;
}

void Global::LockAllVariables() {
    Global::is_locked = true;
}

void Global::UnlockAllVariables() {
    Global::is_locked = false;
}

void Global::increment_seed() {
    seed += 1;
}

void Global::set_n_timesteps(unsigned long n_timesteps) {
    if (is_locked && n_timesteps_init) {
        cerr << "Global variable n_timesteps is already initialized!" << endl;
    } else {
        Global::n_timesteps = n_timesteps;
        Global::n_timesteps_init = true;
    }
}
void Global::set_first_timestep(unsigned long ts) {
    if (first_timestep_init) {
        cerr << "Global variable first_timestep is already initialized! It can only be initialized once, regardless of the lock-state." << endl;
    } else {
        Global::first_timestep = ts;
        Global::first_timestep_init = true;
    }
}
void Global::set_last_timestep(unsigned long ts) {
    if (last_timestep_init) {
        cerr << "Global variable last_timestep is already initialized! It can only be initialized once, regardless of the lock-state." << endl;
    } else {
        Global::last_timestep = ts;
        Global::last_timestep_init = true;
    }
}
void Global::set_n_substations(unsigned long n_substations) {
    if (is_locked && n_substations_init) {
        cerr << "Global variable n_substations is already initialized!" << endl;
    } else {
        Global::n_substations = n_substations;
        Global::n_substations_init = true;
    }
}
void Global::set_n_CUs(unsigned long n_CUs) {
    if (is_locked && n_CUs_init) {
        cerr << "Global variable n_CUs is already initialized!" << endl;
    } else {
        Global::n_CUs = n_CUs;
        Global::n_CUs_init = true;
    }
}
void Global::set_n_pv_profiles(unsigned long value) {
    if (is_locked && n_pv_ts_init) {
        cerr << "Global variable n_pv_profiles is already initialized!" << endl;
    } else {
        Global::n_pv_ts = value;
        Global::n_pv_ts_init = true;
    }
}
void Global::set_n_heatpump_profiles(unsigned long n_hp_ts) {
    if (is_locked && n_hp_ts_init) {
        cerr << "Global variable n_heatpump_profiles is already initialized!" << endl;
    } else {
        Global::n_hp_ts = n_hp_ts;
        Global::n_hp_ts_init = true;
    }
}
void Global::set_n_MUs(unsigned long n_MUs) {
    if (is_locked && n_MUs_init) {
        cerr << "Global variable n_MUs is already initialized!" << endl;
    } else {
        Global::n_MUs = n_MUs;
        Global::n_MUs_init = true;
    }
}
void Global::set_seed(unsigned int value) {
    if (is_locked && seed_set) {
        cerr << "Global variable seed is already initialized!" << endl;
    } else {
        Global::seed = value;
        Global::seed_set = true;
    }
}
void Global::set_compute_weekly_metrics(bool mode) {
    if (is_locked) {
        cerr << "Global variable compute_weekly_metrics is already initialized!" << endl;
    } else {
        Global::compute_weekly_metrics = mode;
    }
}
void Global::set_pvar_vals(bool pvar_val, int pvarID) {
    if (is_locked && pvar_set) {
        cerr << "Values for parameter variation are already set!" << endl;
    } else {
        Global::pvar_selected = pvar_val;
        Global::pvar_id  = pvarID;
        Global::pvar_set = true;
    }
}
void Global::set_repetitions_selected(bool value) {
    if (is_locked && repetitions_selected_set){
        cerr << "Global variable repetitions_selected is already initialized!" << endl;
    } else {
        Global::repetitions_selected = value;
        Global::repetitions_selected_set = true;
    }
}
void Global::set_n_repetitions(unsigned int value) {
    if (is_locked && n_repetitions_set){
        cerr << "Global variable n_repetitions is already initialized!" << endl;
    } else {
        Global::n_repetitions = value;
        Global::n_repetitions_set = true;
    }
}
void Global::set_n_threads(unsigned int value) {
    if (is_locked || n_threads_set){
        cerr << "Global variable n_threads is already initialized!" << endl;
    } else {
        Global::n_threads = value;
        Global::n_threads_set = true;
    }
}
void Global::set_work_stealing(bool value) {
    if (is_locked) {
        cerr << "Global variable work_stealing cannot be set at the moment!" << endl;
    } else {
        Global::work_stealing = value;
    }
}
void Global::set_stop_on_cc_err(bool value) {
    if (is_locked){
        cerr << "Global variable stop_on_cc_err cannot be set at the moment!" << endl;
    } else {
        Global::stop_on_cc_err = value;
    }
}
void Global::set_ts_start_tm(struct tm* ts_start_tm) {
    if (is_locked && ts_start_str_init) {
        cerr << "Global variable ts_start_str is already initialized!" << endl;
    } else {
        Global::ts_start_tm = ts_start_tm;
        Global::ts_start_str_init = true;
    }
}
void Global::set_ts_end_tm(struct tm* ts_end_tm) {
    if (is_locked && ts_end_str_init) {
        cerr << "Global variable ts_end_str is already initialized!" << endl;
    } else {
        Global::ts_end_tm = ts_end_tm;
        Global::ts_end_str_init = true;
    }
}
void Global::set_tsteps_per_hour(int tsteps_per_hour) {
    if (is_locked && tsteps_per_hour_init) {
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
void Global::set_expansion_scenario_id(unsigned long expansion_scenario_id) {
    if (is_locked && expansion_scenario_id_init) {
        cerr << "Global variable expansion_scenario_id is already initialized!" << endl;
    } else {
        Global::expansion_scenario_id = expansion_scenario_id;
        Global::expansion_scenario_id_init = true;
    }
}
void Global::set_break_sac_loop_if_limit_reached(bool value) {
    if (is_locked) {
        cerr << "Global variable break_sac_loop_if_limit_reached is already initialized!" << endl;
    } else {
        Global::break_sac_loop_if_limit_reached = value;
    }
}
void Global::set_exp_pv_sizing_mode(global::PVSizingMode mode) {
    if (is_locked) {
        cerr << "Global variable exp_pv_sizing_mode cannot be modified at the moment!" << endl;
    } else {
        Global::exp_pv_sizing_mode = mode;
    }
}
void Global::set_exp_pv_kWp_static(float value) {
    if (is_locked && exp_pv_kWp_static_init) {
        cerr << "Global variable exp_pv_kWp_static is already initialized!" << endl;
    } else {
        Global::exp_pv_kWp_static      = value;
        Global::exp_pv_kWp_static_init = true;
    }
}
void Global::set_exp_pv_kWp_per_m2(float value) {
    if (is_locked && exp_pv_kWp_per_m2_init) {
        cerr << "Global variable exp_pv_kWp_per_m2 is already initialized!" << endl;
    } else {
        Global::exp_pv_kWp_per_m2      = value;
        Global::exp_pv_kWp_per_m2_init = true;
    }
}
void Global::set_exp_pv_min_kWp_roof_sec(float value) {
    if (is_locked && exp_pv_min_kWp_roof_sec_init) {
        cerr << "Global variable exp_pv_min_kWp_roof_sec is already initialized!" << endl;
    } else {
        Global::exp_pv_min_kWp_roof_sec      = value;
        Global::exp_pv_min_kWp_roof_sec_init = true;
    }
}
void Global::set_exp_pv_max_kWp_roof_sec(float value) {
    if (is_locked && exp_pv_max_kWp_per_sec_init) {
        cerr << "Global variable exp_pv_max_kWp_per_sec is already initialized!" << endl;
    } else {
        Global::exp_pv_max_kWp_per_sec      = value;
        Global::exp_pv_max_kWp_per_sec_init = true;
    }
}
void Global::set_exp_pv_max_kWp_per_unit(float value) {
    if (is_locked && Global::exp_pv_max_kWp_per_unit_init) {
        cerr << "Global variable exp_pv_max_kWp_per_unit is already initialized!" << endl;
    } else {
        Global::exp_pv_max_kWp_per_unit      = value;
        Global::exp_pv_max_kWp_per_unit_init = true;
    }
}
void Global::set_exp_pv_max_kWp_total(float value) {
    if (is_locked && Global::exp_pv_max_kWp_total_init) {
        cerr << "Global variable exp_pv_max_kWp_total is already initialized!" << endl;
    } else {
        Global::exp_pv_max_kWp_total      = value;
        Global::exp_pv_max_kWp_total_init = true;
    }
}
void Global::set_exp_bess_kW(float exp_bess_kW) {
    if (is_locked && exp_bess_kW_init) {
        cerr << "Global variable exp_bess_kW is already initialized!" << endl;
    } else {
        Global::exp_bess_kW = exp_bess_kW;
        Global::exp_bess_kW_init = true;
    }
}
void Global::set_exp_bess_kWh(float exp_bess_kWh) {
    if (is_locked && exp_bess_kWh_init) {
        cerr << "Global variable exp_bess_kWh is already initialized!" << endl;
    } else {
        Global::exp_bess_kWh = exp_bess_kWh;
        Global::exp_bess_kWh_init = true;
    }
}
void Global::set_exp_bess_E_P_ratio(float value) {
    if (is_locked && exp_bess_E_P_ratio_init) {
        cerr << "Global variable exp_bess_E_P_ratio is already initialized!" << endl;
    } else {
        Global::exp_bess_E_P_ratio = value;
        Global::exp_bess_E_P_ratio_init = true;
    }
}
void Global::set_exp_bess_max_capacity(float value) {
    if (is_locked) {
        cerr << "Global variable exp_bess_max_capacity cannot be overwritten at the moment!" << endl;
    } else {
        if (value <= 0 && value != -1.0) { //-1.0 -> variable unset
            cerr << "Global variable exp_bess_max_capacity cannot be set to value " << value << " (allowed range: ]0.0,inf[ )" << endl;
            return;
        } else {
            Global::exp_bess_max_capacity = value;
        }
    }
}
void Global::set_exp_bess_sizingE_boPV(float value) {
    if (is_locked && exp_bess_kWh_init) {
        cerr << "Global variable exp_bess_sizingE_boPV is already initialized!" << endl;
    } else {
        if (value <= 0) {
            cerr << "Global variable exp_bess_sizingE_boPV cannot be set to value " << value << " (allowed range: ]0.0,inf[ )" << endl;
            return;
        } else {
            Global::exp_bess_sizingE_PV_ratio = value;
        }
    }
}
void Global::set_exp_bess_start_soc(float exp_bess_start_soc) {
    if (is_locked && exp_bess_start_soc_init) {
        cerr << "Global variable exp_bess_start_soc is already initialized!" << endl;
    } else {
        if (exp_bess_start_soc > 1.0 || exp_bess_start_soc < 0.0) {
            cerr << "Global variable exp_bess_start_soc cannot be set to value " << exp_bess_start_soc << " (allowed range: [0.0,1.0] )" << endl;
            return;
        }
        Global::exp_bess_start_soc = exp_bess_start_soc;
        Global::exp_bess_start_soc_init = true;
    }
}
void Global::set_exp_bess_effi_in(float value) {
    if (is_locked) {
        cerr << "Global variable exp_bess_effi_in cannot be overwritten at the moment!" << endl;
    } else {
        if (value > 1.0 || value <= 0.0) {
            cerr << "Global variable exp_bess_effi_in cannot be set to value " << value << " (allowed range: ]0.0,1.0] )" << endl;
            return;
        }
        Global::exp_bess_effi_in = value;
    }
}
void Global::set_exp_bess_effi_out(float value) {
    if (is_locked) {
        cerr << "Global variable exp_bess_effi_in cannot be overwritten at the moment!" << endl;
    } else {
        if (value > 1.0 || value <= 0.0) {
            cerr << "Global variable exp_bess_effi_out cannot be set to value " << value << " (allowed range: ]0.0,1.0] )" << endl;
            return;
        }
        Global::exp_bess_effi_out = value;
    }
}
void Global::set_exp_bess_self_ds_ts(float value) {
    if (is_locked) {
        cerr << "Global variable exp_bess_effi_in cannot be overwritten at the moment!" << endl;
    } else {
        if (value > 1.0 || value < 0.0) {
            cerr << "Global variable exp_bess_self_ds_ts cannot be set to value " << value << " (allowed range: [0.0,1.0] )" << endl;
            return;
        }
        Global::exp_bess_self_ds_ts = value;
    }
}
void Global::set_exp_bess_P_for_SOC_0(float value) {
    if (is_locked) {
        cerr << "Global variable exp_bess_P_for_SOC_0 cannot be overwritten at the moment!" << endl;
    } else {
        if (value < 0.0) {
            cerr << "Global variable exp_bess_P_for_SOC_0 cannot be set to value " << value << " (allowed range: [0.0,inf) )" << endl;
            return;
        }
        Global::exp_bess_P_for_SOC_0 = value;
    }
}
void Global::set_exp_bess_P_for_SOC_1(float value) {
    if (is_locked) {
        cerr << "Global variable exp_bess_P_for_SOC_1 cannot be overwritten at the moment!" << endl;
    } else {
        if (value < 0.0) {
            cerr << "Global variable exp_bess_P_for_SOC_1 cannot be set to value " << value << " (allowed range: [0.0,inf) )" << endl;
            return;
        }
        Global::exp_bess_P_for_SOC_1 = value;
    }
}
void Global::set_exp_bess_max_E_total(float value) {
    if (is_locked) {
        cerr << "Global variable exp_bess_max_E_total cannot be overwritten at the moment!" << endl;
    } else {
        if (value < 0.0) {
            cerr << "Global variable exp_bess_max_E_total cannot be set to value " << value << " (allowed range: [0.0,inf) )" << endl;
            return;
        }
        Global::exp_bess_max_E_total = value;
    }
}
void Global::set_exp_bess_max_P_total(float value) {
    if (is_locked) {
        cerr << "Global variable exp_bess_max_P_total cannot be overwritten at the moment!" << endl;
    } else {
        if (value < 0.0) {
            cerr << "Global variable exp_bess_max_P_total cannot be set to value " << value << " (allowed range: [0.0,inf) )" << endl;
            return;
        }
        Global::exp_bess_max_P_total = value;
    }
}
void Global::set_exp_hp_max_n_addition(unsigned long value) {
    if (is_locked) {
        cerr << "Global variable exp_hp_max_n_addition cannot be overwritten at the moment!" << endl;
    } else {
        Global::exp_hp_max_n_addition = value;
        Global::exp_hp_max_n_addition_set = true;
    }
}
void Global::set_exp_ev_max_n_addition(unsigned long value) {
    if (is_locked) {
        cerr << "Global variable exp_ev_max_n_addition cannot be overwritten at the moment!" << endl;
    } else {
        Global::exp_ev_max_n_addition = value;
        Global::exp_ev_max_n_addition_set = true;
    }
}
void Global::set_exp_cs_max_ev_per_cs(unsigned long value) {
    if (is_locked) {
        cerr << "Global variable exp_cs_max_ev_per_cs cannot be overwritten at the moment!" << endl;
    } else {
        Global::exp_cs_max_ev_per_cs = value;
    }
}
void Global::set_open_space_pv_kWp(float open_space_kWp) {
    if (is_locked && open_space_pv_kWp_init) {
        cerr << "Global variable open_space_pv_kWp is already initialized!" << endl;
    } else {
        Global::open_space_pv_kWp = open_space_kWp;
        open_space_pv_kWp_init = true;
    }
}
void Global::set_wind_kWp(float wind_kWp) {
    if (is_locked && wind_kWp_init) {
        cerr << "Global variable wind_kWp is already initialized!" << endl;
    } else {
        Global::wind_kWp = wind_kWp;
        wind_kWp_init = true;
    }
}
void Global::set_feed_in_tariff(float value) {
    if (is_locked && feed_in_tariff_set) {
        cerr << "Global variable feed_in_tariff is already initialized!" << endl;
    } else {
        Global::feed_in_tariff = value;
        feed_in_tariff_set = true;
    }
}
void Global::set_demand_tariff(float value) {
    if (is_locked && demand_tariff_set) {
        cerr << "Global variable demand_tariff is already initialized!" << endl;
    } else {
        Global::demand_tariff = value;
        demand_tariff_set = true;
    }
}
void Global::set_emissions_g_CO2eq_per_kWh(float value) {
    if (is_locked) {
        cerr << "Global variable emissions_per_kWh cannot be set at the moment!" << endl;
    } else {
        Global::emissions_per_kWh = value;
    }
}
void Global::set_npv_discount_rate(float value) {
    if (is_locked && npv_discount_rate_set) {
        cerr << "Global variable npv_discount_rate is already initialized!" << endl;
    } else {
        Global::npv_discount_rate = value;
        npv_discount_rate_set = true;
    }
    // Compute npv_factor_if_const if all required values are defined
    if (npv_time_horizon_set) {
        npv_factor_if_const =
          ( pow( (double) (1 + npv_discount_rate), (double) (npv_time_horizon)) - 1 ) /
          ( pow( (double) (1 + npv_discount_rate), (double) (npv_time_horizon)) * (double) (npv_discount_rate) );
    }
}
void Global::set_npv_time_horizon(unsigned int value) {
    if (is_locked && npv_time_horizon_set) {
        cerr << "Global variable npv_time_horizon is already initialized!" << endl;
    } else {
        Global::npv_time_horizon = value;
        npv_time_horizon_set = true;
    }
    // Compute npv_factor_if_const if all required values are defined
    if (npv_discount_rate_set) {
        npv_factor_if_const =
          ( pow( (double) (1 + npv_discount_rate), (double) (npv_time_horizon)) - 1 ) /
          ( pow( (double) (1 + npv_discount_rate), (double) (npv_time_horizon)) * (double) (npv_discount_rate) );
    }
}
void Global::set_inst_cost_PV_per_kWp(float value) {
    if (is_locked && inst_cost_PV_per_kWp_set) {
        cerr << "Global variable inst_cost_PV_per_kWp is already initialized!" << endl;
    } else {
        Global::inst_cost_PV_per_kWp = value;
        inst_cost_PV_per_kWp_set = true;
    }
}
void Global::set_inst_cost_BS_per_kWh(float value) {
    if (is_locked && inst_cost_BS_per_kWh_set) {
        cerr << "Global variable inst_cost_BS_per_kWh is already initialized!" << endl;
    } else {
        Global::inst_cost_BS_per_kWh = value;
        inst_cost_BS_per_kWh_set = true;
    }
}
void Global::set_hp_flexibility_in_ts(unsigned int value) {
    if (is_locked) {
        cerr << "Global variable hp_flexibility_in_ts cannot be set at the moment!" << endl;
    } else {
        Global::hp_flexibility_in_ts = value;
    }
}
void Global::set_heat_demand_thermalE_to_hpE_conv_f(float value) {
    if (is_locked && heat_demand_thermalE_to_hpE_conv_f_set) {
        cerr << "Global variable heat_demand_thermalE_to_hpE_conv_f is already initialized!" << endl;
    } else {
        Global::heat_demand_thermalE_to_hpE_conv_f = value;
        heat_demand_thermalE_to_hpE_conv_f_set = true;
    }
}
void Global::set_heat_cons_bobv_slope(float value) {
    if (is_locked) {
        cerr << "Global variable heat_cons_bobv_slope cannot be set at the moment!" << endl;
    } else {
        Global::heat_cons_bobv_slope = value;
    }
}
void Global::set_heat_cons_bobv_intercept(float value) {
    if (is_locked) {
        cerr << "Global variable heat_cons_bobv_intercept cannot be set at the moment!" << endl;
    } else {
        Global::heat_cons_bobv_intercept = value;
    }
}
void Global::set_ev_plugin_probability(float value) {
    if (is_locked) {
        cerr << "Global variable ev_plugin_probability cannot be set!" << endl;
    } else {
        if (value > 1.0 || value < 0.0) {
            cerr << "Global variable ev_plugin_probability cannot be set to value " << value << " (allowed range: [0.0,1.0] )" << endl;
            return;
        }
        Global::ev_plugin_probability = value;
    }
}
void Global::set_ev_battery_size_kWh(float value) {
    if (is_locked) {
        cerr << "Global variable ev_battery_size_kWh cannot be set!" << endl;
    } else {
        if (value < 0.0) {
            cerr << "Global variable ev_battery_size_kWh cannot be set to value " << value << " ( allowed range: [0.0,inf[ )" << endl;
            return;
        }
        Global::ev_battery_size_kWh = value;
    }
}
void Global::set_ev_consumption_kWh_km(float value) {
    if (is_locked) {
        cerr << "Global variable ev_consumption_kWh_km cannot be set!" << endl;
    } else {
        if (value < 0.0) {
            cerr << "Global variable ev_consumption_kWh_km cannot be set to value " << value << " ( allowed range: [0.0,inf[ )" << endl;
            return;
        }
        Global::ev_consumption_kWh_km = value;
    }
}
void Global::set_ev_max_charging_power_kW(float value) {
    if (is_locked) {
        cerr << "Global variable ev_max_charging_power_kW cannot be set!" << endl;
    } else {
        if (value < 0.0) {
            cerr << "Global variable ev_max_charging_power_kW cannot be set to value " << value << " ( allowed range: [0.0,inf[ )" << endl;
            return;
        }
        Global::ev_max_charging_power_kW = value;
    }
}
void Global::set_ev_charging_effi(float value) {
    if (is_locked) {
        cerr << "Global variable ev_charging_effi cannot be set!" << endl;
    } else {
        if (value <= 0.0f || value > 1.0f) {
            cerr << "Global variable ev_charging_effi cannot be set to value " << value << " ( allowed range: ]0.0,1.0] )" << endl;
            return;
        }
        Global::ev_charging_effi = value;
    }
}
void Global::set_cs_max_charging_power_kW(float value) {
    if (is_locked) {
        cerr << "Global variable cs_max_charging_power_kW cannot be set!" << endl;
    } else {
        if (value <= 0.0) {
            cerr << "Global variable cs_max_charging_power_kW cannot be set to value " << value << " ( allowed range: ]0.0,inf[ )" << endl;
            return;
        }
        Global::cs_max_charging_power_kW = value;
    }
}
void Global::set_use_emission_time_series_ia(bool use) {
    if (is_locked) {
        cerr << "Global variable use_emission_time_series_ia cannot be set!" << endl;
    } else {
        use_emission_time_series_ia = use;
    }
}
void Global::set_use_prices_time_series_ia(bool use) {
    if (is_locked) {
        cerr << "Global variable use_emission_time_series_ia cannot be set!" << endl;
    } else {
        use_prices_time_series_ia = use;
    }
}
void Global::set_select_only_residential_buildings(bool value) {
    if (is_locked) {
        cerr << "Global variable select_only_residential_buildings cannot be set!" << endl;
    } else {
        select_only_residential_buildings = value;
    }
}
void Global::set_control_horizon_in_ts(unsigned int value) {
    if (is_locked) {
        cerr << "Global variable control_horizon_in_ts cannot be set!" << endl;
    } else {
        if (value < 1) {
            cerr << "Global variable control_horizon_in_ts cannot be set to value " << value << " (allowed range: [1,inf) )" << endl;
            return;
        }
        control_horizon_in_ts = value;
    }
}
void Global::set_control_update_freq_in_ts(unsigned int value) {
    if (is_locked) {
        cerr << "Global variable control_update_freq_in_ts cannot be set!" << endl;
    } else {
        if (value < 1) {
            cerr << "Global variable control_update_freq_in_ts cannot be set to value " << value << " (allowed range: [1,inf) )" << endl;
            return;
        }
        control_update_freq_in_ts = value;
    }
}
void Global::set_input_path(string* path) {
    if (is_locked && input_path_init) {
        cerr << "Input path already set!" << endl;
    } else {
        Global::input_path = *path;
        Global::input_path_init = true;
        // check, if path ends with an "/", add it, if not
        if (Global::input_path.back() != '/') {
            Global::input_path += "/";
        }
    }
}
void Global::set_output_path(string* path) {
    if (is_locked && output_path_init) {
        cerr << "Output path already set!" << endl;
    } else {
        Global::output_path = *path;
        Global::output_path_init = true;
        // check, if path ends with an "/", add it, if not
        if (Global::output_path.back() != '/') {
            Global::output_path += "/";
        }
    }
}
void Global::set_structure_database_name(std::string* fname) {
    if (is_locked) {
        cerr << "System structure database name already set!" << endl;
    } else {
        Global::system_db_name = *fname;
    }
}
void Global::set_ev_data_path(string* path) {
    if (is_locked) {
        cerr << "EV data path is already set!" << endl;
    } else {
        Global::ev_data_path = *path;
        // check, if path ends with an "/", add it, if not
        if (Global::ev_data_path.back() != '/') {
            Global::ev_data_path += "/";
        }
    }
}
void Global::set_output_mode_per_cu(OutputModePerCU mode) {
    if (is_locked && output_mode_per_cu_init) {
        cerr << "Output mode already set!" << endl;
    } else {
        Global::output_mode_per_cu = mode;
        Global::output_mode_per_cu_init = true;
    }
}
void Global::set_exp_profile_mode(ExpansionProfileAllocationMode mode) {
    if (is_locked && exp_profile_mode_init) {
        cerr << "Expansion profile mode already set!" << endl;
    } else {
        Global::exp_profile_mode = mode;
        Global::exp_profile_mode_init = true;
    }
}
void Global::set_cu_selection_mode_fca(global::CUSModeFCA mode) {
    if (is_locked && cu_selection_mode_fca_init) {
        cerr << "CU selection mode for component addition already set!" << endl;
    } else {
        Global::cu_selection_mode_fca = mode;
        Global::cu_selection_mode_fca_init = true;
    }
}
void Global::set_battery_power_computation_mode(global::BatteryPowerComputationMode mode) {
    if (is_locked && bat_power_comp_mode_init) {
        cerr << "Battery power computation mode is already set!" << endl;
    } else {
        Global::bat_power_comp_mode = mode;
        Global::bat_power_comp_mode_init = true;
    }
}
void Global::set_battery_capacity_computation_mode(global::BatteryCapacityComputationMode mode) {
    if (is_locked) {
        cerr << "Battery capacity computation mode is already set!" << endl;
    } else {
        Global::bat_capacity_comp_mode = mode;
    }
}
void Global::set_controller_mode(global::ControllerMode mode) {
    if (is_locked) {
        cerr << "Controller mode cannot be set." << endl;
    } else {
        Global::controller_mode = mode;
    }
}
void Global::set_controller_optimization_target(global::ControllerOptimizationTarget mode) {
    if (is_locked) {
        cerr << "Variable controller_optimization_target cannot be set at the moment." << endl;
    } else {
        controller_optimization_target = mode;
    }
}
void Global::set_cu_list_for_sac_planning(const std::set<unsigned long>* selected_cuIDs) {
    Global::cu_list_for_sac_planning = selected_cuIDs;
}
void Global::set_controller_bs_grid_charging_mode(global::ControllerBSGridChargingMode mode) {
    if (is_locked) {
        cerr << "Variable controller_bs_grid_charging_mode cannot be set at the moment." << endl;
    } else {
        controller_bs_grid_charging_mode = mode;
    }
}
void Global::set_max_parallel_opti_vars(unsigned long value) {
    if (is_locked) {
        cerr << "Variable max_parallel_opti_vars cannot be set at the moment." << endl;
    } else {
        max_parallel_opti_vars = value;
    }
}
void Global::set_annual_heat_demand_limit_fsac(float value) {
    if (is_locked) {
        cerr << "Variables cannot be set currently!" << endl;
    } else {
        if (value < 1.0 && value != -1) {
            cerr << "Global variable annual_heat_demand_limit_fsac cannot be set to value " << value << " (allowed range: -1 or [1,inf) )" << endl;
            return;
        }
        Global::annual_heat_demand_limit_fsac = value;
    }
}
void Global::set_select_buildings_wg_heatd_only(bool value) {
    if (is_locked) {
        cerr << "Variables cannot be set currently!" << endl;
    } else {
        Global::select_buildings_wg_heatd_only = value;
    }
}
void Global::set_create_substation_output(bool value) {
    if (is_locked) {
        cerr << "Variables cannot be set currently!" << endl;
    } else {
        Global::create_substation_output = value;
    }
}
void Global::set_create_control_cmd_output(bool value) {
    if (is_locked) {
        cerr << "Variables cannot be set currently!" << endl;
    } else {
        Global::create_control_cmd_output = value;
    }
}
void Global::set_create_ev_detailed_output(bool value) {
    if (is_locked) {
        cerr << "Variables cannot be set currently!" << endl;
    } else {
        Global::create_ev_detailed_output = value;
    }
}
void Global::set_exp_pv_static_profile_orientation(std::string* value) {
    if (is_locked) {
        cerr << "Variables cannot be set currently!" << endl;
    } else {
        Global::exp_pv_static_profile_orientation = *value;
    }
}
void Global::set_exp_pv_static_profile_idx(int value) {
    if (is_locked) {
        cerr << "Variables cannot be set currently!" << endl;
    } else {
        if (value < 0) {
            cerr << "Global variable exp_pv_static_profile_idx cannot be set to value " << value << " (allowed range: [0,inf) )" << endl;
            return;
        }
        Global::exp_pv_static_profile_idx = value;
    }
}

