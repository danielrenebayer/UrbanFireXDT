#include "surplus_controller.hpp"
#include "units.h"
#include "global.h"
#include <iostream>
#include <algorithm>
#include <random>
#include "simulation_logic.h"

namespace surplus {

// Static member initialization
std::unique_ptr<SurplusController> SurplusController::instance = nullptr;
bool SurplusController::initialized = false;

SurplusController::SurplusController() 
    : last_optimization_ts(0)
    , optimization_frequency_ts(Global::get_surplus_controller_frequency_ts())
    , lookahead_horizon_ts(Global::get_surplus_controller_lookahead_horizon_ts())
    , enabled(Global::get_surplus_controller_enabled())
    , bess_knowledge(Global::get_surplus_controller_BESS_knowledge())
    , thread_manager(nullptr)
    , output_prefix("") {
}

SurplusController& SurplusController::GetInstance() {
    if (!instance) {
        instance = std::unique_ptr<SurplusController>(new SurplusController());
    }
    return *instance;
}

void SurplusController::Initialize(CUControllerThreadGroupManager* thread_manager, const char* output_prefix) {
    if (!initialized) {
        auto& instance = GetInstance();
        instance.thread_manager = thread_manager;
        instance.output_prefix = output_prefix;
        instance.ResetAllData();
        initialized = true;
    }
}

void SurplusController::Cleanup() {
    if (instance) {
        instance.reset();
        initialized = false;
    }
}

bool SurplusController::ShouldRunOptimization(unsigned long current_ts) const {
    if (!enabled) {
        return false;
    }
    
    // Run optimization at the beginning and then every optimization_frequency_ts steps
    if (last_optimization_ts == 0) {
        return true;
    }
    
    return (current_ts - last_optimization_ts) >= optimization_frequency_ts;
}

bool SurplusController::ExecuteOptimization(unsigned long ts_horizon_start) {
    if (!enabled) {
        return true; 
    }

    unsigned long ts_horizon_end = ts_horizon_start + lookahead_horizon_ts - 1;
    if (ts_horizon_end > Global::get_last_timestep()) {
        ts_horizon_end = Global::get_last_timestep();
    }
    unsigned long ts_optimization_end = ts_horizon_start + optimization_frequency_ts - 1;
    if (ts_optimization_end > ts_horizon_end) {
        ts_optimization_end = ts_horizon_end;
    }
    
    std::cout << "SurplusController: Executing global optimization for timestep " << ts_horizon_start <<  " to " << ts_horizon_end << std::endl;

    const std::vector<ControlUnit*>& units = ControlUnit::GetArrayOfInstances();

    // A: Fetch nececarry data
    std::unordered_map<unsigned long, std::vector<double>> grid_demand_kWh; ///< Grid demand per unit ID for all timesteps in horizon
    std::unordered_map<unsigned long, std::vector<double>> bs_stored_energy_kWh; ///< BESS stored energy per unit ID for all timesteps in horizon. Assume initial SoC 0%
    std::unordered_map<unsigned long, std::vector<double>> bs_power_kW; ///< BESS power per unit ID for all timesteps in horizon
    std::unordered_map<unsigned long, double> bs_capacity_kWh; ///< BESS total capacity per unit ID
    std::unordered_map<unsigned long, double> bs_max_charge_kWh; ///< BESS max charge energy in one timestep per unit ID
    std::vector<unsigned long> units_with_bs; ///< List of unit IDs that have a BESS

    double self_discharge_rate = Global::get_exp_bess_self_ds_ts();
    double efficiency_in = Global::get_exp_bess_effi_in();
    double efficiency_out = Global::get_exp_bess_effi_out();

    for (ControlUnit* unit : units) {

        // only process if unit has simulated BESS 
        if (!unit->has_bs_sim_added()){
            continue;
        }
        // add unit to list
        auto unit_id = unit->get_unitID();
        units_with_bs.push_back(unit_id);
        
        // Collect static BESS data
        if(!bess_knowledge){
            // If no BESS knowledge, assume we start at 0% SoC and always have 0 kW power
            bs_stored_energy_kWh[unit_id] = std::vector<double>(ts_horizon_end - ts_horizon_start + 1, 0.0);
            bs_power_kW[unit_id] = std::vector<double>(ts_horizon_end - ts_horizon_start + 1, 0.0);
        }
        if(unit->has_bs_sim_added()){
            bs_capacity_kWh[unit_id] = unit->get_sim_comp_bs_E_kWh();
            bs_max_charge_kWh[unit_id] = unit->get_sim_comp_bs_P_kW() * Global::get_time_step_size_in_h();
        }
        
        // Initialize request vectors if not yet done, or resize if length changed (e.g., due to parameter variation)
        if (unit_charge_requests_kW[unit_id].empty() || unit_charge_requests_kW[unit_id].size() != optimization_frequency_ts) {
            unit_charge_requests_kW[unit_id].resize(optimization_frequency_ts, 0.0);
        }
        // Discharge requests need full horizon length (to keep future discharge plans from previous optimizations)
        if (unit_discharge_requests_kW[unit_id].empty() || unit_discharge_requests_kW[unit_id].size() != lookahead_horizon_ts) {
            unit_discharge_requests_kW[unit_id].resize(lookahead_horizon_ts, 0.0);
        }
    }

    // Lookahead simulation
    auto lookahead_result = LookaheadSimulation(ts_horizon_start, ts_horizon_end);
    std::vector<double> total_load_kWh = lookahead_result.total_load;
    future_surplus_log = total_load_kWh; // Log for analysis
    if(bess_knowledge){
        bs_stored_energy_kWh = lookahead_result.bs_stored_energy_kWh;
        bs_power_kW = lookahead_result.bs_power_kW;
    }
    grid_demand_kWh = lookahead_result.grid_demand_kWh;

    // B: Calculate surplus allocation
    for (unsigned long ts = ts_horizon_start; ts <= ts_optimization_end; ++ts) {
        if (total_load_kWh[ts - ts_horizon_start] >= 0.0) {
            continue; // No surplus to allocate in this timestep, continue with next timestep
        }
        
        // Iteratively allocate surplus energy to battery for replacing future demand 
        // Iterate future timesteps first to prioritize closer demands
        for (unsigned long future_ts = ts+1; future_ts <= ts_horizon_end; ++future_ts) {
            
            // Break if all surplus has been consumed
            if (total_load_kWh[ts - ts_horizon_start] >= 0.0) {
                break; // No more surplus to allocate in this timestep
            }
            
            // Skip if this future timestep is already in surplus
            if (total_load_kWh[future_ts - ts_horizon_start] <= 0.0) {
                continue; // this future timestep is already in surplus
            }
            
            // Randomize units_with_bs order for fairness at each future timestep
            static std::mt19937 rng(Global::is_seed_set() ? Global::get_seed() : std::random_device{}());
            std::shuffle(units_with_bs.begin(), units_with_bs.end(), rng);
            
            for (auto unit_id : units_with_bs) {

                // Breaking conditions for allocation loop (End this timestep - no more surplus)
                if (total_load_kWh[ts - ts_horizon_start] >= 0.0) {
                    break; // No more surplus to allocate in this timestep
                }
                if (total_load_kWh[future_ts - ts_horizon_start] <= 0.0) {
                    break; // this future timestep is already in surplus
                }

                // Continue conditions (Skip this unit for this future_ts)
                if (bs_power_kW[unit_id][ts - ts_horizon_start] >= bs_max_charge_kWh[unit_id] / Global::get_time_step_size_in_h()) {
                    continue; // No more power capacity left to charge in this timestep
                }
                if (grid_demand_kWh[unit_id][future_ts - ts_horizon_start] == 0.0) {
                    continue; // no demand to prevent at this future timestep
                }

                // Grid demand this future timestep (what should be replaced by charging now)
                // This var is further limited to all constraints. It always represents energy in kWh on grid side!
                double demand_this_ts = grid_demand_kWh[unit_id][future_ts - ts_horizon_start];

                // Limit demand to available BESS power at discharging timestep
                // The unit should not be both drawing from grid and charging BESS at the same time (would be grid charging, which is not allowed)
                if (bs_power_kW[unit_id][future_ts - ts_horizon_start] > 0.0) {
                    // auto debug_bs_power = bs_power_kW[unit_id][future_ts - ts_horizon_start];
                    // auto debug_grid_demand = grid_demand_kWh[unit_id][future_ts - ts_horizon_start];
                    std::cerr << "SurplusController Optimization Error: Unit " << unit_id << " has grid demand and charging BESS at ts " << future_ts << std::endl;
                    return false;
                }
                double max_demand_from_discharge_power = bs_max_charge_kWh[unit_id] + bs_power_kW[unit_id][future_ts - ts_horizon_start] * Global::get_time_step_size_in_h();
                demand_this_ts = std::min(demand_this_ts, max_demand_from_discharge_power);

                // Consider self-discharge: battery must store more to compensate for self-discharge over time
                demand_this_ts = demand_this_ts / pow(1-self_discharge_rate,future_ts - ts);

                // This accounts for losses during charging (efficiency_in) and discharging (efficiency_out)
                // Efficiency losses are only considered for what is charged from grid into the battery. If the battery initially discharges, no efficiency loss is considered for the part that just prevents this discharge.
                // demand_this_ts is now slittet into the two parts demand_prevent_discharge_kWh and demand_stored_kWh seperately, both in grid energy terms. 
                // Both energy levels are from now on stored into this BESS (either charged or prevented discharge), but only demand_stored_kWh causes efficiency losses.
                double demand_prevent_discharge_kWh = 0.0;
                double demand_stored_kWh = 0.0;
                if (bs_power_kW[unit_id][ts - ts_horizon_start] < 0.0) {
                    demand_prevent_discharge_kWh = std::min((-bs_power_kW[unit_id][ts - ts_horizon_start] * Global::get_time_step_size_in_h()), demand_this_ts);
                    demand_stored_kWh = demand_this_ts - demand_prevent_discharge_kWh;
                }
                else{
                    demand_stored_kWh = demand_this_ts;
                }
                // Efficiency losses only apply to the part that is actually charged from grid
                demand_stored_kWh = demand_stored_kWh / (efficiency_in * efficiency_out);

                // Limit demand to available BESS power at charging timestep, subtract if currently discharging
                demand_stored_kWh = std::min({demand_stored_kWh, bs_max_charge_kWh[unit_id], bs_max_charge_kWh[unit_id] - bs_power_kW[unit_id][ts - ts_horizon_start] * Global::get_time_step_size_in_h()});

                // Limit demand to current surplus
                demand_prevent_discharge_kWh = std::min(demand_prevent_discharge_kWh, -total_load_kWh[ts - ts_horizon_start]);
                demand_stored_kWh = std::min(demand_stored_kWh, -total_load_kWh[ts - ts_horizon_start] - demand_prevent_discharge_kWh);

                // Limit demand, so that no surplus is created / increased at future timestep
                double future_total_load_if_full_demand = total_load_kWh[future_ts - ts_horizon_start] - (demand_prevent_discharge_kWh + demand_stored_kWh * efficiency_in * efficiency_out) * pow(1-self_discharge_rate, future_ts - ts);
                if (future_total_load_if_full_demand < 0.0) {
                    double excess_energy_kWh = -future_total_load_if_full_demand;
                    // Reduce demand_stored first
                    double stored_energy_at_grid_kWh = demand_stored_kWh * efficiency_in * efficiency_out * pow(1-self_discharge_rate, future_ts - ts);
                    if (stored_energy_at_grid_kWh >= excess_energy_kWh) {
                        // Can handle by reducing stored only
                        double reduced_stored_energy_at_grid_kWh = stored_energy_at_grid_kWh - excess_energy_kWh;
                        demand_stored_kWh = reduced_stored_energy_at_grid_kWh / (efficiency_in * efficiency_out * pow(1-self_discharge_rate, future_ts - ts));
                    } else {
                        // Need to reduce both: stored to zero, and reduce prevent_discharge
                        demand_stored_kWh = 0.0;
                        double remaining_excess_kWh = excess_energy_kWh - stored_energy_at_grid_kWh;
                        demand_prevent_discharge_kWh = std::max(0.0, demand_prevent_discharge_kWh - remaining_excess_kWh / pow(1-self_discharge_rate, future_ts - ts));
                    }
                }

                // Limit demand to available capacity in BESS at all timesteps between charging and discharging
                for (unsigned long intermediate_ts = ts; intermediate_ts < future_ts; ++intermediate_ts) {
                    double intermediate_left_capacity_kWh = bs_capacity_kWh[unit_id] - bs_stored_energy_kWh[unit_id][intermediate_ts - ts_horizon_start];
                    // Calculate battery energy from both parts: prevent needs /eff_out (stays in battery), stored needs *eff_in (charged from grid)
                    double battery_energy_at_ts = demand_prevent_discharge_kWh / efficiency_out + demand_stored_kWh * efficiency_in;
                    double energy_at_intermediate = battery_energy_at_ts * pow(1-self_discharge_rate, intermediate_ts - ts);
                    
                    if (energy_at_intermediate > intermediate_left_capacity_kWh) {
                        // Reduce demand to fit capacity - prioritize reducing newely stored energy
                        double max_battery_energy_at_ts = intermediate_left_capacity_kWh / pow(1-self_discharge_rate, intermediate_ts - ts);
                        double excess_energy = battery_energy_at_ts - max_battery_energy_at_ts;
                        
                        // First reduce demand_stored 
                        double stored_battery_energy = demand_stored_kWh * efficiency_in;
                        if (stored_battery_energy >= excess_energy) {
                            // Can handle by reducing stored only
                            demand_stored_kWh = (stored_battery_energy - excess_energy) / efficiency_in;
                        } else {
                            // Need to reduce both: stored to zero, and reduce prevent_discharge
                            demand_stored_kWh = 0.0;
                            demand_prevent_discharge_kWh = max_battery_energy_at_ts * efficiency_out;
                        }
                        demand_this_ts = demand_prevent_discharge_kWh + demand_stored_kWh;
                    }
                }

                // Sum back to total demand of this future timestep
                demand_this_ts = demand_prevent_discharge_kWh + demand_stored_kWh;

                // Update data structures to propagate changes
                if (demand_this_ts > 0.0) {
                    // Save charge request if in this control period (convert from kWh to kW)
                    if (ts <= ts_optimization_end) {
                        unit_charge_requests_kW[unit_id][ts - ts_horizon_start] += demand_this_ts / Global::get_time_step_size_in_h();
                    }
                    // Storage energy update: compute battery energy from both parts
                    double battery_energy_charged_kWh = demand_prevent_discharge_kWh / efficiency_out + demand_stored_kWh * efficiency_in;
                    for (unsigned long intermediate_ts = ts; intermediate_ts < future_ts; ++intermediate_ts) {
                        // Apply self-discharge to battery energy
                        double energy_at_intermediate_ts = battery_energy_charged_kWh * pow(1-self_discharge_rate, intermediate_ts - ts);
                        bs_stored_energy_kWh[unit_id][intermediate_ts - ts_horizon_start] += energy_at_intermediate_ts;
                    }
                    
                    // Calculate what comes out of battery at discharge timestep (with self-discharge and discharge efficiency)
                    double energy_discharged_to_grid_kWh = battery_energy_charged_kWh * pow(1-self_discharge_rate, future_ts - ts) * efficiency_out;

                    // Save discharge request for the full horizon (convert from kWh to kW)
                    unit_discharge_requests_kW[unit_id][future_ts - ts_horizon_start] += energy_discharged_to_grid_kWh / Global::get_time_step_size_in_h();
                    
                    // Storage power update at charging timestep (grid side)
                    bs_power_kW[unit_id][ts - ts_horizon_start] += demand_this_ts / Global::get_time_step_size_in_h();
                    // Storage power update at discharging timestep (grid side)
                    bs_power_kW[unit_id][future_ts - ts_horizon_start] -= energy_discharged_to_grid_kWh / Global::get_time_step_size_in_h();
                    // Update total load at charging timestep
                    total_load_kWh[ts - ts_horizon_start] += demand_this_ts;
                    // Update total load at discharging timestep
                    total_load_kWh[future_ts - ts_horizon_start] -= energy_discharged_to_grid_kWh;   
                    // Update grid demand at discharging timestep
                    grid_demand_kWh[unit_id][future_ts - ts_horizon_start] -= energy_discharged_to_grid_kWh;                
                }
            }
        }
    }

    // Validate if everything is within limits, TODO: Can propabely be removed later, just for debugging purposes
    for (auto unit_id : units_with_bs) {

        // auto debug_power = bs_power_kW[unit_id];
        // auto debug_storage = bs_stored_energy_kWh[unit_id];
        // auto debug_max_power = bs_max_charge_kWh[unit_id];
        // auto debug_capacity = bs_capacity_kWh[unit_id];

        for (unsigned long ts = ts_horizon_start; ts <= ts_horizon_end; ++ts) {
            // Check power limits
            if (bs_power_kW[unit_id][ts - ts_horizon_start] > bs_max_charge_kWh[unit_id] / Global::get_time_step_size_in_h() + 1e-6) {
                std::cerr << "SurplusController Optimization Error: Exceeded max charge power for unit " << unit_id << " at timestep " << ts << std::endl;
                return false;
            }
            if (bs_power_kW[unit_id][ts - ts_horizon_start] < -bs_max_charge_kWh[unit_id] / Global::get_time_step_size_in_h() - 1e-6) {
                std::cerr << "SurplusController Optimization Error: Exceeded max discharge power for unit " << unit_id << " at timestep " << ts << std::endl;
                return false;
            }
            // Check capacity limits
            if (bs_stored_energy_kWh[unit_id][ts - ts_horizon_start] > bs_capacity_kWh[unit_id] + 1e-6) {
                std::cerr << "SurplusController Optimization Error: Exceeded BESS capacity for unit " << unit_id << " at timestep " << ts << std::endl;
                return false;
            }
            if (bs_stored_energy_kWh[unit_id][ts - ts_horizon_start] < -1e-6) {
                std::cerr << "SurplusController Optimization Error: Below zero BESS storage for unit " << unit_id << " at timestep " << ts << std::endl;
                return false;
            }
        }
    }

    // ... more validation
    for (unsigned long ts = ts_horizon_start; ts <= ts_optimization_end; ++ts) {
        // Check that a timestep either has only charge requests or only discharge requests over all units
        bool has_charge = false;
        bool has_discharge = false;
        for (auto unit_id : units_with_bs) {
            if (unit_charge_requests_kW[unit_id][ts - ts_horizon_start] > 1e-6) {
                has_charge = true;
            }
            if (unit_discharge_requests_kW[unit_id][ts - ts_horizon_start] > 1e-6) {
                has_discharge = true;
            }
        }
        // Allow simultaneous charge and discharge only if no BESS knowledge and we are beyond optimization horizon 
        // Otherwise in non-BESS-knowledge mode we might have a charge/discharge request from previous optimization run, but have to rule against it in the current optimization because the Units didn't behave as expected
        if (has_charge && has_discharge) {
            if(bess_knowledge || ts > ts_optimization_end){
                std::cerr << "SurplusController Optimization Error: Both charge and discharge request at timestep " << ts << std::endl;
                return false;
            }
            else{
                std::cout << "SurplusController Warning: Both charge and discharge request at timestep " << ts << " in no-BESS-knowledge mode." << std::endl;
            }
        }
    }

    // Save optimization timestep for keeping track when to run next
    last_optimization_ts = ts_horizon_start;
    return true;
}

double SurplusController::GetChargeRequest(unsigned long unit_id) const {
    auto it = unit_charge_requests_kW.find(unit_id);
    if (it == unit_charge_requests_kW.end() || it->second.empty()) {
        return 0.0; // No charge request available
    }

    return it->second[0];
}

double SurplusController::GetDischargeRequest(unsigned long unit_id) const {
    auto it = unit_discharge_requests_kW.find(unit_id);
    if (it == unit_discharge_requests_kW.end() || it->second.empty()) {
        return 0.0; // No discharge request available
    }

    return it->second[0];
}

void SurplusController::ShiftTimeSeriesData() {
    // Shift all time series data by one step
    for (auto& pair : unit_charge_requests_kW) {
        auto& time_series = pair.second;
        
        // Shift data and pad with zeros
        std::move(time_series.begin() + 1, time_series.end(), time_series.begin());
        
        // Fill last position with zero
        if (!time_series.empty()) {
            time_series.back() = 0.0;
        }
    }
    for (auto& pair : unit_discharge_requests_kW) {
        auto& time_series = pair.second;
        
        // Shift data and pad with zeros
        std::move(time_series.begin() + 1, time_series.end(), time_series.begin());
        
        // Fill last position with zero
        if (!time_series.empty()) {
            time_series.back() = 0.0;
        }
    }
}

void SurplusController::ResetAllData() {
    unit_charge_requests_kW.clear();
    unit_discharge_requests_kW.clear();
    last_optimization_ts = 0;
    
    // Refresh configuration parameters from Global (important for parameter variations)
    optimization_frequency_ts = Global::get_surplus_controller_frequency_ts();
    lookahead_horizon_ts = Global::get_surplus_controller_lookahead_horizon_ts();
    enabled = Global::get_surplus_controller_enabled();
    bess_knowledge = Global::get_surplus_controller_BESS_knowledge();
}

LookaheadResult SurplusController::LookaheadSimulation(unsigned long ts_horizon_start, unsigned long ts_horizon_end) {
    // Save current state
    auto state = ControlUnit::SaveAllInternalStates();    
    auto discharge_requests_copy = unit_discharge_requests_kW;  // Save discharge requests, as they might extend into the next lookahead horizon if (horizon > optimization frequency)


    // Pre-run the simulation for the length of the optimization horizon to get forecasted surplus and loads of each control unit
    std::vector<double> total_load(ts_horizon_end - ts_horizon_start + 1, 0.0);
    std::unordered_map<unsigned long, std::vector<double>> grid_demand_kWh;
    std::unordered_map<unsigned long, std::vector<double>> bs_stored_energy_kWh;
    std::unordered_map<unsigned long, std::vector<double>> bs_power_kW;
    const double totalBatteryCapacity_kWh = ControlUnit::GetAllSimCompBatteriesCapacity_kWh();
    
    for (unsigned long ts = ts_horizon_start; ts <= ts_horizon_end; ++ts) {
        double total_load_kW = 0.0;
        // Run one step
        simulation::oneStep(ts, totalBatteryCapacity_kWh, thread_manager, output_prefix, NULL, &total_load_kW, false);
        
        total_load[ts - ts_horizon_start] = total_load_kW * Global::get_time_step_size_in_h();

        // Extract control unit specific parameters
        const std::vector<ControlUnit*>& units = ControlUnit::GetArrayOfInstances();
        for (ControlUnit* unit : units) {
            // unit specific data only needed for units with simulated BESS
            if(unit->has_bs_sim_added()){
                auto unit_id = unit->get_unitID();
                // grid demand
                auto load = unit->get_current_load_vSMeter_kW() * Global::get_time_step_size_in_h();
                // Only positive demand counts toward grid demand
                grid_demand_kWh[unit_id].push_back((load > 1e-6) ? load : 0.0);
                // BESS stored energy
                auto bs = unit->get_component_BS();
                bs_stored_energy_kWh[unit_id].push_back(bs->get_currentCharge_kWh());
                // BESS power
                bs_power_kW[unit_id].push_back(bs->get_currentLoad_kW());
            }
        }

        // Shift Discharge Requests in surplus controller, as there might be discharge requests scheduled from previous optimizations
        ShiftTimeSeriesData();
    }

    // Restore state
    ControlUnit::RestoreAllInternalStates(state);
    unit_discharge_requests_kW = discharge_requests_copy;

    // Return surplus forecast and grid demand
    return {total_load, grid_demand_kWh, bs_stored_energy_kWh, bs_power_kW};
};

// Static convenience methods for ControlUnit access
double SurplusController::GetChargeRequestForUnit(unsigned long unit_id) {
    return GetInstance().GetChargeRequest(unit_id);
}

double SurplusController::GetDischargeRequestForUnit(unsigned long unit_id) {
    return GetInstance().GetDischargeRequest(unit_id);
}

double SurplusController::GetScheduledSurplusToBESS() {
    auto& instance = GetInstance();
    if (!instance.enabled) {
        return 0.0;
    }
    double total_surplus = 0.0;
    for (const auto& pair : instance.unit_charge_requests_kW) {
        if (!pair.second.empty()) {
            total_surplus += pair.second[0];
        }
    }
    return total_surplus;
}

double SurplusController::GetScheduledSurplusToUnit() {
    auto& instance = GetInstance();
    if (!instance.enabled) {
        return 0.0;
    }
    double total_discharge = 0.0;
    for (const auto& pair : instance.unit_discharge_requests_kW) {
        if (!pair.second.empty()) {
            total_discharge += pair.second[0];
        }
    }
    return total_discharge;
}

double SurplusController::GetActualSurplusToBESS() {
    auto& instance = GetInstance();
    if (!instance.enabled) {
        return 0.0;
    }
    double total_actual_surplus = 0.0;
    const std::vector<ControlUnit*>& units = ControlUnit::GetArrayOfInstances();
    for (ControlUnit* unit : units) {
        if (!unit->has_bs_sim_added()){
            continue;
        }
        auto bs = unit->get_component_BS();
        double actual_surplus_charge = bs->get_currentLoad_kW() - unit->get_initial_charge_request_kW();
        total_actual_surplus += actual_surplus_charge;
    }
    return total_actual_surplus;
}

double SurplusController::GetBESSChargeRequest() {
    double sum = 0.0;
    const std::vector<ControlUnit*>& units = ControlUnit::GetArrayOfInstances();
    for (ControlUnit* unit : units) {
        if (!unit->has_bs_sim_added()){
            continue;
        }
        auto bs = unit->get_component_BS();
        sum += bs->get_chargeRequest();
    }
    return sum;
}

double SurplusController::GetBESSLoad() {
    double sum = 0.0;
    const std::vector<ControlUnit*>& units = ControlUnit::GetArrayOfInstances();
    for (ControlUnit* unit : units) {
        if (!unit->has_bs_sim_added()){
            continue;
        }
        auto bs = unit->get_component_BS();
        sum += bs->get_currentLoad_kW();
    }
    return sum;
}

double SurplusController::GetBESSSurplusEnergy() {
    double sum = 0.0;
    const std::vector<ControlUnit*>& units = ControlUnit::GetArrayOfInstances();
    for (ControlUnit* unit : units) {
        if (!unit->has_bs_sim_added()){
            continue;
        }
        auto bs = unit->get_component_BS();
        sum += bs->get_currentE_from_surplus();
    }
    return sum;
}

double SurplusController::GetFutureSurplusLog(unsigned long ts){
    auto& instance = GetInstance();
    if (!instance.enabled) {
        return 0.0;
    }
    return instance.future_surplus_log[ts-instance.last_optimization_ts];
}

}// namespace surplus