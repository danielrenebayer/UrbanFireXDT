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
    
    std::cout << "SurplusController: Executing global optimization for timestep " << ts_horizon_start <<  " to " << ts_horizon_end << std::endl;

    const std::vector<ControlUnit*>& units = ControlUnit::GetArrayOfInstances();

    // A: Fetch nececarry data
    std::unordered_map<unsigned long, std::vector<double>> grid_demand_kWh; ///< Grid demand per unit ID for all timesteps in horizon
    std::unordered_map<unsigned long, std::vector<double>> bs_stored_energy_kWh; ///< BESS stored energy per unit ID for all timesteps in horizon. Assume initial SoC 0%
    std::unordered_map<unsigned long, std::vector<double>> bs_power_kW; ///< BESS power per unit ID for all timesteps in horizon
    std::unordered_map<unsigned long, double> bs_capacity_kWh; ///< BESS total capacity per unit ID
    std::unordered_map<unsigned long, double> bs_max_charge_kWh; ///< BESS max charge energy in one timestep per unit ID
    std::vector<unsigned long> units_with_bs; ///< List of unit IDs that have a BESS

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
        
        // Reset charge requests
        unit_charge_requests[unit_id].assign(optimization_frequency_ts, 0.0);
    }

    // Lookahead simulation
    auto lookahead_result = LookaheadSimulation(ts_horizon_start, ts_horizon_end);
    std::vector<double> total_load = lookahead_result.future_surplus;
    future_surplus_log = total_load;
    if(bess_knowledge){
        bs_stored_energy_kWh = lookahead_result.bs_stored_energy_kWh;
        bs_power_kW = lookahead_result.bs_power_kW;
    }
    grid_demand_kWh = lookahead_result.grid_demand_kWh;

    // B: Calculate surplus allocation
    for (unsigned long ts = ts_horizon_start; ts <= ts_horizon_end; ++ts) {
        if (total_load[ts - ts_horizon_start] >= 0.0) {
            continue; // No surplus to allocate in this timestep, continue with next timestep
        }
        // Randomize units_with_bs order for fairness
        static std::mt19937 rng(Global::is_seed_set() ? Global::get_seed() : std::random_device{}());
        std::shuffle(units_with_bs.begin(), units_with_bs.end(), rng);
        for (auto unit_id : units_with_bs) {
            if (total_load[ts - ts_horizon_start] >= 0.0) {
                break; // Surplus fully allocated, continue with next timestep
            }
            auto current_left_capacity = bs_capacity_kWh[unit_id] - bs_stored_energy_kWh[unit_id][ts - ts_horizon_start]; 
            double future_demand = 0.0;
            for (unsigned long future_ts = ts+1; future_ts <= ts_horizon_end; ++future_ts) {
            // for (unsigned long future_ts = ts; future_ts <= ts_horizon_end; ++future_ts) {
                future_demand += grid_demand_kWh[unit_id][future_ts - ts_horizon_start];
            }
            // Determine how much can be charged this timestep based on all constraints
            double left_power = bs_max_charge_kWh[unit_id] - bs_power_kW[unit_id][ts - ts_horizon_start] * Global::get_time_step_size_in_h();
            double charge_possible = std::min({future_demand, current_left_capacity, left_power, -total_load[ts - ts_horizon_start]});
            // Update data structures
            if (charge_possible > 0.0) {
                // Save charge request if in this control period
                if (ts <= (ts_horizon_start + optimization_frequency_ts - 1)) {
                    unit_charge_requests[unit_id][ts - ts_horizon_start] += charge_possible;
                }
                bs_stored_energy_kWh[unit_id][ts - ts_horizon_start] += charge_possible;
                bs_power_kW[unit_id][ts - ts_horizon_start] += charge_possible / Global::get_time_step_size_in_h();
                total_load[ts - ts_horizon_start] += charge_possible;
                // propagate to grid demand and storage in the following hours of the day
                for (unsigned long future_ts = ts+1; future_ts <= ts_horizon_end; ++future_ts) {
                // for (unsigned long future_ts = ts; future_ts <= ts_horizon_end; ++future_ts) {
                    if (charge_possible <= 0.0) {
                        break; // stop when all stored energy is allocated to later grid demand
                    }
                    double left_power = bs_max_charge_kWh[unit_id] + bs_power_kW[unit_id][future_ts - ts_horizon_start] * Global::get_time_step_size_in_h();
                    double discharge_possible = std::min({charge_possible, grid_demand_kWh[unit_id][future_ts - ts_horizon_start], left_power});
                    grid_demand_kWh[unit_id][future_ts - ts_horizon_start] -= discharge_possible;
                    charge_possible -= discharge_possible;
                    bs_stored_energy_kWh[unit_id][future_ts - ts_horizon_start] += charge_possible;
                    bs_power_kW[unit_id][future_ts - ts_horizon_start] -= discharge_possible / Global::get_time_step_size_in_h();
                    total_load[future_ts - ts_horizon_start] -= discharge_possible;
                }
            }
        }
    }

    // Save optimization timestep for keeping track when to run next
    last_optimization_ts = ts_horizon_start;
    return true;
}

double SurplusController::GetChargeRequest(unsigned long unit_id) const {
    auto it = unit_charge_requests.find(unit_id);
    if (it == unit_charge_requests.end() || it->second.empty()) {
        return 0.0; // No charge request available
    }

    return it->second[0];
}

void SurplusController::ShiftTimeSeriesData() {
    // Shift all time series data by one step
    for (auto& pair : unit_charge_requests) {
        auto& time_series = pair.second;
        
        // Shift data and pad with zeros
        std::move(time_series.begin() + 1, time_series.end(), time_series.begin());
        
        // Fill last position with zero
        time_series.back() = 0.0;
    }
}

void SurplusController::ResetAllData() {
    unit_charge_requests.clear();
    last_optimization_ts = 0;
}

LookaheadResult SurplusController::LookaheadSimulation(unsigned long ts_horizon_start, unsigned long ts_horizon_end) {
    // Save current state
    auto state = ControlUnit::SaveAllInternalStates();    

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
        
        total_load[ts - ts_horizon_start] = total_load_kW;

        // Extract control unit specific parameters
        const std::vector<ControlUnit*>& units = ControlUnit::GetArrayOfInstances();
        for (ControlUnit* unit : units) {
            // unit specific data only needed for units with simulated BESS
            if(unit->has_bs_sim_added()){
                auto unit_id = unit->get_unitID();
                // grid demand
                auto load = unit->get_current_load_vSMeter_kW() * Global::get_time_step_size_in_h();
                // Only positive demand counts toward grid demand
                grid_demand_kWh[unit_id].push_back((load > 0.0) ? load : 0.0);
                // BESS stored energy
                auto bs = unit->get_component_BS();
                bs_stored_energy_kWh[unit_id].push_back(bs->get_currentCharge_kWh());
                // BESS power
                bs_power_kW[unit_id].push_back(bs->get_currentLoad_kW());
            }
        }
    }

    // Restore state
    ControlUnit::RestoreAllInternalStates(state);

    // Return surplus forecast and grid demand
    return {total_load, grid_demand_kWh, bs_stored_energy_kWh, bs_power_kW};
};

// Static convenience methods for ControlUnit access
double SurplusController::GetChargeRequestForUnit(unsigned long unit_id) {
    return GetInstance().GetChargeRequest(unit_id);
}

double SurplusController::GetScheduledSurplusToBESS() {
    auto& instance = GetInstance();
    if (!instance.enabled) {
        return 0.0;
    }
    double total_surplus = 0.0;
    for (const auto& pair : instance.unit_charge_requests) {
        if (!pair.second.empty()) {
            total_surplus += pair.second[0];
        }
    }
    return total_surplus;
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

double SurplusController::GetFutureSurplusLog(unsigned long ts){
    auto& instance = GetInstance();
    if (!instance.enabled) {
        return 0.0;
    }
    return instance.future_surplus_log[ts-instance.last_optimization_ts];
}

}// namespace surplus