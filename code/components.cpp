
#include "components.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <ranges>
#include <string>
#include <vector>

#include "cache_helper.hpp"
#include "vehicles.h"
#include "global.h"
#include "output.h"


// ----------------------------- //
//      Implementation of        //
//   BaseComponentSemiFlexible   //
// ----------------------------- //

BaseComponentSemiFlexible::BaseComponentSemiFlexible() {
    // initialization of the cached vector
    set_horizon_in_ts( Global::get_control_horizon_in_ts() );
}

BaseComponentSemiFlexible::~BaseComponentSemiFlexible() {
}

void BaseComponentSemiFlexible::set_horizon_in_ts(unsigned int new_horizon) {
    future_maxE_storage.clear();
    future_minE_storage.clear();

    future_maxE_storage.resize(new_horizon, 0.0);
    future_minE_storage.resize(new_horizon, 0.0);
}



// ----------------------------- //
//      Implementation of        //
//        RoofSectionPV          //
// ----------------------------- //

std::map<std::string, size_t> RoofSectionPV::next_pv_idx;
std::map<std::string, std::default_random_engine>            RoofSectionPV::random_generators;
std::map<std::string, std::uniform_int_distribution<size_t>> RoofSectionPV::distributions;

RoofSectionPV::RoofSectionPV(size_t locationID, float this_section_kWp, const std::string& orientation)
    : this_section_kWp(this_section_kWp), orientation(orientation)
{
    // 1)
    // take the correct profile
    profile_index = 0;
    if (Global::get_exp_pv_sizing_mode() == global::PVSizingMode::StaticPVSize &&
        Global::get_exp_pv_static_profile_orientation() != "" &&
        Global::get_exp_pv_static_profile_idx() >= 0)
    {
        profile_index = Global::get_exp_pv_static_profile_idx();
    }
    else if (Global::get_exp_profile_mode() == global::ExpansionProfileAllocationMode::AsInData) {
        // case: selection as it appears in data
        // look up next index for a given orientation
        // if this orientation is not in the map up to now, take the last available index and add this orientation to the map
        if (next_pv_idx.contains(orientation)) {
            size_t& next_pv_idx_so = next_pv_idx[orientation];
            profile_index = next_pv_idx_so;
            next_pv_idx_so++; // increment next index by one until end is reached
            if (next_pv_idx_so >= global::pv_profiles_information[orientation])
                next_pv_idx_so = 0;
        } else {
            profile_index = global::pv_profiles_information[orientation] - 1; // take last index available
            next_pv_idx[orientation] = 0; // set this to first, as this is the following on to the last one
        }
    } else {
        // case random selection
        std::optional<size_t> cached_profile_idx = PVProfileIDCache::GetInstance().readCache(locationID, orientation);
        if (cached_profile_idx.has_value()) {
            //If there is a cache file for the profile allocation -> use it!
            profile_index = *cached_profile_idx;
        } else {
            // Otherwise, sample a new number
            if (!random_generators.contains(orientation)) {
                random_generators[orientation] = std::default_random_engine{};
                if (Global::is_seed_set()) random_generators[orientation].seed(Global::get_seed());
                distributions[orientation] = std::uniform_int_distribution<size_t>(0, global::pv_profiles_information[orientation]-1);
            }
            // randomly select new index
            profile_index = distributions[orientation](random_generators[orientation]);
            // update cache
            PVProfileIDCache::GetInstance().updateCache(locationID, orientation, profile_index);
        }
    }
    // 2)
    // select the correct profile
    if (global::pv_profiles_per_ori[orientation].size() <= 0) {
        std::cerr << "Error: There is no feed-in profile given for the orientation " << orientation << std::endl;
        throw runtime_error("There is no feed-in profile given for a selected orientation!");
    }
    if (global::pv_profiles_per_ori[orientation].size() <= profile_index) {
        std::cerr << "Error: There is no feed-in profile given for the orientation " << orientation << " and index " << profile_index << std::endl;
        throw runtime_error("There is no feed-in profile given for a selected orientation and the given index!");
    }
    // 3) initialize the pointer to the profile data
    profile_data = global::pv_profiles_per_ori.at(orientation)[profile_index];
}

float RoofSectionPV::get_generation_at_ts_kW(unsigned long ts) const {
    unsigned long tsID = ts - 1;
    return this_section_kWp * profile_data[tsID];
}

void RoofSectionPV::update_section_kWp(float new_kWp) {
    if (new_kWp < 0.0f) new_kWp = 0.0f;
    this_section_kWp = new_kWp;
}





// ----------------------------- //
//      Implementation of        //
//         ComponentPV           //
// ----------------------------- //

ComponentPV::ComponentPV(float kWp, unsigned long locationID)
    : total_kWp(kWp)
{
    /*
     *
     * Constructor for ComponentPV
     * in the case of a static kWp computation
     *
     */
    currentGeneration_kW = 0;
    generation_cumsum_total_kWh = 0.0;
    generation_cumsum_cweek_kWh = 0.0;

  if (Global::get_exp_pv_static_profile_orientation() == "") {
    // Case 1: if no static profile is selected, use existing roof data
    //
    // attach roof sections as defined in data
    auto& roof_section_vec = global::roof_section_orientations[locationID];
    //unsigned long number_of_sections = roof_section_vec.size();
    float complete_roof_area = accumulate(
                begin(roof_section_vec),
                end(roof_section_vec),
                0.0f,
                [](const float prev, const auto& elem){ return prev + elem.first; });;
    // iterate over all roof sections
    // and get all orientations and areas per roof section
    for (auto& section_tuple : roof_section_vec) {
        float  section_roof_area   = section_tuple.first;
        std::string section_orientation = section_tuple.second;
        // 1)
        // compute the share of this section among all existing sections
        float share_of_total_area = section_roof_area / complete_roof_area;
        // 2)
        // create and add section to list
        float section_kWp = share_of_total_area * kWp;
        roof_sections.emplace_back(locationID, section_kWp, section_orientation);
    }
  } else {
    // Case 2: Only use one component facing the given orientation
    string ori = Global::get_exp_pv_static_profile_orientation();
    roof_sections.emplace_back(locationID, kWp, ori);
  }

    // initialization of the cached vector
    set_horizon_in_ts( Global::get_control_horizon_in_ts() );
}

ComponentPV::ComponentPV(float kWp_per_m2, float min_kWp_sec, float max_kWp_sec, float max_kWp_unit, unsigned long locationID)
{
    /*
     *
     * Constructor for ComponentPV
     * in the case of a dynamic kWp computation
     *
     */

    // general checks
    if (kWp_per_m2 <= 0.0) {
        throw runtime_error("Error in ComponentPV constructor: kWp_per_m2 <= 0!");
    }
    if (max_kWp_sec > 0 && min_kWp_sec > max_kWp_sec) {
        throw runtime_error("Error in ComponentPV constructor: min_kWp_sec > max_kWp_sec");
    }

    currentGeneration_kW = 0;
    generation_cumsum_total_kWh = 0.0;
    generation_cumsum_cweek_kWh = 0.0;
    total_kWp = 0;
    // attach roof sections as defined in data
    // // float min_roof_area = min_kWp / kWp_per_m2;
    // // float max_roof_area = max_kWp / kWp_per_m2;
    auto& roof_section_vec = global::roof_section_orientations[locationID];
    std::vector<std::pair<float,std::string>> vec_of_sections; // helper vector required to check, if Global::exp_pv_max_kWp_per_unit is reached
    // iterate over all roof sections
    // 1)
    // get all orientations and areas per roof section
    for (auto& section_tuple : roof_section_vec) {
        float       section_roof_area   = section_tuple.first;
        std::string section_orientation = section_tuple.second;
        // 1)
        // compute the kWp of this area and decide whether to use it or not
        float section_kWp = section_roof_area * kWp_per_m2;
        if (section_kWp < min_kWp_sec) {
            // ignore this section
            continue;
        }
        // if max_kWp == 0 -> ignore this value
        if (max_kWp_sec > 0 && section_kWp > max_kWp_sec) {
            section_kWp = max_kWp_sec;
        }
        total_kWp += section_kWp;
        // 2)
        // add section to list
        vec_of_sections.emplace_back(section_kWp, section_orientation);
    }
    // 2)
    // only in the case of Global::exp_pv_max_kWp_per_unit is set (i.e. >= 0.0)
    if (max_kWp_unit >= 0.0) {
        // is total_kWp > max_kWp_unit (which is most of the time = Global::get_exp_pv_max_kWp_per_unit(), except of parameter variations)
        // if yes, lower installed power with the same percantage over all sections
        if ( total_kWp > max_kWp_unit ) {
            float reduction_factor = max_kWp_unit / total_kWp;
            total_kWp = max_kWp_unit;
            for (auto& section_tuple : vec_of_sections) {
                section_tuple.first = section_tuple.first * reduction_factor;
            }
        }
    }
    // 3)
    // instanziate objects of roof sections
    for (auto& section_tuple : vec_of_sections) {
        float       section_kWp         = section_tuple.first;
        std::string section_orientation = section_tuple.second;
        roof_sections.emplace_back(locationID, section_kWp, section_orientation);
    }

    // initialization of the cached vector
    set_horizon_in_ts( Global::get_control_horizon_in_ts() );
}

float ComponentPV::get_generation_at_ts_kW(unsigned long ts) const {
    if (ts <= 0 || ts > Global::get_n_timesteps()) {
        return 0.0;
    }
    float generation_kW = 0.0;
    for (const RoofSectionPV& section : roof_sections)
        generation_kW += section.get_generation_at_ts_kW(ts);
    return generation_kW;
}

string* ComponentPV::get_section_string(const string& prefix_per_line) {
    string* return_str = new string();
    unsigned long counter = 1;
    string local_line_prefix { prefix_per_line };
    if (local_line_prefix.length() > 0)
        local_line_prefix.append(",");
    for (RoofSectionPV& section : roof_sections) {
        return_str->append( local_line_prefix );
        return_str->append( to_string( counter ) );
        return_str->append(",");
        return_str->append( to_string( section.get_section_kWp() )   );
        return_str->append(",");
        return_str->append(            section.get_orientation()     );
        return_str->append(",");
        return_str->append( to_string( section.get_profile_index() ) );
        return_str->append("\n");
        counter++;
    }
    return return_str;
}

std::list<std::vector<double>> ComponentPV::get_total_generation_by_section_kW() const {
    const unsigned long ts_start = Global::get_first_timestep();
    const unsigned long ts_end   = Global::get_last_timestep();
    const unsigned long horizon_len = ts_end - ts_start + 1;
    std::list<std::vector<double>> output_list;

    for (const RoofSectionPV& section : roof_sections) {
        std::vector<double>& this_sec_vec = output_list.emplace_back();
        this_sec_vec.reserve(horizon_len);

        for (unsigned long ts = ts_start; ts <= ts_end; ts++) {
            this_sec_vec.push_back( section.get_generation_at_ts_kW(ts) );
        }
    }

    return output_list;
}

std::list<double> ComponentPV::get_kWp_per_section() const {
    std::list<double> output_list;

    for (const RoofSectionPV& section : roof_sections) {
        output_list.push_back( section.get_section_kWp() );
    }

    return output_list;
}

void ComponentPV::calculateCurrentFeedin(unsigned long ts) {
    currentGeneration_kW = 0.0;
    for (RoofSectionPV& section : roof_sections)
        currentGeneration_kW += section.get_generation_at_ts_kW(ts);
    // compute total generation
    double e = Global::get_time_step_size_in_h() * currentGeneration_kW;
    generation_cumsum_total_kWh += e;
    generation_cumsum_cweek_kWh += e;
    // update future generation vector
    for (size_t tOffset = 0; tOffset < Global::get_control_horizon_in_ts(); tOffset++) {
        future_generation_kW[tOffset] = get_generation_at_ts_kW(ts + tOffset); // get_generation_at_ts_kW() returns 0.0 if outside data
    }
}

void ComponentPV::resetWeeklyCounter() {
    generation_cumsum_cweek_kWh = 0.0;
}

void ComponentPV::resetInternalState() {
    generation_cumsum_total_kWh = 0.0;
    generation_cumsum_cweek_kWh = 0.0;
    currentGeneration_kW = 0.0;
}

ComponentStateVariant ComponentPV::saveInternalState() const {
    PVComponentState state;
    state.currentGeneration_kW = currentGeneration_kW;
    state.generation_cumsum_total_kWh = generation_cumsum_total_kWh;
    state.generation_cumsum_cweek_kWh = generation_cumsum_cweek_kWh;
    return state;
}

void ComponentPV::restoreInternalState(const ComponentStateVariant& state) {
    try {
        const PVComponentState& pv_state = std::get<PVComponentState>(state);
        currentGeneration_kW = pv_state.currentGeneration_kW;
        generation_cumsum_total_kWh = pv_state.generation_cumsum_total_kWh;
        generation_cumsum_cweek_kWh = pv_state.generation_cumsum_cweek_kWh;
    } catch (const std::bad_variant_access& e) {
        throw std::runtime_error("ComponentPV::restoreInternalState(): Invalid state type provided");
    }
}

void ComponentPV::set_horizon_in_ts(unsigned int new_horizon) {
    future_generation_kW.clear();
    future_generation_kW.resize(new_horizon, 0.0);
}

void ComponentPV::set_kWp_per_section(const std::list<double>& new_values) {
    total_kWp = 0.0;
    auto it = new_values.begin();
    for (RoofSectionPV& section : roof_sections) {
        // safety checks
        if (it == new_values.end()) {
            throw std::runtime_error("ComponentPV::set_kWp_per_section(): not enough values provided for all roof sections");
        }
        // set the new value
        float v = static_cast<float>(*it++);
        if (v < 0.0) v = 0.0;
        section.update_section_kWp(v);
        total_kWp += v;
    }
    // final safety checks
    if (it != new_values.end()) {
        throw std::runtime_error("ComponentPV::set_kWp_per_section(): too many values provided compared to roof sections");
    }
}




// ----------------------------- //
//      Implementation of        //
//         ComponentBS           //
// ----------------------------- //

ComponentBS::ComponentBS(
    double maxE_kWh,
    double maxP_kW,
    double E_over_P_ratio,
    double discharge_rate_per_step,
    double efficiency_in,
    double efficiency_out,
    double initial_SoC
) : maxE_kWh(maxE_kWh),
    discharge_rate_per_step(discharge_rate_per_step), efficiency_in(efficiency_in),
    efficiency_out(efficiency_out), initial_SoC(initial_SoC)
{
    SOC               = 0;
    currentE_kWh      = 0;
    currentP_kW       = 0;
    currentE_from_grid_kWh = 0.0;
    currentE_from_surplus_kWh = 0.0;
    currentP_from_grid_kW  = 0.0;
    charge_request_kW = 0;
    cweek_E_withdrawn_kWh = 0.0;
    total_E_withdrawn_kWh = 0.0;
    cweek_E_withdrawn_from_grid_kWh = 0.0;
    total_E_withdrawn_from_grid_kWh = 0.0;
    n_ts_SOC_empty    = 0;
    n_ts_SOC_full     = 0;

    if (Global::get_battery_power_computation_mode() == global::BatteryPowerComputationMode::UseEOverPRatio) {
        this->E_over_P_ratio = E_over_P_ratio;
        this->maxP_kW = maxE_kWh / E_over_P_ratio;
    } else {
        this->E_over_P_ratio = maxE_kWh / maxP_kW;
        this->maxP_kW = maxP_kW;
    }

    if (initial_SoC > 0) {
        SOC = initial_SoC;
        currentE_kWh = maxE_kWh * initial_SoC;
    }

}

ComponentBS::ComponentBS(
    double maxE_kWh,
    double discharge_rate_per_step,
    double efficiency_in,
    double initial_SoC
) : maxE_kWh(maxE_kWh),
    discharge_rate_per_step(discharge_rate_per_step),
    efficiency_in(efficiency_in),
    efficiency_out(1.0),
    initial_SoC(initial_SoC)
{
    SOC               = 0;
    currentE_kWh      = 0;
    currentP_kW       = 0;
    currentE_from_grid_kWh = 0.0;
    currentE_from_surplus_kWh = 0.0;
    currentP_from_grid_kW  = 0.0;
    charge_request_kW = 0;
    cweek_E_withdrawn_kWh = 0.0;
    total_E_withdrawn_kWh = 0.0;
    cweek_E_withdrawn_from_grid_kWh = 0.0;
    total_E_withdrawn_from_grid_kWh = 0.0;
    n_ts_SOC_empty    = 0;
    n_ts_SOC_full     = 0;

    maxP_kW = maxE_kWh / Global::get_time_step_size_in_h(); // This is required so that the battery can be emptied in one single step at the arrival

    if (initial_SoC > 0) {
        SOC = initial_SoC;
        currentE_kWh = maxE_kWh * initial_SoC;
    }
}

void ComponentBS::set_grid_charged_amount(double grid_charged_kW) {
    // a guard, that it is only called when the BS is charged
    if (currentP_kW > 0) {
#ifdef DEBUG
        // a check, if grid_charged_kW <= currentP_kW
        if (grid_charged_kW > currentP_kW) {
            throw std::runtime_error("Error in ComponentBS: grid_charged_kW > currentP_kW");
        }
#endif
        currentE_from_grid_kWh += grid_charged_kW * Global::get_time_step_size_in_h() * efficiency_in;
    }
}

void ComponentBS::set_surplus_charged_amount(double surplus_charged_kW) {
    currentE_from_surplus_kWh += surplus_charged_kW * Global::get_time_step_size_in_h() * efficiency_in;
    if (currentE_from_surplus_kWh > currentE_kWh) {
        currentE_from_surplus_kWh = currentE_kWh;
    }
}


void ComponentBS::reset_surplus_charged_amount() {
    currentE_from_surplus_kWh = 0.0;
}

void ComponentBS::set_SOE_without_computations(double new_SOE_kWh) {
    currentE_kWh = new_SOE_kWh;
}

void ComponentBS::set_maxE_kWh(double value) {
    maxE_kWh = value;
    if (Global::get_battery_power_computation_mode() == global::BatteryPowerComputationMode::UseEOverPRatio)
        this->maxP_kW = maxE_kWh / E_over_P_ratio;
}

void ComponentBS::set_maxP_kW(double value) {
    if (Global::get_battery_power_computation_mode() == global::BatteryPowerComputationMode::AsDefinedByConfigVar) {
        maxP_kW  = value;
        E_over_P_ratio = maxE_kWh / maxP_kW;
    }
}

void ComponentBS::set_maxP_by_EPRatio(double EP_ratio) {
    if (Global::get_battery_power_computation_mode() == global::BatteryPowerComputationMode::UseEOverPRatio) {
        this->E_over_P_ratio = EP_ratio;
        this->maxP_kW = this->maxE_kWh / EP_ratio;
    }
}

void ComponentBS::set_efficiency_in(double value) {
    efficiency_in = value;
}

void ComponentBS::set_efficiency_out(double value) {
    efficiency_out = value;
}

void ComponentBS::set_self_discharge_rate(double value) {
    discharge_rate_per_step = value;
}

double const ComponentBS::validateNoSurplusChargeRequest(double charge_request_kW) {
    double timestep_size_in_h = Global::get_time_step_size_in_h();

    double currentE_after_self_discharge_kWh = currentE_kWh - discharge_rate_per_step * currentE_kWh;
    double currentE_from_surplus_after_self_discharge_kWh = currentE_from_surplus_kWh - discharge_rate_per_step * currentE_from_surplus_kWh;

    if (charge_request_kW > 0) {
        // Charging: limit to maxP_kW
        if (charge_request_kW > maxP_kW)
            charge_request_kW = maxP_kW;
        
        // Check if charging would exceed battery capacity
        double new_charge_kWh = currentE_after_self_discharge_kWh + timestep_size_in_h * charge_request_kW * efficiency_in;
        if (new_charge_kWh > maxE_kWh) {
            // Adjust charge request to exactly fill the battery
            charge_request_kW = (maxE_kWh - currentE_after_self_discharge_kWh) / timestep_size_in_h / efficiency_in;
        }
    } else if (charge_request_kW < 0) {
        // Discharging: limit to maxP_kW
        if (-charge_request_kW > maxP_kW)
            charge_request_kW = -maxP_kW;
        
        // Calculate available non-surplus energy (only return charged energy from non-surplus sources)
        double available_non_surplus_energy_kWh = currentE_after_self_discharge_kWh - currentE_from_surplus_after_self_discharge_kWh;
        
        // Check if discharging would deplete non-surplus energy below zero
        double new_charge_kWh = currentE_after_self_discharge_kWh + timestep_size_in_h * charge_request_kW / efficiency_out;
        if (new_charge_kWh < currentE_from_surplus_after_self_discharge_kWh) {
            // Adjust discharge request to only use available non-surplus energy
            charge_request_kW = -available_non_surplus_energy_kWh / timestep_size_in_h * efficiency_out;
        }
    }

    return charge_request_kW;
}

void ComponentBS::calculateActions() {
    double timestep_size_in_h = Global::get_time_step_size_in_h();
    double new_charge_kWh;

    currentP_kW = 0;

    // Calculate Self-discharge
    currentE_kWh -= discharge_rate_per_step * currentE_kWh;
    currentE_from_grid_kWh -= discharge_rate_per_step * currentE_from_grid_kWh;
    currentE_from_surplus_kWh -= discharge_rate_per_step * currentE_from_surplus_kWh;

    // Charging and discharging
    if (charge_request_kW > 0) {
        // charging requested
        if (charge_request_kW > maxP_kW)
            charge_request_kW = maxP_kW;
        new_charge_kWh = currentE_kWh + timestep_size_in_h*charge_request_kW*efficiency_in;
        if (new_charge_kWh > maxE_kWh)
            new_charge_kWh = maxE_kWh;
        currentP_kW  = (new_charge_kWh - currentE_kWh)/timestep_size_in_h/efficiency_in;
        currentE_kWh = new_charge_kWh;
    } else if (charge_request_kW < 0) {
        // discharging requested
        if (-charge_request_kW > maxP_kW)
            charge_request_kW = -maxP_kW;
        new_charge_kWh = currentE_kWh + timestep_size_in_h*charge_request_kW/efficiency_out;
        if (new_charge_kWh < 0)
            new_charge_kWh = 0;
        double energy_taken_kWh = new_charge_kWh - currentE_kWh;
        currentP_kW  = energy_taken_kWh / timestep_size_in_h * efficiency_out;
        currentE_kWh = new_charge_kWh;
        // add withrawn energy to summation variable (mind energy_taken_kWh < 0)
        total_E_withdrawn_kWh -= energy_taken_kWh;
        cweek_E_withdrawn_kWh -= energy_taken_kWh;
        // update currentE_from_grid_kWh as well ...
        currentP_from_grid_kW = 0.0;
        if (currentE_kWh < currentE_from_grid_kWh) {
            double removed_amount_from_grid_kWh = currentE_from_grid_kWh - currentE_kWh;
            currentP_from_grid_kW = removed_amount_from_grid_kWh / timestep_size_in_h * efficiency_out;
            currentE_from_grid_kWh = currentE_kWh; // set to current SOE
            // add to summation variables
            cweek_E_withdrawn_from_grid_kWh += removed_amount_from_grid_kWh;
            total_E_withdrawn_from_grid_kWh += removed_amount_from_grid_kWh;
        }
        // ... and surplus energy as well
        if (currentE_kWh < currentE_from_surplus_kWh) {
            currentE_from_surplus_kWh = currentE_kWh; // set to current SOE
        }
    }

    // calculate new SOC value
    SOC = currentE_kWh / maxE_kWh;

    // calculate n_ts_SOC_empty or n_ts_SOC_full
    if (currentE_kWh <= 0)
        n_ts_SOC_empty++;
    else if (currentE_kWh >= maxE_kWh)
        n_ts_SOC_full++;
}

void ComponentBS::resetWeeklyCounter() {
    cweek_E_withdrawn_kWh = 0.0;
    cweek_E_withdrawn_from_grid_kWh = 0.0;
}

void ComponentBS::resetInternalState() {
    //
    // This method resets the internal SOC to the initial level
    //
    SOC = initial_SoC;
    currentE_kWh = maxE_kWh * initial_SoC;
    currentE_from_grid_kWh  = 0.0;
    currentP_kW           = 0.0;
    currentE_from_surplus_kWh = 0.0;
    currentP_from_grid_kW   = 0.0;
    cweek_E_withdrawn_kWh = 0.0;
    total_E_withdrawn_kWh = 0.0;
    cweek_E_withdrawn_from_grid_kWh = 0.0;
    total_E_withdrawn_from_grid_kWh = 0.0;
    n_ts_SOC_empty = 0;
    n_ts_SOC_full  = 0;
}

ComponentStateVariant ComponentBS::saveInternalState() const {
    BSComponentState state;
    state.SOC = SOC;
    state.currentE_kWh = currentE_kWh;
    state.currentP_kW = currentP_kW;
    state.currentE_from_grid_kWh = currentE_from_grid_kWh;
    state.currentE_from_surplus_kWh = currentE_from_surplus_kWh;
    state.currentP_from_grid_kW = currentP_from_grid_kW;
    state.cweek_E_withdrawn_kWh = cweek_E_withdrawn_kWh;
    state.total_E_withdrawn_kWh = total_E_withdrawn_kWh;
    state.cweek_E_withdrawn_from_grid_kWh = cweek_E_withdrawn_from_grid_kWh;
    state.total_E_withdrawn_from_grid_kWh = total_E_withdrawn_from_grid_kWh;
    state.n_ts_SOC_empty = n_ts_SOC_empty;
    state.n_ts_SOC_full = n_ts_SOC_full;
    return state;
}

void ComponentBS::restoreInternalState(const ComponentStateVariant& state) {
    try {
        const BSComponentState& bs_state = std::get<BSComponentState>(state);
        SOC = bs_state.SOC;
        currentE_kWh = bs_state.currentE_kWh;
        currentP_kW = bs_state.currentP_kW;
        currentE_from_grid_kWh = bs_state.currentE_from_grid_kWh;
        currentE_from_surplus_kWh = bs_state.currentE_from_surplus_kWh;
        currentP_from_grid_kW = bs_state.currentP_from_grid_kW;
        cweek_E_withdrawn_kWh = bs_state.cweek_E_withdrawn_kWh;
        total_E_withdrawn_kWh = bs_state.total_E_withdrawn_kWh;
        cweek_E_withdrawn_from_grid_kWh = bs_state.cweek_E_withdrawn_from_grid_kWh;
        total_E_withdrawn_from_grid_kWh = bs_state.total_E_withdrawn_from_grid_kWh;
        n_ts_SOC_empty = bs_state.n_ts_SOC_empty;
        n_ts_SOC_full = bs_state.n_ts_SOC_full;
    } catch (const std::bad_variant_access& e) {
        throw std::runtime_error("ComponentBS::restoreInternalState(): Invalid state type provided");
    }
}





// ----------------------------- //
//      Implementation of        //
//         ComponentHP           //
// ----------------------------- //

size_t ComponentHP::next_hp_idx = 0;
bool   ComponentHP::random_generator_init = false;
std::default_random_engine*            ComponentHP::random_generator = NULL;
std::uniform_int_distribution<size_t>* ComponentHP::distribution     = NULL;

ComponentHP::ComponentHP(const ControlUnit* connected_unit, float annual_econs_kWh)
    : connected_unit(connected_unit),
      yearly_electricity_consumption_kWh(annual_econs_kWh),
      scaling_factor(annual_econs_kWh/1000.0f)
{
    // initialize the unshiftable load storage
    future_maxP_storage.clear();
    future_minP_storage.clear();
    future_maxP_storage.resize( Global::get_control_horizon_in_ts(), 0.0);
    future_minP_storage.resize( Global::get_control_horizon_in_ts(), 0.0);

    // select heat pump profile static or random
    size_t this_hp_profile_idx;
    if (Global::get_exp_profile_mode() == global::ExpansionProfileAllocationMode::AsInData) {
        this_hp_profile_idx = next_hp_idx;
        // increment next index by one
        next_hp_idx++;
        if (next_hp_idx >= Global::get_n_heatpump_profiles())
            next_hp_idx = 0;
    } else {
        std::optional<size_t> cached_profile_idx = HPProfileIDCache::GetInstance().readCache(connected_unit->get_unitID());
        if (cached_profile_idx.has_value()) {
            //If there is a cache file for the profile allocation -> use it!
            this_hp_profile_idx = *cached_profile_idx;
        } else {
            // Otherwise, sample a new number
            if (!random_generator_init)
                ComponentHP::InitializeRandomGenerator();
            // randomly select new index
            this_hp_profile_idx = (*distribution)(*random_generator);
            // update cache
            HPProfileIDCache::GetInstance().updateCache(connected_unit->get_unitID(), this_hp_profile_idx);
        }
    }
    //
    // reference the profile
    profile_data = global::hp_profiles[this_hp_profile_idx];
    profile_cumsum = global::hp_profiles_cumsum[this_hp_profile_idx];
    // further initialization
    currentDemand_kW = 0;
    total_consumption_kWh = 0.0;
    cweek_consumption_kWh = 0.0;
    // computation of rated power (without AUX heating mode)
    // = max of shiftable demand time series
    // TODO: Update computation: Detect AUX heating mode
    rated_power_kW = 0.0;
    for (unsigned long tsID = 0; tsID < Global::get_n_timesteps(); tsID++) {
        float np = profile_data[tsID] * scaling_factor;
        if (np > rated_power_kW)
            rated_power_kW = np;
    }
#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    state_s1 = false;
#endif
}

void ComponentHP::set_horizon_in_ts(unsigned int new_horizon) {
    BaseComponentSemiFlexible::set_horizon_in_ts(new_horizon);
    future_maxP_storage.clear();
    future_minP_storage.clear();
    future_maxP_storage.resize(new_horizon, 0.0);
    future_minP_storage.resize(new_horizon, 0.0);
}

void ComponentHP::computeNextInternalState(unsigned long ts) {
#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    if (!state_s1) {
        state_s1 = true;
    } else {
        throw std::runtime_error("Method ComponentHP::computeNextInternalState() cannot be called at the moment! Call ComponentHP::setDemandToProfileData() or ComponentHP::setDemandToGivenValue() first.");
    }
#endif
    //
    unsigned long tsID = ts - 1;
    // compute variables for future min/max consumption and information on not-shiftable consumption
    double last_known_maxE_val = 0.0;
    double last_known_minE_val = 0.0;
    for (size_t tOffset = 0; tOffset < Global::get_control_horizon_in_ts(); tOffset++) {
        // power
        // TODO: detect AUX heating and increase hp min and max power in these steps
        future_maxP_storage[tOffset] = rated_power_kW;
        future_minP_storage[tOffset] = 0.0;
        // energy
        // for min profile -> right shift
        if (tsID + tOffset < Global::get_hp_flexibility_in_ts()) {
            future_minE_storage[tOffset] = 0.0;
        } else {
            size_t minProfilePos = tsID + tOffset - Global::get_hp_flexibility_in_ts();
            if (minProfilePos < Global::get_n_timesteps()) {
                double new_val = profile_cumsum[minProfilePos] * scaling_factor - total_consumption_kWh;
                if (new_val < 0)
                    new_val = 0.0;
                future_minE_storage[tOffset] = new_val;
                last_known_minE_val = new_val;
                // only required when reaching end of simulation (i.e., maxProfilePos >= Global::get_n_timesteps())
                // check if last_known_minE_val <= last_known_maxE_val holds
                if (last_known_maxE_val < last_known_minE_val)
                    last_known_maxE_val = last_known_minE_val;
            } else {
                future_minE_storage[tOffset] = last_known_minE_val;
            }
        }
        // for max profile -> left shift
        size_t maxProfilePos = tsID + tOffset + Global::get_hp_flexibility_in_ts();
        if (maxProfilePos < Global::get_n_timesteps() && tOffset < Global::get_control_horizon_in_ts() - Global::get_hp_flexibility_in_ts()) {
            double new_val = profile_cumsum[maxProfilePos] * scaling_factor - total_consumption_kWh;
            if (new_val < 0)
                new_val = 0.0;
            future_maxE_storage[tOffset] = new_val;
            last_known_maxE_val = new_val;
        } else {
            future_maxE_storage[tOffset] = last_known_maxE_val;
        }
    }
    // in the end, min and max profiles must come to the same point
    unsigned long last_tOffset = Global::get_control_horizon_in_ts() - 1;
    if (last_tOffset + ts >= Global::get_n_timesteps()) {
        future_minE_storage[last_tOffset] = future_maxE_storage[last_tOffset];
    }
}

void ComponentHP::setDemandToProfileData(unsigned long ts) {
#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    if (state_s1) {
        state_s1 = false;
    } else {
        throw std::runtime_error("Method ComponentHP::setDemandToProfileData() cannot be called at the moment! Call ComponentHP::computeNextInternalState() first.");
    }
#endif
    unsigned long tsID = ts - 1;
    currentDemand_kW = profile_data[tsID] * scaling_factor;
    double e = currentDemand_kW * Global::get_time_step_size_in_h();
    total_consumption_kWh += e;
    cweek_consumption_kWh += e;
}

#define epsilon_hp 0.0001
bool ComponentHP::setDemandToGivenValue(double new_demand_kW) {
#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    if (state_s1) {
        state_s1 = false;
    } else {
        throw std::runtime_error("Method ComponentHP::setDemandToGivenValue() cannot be called at the moment! Call ComponentHP::computeNextInternalState() first.");
    }
#endif
    bool no_error = true;
    double e = new_demand_kW * Global::get_time_step_size_in_h(); // amount of energy consumed in this time step
    //double new_total_e = total_consumption_kWh + e;
    // check, if power is within the min/max power bands
    // Output warnings only, if demand is below bound + epsilon
    if (new_demand_kW > future_maxP_storage[0]) {
        if (new_demand_kW > future_maxP_storage[0] + epsilon_hp) {
            no_error = false;
            std::cerr << "Warning in Component HP: Set demand is violating the upper power bound!\n";
            std::cerr << "     new P = " << std::fixed << std::setprecision(3) << new_demand_kW << ", but upper limit = " << std::fixed << std::setprecision(3) << future_maxP_storage[0] << "\n";
        }
        //
        new_demand_kW = future_maxP_storage[0];
        //
    } else if (new_demand_kW < future_minP_storage[0]) {
        if (new_demand_kW < future_minP_storage[0] - epsilon_hp) {
            no_error = false;
            std::cerr << "Warning in Component HP: Set demand is violating the lower power bound!\n";
            std::cerr << "     new P = " << std::fixed << std::setprecision(3) << new_demand_kW << ", but lower limit = " << std::fixed << std::setprecision(3) << future_minP_storage[0] << "\n";
        }
        //
        new_demand_kW = future_minP_storage[0];
    }
    // check, if the new demand is within the min/max energy consumption bands
    if (e > future_maxE_storage[0]) {
        if (e > future_maxE_storage[0] + epsilon_hp) {
            no_error = false;
            std::cerr << "Warning in Component HP: Set demand (" << std::fixed << std::setprecision(3) << new_demand_kW << " kW) is violating the upper energy consumption bound!\n";
            std::cerr << "     new E = " << std::fixed << std::setprecision(3) << e << ", but upper limit = " << std::fixed << std::setprecision(3) << future_maxE_storage[0] << "\n";
        }
        // set new_total_e to the upper limit
        e = future_maxE_storage[0];
        new_demand_kW = e / Global::get_time_step_size_in_h();
        //std::cerr << "Setting the new demand to " << std::fixed << std::setprecision(3) << new_demand_kW << " kW.\n";
    } else if (e < future_minE_storage[0]) {
        if (e < future_minE_storage[0] - epsilon_hp) {
            no_error = false;
            std::cerr << "Warning in Component HP: Set demand (" << std::fixed << std::setprecision(3) << new_demand_kW << " kW) is violating the lower energy consumption bound!\n";
            std::cerr << "     new E = " << std::fixed << std::setprecision(3) << e << ", but lower limit = " << std::fixed << std::setprecision(3) << future_minE_storage[0] << "\n";
        }
        // set new_total_e to the lower limit
        e = future_minE_storage[0];
        new_demand_kW = e / Global::get_time_step_size_in_h();
        //std::cerr << "Setting the new demand to " << std::fixed << std::setprecision(3) << new_demand_kW << " kW.\n";
    }
    // update current demand
    currentDemand_kW = new_demand_kW;
    // update cumulative variables
    total_consumption_kWh += e;
    cweek_consumption_kWh += e;

    return no_error;
}

void ComponentHP::resetWeeklyCounter() {
    cweek_consumption_kWh = 0.0;
}

void ComponentHP::resetInternalState() {
    currentDemand_kW = 0.0;
    total_consumption_kWh = 0.0;
    cweek_consumption_kWh = 0.0;
}

ComponentStateVariant ComponentHP::saveInternalState() const {
    HPComponentState state;
    state.currentDemand_kW = currentDemand_kW;
    state.total_consumption_kWh = total_consumption_kWh;
    state.cweek_consumption_kWh = cweek_consumption_kWh;
    return state;
}

void ComponentHP::restoreInternalState(const ComponentStateVariant& state) {
    try {
        const HPComponentState& hp_state = std::get<HPComponentState>(state);
        currentDemand_kW = hp_state.currentDemand_kW;
        total_consumption_kWh = hp_state.total_consumption_kWh;
        cweek_consumption_kWh = hp_state.cweek_consumption_kWh;
    } catch (const std::bad_variant_access& e) {
        throw std::runtime_error("ComponentHP::restoreInternalState(): Invalid state type provided");
    }
}

void ComponentHP::InitializeRandomGenerator() {
    if (random_generator_init) {
        throw runtime_error("Error: Static variables of ComponentHP are already initialized.");
    } else {
        random_generator_init = true;
        random_generator = new std::default_random_engine();
        if (Global::is_seed_set()) {
            random_generator->seed(Global::get_seed());
        }
        distribution     = new std::uniform_int_distribution<size_t>(0, Global::get_n_heatpump_profiles()-1);
    }
}

void ComponentHP::VacuumStaticVariables() {
    if (random_generator_init) {
        random_generator_init = false;
        delete random_generator;
        delete distribution;
    }
}




// ----------------------------- //
//      Implementation of        //
//         ComponentCS           //
// ----------------------------- //

ComponentCS::ComponentCS(ControlUnit* calling_control_unit, unsigned int number_of_flats) :
    installation_place(calling_control_unit), n_chargers(number_of_flats)
{
    enabled = false;
    current_demand_kW = 0.0;
    total_consumption_kWh  = 0.0;
    cweek_consumption_kWh  = 0.0;

    // Computation of maximum charging power
    max_charging_power = 11.0 * number_of_flats;
    if (Global::get_cs_max_charging_power_kW() > 0.0 &&
        Global::get_cs_max_charging_power_kW() < max_charging_power
    ) {
        max_charging_power = Global::get_cs_max_charging_power_kW();
    }

#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    is_callable_setCarStatesForTimeStep = true;
    is_callable_set_charging_value      = false;
#endif
}

ComponentCS::~ComponentCS() {
    for (EVFSM* ev : listOfEVs)
        delete ev;
}

double ComponentCS::get_max_P_kW() const {
    if (enabled)
        return max_charging_power;
    return 0.0;
}

unsigned long ComponentCS::get_n_EVs_pc() const {
    return std::ranges::count_if(listOfEVs, [](EVFSM* ev){return ev->get_current_state() == EVState::ConnectedAtHome;});
}

unsigned long ComponentCS::get_n_EVs_pnc() const {
    return std::ranges::count_if(listOfEVs.begin(), listOfEVs.end(), [](EVFSM* ev){return ev->get_current_state() == EVState::DisconnectedAtHome;});
}

unsigned long ComponentCS::get_n_EVs() const {
    if (enabled)
        return listOfEVs.size();
    return 0;
}

unsigned long ComponentCS::get_control_unit_id() const {
    return installation_place->get_unitID();
}

unsigned long ComponentCS::get_possible_n_EVs() const {
    return listOfEVs.size();
}

void ComponentCS::enable_station() {
    enabled = true;
}

void ComponentCS::disable_station() {
    enabled = false;
}

void ComponentCS::resetWeeklyCounter() {
    cweek_consumption_kWh = 0.0;
}

void ComponentCS::resetInternalState() {
    current_demand_kW = 0.0;
    total_consumption_kWh      = 0.0;
    cweek_consumption_kWh      = 0.0;
    //charging_order_req.clear();
    //charging_order_pos.clear();
    //
    for (EVFSM* ev : listOfEVs) {
        ev->resetInternalState();
    }
}

ComponentStateVariant ComponentCS::saveInternalState() const {
    CSComponentState state;
    state.current_demand_kW = current_demand_kW;
    state.total_consumption_kWh = total_consumption_kWh;
    state.cweek_consumption_kWh = cweek_consumption_kWh;
    
    // Save states of all connected EVs
    for (const EVFSM* ev : listOfEVs) {
        ComponentStateVariant ev_state_variant = ev->saveInternalState();
        EVFSMComponentState ev_state = std::get<EVFSMComponentState>(ev_state_variant);
        // Use the carID from the EV state as the map key
        unsigned long carID = ev_state.carID;
        state.ev_states_by_id[carID] = ev_state;
    }
    
    return state;
}

void ComponentCS::restoreInternalState(const ComponentStateVariant& state) {
    try {
        const CSComponentState& cs_state = std::get<CSComponentState>(state);
        current_demand_kW = cs_state.current_demand_kW;
        total_consumption_kWh = cs_state.total_consumption_kWh;
        cweek_consumption_kWh = cs_state.cweek_consumption_kWh;
        
        // Restore states of all connected EVs
        for (EVFSM* ev : listOfEVs) {
            unsigned long carID = ev->get_carID();
            
            auto it = cs_state.ev_states_by_id.find(carID);
            if (it != cs_state.ev_states_by_id.end()) {
                ComponentStateVariant ev_state_variant = it->second;
                ev->restoreInternalState(ev_state_variant);
            } else {
                throw std::runtime_error("ComponentCS::restoreInternalState(): No saved state found for carID " + std::to_string(carID));
            }
        }
        
    } catch (const std::bad_variant_access& e) {
        throw std::runtime_error("ComponentCS::restoreInternalState(): Invalid state type provided");
    }
}

void ComponentCS::add_ev(unsigned long carID) {
    listOfEVs.push_back(new EVFSM(carID, this));
    //
    // Reinitialize the vectors for future min/max storage
    future_minE_storage.assign(listOfEVs.size(), NULL);
    future_maxE_storage.assign(listOfEVs.size(), NULL);
    future_maxP_storage.assign(listOfEVs.size(), NULL);
}

void ComponentCS::set_horizon_in_ts(unsigned int new_horizon) {
    for (EVFSM* ev : listOfEVs) {
        ev->set_horizon_in_ts(new_horizon);
    }
}

void ComponentCS::preprocess_ev_data() {
    for (EVFSM* ev : listOfEVs) {
        ev->preprocessTourInformation();
    }
}

void ComponentCS::setCarStatesForTimeStep(unsigned long ts) {

#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    if (is_callable_setCarStatesForTimeStep) {
        is_callable_setCarStatesForTimeStep = false;
        is_callable_set_charging_value = true;
    } else {
        throw std::runtime_error("Method ComponentCS::setCarStatesForTimeStep() cannot be called at the moment!");
    }
#endif

    // 1. set new car states
    for (EVFSM* ev : listOfEVs) {
        ev->setCarStateForTimeStep(ts);
    }

    // 2. compute new vectors storing future min/max values
    for (unsigned long evIdx = 0; evIdx < listOfEVs.size(); evIdx++) {
        future_minE_storage[evIdx] = listOfEVs[evIdx]->get_future_min_consumption_kWh();
        future_maxE_storage[evIdx] = listOfEVs[evIdx]->get_future_max_consumption_kWh();
        future_maxP_storage[evIdx] = listOfEVs[evIdx]->get_future_max_power_kW();
    }
}

void ComponentCS::setDemandToProfileData(unsigned long ts) {

#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    if (is_callable_set_charging_value) {
        is_callable_set_charging_value = false;
        is_callable_setCarStatesForTimeStep = true;
    } else {
        throw std::runtime_error("Method ComponentCS::setDemandToProfileData() cannot be called at the moment!");
    }
#endif

    current_demand_kW = 0.0;
    // if no upper limit for charging power is set -> use EVFSM::setDemandToProfileData()
    // else use EVFSM::setDemandToGivenValue()
    if (Global::get_cs_max_charging_power_kW() <= 0.0) {
        for (EVFSM* ev : listOfEVs) {
            ev->setDemandToProfileData(ts);
            current_demand_kW += ev->get_currentDemand_kW();
        }
    } else {
        double remaining_power_kW = Global::get_cs_max_charging_power_kW();
        std::vector<double> new_power_per_EV_kW(listOfEVs.size(), 0.0);
        // Loop over all EVs and get min demand
        for (size_t ev_idx = 0; ev_idx < listOfEVs.size(); ev_idx++) {
            EVFSM* ev = listOfEVs[ev_idx];
            double min_demand = ev->get_future_min_consumption_kWh()->at(0) / Global::get_time_step_size_in_h(); // value at ts 0 always exists, as Global::control_horizon_in_ts >= 1 always holds!
            new_power_per_EV_kW[ev_idx] = min_demand;
            remaining_power_kW -= min_demand;
        }
        // Loop over all EVs again and assign max demand, if remaining power > 0
        for (size_t ev_idx = 0; ev_idx < listOfEVs.size(); ev_idx++) {
            EVFSM* ev = listOfEVs[ev_idx];
            double max_demand = ev->get_future_max_consumption_kWh()->at(0) / Global::get_time_step_size_in_h();
            if (remaining_power_kW > 0.0) {
                if (remaining_power_kW > max_demand) {
                    new_power_per_EV_kW[ev_idx] = max_demand;
                    remaining_power_kW -= max_demand;
                } else {
                    new_power_per_EV_kW[ev_idx] = remaining_power_kW;
                    remaining_power_kW = 0.0;
                    break;
                }
            } else {
                break;
            }
        }
        // Set actual power of all EVs
        for (size_t ev_idx = 0; ev_idx < listOfEVs.size(); ev_idx++) {
            EVFSM* ev = listOfEVs[ev_idx];
            ev->setDemandToGivenValue(new_power_per_EV_kW[ev_idx]);
            current_demand_kW += ev->get_currentDemand_kW();
        }
    }
    // compute current demand and cumulative sum
    double e = current_demand_kW * Global::get_time_step_size_in_h();
    total_consumption_kWh += e;
    cweek_consumption_kWh += e;

}

bool ComponentCS::setDemandToGivenValues(std::vector<double>& charging_power_per_EV_kW) {

#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    if (is_callable_set_charging_value) {
        is_callable_set_charging_value = false;
        is_callable_setCarStatesForTimeStep = true;
    } else {
        throw std::runtime_error("Method ComponentCS::setDemandToGivenValues() cannot be called at the moment!");
    }
#endif
    bool no_error = true;

    current_demand_kW = 0.0;
    for (size_t ev_idx = 0; ev_idx < listOfEVs.size() && ev_idx < charging_power_per_EV_kW.size(); ev_idx++) {
        bool ret_val = listOfEVs[ev_idx]->setDemandToGivenValue( charging_power_per_EV_kW[ev_idx] );
        current_demand_kW += listOfEVs[ev_idx]->get_currentDemand_kW();
        // Check if an 'error' occured, i.e., if the bounds are violated
        if (!ret_val) {
            no_error = false;
        }
    }
    // compute current demand and cumulative sum
    double e = current_demand_kW * Global::get_time_step_size_in_h();
    total_consumption_kWh += e;
    cweek_consumption_kWh += e;

    return no_error;
}




// ----------------------------- //
//      Implementation of        //
//            EVFSM              //
// ----------------------------- //

std::map<unsigned long, EVFSM*> EVFSM::list_of_cars;
std::default_random_engine*            EVFSM::random_generator = new std::default_random_engine();
std::uniform_real_distribution<float>* EVFSM::distribution     = new std::uniform_real_distribution<float>(0, 1);
const std::string EVFSM::MetricsStringHeaderAnnual = "CarID,HomeControlUnitID,Driving distance [km],E used for driving [kWh],Home-charged E [kWh],Home-discharged E [kWh],n ts home-connected";

EVFSM::EVFSM(unsigned long carID, ComponentCS* homeStation) :
    carID(carID),
    econs_kWh_per_km(Global::get_ev_consumption_kWh_km()),
    homeStation(homeStation)
{
    // Register this object in the global list of cars
    if (EVFSM::list_of_cars.contains(carID)) {
        throw runtime_error("Error: There is already an instance of class EVFSM with carID = " + std::to_string(carID));
    }
    EVFSM::list_of_cars.emplace(carID, this);
    //
    // Create battery (using the secondary constructor)
    battery = new ComponentBS(Global::get_ev_battery_size_kWh(), 0.0, 1.0, 1.0);
    // Initialize variables for the state
    current_state      = EVState::ConnectedAtHome;
    current_ts         = 0;
    //current_state_icah = EVStateIfConnAtHome::ChargingPossible;
    energy_demand_per_tour_ts = 0.0;
    current_P_kW              = 0.0;
    sum_of_driving_distance_km    = 0.0;
    sum_of_E_used_for_driving_kWh = 0.0;
    sum_of_E_charged_home_kWh     = 0.0;
    sum_of_E_discharged_home_kWh  = 0.0;
    sum_of_ts_EV_is_connected     = 0;
    // Create the list of tours
    list_of_tours_pd.reserve(6);
    for (unsigned int day = 0; day < 7; day++) {
        list_of_tours_pd.push_back( new std::vector<WeeklyVehicleTour> () );
    }
    // Initialize empty vector for storing future max power
    future_maxP_storage.assign(Global::get_control_horizon_in_ts(), 0.0);

#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    state_s1 = true;
    state_s2 = false;
    state_s3 = false;
#endif
}

EVFSM::~EVFSM() {
    delete battery;
    for (auto day_tours : list_of_tours_pd) {
        delete day_tours;
    }
    list_of_tours_pd.clear();
    list_of_all_tours.clear();
    // remove this car from the list_of_cars
    EVFSM::list_of_cars.erase(this->carID);
}

std::string* EVFSM::get_metrics_string_annual() {
    std::string* retstr = new string;
    *retstr += std::to_string(carID) + ",";
    *retstr += std::to_string(homeStation->get_control_unit_id()) + ",";
    *retstr += std::to_string(sum_of_driving_distance_km)     + ",";
    *retstr += std::to_string(sum_of_E_used_for_driving_kWh)  + ",";
    *retstr += std::to_string(sum_of_E_charged_home_kWh)      + ",";
    *retstr += std::to_string(sum_of_E_discharged_home_kWh)   + ",";
    *retstr += std::to_string(sum_of_ts_EV_is_connected);
    return retstr;
}

void EVFSM::set_horizon_in_ts(unsigned int new_horizon) {
    BaseComponentSemiFlexible::set_horizon_in_ts(new_horizon);
    future_maxP_storage.clear();
    future_maxP_storage.assign(new_horizon, 0.0);
}

void EVFSM::add_weekly_tour(
    short weekday,
    unsigned long departure_ts_of_day,
    unsigned long ts_duration,
    double tour_length_km,
    bool with_work)
{
#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    if (!state_s1) {
        throw std::runtime_error("Method EVFSM::add_weekly_tour() cannot be called after EVFSM::preprocessTourInformation() has been called.");
    }
#endif

    if (weekday > 6)
        throw std::runtime_error("Error when adding a new vehicle tour for carID " + std::to_string(carID) + ": A weekday > 6 is not possible!");
    
    if (ts_duration == 0)
        ts_duration = 1;

    unsigned long this_new_tour_id = list_of_all_tours.size();
    if (list_of_all_tours.size() > 0) {
        WeeklyVehicleTour* last_know_tour = list_of_all_tours.back();
        // Check, if tours are added in the wrong ordering
        if (last_know_tour->day_of_week > weekday) {
            throw std::runtime_error("Error when adding a new vehicle tour for carID " + std::to_string(carID) + ": Weekday of new tour is bevore the latest added tour!");
        } else if (last_know_tour->day_of_week == weekday && last_know_tour->departure_ts_of_day > departure_ts_of_day) {
            throw std::runtime_error("Error when adding a new vehicle tour for carID " + std::to_string(carID) + ": Time step of departure of new tour is bevore the latest added tour!");
        }
        // Check, if tours are overlapping
        // A) For the previous tour
        float atime_of_week_prev = 24 * Global::get_time_step_size_in_h() * last_know_tour->day_of_week + (float) (last_know_tour->departure_ts_of_day) + (float) (last_know_tour->ts_duration); // arrival time of week of the previous trip
        float dtime_of_week_new  = 24 * Global::get_time_step_size_in_h() * weekday + (float) (departure_ts_of_day); // departure time of week of the trip to add
        if (atime_of_week_prev > dtime_of_week_new) {
            std::cerr << "Warning in carID = " << carID << ": A tour is overlapping with its previous tour (weekday=" << last_know_tour->day_of_week << ", dep. ts=" << last_know_tour->departure_ts_of_day << ", ts. dur=" << last_know_tour->ts_duration << "). ";
            std::cerr << "Ignoring the second tour (weekday=" << weekday << ", dep. ts=" << departure_ts_of_day << ").\n";
            return;
        }
        // B) For the first tour in the next week
        WeeklyVehicleTour* first_known_tour = list_of_all_tours.front();
        unsigned long ts_of_a_week = (unsigned long) std::floor(((double) (7 * 24) / (double) Global::get_time_step_size_in_h()));
        float atime_of_week_new  = 24 * Global::get_time_step_size_in_h() * weekday + (float) (departure_ts_of_day) + (float) (ts_duration); // arrival time of week of the trip to add
        if ((unsigned long) atime_of_week_new > ts_of_a_week) {
            float atime_of_week_first = Global::get_time_step_size_in_h() * first_known_tour->day_of_week + (float) (first_known_tour->departure_ts_of_day);
            if ((unsigned long) (atime_of_week_new) % ts_of_a_week > (unsigned long) atime_of_week_first) {
                std::cerr << "Warning in carID = " << carID << ": A tour is overlapping with the first known tour (weekday=" << first_known_tour->day_of_week << ", dep. ts=" << first_known_tour->departure_ts_of_day << "). ";
                std::cerr << "Ignoring the second tour (weekday=" << weekday << ", dep. ts=" << departure_ts_of_day << ", ts. dur=" << ts_duration << ").\n";
                return;
            }
        }
        // Set next tour ID for the last know tour
        last_know_tour->next_tour_id = this_new_tour_id;
    }
    // append new tour
    WeeklyVehicleTour& new_tour = list_of_tours_pd[weekday]->emplace_back(0, weekday, departure_ts_of_day, ts_duration, tour_length_km, with_work);
    list_of_all_tours.push_back( &new_tour );
}

void EVFSM::preprocessTourInformation() {
#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    if (!state_s1) {
        throw std::runtime_error("Method EVFSM::preprocessTourInformation() must be called once only after the last call of EVFSM::add_weekly_tour().");
    }
    state_s1 = false;
    state_s2 = true;
#endif
    resetInternalState();
    // empty precomputed vectors
    prec_vec_of_states.clear();
    prec_vec_of_minE.clear();
    prec_vec_of_maxE.clear();
    prec_vec_of_driving_distance_km.clear();
    prec_vec_of_maxP_kW.clear();
    prec_vec_of_curr_BS_E_cons_kWh.clear();
    prec_vec_of_states.resize(Global::get_n_timesteps(), EVState::ConnectedAtHome);
    prec_vec_of_minE.resize(Global::get_n_timesteps(), 0.0); // initialize vectors with 0 by default
    prec_vec_of_maxE.resize(Global::get_n_timesteps(), 0.0);
    prec_vec_of_driving_distance_km.resize(Global::get_n_timesteps(), 0.0);
    prec_vec_of_maxP_kW.resize(Global::get_n_timesteps(), 0.0);
    prec_vec_of_curr_BS_E_cons_kWh.resize(Global::get_n_timesteps(), 0.0);
    // list of all tours / complete tour plan
    std::vector<struct SingleVehicleTour> complete_tour_plan;
    // loop over all simulation time steps twice
    // -- first, to generate the tour plan for the complete simulation time span
    unsigned long current_tour_ts_start =    0; // Time step of current tour start (only valid if a tour is ongoing)
    unsigned long ts_since_departure    =    0; // Number of time steps passed until the departure of the current tour (only valid if a tour is ongoing)
    WeeklyVehicleTour* current_wTour    = NULL; ///< Reference to the current weekly tour
    for (unsigned long ts = Global::get_first_timestep(); ts <= Global::get_last_timestep(); ts++) {
        // get current time
        struct tm* current_tm_l = global::time_localtime_l->at(ts - 1);
        int dayOfWeek_l = current_tm_l->tm_wday; // get day of week in the format 0->Sunday, 6->Saturday
        int hourOfDay_l = current_tm_l->tm_hour;
        // Convert format: 0->Monday, 6->Sunday
        dayOfWeek_l = dayOfWeek_l - 1;
        if (dayOfWeek_l < 0)
            dayOfWeek_l += 7;
        // if there is a currently active tour: check, if it reached home
        if (current_wTour != NULL) {
            ts_since_departure += 1;
            // check, if this tour is finished?
            if (ts_since_departure >= current_wTour->ts_duration) {
                double this_tour_econs_kWh = current_wTour->tour_length_km * econs_kWh_per_km;
                complete_tour_plan.emplace_back(current_tour_ts_start, ts, this_tour_econs_kWh, current_wTour);
                // remove current tour
                current_wTour = NULL;
            }
        }
        // does a new tour start?
        if (current_wTour == NULL) {
            // loop over all tours on this day, is there a tour starting right now?
            for ( WeeklyVehicleTour &vt : *(list_of_tours_pd)[dayOfWeek_l] ) {
                if (vt.departure_ts_of_day == (unsigned int) hourOfDay_l) { // mind the shift: Car tours are left-aligned, hours are right aligned
                    // TODO: Handle other time resolution than an hourly resolution
                    current_wTour = &vt;
                    ts_since_departure = 0;
                    current_tour_ts_start = ts;
                    current_state = EVState::Driving;
                    break; // inner loop
                }
            }
        }
    }
    // -- second, to compute upper and lower cumluative charging electricity consumption
    current_state = EVState::ConnectedAtHome; // the EV is at home (connected) at the beginning of the simulation run 
    double last_minE_value = 0.0; // actual minE value
    double last_maxE_value = 0.0; // actual maxE value ( <= cumsum_E_driving_kWh )
    double cumsum_E_driving_kWh = 0.0; // cumulative sum of electricity used for driving
    double cumsum_driving_distance_km = 0.0;
    const double E_CHARGABLE_PER_TS_kWh = Global::get_time_step_size_in_h() * Global::get_ev_max_charging_power_kW();
    struct SingleVehicleTour* current_sTour = NULL; // dragging pointer to current single tour
    //struct SingleVehicleTour* next_sTour    = NULL; // dragging pointer to next single tour
    auto next_sTour = complete_tour_plan.begin();
    // loop over all tours that have finished until the beginning of the simulation
    for (unsigned long ts = 1; ts < Global::get_first_timestep(); ts++) {
        if (current_sTour != NULL && current_sTour->ts_arrival >= ts) {
            current_sTour = NULL;
            current_state = EVState::ConnectedAtHome;
        }
        if (next_sTour != complete_tour_plan.end() && next_sTour->ts_start >= ts) {
            current_sTour = &(*next_sTour);
            current_state = EVState::Driving;
            next_sTour++;
        }
        prec_vec_of_states.at(ts-1) = current_state;
    }
    // main loop for simulation time span
    bool ev_fully_charged_at_next_dep = false; // True, if the EV must be fully charged at the next departure (or at least charged in every time step)!
    unsigned long ts_since_last_connection = 0;
    // double current_min_cumsum_SOE_kWh = battery->get_maxE_kWh(); // Start with a full battery
    for (unsigned long ts = Global::get_first_timestep(); ts <= Global::get_last_timestep(); ts++) {
        unsigned long tsID = ts - 1;
        ts_since_last_connection++;
        // is the currently active tour finished?
        if (current_sTour != NULL && current_sTour->ts_arrival <= ts) {
            // arrival of current tour
            cumsum_E_driving_kWh += current_sTour->energy_consumption_kWh;
            cumsum_driving_distance_km += current_sTour->weekly_tour->tour_length_km;
            // Computation of EVState: Connected or not connected at home?
            if (Global::get_ev_plugin_probability() >= 1.0) {
                current_state = EVState::ConnectedAtHome;
                ts_since_last_connection = 0;
            } else {
                // Pluggin-in required, if: SOC <= 0.35, or if SOC is too low for next tour, or if sampling says so
                if (
                    battery->get_SOC() <= 0.35 ||
                    (*distribution)(*random_generator) <= Global::get_ev_plugin_probability()
                ) {
                    current_state = EVState::ConnectedAtHome;
                    ts_since_last_connection = 0;
                    // Attention: in this case (i.e., plugin probability < 1.0 and if connected), EV MUST always be fully charged (or at least as much as possible in the time at home)
                    ev_fully_charged_at_next_dep = true;
                } else {
                    current_state = EVState::DisconnectedAtHome;
                }
            }
            // remove current tour
            current_sTour = NULL;
            energy_demand_per_tour_ts = 0.0;
        }
        // is there no current tour, and is there still a next tour, and if yes, does this tour start right now?
        if (current_sTour == NULL && next_sTour != complete_tour_plan.end() && next_sTour->ts_start <= ts) {
            current_sTour = &(*next_sTour);
            current_state = EVState::Driving;
            // compute (mean) energy demand per tour time step
            energy_demand_per_tour_ts = (current_sTour->weekly_tour->tour_length_km * econs_kWh_per_km) / (double) (current_sTour->weekly_tour->ts_duration);
            // Compute new last_minE_value
            const double PREV_last_minE_value = last_minE_value;
            if (ev_fully_charged_at_next_dep) {
                // BS must be full (or set to max. chargeable amount) on departure in this case
                last_minE_value = last_maxE_value;
            } else {
                last_minE_value = last_maxE_value - battery->get_maxE_kWh();
                if (last_minE_value < 0)
                    last_minE_value = 0.0;
            }
            // check for infeasability (i.e., car is consuming more energy than can be charged)
            if (last_maxE_value < last_minE_value) {
                std::cerr << "Error in processing data for CarID " << carID << ". last_maxE_value < last_minE_value at time step " << ts << ". Optimization model will get infeasible!" << std::endl;
            }
            // set new EV battery state of charge
            double must_be_charged_kWh = last_minE_value - PREV_last_minE_value;
            battery->set_SOE_without_computations(must_be_charged_kWh);
            // set next tour iterator to next tour and unset boolean variables
            next_sTour++;
            ev_fully_charged_at_next_dep = false;
        }
        // Update the current EV battery state if driving
        if (current_state == EVState::Driving) {
            // remove consumed energy per step from the battery
            battery->set_chargeRequest(-energy_demand_per_tour_ts);
            battery->calculateActions();
            prec_vec_of_curr_BS_E_cons_kWh[tsID] = energy_demand_per_tour_ts;
        }
        // Update current maxE values (if connected at home, else use old values)
        // and set current available max charging power
        if (current_state == EVState::ConnectedAtHome) {
            last_maxE_value += E_CHARGABLE_PER_TS_kWh;
            if (last_maxE_value > cumsum_E_driving_kWh) {
                last_maxE_value = cumsum_E_driving_kWh;
            }
            prec_vec_of_maxP_kW[tsID] = Global::get_ev_max_charging_power_kW();
        } // else: prec_vec_of_maxP_kW[tsID] = 0 by default
        // Write current state values to the pre-computed dict
        prec_vec_of_states[tsID] = current_state;
        prec_vec_of_minE[tsID] = last_minE_value;
        prec_vec_of_maxE[tsID] = last_maxE_value;
        prec_vec_of_driving_distance_km[tsID] = cumsum_driving_distance_km;
    }
    // Update prec_vec_of_minE / maxE until the end of (without sampling what will happen)
    for (unsigned long ts = Global::get_last_timestep() + 1; ts < Global::get_n_timesteps(); ts++) {
        prec_vec_of_minE.at(ts-1) = last_minE_value;
        prec_vec_of_maxE.at(ts-1) = last_maxE_value;
        prec_vec_of_states.at(ts-1) = current_state;
        prec_vec_of_driving_distance_km.at(ts-1) = cumsum_driving_distance_km;
    }
    // reset internal state again for the main simulation run
    resetInternalState();
}

void EVFSM::resetInternalState() {
    battery->resetInternalState();
    current_state      = EVState::ConnectedAtHome;
    //current_state_icah = EVStateIfConnAtHome::ChargingPossible;
    energy_demand_per_tour_ts = 0.0;
    current_P_kW              = 0.0;
    sum_of_driving_distance_km    = 0.0;
    sum_of_E_used_for_driving_kWh = 0.0;
    sum_of_E_charged_home_kWh     = 0.0;
    sum_of_E_discharged_home_kWh  = 0.0;
    sum_of_ts_EV_is_connected     = 0;
}

ComponentStateVariant EVFSM::saveInternalState() const {
    EVFSMComponentState state;
    state.carID = carID;
    state.current_state = current_state;
    state.current_ts = current_ts;
    state.energy_demand_per_tour_ts = energy_demand_per_tour_ts;
    state.current_P_kW = current_P_kW;
    state.sum_of_driving_distance_km = sum_of_driving_distance_km;
    state.sum_of_E_used_for_driving_kWh = sum_of_E_used_for_driving_kWh;
    state.sum_of_E_charged_home_kWh = sum_of_E_charged_home_kWh;
    state.sum_of_E_discharged_home_kWh = sum_of_E_discharged_home_kWh;
    state.sum_of_ts_EV_is_connected = sum_of_ts_EV_is_connected;
    
    // Save the nested battery state
    ComponentStateVariant battery_state_variant = battery->saveInternalState();
    state.battery_state = std::get<BSComponentState>(battery_state_variant);
    
    return state;
}

void EVFSM::restoreInternalState(const ComponentStateVariant& state) {
    try {
        const EVFSMComponentState& ev_state = std::get<EVFSMComponentState>(state);
        
        // Verify that the carID matches
        if (ev_state.carID != carID) {
            throw std::runtime_error("EVFSM::restoreInternalState(): CarID mismatch - expected " + std::to_string(carID) + ", got " + std::to_string(ev_state.carID));
        }
        
        current_state = ev_state.current_state;
        current_ts = ev_state.current_ts;
        energy_demand_per_tour_ts = ev_state.energy_demand_per_tour_ts;
        current_P_kW = ev_state.current_P_kW;
        sum_of_driving_distance_km = ev_state.sum_of_driving_distance_km;
        sum_of_E_used_for_driving_kWh = ev_state.sum_of_E_used_for_driving_kWh;
        sum_of_E_charged_home_kWh = ev_state.sum_of_E_charged_home_kWh;
        sum_of_E_discharged_home_kWh = ev_state.sum_of_E_discharged_home_kWh;
        sum_of_ts_EV_is_connected = ev_state.sum_of_ts_EV_is_connected;
        
        // Restore the nested battery state
        ComponentStateVariant battery_state_variant = ev_state.battery_state;
        battery->restoreInternalState(battery_state_variant);
        
    } catch (const std::bad_variant_access& e) {
        throw std::runtime_error("EVFSM::restoreInternalState(): Invalid state type provided");
    }
}

void EVFSM::setCarStateForTimeStep(unsigned long ts) {
#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    if (!state_s2) {
        throw std::runtime_error("Method EVFSM::setCarStateForTimeStep() cannot be called at the moment!");
    }
    state_s2 = false;
    state_s3 = true;
#endif
    current_ts = ts;
    unsigned long tsID = ts - 1;
    // Read values from precomputed vectors
    current_state = prec_vec_of_states[tsID];
    sum_of_E_used_for_driving_kWh = prec_vec_of_maxE[tsID];
    sum_of_driving_distance_km = prec_vec_of_driving_distance_km[tsID];
    current_P_kW = 0.0;
    // set inherited variables for semi-flexible component
    double last_known_maxE_val = 0.0;
    double last_known_minE_val = 0.0;
    double last_known_maxP_val = 0.0;
    for (size_t tOffset = 0; tOffset < Global::get_control_horizon_in_ts(); tOffset++) {
        size_t tsIDaOffset = tsID + tOffset;
        if (tsIDaOffset < Global::get_n_timesteps()) {
            last_known_maxE_val = prec_vec_of_maxE[tsIDaOffset] - sum_of_E_charged_home_kWh;
            last_known_minE_val = prec_vec_of_minE[tsIDaOffset] - sum_of_E_charged_home_kWh;
            if (last_known_maxE_val < 0) // this might happen due to precission/rounding errors
                last_known_maxE_val = 0.0;
            if (last_known_minE_val < 0)
                last_known_minE_val = 0.0;
            last_known_maxP_val = prec_vec_of_maxP_kW[tsIDaOffset];
        }
        future_maxE_storage[tOffset] = last_known_maxE_val;
        future_minE_storage[tOffset] = last_known_minE_val;
        future_maxP_storage[tOffset] = last_known_maxP_val;
    }
    // in the end (of the complete simulation span), min and max profiles must come to the same point
    unsigned long last_tOffset = Global::get_control_horizon_in_ts() - 1;
    if (last_tOffset + ts >= Global::get_n_timesteps()) {
        future_minE_storage[last_tOffset] = future_maxE_storage[last_tOffset];
    }
    // Remove energy from the battery if car is driving
    if (current_state == EVState::Driving) {
        battery->set_chargeRequest(-prec_vec_of_curr_BS_E_cons_kWh[tsID]);
        battery->calculateActions();
    } else if (current_state == EVState::ConnectedAtHome) {
        // Increment counters
        sum_of_ts_EV_is_connected += 1;
    }
}

void EVFSM::setDemandToProfileData(unsigned long ts) {
#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    if (!state_s3) {
        throw std::runtime_error("Method EVFSM::setDemandToProfileData() cannot be called at the moment!");
    }
    state_s2 = true;
    state_s3 = false;
#endif
    // immediate charging
    if (current_state == EVState::ConnectedAtHome && battery->get_SOC() < 1.0) {
        battery->set_chargeRequest( Global::get_ev_max_charging_power_kW() );
        battery->calculateActions();
        // update the amount of charged / discharged electricity
        current_P_kW = battery->get_currentLoad_kW();
        if ( current_P_kW > 0 ) {
            sum_of_E_charged_home_kWh    += current_P_kW * Global::get_time_step_size_in_h();
        } else if ( current_P_kW < 0 ) {
            sum_of_E_discharged_home_kWh -= current_P_kW * Global::get_time_step_size_in_h();
        }
    } else {
        current_P_kW = 0.0;
    }

    if (Global::get_create_ev_detailed_output()) {
        output::outputEVStateDetails(ts, carID, current_state, (float) current_P_kW, (float) sum_of_E_charged_home_kWh, (float) prec_vec_of_minE[current_ts-1], (float) prec_vec_of_maxE[current_ts-1], (float) battery->get_currentCharge_kWh());
    }
}

// define a small value added to the min/max cumsum energy consumption to remove misplaced warnings due to rounding errors
#define epsilon_ev 0.001

bool EVFSM::setDemandToGivenValue(double new_demand_kW) {
#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    if (!state_s3) {
        throw std::runtime_error("Method EVFSM::setDemandToGivenValue() cannot be called at the moment!");
    }
    state_s2 = true;
    state_s3 = false;
#endif
    bool no_error = true;
    if (new_demand_kW < 0) {
        if (new_demand_kW < -epsilon_ev) {
            std::cerr << "Warning: Bidirection EV charging (" << std::fixed << std::setprecision(3) << new_demand_kW << " kW) currently not supported for EVFSM::setDemandToGivenValue()." << std::endl;
            return false;
        }
        return true; // ignore small deviations below epsilon_ev
    }
    double e = new_demand_kW * Global::get_time_step_size_in_h();
    double new_total_e = sum_of_E_charged_home_kWh + e;
    // check, if the new demand is within the min/max bands
    // Output warnings only if absolute bound violation > epsilon_ev
    if (new_total_e > prec_vec_of_maxE[current_ts-1]) {
        if (new_total_e > prec_vec_of_maxE[current_ts-1] + epsilon_ev) {
            no_error = false;
            std::cerr << "Warning in EVFSM: Set demand (" << std::fixed << std::setprecision(3) << new_demand_kW << " kW) is violating the upper bound for carID " << carID << " at time step " << current_ts << "!";
        }
        e = prec_vec_of_maxE[current_ts-1] - sum_of_E_charged_home_kWh;
        new_demand_kW = e / Global::get_time_step_size_in_h();
        //std::cerr << "Setting the new demand to " << std::fixed << std::setprecision(3) << new_demand_kW << " kW.\n";
    } else if (new_total_e < prec_vec_of_minE[current_ts-1]) {
        if (new_total_e < prec_vec_of_minE[current_ts-1] - epsilon_ev) {
            no_error = false;
            std::cerr << "Warning in EVFSM: Set demand (" << std::fixed << std::setprecision(3) << new_demand_kW << " kW) is violating the lower bound for carID " << carID << " at time step " << current_ts << "!";
        }
        e = prec_vec_of_minE[current_ts-1] - sum_of_E_charged_home_kWh;
        new_demand_kW = e / Global::get_time_step_size_in_h();
        //std::cerr << "Setting the new demand to " << std::fixed << std::setprecision(3) << new_demand_kW << " kW.\n";
    }
    // update current demand
    battery->set_chargeRequest( new_demand_kW );
    battery->calculateActions();
    current_P_kW = new_demand_kW;
    // update cumulative variables
    sum_of_E_charged_home_kWh += e;

    if (Global::get_create_ev_detailed_output()) {
        output::outputEVStateDetails(current_ts, carID, current_state, (float) current_P_kW, (float) sum_of_E_charged_home_kWh, (float) prec_vec_of_minE[current_ts-1], (float) prec_vec_of_maxE[current_ts-1], (float) battery->get_currentCharge_kWh());
    }
    return no_error;
}

void EVFSM::AddWeeklyTour(
    unsigned long carID,              short weekday,
    unsigned long departure_ts_of_day,unsigned long ts_duration,
    double tour_length_km,            bool with_work)
{
    EVFSM::list_of_cars.at(carID)->add_weekly_tour(weekday, departure_ts_of_day, ts_duration, tour_length_km, with_work);
}

void EVFSM::VacuumStaticVariables() {
    EVFSM::list_of_cars.clear();
}

void EVFSM::SetSeed(unsigned int seed) {
    random_generator->seed(seed);
}

