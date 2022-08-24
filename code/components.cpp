
#include "components.h"

#include <iostream>
#include <numeric>
#include <string>

// ----------------------------- //
//      Implementation of        //
//        RoofSectionPV          //
// ----------------------------- //

RoofSectionPV::RoofSectionPV(float this_section_kWp, std::string orientation, size_t profile_index)
    : this_section_kWp(this_section_kWp)
{
    // select the correct profile
    if (global::pv_profiles_per_ori[orientation].size() <= 0) {
        std::cerr << "Error: There is no feed-in profile given for the orientation " << orientation << std::endl;
        throw runtime_error("There is no feed-in profile given for a selected orientation!");
    }
    profile_data = global::pv_profiles_per_ori.at(orientation)[profile_index];
}

float RoofSectionPV::get_currentFeedin_kW(int ts) {
    int tsID = ts - 1;
    return this_section_kWp * profile_data[tsID];
}




// ----------------------------- //
//      Implementation of        //
//         ComponentPV           //
// ----------------------------- //

std::map<std::string, size_t> ComponentPV::next_pv_idx;

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
    // attach roof sections as defined in data
    auto& roof_section_vec = global::roof_section_orientations[locationID];
    //unsigned long number_of_sections = roof_section_vec.size();
    float complete_roof_area = accumulate(
                begin(roof_section_vec),
                end(roof_section_vec),
                0,
                [](const size_t prev, const auto& elem){ return prev + elem.first; });;
    // iterate over all roof sections
    // and get all orientations and areas per roof section
    for (auto& section_tuple : roof_section_vec) {
        float  section_roof_area   = section_tuple.first;
        std::string section_orientation = section_tuple.second;
        // 1)
        // take the correct profile
        size_t pv_profile_idx = 0;
        // TODO: select randomly as an alternative
        if (next_pv_idx.contains(section_orientation)) {
            size_t& next_pv_idx_so = next_pv_idx[section_orientation];
            pv_profile_idx = next_pv_idx_so;
            next_pv_idx_so++; // increment next index by one until end is reached
            if (next_pv_idx_so >= global::pv_profiles_information[section_orientation])
                next_pv_idx_so = 0;
        }
        // 2)
        // compute the share of this section among all existing sections
        float share_of_total_area = section_roof_area / complete_roof_area;
        // 3)
        // create and add section to list
        float section_kWp = share_of_total_area * kWp;
        roof_sections.emplace_back(section_kWp, section_orientation, pv_profile_idx);
    }
}

ComponentPV::ComponentPV(float kWp_per_m2, float min_kWp, float max_kWp, unsigned long locationID)
{
    /*
     *
     * Constructor for ComponentPV
     * in the case of a dynamic kWp computation
     *
     */

    currentGeneration_kW = 0;
    total_kWp = 0;
    // attach roof sections as defined in data
    // // float min_roof_area = min_kWp / kWp_per_m2;
    // // float max_roof_area = max_kWp / kWp_per_m2;
    auto& roof_section_vec = global::roof_section_orientations[locationID];
    // iterate over all roof sections
    // and get all orientations and areas per roof section
    for (auto& section_tuple : roof_section_vec) {
        float       section_roof_area   = section_tuple.first;
        std::string section_orientation = section_tuple.second;
        // 1)
        // compute the kWp of this area and decide whether to use it or not
        float section_kWp = section_roof_area * kWp_per_m2;
        if (section_kWp < min_kWp) {
            // ignore this section
            continue;
        }
        if (section_kWp > max_kWp) {
            section_kWp = max_kWp;
        }
        total_kWp += section_kWp;
        // 2)
        // take the correct profile
        size_t pv_profile_idx = 0;
        // TODO: select randomly as an alternative
        if (next_pv_idx.contains(section_orientation)) {
            size_t& next_pv_idx_so = next_pv_idx[section_orientation];
            pv_profile_idx = next_pv_idx_so;
            next_pv_idx_so++; // increment next index by one until end is reached
            if (next_pv_idx_so >= global::pv_profiles_information[section_orientation])
                next_pv_idx_so = 0;
        }
        // 3)
        // create and add section to list
        roof_sections.emplace_back(section_kWp, section_orientation, pv_profile_idx);
    }
}

void ComponentPV::calculateCurrentFeedin(int ts) {
    currentGeneration_kW = 0.0;
    for (RoofSectionPV& section : roof_sections)
        currentGeneration_kW += section.get_currentFeedin_kW(ts);
}





// ----------------------------- //
//      Implementation of        //
//         ComponentBS           //
// ----------------------------- //

ComponentBS::ComponentBS(float maxE_kWh, float maxP_kW,
        float discharge_rate_per_step, float efficiency, float initial_SoC) : maxE_kWh(maxE_kWh),maxP_kW(maxP_kW),
        discharge_rate_per_step(discharge_rate_per_step), efficiency(efficiency),
        initial_SoC(initial_SoC) {

    SOC               = 0;
    currentE_kWh      = 0;
    currentP_kW       = 0;
    charge_request_kW = 0;

    if (initial_SoC > 0) {
        SOC = initial_SoC;
        currentE_kWh = maxE_kWh * initial_SoC;
    }

}

void ComponentBS::calculateActions() {
    float timestep_size_in_h = Global::get_time_step_size_in_h();
    float new_charge_kWh;

    currentP_kW = 0;

    // Calculate Self-discharge
    currentE_kWh -= discharge_rate_per_step * currentE_kWh;
    // Calculate efficiency
    charge_request_kW = charge_request_kW * efficiency;

    // Charging and discharging
    if (charge_request_kW > 0) {
        // charging requested
        if (charge_request_kW > maxP_kW)
            charge_request_kW = maxP_kW;
        new_charge_kWh = currentE_kWh + timestep_size_in_h*charge_request_kW;
        if (new_charge_kWh > maxE_kWh)
            new_charge_kWh = maxE_kWh;
        currentP_kW  = (new_charge_kWh - currentE_kWh)/timestep_size_in_h;
        currentE_kWh = new_charge_kWh;
    } else if (charge_request_kW < 0) {
        // discharging requested
        if (-charge_request_kW > maxP_kW)
            charge_request_kW = -maxP_kW;
        new_charge_kWh = currentE_kWh + timestep_size_in_h*charge_request_kW;
        if (new_charge_kWh < 0)
            new_charge_kWh = 0;
        currentP_kW  = (new_charge_kWh - currentE_kWh)/timestep_size_in_h;
        currentE_kWh = new_charge_kWh;
    }

    // calculate new SOC value
    SOC = currentE_kWh / maxE_kWh;
}

void ComponentBS::resetInternalState() {
    //
    // This method resets the internal SOC to the initial level
    //
    SOC = initial_SoC;
    currentE_kWh = maxE_kWh * initial_SoC;
}





// ----------------------------- //
//      Implementation of        //
//         ComponentHP           //
// ----------------------------- //

ComponentHP::ComponentHP(size_t profile_index, float yearly_econs_kWh)
    : yearly_electricity_consumption_kWh(yearly_econs_kWh),
      scaling_factor(yearly_econs_kWh/1000.0)
{
    // reference the profile
    profile_data = global::hp_profiles[profile_index];
    // further initialization
    currentDemand_kW = 0;
}

void ComponentHP::calculateCurrentFeedin(int ts) {
    int tsID = ts - 1;
    currentDemand_kW = profile_data[tsID] * scaling_factor;
}
