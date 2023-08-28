
#include "components.h"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include <ranges>
#include <string>
#include <vector>

#include "vehicles.h"

// ----------------------------- //
//      Implementation of        //
//        RoofSectionPV          //
// ----------------------------- //

std::map<std::string, size_t> RoofSectionPV::next_pv_idx;
std::map<std::string, std::default_random_engine>            RoofSectionPV::random_generators;
std::map<std::string, std::uniform_int_distribution<size_t>> RoofSectionPV::distributions;

RoofSectionPV::RoofSectionPV(float this_section_kWp, std::string& orientation)
    : this_section_kWp(this_section_kWp), orientation(orientation)
{
    // 1)
    // take the correct profile
    profile_index = 0;
    if (Global::get_exp_pv_static_mode() &&
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
        if (!random_generators.contains(orientation)) {
            random_generators[orientation] = std::default_random_engine{};
            distributions[orientation] = std::uniform_int_distribution<size_t>(0, global::pv_profiles_information[orientation]-1);
        }
        profile_index = distributions[orientation](random_generators[orientation]);
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

float RoofSectionPV::get_currentFeedin_kW(unsigned long ts) {
    unsigned long tsID = ts - 1;
    return this_section_kWp * profile_data[tsID];
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
    total_generation_kWh = 0.0;

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
        roof_sections.emplace_back(section_kWp, section_orientation);
    }
  } else {
    // Case 2: Only use one component facing the given orientation
    string ori = Global::get_exp_pv_static_profile_orientation();
    roof_sections.emplace_back(kWp, ori);
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

    // general checks
    if (kWp_per_m2 <= 0.0) {
        throw runtime_error("Error in ComponentPV constructor: kWp_per_m2 <= 0!");
    }
    if (max_kWp > 0 && min_kWp > max_kWp) {
        throw runtime_error("Error in ComponentPV constructor: min_kWp > max_kWp");
    }

    currentGeneration_kW = 0;
    total_generation_kWh = 0.0;
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
        if (section_kWp < min_kWp) {
            // ignore this section
            continue;
        }
        // if max_kWp == 0 -> ignore this value
        if (max_kWp > 0 && section_kWp > max_kWp) {
            section_kWp = max_kWp;
        }
        total_kWp += section_kWp;
        // 2)
        // add section to list
        vec_of_sections.emplace_back(section_kWp, section_orientation);
    }
    // 2)
    // only in the case of Global::exp_pv_max_kWp_per_unit is set (i.e. >= 0.0)
    if (Global::get_exp_pv_max_kWp_per_unit() >= 0.0) {
        // is total_kWp > Global::get_exp_pv_max_kWp_per_unit()
        // if yes, lower installed power with the same percantage over all sections
        if ( total_kWp > Global::get_exp_pv_max_kWp_per_unit() ) {
            float reduction_factor = Global::get_exp_pv_max_kWp_per_unit() / total_kWp;
            total_kWp = Global::get_exp_pv_max_kWp_per_unit();
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
        roof_sections.emplace_back(section_kWp, section_orientation);
    }
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

void ComponentPV::calculateCurrentFeedin(unsigned long ts) {
    currentGeneration_kW = 0.0;
    for (RoofSectionPV& section : roof_sections)
        currentGeneration_kW += section.get_currentFeedin_kW(ts);
    // compute total generation
    total_generation_kWh += Global::get_time_step_size_in_h() * currentGeneration_kW;
}

void ComponentPV::resetInternalState() {
    total_generation_kWh = 0.0;
    currentGeneration_kW = 0.0;
}





// ----------------------------- //
//      Implementation of        //
//         ComponentBS           //
// ----------------------------- //

ComponentBS::ComponentBS(
    float maxE_kWh,
    float maxP_kW,
    float E_over_P_ratio,
    float discharge_rate_per_step,
    float efficiency_in,
    float efficiency_out,
    float initial_SoC
) : maxE_kWh(maxE_kWh),
    discharge_rate_per_step(discharge_rate_per_step), efficiency_in(efficiency_in),
    efficiency_out(efficiency_out), initial_SoC(initial_SoC)
{
    SOC               = 0;
    currentE_kWh      = 0;
    currentP_kW       = 0;
    charge_request_kW = 0;
    total_E_withdrawn_kWh = 0.0;
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
    float maxE_kWh,
    float discharge_rate_per_step,
    float efficiency_in,
    float initial_SoC
) : maxE_kWh(maxE_kWh),
    discharge_rate_per_step(discharge_rate_per_step),
    efficiency_in(efficiency_in),
    efficiency_out(1.0),
    initial_SoC(initial_SoC)
{
    SOC               = 0;
    currentE_kWh      = 0;
    currentP_kW       = maxE_kWh; // This is required so that the battery can be emptied in one single step at the arrival
    charge_request_kW = 0;
    total_E_withdrawn_kWh = 0.0;
    n_ts_SOC_empty    = 0;
    n_ts_SOC_full     = 0;

    if (initial_SoC > 0) {
        SOC = initial_SoC;
        currentE_kWh = maxE_kWh * initial_SoC;
    }
}

void ComponentBS::set_maxE_kWh(float value) {
    maxE_kWh = value;
    if (Global::get_battery_power_computation_mode() == global::BatteryPowerComputationMode::UseEOverPRatio)
        this->maxP_kW = maxE_kWh / E_over_P_ratio;
}

void ComponentBS::set_maxP_kW(float value) {
    if (Global::get_battery_power_computation_mode() == global::BatteryPowerComputationMode::AsDefinedByConfigVar) {
        maxP_kW  = value;
        E_over_P_ratio = maxE_kWh / maxP_kW;
    }
}

void ComponentBS::set_maxP_by_EPRatio(float EP_ratio) {
    if (Global::get_battery_power_computation_mode() == global::BatteryPowerComputationMode::UseEOverPRatio) {
        this->E_over_P_ratio = EP_ratio;
        this->maxP_kW = this->maxE_kWh / EP_ratio;
    }
}

void ComponentBS::calculateActions() {
    float timestep_size_in_h = Global::get_time_step_size_in_h();
    float new_charge_kWh;

    currentP_kW = 0;

    // Calculate Self-discharge
    currentE_kWh -= discharge_rate_per_step * currentE_kWh;

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
        float energy_taken_kWh = new_charge_kWh - currentE_kWh;
        currentP_kW  = energy_taken_kWh / timestep_size_in_h * efficiency_out;
        currentE_kWh = new_charge_kWh;
        // add withrawn energy to summation variable (mind energy_taken_kWh < 0)
        total_E_withdrawn_kWh -= energy_taken_kWh;
    }

    // calculate new SOC value
    SOC = currentE_kWh / maxE_kWh;

    // calculate n_ts_SOC_empty or n_ts_SOC_full
    if (currentE_kWh <= 0)
        n_ts_SOC_empty++;
    else if (currentE_kWh >= maxE_kWh)
        n_ts_SOC_full++;
}

void ComponentBS::resetInternalState() {
    //
    // This method resets the internal SOC to the initial level
    //
    SOC = initial_SoC;
    currentE_kWh = maxE_kWh * initial_SoC;
    total_E_withdrawn_kWh = 0.0;
    n_ts_SOC_empty = 0;
    n_ts_SOC_full  = 0;
}





// ----------------------------- //
//      Implementation of        //
//         ComponentHP           //
// ----------------------------- //

size_t ComponentHP::next_hp_idx = 0;
bool   ComponentHP::random_generator_init = false;
std::default_random_engine*            ComponentHP::random_generator = NULL;
std::uniform_int_distribution<size_t>* ComponentHP::distribution     = NULL;

ComponentHP::ComponentHP(float yearly_econs_kWh)
    : yearly_electricity_consumption_kWh(yearly_econs_kWh),
      scaling_factor(yearly_econs_kWh/1000.0f)
{
    // select heat pump profile static or random
    size_t this_hp_profile_idx;
    if (Global::get_exp_profile_mode() == global::ExpansionProfileAllocationMode::AsInData) {
        this_hp_profile_idx = next_hp_idx;
        // increment next index by one
        next_hp_idx++;
        if (next_hp_idx >= Global::get_n_heatpump_profiles())
            next_hp_idx = 0;
    } else {
        if (!random_generator_init)
            ComponentHP::InitializeRandomGenerator();
        // randomly select new index
        this_hp_profile_idx = (*distribution)(*random_generator);
    }
    //
    // reference the profile
    profile_data = global::hp_profiles[this_hp_profile_idx];
    // further initialization
    currentDemand_kW = 0;
    total_demand_kWh = 0.0;
}

void ComponentHP::calculateCurrentFeedin(unsigned long ts) {
    unsigned long tsID = ts - 1;
    currentDemand_kW = profile_data[tsID] * scaling_factor;
}

void ComponentHP::resetInternalState() {
    currentDemand_kW = 0.0;
    total_demand_kWh = 0.0;
}

void ComponentHP::InitializeRandomGenerator() {
    if (random_generator_init) {
        throw runtime_error("Error: Static variables of ComponentHP are already initialized.");
    } else {
        random_generator_init = true;
        random_generator = new std::default_random_engine();
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

ComponentCS::ComponentCS() : max_charging_power(11) { // TODO: make max charging power configurable, or at least dependent of the number of flats in a building
    enabled = false;
    current_demand_kW = 0.0;
    total_demand_kWh  = 0.0;
    charging_power_required_kW = 0.0;
    charging_power_possible_kW = 0.0;
}

ComponentCS::~ComponentCS() {
    for (EVFSM* ev : listOfEVs)
        delete ev;
}

float ComponentCS::get_max_P_kW() const {
    if (enabled)
        return max_charging_power;
    return 0.0;
}

double ComponentCS::get_max_curr_charging_power_kW() const {
    if (charging_power_required_kW + charging_power_possible_kW >= max_charging_power) {
        return max_charging_power;
    } else {
        return charging_power_required_kW + charging_power_possible_kW;
    }
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

void ComponentCS::enable_station() {
    enabled = true;
}

void ComponentCS::disable_station() {
    enabled = false;
}

void ComponentCS::resetInternalState() {
    current_demand_kW = 0.0;
    total_demand_kWh  = 0.0;
    charging_power_required_kW = 0.0;
    charging_power_possible_kW = 0.0;
    charging_order_req.clear();
    charging_order_pos.clear();
    //
    for (EVFSM* ev : listOfEVs) {
        ev->resetInternalState();
    }
}

void ComponentCS::add_ev(unsigned long carID) {
    listOfEVs.push_back(new EVFSM(carID, this));
}

void ComponentCS::setCarStatesForTimeStep(unsigned long ts, int dayOfWeek_l, int hourOfDay_l) {
    // 1. set new car states
    for (EVFSM* ev : listOfEVs) {
        ev->setCarStateForTimeStep(ts, dayOfWeek_l, hourOfDay_l);
    }
    // 2. iterate over all cars currently at home -> how much energy do they require to charge?
    charging_power_required_kW = 0.0;
    charging_power_possible_kW = 0.0;
    charging_order_req.clear();
    charging_order_pos.clear();
    // 2a. The cars that require charging
    for (EVFSM* ev : listOfEVs
                   | std::views::filter(
                        [](EVFSM* ev) {
                            return ev->get_current_state() == EVState::ConnectedAtHome &&
                                   ev->get_current_state_icah() == EVStateIfConnAtHome::ChargingRequired;
                        }))
    {
        charging_power_required_kW += ev->get_max_curr_charging_power_kW();
        charging_order_req.push_back(ev);
    }
    // 2b. The cars that can charge
    for (EVFSM* ev : listOfEVs
                   | std::views::filter(
                        [](EVFSM* ev) {
                            return ev->get_current_state() == EVState::ConnectedAtHome &&
                                   ( ev->get_current_state_icah() == EVStateIfConnAtHome::ChargingPossible ||
                                     ev->get_current_state_icah() == EVStateIfConnAtHome::BothPossible);
                        }))
    {
        charging_power_possible_kW += ev->get_max_curr_charging_power_kW();
        charging_order_pos.push_back(ev);
    }
}

void ComponentCS::set_charging_value(double power_kW) {
    current_demand_kW = 0.0;
    double remaining_power_kW = max_charging_power - current_demand_kW; // The remaining power
    for (EVFSM* ev : charging_order_req) {
        // max. 11 kW per car
        // Set charging request
        if (remaining_power_kW >= 11.0) {
            ev->set_current_charging_power(11.0);
        } else {
            ev->set_current_charging_power(remaining_power_kW);
        }
        // get actual charging power
        remaining_power_kW -= ev->get_current_charging_power();
    }
    if (remaining_power_kW > 0) {
        for (EVFSM* ev : charging_order_pos) {
            // max. 11 kW per car
            // Set charging request
            if (remaining_power_kW >= 11.0) {
                ev->set_current_charging_power(11.0);
            } else {
                ev->set_current_charging_power(remaining_power_kW);
            }
            // get actual charging power
            remaining_power_kW -= ev->get_current_charging_power();
        }
    }
    // compute cumulative sum
    total_demand_kWh += current_demand_kW * Global::get_time_step_size_in_h();
}
