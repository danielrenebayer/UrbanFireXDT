
#include "components.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <ranges>
#include <string>
#include <vector>

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
            if (Global::is_seed_set()) random_generators[orientation].seed(Global::get_seed());
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

float RoofSectionPV::get_generation_at_ts_kW(unsigned long ts) const {
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
        roof_sections.emplace_back(section_kWp, section_orientation);
    }
  } else {
    // Case 2: Only use one component facing the given orientation
    string ori = Global::get_exp_pv_static_profile_orientation();
    roof_sections.emplace_back(kWp, ori);
  }

    // initialization of the cached vector
    set_horizon_in_ts( Global::get_control_horizon_in_ts() );
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

void ComponentPV::set_horizon_in_ts(unsigned int new_horizon) {
    future_generation_kW.clear();
    future_generation_kW.resize(new_horizon, 0.0);
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
    currentP_kW       = 0;
    charge_request_kW = 0;
    cweek_E_withdrawn_kWh = 0.0;
    total_E_withdrawn_kWh = 0.0;
    n_ts_SOC_empty    = 0;
    n_ts_SOC_full     = 0;

    maxP_kW = maxE_kWh / Global::get_time_step_size_in_h(); // This is required so that the battery can be emptied in one single step at the arrival

    if (initial_SoC > 0) {
        SOC = initial_SoC;
        currentE_kWh = maxE_kWh * initial_SoC;
    }
}

void ComponentBS::set_SOE_without_computations(float new_SOE_kWh) {
    currentE_kWh = new_SOE_kWh;
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
        cweek_E_withdrawn_kWh -= energy_taken_kWh;
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
}

void ComponentBS::resetInternalState() {
    //
    // This method resets the internal SOC to the initial level
    //
    SOC = initial_SoC;
    currentE_kWh = maxE_kWh * initial_SoC;
    cweek_E_withdrawn_kWh = 0.0;
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
    // initialize the unshiftable load storage
    future_unshiftable_storage.clear();
    future_unshiftable_storage.resize( Global::get_control_horizon_in_ts(), 0.0);

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
    profile_data_shiftable = global::hp_profiles_shiftable[this_hp_profile_idx];
    profile_data_not_shift = global::hp_profiles_not_shift[this_hp_profile_idx];
    profile_shiftable_cumsum = global::hp_profiles_s_cumsum[this_hp_profile_idx];
    // further initialization
    currentDemand_kW = 0;
    total_consumption_kWh = 0.0;
    cweek_consumption_kWh = 0.0;
    // computation of rated power (without AUX heating mode)
    // = max of shiftable demand time series
    rated_power_kW = 0.0;
    for (unsigned long tsID = 0; tsID < Global::get_n_timesteps(); tsID++) {
        float np = profile_data_shiftable[tsID] * scaling_factor;
        if (np > rated_power_kW)
            rated_power_kW = np;
    }
#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    state_s1 = false;
#endif
}

void ComponentHP::set_horizon_in_ts(unsigned int new_horizon) {
    BaseComponentSemiFlexible::set_horizon_in_ts(new_horizon);
    future_unshiftable_storage.clear();
    future_unshiftable_storage.resize(new_horizon, 0.0);
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
        // not-shiftable part
        future_unshiftable_storage[tOffset] = profile_data_not_shift[tsID + tOffset] ? tsID + tOffset < Global::get_n_timesteps() : 0.0;
        // shiftable part
        // for min profile -> left shift
        size_t maxProfilePos = tsID + tOffset + Global::get_hp_flexibility_in_ts();
        if (maxProfilePos < Global::get_n_timesteps()) {
            double new_val = profile_shiftable_cumsum[maxProfilePos] - total_consumption_kWh;
            if (new_val < 0)
                new_val = 0.0;
            future_maxE_storage[tOffset] = new_val;
            last_known_maxE_val = new_val;
        } else {
            future_maxE_storage[tOffset] = last_known_maxE_val;
        }
        // for min profile -> right shift
        if (tsID + tOffset < Global::get_hp_flexibility_in_ts()) {
            future_minE_storage[tOffset] = 0.0;
        } else {
            size_t minProfilePos = tsID + tOffset - Global::get_hp_flexibility_in_ts();
            if (minProfilePos < Global::get_n_timesteps()) {
                double new_val = profile_shiftable_cumsum[minProfilePos] - total_consumption_kWh;
                if (new_val < 0)
                    new_val = 0.0;
                future_minE_storage[tOffset] = new_val;
                last_known_minE_val = new_val;
            } else {
                future_minE_storage[tOffset] = last_known_minE_val;
            }
        }
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
    currentDemand_kW = ( profile_data_shiftable[tsID] + profile_data_not_shift[tsID] ) * scaling_factor;
    double e = currentDemand_kW * Global::get_time_step_size_in_h();
    total_consumption_kWh += e;
    cweek_consumption_kWh += e;
}

void ComponentHP::setDemandToGivenValue(float new_demand_kW) {
#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    if (state_s1) {
        state_s1 = false;
    } else {
        throw std::runtime_error("Method ComponentHP::setDemandToGivenValue() cannot be called at the moment! Call ComponentHP::computeNextInternalState() first.");
    }
#endif
    double e = new_demand_kW * Global::get_time_step_size_in_h();
    double new_total_e = total_consumption_kWh + e;
    // check, if the new demand is within the min/max bands
    double unshiftable_e = future_unshiftable_storage[0] * Global::get_time_step_size_in_h();
    if (new_total_e - unshiftable_e > future_maxE_storage[0] ||
        new_total_e - unshiftable_e < future_minE_storage[0]
    ) {
        std::cerr << "Warning in Component HP: Heat consumption violating bounds!" << std::endl;
    }
    // update current demand
    currentDemand_kW = new_demand_kW;
    // update cumulative variables
    total_consumption_kWh += e;
    cweek_consumption_kWh += e;
}

void ComponentHP::resetWeeklyCounter() {
    cweek_consumption_kWh = 0.0;
}

void ComponentHP::resetInternalState() {
    currentDemand_kW = 0.0;
    total_consumption_kWh = 0.0;
    cweek_consumption_kWh = 0.0;
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
    installation_place(calling_control_unit)
{
    enabled = false;
    current_demand_kW = 0.0;
    total_consumption_kWh  = 0.0;
    cweek_consumption_kWh  = 0.0;

    // Computation of maximum charging power
    max_charging_power = 11.0 * number_of_flats;

#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    is_callable_setCarStatesForTimeStep = true;
    is_callable_set_charging_value      = false;
#endif
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

void ComponentCS::add_ev(unsigned long carID) {
    listOfEVs.push_back(new EVFSM(carID, this));
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
        float remaining_power_kW = Global::get_cs_max_charging_power_kW();
        std::vector<double> new_power_per_EV_kW(listOfEVs.size(), 0.0);
        // Loop over all EVs and get min demand
        for (size_t ev_idx = 0; ev_idx < listOfEVs.size(); ev_idx++) {
            EVFSM* ev = listOfEVs[ev_idx];
            float min_demand = ev->get_future_min_consumption_kWh()->at(0) / Global::get_time_step_size_in_h(); // value at ts 0 always exists, as Global::control_horizon_in_ts >= 1 always holds!
            new_power_per_EV_kW[ev_idx] = min_demand;
            remaining_power_kW -= min_demand;
        }
        // Loop over all EVs again and assign max demand, if remaining power > 0
        for (size_t ev_idx = 0; ev_idx < listOfEVs.size(); ev_idx++) {
            EVFSM* ev = listOfEVs[ev_idx];
            float max_demand = ev->get_future_max_consumption_kWh()->at(0) / Global::get_time_step_size_in_h();
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

void ComponentCS::setDemandToGivenValues(std::vector<float>& charging_power_per_EV_kW) {

#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    if (is_callable_set_charging_value) {
        is_callable_set_charging_value = false;
        is_callable_setCarStatesForTimeStep = true;
    } else {
        throw std::runtime_error("Method ComponentCS::setDemandToGivenValues() cannot be called at the moment!");
    }
#endif

    current_demand_kW = 0.0;
    for (size_t ev_idx = 0; ev_idx < listOfEVs.size() && ev_idx < charging_power_per_EV_kW.size(); ev_idx++) {
        listOfEVs[ev_idx]->setDemandToGivenValue( charging_power_per_EV_kW[ev_idx] );
        current_demand_kW += listOfEVs[ev_idx]->get_currentDemand_kW();
    }
    // compute current demand and cumulative sum
    double e = current_demand_kW * Global::get_time_step_size_in_h();
    total_consumption_kWh += e;
    cweek_consumption_kWh += e;

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
    battery = new ComponentBS(Global::get_ev_battery_size_kWh(), 0.0f, Global::get_ev_charging_effi(), 1.0f);
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
        throw runtime_error("Error when adding a new vehicle tour: A weekday > 6 is not possible!");
    
    if (ts_duration == 0)
        ts_duration = 1;

    unsigned long this_new_tour_id = list_of_all_tours.size();
    if (list_of_all_tours.size() > 0) {
        WeeklyVehicleTour* last_know_tour = list_of_all_tours.back();
        // Check, if tours are added in the wrong ordering
        if (last_know_tour->day_of_week > weekday) {
            throw runtime_error("Error when adding a new vehicle tour: Weekday of new tour is bevore the latest added tour!");
        } else if (last_know_tour->day_of_week == weekday && last_know_tour->departure_ts_of_day > departure_ts_of_day) {
            throw runtime_error("Error when adding a new vehicle tour: Time step of departure of new tour is bevore the latest added tour!");
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
    double last_minE_value = 0.0;
    double last_maxE_value = 0.0; // = cumulative sum of electricity used for driving
    double cumsum_driving_distance_km = 0.0;
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
    double current_min_cumsum_SOE_kWh = battery->get_maxE_kWh(); // Start with a full battery
    for (unsigned long ts = Global::get_first_timestep(); ts <= Global::get_last_timestep(); ts++) {
        unsigned long tsID = ts - 1;
        ts_since_last_connection++;
        // is the currently active tour finished?
        if (current_sTour != NULL && current_sTour->ts_arrival <= ts) {
            // arrival of current tour
            last_maxE_value += current_sTour->energy_consumption_kWh;
            cumsum_driving_distance_km += current_sTour->weekly_tour->tour_length_km;
            // Computation of EVState: Connected or not connected at home?
            if (Global::get_ev_plugin_probability() >= 1.0) {
                current_state = EVState::ConnectedAtHome;
                ts_since_last_connection = 0;
            } else {
                // Pluggin-in required, if: SOC <= 0.35, or if SOC is too low for next tour, or if sampling says so
                if (
                    battery->get_SOC() <= 0.35 ||
                    // TODO: SOC is too low for the next tour ... is this really required ???
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
            energy_demand_per_tour_ts = (float) (current_sTour->weekly_tour->tour_length_km * econs_kWh_per_km) / (float) (current_sTour->weekly_tour->ts_duration);
            // Compute new min SOE to fulfil this new tour, i.e., max(0.35*BAT_CAPACITY, THIS_TOUR_CONSUMPTION)
            double min_req_SOE = 0.35 * battery->get_maxE_kWh();
            if (min_req_SOE < current_sTour->energy_consumption_kWh) {
                min_req_SOE = current_sTour->energy_consumption_kWh;
            }
            if (ev_fully_charged_at_next_dep) {
                min_req_SOE = battery->get_maxE_kWh();
                // check, if there is enough time to charge this amount, otherwise reduce it!
                double max_possible_kWh = ts_since_last_connection * Global::get_ev_max_charging_power_kW() * Global::get_time_step_size_in_h();
                if (max_possible_kWh < min_req_SOE) {
                    min_req_SOE = max_possible_kWh;
                }
            }
            current_min_cumsum_SOE_kWh -= min_req_SOE;
            // if min cumsum SOE would be smaller 0.0, we need to charge ...
            if (current_min_cumsum_SOE_kWh < 0.0) {
                battery->set_SOE_without_computations(-current_min_cumsum_SOE_kWh);
                last_minE_value += -current_min_cumsum_SOE_kWh;
                current_min_cumsum_SOE_kWh = 0.0;
            }
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
        // Set current available max charging power
        if (current_state == EVState::ConnectedAtHome) {
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
    for (size_t tOffset = 0; tOffset < Global::get_control_horizon_in_ts(); tOffset++) {
        size_t tsIDaOffset = tsID + tOffset;
        if (tsIDaOffset < Global::get_n_timesteps()) {
            last_known_maxE_val = prec_vec_of_maxE[tsIDaOffset] - sum_of_E_charged_home_kWh;
            last_known_minE_val = prec_vec_of_minE[tsIDaOffset] - sum_of_E_charged_home_kWh;
            if (last_known_maxE_val < 0) // TODO: How can this happen?
                last_known_maxE_val = 0.0;
            if (last_known_minE_val < 0)
                last_known_minE_val = 0.0;
        }
        future_maxE_storage[tOffset] = last_known_maxE_val;
        future_minE_storage[tOffset] = last_known_minE_val;
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
        output::outputEVStateDetails(ts, carID, current_state, current_P_kW, (float) sum_of_E_charged_home_kWh, (float) prec_vec_of_minE[current_ts-1], (float) prec_vec_of_maxE[current_ts-1], (float) battery->get_currentCharge_kWh());
    }
}

// define a small value added to the min/max cumsum energy consumption to remove misplaced warnings due to rounding errors
#define epsilon 0.001

void EVFSM::setDemandToGivenValue(float new_demand_kW) {
#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
    if (!state_s3) {
        throw std::runtime_error("Method EVFSM::setDemandToGivenValue() cannot be called at the moment!");
    }
    state_s2 = true;
    state_s3 = false;
#endif
    if (new_demand_kW < 0) {
        std::cerr << "Warning: Bidirection EV charging currently not supported for EVFSM::setDemandToGivenValue()." << std::endl;
        return;
    }
    double e = new_demand_kW * Global::get_time_step_size_in_h();
    double new_total_e = sum_of_E_charged_home_kWh + e;
    // check, if the new demand is within the min/max bands
    if (new_total_e > prec_vec_of_maxE[current_ts-1] + epsilon ||
        new_total_e < prec_vec_of_minE[current_ts-1] - epsilon
    ) {
        std::cerr << "Warning in EVFSM: Required demand is violating bounds for carID " << carID << " at time step " << current_ts << "!" << std::endl;
    }
    // update current demand
    battery->set_chargeRequest( new_demand_kW );
    battery->calculateActions();
    current_P_kW = new_demand_kW;
    // update cumulative variables
    sum_of_E_charged_home_kWh += e;

    if (Global::get_create_ev_detailed_output()) {
        output::outputEVStateDetails(current_ts, carID, current_state, current_P_kW, (float) sum_of_E_charged_home_kWh, (float) prec_vec_of_minE[current_ts-1], (float) prec_vec_of_maxE[current_ts-1], (float) battery->get_currentCharge_kWh());
    }
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

