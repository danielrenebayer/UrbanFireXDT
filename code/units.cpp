
#include "units.h"

#include <algorithm> /* using min() */
#include <atomic>
#include <condition_variable>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <list>
#include <mutex>
#include <optional>
#include <vector>

#include <cmath> /* using: pow() */

#include "global.h"
#include "output.h"
#include "sac_planning.h"
#include "surplus_controller.hpp"

#include "optimization_unit_general.hpp"
#ifdef USE_OR_TOOLS
    #include "optimization_unit_or_tools.hpp"
#elif USE_GUROBI
    #include "optimization_unit_gurobi.hpp"
#endif

using namespace std;


// ----------------------------- //
//      Implementation of        //
//          Substation           //
// ----------------------------- //

size_t Substation::st__n_substations = 0;
std::vector<Substation*> Substation::st__substation_list;
std::map<unsigned long, unsigned long> Substation::public_to_internal_id;

bool Substation::InstantiateNewSubstation(unsigned long public_id, std::string* name) {
    // check, if public_id is known
    if (public_to_internal_id.contains(public_id))
        return false;
    // define new internal id
    unsigned long new_internal_id = st__n_substations;
    // increment the list
    st__n_substations++;
    // register public id
    public_to_internal_id.insert(std::pair<unsigned long, unsigned long>(public_id, new_internal_id));
    // create instance
    Substation* new_obj = new Substation(new_internal_id, public_id, name);
    // add to list of instances
    st__substation_list.push_back(new_obj);
    //
    return true;
}

Substation::Substation(unsigned long internal_id, unsigned long public_id, std::string* name)
	: internal_id(internal_id), public_id(public_id), name(name)
{
	/*
	 * Constructor implementation of class Substation.
	 */
    // initialize instance variables
    connected_units = new list<ControlUnit*>();
    station_load    = 0.0;
    resident_load   = 0.0;
    resident_demand = 0.0;
    total_demand    = 0.0;
    total_BESS_demand = 0.0;
    pv_gen_total_kW       = 0.0;
    bs_gen_total_kW       = 0.0;
    chp_gen_total_kW      = 0.0;
    wind_gen_total_kW     = 0.0;
    unknown_gen_total_kW  = 0.0;
    pv_gen_expo_kW        = 0.0;
    bs_gen_expo_kW        = 0.0;
    chp_gen_expo_kW       = 0.0;
    wind_gen_expo_kW      = 0.0;
    unknown_gen_expo_kW   = 0.0;
}

Substation::~Substation() {
	delete name;
	delete connected_units;
}

void Substation::add_unit(ControlUnit* unit) {
	connected_units->push_back(unit);
}

void Substation::calc_load() {
    station_load    = 0.0;
    resident_load   = 0.0; // Load only of residential buildings
    resident_demand = 0.0; // Residential load, only demand
    total_demand    = 0.0;
    total_BESS_demand = 0.0;
    pv_gen_total_kW       = 0.0;
    bs_gen_total_kW       = 0.0;
    chp_gen_total_kW      = 0.0;
    wind_gen_total_kW     = 0.0;
    unknown_gen_total_kW  = 0.0;
    pv_gen_expo_kW        = 0.0;
    bs_gen_expo_kW        = 0.0;
    chp_gen_expo_kW       = 0.0;
    wind_gen_expo_kW      = 0.0;
    unknown_gen_expo_kW   = 0.0;
    for (ControlUnit* cu : *connected_units) {
        station_load += cu->get_current_load_vSMeter_kW();
        total_demand += cu->get_current_demand_wo_BS_or_gen_kW();
        total_BESS_demand    += cu->get_current_BS_demand_kW();
        pv_gen_total_kW      += cu->get_current_PV_generation_kW();
        bs_gen_total_kW      += cu->get_current_BS_generation_kW();
        chp_gen_total_kW     += cu->get_current_CHP_generation_kW();
        wind_gen_total_kW    += cu->get_current_wind_generation_kW();
        unknown_gen_total_kW += cu->get_current_unknown_generation_kW();
        if (cu->is_residential()) {
            resident_load += cu->get_current_load_vSMeter_kW();
            if (cu->get_current_load_vSMeter_kW() > 0) {
                resident_demand += cu->get_current_load_vSMeter_kW();
            }
        }
        // In case of a feed-in, detect the source of the feed-in
        if (cu->get_current_load_vSMeter_kW() < 0.0) {
            pv_gen_expo_kW       += cu->get_current_PV_feedin_to_grid_kW();
            bs_gen_expo_kW       += cu->get_current_BS_feedin_to_grid_kW();
            chp_gen_expo_kW      += cu->get_current_CHP_feedin_to_grid_kW();
            wind_gen_expo_kW     += cu->get_current_wind_feedin_to_grid_kW();
            unknown_gen_expo_kW  += cu->get_current_unknown_feedin_to_grid_kW();
        }
    }
}

void Substation::InitializeStaticVariables(unsigned long n_substations) {
    st__substation_list.reserve(n_substations);
}

void Substation::VacuumInstancesAndStaticVariables() {
	// delete all instances
    for (size_t i = 0; i < st__n_substations; i++) {
        delete st__substation_list[i];
        st__substation_list[i] = NULL;
    }
	st__n_substations = 0;
    st__substation_list.clear();
    public_to_internal_id.clear();
}

Substation* Substation::GetInstancePublicID(unsigned long public_id) {
    try {
        size_t internal_id = public_to_internal_id.at(public_id); // throws std::out_of_range if public_id is unknown
		return st__substation_list[internal_id];
    } catch (std::out_of_range const&) {
		return NULL;
    }
}




// ----------------------------- //
//      Implementation of        //
//         ControlUnit           //
// ----------------------------- //

size_t ControlUnit::st__n_CUs = 0;
std::vector<ControlUnit*> ControlUnit::st__cu_list;
std::map<unsigned long, unsigned long> ControlUnit::public_to_internal_id;
std::vector<double>* ControlUnit::st__empty_vector_for_time_horizon = NULL;
const std::string ControlUnit::MetricsStringHeaderAnnual = "UnitID,SCR,SSR,NPV,ALR,BDR,RBC,Sum of demand [kWh],Sum of MU demand [kWh],Sum of self-consumed e. [kWh],Sum of PV-generated e. [kWh],Sum of grid feed-in [kWh],Sum of grid demand [kWh],BS EFC,BS n_ts_empty,BS n_ts_full,BS total E withdrawn [kWh],BS total E withdrawn grid-originated [kWh],Sum of HP demand [kWh],Sum of CS demand [kWh],Peak grid demand [kW],Emissions cbgd [kg CO2eq],Avoided emissions [kg CO2eq],Actual electricity costs [CU],Electricity cons. costs [CU],Avoided electricity cons. costs [CU],Feed-in revenue [CU],Sim. PV max P [kWp],Sim. BS P [kW],Sim. BS E [kWh],n EVs,Sim. CS max P [kW],Simulated PV,Simulated BS,Simulated HP,Simulated CS,n errors in cntrl,n error cntrl cmd appl";
const std::string ControlUnit::MetricsStringHeaderWeekly = "UnitID,Week number,SCR,SSR,Sum of demand [kWh],Sum of MU demand [kWh],Sum of self-consumed e. [kWh],Sum of PV-generated e. [kWh],Sum of grid feed-in [kWh],Sum of grid demand [kWh],BS EFC,BS total E withdrawn [kWh],BS total E withdrawn grid-originated [kWh],Sum of HP demand [kWh],Sum of CS demand [kWh],Peak grid demand [kW],Emissions cbgd [kg CO2eq],Avoided emissions [kg CO2eq],Actual electricity costs [CU],Electricity cons. costs [CU],Avoided electricity cons. costs [CU],Feed-in revenue [CU],Sim. PV max P [kWp],Sim. BS P [kW],Sim. BS E [kWh],n EVs,Sim. CS max P [kW],Simulated PV,Simulated BS,Simulated HP,Simulated CS,n errors in cntrl,n error cntrl cmd appl";
std::atomic<unsigned long> ControlUnit::optimization_call_counter = 0;
std::atomic<unsigned long> ControlUnit::n_curr_optimization_vars = 0;
std::mutex ControlUnit::n_curr_optimization_vars_mutex;
std::condition_variable ControlUnit::curr_opti_cv;

bool ControlUnit::InstantiateNewControlUnit(unsigned long public_unitID, unsigned long substation_id, unsigned long locationID, bool residential, unsigned int n_flats) {
    // check, if public_id is known
    if (public_to_internal_id.contains(public_unitID))
        return false;
    // define new internal id
    unsigned long new_internal_id = st__n_CUs;
    // increment the list
    st__n_CUs++;
    // register public id
    public_to_internal_id.insert(std::pair<unsigned long, unsigned long>(public_unitID, new_internal_id));
    // create instance
    ControlUnit* new_obj = new ControlUnit(new_internal_id, public_unitID, substation_id, locationID, residential, n_flats);
    // add to list of instances
    st__cu_list.push_back(new_obj);
    //
    return true;
}

ControlUnit::ControlUnit(unsigned long internalID, unsigned long publicID, unsigned long substation_id, unsigned long locationID, bool residential, unsigned int n_flats)
    : internal_id(internalID), unitID(publicID), higher_level_subst(Substation::GetInstancePublicID(substation_id)), locationID(locationID), residential(residential),
      generate_output_for_this_unit( global::unitIDs_selected_for_output.size() == 0 || std::find( global::unitIDs_selected_for_output.begin(), global::unitIDs_selected_for_output.end(), publicID ) != global::unitIDs_selected_for_output.end() )
{
	//
	// initialize instance variables
	connected_units = new list<MeasurementUnit*>();
	has_sim_pv      = false;
	has_sim_bs      = false;
	has_sim_hp      = false;
    has_sim_cs      = false;
	sim_comp_pv     = NULL;
	sim_comp_bs     = NULL;
	sim_comp_hp     = NULL;
	current_load_vSM_kW   = 0;
	self_produced_load_kW = 0;
    current_total_consumption_kW = 0.0;
    current_PV_feedin_to_grid_kW      = 0.0;
    current_BS_feedin_to_grid_kW      = 0.0;
    current_CHP_feedin_to_grid_kW     = 0.0;
    current_wind_feedin_to_grid_kW    = 0.0;
    current_unknown_feedin_to_grid_kW = 0.0;

    output_obj      = NULL;
    is_expandable_with_pv_cache          = false;
    is_expandable_with_pv_cache_computed = false;
    is_expandable_with_hp_cache          = false;
    is_expandable_with_hp_cache_computed = false;

#ifdef PYTHON_MODULE
    py_control_commands_obtained = false;
    py_cmd_p_bs_kW = 0.0;
    py_cmd_p_hp_kW = 0.0;
    py_cmd_p_ev_kW.clear();
#endif

	//
	// add this control unit to the list of
	// connected units in the connected substation
    if (higher_level_subst == NULL) {
        cerr << "Substation with ID " << substation_id << " is not available!" << endl;
        throw runtime_error("Invalid substation ID given to constructor of ControlUnit class!");
    }
	higher_level_subst->add_unit(this);

    //
    // For the computation of the evaluation metrics (like SSC/SSR/NPV)
    // one needs to store some information about every time step for the current
    // parameter variation setting.
    sum_of_consumption_kWh    = 0.0;
    sum_of_self_cons_kWh      = 0.0;
    sum_of_mu_cons_kWh        = 0.0;
    sum_of_feed_into_grid_kWh = 0.0;
    sum_of_grid_demand_kWh    = 0.0;
    sum_of_rem_pow_costs_EUR       = 0.0;
    sum_of_saved_pow_costs_EUR     = 0.0;
    sum_of_feedin_revenue_EUR      = 0.0;
    sum_of_emissions_cbgd_kg_CO2eq = 0.0;
    sum_of_emissions_avoi_kg_CO2eq = 0.0;
    sum_of_errors_in_cntrl         = 0;
    sum_of_errors_in_cntrl_cmd_appl= 0;
    peak_grid_demand_kW            = 0.0;
    // weekly counters
    sum_of_cweek_consumption_kWh         = 0.0;
    sum_of_cweek_self_cons_kWh           = 0.0;
    sum_of_cweek_mu_cons_kWh             = 0.0;
    sum_of_cweek_feed_into_grid_kWh      = 0.0;
    sum_of_cweek_grid_demand_kWh         = 0.0;
    sum_of_cweek_rem_pow_costs_EUR       = 0.0;
    sum_of_cweek_saved_pow_costs_EUR     = 0.0;
    sum_of_cweek_feedin_revenue_EUR      = 0.0;
    sum_of_cweek_emissions_cbgd_kg_CO2eq = 0.0;
    sum_of_cweek_emissions_avoi_kg_CO2eq = 0.0;
    sum_of_cweek_errors_in_cntrl         = 0;
    sum_of_cweek_errors_in_cntrl_cmd_appl= 0;
    cweek_peak_grid_demand_kW            = 0.0;

    is_sim_expanded = false;

    // Generate new instance for the EV charging station (regardless if it is required or not)
    sim_comp_cs = new ComponentCS(this, n_flats);

    ts_since_last_opti_run = Global::get_control_update_freq_in_ts();

    worker_thread = NULL;

    // Initialize the controller to NULL by default
    optimized_controller = NULL;
}

ControlUnit::~ControlUnit() {
	delete connected_units;
	if (has_sim_pv) delete sim_comp_pv;
	if (has_sim_bs) delete sim_comp_bs;
	if (has_sim_hp) delete sim_comp_hp;
    delete sim_comp_cs;
    if (optimized_controller != NULL)
        delete optimized_controller;
}

void ControlUnit::add_unit(MeasurementUnit* unit) {
	connected_units->push_back(unit);
}

bool ControlUnit::has_electricity_demand() {
	for (MeasurementUnit* mu : *connected_units) {
		if (mu->has_demand())
			return true;
	}
	return false;
}

bool ControlUnit::has_pv() {
	if (has_sim_pv)
		return true;
	for (MeasurementUnit* mu : *connected_units) {
		if (mu->has_pv())
			return true;
	}
	return false;
}

bool ControlUnit::has_bs() {
	if (has_sim_bs)
		return true;
	for (MeasurementUnit* mu : *connected_units) {
		if (mu->has_bs())
			return true;
	}
	return false;
}

bool ControlUnit::has_hp() {
	if (has_sim_hp)
		return true;
	for (MeasurementUnit* mu : *connected_units) {
		if (mu->has_hp())
			return true;
	}
	return false;
}

bool ControlUnit::has_cs() {
	if (has_sim_cs)
		return true;
	for (MeasurementUnit* mu : *connected_units) {
		if (mu->has_evchst())
			return true;
	}
	return false;
}

bool ControlUnit::has_chp() {
    for (MeasurementUnit* mu : *connected_units) {
        if (mu->has_chp())
            return true;
    }
    return false;
}

bool ControlUnit::has_wind() {
    for (MeasurementUnit* mu : *connected_units) {
        if (mu->has_wind())
            return true;
    }
    return false;
}

bool ControlUnit::has_bs_sim_added() {
	return has_sim_bs;
}

bool ControlUnit::is_expandable_with_pv() {
    if (!is_expandable_with_pv_cache_computed) {
        // compute value first, if it has not been computed
        is_expandable_with_pv_cache_computed = true;
        // computation itselfe
        bool geodata_flag_set  = global::locations_with_geodata.contains(locationID);
        bool roofdata_flag_set = global::roof_section_orientations.contains(locationID);
        if (geodata_flag_set && roofdata_flag_set) {
            is_expandable_with_pv_cache = true;
        } else if (geodata_flag_set && !roofdata_flag_set) {
            is_expandable_with_pv_cache = false;
            cerr << "Warning: Geodata is claimed to be available for control unit with ID " << unitID << ", but still no roof data has been found!" << endl;
        } else if (!geodata_flag_set && roofdata_flag_set) {
            is_expandable_with_pv_cache = true;
            cerr << "Warning: Roofdata is available for control unit with ID " << unitID << ", even though it is claimed that there is no geodata for it!" << endl;
        } else {
            is_expandable_with_pv_cache = false;
        }
    }

    return is_expandable_with_pv_cache;
}

bool ControlUnit::is_expandable_with_hp() {
    if (!is_expandable_with_hp_cache_computed) {
        // compute value first, if it has not been computed
        is_expandable_with_hp_cache_computed = true;
        // defaults to true
        is_expandable_with_hp_cache = true;
        // check, if annual heat consumption > 0, otherwise, no HP can be added!
        // if heat consumption above the global limit for HP addition, also do not include this
        if (
            get_annual_heat_demand_th_kWh() <= 0.0 ||
            (
                Global::get_annual_heat_demand_limit_fsac() >= 1.0 &&
                get_annual_heat_demand_th_kWh() > Global::get_annual_heat_demand_limit_fsac()
            )
        ) {
            is_expandable_with_hp_cache = false;
        }
        // only select units with known heat consumption?
        if (Global::get_select_buildings_wg_heatd_only() && !heat_demand_given_in_data())
            is_expandable_with_hp_cache = false;
        // if there is an chp, exclude this control unit from HP expansion
        if (has_chp()) {
            is_expandable_with_hp_cache = false;
        }
    }

    return is_expandable_with_hp_cache;
}

bool ControlUnit::heat_demand_given_in_data() {
    return global::annual_heat_demand_kWh[locationID] > 0;
}

int ControlUnit::get_exp_combi_bit_repr() {
    return get_exp_combi_bit_repr_from_MUs() | get_exp_combi_bit_repr_sim_added();
}

int ControlUnit::get_exp_combi_bit_repr_from_MUs() {
    //
    // This function returns the expansion combination in bitwise representation
    // as given by the sum of the connected measurement units
    //
    int combination = 0;
    for (MeasurementUnit* mu : *connected_units) {
        combination = combination | mu->get_expansion_combination();
    }
    return combination;
}

int ControlUnit::get_exp_combi_bit_repr_sim_added() {
    //
    // This function returns the expansion combination in bitwise representation
    // that is added by the simulation
    //
    int combination = 0;
    if (has_sim_pv)
        combination = combination | expansion::MaskPV;
    if (has_sim_bs)
        combination = combination | expansion::MaskBS;
    if (has_sim_hp)
        combination = combination | expansion::MaskHP;
    if (has_sim_cs)
        combination = combination | expansion::MaskWB;
    return combination;
}

double ControlUnit::get_current_BS_demand_kW() const {
    if (has_sim_bs) {
        double cl = sim_comp_bs->get_currentLoad_kW();
        if (cl > 0.0)
            return cl;
    }
    return 0.0;
}

double ControlUnit::get_current_PV_generation_kW() const {
    if (has_sim_pv) {
        return sim_comp_pv->get_currentGeneration_kW();
    }
    return 0.0;
}

double ControlUnit::get_current_BS_generation_kW() const {
    if (has_sim_bs) {
        double cl = sim_comp_bs->get_currentLoad_kW();
        if (cl < 0.0)
            return -cl;
    }
    return 0.0;
}

double ControlUnit::get_current_CHP_generation_kW() const {
    double current_generation_kW = 0.0;
    for (MeasurementUnit* mu : *connected_units) {
        if (mu->has_chp()) {
            double cf = mu->get_current_ts_rsm_value();
            if (cf < 0.0) {
                current_generation_kW += -cf;
            }
        }
    }
    return current_generation_kW;
}
double ControlUnit::get_current_wind_generation_kW() const {
    double current_generation_kW = 0.0;
    for (MeasurementUnit* mu : *connected_units) {
        if (mu->has_wind()) {
            double cf = mu->get_current_ts_rsm_value();
            if (cf < 0.0) {
                current_generation_kW += -cf;
            }
        }
    }
    return current_generation_kW;
}

double ControlUnit::get_current_unknown_generation_kW() const {
    return 0.0;
}

double ControlUnit::get_mean_annual_MU_el_demand_kWh() const {
    double total_demand = 0.0;
    for (MeasurementUnit* mu : *connected_units) {
        if (mu->has_demand()) {
            total_demand += mu->get_total_demand_kWh();
        }
    }
    // compute the complete time span length of the input data (in hours)
    double data_len_in_h = (double) (Global::get_n_timesteps()) * Global::get_time_step_size_in_h();
    double data_len_in_years = data_len_in_h / ( 365 * 24 );
    return total_demand / data_len_in_years;
}

double ControlUnit::get_sim_comp_pv_kWp() {
    if (has_sim_pv)
        return sim_comp_pv->get_kWp();
    return 0;
}

double ControlUnit::get_sim_comp_bs_P_kW() {
    if (has_sim_bs)
        return sim_comp_bs->get_maxP_kW();
    return 0;
}

double ControlUnit::get_sim_comp_bs_E_kWh() {
    if (has_sim_bs)
        return sim_comp_bs->get_maxE_kWh();
    return 0;
}

float ControlUnit::get_annual_heat_demand_th_kWh() {
    float heat_demand = global::annual_heat_demand_kWh[locationID];
    if (heat_demand <= 0) {
        float building_volume = global::building_volumes_m3[locationID];
        if (building_volume > 0.0) {
            // heat demand not given by gas or other data, so we need to approximate it based on a linear regression
            return Global::get_heat_cons_bobv_slope() * building_volume + Global::get_heat_cons_bobv_intercept();
        } else {
            return 0.0;
        }
    }
    return heat_demand;
}

float ControlUnit::get_annual_hp_el_cons_kWh() {
    return get_annual_heat_demand_th_kWh() / Global::get_heat_demand_thermalE_to_hpE_conv_f();
}

double ControlUnit::get_sim_comp_cs_max_P_kW() const {
    return sim_comp_cs->get_max_P_kW();
}

size_t ControlUnit::get_sim_comp_cs_n_EVs() const {
    return sim_comp_cs->get_n_EVs();
}

size_t ControlUnit::get_sim_comp_cs_possible_n_EVs() const {
    return sim_comp_cs->get_possible_n_EVs();
}

double ControlUnit::get_SSR() {
    double SSR = 0.0;
    if (sum_of_consumption_kWh > 0)
        SSR = sum_of_self_cons_kWh / sum_of_consumption_kWh;
    return SSR;
}

double ControlUnit::get_SCR() {
    double SCR = 0.0;
    if (has_sim_pv && sim_comp_pv->get_total_generation_kWh() > 0)
        SCR = sum_of_self_cons_kWh / sim_comp_pv->get_total_generation_kWh();
    return SCR;
}

double ControlUnit::get_NPV() {
    // 0. return 0 if no PV or battery is available
    if (!(has_sim_pv || has_sim_bs))
        return 0;
    // 1. compute investment cost
    double investemnt_costs = 0.0;
    if (has_sim_pv)
        investemnt_costs += sim_comp_pv->get_kWp() * Global::get_inst_cost_PV_per_kWp();
    if (has_sim_bs)
        investemnt_costs += sim_comp_bs->get_maxE_kWh() * Global::get_inst_cost_BS_per_kWh();
    if (has_sim_hp)
        investemnt_costs += (double) sim_comp_hp->get_rated_power_without_AUX() * Global::get_inst_cost_HP_per_kW();
    if (has_sim_cs)
        investemnt_costs += (double) sim_comp_cs->get_n_chargers() * Global::get_inst_cost_CS_per_unit();
    // 2. compute savings per year
    //const double& saved_energy_kWh = sum_of_self_cons_kWh; // energy in kWh for which no grid demand was required
    double savings_per_year = sum_of_saved_pow_costs_EUR + sum_of_feedin_revenue_EUR;
    // 3. compute total savings
    double total_savings = savings_per_year * Global::get_npv_factor_if_const(); // new computation rule using the present value of annuity formula
    /*for (unsigned int t = 1; t <= Global::get_npv_time_horizon(); t++) {
        total_savings += savings_per_year / pow( (double)(1 + Global::get_npv_discount_rate()), (double) t );
    }*/
    return - investemnt_costs + total_savings;
}

string* ControlUnit::get_metrics_string_annual() {
        double SCR = get_SCR();
        double SSR = get_SSR();
        double NPV = get_NPV();
        double sum_of_PV_generated_kWh = 0.0;
        if (has_sim_pv) {
            sum_of_PV_generated_kWh = sim_comp_pv->get_total_generation_kWh();
        }
        double bat_EFC = 0.0;
        double bat_E_withdrawn = 0.0;
        unsigned long bat_n_ts_empyt = 0;
        unsigned long bat_n_ts_full  = 0;
        if (has_sim_bs) {
            bat_EFC = sim_comp_bs->get_current_EFC();
            bat_E_withdrawn= sim_comp_bs->get_total_withdrawn_E_kWh();
            bat_n_ts_empyt = sim_comp_bs->get_n_ts_empty();
            bat_n_ts_full  = sim_comp_bs->get_n_ts_full();
        }
        // computation of the array to load ratio (see Nyholm et al. 2016 in Applied Energy, https://doi.org/10.1016/j.apenergy.2016.08.172)
        double ALR = 0.0;
        if (has_sim_pv) {
            ALR = sim_comp_pv->get_kWp() / ( sum_of_consumption_kWh / (double) Global::get_n_timesteps() / Global::get_time_step_size_in_h() );
        }
        // compuation of the Battery-to-Demand ratio
        // and the Relative Battery Capacity (RBC)
        double BDR = 0.0;
        double RBC = 0.0;
        if (has_sim_bs) {
            BDR = sim_comp_bs->get_maxE_kWh() / ( sum_of_consumption_kWh / (double) Global::get_n_timesteps() / Global::get_time_step_size_in_h() );
            if (has_sim_pv) {
                RBC = 1000.0 * sim_comp_bs->get_maxE_kWh() / sim_comp_pv->get_total_generation_kWh();
            }
        }
        // See ControlUnit::MetricsStringHeaderAnnual for the header definition
        //
        string* retstr = new string;
        *retstr += to_string(unitID) + ",";
        *retstr += to_string(SCR) + ",";
        *retstr += to_string(SSR) + ",";
        *retstr += to_string(NPV) + ",";
        *retstr += to_string(ALR) + ",";
        *retstr += to_string(BDR) + ",";
        *retstr += to_string(RBC) + ",";
        *retstr += to_string(sum_of_consumption_kWh) + ",";
        *retstr += to_string(sum_of_mu_cons_kWh)     + ",";
        *retstr += to_string(sum_of_self_cons_kWh)   + ",";
        *retstr += to_string(sum_of_PV_generated_kWh)+ ",";
        *retstr += to_string(sum_of_feed_into_grid_kWh)+ ",";
        *retstr += to_string(sum_of_grid_demand_kWh)+ ",";
        *retstr += to_string(bat_EFC) + ",";
        *retstr += to_string(bat_n_ts_empyt) + ",";
        *retstr += to_string(bat_n_ts_full) + ",";
        *retstr += to_string(bat_E_withdrawn) + ",";
        *retstr += to_string((has_sim_bs) ? sim_comp_bs->get_total_withdrawn_E_gridOnly_kWh() : 0.0) + ",";
        *retstr += to_string((has_sim_hp) ? sim_comp_hp->get_total_consumption_kWh() : 0.0) + ",";
        *retstr += to_string((has_sim_cs) ? sim_comp_cs->get_total_consumption_kWh() : 0.0) + ",";
        *retstr += to_string(peak_grid_demand_kW) + ",";
        *retstr += to_string(sum_of_emissions_cbgd_kg_CO2eq) + ",";
        *retstr += to_string(sum_of_emissions_avoi_kg_CO2eq) + ",";
        *retstr += to_string(sum_of_rem_pow_costs_EUR - sum_of_feedin_revenue_EUR) + ",";
        *retstr += to_string(sum_of_rem_pow_costs_EUR)       + ",";
        *retstr += to_string(sum_of_saved_pow_costs_EUR)     + ",";
        *retstr += to_string(sum_of_feedin_revenue_EUR)      + ",";
        *retstr += to_string(get_sim_comp_pv_kWp()) + ",";
        *retstr += to_string(get_sim_comp_bs_P_kW()) + ",";
        *retstr += to_string(get_sim_comp_bs_E_kWh()) + ",";
        *retstr += to_string(get_sim_comp_cs_n_EVs()) + ",";
        *retstr += to_string(get_sim_comp_cs_max_P_kW()) + ",";
        *retstr += to_string(has_sim_pv ? 1 : is_expandable_with_pv() ? 0 : -1) + ",";
        *retstr += to_string(has_sim_bs) + ",";
        *retstr += to_string(has_sim_hp ? 1 : is_expandable_with_hp() ? 0 : -1) + ",";
        *retstr += to_string(has_sim_cs) + ",";
        *retstr += to_string(sum_of_errors_in_cntrl) + ",";
        *retstr += to_string(sum_of_errors_in_cntrl_cmd_appl);
        return retstr;
}

string* ControlUnit::get_metrics_string_weekly_wr(unsigned long week_number) {
        double SCR_cweek = 0.0;
        if (has_sim_pv && sim_comp_pv->get_cweek_generation_kWh() > 0) {
            SCR_cweek = sum_of_cweek_self_cons_kWh / sim_comp_pv->get_cweek_generation_kWh();
        }
        double SSR_cweek = 0.0;
        if (sum_of_cweek_consumption_kWh > 0) {
            SSR_cweek = sum_of_cweek_self_cons_kWh / sum_of_cweek_consumption_kWh;
        }
        double sum_of_cweek_PV_generated_kWh = 0.0;
        if (has_sim_pv) {
            sum_of_cweek_PV_generated_kWh = sim_comp_pv->get_cweek_generation_kWh();
            sim_comp_pv->resetWeeklyCounter();
        }
        double bat_EFC_cweek = 0.0;
        double bat_E_cweek_withdrawn = 0.0;
        double bat_E_cweek_gridOnly_withdrawn = 0.0;
        if (has_sim_bs) {
            bat_EFC_cweek = sim_comp_bs->get_cweek_EFC();
            bat_E_cweek_withdrawn= sim_comp_bs->get_cweek_withdrawn_E_kWh();
            bat_E_cweek_gridOnly_withdrawn = sim_comp_bs->get_cweek_withdrawn_E_gridOnly_kWh();
            sim_comp_bs->resetWeeklyCounter();
        }
        // See ControlUnit::MetricsStringHeader for the header definition
        //
        string* retstr = new string;
        *retstr += to_string(unitID)      + ",";
        *retstr += to_string(week_number) + ",";
        *retstr += to_string(SCR_cweek)   + ",";
        *retstr += to_string(SSR_cweek)   + ",";
        *retstr += to_string(sum_of_cweek_consumption_kWh)   + ",";
        *retstr += to_string(sum_of_cweek_mu_cons_kWh)       + ",";
        *retstr += to_string(sum_of_cweek_self_cons_kWh)     + ",";
        *retstr += to_string(sum_of_cweek_PV_generated_kWh)  + ",";
        *retstr += to_string(sum_of_cweek_feed_into_grid_kWh)+ ",";
        *retstr += to_string(sum_of_cweek_grid_demand_kWh)   + ",";
        *retstr += to_string(bat_EFC_cweek)         + ",";
        *retstr += to_string(bat_E_cweek_withdrawn) + ",";
        *retstr += to_string(bat_E_cweek_gridOnly_withdrawn) + ",";
        *retstr += to_string((has_sim_hp) ? sim_comp_hp->get_cweek_consumption_kWh() : 0.0) + ",";
        *retstr += to_string((has_sim_cs) ? sim_comp_cs->get_cweek_consumption_kWh() : 0.0) + ",";
        *retstr += to_string(cweek_peak_grid_demand_kW) + ",";
        *retstr += to_string(sum_of_cweek_emissions_cbgd_kg_CO2eq) + ",";
        *retstr += to_string(sum_of_cweek_emissions_avoi_kg_CO2eq) + ",";
        *retstr += to_string(sum_of_cweek_rem_pow_costs_EUR - sum_of_cweek_feedin_revenue_EUR) + ",";
        *retstr += to_string(sum_of_cweek_rem_pow_costs_EUR)       + ",";
        *retstr += to_string(sum_of_cweek_saved_pow_costs_EUR)     + ",";
        *retstr += to_string(sum_of_cweek_feedin_revenue_EUR)      + ",";
        *retstr += to_string(get_sim_comp_pv_kWp())   + ",";
        *retstr += to_string(get_sim_comp_bs_P_kW())  + ",";
        *retstr += to_string(get_sim_comp_bs_E_kWh()) + ",";
        *retstr += to_string(get_sim_comp_cs_n_EVs()) + ",";
        *retstr += to_string(get_sim_comp_cs_max_P_kW()) + ",";
        *retstr += to_string(has_sim_pv) + ",";
        *retstr += to_string(has_sim_bs) + ",";
        *retstr += to_string(has_sim_hp) + ",";
        *retstr += to_string(has_sim_cs) + ",";
        *retstr += to_string(sum_of_cweek_errors_in_cntrl) + ",";
        *retstr += to_string(sum_of_cweek_errors_in_cntrl_cmd_appl);
        // reset weekly counters
        sum_of_cweek_consumption_kWh         = 0.0;
        sum_of_cweek_self_cons_kWh           = 0.0;
        sum_of_cweek_mu_cons_kWh             = 0.0;
        sum_of_cweek_feed_into_grid_kWh      = 0.0;
        sum_of_cweek_grid_demand_kWh         = 0.0;
        sum_of_cweek_rem_pow_costs_EUR       = 0.0;
        sum_of_cweek_saved_pow_costs_EUR     = 0.0;
        sum_of_cweek_feedin_revenue_EUR      = 0.0;
        sum_of_cweek_emissions_cbgd_kg_CO2eq = 0.0;
        sum_of_cweek_emissions_avoi_kg_CO2eq = 0.0;
        sum_of_cweek_errors_in_cntrl         = 0;
        sum_of_cweek_errors_in_cntrl_cmd_appl= 0;
        cweek_peak_grid_demand_kW            = 0.0;
        // reset HP and CS (PV and BS are resetted above)
        if (has_sim_hp) sim_comp_hp->resetWeeklyCounter();
        if (has_sim_cs) sim_comp_cs->resetWeeklyCounter();

        return retstr;
}

string* ControlUnit::get_pv_section_string() {
    if (has_sim_pv) {
        return sim_comp_pv->get_section_string(to_string(unitID));
    } else {
        return new string();
    }
}

void ControlUnit::add_exp_pv() {
    if (has_pv())
        cerr << "Warning: Control unit with location id " << locationID << " already has a PV installation!" << endl;
    if (!has_sim_pv) {
        if (has_sim_bs) {
            std::string error_msg = "Error in ControlUnit::add_exp_pv(): Control unit with ID " + std::to_string( unitID ) + " and location ID " + std::to_string( locationID ) + " has already a simulated BS. This method (adding PV) MUST be called before adding BS!";
            cerr << error_msg << endl;
            throw std::logic_error(error_msg);
        }
        // if no error occurs ...
        has_sim_pv  = true;
        if (Global::get_exp_pv_sizing_mode() == global::PVSizingMode::StaticPVSize) {
            sim_comp_pv = new ComponentPV(Global::get_exp_pv_kWp_static(), locationID);
        } else if (
            Global::get_exp_pv_sizing_mode() == global::PVSizingMode::MaxAvailableRoofArea ||
            Global::get_exp_pv_sizing_mode() == global::PVSizingMode::Optimized
        ) {
            sim_comp_pv = new ComponentPV(Global::get_exp_pv_kWp_per_m2(),
                                          Global::get_exp_pv_min_kWp_roof_sec(),
                                          Global::get_exp_pv_max_kWp_roof_sec(),
                                          Global::get_exp_pv_max_kWp_per_unit(),
                                          locationID);
            // Hint for the case of Global::get_exp_pv_sizing_mode() == global::PVSizingMode::Optimized
            //    --> Just instanciate a arbitrarily sized PV installation - the final sizing will be done when the optimization is called for the first time
        }
    }
}

void ControlUnit::add_exp_bs() {
    if (has_bs())
        cerr << "Warning: Control unit with location id " << locationID << " already has a battery!" << endl;
    if (!has_sim_bs) {
        // compute the battery capacity
        double new_battery_capacity_kWh = 0.0;
        if ( !has_sim_pv ||
             Global::get_battery_capacity_computation_mode() == global::BatteryCapacityComputationMode::Constant
           )
        {
            // if there is no battery -> always use static values
            if (!has_sim_pv && Global::get_battery_capacity_computation_mode() != global::BatteryCapacityComputationMode::Constant) {
                std::cerr << "Warning in Control Unit " + std::to_string(unitID) + ": No (simulated) PV installation available, but BS should be added and BS capacity computation not constant!" << std::endl;
            }
            new_battery_capacity_kWh = Global::get_exp_bess_kWh();
        } else if (Global::get_battery_capacity_computation_mode() == global::BatteryCapacityComputationMode::BasedOnNominalPVPower) {
            double pv_kWp = sim_comp_pv->get_kWp();
            new_battery_capacity_kWh = pv_kWp * Global::get_exp_bess_sizingE_boPV();
        } else if (Global::get_battery_capacity_computation_mode() == global::BatteryCapacityComputationMode::BasedOnAnnualConsumption) {
            // round on two digits
            new_battery_capacity_kWh = round( (get_mean_annual_MU_el_demand_kWh() / 1000) * 100 ) / 100.0;
        } else if (Global::get_battery_capacity_computation_mode() == global::BatteryCapacityComputationMode::BasedOnAnnualConsumptionWithHeatPump) {
            double ann_cons_kWh = get_mean_annual_MU_el_demand_kWh();
            if (has_sim_hp) {
                ann_cons_kWh += get_annual_hp_el_cons_kWh();
            }
            // round on two digits
            new_battery_capacity_kWh = round( (ann_cons_kWh / 1000) * 100 ) / 100.0;
        } else if (Global::get_battery_capacity_computation_mode() == global::BatteryCapacityComputationMode::Optimized) {
            // Just set the battery capacity to 0.0 - the actual value will be set later, when the optimization is called (for the first time)
            new_battery_capacity_kWh = 0.0;
        }
        // respect maximum addition
        if (Global::get_exp_bess_max_capacity() > 0.0) {
            if (Global::get_exp_bess_max_capacity() < new_battery_capacity_kWh) {
                new_battery_capacity_kWh = Global::get_exp_bess_max_capacity();
            }
        }
        // initialize the object
        has_sim_bs  = true;
        sim_comp_bs = new ComponentBS(new_battery_capacity_kWh /*maxE*/,
                                      Global::get_exp_bess_kW()  /*maxP*/,
                                      Global::get_exp_bess_E_P_ratio() /* E over P ratio */,
                                      Global::get_exp_bess_self_ds_ts() /* discharge rate per step */,
                                      Global::get_exp_bess_effi_in() /* efficiency in */,
                                      Global::get_exp_bess_effi_out() /* efficiency out */,
                                      Global::get_exp_bess_start_soc() /* inital SOC */);
    }
}

void ControlUnit::add_exp_hp() {
    if (has_hp())
        cerr << "Warning: Control unit with location id " << locationID << " already has a heat pump!" << endl;
    if (!has_sim_hp) {
        //
        // create and link component
        has_sim_hp  = true;
        sim_comp_hp = new ComponentHP(this, this->get_annual_hp_el_cons_kWh());
    }
}

void ControlUnit::add_exp_cs() {
    if (has_cs())
        cerr << "Warning: Control unit with location id " << locationID << " already has an EV charging station!" << endl;
    if (!has_sim_cs) {
        has_sim_cs  = true;
        sim_comp_cs->enable_station();
    }
}

void ControlUnit::preprocess_ev_data() {
    sim_comp_cs->preprocess_ev_data();
}

void ControlUnit::add_ev(unsigned long carID) {
    sim_comp_cs->add_ev(carID);
#ifdef PYTHON_MODULE
    py_cmd_p_ev_kW.push_back(0.0);
#endif
}

void ControlUnit::set_output_object(CUOutput* output_obj) {
    this->output_obj = output_obj;
}

void ControlUnit::set_exp_pv_params_A(float value) {
    // do something only if this unit is selected for PV expansion
    if (has_sim_pv) {
        if (Global::get_exp_pv_sizing_mode() == global::PVSizingMode::StaticPVSize) {
            delete sim_comp_pv;
            sim_comp_pv = new ComponentPV(value, locationID);
        } else {
            throw runtime_error("Error: ControlUnit::set_exp_pv_params_A() has been called even though PV static mode is not set!");
        }
    }
}

void ControlUnit::set_exp_pv_params_B(float kWp_per_m2, float min_kWp_sec, float max_kWp_sec, float max_kWp_unit) {
    // do something only if this unit is selected for PV expansion
    if (has_sim_pv) {
        if (Global::get_exp_pv_sizing_mode() == global::PVSizingMode::StaticPVSize) {
            throw runtime_error("Error: ControlUnit::set_exp_pv_params_B() has been called even though PV static mode is set!");
        } else {
            delete sim_comp_pv;
            sim_comp_pv = new ComponentPV(kWp_per_m2, min_kWp_sec, max_kWp_sec, max_kWp_unit, locationID);
        }
    }
}

void ControlUnit::set_exp_bs_maxE_kWh(float value) {
    if (has_sim_bs)
        sim_comp_bs->set_maxE_kWh(value);
}

void ControlUnit::set_exp_bs_maxP_kW(float value) {
    if (has_sim_bs)
        sim_comp_bs->set_maxP_kW(value);
}

void ControlUnit::set_exp_bs_E_P_ratio(float value) {
    if (has_sim_bs)
        sim_comp_bs->set_maxP_by_EPRatio(value);
}

void ControlUnit::change_control_horizon_in_ts(unsigned int new_horizon) {
    if (optimized_controller != NULL)
        optimized_controller->reset(new_horizon);
    if (has_sim_hp)
        sim_comp_hp->set_horizon_in_ts(new_horizon);
    if (has_sim_cs)
        sim_comp_cs->set_horizon_in_ts(new_horizon);
}

void ControlUnit::remove_sim_added_pv() {
    if (has_sim_pv) {
        has_sim_pv = false;
        delete sim_comp_pv;
        sim_comp_pv = NULL;
    } else {
        cerr << "Warning: A PV installation should be removed from one CU, but no one is sim. added!" << endl;
    }
}

void ControlUnit::remove_sim_added_bs() {
    if (has_sim_bs) {
        has_sim_bs = false;
        delete sim_comp_bs;
        sim_comp_bs = NULL;
    } else {
        cerr << "Warning: A battery storage should be removed from one CU, but no one is sim. added!" << endl;
    }
}

void ControlUnit::remove_sim_added_components() {
    if (has_sim_pv) {
        has_sim_pv = false;
        delete sim_comp_pv;
        sim_comp_pv = NULL;
    }
	if (has_sim_bs) {
        has_sim_bs = false;
        delete sim_comp_bs;
        sim_comp_bs = NULL;
    }
	if (has_sim_hp) {
        has_sim_hp = false;
        delete sim_comp_hp;
        sim_comp_hp = NULL;
    }
	if (has_sim_cs) {
        has_sim_cs = false;
        sim_comp_cs->disable_station();
    }
    is_sim_expanded = false;
}

void ControlUnit::reset_internal_state() {
    current_load_vSM_kW = 0.0;
    self_produced_load_kW = 0.0;
    current_total_consumption_kW = 0.0;
    current_PV_feedin_to_grid_kW      = 0.0;
    current_BS_feedin_to_grid_kW      = 0.0;
    current_CHP_feedin_to_grid_kW     = 0.0;
    current_wind_feedin_to_grid_kW    = 0.0;
    current_unknown_feedin_to_grid_kW = 0.0;
    //
#ifdef PYTHON_MODULE
    py_control_commands_obtained = false;
    py_cmd_p_bs_kW = 0.0;
    py_cmd_p_hp_kW = 0.0;
#endif
    //
    sum_of_consumption_kWh    = 0.0;
    sum_of_self_cons_kWh      = 0.0;
    sum_of_mu_cons_kWh        = 0.0;
    sum_of_feed_into_grid_kWh = 0.0;
    sum_of_grid_demand_kWh    = 0.0;
    sum_of_rem_pow_costs_EUR       = 0.0;
    sum_of_saved_pow_costs_EUR     = 0.0;
    sum_of_feedin_revenue_EUR      = 0.0;
    sum_of_emissions_cbgd_kg_CO2eq = 0.0;
    sum_of_emissions_avoi_kg_CO2eq = 0.0;
    sum_of_errors_in_cntrl         = 0;
    sum_of_errors_in_cntrl_cmd_appl= 0;
    peak_grid_demand_kW            = 0.0;
    //
    sum_of_cweek_consumption_kWh         = 0.0;
    sum_of_cweek_self_cons_kWh           = 0.0;
    sum_of_cweek_mu_cons_kWh             = 0.0;
    sum_of_cweek_feed_into_grid_kWh      = 0.0;
    sum_of_cweek_grid_demand_kWh         = 0.0;
    sum_of_cweek_rem_pow_costs_EUR       = 0.0;
    sum_of_cweek_saved_pow_costs_EUR     = 0.0;
    sum_of_cweek_feedin_revenue_EUR      = 0.0;
    sum_of_cweek_emissions_cbgd_kg_CO2eq = 0.0;
    sum_of_cweek_emissions_avoi_kg_CO2eq = 0.0;
    sum_of_cweek_errors_in_cntrl         = 0;
    sum_of_cweek_errors_in_cntrl_cmd_appl= 0;
    cweek_peak_grid_demand_kW            = 0.0;
    //
    if (has_sim_pv) {
        sim_comp_pv->resetInternalState();
    }
    if (has_sim_bs) {
        sim_comp_bs->resetInternalState();
    }
    if (has_sim_hp) {
        sim_comp_hp->resetInternalState();
    }
    sim_comp_cs->resetInternalState();
    //
    ts_since_last_opti_run = Global::get_control_update_freq_in_ts();
}


float ControlUnit::get_rsm_load_at_ts(unsigned long ts) const{
    if (ts <= 0 || ts > Global::get_n_timesteps())
        return 0.0;
    float total_load = 0.0;
    for (MeasurementUnit* mu : *connected_units) {
        total_load += mu->get_rsm_load_at_ts(ts);
    }
    return total_load;
}

float ControlUnit::get_rsm_demand_at_ts(unsigned long ts) const{
    if (ts <= 0 || ts > Global::get_n_timesteps())
        return 0.0;
    float total_load = get_rsm_load_at_ts(ts);
    if (total_load < 0.0)
        return 0.0;
    return total_load;
}

float ControlUnit::get_rsm_feedin_at_ts(unsigned long ts) const{
    if (ts <= 0 || ts > Global::get_n_timesteps())
        return 0.0;
    float total_load = get_rsm_load_at_ts(ts);
    if (total_load > 0.0)
        return 0.0;
    return -total_load;
}

bool ControlUnit::compute_next_value(unsigned long ts) {
    //
    // This function computes the next value
    // for this complete control unit.
    // It also calls the methods of all connected components
    // to calculate their actions for the next step.
    //

    // Part A: Update all components and set their internal variables to the new time step
    // Component PV
    if (has_sim_pv) {
        sim_comp_pv->calculateCurrentFeedin(ts);
    }
    // Component HP
    if (has_sim_hp) {
        sim_comp_hp->computeNextInternalState(ts);
    }
    // EV Charging station (propagates the new time step to the connected vehicles)
    if (sim_comp_cs->is_enabled()) {
        sim_comp_cs->setCarStatesForTimeStep(ts);
    }

    unsigned long n_cars = 0;
    if (has_sim_cs)
        n_cars = sim_comp_cs->get_n_EVs();
    
    bool error_happend_during_cntrl_cmd_appl = false;

    // Part B: Make decisions using an optimization if selected
    if (Global::get_controller_mode() == global::ControllerMode::OptimizedWithPerfectForecast &&
        ( has_sim_pv || has_sim_bs || has_sim_hp || has_sim_cs )

#ifdef PYTHON_MODULE
        // Jump optimization if control already given externally
        && !py_control_commands_obtained
#endif

       )
    {
        // Check number of time steps passed since last optimization run
        ts_since_last_opti_run++;
        if (ts_since_last_opti_run >= Global::get_control_update_freq_in_ts()) {
            ts_since_last_opti_run = 0;
            // run optimization in this case
            if (optimized_controller == NULL) {
#ifdef USE_OR_TOOLS
                optimized_controller = new ORToolsLPController(unitID, Global::get_control_horizon_in_ts(), n_cars);
#elif USE_GUROBI
                optimized_controller = new GurobiLPController(unitID, Global::get_control_horizon_in_ts(), n_cars);
#else
                throw std::runtime_error("No optimization backend selected during compile time!");
#endif
            }
            bool optimize_PV_size = false; // this value is set some lines below, if true required
            bool optimize_BS_size = has_sim_bs && Global::get_battery_capacity_computation_mode() == global::BatteryCapacityComputationMode::Optimized;
            //
            // 3. Get the not-shiftable and shiftable loads over the prediction horizon
            // - measurement units
            std::vector<float> future_resid_demand_kW(Global::get_control_horizon_in_ts());
            for (size_t fts = 0; fts < Global::get_control_horizon_in_ts(); fts++) {
                future_resid_demand_kW[fts] = 0.0;
                for (MeasurementUnit* mu : *connected_units)
                    future_resid_demand_kW[fts] += mu->get_rsm_demand_at_ts(ts + fts);
            }
            // - PV generation: sim_comp_pv->get_generation_at_ts_kW( ... future steps over horizon ... )
            const std::vector<double>* future_pv_generation_kW = st__empty_vector_for_time_horizon;
            if (has_sim_pv) {
                future_pv_generation_kW = sim_comp_pv->get_future_generation_kW();
            }
            // - -> get PV generation over the complete time span if Global::get_exp_pv_sizing_mode() == global::PVSizingMode::Optimized
            std::optional< std::list<std::vector<double>> > total_PV_generation_per_section_kW;
            std::optional< std::list<double>              > max_PV_power_per_section_kWp;
            if (has_sim_pv && Global::get_exp_pv_sizing_mode() == global::PVSizingMode::Optimized) {
                total_PV_generation_per_section_kW = sim_comp_pv->get_total_generation_by_section_kW();
                max_PV_power_per_section_kWp       = sim_comp_pv->get_kWp_per_section();
                optimize_PV_size = true;
            }
            // - heat pump
            const std::vector<double>* future_hp_shiftable_maxP = st__empty_vector_for_time_horizon;
            const std::vector<double>* future_hp_shiftable_minP = st__empty_vector_for_time_horizon;
            const std::vector<double>* future_hp_shiftable_maxE = st__empty_vector_for_time_horizon;
            const std::vector<double>* future_hp_shiftable_minE = st__empty_vector_for_time_horizon;
            //float max_p_hp_kW = 0.0;
            if (has_sim_hp) {
                future_hp_shiftable_maxP = sim_comp_hp->get_future_max_power_kW();
                future_hp_shiftable_minP = sim_comp_hp->get_future_min_power_kW();
                future_hp_shiftable_maxE = sim_comp_hp->get_future_max_consumption_kWh();
                future_hp_shiftable_minE = sim_comp_hp->get_future_min_consumption_kWh();
                //max_p_hp_kW = sim_comp_hp->get_rated_power_without_AUX();
            }
            // - charging station
            double max_p_cs_kW = 0.0;
            const std::vector<const std::vector<double>*>* future_ev_shiftable_maxE = NULL;
            const std::vector<const std::vector<double>*>* future_ev_shiftable_minE = NULL;
            const std::vector<const std::vector<double>*>* future_ev_maxP           = NULL;
            if (has_sim_cs) {
                future_ev_shiftable_maxE = sim_comp_cs->get_future_max_consumption_kWh();
                future_ev_shiftable_minE = sim_comp_cs->get_future_min_consumption_kWh();
                future_ev_maxP = sim_comp_cs->get_future_max_power_kW();
                max_p_cs_kW = sim_comp_cs->get_max_P_kW();
            }
            // - current battery SOC
            double current_bs_charge_kWh = 0.0;
            double max_e_bs_kWh          = 0.0;
            double max_p_bs_kW           = 0.0;
            if (has_sim_bs) {
                current_bs_charge_kWh = sim_comp_bs->get_currentCharge_kWh();
                max_e_bs_kWh          = sim_comp_bs->get_maxE_kWh();
                max_p_bs_kW           = sim_comp_bs->get_maxP_kW();
            }
            //
            // 4a. Check if number of optimization vars would not exceed the current limit
            unsigned long n_opti_vars = 0;
            if (Global::get_max_parallel_opti_vars() > 0) {
                // Compute the number of optimization variables for this unit
                n_opti_vars  = Global::get_control_horizon_in_ts() * 10 + 1;
                n_opti_vars += Global::get_control_horizon_in_ts() * n_cars * 2;
                if (Global::get_controller_bs_grid_charging_mode() == global::ControllerBSGridChargingMode::NoGridCharging) {
                    n_opti_vars += Global::get_control_horizon_in_ts() * 4;
                } else if (Global::get_controller_bs_grid_charging_mode() == global::ControllerBSGridChargingMode::OnlyGridCharging) {
                    n_opti_vars += Global::get_control_horizon_in_ts() * 6;
                }
                if (Global::get_controller_optimization_target() == global::ControllerOptimizationTarget::PeakLoad) {
                    n_opti_vars += 1;
                }
                // Check only if a limit is set
                {
                    std::unique_lock<std::mutex> lock(n_curr_optimization_vars_mutex);
                    // If n_curr_optimization_vars == 0, we execute the optimization regardless of n_opti_vars,
                    // because we would end up in a dead lock otherwise
                    if (n_curr_optimization_vars == 0) {
                        n_curr_optimization_vars += n_opti_vars;
                    } else {
                        // wait until it is possible to continue
                        curr_opti_cv.wait(lock, [&]() {
                            return n_curr_optimization_vars + n_opti_vars <= Global::get_max_parallel_opti_vars() ||
                                   n_curr_optimization_vars == 0;
                        } );
                        n_curr_optimization_vars += n_opti_vars;
                    }
                }
            }
            //
            // 4b. Run the optimization and get the results (and apply optimal PV/BS size if required)
            optimization_call_counter++;
            bool opti_ret_val = optimized_controller->updateController(
                ts, max_p_bs_kW, max_e_bs_kWh, max_p_cs_kW, current_bs_charge_kWh,
                future_resid_demand_kW,   *future_pv_generation_kW,
                *future_hp_shiftable_maxP, *future_hp_shiftable_minP,
                *future_hp_shiftable_maxE, *future_hp_shiftable_minE,
                 future_ev_shiftable_maxE,  future_ev_shiftable_minE,  future_ev_maxP,
                 optimize_PV_size, optimize_BS_size,
                 optimize_PV_size ? &*total_PV_generation_per_section_kW : NULL,
                 optimize_PV_size ? &*max_PV_power_per_section_kWp : NULL
                );
            // Reduce the number of variables for currently running optimizations again
            if (Global::get_max_parallel_opti_vars() > 0) {
                //std::unique_lock<std::mutex> lock(n_curr_optimization_vars_mutex);
                n_curr_optimization_vars -= n_opti_vars;
                //lock.unlock(); // locking / unlocking not required here!
                curr_opti_cv.notify_all();
            }
            // Create output with optimization details if selected
            if (Global::get_create_control_cmd_output() && generate_output_for_this_unit) {
                output::outputControlCommandDetails(
                    ts, unitID, opti_ret_val, max_p_bs_kW, max_e_bs_kWh, max_p_cs_kW, current_bs_charge_kWh,
                    future_resid_demand_kW,   *future_pv_generation_kW,
                    *future_hp_shiftable_maxP, *future_hp_shiftable_minP,
                    *future_hp_shiftable_maxE, *future_hp_shiftable_minE,
                    future_ev_shiftable_maxE,  future_ev_shiftable_minE,  future_ev_maxP,
                    optimized_controller->get_future_bs_power_kW(),
                    optimized_controller->get_future_hp_power_kW(),
                    optimized_controller->get_future_ev_power_kW()
                );
            }
            // Check optimization return value
            if (!opti_ret_val) {
                std::cerr << "Optimization error for unit with ID " << unitID << " at time step " << ts << ".\n";
                sum_of_errors_in_cntrl++;
                sum_of_cweek_errors_in_cntrl++;
                if (Global::get_stop_on_cc_err()) {
                    // return error ONLY if selected as a program option
                    return false;
                } else {
                    // Use value of the last optimization run if there is no new value and shift vector by one position.
                    // If the control commands received in this way have become invalid, the components themselves detect this and set the demand to a permissible value.
                    optimized_controller->shiftVectorsByOnePlace();
                }
            }
            // Process result of PV and/or BS sizing
            if (optimize_PV_size) {
                auto new_values = optimized_controller->get_optimal_PV_size_per_section_kW();
                sim_comp_pv->set_kWp_per_section( new_values );
            }
            if (optimize_BS_size) {
                sim_comp_bs->set_maxE_kWh( optimized_controller->get_optimal_BS_size_kWh() );
            }

        } else {
            // else (i.e., no opti executed in this position):
            // just move the cache by one
            optimized_controller->shiftVectorsByOnePlace();
        }

        //
        // 5. Propagate the results to the components
        if (has_sim_bs) {
            // TODO: Instead of this, implement minimum charge_from_grid_into_BS rule into the optimization model directly
            // TODO: New scheduling ideas from rule_based not implemented here yet!
            if (Global::get_surplus_controller_enabled()) {
                sim_comp_bs->set_chargeRequest( optimized_controller->get_future_bs_power_kW()[0] + surplus::SurplusController::GetChargeRequestForUnit(this->get_unitID()) );
            } else {
                sim_comp_bs->set_chargeRequest( optimized_controller->get_future_bs_power_kW()[0] );
            }
            
        }
        if (has_sim_hp) {
            bool ret_val = sim_comp_hp->setDemandToGivenValue( optimized_controller->get_future_hp_power_kW()[0] );
            if (!ret_val) error_happend_during_cntrl_cmd_appl = true;
        }
        if (has_sim_cs) {
            const std::vector<std::vector<double>>& ev_power_kW = optimized_controller->get_future_ev_power_kW();
            std::vector<double> ev_power_kW_this_ts;
            for (unsigned int evIdx = 0; evIdx < n_cars; evIdx++)
                ev_power_kW_this_ts.push_back( ev_power_kW[evIdx][0] );
            bool ret_val = sim_comp_cs->setDemandToGivenValues( ev_power_kW_this_ts );
            if (!ret_val) error_happend_during_cntrl_cmd_appl = true;
        }

    // In case of presence of control commands obtained from the python module: apply those
#ifdef PYTHON_MODULE
    } else if (py_control_commands_obtained) {
        if (has_sim_bs) {
            sim_comp_bs->set_chargeRequest( py_cmd_p_bs_kW );
        }
        if (has_sim_hp) {
            bool ret_val = sim_comp_hp->setDemandToGivenValue( py_cmd_p_hp_kW );
            if (!ret_val) error_happend_during_cntrl_cmd_appl = true;
        }
        if (has_sim_cs) {
            bool ret_val = sim_comp_cs->setDemandToGivenValues( py_cmd_p_ev_kW );
            if (!ret_val) error_happend_during_cntrl_cmd_appl = true;
        }
#endif

    } else {
        // If the rule-based strategy is selected, just use the current profile for the new demand
        // The rule-based strategy is also the fall back option
        if (has_sim_hp)
            sim_comp_hp->setDemandToProfileData(ts);
        if (has_sim_cs)
            sim_comp_cs->setDemandToProfileData(ts);
    }

    // Process information on errors during command application
    if (error_happend_during_cntrl_cmd_appl) {
        sum_of_errors_in_cntrl_cmd_appl++;
        sum_of_cweek_errors_in_cntrl_cmd_appl++;
    }


    // Part C: Compute the load at the virtual smart meter

    current_total_consumption_kW = 0.0;
    double current_load_all_rSMs_kW = 0.0;
    double load_pv = 0.0;
    double load_hp = 0.0;
    double load_cs = 0.0;
    double load_bs = 0.0;
    double bs_SOC  = 0.0;
    unsigned long n_cars_pc  = 0; // Number of cars parking at home, and connected
    unsigned long n_cars_pnc = 0; // Number of cars parking at home, but not connected

    //
    // 1. get sum of all real smart meter values
    for (MeasurementUnit* mu : *connected_units) {
        if (!mu->compute_next_value(ts))
            return false;
        current_load_all_rSMs_kW += mu->get_current_ts_rsm_value();
    }
    current_load_vSM_kW = current_load_all_rSMs_kW;
    current_total_consumption_kW = current_load_all_rSMs_kW;
    if (current_total_consumption_kW < 0) current_total_consumption_kW = 0;
    // // //float load_bevore_local_pv_bess = current_load_vSM_kW;
    //
    // 2. get PV feedin
    if (has_sim_pv) {
        load_pv = sim_comp_pv->get_currentGeneration_kW();
        current_load_vSM_kW -= load_pv;
    }
    //
    // 3. get the demand of the heat pump
    if (has_sim_hp) {
        load_hp = sim_comp_hp->get_currentDemand_kW();
        current_load_vSM_kW += load_hp;
        current_total_consumption_kW += load_hp;
    }
    //
    // 4. get the effect of the EV charging station
    if (sim_comp_cs->is_enabled()) {
        load_cs = sim_comp_cs->get_currentDemand_kW();
        current_load_vSM_kW += load_cs;
        if (load_cs > 0)
            current_total_consumption_kW += load_cs;
        // TODO in bidirectional charging case: If load_cs < 0 (ie the battery is discharged), we need to account that as well!
        // information for output
        n_cars_pc  = sim_comp_cs->get_n_EVs_pc();
        n_cars_pnc = sim_comp_cs->get_n_EVs_pnc();
    }
    // if load_evchst > 0: current_total_consumption_kW += load_evchst; // ... only problem: EV can potentially feed-in energy taken from somewhere else
    //
    // 5. send situation to battery storage and get its resulting action
    double charge_request_kW = 0.0;
    if (has_sim_bs) {
        if (
            Global::get_controller_mode() == global::ControllerMode::RuleBased
#ifdef PYTHON_MODULE
            && !py_control_commands_obtained
#endif
        ) {
            if (Global::get_surplus_controller_enabled()) {
                // Checking non-linear behavior where battery cannot discharge full requested power
                initial_charge_request_kW = sim_comp_bs->validateChargeRequest(-current_load_vSM_kW);
                charge_request_kW = initial_charge_request_kW + surplus::SurplusController::GetChargeRequestForUnit(this->get_unitID());
            } else {
                initial_charge_request_kW = -current_load_vSM_kW;
                charge_request_kW = -current_load_vSM_kW;
            }
            sim_comp_bs->set_chargeRequest( charge_request_kW );
        } // no else case here ... in the case of an optimized controller the action has been set before
        sim_comp_bs->calculateActions();
        load_bs = sim_comp_bs->get_currentLoad_kW();
        current_load_vSM_kW += load_bs;

        bs_SOC = sim_comp_bs->get_SOC();

    }

    //
    // compute self-produced load, that is directly consumed
    if (
#ifndef PYTHON_MODULE
         (
          Global::get_controller_mode() == global::ControllerMode::RuleBased ||
          ( Global::get_controller_mode() == global::ControllerMode::OptimizedWithPerfectForecast && Global::get_controller_bs_grid_charging_mode() == global::ControllerBSGridChargingMode::NoGridCharging )
         ) && !Global::get_surplus_controller_enabled()
#else
          false
#endif
    ) {
        // in this case, load_pv - load_bs >= 0 always holds --> calculation is easy as BS can ouly be charged from the PV, and thus stores always 100% PV generated energy
        self_produced_load_kW = std::min(current_total_consumption_kW, load_pv - load_bs);
    } else {
        // separate two cases: BS charging or discharging
        if (load_bs > 0) {
            // - case: BS charging
            self_produced_load_kW = std::min( load_pv, current_total_consumption_kW);
            // analyze the source of the charged energy and send this information to the BS
            if (has_sim_bs) {
                double bs_charging_from_pv = std::min( load_bs, load_pv - self_produced_load_kW);
                sim_comp_bs->set_grid_charged_amount(load_bs - bs_charging_from_pv); // send this information to the battery storage
            }
        } else {
            // - case: BS discharging
            // get amount of load_bs that has been taken from PV feed-in
            double bs_generation_from_pv = 0.0;
            if (has_sim_bs && load_bs < 0.0) {
                bs_generation_from_pv = - ( sim_comp_bs->get_gridOnly_discharge_kW() + load_bs ); // Mind: load_bs < 0, but sim_comp_bs->get_gridOnly_discharge_kW() >= 0
            }
            // calculate self-produced amount of locally consumed energy
            self_produced_load_kW = std::min(load_pv + bs_generation_from_pv, current_total_consumption_kW);
        }
    }

    // Compute current energy flows out of the power flows
    double grid_feedin_kWh = 0.0;
    if (current_load_vSM_kW < 0) {
        grid_feedin_kWh = -current_load_vSM_kW * Global::get_time_step_size_in_h();
    }
    double grid_demand_kWh = 0.0;
    if (current_load_vSM_kW > 0) {
        grid_demand_kWh =  current_load_vSM_kW  * Global::get_time_step_size_in_h();
    }
    double self_cons_kWh       = Global::get_time_step_size_in_h() * self_produced_load_kW;
    double cons_kWh            = Global::get_time_step_size_in_h() * current_total_consumption_kW;
    // add values to summation variables
    sum_of_consumption_kWh       += cons_kWh;
    sum_of_cweek_consumption_kWh += cons_kWh;
    sum_of_self_cons_kWh         += self_cons_kWh;
    sum_of_cweek_self_cons_kWh   += self_cons_kWh;
    if (current_load_all_rSMs_kW > 0) {
        double mu_cons = Global::get_time_step_size_in_h() * current_load_all_rSMs_kW;
        sum_of_mu_cons_kWh       += mu_cons;
        sum_of_cweek_mu_cons_kWh += mu_cons;
    }
    sum_of_feed_into_grid_kWh += grid_feedin_kWh;
    sum_of_grid_demand_kWh    += grid_demand_kWh;
    sum_of_cweek_feed_into_grid_kWh += grid_feedin_kWh;
    sum_of_cweek_grid_demand_kWh    += grid_demand_kWh;
    //
    if (peak_grid_demand_kW < current_load_vSM_kW)
        peak_grid_demand_kW = current_load_vSM_kW;
    // add values to summation variables where potentially the global time series have to be conducted
    //  a) price computation
    double p, p1, p2;
    p = grid_feedin_kWh * Global::get_feed_in_tariff();
    sum_of_feedin_revenue_EUR       += p;
    sum_of_cweek_feedin_revenue_EUR += p;
    if (global::eprices_local_ts != NULL) {
        p1 = grid_demand_kWh * global::eprices_local_ts[ts - 1];
        p2 = self_cons_kWh   * global::eprices_local_ts[ts - 1];
    } else {
        p1 = grid_demand_kWh * Global::get_demand_tariff();
        p2 = self_cons_kWh   * Global::get_demand_tariff();
    }
    sum_of_rem_pow_costs_EUR         += p1;
    sum_of_saved_pow_costs_EUR       += p2;
    sum_of_cweek_rem_pow_costs_EUR   += p1;
    sum_of_cweek_saved_pow_costs_EUR += p2;
    //  b) emission computations
    double r1, r2;
    if (global::emission_ts != NULL) {
        r1 = grid_demand_kWh * global::emission_ts[ts - 1] / 1000;
        r2 = self_cons_kWh   * global::emission_ts[ts - 1] / 1000;
    } else {
        r1 = grid_demand_kWh * Global::get_emissions_g_CO2eq_per_kWh() / 1000;
        r2 = self_cons_kWh   * Global::get_emissions_g_CO2eq_per_kWh() / 1000;
    }
    sum_of_emissions_cbgd_kg_CO2eq       += r1;
    sum_of_emissions_avoi_kg_CO2eq       += r2;
    sum_of_cweek_emissions_cbgd_kg_CO2eq += r1;
    sum_of_cweek_emissions_avoi_kg_CO2eq += r2;
    //
    if (cweek_peak_grid_demand_kW < current_load_vSM_kW)
        cweek_peak_grid_demand_kW = current_load_vSM_kW;

    // Part C.2: In case of a grid feed-in, determine the source of the feed-in
    current_PV_feedin_to_grid_kW      = 0.0;
    current_BS_feedin_to_grid_kW      = 0.0;
    current_CHP_feedin_to_grid_kW     = 0.0;
    current_wind_feedin_to_grid_kW    = 0.0;
    current_unknown_feedin_to_grid_kW = 0.0;
    if (current_load_vSM_kW < 0.0) {
        // Attribute the feedin iteratively to the components, until the remaining power is zero
        double remaining_feedin_to_attribute_kW = -current_load_vSM_kW;
        //
        if (has_sim_pv && sim_comp_pv->get_currentGeneration_kW() > 0.0) {
            current_PV_feedin_to_grid_kW = sim_comp_pv->get_currentGeneration_kW();
            if (current_PV_feedin_to_grid_kW > remaining_feedin_to_attribute_kW)
                current_PV_feedin_to_grid_kW = remaining_feedin_to_attribute_kW;
            remaining_feedin_to_attribute_kW -= current_PV_feedin_to_grid_kW;
        }
        if (remaining_feedin_to_attribute_kW > 0.0 &&
            !has_sim_pv && has_pv() )
        {
            // No simulated PV, but still a PV existing
            // --> Attribute completly to the PV installation
            current_PV_feedin_to_grid_kW = remaining_feedin_to_attribute_kW;
            remaining_feedin_to_attribute_kW = 0.0;
        }
        //
        if (remaining_feedin_to_attribute_kW > 0.0 &&
            has_sim_bs && sim_comp_bs->get_currentLoad_kW() < 0.0)
        {
            current_BS_feedin_to_grid_kW = sim_comp_bs->get_currentLoad_kW();
            if (current_BS_feedin_to_grid_kW > remaining_feedin_to_attribute_kW)
                current_BS_feedin_to_grid_kW = remaining_feedin_to_attribute_kW;
            remaining_feedin_to_attribute_kW -= current_BS_feedin_to_grid_kW;
        }
        //
        if (remaining_feedin_to_attribute_kW > 0.0 &&
            has_chp())
        {
            // If CHP is present in the input data -> attribute all remaining feed-in to the CHP
            current_CHP_feedin_to_grid_kW = remaining_feedin_to_attribute_kW;
            remaining_feedin_to_attribute_kW = 0.0;
        }
        //
        if (remaining_feedin_to_attribute_kW > 0.0 &&
            has_wind())
        {
            // If wind turbine is present in the input data -> attribute all remaining feed-in to the wind turbine
            current_wind_feedin_to_grid_kW = remaining_feedin_to_attribute_kW;
            remaining_feedin_to_attribute_kW = 0.0;
        }
        // Attribute the final feedin to unknown
        if (remaining_feedin_to_attribute_kW > 0.0)
        {
            current_unknown_feedin_to_grid_kW = remaining_feedin_to_attribute_kW;
            // obsolet: remaining_feedin_to_attribute_kW = 0.0;
        }
    }

    //
    // output current status
    if (output_obj != NULL && generate_output_for_this_unit)
        output_obj->output_for_one_cu(
            unitID, ts,
            current_load_vSM_kW,
            current_load_all_rSMs_kW,
            self_produced_load_kW,
            load_pv, bs_SOC, load_bs,
            load_hp, load_cs,
            n_cars_pc, n_cars_pnc);

#ifdef PYTHON_MODULE
    py_control_commands_obtained = false;
#endif

    return true;
}

#ifdef PYTHON_MODULE
void ControlUnit::send_control_commands_from_py_interface(double p_bs_kW, double p_hp_kW, const std::vector<double>& p_ev_kW) {
    py_control_commands_obtained = true;
    py_cmd_p_bs_kW = p_bs_kW;
    py_cmd_p_hp_kW = p_hp_kW;
    for (unsigned long evIdx = 0; evIdx < py_cmd_p_ev_kW.size() && evIdx < p_ev_kW.size(); evIdx++) {
        py_cmd_p_ev_kW[evIdx] = p_ev_kW[evIdx];
    }
}
#endif

void ControlUnit::InitializeStaticVariables(unsigned long n_CUs) {
    // reserve enough space in vector of CUs
    st__cu_list.reserve(n_CUs);
    // define static empty vector for the optimization
    if (st__empty_vector_for_time_horizon == NULL) {
        st__empty_vector_for_time_horizon = new std::vector<double>;
        st__empty_vector_for_time_horizon->clear();
        st__empty_vector_for_time_horizon->resize( Global::get_control_horizon_in_ts(), 0.0 );
    }
    // Initialize static members
    if (Global::get_controller_mode() != global::ControllerMode::RuleBased) {
#ifdef USE_OR_TOOLS
#elif USE_GUROBI
        GurobiLPController::InitializeGurobiEnvironment();
#endif
    }
}

void ControlUnit::VacuumInstancesAndStaticVariables() {
    // delete all instances
    for (unsigned long i = 0; i < st__n_CUs; i++) {
        delete st__cu_list[i];
        st__cu_list[i] = NULL;
    }
    st__n_CUs = 0;
    // static variables
    if (st__empty_vector_for_time_horizon != NULL) {
        delete st__empty_vector_for_time_horizon;
        st__empty_vector_for_time_horizon = NULL;
    }
    public_to_internal_id.clear();
    if (Global::get_controller_mode() != global::ControllerMode::RuleBased) {
#ifdef USE_OR_TOOLS
#elif USE_GUROBI
        GurobiLPController::VaccumAllStaticVariables();
#endif
    }
}

ControlUnit* ControlUnit::GetInstancePublicID(unsigned long public_unitID) {
    try {
        size_t internal_id = public_to_internal_id.at(public_unitID); // throws std::out_of_range if public_id is unknown
		return st__cu_list[internal_id];
    } catch (std::out_of_range const&) {
		return NULL;
    }
}

ControlUnit* ControlUnit::GetInstancePublicIDWE(unsigned long public_unitID) {
    try {
        size_t internal_id = public_to_internal_id.at(public_unitID); // throws std::out_of_range if public_id is unknown
		return st__cu_list[internal_id];
    } catch (std::out_of_range const&) {
		throw runtime_error("ControlUnitID " + std::to_string(public_unitID) + " is out of range!");
    }
}

//ControlUnit* ControlUnit::GetInstanceAtLocationID(unsigned long locationID) {
//    return location_to_cu_map[locationID];
//}

void ControlUnit::PreprocessEVData() {
    for (ControlUnit* cu : st__cu_list) {
        cu->preprocess_ev_data();
    }
}

void ControlUnit::ResetAllInternalStates() {
    //
    // Call this function to reset all internal states of
    // the unit and all connected components (as long as they
    // have interal states, like a battery has the current SOC
    // as internal state).
    //
    for (unsigned long i = 0; i < st__n_CUs; i++) {
        st__cu_list[i]->reset_internal_state();
    }
}

ControlUnitState ControlUnit::save_internal_state() const {
    ControlUnitState state;
    state.unitID = unitID;
    
    // Save component states
    if (has_sim_pv) {
        state.PVState = sim_comp_pv->saveInternalState();
    }
    if (has_sim_bs) {
        state.BSState = sim_comp_bs->saveInternalState();
    }
    if (has_sim_hp) {
        state.HPState = sim_comp_hp->saveInternalState();
    }
    if (has_sim_cs) {
        state.CSState = sim_comp_cs->saveInternalState();
    }
    
    // Save optimization state
    if (optimized_controller != nullptr) {
        state.optimizationState = optimized_controller->saveInternalState();
    }
    
    // Save summation variables from beginning of simulation run
    state.sum_of_consumption_kWh = sum_of_consumption_kWh;
    state.sum_of_self_cons_kWh = sum_of_self_cons_kWh;
    state.sum_of_mu_cons_kWh = sum_of_mu_cons_kWh;
    state.sum_of_feed_into_grid_kWh = sum_of_feed_into_grid_kWh;
    state.sum_of_grid_demand_kWh = sum_of_grid_demand_kWh;
    state.sum_of_rem_pow_costs_EUR = sum_of_rem_pow_costs_EUR;
    state.sum_of_saved_pow_costs_EUR = sum_of_saved_pow_costs_EUR;
    state.sum_of_feedin_revenue_EUR = sum_of_feedin_revenue_EUR;
    state.sum_of_emissions_cbgd_kg_CO2eq = sum_of_emissions_cbgd_kg_CO2eq;
    state.sum_of_emissions_avoi_kg_CO2eq = sum_of_emissions_avoi_kg_CO2eq;
    state.sum_of_errors_in_cntrl = sum_of_errors_in_cntrl;
    state.sum_of_errors_in_cntrl_cmd_appl = sum_of_errors_in_cntrl_cmd_appl;
    state.peak_grid_demand_kW = peak_grid_demand_kW;
    
    // Save weekly summation variables
    state.sum_of_cweek_consumption_kWh = sum_of_cweek_consumption_kWh;
    state.sum_of_cweek_self_cons_kWh = sum_of_cweek_self_cons_kWh;
    state.sum_of_cweek_mu_cons_kWh = sum_of_cweek_mu_cons_kWh;
    state.sum_of_cweek_feed_into_grid_kWh = sum_of_cweek_feed_into_grid_kWh;
    state.sum_of_cweek_grid_demand_kWh = sum_of_cweek_grid_demand_kWh;
    state.sum_of_cweek_rem_pow_costs_EUR = sum_of_cweek_rem_pow_costs_EUR;
    state.sum_of_cweek_saved_pow_costs_EUR = sum_of_cweek_saved_pow_costs_EUR;
    state.sum_of_cweek_feedin_revenue_EUR = sum_of_cweek_feedin_revenue_EUR;
    state.sum_of_cweek_emissions_cbgd_kg_CO2eq = sum_of_cweek_emissions_cbgd_kg_CO2eq;
    state.sum_of_cweek_emissions_avoi_kg_CO2eq = sum_of_cweek_emissions_avoi_kg_CO2eq;
    state.sum_of_cweek_errors_in_cntrl = sum_of_cweek_errors_in_cntrl;
    state.sum_of_cweek_errors_in_cntrl_cmd_appl = sum_of_cweek_errors_in_cntrl_cmd_appl;
    state.cweek_peak_grid_demand_kW = cweek_peak_grid_demand_kW;
    
    // Save optimization cache
    state.ts_since_last_opti_run = ts_since_last_opti_run;
    
    // Save current state variables
    state.current_load_vSM_kW = current_load_vSM_kW;
    state.self_produced_load_kW = self_produced_load_kW;
    state.current_total_consumption_kW = current_total_consumption_kW;
    state.current_PV_feedin_to_grid_kW = current_PV_feedin_to_grid_kW;
    state.current_BS_feedin_to_grid_kW = current_BS_feedin_to_grid_kW;
    state.current_CHP_feedin_to_grid_kW = current_CHP_feedin_to_grid_kW;
    state.current_wind_feedin_to_grid_kW = current_wind_feedin_to_grid_kW;
    state.current_unknown_feedin_to_grid_kW = current_unknown_feedin_to_grid_kW;
    
    return state;
}

void ControlUnit::restore_internal_state(const ControlUnitState& state) {
    // Validate unitID
    if (state.unitID != unitID) {
        throw std::runtime_error("ControlUnit::restore_internal_state: unitID mismatch. Expected " + 
                                 std::to_string(unitID) + " but got " + std::to_string(state.unitID));
    }
    
    // Restore component states
    if (has_sim_pv) {
        if (!std::holds_alternative<PVComponentState>(state.PVState)) {
            throw std::runtime_error("ControlUnit::restore_internal_state: PV component exists but state does not contain PVComponentState for unitID " + 
                                     std::to_string(unitID));
        }
        sim_comp_pv->restoreInternalState(state.PVState);
    }
    if (has_sim_bs) {
        if (!std::holds_alternative<BSComponentState>(state.BSState)) {
            throw std::runtime_error("ControlUnit::restore_internal_state: BS component exists but state does not contain BSComponentState for unitID " + 
                                     std::to_string(unitID));
        }
        sim_comp_bs->restoreInternalState(state.BSState);
    }
    if (has_sim_hp) {
        if (!std::holds_alternative<HPComponentState>(state.HPState)) {
            throw std::runtime_error("ControlUnit::restore_internal_state: HP component exists but state does not contain HPComponentState for unitID " + 
                                     std::to_string(unitID));
        }
        sim_comp_hp->restoreInternalState(state.HPState);
    }
    if (has_sim_cs) {
        if (!std::holds_alternative<CSComponentState>(state.CSState)) {
            throw std::runtime_error("ControlUnit::restore_internal_state: CS component exists but state does not contain CSComponentState for unitID " + 
                                     std::to_string(unitID));
        }
        sim_comp_cs->restoreInternalState(state.CSState);
    }
    
    // Restore optimization state
    if (optimized_controller != nullptr) {
        optimized_controller->restoreInternalState(state.optimizationState);
    }
    
    // Restore summation variables from beginning of simulation run
    sum_of_consumption_kWh = state.sum_of_consumption_kWh;
    sum_of_self_cons_kWh = state.sum_of_self_cons_kWh;
    sum_of_mu_cons_kWh = state.sum_of_mu_cons_kWh;
    sum_of_feed_into_grid_kWh = state.sum_of_feed_into_grid_kWh;
    sum_of_grid_demand_kWh = state.sum_of_grid_demand_kWh;
    sum_of_rem_pow_costs_EUR = state.sum_of_rem_pow_costs_EUR;
    sum_of_saved_pow_costs_EUR = state.sum_of_saved_pow_costs_EUR;
    sum_of_feedin_revenue_EUR = state.sum_of_feedin_revenue_EUR;
    sum_of_emissions_cbgd_kg_CO2eq = state.sum_of_emissions_cbgd_kg_CO2eq;
    sum_of_emissions_avoi_kg_CO2eq = state.sum_of_emissions_avoi_kg_CO2eq;
    sum_of_errors_in_cntrl = state.sum_of_errors_in_cntrl;
    sum_of_errors_in_cntrl_cmd_appl = state.sum_of_errors_in_cntrl_cmd_appl;
    peak_grid_demand_kW = state.peak_grid_demand_kW;
    
    // Restore weekly summation variables
    sum_of_cweek_consumption_kWh = state.sum_of_cweek_consumption_kWh;
    sum_of_cweek_self_cons_kWh = state.sum_of_cweek_self_cons_kWh;
    sum_of_cweek_mu_cons_kWh = state.sum_of_cweek_mu_cons_kWh;
    sum_of_cweek_feed_into_grid_kWh = state.sum_of_cweek_feed_into_grid_kWh;
    sum_of_cweek_grid_demand_kWh = state.sum_of_cweek_grid_demand_kWh;
    sum_of_cweek_rem_pow_costs_EUR = state.sum_of_cweek_rem_pow_costs_EUR;
    sum_of_cweek_saved_pow_costs_EUR = state.sum_of_cweek_saved_pow_costs_EUR;
    sum_of_cweek_feedin_revenue_EUR = state.sum_of_cweek_feedin_revenue_EUR;
    sum_of_cweek_emissions_cbgd_kg_CO2eq = state.sum_of_cweek_emissions_cbgd_kg_CO2eq;
    sum_of_cweek_emissions_avoi_kg_CO2eq = state.sum_of_cweek_emissions_avoi_kg_CO2eq;
    sum_of_cweek_errors_in_cntrl = state.sum_of_cweek_errors_in_cntrl;
    sum_of_cweek_errors_in_cntrl_cmd_appl = state.sum_of_cweek_errors_in_cntrl_cmd_appl;
    cweek_peak_grid_demand_kW = state.cweek_peak_grid_demand_kW;
    
    // Restore optimization cache
    ts_since_last_opti_run = state.ts_since_last_opti_run;
    
    // Restore current state variables
    current_load_vSM_kW = state.current_load_vSM_kW;
    self_produced_load_kW = state.self_produced_load_kW;
    current_total_consumption_kW = state.current_total_consumption_kW;
    current_PV_feedin_to_grid_kW = state.current_PV_feedin_to_grid_kW;
    current_BS_feedin_to_grid_kW = state.current_BS_feedin_to_grid_kW;
    current_CHP_feedin_to_grid_kW = state.current_CHP_feedin_to_grid_kW;
    current_wind_feedin_to_grid_kW = state.current_wind_feedin_to_grid_kW;
    current_unknown_feedin_to_grid_kW = state.current_unknown_feedin_to_grid_kW;
}

std::map<unsigned long, ControlUnitState> ControlUnit::SaveAllInternalStates() {
    std::map<unsigned long, ControlUnitState> states;
    
    for (unsigned long i = 0; i < st__n_CUs; i++) {
        ControlUnit* cu = st__cu_list[i];
        ControlUnitState state = cu->save_internal_state();
        states[cu->get_unitID()] = state;
    }
    
    return states;
}

void ControlUnit::RestoreAllInternalStates(const std::map<unsigned long, ControlUnitState> &states) {
    for (unsigned long i = 0; i < st__n_CUs; i++) {
        ControlUnit* cu = st__cu_list[i];
        unsigned long unitID = cu->get_unitID();
        auto it = states.find(unitID);
        if (it == states.end()) {
            throw std::runtime_error("ControlUnit::RestoreAllInternalStates: No state found for unitID " + 
                                     std::to_string(unitID));
        }
        cu->restore_internal_state(it->second);
    }
}

void ControlUnit::RemoveAllSimAddedComponents() {
    for (unsigned long i = 0; i < st__n_CUs; i++) {
        st__cu_list[i]->remove_sim_added_components();
    }
}

void ControlUnit::ChangeControlHorizonInTS(unsigned int new_horizon) {
    // Update global default variables
    if (st__empty_vector_for_time_horizon != NULL)
        delete st__empty_vector_for_time_horizon;
    st__empty_vector_for_time_horizon = new std::vector<double>;
    st__empty_vector_for_time_horizon->clear();
    st__empty_vector_for_time_horizon->resize( new_horizon, 0.0 );
    // Call method for every register object
    for (ControlUnit* cu : st__cu_list) {
        cu->change_control_horizon_in_ts(new_horizon);
    }
}

size_t ControlUnit::GetNumberOfCUsWithSimCompPV() {
    size_t counter = 0;
    for (unsigned long i = 0; i < st__n_CUs; i++) {
        if (st__cu_list[i]->has_sim_pv)
            counter++;
    }
    return counter;
}

size_t ControlUnit::GetNumberOfCUsWithSimCompHP() {
    size_t counter = 0;
    for (unsigned long i = 0; i < st__n_CUs; i++) {
        if (st__cu_list[i]->has_sim_hp)
            counter++;
    }
    return counter;
}

size_t ControlUnit::GetNumberOfCUsWithSimCompEV() {
    size_t counter = 0;
    for (unsigned long i = 0; i < st__n_CUs; i++) {
        if (st__cu_list[i]->has_sim_cs)
            counter++;
    }
    return counter;
}

double ControlUnit::GetAllSimCompBatteriesCharge_kWh() {
    double sum_kWh = 0.0;
    for (unsigned long i = 0; i < st__n_CUs; i++) {
        if (st__cu_list[i]->has_sim_bs)
            sum_kWh += st__cu_list[i]->sim_comp_bs->get_currentCharge_kWh();
    }
    return sum_kWh;
}

double ControlUnit::GetAllSimCompBatteriesCapacity_kWh() {
    double sum_kWh = 0.0;
    for (unsigned long i = 0; i < st__n_CUs; i++) {
        if (st__cu_list[i]->has_sim_bs)
            sum_kWh += st__cu_list[i]->sim_comp_bs->get_maxE_kWh();
    }
    return sum_kWh;
}





// ----------------------------- //
//      Implementation of        //
//       MeasurementUnit         //
// ----------------------------- //

size_t MeasurementUnit::st__n_MUs = 0;
std::vector<MeasurementUnit*> MeasurementUnit::st__mu_list;
std::map<unsigned long, unsigned long> MeasurementUnit::public_to_internal_id;

bool MeasurementUnit::InstantiateNewMeasurementUnit(size_t meUID, size_t public_unitID, std::string * meterPointName, size_t locID, 
                bool has_demand, bool has_feedin, bool has_pv_resid, bool has_pv_opens,
                bool has_bess,   bool has_hp,     bool has_wind,     bool has_evchst,
                bool has_chp, const std::string& data_source_path)
{
    // check, if public_id is known
    if (public_to_internal_id.contains(meUID))
        return false;
    // define new internal id
    unsigned long new_internal_id = st__n_MUs;
    // increment the list
    st__n_MUs++;
    // register public id
    public_to_internal_id.insert(std::pair<unsigned long, unsigned long>(meUID, new_internal_id));
    // create instance
    MeasurementUnit* new_obj = new MeasurementUnit(new_internal_id, meUID, public_unitID, meterPointName, locID, has_demand, has_feedin, has_pv_resid, has_pv_opens, has_bess, has_hp, has_wind, has_evchst, has_chp, data_source_path);
    // add to list of instances
    st__mu_list.push_back(new_obj);
    //
    return true;
}

MeasurementUnit::MeasurementUnit(size_t internalID, size_t meUID, size_t public_unitID, string * meterPointName, size_t locID,
                                 bool has_demand, bool has_feedin, bool has_pv_resid, bool has_pv_opens,
                                 bool has_bess,   bool has_hp,     bool has_wind,     bool has_evchst,
                                 bool has_chp, const std::string& data_source_path) :
    internal_id(internalID),
    meUID(meUID),
    higher_level_cu(ControlUnit::GetInstancePublicID(public_unitID)),
    meterPointName(meterPointName), locationID(locID),
    data_source_path(data_source_path)
{
    //
    // initialize instance variables
    current_load_rsm_kW = 0;
    rsm_has_demand = has_demand;
    rsm_has_feedin = has_feedin;
    rsm_with_pv_residential = has_pv_resid;
    rsm_with_pv_open_space  = has_pv_opens;
    rsm_with_bess  = has_bess;
    rsm_with_hp    = has_hp;
    rsm_with_evchst= has_evchst;
    rsm_with_wind  = has_wind;
    rsm_with_chp   = has_chp;
    data_loaded       = false;
    data_timestepID   = NULL;
    data_value_demand = NULL;
    data_value_feedin = NULL;
    //data_status_demand=NULL;
    //data_status_feedin=NULL;
    expansion_combination = expansion::genExpCombiAsBitRepr(has_pv_resid||has_pv_opens, has_bess, has_hp, has_evchst);

    if (higher_level_cu == NULL) {
        std::cerr << "Control unit with ID " << public_unitID << " not found! Cannot create measurement unit with meUID " << meUID << "." << std::endl;
        throw runtime_error("ControlUnitID " + std::to_string(public_unitID) + " is out of range!"); 
    }

    //
    // add this measurement unit to the list of
    // connected units in the control unit
    higher_level_cu->add_unit(this);
}

MeasurementUnit::~MeasurementUnit() {
    delete[] data_timestepID;
    delete[] data_value_demand;
    delete[] data_value_feedin;
    //delete[] data_status_demand;
    //delete[] data_status_feedin;

    delete meterPointName;

    data_timestepID = NULL;
    data_value_demand = NULL;
    data_value_feedin = NULL;
    //data_status_demand = NULL;
    //data_status_feedin = NULL;
}

bool MeasurementUnit::load_data() {
    //
    // initialize the arrays
    //
    data_timestepID    = new int[Global::get_n_timesteps()];    //int32    -> 4 bytes
    data_value_demand  = new float[Global::get_n_timesteps()];  //float32  -> 4 bytes
    data_value_feedin  = new float[Global::get_n_timesteps()];  //float32  -> 4 bytes
    //data_status_demand = new char[Global::get_n_timesteps()]; //int8     -> 1 byte
    //data_status_feedin = new char[Global::get_n_timesteps()]; //int8     -> 1 bytes

    //
    // open and parse the csv file
    //
    ifstream smeter_input;
    smeter_input.open(data_source_path);
    if (!smeter_input.good()) {
        cout << "Error when connecting to the smart meter file with path " << data_source_path << endl;
        return false;
    } else {
        string currLineString;
        // check first line (i.e. the header): do the columns equal what they should be?
        getline( smeter_input, currLineString );
        if ( std::count(currLineString.begin(), currLineString.end(), ',') != 4 ) {
            #ifdef DEBUG
            cerr << "DEBUG MESSAGE: std::count(currLineString.begin(), currLineString.end(), ',') = " << std::count(currLineString.begin(), currLineString.end(), ',') << endl;
            #endif
            cerr << "One smart meter data file header has not exactly 5 columns!" << endl;
            cerr << "The data file must have the following columns (mind the order):" << endl;
            cerr << "TimestepID,PowerDemand_kW,PowerDemandStatus,PowerFeedin_kW,PowerFeedinStatus" << endl;
            return false; // destructor of smeter_input closes it by default
        }
        stringstream headerSStream(currLineString);
        string headerSplitted[5];
        for (uint col = 0; col < 5; col++) {
            std::getline( headerSStream, headerSplitted[col], ',' );
        }
        if (
            headerSplitted[0] != "TimestepID"        ||
            headerSplitted[1] != "PowerDemand_kW"    ||
            headerSplitted[2] != "PowerDemandStatus" ||
            headerSplitted[3] != "PowerFeedin_kW"    ||
            headerSplitted[4] != "PowerFeedinStatus"
        ) {
            cerr << "There is one smart meter data file with a wrong header!" << endl;
            cerr << "The data file must have the following columns (mind the order):" << endl;
            cerr << "TimestepID,PowerDemand_kW,PowerDemandStatus,PowerFeedin_kW,PowerFeedinStatus" << endl;
            return false; // destructor of smeter_input closes it by default
        }
        //
        // main loop for collecting all data
        for (size_t r = 0; r < Global::get_n_timesteps(); r++) {
            // iterate over every row
            getline( smeter_input, currLineString );
            stringstream currLineStream( currLineString );
            string currLineSplitted[5];
            // Check if current line has more or less than 5 times an ',' -> throw an error!
            if ( std::count(currLineString.begin(), currLineString.end(), ',') != 4 ) {
                cerr << "One smart meter data file holds a row with not exactly 5 values!" << endl;
                return false; // destructor of smeter_input closes it by default
            }
            // read the values
            for (int col = 0; col < 5; col++) {
                // split this row on the ","
                getline( currLineStream, currLineSplitted[col], ',' );
            }
            // convert individual strings to int / float / char
            data_timestepID[   r ]  = stoi( currLineSplitted[0] );
            data_value_demand[ r ]  = stof( currLineSplitted[1] );
            data_value_feedin[ r ]  = stof( currLineSplitted[3] );
        }
        smeter_input.close();
    }
    
    #ifdef DEBUG_EXTRA_OUTPUT
    for (int i = 0; i < 10; i++) {
        printf("%6i %6.4f %6.4f \n",
            data_timestepID[i],
            data_value_demand[i],
            data_value_feedin[i]);
    }
    cout << endl;
    #endif

    data_loaded = true;
    return true;

}

double MeasurementUnit::get_total_demand_kWh() const {
    double cumsum = 0.0;
    for (unsigned long ts = 1; ts <= Global::get_n_timesteps(); ts++) {
        cumsum += data_value_demand[ts-1];
    }
    return cumsum;
}

bool MeasurementUnit::compute_next_value(unsigned long ts) {
    if (ts <= 0 || ts > Global::get_n_timesteps()) {
        current_load_rsm_kW = 0.0;
        return false;
    }
    unsigned long tsID = ts - 1;
    current_load_rsm_kW = data_value_demand[tsID] - data_value_feedin[tsID];
    return true;
}

float MeasurementUnit::get_rsm_load_at_ts(unsigned long ts) const {
    if (ts <= 0 || ts > Global::get_n_timesteps())
        return 0.0;
    return data_value_demand[ts - 1] - data_value_feedin[ts - 1];
}

float MeasurementUnit::get_rsm_demand_at_ts(unsigned long ts) const {
    if (ts <= 0 || ts > Global::get_n_timesteps())
        return 0.0;
    return data_value_demand[ts - 1];
}

float MeasurementUnit::get_rsm_feedin_at_ts(unsigned long ts) const {
    if (ts <= 0 || ts > Global::get_n_timesteps())
        return 0.0;
    return data_value_feedin[ts - 1];
}

void MeasurementUnit::InitializeStaticVariables(unsigned long n_MUs) {
    st__mu_list.reserve(n_MUs);
}

void MeasurementUnit::VacuumInstancesAndStaticVariables() {
    // delete all instances
    for (size_t i = 0; i < st__n_MUs; i++) {
        delete st__mu_list[i];
        st__mu_list[i] = NULL;
    }
    st__n_MUs = 0;
}

bool MeasurementUnit::LoadDataForAllInstances() {
    if (Global::get_n_threads() >= 2) {
        size_t n_threads = Global::get_n_threads();
        // create a vector of lists where the objects are stored that should be handeled by one thread
        std::vector<std::list<MeasurementUnit*>> object_allocation_list(n_threads);
        size_t counter  = 0;
        for (MeasurementUnit* mu : st__mu_list) {
            object_allocation_list[counter % n_threads].push_back(mu);
            counter += 1;
        }
        // create an atomic flag to indicate if an error has occured
        std::atomic<bool> error_occured(false);
        // start the threads
        std::vector<std::thread> thread_list;
        for (size_t thread_id = 0; thread_id < n_threads; thread_id++) {
            thread_list.emplace_back(
                [&, thread_id] {
                    for (MeasurementUnit* mu : object_allocation_list[thread_id]) {
                        if (! (mu->load_data()) ) {
                            cerr << "Error when loading data for measurement unit with id " << mu->get_meUID() << endl;
                            error_occured.store(true);
                            return;
                        }
                    }
                }
            );
        }
        // wait until all threads are finished
        for (std::thread& t : thread_list) {
            if (t.joinable())
                t.join();
        }
        // check for errors
        if (error_occured.load()) {
            return false;
        }
    } else {
        for (MeasurementUnit* mu : st__mu_list) {
            if (! (mu->load_data()) ) {
                cerr << "Error when loading data for measurement unit with id " << mu->get_meUID() << endl;
                return false;
            }
        }
    }
    return true;
}

MeasurementUnit* MeasurementUnit::GetInstancePublicID(unsigned long meUID) {
    try {
        size_t internal_id = public_to_internal_id.at(meUID); // throws std::out_of_range if public_id is unknown
		return st__mu_list[internal_id];
    } catch (std::out_of_range const&) {
		return NULL;
    }
}



// ----------------------------- //
//      Implementation of        //
//      OpenSpacePVOrWind        //
// ----------------------------- //
OpenSpacePVOrWind::OpenSpacePVOrWind(float kWp, OpenSpacePVOrWindType type)
: kWp(kWp) {
    // select the correct profile array
    if (type == OpenSpacePVOrWindType::PV) {
        if (global::pv_profiles_per_ori["S"].size() <= 0) {
            cerr << "Error: There is no south feed-in profile given!" << endl;
            cerr << "It is required for open-space installations." << endl;
            throw runtime_error("There is no south feed-in profile given!");
        }
        profile_data = global::pv_profiles_per_ori["S"][0];
    } else {
        profile_data = global::wind_profile;
    }
}

bool OpenSpacePVOrWind::compute_next_value(unsigned long ts) {
    if (ts <= 0 || ts > Global::get_n_timesteps()) {
        current_feedin_kW = 0.0;
        return false;
    }
    unsigned long tsID = ts - 1;
    current_feedin_kW = profile_data[tsID] * kWp;
    return true;
}
