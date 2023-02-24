
#include "units.h"

#include <algorithm> /* using min() */
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <list>

#include <cmath> /* using: pow() */

#include "global.h"
#include "output.h"
#include "sac_planning.h"

using namespace std;


// ----------------------------- //
//      Implementation of        //
//          Substation           //
// ----------------------------- //

bool Substation::st__substation_list_init    = false;
size_t Substation::st__n_substations            = 0;
size_t Substation::st__new_substation_position  = 0;
Substation** Substation::st__substation_list = NULL;

Substation::Substation(unsigned long id, string* name)
	: id(id), name(name)
{
	/*
	 * Constructor implementation of class Substation.
	 *
	 * It constructs a new instance and checks if id is the next id
	 */
	// initialize instance variables
	connected_units = new list<ControlUnit*>();

	//
	// add to class variables
	if (!st__substation_list_init) {
		cerr << "Error: global list of substations has not been initialized" << endl;
		throw runtime_error("Global list of substations has not been initialized!");
		return;
	}
	// Attention: argument id starts with 1, not 0
	if (id <= 0 || id > st__n_substations) {
		cerr << "Error when creating a substation: id <= 0 or id > n_substations" << endl;
		throw runtime_error("id <= 0 or id > n_substations");
		return;
	}
	if (id - 1 != st__new_substation_position) {
		cerr << "Error when creating a substation: id - 1 != st__new_substation_position" << endl;
		cerr << "A reason might be, that Substation IDs are not ordered sequentially!" << endl;
		throw runtime_error("Substation IDs are not ordered sequentially!");
		return;
	}
	st__substation_list[st__new_substation_position] = this;
	st__new_substation_position++;
}

Substation::~Substation() {
	delete name;
	delete connected_units;
}

void Substation::add_unit(ControlUnit* unit) {
	connected_units->push_back(unit);
}

float Substation::calc_load() {
	float total_load = 0.0;
	for (ControlUnit* cu : *connected_units) {
		total_load += cu->get_current_load_vSMeter_kW();
	}
	return total_load;
}

void Substation::InitializeStaticVariables(unsigned long n_substations) {
	st__substation_list_init = true;
	st__substation_list = new Substation*[n_substations];
	st__n_substations   = n_substations;
}

void Substation::VacuumInstancesAndStaticVariables() {
	// delete all instances
	for (unsigned long i = 0; i < st__new_substation_position; i++)
		delete st__substation_list[i];
	// delete static variables and arrays
	delete[] st__substation_list;
	st__substation_list = NULL;
	st__substation_list_init = false;
}

inline Substation* Substation::GetInstance(unsigned long id) {
	if (id > 0 && id <= st__n_substations)
		return st__substation_list[id - 1];
	else
		return NULL;
}




// ----------------------------- //
//      Implementation of        //
//         ControlUnit           //
// ----------------------------- //

bool ControlUnit::st__cu_list_init     = false;
size_t ControlUnit::st__n_CUs             = 0;
size_t ControlUnit::st__new_CU_position   = 0;
ControlUnit** ControlUnit::st__cu_list = NULL;

ControlUnit::ControlUnit(unsigned long unitID, unsigned long substation_id, unsigned long locationID)
    : unitID(unitID), higher_level_subst(Substation::GetInstance(substation_id)), locationID(locationID)
{
	//
	// initialize instance variables
	connected_units = new list<MeasurementUnit*>();
	has_sim_pv      = false;
	has_sim_bs      = false;
	has_sim_hp      = false;
	has_sim_wb      = false;
	sim_comp_pv     = NULL;
	sim_comp_bs     = NULL;
	sim_comp_hp     = NULL;
	sim_comp_wb     = NULL;
	current_load_vSM_kW   = 0;
	self_produced_load_kW = 0;
    output_obj      = NULL;

	//
	// add to class variables
	if (!st__cu_list_init) {
		cerr << "Error: Global list of control units has not been initialized" << endl;
		throw runtime_error("Global list of control units has not been initialized!");
		return;
	}
	// Attention: argument unitID starts with 1, not 0
	if (unitID <= 0 || unitID > st__n_CUs) {
		cerr << "Error when creating a control unit: UnitID <= 0 or UnitID > n_CUs" << endl;
		throw runtime_error("UnitID <= 0 or UnitID > n_CUs");
		return;
	}
	if (unitID - 1 != st__new_CU_position) {
		cerr << "Error when creating a control unit: id - 1 != st__new_CU_position" << endl;
		cerr << "A reason might be, that Control Unit IDs are not ordered sequentially!" << endl;
		throw runtime_error("Control Unit IDs are not ordered sequentially!");
		return;
	}
	st__cu_list[st__new_CU_position] = this;
	st__new_CU_position++;

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
    sum_of_PV_generated_kWh   = 0.0;
    sum_of_feed_into_grid_kWh = 0.0;
    /*
    if (Global::get_comp_eval_metrics()) {
        create_history_output = true;
        history_self_prod_load_kW       = new float[Global::get_n_timesteps()];
        history_pv_generation_kW        = new float[Global::get_n_timesteps()];
        history_avg_consumption_load_kW = new float[Global::get_n_timesteps()];
    } else {
        create_history_output = false;
        history_self_prod_load_kW       = NULL;
        history_pv_generation_kW        = NULL;
        history_avg_consumption_load_kW = NULL;
    }*/
}

ControlUnit::~ControlUnit() {
	delete connected_units;
	if (has_sim_pv) delete sim_comp_pv;
	if (has_sim_bs) delete sim_comp_bs;
	if (has_sim_hp) delete sim_comp_hp;
	if (has_sim_wb) delete sim_comp_wb;
    /*
    if (create_history_output) {
        delete[] history_self_prod_load_kW;
        delete[] history_pv_generation_kW;
        delete[] history_avg_consumption_load_kW;
    }*/
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

bool ControlUnit::has_wb() {
	if (has_sim_wb)
		return true;
	for (MeasurementUnit* mu : *connected_units) {
		if (mu->has_wb())
			return true;
	}
	return false;
}

bool ControlUnit::has_bs_sim_added() {
	return has_sim_bs;
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
    if (has_sim_wb)
        combination = combination | expansion::MaskWB;
    return combination;
}

float ControlUnit::get_sim_comp_pv_kWp() {
    if (has_sim_pv)
        return sim_comp_pv->get_kWp();
    return 0;
}

float ControlUnit::get_sim_comp_bs_P_kW() {
    if (has_sim_bs)
        return sim_comp_bs->get_maxP_kW();
    return 0;
}

float ControlUnit::get_sim_comp_bs_E_kWh() {
    if (has_sim_bs)
        return sim_comp_bs->get_maxE_kWh();
    return 0;
}

float ControlUnit::get_annual_hp_el_cons() {
    return global::yearly_hp_energy_demand_kWh[locationID];
}

double ControlUnit::get_SSR() {
    double SSR = 0.0;
    if (sum_of_consumption_kWh > 0)
        SSR = sum_of_self_cons_kWh / sum_of_consumption_kWh;
    return SSR;
}

double ControlUnit::get_SCR() {
    double SCR = 0.0;
    if (sum_of_PV_generated_kWh > 0)
        SCR = sum_of_self_cons_kWh / sum_of_PV_generated_kWh;
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
    // 2. compute savings per year
    const double& saved_energy_kWh = sum_of_self_cons_kWh; // energy in kWh for which no grid demand was required
    double savings_per_year = 
        saved_energy_kWh * Global::get_demand_tariff() +
        sum_of_feed_into_grid_kWh * Global::get_feed_in_tariff();
    // 3. compute total savings
    double total_savings = 0.0;
    for (unsigned int t = 1; t <= Global::get_npv_time_horizon(); t++) {
        total_savings += savings_per_year / pow( (double)(1 + Global::get_npv_discount_rate()), (double) t );
    }
    return - investemnt_costs + total_savings;
}

string* ControlUnit::get_metrics_string() {
        double SCR = get_SCR();
        double SSR = get_SSR();
        double NPV = get_NPV();
        double bat_EFC = 0.0;
        if (has_sim_bs)
            bat_EFC = sim_comp_bs->get_current_EFC();
        string* retstr = new string;
        *retstr += to_string(unitID) + ",";
        *retstr += to_string(SCR) + ",";
        *retstr += to_string(SSR) + ",";
        *retstr += to_string(NPV) + ",";
        *retstr += to_string(sum_of_consumption_kWh) + ",";
        *retstr += to_string(sum_of_self_cons_kWh)   + ",";
        *retstr += to_string(sum_of_PV_generated_kWh)+ ",";
        *retstr += to_string(sum_of_feed_into_grid_kWh)+ ",";
        *retstr += to_string(bat_EFC);
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
        has_sim_pv  = true;
        if (Global::get_exp_pv_static_mode()) {
            sim_comp_pv = new ComponentPV(Global::get_exp_pv_kWp_static(), locationID);
        } else {
            sim_comp_pv = new ComponentPV(Global::get_exp_pv_kWp_per_m2(),
                                          Global::get_exp_pv_min_kWp_roof_sec(),
                                          Global::get_exp_pv_max_kWp_roof_sec(),
                                          locationID);
        }
    }
}

void ControlUnit::add_exp_bs() {
    if (has_bs())
        cerr << "Warning: Control unit with location id " << locationID << " already has a battery!" << endl;
    if (!has_sim_bs) {
        has_sim_bs  = true;
        sim_comp_bs = new ComponentBS(Global::get_exp_bess_kWh() /*maxE*/,
                                      Global::get_exp_bess_kW()  /*maxP*/,
                                      Global::get_exp_bess_E_P_ratio() /* E over P ratio */,
                                      0.0 /* discharge rate per step */,
                                      1.0 /* efficiency */,
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
        sim_comp_hp = new ComponentHP(this->get_annual_hp_el_cons());
    }
}

void ControlUnit::add_exp_wb() {
    if (has_wb())
        cerr << "Warning: Control unit with location id " << locationID << " already has a wallbox!" << endl;
    if (!has_sim_wb) {
        has_sim_wb  = true;
        sim_comp_wb = new ComponentWB();
    }
}

void ControlUnit::set_output_object(CUOutput* output_obj) {
    this->output_obj = output_obj;
}

void ControlUnit::set_exp_pv_params_A(float value) {
    // do something only if this unit is selected for PV expansion
    if (has_sim_pv) {
        if (Global::get_exp_pv_static_mode()) {
            delete sim_comp_pv;
            sim_comp_pv = new ComponentPV(value, locationID);
        } else {
            throw runtime_error("Error: ControlUnit::set_exp_pv_params_A() has been callen even though PV static mode is set!");
        }
    }
}

void ControlUnit::set_exp_pv_params_B(float kWp_per_m2, float min_kWp, float max_kWp) {
    // do something only if this unit is selected for PV expansion
    if (has_sim_pv) {
        if (Global::get_exp_pv_static_mode()) {
            throw runtime_error("Error: ControlUnit::set_exp_pv_params_B() has been callen even though PV static mode is not set!");
        } else {
            delete sim_comp_pv;
            sim_comp_pv = new ComponentPV(kWp_per_m2, min_kWp, max_kWp, locationID);
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
	if (has_sim_wb) {
        has_sim_wb = false;
        delete sim_comp_wb;
        sim_comp_wb = NULL;
    }
}

bool ControlUnit::compute_next_value(unsigned long ts) {
    //
    // This function computes the next value
    // for this complete control unit.
    // It also calls the methods of all connected components
    // to calculate their actions for the next step.
    //

    float current_load_all_rSMs_kW = 0.0;
    float total_consumption = 0.0;
    float load_pv = 0.0;
    float load_hp = 0.0;
    float load_wb = 0.0;
    float load_bs = 0.0;
    float bs_SOC  = 0.0;

    //
    // 1. get sum of all real smart meter values
    for (MeasurementUnit* mu : *connected_units) {
        if (!mu->compute_next_value(ts))
            return false;
        current_load_all_rSMs_kW += mu->get_current_ts_rsm_value();
    }
    current_load_vSM_kW = current_load_all_rSMs_kW;
    total_consumption   = current_load_all_rSMs_kW;
    if (total_consumption < 0) total_consumption = 0;
    // // //float load_bevore_local_pv_bess = current_load_vSM_kW;
    //
    // 2. get PV feedin
    if (has_sim_pv) {
        sim_comp_pv->calculateCurrentFeedin(ts);
        load_pv = sim_comp_pv->get_currentGeneration_kW();
        current_load_vSM_kW -= load_pv;
    }
    //
    // 3. get the demand of the heat pump
    if (has_sim_hp) {
        sim_comp_hp->calculateCurrentFeedin(ts);
        load_hp = sim_comp_hp->get_currentDemand_kW();
        current_load_vSM_kW += load_hp;
        total_consumption   += load_hp;
    }
    //
    // 4. get the effect of the e-car / wallbox
    // TODO
    // if load_wb > 0: total_consumption += load_wb; // ... only problem: wallbox feeds in energy taken from somewhere else
    //
    // 5. send situation to battery storage and get its resulting action
    if (has_sim_bs) {
        sim_comp_bs->set_chargeRequest( -current_load_vSM_kW );
        sim_comp_bs->calculateActions();
        load_bs = sim_comp_bs->get_currentLoad_kW();
        current_load_vSM_kW += load_bs;

        bs_SOC = sim_comp_bs->get_SOC();
    }

    //
    // compute self-produced load, that is directly consumed
    if (current_load_all_rSMs_kW < 0) {
        self_produced_load_kW = std::min(load_hp, load_pv-load_bs);
    } else {
        self_produced_load_kW = std::min(current_load_all_rSMs_kW + load_hp, load_pv-load_bs);
    }

    double grid_feedin_kW = 0.0;
    if (current_load_vSM_kW < 0) {
        grid_feedin_kW = -current_load_vSM_kW;
    }

    //
    // add values to summation variables
    sum_of_consumption_kWh    += Global::get_time_step_size_in_h() * total_consumption;
    sum_of_self_cons_kWh      += Global::get_time_step_size_in_h() * self_produced_load_kW;
    sum_of_PV_generated_kWh   += Global::get_time_step_size_in_h() * load_pv;
    sum_of_feed_into_grid_kWh += Global::get_time_step_size_in_h() * grid_feedin_kW;
    /*if (create_history_output) {
        history_self_prod_load_kW[ts - 1]     = ;
        history_pv_generation_kW[ ts - 1]     = load_pv;
        history_avg_consumption_load_kW[ts-1] = total_consumption;
    }*/

    //
    // output current status
    if (output_obj != NULL)
        output_obj->output_for_one_cu(
            unitID, ts,
            current_load_vSM_kW,
            current_load_all_rSMs_kW,
            self_produced_load_kW,
            load_pv, bs_SOC, load_bs,
            load_hp, load_wb);

    return true;
}

void ControlUnit::InitializeStaticVariables(unsigned long n_CUs) {
    st__cu_list_init = true;
    st__cu_list = new ControlUnit*[n_CUs];
    st__n_CUs   = n_CUs;
}

void ControlUnit::VacuumInstancesAndStaticVariables() {
    // delete all instances
    for (unsigned long i = 0; i < st__new_CU_position; i++)
        delete st__cu_list[i];
    // delete static variables and arrays
    delete[] st__cu_list;
    st__cu_list = NULL;
    st__cu_list_init = false;
}

inline ControlUnit* ControlUnit::GetInstance(unsigned long unitID) {
    if (unitID > 0 && unitID <= st__n_CUs)
        return st__cu_list[unitID - 1];
    else
        return NULL;
}

void ControlUnit::ResetAllInternalStates() {
    //
    // Call this function to reset all internal states of
    // the unit and all connected components (as long as they
    // have interal states, like a battery has the current SOC
    // as internal state).
    //
    for (unsigned long i = 0; i < st__n_CUs; i++) {
        ControlUnit* e_i = st__cu_list[i];
        //
        e_i->current_load_vSM_kW = 0.0;
        e_i->self_produced_load_kW = 0.0;
        e_i->sum_of_consumption_kWh    = 0.0;
        e_i->sum_of_self_cons_kWh      = 0.0;
        e_i->sum_of_PV_generated_kWh   = 0.0;
        e_i->sum_of_feed_into_grid_kWh = 0.0;
        //
        if (e_i->has_sim_bs) {
            e_i->sim_comp_bs->resetInternalState();
        }
    }
}

void ControlUnit::RemoveAllSimAddedComponents() {
    for (unsigned long i = 0; i < st__n_CUs; i++) {
        st__cu_list[i]->remove_sim_added_components();
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





// ----------------------------- //
//      Implementation of        //
//       MeasurementUnit         //
// ----------------------------- //

bool MeasurementUnit::st__mu_list_init   = false;
size_t MeasurementUnit::st__n_MUs           = 0;
size_t MeasurementUnit::st__new_MU_position = 0;
MeasurementUnit** MeasurementUnit::st__mu_list = NULL;

MeasurementUnit::MeasurementUnit(size_t meUID, size_t unitID, string * meterPointName, size_t locID,
                                 bool has_demand, bool has_feedin, bool has_pv_resid, bool has_pv_opens,
                                 bool has_bess,   bool has_hp,     bool has_wb,       bool has_chp) :
    meUID(meUID),
    higher_level_cu(ControlUnit::GetInstance(unitID)),
    meterPointName(meterPointName), locationID(locID) {
    //
    // initialize instance variables
    current_load_rsm_kW = 0;
    rsm_has_demand = has_demand;
    rsm_has_feedin = has_feedin;
    rsm_with_pv_residential = has_pv_resid;
    rsm_with_pv_open_space  = has_pv_opens;
    rsm_with_bess  = has_bess;
    rsm_with_hp    = has_hp;
    rsm_with_wb    = has_wb;
    rsm_with_chp   = has_chp;
    data_loaded       = false;
    data_timestepID   = NULL;
    data_value_demand = NULL;
    data_value_feedin = NULL;
    //data_status_demand=NULL;
    //data_status_feedin=NULL;
    expansion_combination = expansion::genExpCombiAsBitRepr(has_pv_resid||has_pv_opens, has_bess, has_hp, has_wb);

    //
    // add to class variables
    if (!st__mu_list_init) {
        cerr << "Error: global list of measurement units has not been initialized" << endl;
        throw runtime_error("Global list of measurement units has not been initialized!");
        return;
    }
    // Attention: argument meUID starts with 1, not 0
    if (meUID <= 0 || meUID > st__n_MUs) {
        cerr << "Error when creating a measurement unit: meUID <= 0 or meUID > n_MUs" << endl;
        throw runtime_error("meUID <= 0 or meUID > n_MUs");
        return;
    }
    if (meUID - 1 != st__new_MU_position) {
        cerr << "Error when creating a measurement unit: meUID - 1 != st__new_MU_position" << endl;
        cerr << "A reason might be, that Measurement Unit IDs are not ordered sequentially!" << endl;
        throw runtime_error("Measurement Unit IDs are not ordered sequentially!");
        return;
    }
    st__mu_list[st__new_MU_position] = this;
    st__new_MU_position++;

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

bool MeasurementUnit::load_data(const char * filepath) {
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
    smeter_input.open(filepath);
    if (!smeter_input.good()) {
        cout << "Error when connecting to the smart meter file with path " << filepath << endl;
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
            cerr << "TimestepID,Value_Demand,Status_Demand,Value_Feedin,Status_Feedin" << endl;
            return false; // destructor of smeter_input closes it by default
        }
        stringstream headerSStream(currLineString);
        string headerSplitted[5];
        for (uint col = 0; col < 5; col++) {
            std::getline( headerSStream, headerSplitted[col], ',' );
        }
        if (
            headerSplitted[0] != "TimestepID"    ||
            headerSplitted[1] != "Value_Demand"  ||
            headerSplitted[2] != "Status_Demand" ||
            headerSplitted[3] != "Value_Feedin"  ||
            headerSplitted[4] != "Status_Feedin"
        ) {
            cerr << "There is one smart meter data file with a wrong header!" << endl;
            cerr << "The data file must have the following columns (mind the order):" << endl;
            cerr << "TimestepID,Value_Demand,Status_Demand,Value_Feedin,Status_Feedin" << endl;
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

inline const std::string * MeasurementUnit::get_meterPointName() const {
    return meterPointName;
}

inline unsigned long MeasurementUnit::get_meUID() const{
    return meUID;
}

inline unsigned long MeasurementUnit::get_locationID() const {
    return locationID;
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

void MeasurementUnit::InitializeStaticVariables(unsigned long n_MUs) {
    st__mu_list_init = true;
    st__mu_list = new MeasurementUnit*[n_MUs];
    st__n_MUs   = n_MUs;
}

void MeasurementUnit::VacuumInstancesAndStaticVariables() {
    // delete all instances
    for (size_t i = 0; i < st__new_MU_position; i++)
        delete st__mu_list[i];
    // delete static variables and arrays
    delete[] st__mu_list;
    st__mu_list = NULL;
    st__mu_list_init = false;
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
