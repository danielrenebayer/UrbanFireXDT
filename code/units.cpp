
#include "units.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>

#include "global.h"
#include "simulation_setup.h"

using namespace std;


// ----------------------------- //
//      Implementation of        //
//          Substation           //
// ----------------------------- //

bool Substation::st__substation_list_init    = false;
int Substation::st__n_substations            = 0;
int Substation::st__new_substation_position  = 0;
Substation** Substation::st__substation_list = NULL;

Substation::Substation(int id, string* name)
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
		throw "Global list of substations has not been initialized!";
		return;
	}
	// Attention: argument id starts with 1, not 0
	if (id <= 0 || id > st__n_substations) {
		cerr << "Error when creating a substation: id <= 0 or id > n_substations" << endl;
		throw "id <= 0 or id > n_substations";
		return;
	}
	if (id - 1 != st__new_substation_position) {
		cerr << "Error when creating a substation: id - 1 != st__new_substation_position" << endl;
		cerr << "A reason might be, that Substation IDs are not ordered sequentially!" << endl;
		throw "Substation IDs are not ordered sequentially!";
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

void Substation::InitializeStaticVariables(int n_substations) {
	st__substation_list_init = true;
	st__substation_list = new Substation*[n_substations];
	st__n_substations   = n_substations;
}

void Substation::VacuumInstancesAndStaticVariables() {
	// delete all instances
	for (int i = 0; i < st__new_substation_position; i++)
		delete st__substation_list[i];
	// delete static variables and arrays
	delete[] st__substation_list;
	st__substation_list = NULL;
	st__substation_list_init = false;
}

inline Substation* Substation::GetInstance(int id) {
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
int ControlUnit::st__n_CUs             = 0;
int ControlUnit::st__new_CU_position   = 0;
ControlUnit** ControlUnit::st__cu_list = NULL;

ControlUnit::ControlUnit(int unitID, int substation_id)
    : unitID(unitID), higher_level_subst(Substation::GetInstance(substation_id))
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

	//
	// add to class variables
	if (!st__cu_list_init) {
		cerr << "Error: Global list of control units has not been initialized" << endl;
		throw "Global list of control units has not been initialized!";
		return;
	}
	// Attention: argument unitID starts with 1, not 0
	if (unitID <= 0 || unitID > st__n_CUs) {
		cerr << "Error when creating a control unit: UnitID <= 0 or UnitID > n_CUs" << endl;
		throw "UnitID <= 0 or UnitID > n_CUs";
		return;
	}
	if (unitID - 1 != st__new_CU_position) {
		cerr << "Error when creating a substation: id - 1 != st__new_substation_position" << endl;
		cerr << "A reason might be, that Control Unit IDs are not ordered sequentially!" << endl;
		throw "Control Unit IDs are not ordered sequentially!";
		return;
	}
	st__cu_list[st__new_CU_position] = this;
	st__new_CU_position++;
}

ControlUnit::~ControlUnit() {
	delete connected_units;
	if (has_sim_pv) delete sim_comp_pv;
	if (has_sim_bs) delete sim_comp_bs;
	if (has_sim_hp) delete sim_comp_hp;
	if (has_sim_wb) delete sim_comp_wb;
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

int ControlUnit::get_exp_combi_bit_repr() {
	int combination = 0;
	for (MeasurementUnit* mu : *connected_units) {
		combination = combination | mu->get_expansion_combination();
	}
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

void ControlUnit::add_exp_pv() {
	if (!has_sim_pv) {
		has_sim_pv  = true;
		sim_comp_pv = new ComponentPV(Global::get_exp_pv_kWp());
	}
}

void ControlUnit::add_exp_bs() {
	if (!has_sim_bs) {
		has_sim_bs  = true;
		sim_comp_bs = new ComponentBS(Global::get_exp_bess_kWh(), Global::get_exp_bess_kW(), 0.0, 1.0, Global::get_exp_bess_start_soc());
	}
}

void ControlUnit::add_exp_hp() {
	if (!has_sim_hp) {
		has_sim_hp  = true;
		sim_comp_hp = new ComponentHP();
	}
}

void ControlUnit::add_exp_wb() {
	if (!has_sim_wb) {
		has_sim_wb  = true;
		sim_comp_wb = new ComponentWB();
	}
}

bool ControlUnit::compute_next_value(int ts) {
	//
	// This function computes the next value
	// for this complete control unit.
	// It also calls the methods of all connected components
	// to calculate their actions for the next step.
	//

	//
	// 1. get sum of all real smart meter values
	current_load_vSM_kW = 0.0;
	for (MeasurementUnit* mu : *connected_units) {
		if (!mu->compute_next_value(ts))
			return false;
		current_load_vSM_kW += mu->get_current_ts_rsm_value();
	}
	float load_bevore_local_pv_bess = current_load_vSM_kW;
	//
	// 2. get PV feedin
	if (has_sim_pv) {
		sim_comp_pv->calculateCurrentFeedin(ts);
		current_load_vSM_kW -= sim_comp_pv->get_currentGeneration_kW();
	}
	//
	// 3. send situation to battery storage and get its resulting action
	if (has_sim_bs) {
		sim_comp_bs->set_chargeRequest( -current_load_vSM_kW );
		sim_comp_bs->calculateActions();
		current_load_vSM_kW += sim_comp_bs->get_currentLoad_kW();
	}
	// TODO: implement 4 and 5: sector coupling
	// 4. the effect of the heat pump
	// 5. the effect of the car

	//
	// compute self-produced load, that is directly consumed
	self_produced_load_kW = 0;
	if (load_bevore_local_pv_bess > 0) {
		if (current_load_vSM_kW > 0) {
			// less local production than than consumption -> additional supply from the grid
			self_produced_load_kW = load_bevore_local_pv_bess - current_load_vSM_kW;
			// only check
			if (self_produced_load_kW < 0)
				cerr << "Error: load_self_produced < 0!" << endl;
		} else {
			// more (or equal) local production than consumption -> additional feedin (or nothing)
			self_produced_load_kW = load_bevore_local_pv_bess;
		}
	}

	return true;
}

void ControlUnit::InitializeStaticVariables(int n_CUs) {
	st__cu_list_init = true;
	st__cu_list = new ControlUnit*[n_CUs];
	st__n_CUs   = n_CUs;
}

void ControlUnit::VacuumInstancesAndStaticVariables() {
	// delete all instances
	for (int i = 0; i < st__new_CU_position; i++)
		delete st__cu_list[i];
	// delete static variables and arrays
	delete[] st__cu_list;
	st__cu_list = NULL;
	st__cu_list_init = false;
}

inline ControlUnit* ControlUnit::GetInstance(int unitID) {
	if (unitID > 0 && unitID <= st__n_CUs)
		return st__cu_list[unitID - 1];
	else
		return NULL;
}





// ----------------------------- //
//      Implementation of        //
//       MeasurementUnit         //
// ----------------------------- //

bool MeasurementUnit::st__mu_list_init   = false;
int MeasurementUnit::st__n_MUs           = 0;
int MeasurementUnit::st__new_MU_position = 0;
MeasurementUnit** MeasurementUnit::st__mu_list = NULL;

MeasurementUnit::MeasurementUnit(int meloID, int unitID, string * melo, int locID,
								 bool has_demand, bool has_feedin, bool has_pv_resid, bool has_pv_opens,
								 bool has_bess,   bool has_hp,     bool has_wb,       bool has_chp) :
	meloID(meloID),
	higher_level_cu(ControlUnit::GetInstance(unitID)),
	melo(melo), locationID(locID) {
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
		throw "Global list of measurement units has not been initialized!";
		return;
	}
	// Attention: argument meloID starts with 1, not 0
	if (meloID <= 0 || meloID > st__n_MUs) {
		cerr << "Error when creating a measurement unit: meloID <= 0 or meloID > n_MUs" << endl;
		throw "meloID <= 0 or meloID > n_MUs";
		return;
	}
	if (meloID - 1 != st__new_MU_position) {
		cerr << "Error when creating a measurement unit: meloID - 1 != st__new_MU_position" << endl;
		cerr << "A reason might be, that Measurement Unit IDs are not ordered sequentially!" << endl;
		throw "Measurement Unit IDs are not ordered sequentially!";
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

    delete melo;

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
		getline( smeter_input, currLineString ); // jump first line, as this is the header
		for (int r = 0; r < Global::get_n_timesteps(); r++) {
			// iterate over every row
			getline( smeter_input, currLineString );
			stringstream currLineStream( currLineString );
			string currLineSplitted[5];
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

inline const std::string * MeasurementUnit::get_melo() const {
    return melo;
}

inline const int MeasurementUnit::get_meloID() const{
    return meloID;
}

inline const int MeasurementUnit::get_locationID() const {
    return locationID;
}

bool MeasurementUnit::compute_next_value(int ts) {
	if (ts <= 0 || ts > Global::get_n_timesteps()) {
		current_load_rsm_kW = 0.0;
		return false;
	}
	int tsID = ts - 1;
	current_load_rsm_kW = data_value_demand[tsID] - data_value_feedin[tsID];
	return true;
}

void MeasurementUnit::InitializeStaticVariables(int n_MUs) {
	st__mu_list_init = true;
	st__mu_list = new MeasurementUnit*[n_MUs];
	st__n_MUs   = n_MUs;
}

void MeasurementUnit::VacuumInstancesAndStaticVariables() {
	// delete all instances
	for (int i = 0; i < st__new_MU_position; i++)
		delete st__mu_list[i];
	// delete static variables and arrays
	delete[] st__mu_list;
	st__mu_list = NULL;
	st__mu_list_init = false;
}
