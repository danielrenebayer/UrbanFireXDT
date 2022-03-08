
#include "units.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>

#include "global.h"

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
}

void ControlUnit::add_unit(MeasurementUnit* unit) {
	connected_units->push_back(unit);
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




// ----------------------------- //
//      Implementation of        //
//       MeasurementUnit         //
// ----------------------------- //

bool MeasurementUnit::st__mu_list_init   = false;
int MeasurementUnit::st__n_MUs           = 0;
int MeasurementUnit::st__new_MU_position = 0;
MeasurementUnit** MeasurementUnit::st__mu_list = NULL;

MeasurementUnit::MeasurementUnit(int meloID, string * melo, int locID) :
    locationID(locID), meloID(meloID), melo(melo) {
	//
	// initialize instance variables
    current_load_rsm_kW = 0;
    rsm_has_demand = false;
    rsm_has_feedin = false;
    rsm_with_pv = false;
    rsm_with_bess = false;
    rsm_with_hp = false;
    rsm_with_wb = false;
    data_loaded       = false;
    data_timestepID   = NULL;
	data_value_demand = NULL;
	data_value_feedin = NULL;
    //data_status_demand=NULL;
    //data_status_feedin=NULL;

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
    melo = NULL;
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
	
	for (int i = 0; i < 10; i++) {
		printf("%6i %6.4f %6.4f \n",
			data_timestepID[i],
			data_value_demand[i],
			data_value_feedin[i]);
	}
	cout << endl;

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

inline bool MeasurementUnit::has_feedin() {
    return rsm_has_feedin;
}

inline bool MeasurementUnit::has_demand() {
    return rsm_has_demand;
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
