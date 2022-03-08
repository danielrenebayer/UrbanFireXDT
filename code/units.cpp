
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
Substation::Substation(int id, string* name)
	: id(id), name(name)
{
	connected_units = new list<ControlUnit*>();
}

Substation::~Substation() {
	delete name;
	delete connected_units;
}

void Substation::add_unit(ControlUnit* unit) {
	connected_units->push_back(unit);
}

void Substation::PreInitializeStaticVariables() {
	st__substation_list_init = false;
	st__substation_list = NULL;
}

void Substation::InitializeStaticVariables(int n_substations) {
	st__substation_list_init = true;
	st__substation_list = new Substation*[n_substations];
}

void Substation::VacuumStaticVariables() {
	delete[] st__substation_list;
	st__substation_list = NULL;
	st__substation_list_init = false;
}




// ----------------------------- //
//      Implementation of        //
//         ControlUnit           //
// ----------------------------- //
ControlUnit::ControlUnit(int unitID, int substation_id)
    : unitID(unitID), higher_level_subst(substation_id) // TODO das geht so nicht!
{
	connected_units = new list<MeasurementUnit*>();
}

ControlUnit::~ControlUnit() {
	delete connected_units;
}

void ControlUnit::add_unit(MeasurementUnit* unit) {
	connected_units->push_back(unit);
}

void ControlUnit::PreInitializeStaticVariables() {
	st__cu_list_init = false;
	st__cu_list = NULL;
}

void ControlUnit::InitializeStaticVariables(int n_CUs) {
	st__cu_list_init = true;
	st__cu_list = new ControlUnit*[n_CUs];
}

void ControlUnit::VacuumStaticVariables() {
	delete[] st__cu_list;
	st__cu_list = NULL;
	st__cu_list_init = false;
}




// ----------------------------- //
//      Implementation of        //
//       MeasurementUnit         //
// ----------------------------- //

MeasurementUnit::MeasurementUnit(int meloID, string * melo, int locID) :
    locationID(locID), meloID(meloID), melo(melo) {
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
    data_timestepID    = new int[global::n_timesteps];    //int32    -> 4 bytes
	data_value_demand  = new float[global::n_timesteps];  //float32  -> 4 bytes
	data_value_feedin  = new float[global::n_timesteps];  //float32  -> 4 bytes
    //data_status_demand = new char[global::n_timesteps]; //int8     -> 1 byte
	//data_status_feedin = new char[global::n_timesteps]; //int8     -> 1 bytes

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
		for (int r = 0; r < global::n_timesteps; r++) {
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

void MeasurementUnit::PreInitializeStaticVariables() {
	st__mu_list_init = false;
	st__mu_list = NULL;
}

void MeasurementUnit::InitializeStaticVariables(int n_MUs) {
	st__mu_list_init = true;
	st__mu_list = new MeasurementUnit*[n_MUs];
}

void MeasurementUnit::VacuumStaticVariables() {
	delete[] st__mu_list;
	st__mu_list = NULL;
	st__mu_list_init = false;
}
