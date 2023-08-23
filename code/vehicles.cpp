#include "vehicles.h"


// ----------------------------- //
//      Implementation of        //
//            EVFSM              //
// ----------------------------- //

std::map<unsigned long, EVFSM*> EVFSM::list_of_cars;

EVFSM::EVFSM(unsigned long carID, ComponentCS* homeStation) : carID(carID), homeStation(homeStation) {
    // Register this object in the global list of cars
    if (EVFSM::list_of_cars.contains(carID)) {
        throw runtime_error("Error: There is already an instance of class EVFSM with carID = " + std::to_string(carID));
    }
    EVFSM::list_of_cars.emplace(carID, this);
    //
    // Create battery (using the secondary constructor)
    battery = new ComponentBS(30, 0.0, 1.0, 0.25);
    // Initialize variables for the state
    current_state = EVState::AtHome;
    // Create the list of tours
    list_of_tours.reserve(6);
    for (unsigned int day = 0; day < 7; day++) {
        list_of_tours.push_back( new std::vector<VehicleTour> () );
    }
}

EVFSM::~EVFSM() {
    delete battery;
    for (auto day_tours : list_of_tours) {
        delete day_tours;
    }
    list_of_tours.clear();
}

void EVFSM::add_weekly_tour(
    unsigned short weekday,
    unsigned int departure_ts_of_day,
    unsigned int ts_duration,
    double trip_length_km,
    bool with_work)
{
    if (weekday > 6)
        throw runtime_error("A weekday > 6 is not possible!");
    
    list_of_tours[weekday]->emplace_back(departure_ts_of_day, ts_duration, trip_length_km, with_work);
}

void EVFSM::AddWeeklyTour(
    unsigned long carID,              unsigned short weekday,
    unsigned int departure_ts_of_day, unsigned int ts_duration,
    double trip_length_km,            bool with_work)
{
    EVFSM::list_of_cars.at(carID)->add_weekly_tour(weekday, departure_ts_of_day, ts_duration, trip_length_km, with_work);
}

void EVFSM::VaccuumStaticVariables() {
    EVFSM::list_of_cars.clear();
}
