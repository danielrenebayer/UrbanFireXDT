#include "vehicles.h"

#include <cmath>
#include <iostream>
#include <random>

#include "components.h"
#include "global.h"


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
    current_state_icah = EVStateIfConnAtHome::ChargingPossible;
    current_tour       = NULL;
    next_tour          = NULL;
    ts_since_departure = 0;
    energy_demand_per_tour_ts = 0.0;
    max_curr_available_p_kW   = 0.0;
    sum_of_driving_distance_km    = 0.0;
    sum_of_E_used_for_driving_kWh = 0.0;
    sum_of_E_charged_home_kWh     = 0.0;
    sum_of_E_discharged_home_kWh  = 0.0;
    sum_of_ts_EV_is_connected_kWh = 0;
    // Create the list of tours
    list_of_tours_pd.reserve(6);
    for (unsigned int day = 0; day < 7; day++) {
        list_of_tours_pd.push_back( new std::vector<VehicleTour> () );
    }
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

float EVFSM::get_current_charging_power() const {
    return battery->get_currentLoad_kW();
}

std::string* EVFSM::get_metrics_string_annual() {
    std::string* retstr = new string;
    *retstr += std::to_string(carID) + ",";
    *retstr += std::to_string(homeStation->get_control_unit_id()) + ",";
    *retstr += std::to_string(sum_of_driving_distance_km)     + ",";
    *retstr += std::to_string(sum_of_E_used_for_driving_kWh)  + ",";
    *retstr += std::to_string(sum_of_E_charged_home_kWh)      + ",";
    *retstr += std::to_string(sum_of_E_discharged_home_kWh)   + ",";
    *retstr += std::to_string(sum_of_ts_EV_is_connected_kWh);
    return retstr;
}

void EVFSM::add_weekly_tour(
    unsigned short weekday,
    unsigned int departure_ts_of_day,
    unsigned int ts_duration,
    double tour_length_km,
    bool with_work)
{
    if (weekday > 6)
        throw runtime_error("Error when adding a new vehicle tour: A weekday > 6 is not possible!");
    
    unsigned long this_new_tour_id = list_of_all_tours.size();
    if (list_of_all_tours.size() > 0) {
        VehicleTour* last_know_tour = list_of_all_tours.back();
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
        VehicleTour* first_known_tour = list_of_all_tours.front();
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
    VehicleTour& new_tour = list_of_tours_pd[weekday]->emplace_back(0, weekday, departure_ts_of_day, ts_duration, tour_length_km, with_work);
    list_of_all_tours.push_back( &new_tour );
}

void EVFSM::resetInternalState() {
    battery->resetInternalState();
    current_state      = EVState::ConnectedAtHome;
    current_state_icah = EVStateIfConnAtHome::ChargingPossible;
    current_tour       = NULL;
    next_tour          = NULL;
    ts_since_departure = 0;
    energy_demand_per_tour_ts = 0.0;
    max_curr_available_p_kW   = 0.0;
    sum_of_driving_distance_km    = 0.0;
    sum_of_E_used_for_driving_kWh = 0.0;
    sum_of_E_charged_home_kWh     = 0.0;
    sum_of_E_discharged_home_kWh  = 0.0;
    sum_of_ts_EV_is_connected_kWh = 0;
}

void EVFSM::setCarStateForTimeStep(unsigned long ts, unsigned int dayOfWeek_l, unsigned int hourOfDay_l) {
    // if there is a currently active tour: check, if it reached home
    if (current_tour != NULL) {
        ts_since_departure += 1;
        // remove consumed energy per step from the battery
        battery->set_chargeRequest(-energy_demand_per_tour_ts);
        battery->calculateActions();
        // check, if this tour is finished?
        if (ts_since_departure >= current_tour->ts_duration) {
            // add tour length of the finished tour to the total driving distance
            sum_of_driving_distance_km += current_tour->tour_length_km;
            // add consumed electricty for driving
            sum_of_E_used_for_driving_kWh += current_tour->tour_length_km * econs_kWh_per_km;
            // remove current tour
            current_tour = NULL;
            // Calculate if car will be connected or not (depending on the current SOC and its connection probability)
            //   Always connect if SOC <= 0.35
            if (battery->get_SOC() <= 0.35) {
                current_state = EVState::ConnectedAtHome;
            } else {
                // sample randomly if car will be conntected or not
                if (Global::get_ev_plugin_probability() >= 1.0 || (*distribution)(*random_generator) <= Global::get_ev_plugin_probability())
                    current_state = EVState::ConnectedAtHome;
                else
                    current_state = EVState::DisconnectedAtHome;
            }
        }
    }
    // does a new tour start?
    if (current_tour == NULL) {
        // loop over all tours on this day, is there a tour starting right now?
        for ( VehicleTour &vt : *(list_of_tours_pd)[dayOfWeek_l] ) {
            if (vt.departure_ts_of_day == hourOfDay_l) { // mind the shift: Car tours are left-aligned, hours are right aligned
                current_tour = &vt;
                ts_since_departure = 0;
                // compute (mean) energy demand per tour time step
                energy_demand_per_tour_ts = (float) (vt.tour_length_km * econs_kWh_per_km) / (float) (vt.ts_duration);
                next_tour = list_of_all_tours[vt.next_tour_id];
                current_state = EVState::Driving;
                break;
            }
        }
    }
    //
    // if connected, determine the internal state
    if (current_state == EVState::ConnectedAtHome) {
        sum_of_ts_EV_is_connected_kWh += 1;
        // the goal is: the car needs
        //   1) either at least 35 percent SOC when leaving,
        //   2) or at least as much energy in the battery required for the next trip + 20 percent points
        if (battery->get_SOC() < 1.0) {
            current_state_icah = EVStateIfConnAtHome::ChargingPossible;
            // TODO: BothPossible could also be the case in this situation!
            // look at next_tour for this information
            //
            // Compute the maximal chargable power at the current time step
            float free_space_kWh = battery->get_maxE_kWh() - battery->get_currentCharge_kWh();
            max_curr_available_p_kW = free_space_kWh / Global::get_time_step_size_in_h();
        } else {
            current_state_icah = EVStateIfConnAtHome::DischargingPossible;
            max_curr_available_p_kW = 0.0;
        }
    } else {
        max_curr_available_p_kW = 0.0;
    }
}

void EVFSM::set_current_charging_power(float power_kW) {
    battery->set_chargeRequest(power_kW);
    battery->calculateActions();
    // update the amount of charged / discharged electricity
    float cload = battery->get_currentLoad_kW();
    if ( cload > 0 ) {
        sum_of_E_charged_home_kWh    += cload * Global::get_time_step_size_in_h();
    } else if ( cload < 0 ) {
        sum_of_E_discharged_home_kWh -= cload * Global::get_time_step_size_in_h();
    }
}

void EVFSM::AddWeeklyTour(
    unsigned long carID,              unsigned short weekday,
    unsigned int departure_ts_of_day, unsigned int ts_duration,
    double tour_length_km,            bool with_work)
{
    EVFSM::list_of_cars.at(carID)->add_weekly_tour(weekday, departure_ts_of_day, ts_duration, tour_length_km, with_work);
}

void EVFSM::VaccuumStaticVariables() {
    EVFSM::list_of_cars.clear();
}

void EVFSM::SetSeed(unsigned int seed) {
    random_generator->seed(seed);
}
