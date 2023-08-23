/*
 * vehicles.h
 *
 * It contains all classes and structures for
 * modeling and representing
 *  1) the (electric) vehicle driving profiles
 *  2) the ev itself modeled as a finite-state machine
 */

#ifndef VEHICLES_H
#define VEHICLES_H

#include <map>
#include <vector>

// The following classes are defined in this header file:
class EVFSM;

#include "components.h"


/*!
 * This struct represents one vehicle tour.
 * All time related parameters have to be aligned with the time step size of the simulation.
 */
struct VehicleTour {
    unsigned int departure_ts_of_day;
    unsigned int ts_duration;
    double trip_length_km;
    bool   with_work;
};

/*!
 * This struct represents the state of the EV, that is modeled as a finite-state machine.
 */
enum struct EVState {
    ChargingAtHome,
    AtHome,
    Driving,
    ChargingOnTheWay
};

/*!
 * This class represents a EV (electric vehicle) modeled as a finite-state machine.
 *
 * It is a finite-state machine with the states:
 *  - charging at home
 *  - at home (bu not charging)
 *  - driving (i.e., away from home)
 *  - charging somewhere on the way
 * 
 * Moreover, it has a battery component (the same as it is used for simulating the residential batteries).
 * 
 * Finally, it holds an attached driving profile.
 */
class EVFSM {

    public:
        EVFSM(unsigned long carID, ComponentCS* homeStation);
        ~EVFSM();
        // getters
        EVState get_current_state() const { return current_state; }
        // modifiers
        void add_weekly_tour(unsigned short weekday, unsigned int departure_ts_of_day, unsigned int ts_duration, double trip_length_km, bool with_work); ///< This method adds a home-centered car tour to the current car. All parameters that represent a time must have the same alignment as the global time information.
        // static methods
        static void AddWeeklyTour(unsigned long carID, unsigned short weekday, unsigned int departure_ts_of_day, unsigned int ts_duration, double trip_length_km, bool with_work); ///< This class method adds a home-centered car tour to the car with ID carID. All parameters that represent a time must have the same alignment as the global time information.
        static void VaccuumStaticVariables();

    private:
        // constant members
        const unsigned long carID;
        ComponentCS const* homeStation;
        // variable members
        EVState current_state;
        std::vector<std::vector<VehicleTour>*> list_of_tours;
        // attached members
        ComponentBS* battery;
        //
        // class members
        static std::map<unsigned long, EVFSM*> list_of_cars;
};


#endif

