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
#include <random>
#include <vector>

// The following classes are defined in this header file:
class EVFSM;

#include "components.h"


/*!
 * This struct represents one vehicle tour.
 * All time related parameters have to be aligned with the time step size of the simulation.
 */
struct VehicleTour {
    unsigned long next_tour_id;
    unsigned short day_of_week;
    unsigned int departure_ts_of_day;
    unsigned int ts_duration;
    double trip_length_km;
    bool   with_work;
};

/*!
 * This struct represents the state of the EV, that is modeled as a finite-state machine.
 */
enum struct EVState {
    ConnectedAtHome,
    DisconnectedAtHome,
    Driving,
    ChargingOnTheWay
};

/*!
 * This struct represents the state of the EV if it is connectet to the CS at home, that is modeled as a finite-state machine.
 */
enum struct EVStateIfConnAtHome {
    ChargingRequired,
    ChargingPossible,
    DischargingPossible,
    BothPossible
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
        EVStateIfConnAtHome get_current_state_icah() const { return current_state_icah; } ///< Returns the current state of the EV iff it is connected at home
        double get_current_charging_power() const; ///< Gets the current charging power in kW; Only valid after calling `set_current_charging_power()`; This value can be smaller than the power set, as it could be that the battery cannot fulfil the request completly
        // TODO: Do we need both methods?
        //double get_min_curr_charging_power_kW() const; ///< Returns the minimal charging power at the current time step. The charging station requires at least this portion to fulfil 
        double get_max_curr_charging_power_kW() const { return max_curr_available_p_kW; } ///< Returns the maximal available charing power at the current time step.
        // modifiers (on structural level of the simulation)
        void add_weekly_tour(unsigned short weekday, unsigned int departure_ts_of_day, unsigned int ts_duration, double trip_length_km, bool with_work); ///< This method adds a home-centered car tour to the current car. All parameters that represent a time must have the same alignment as the global time information.
        void resetInternalState(); ///< Resets the internal state
        // modifiers (in the course of simulation time)
        void setCarStateForTimeStep(unsigned long ts, int dayOfWeek_l, int hourOfDay_l);
        void set_current_charging_power(double power_kW); ///< Sets the current charging power in kW
        // static methods
        static void AddWeeklyTour(unsigned long carID, unsigned short weekday, unsigned int departure_ts_of_day, unsigned int ts_duration, double trip_length_km, bool with_work); ///< This class method adds a home-centered car tour to the car with ID carID. All parameters that represent a time must have the same alignment as the global time information.
        static void VaccuumStaticVariables();
        static void SetSeed(unsigned int seed); ///< Sets the seed for the EVFSM-class random number generator

    private:
        // constant members
        const unsigned long carID;       ///< ID of the car. The ID does not necessarily be successive but it must be unique.
        const float econs_kWh_per_km;    ///< Energy consumption of the current EV in kWh/km
        ComponentCS const* homeStation;  ///< Reference to the home station of the EV
        // variable members (constant during a simulation run)
        std::vector<std::vector<VehicleTour>*> list_of_tours_pd; ///< Vector of vector of tours, one individual vector per week day (0 -> monday, 6 -> friday)
        std::vector<VehicleTour*> list_of_all_tours; ///< Vector of all tours. Is the same as `list_of_tours_pd` in a flattened order, i.e., witout day information.
        // attached members
        ComponentBS* battery;
        // variable members, variable during a simulation run
        EVState current_state;           ///< Internal current state of the EV
        EVStateIfConnAtHome current_state_icah; ///< Internal current state of the EV iff it is connected at home
        VehicleTour* current_tour;       ///< Reference to the current tour
        VehicleTour* next_tour;          ///< Reference to the next tour
        unsigned int ts_since_departure; ///< Number of time steps passed until the departure of the current tour (Only valid if a tour is ongoing)
        float energy_demand_per_tour_ts; ///< The mean energy demand per tour time step. This is the demand of the total tour divided by the number of time steps of the tour -> We assume a linear decay of the battery SOC, ignoring stops
        double max_curr_available_p_kW;  ///< The currently maximal available charging power for this vehicle
        //
        // class members
        static std::map<unsigned long, EVFSM*> list_of_cars;
        static std::default_random_engine*            random_generator; ///< Generator required for random sampling
        static std::uniform_real_distribution<float>* distribution    ; ///< Required for random sampling
};


#endif

