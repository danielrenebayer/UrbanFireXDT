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

/*!
 * This struct represents one weekly-repeating vehicle tour.
 * All time related parameters have to be aligned with the time step size of the simulation.
 */
struct WeeklyVehicleTour {
    unsigned long next_tour_id; ///< 0, if there is no new tour
    short day_of_week; ///< The day of the week a tour starts. Note: 0 is Monday (like in Python, not like tm_wday where 0 is a Sunday)
    unsigned long departure_ts_of_day; ///< Number of time steps after midnight where the tour starts
    unsigned long ts_duration; ///< Tour duration in time steps, arrival at the beginning of a new time step (i.e., ts_duration = 1 --> arrival after one time step)
    double tour_length_km; ///< Total length of the tour in km
    bool   with_work;
};

/*!
 * Representation of one individual tour that happens only once in a simulation span.
 */
struct SingleVehicleTour {
    unsigned long ts_start; ///< Time step, when tour starts
    unsigned long ts_arrival; ///< Time step, when tour ends (i.e., arrival of the car at the beginning of the time step)
    double energy_consumption_kWh; ///< Energy consumption for this trip
    struct WeeklyVehicleTour* weekly_tour; ///< Reference to the corresponding weekly tour

    SingleVehicleTour(unsigned long ts_start, unsigned long ts_arrival, double econs, struct WeeklyVehicleTour* weekly_tour) : 
    ts_start(ts_start), ts_arrival(ts_arrival), energy_consumption_kWh(econs), weekly_tour(weekly_tour) {}
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

#endif

