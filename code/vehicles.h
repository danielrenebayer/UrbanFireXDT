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
    unsigned long next_tour_id;
    unsigned short day_of_week;
    unsigned int departure_ts_of_day;
    unsigned int ts_duration;
    double tour_length_km;
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

