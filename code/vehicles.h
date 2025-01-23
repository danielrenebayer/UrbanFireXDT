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
 * This struct represents one vehicle tour.
 * All time related parameters have to be aligned with the time step size of the simulation.
 */
struct VehicleTour {
    unsigned long next_tour_id;
    unsigned short day_of_week;
    unsigned int departure_ts_of_day;
    unsigned int ts_duration;
    double tour_length_km;
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

#endif

