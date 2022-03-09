
#include "components.h"

// ----------------------------- //
//      Implementation of        //
//         ComponentPV           //
// ----------------------------- //

ComponentPV::ComponentPV(float kWp) : kWp(kWp) {
    currentGeneration_kW = 0;
}

inline float ComponentPV::get_kWp() const {
    return kWp;
}
float ComponentPV::get_currentGeneration_kW() {
    // TODO
    return 0;
}
void ComponentPV::calculateCurrentFeedin() {
    // TODO
    // currentGeneration_kW = ...
    // bisher: current_generation_kW = main.pv_wind_data_storage.currentPVPenetration*kWp;
}





// ----------------------------- //
//      Implementation of        //
//         ComponentBS           //
// ----------------------------- //

ComponentBS::ComponentBS(float maxE_kWh, float maxP_kW,
        float discharge_rate_per_step, float efficiency) : maxE_kWh(maxE_kWh),maxP_kW(maxP_kW),
        discharge_rate_per_step(discharge_rate_per_step), efficiency(efficiency) {

    SOC               = 0;
    currentE_kWh      = 0;
    currentP_kW       = 0;
    charge_request_kW = 0;

}
inline float ComponentBS::get_SOC() const {
    return SOC;
}
inline float ComponentBS::get_currentCharge_kWh() const {
    return currentE_kWh;
}
inline float ComponentBS::get_currentLoad_kW() const {
    return currentP_kW;
}
inline void  ComponentBS::set_chargeRequest(float requested_charge_kW) {
    charge_request_kW = requested_charge_kW;
}
void ComponentBS::calculateActions() {
    // TODO
    /* bisher:
    double timestep_size_in_h = 1;//main.simulationController.timestepsize_in_h;
    double new_charge_kWh;

    current_P_load_kW = 0;

    // Self-discharge
    current_E_charge_kWh -= discharge_rate_per_step * current_E_charge_kWh;
    charge_request_kW     = charge_request_kW * efficiency;

    if (status_output) {
        System.out.println("Bat. state start: charge request: " + Double.toString(charge_request_kW) + 
        ", curr. charge: " + Double.toString(current_E_charge_kWh));
    }

    // Charging and discharging
    if (!charge_request_set) {
        System.out.println("Important warning: charge_request not set!");
    } else {
        if (charge_request_kW > 0) {
            // charging requested
            if (charge_request_kW > maxP_kW)
                charge_request_kW = maxP_kW;
            new_charge_kWh = current_E_charge_kWh + timestep_size_in_h*charge_request_kW;
            if (new_charge_kWh > maxE_kWh)
                new_charge_kWh = maxE_kWh;
            current_P_load_kW = (new_charge_kWh - current_E_charge_kWh)/timestep_size_in_h;
            current_E_charge_kWh = new_charge_kWh;
        } else if (charge_request_kW < 0) {
            // discharging requested
            if (-charge_request_kW > maxP_kW)
                charge_request_kW = -maxP_kW;
            new_charge_kWh = current_E_charge_kWh + timestep_size_in_h*charge_request_kW;
            if (new_charge_kWh < 0)
                new_charge_kWh = 0;
            current_P_load_kW    = (new_charge_kWh - current_E_charge_kWh)/timestep_size_in_h;
            current_E_charge_kWh = new_charge_kWh;
        }
        charge_request_set = false;
    }

    if (status_output) {
        System.out.println("Bat. state end: charge request: curr. charge: " + Double.toString(current_E_charge_kWh)
        + ", current_P_load_kW: " + Double.toString(current_P_load_kW));
    }

    charge_actual_set = true;
    SOC = current_E_charge_kWh / maxE_kWh;
    */
}

