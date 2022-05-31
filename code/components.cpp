
#include "components.h"

// ----------------------------- //
//      Implementation of        //
//         ComponentPV           //
// ----------------------------- //

ComponentPV::ComponentPV(float kWp) : kWp(kWp) {
    currentGeneration_kW = 0;
}

void ComponentPV::calculateCurrentFeedin(int ts) {
    int tsID = ts - 1;
    currentGeneration_kW = kWp * global::pv_profile[tsID];
}





// ----------------------------- //
//      Implementation of        //
//         ComponentBS           //
// ----------------------------- //

ComponentBS::ComponentBS(float maxE_kWh, float maxP_kW,
        float discharge_rate_per_step, float efficiency, float initial_SoC) : maxE_kWh(maxE_kWh),maxP_kW(maxP_kW),
        discharge_rate_per_step(discharge_rate_per_step), efficiency(efficiency),
        initial_SoC(initial_SoC) {

    SOC               = 0;
    currentE_kWh      = 0;
    currentP_kW       = 0;
    charge_request_kW = 0;

    if (initial_SoC > 0) {
        SOC = initial_SoC;
        currentE_kWh = maxE_kWh * initial_SoC;
    }

}

void ComponentBS::calculateActions() {
    float timestep_size_in_h = 1;//main.simulationController.timestepsize_in_h;
    float new_charge_kWh;

    currentP_kW = 0;

    // Calculate Self-discharge
    currentE_kWh -= discharge_rate_per_step * currentE_kWh;
    // Calculate efficiency
    charge_request_kW = charge_request_kW * efficiency;

    // Charging and discharging
    if (charge_request_kW > 0) {
        // charging requested
        if (charge_request_kW > maxP_kW)
            charge_request_kW = maxP_kW;
        new_charge_kWh = currentE_kWh + timestep_size_in_h*charge_request_kW;
        if (new_charge_kWh > maxE_kWh)
            new_charge_kWh = maxE_kWh;
        currentP_kW  = (new_charge_kWh - currentE_kWh)/timestep_size_in_h;
        currentE_kWh = new_charge_kWh;
    } else if (charge_request_kW < 0) {
        // discharging requested
        if (-charge_request_kW > maxP_kW)
            charge_request_kW = -maxP_kW;
        new_charge_kWh = currentE_kWh + timestep_size_in_h*charge_request_kW;
        if (new_charge_kWh < 0)
            new_charge_kWh = 0;
        currentP_kW  = (new_charge_kWh - currentE_kWh)/timestep_size_in_h;
        currentE_kWh = new_charge_kWh;
    }

    // calculate new SOC value
    SOC = currentE_kWh / maxE_kWh;
}

void ComponentBS::resetInternalState() {
    //
    // This method resets the internal SOC to the initial level
    //
    SOC = initial_SoC;
    currentE_kWh = maxE_kWh * initial_SoC;
}

