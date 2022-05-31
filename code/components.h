/*
 * components.h
 *
 * It contains all individual components
 *
 */

#include "global.h"

#ifndef COMPONENTS_H
#define COMPONENTS_H

class ComponentPV {
    public:
        ComponentPV(float kWp);
        // getter methods
        float get_kWp() const            { return kWp; }
        float get_currentGeneration_kW() { return currentGeneration_kW; }
        // update / action methods
        void  calculateCurrentFeedin(int ts);
        void  set_kWp(float value)       { kWp = value; }
    private:
        // semi-constant member variables, i.e. they might change for parameter variations
        float kWp;
        // member variables that can change over time
        float currentGeneration_kW;
        // TODO: Ausrichtung beachten -> dann current Feedin ueber globalstrahlung und sonnenstand ausrechnen
};

class ComponentBS {
    public:
        ComponentBS(float maxE_kWh, float maxP_kW, float discharge_rate_per_step, float efficiency, float initial_SoC);
        // getter methods
        float get_SOC() const               { return SOC; }
        float get_currentCharge_kWh() const { return currentE_kWh; }
        float get_currentLoad_kW() const    { return currentP_kW;  }
        const float get_maxE_kWh() const    { return maxE_kWh;     }
        const float get_maxP_kW()  const    { return maxP_kW;      }
        // setter methods
        void  set_chargeRequest(float requested_charge_kW) { charge_request_kW = requested_charge_kW; }
        void  set_maxE_kWh(float value)     { maxE_kWh = value; }
        void  set_maxP_kW (float value)     { maxP_kW  = value; }
        // update / action methods
        void calculateActions();
        void resetInternalState();
    private:
        // semi-constant member variables, i.e. they might change for parameter variations
        float maxE_kWh;
        float maxP_kW;
        // constant member variables (other languages might call this 'final')
        const float discharge_rate_per_step;
        const float efficiency;
        const float initial_SoC;
        // member variables that can change over time
        float SOC;
        float currentE_kWh;
        float currentP_kW;
        float charge_request_kW;
};

class ComponentHP {
    // TODO: Implement Heat Pump
};

class ComponentWB {
    // TODO: Implement Wallboxes
};

#endif
