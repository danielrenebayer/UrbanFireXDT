/*
 * components.h
 *
 * It contains all individual components
 *
 */

#ifndef COMPONENTS_H
#define COMPONENTS_H

#include <map>
#include <random>
#include <string>
#include <vector>

// The following classes are defined in this header file:
class RoofSectionPV;
class ComponentPV;
class ComponentBS;
class ComponentHP;
class ComponentCS;

#include "global.h"
#include "vehicles.h"


/*!
 * An abstract class acting as a super-class for all existing components.
 */
class BaseComponent {
    public:
		virtual ~BaseComponent() {}
        // modifieres (on structural level of the simulation)
        virtual void resetInternalState() = 0; ///< Resets the internal simulation state (like counters), but retains the structural state like attached sub-components
        // modifiers (in the course of simulation time)
        //
};

/*!
 * An abstract class acting as super-class for all components that can be controlled, but still have to fulfil a task.
 * An example would be a heat pump or a charging station. Both offer flexibility to some extend, but their main purpose
 * is to generate heat or charge a EV. Thus their control is semi-flexible, as the main purpose still has to be fulfilled.
 * 
 * Instead, a battery storage is a fully-flexible controlable element, as it has no other task to fulfil.
 */
class BaseComponentSemiFlexible : public BaseComponent {
    public:
		virtual ~BaseComponentSemiFlexible() {}
        virtual float get_currentDemand_kW() const = 0;
};


class RoofSectionPV {
    public:
        RoofSectionPV(float this_section_kWp, std::string& orientation);
        float get_currentFeedin_kW(unsigned long ts);
        float get_section_kWp()    const { return this_section_kWp; }
        const std::string& get_orientation()   const { return orientation;      }
        size_t get_profile_index() const { return profile_index;    }
    private:
        // constant member variables
        const float* profile_data; ///< Reference to the array of size Global::get_n_timesteps(), where the profile is stored. Should be a part of global::pv_profiles_data
        const float  this_section_kWp;
        const std::string orientation;
        size_t profile_index;
        //
        // static data for selecting the next time series for expansion
        static std::map<std::string, size_t> next_pv_idx; ///< Next index per orientation (given as string)
        // static data for selecting the next time series randomly
        static std::map<std::string, std::default_random_engine>            random_generators;
        static std::map<std::string, std::uniform_int_distribution<size_t>> distributions;
};

class ComponentPV : public BaseComponent {
    public:
        ComponentPV(float kWp_static, unsigned long locationID); ///< Constructor in the case of static kWp computation
        ComponentPV(float kWp_per_m2, float min_kWp, float max_kWp, unsigned long locationID); ///< Constructor in the case of dynamic kWp computation. If max_kWp is set to a value <= 0, it will be ignored
        // getter methods
        float get_kWp() const            { return total_kWp; }
        float get_currentGeneration_kW() { return currentGeneration_kW; }
        double get_total_generation_kWh(){ return total_generation_kWh; } ///< Returns the total produced energy in kWh from the start of the simulation run until the current time step
        std::string* get_section_string(const std::string& prefix_per_line); ///< Returns a string listing information about all existing, simulated roof sections - one line per section
        // update / action methods during simulation run
        void  calculateCurrentFeedin(unsigned long ts);
        // methods for simulation reset
        void resetInternalState();
        // modification methods for structural modifications
        //void  set_kWp(float value);
    private:
        // constant members
        std::vector<RoofSectionPV> roof_sections;
        // semi-constant member variables, i.e. they might change for parameter variations
        float total_kWp;
        // member variables that can change over time
        float currentGeneration_kW;
        double total_generation_kWh; ///< Total generation over the complete simulation run
};

class ComponentBS : public BaseComponent {
    public:
        ComponentBS(float maxE_kWh, float maxP_kW, float E_over_P_ratio, float discharge_rate_per_step, float efficiency_in, float efficiency_out, float initial_SoC); ///< Default constructor if the battery should represent a large-scale or household-style battery system
        ComponentBS(float maxE_kWh, float discharge_rate_per_step, float efficiency_in, float initial_SoC); ///< Constructor to use if it is the battery of an EV
        // getter methods
        float get_SOC() const               { return SOC; }          ///< Returns the current state of charge of the battery
        float get_currentCharge_kWh() const { return currentE_kWh; } ///< Returns the charge of the battery (i.e. how much kWh are stored inside)
        float get_currentLoad_kW() const    { return currentP_kW;  } ///< Returns the load of the battery (from outside perspective)
        float get_maxE_kWh()       const    { return maxE_kWh;     }
        float get_maxP_kW()        const    { return maxP_kW;      }
        double get_current_EFC()   const    { return (maxE_kWh > 0) ? total_E_withdrawn_kWh / maxE_kWh : 0.0; } ///< Returns the equivalent full cycles (EFC)
        unsigned long get_n_ts_empty() const { return n_ts_SOC_empty; } ///< Returns the number of time steps where the battery is empty
        unsigned long get_n_ts_full()  const { return n_ts_SOC_full;  } ///< Returns the number of time steps where the battery is fully charged
        double get_total_withdrawn_E_kWh() const { return total_E_withdrawn_kWh; } ///< Returns the total energy that is taken from the battery from the beginning of the simulation run until now
        // setter methods
        void  set_chargeRequest(float requested_charge_kW) { charge_request_kW = requested_charge_kW; }
        void  set_maxE_kWh(float value);
        void  set_maxP_kW (float value);
        void  set_maxP_by_EPRatio(float EP_ratio);
        // update / action methods
        void calculateActions();
        void resetInternalState();
    private:
        // semi-constant member variables, i.e. they might change for parameter variations
        float maxE_kWh;
        float maxP_kW;
        float E_over_P_ratio;
        // constant member variables (other languages might call this 'final')
        const float discharge_rate_per_step;
        const float efficiency_in;
        const float efficiency_out;
        const float initial_SoC;
        // member variables that can change over time
        float SOC;
        float currentE_kWh; ///< current energy inside the battery
        float currentP_kW;  ///< current power from the perspective of outside (internally efficiency has to be added)
        float charge_request_kW;
        double total_E_withdrawn_kWh; ///< summation variable for EFC computation
        unsigned long n_ts_SOC_empty; ///< Number of time steps where battery is empty
        unsigned long n_ts_SOC_full;  ///< Number of time steps where battery is full
};

class ComponentHP : public BaseComponentSemiFlexible {
    public:
        ComponentHP(float yearly_econs_kWh);
        // getter methods
        using BaseComponentSemiFlexible::get_currentDemand_kW;
        float  get_currentDemand_kW() const { return currentDemand_kW; }
        double get_total_demand_kWh() const { return total_demand_kWh; }
        // update / action methods
        void calculateCurrentFeedin(unsigned long ts);
        void resetInternalState();
        //
        // static methods for initializing the random generators
        static void InitializeRandomGenerator();
        static void VacuumStaticVariables();
    private:
        // constant member variables
        const float yearly_electricity_consumption_kWh;
        const float scaling_factor; ///< Factor to scale the profile to fit the yearly_electricity_consumption_kWh
        const float* profile_data; ///< Reference to the profile, should be one of global::hp_profiles
        // member variables that can change over time
        float currentDemand_kW;
        double total_demand_kWh; ///< Total demand since the beginning of the simulation period
        //
        // static data for selecting the next time series for expansion
        static size_t next_hp_idx;
        static bool random_generator_init;
        static std::default_random_engine* random_generator;
        static std::uniform_int_distribution<size_t>* distribution;
};

/**
 * This class represents a charging station.
 * Every control unit holds a charging station by default, which is disabled.
 * Enabling a station in the simulation represents to build such a station in
 * the reality.
 */
class ComponentCS : public BaseComponentSemiFlexible {
    public:
        ComponentCS();
        ~ComponentCS();
        // getters
        bool is_enabled() const { return enabled; }
        float  get_max_P_kW() const; ///< Returns the maximum available charging power for this station in kW
        using BaseComponentSemiFlexible::get_currentDemand_kW;
        float  get_currentDemand_kW() const { return current_demand_kW; } ///< Returns the current power in kW for the current time step. Only valid after the call of set_charging_value(). This value might differ from the request set using set_charging_value().
        double get_total_demand_kWh() const { return total_demand_kWh;  } ///< Returns the total consumed energy in kWh after the current time step. Only valid after the call of set_charging_value()
        //double get_min_curr_charging_power_kW() const; ///< Returns the minimal charging power at the current time step. The charging station requires at least this portion to fulfil all charging demands.
        float  get_max_curr_charging_power_kW() const; ///< The maximal charging power that could be charged into the currently parking cars at the station.
        unsigned long get_n_EVs_pc()  const; ///< Returns the number of EVs that are currently at home AND connected with the station
        unsigned long get_n_EVs_pnc() const; ///< Returns the number of EVs that are currently at home AND NOT connected with the station
        unsigned long get_n_EVs()     const; ///< Returns the number of connected EVs if the component is enabled, otherwise 0 is returned.
        unsigned long get_possible_n_EVs() const; ///< Returns the number of possible connected EVs if the component would be enabled
        // modifieres (on structural level of the simulation)
        void enable_station();
        void disable_station();
        void resetInternalState();
        void add_ev(unsigned long carID);
        // modifiers (in the course of simulation time)
        void setCarStatesForTimeStep(unsigned long ts, unsigned int dayOfWeek_l, unsigned int hourOfDay_l);
        void set_charging_value(float power_kW); ///< Sets the charging value in kW that should be charged into the connected EVs
    private:
        // constant members
        const float max_charging_power;
        // variable members, constant during one simulation run
        bool enabled;
        std::vector<EVFSM*> listOfEVs;
        // variable members, variable during simulation run
        float  current_demand_kW;
        double total_demand_kWh;
        float  charging_power_required_kW;
        float  charging_power_possible_kW;
        std::list<EVFSM*> charging_order_req; ///< List of cars sorted after their charging urgency for dispatching the charging process; This list contains all cars that require charging
        std::list<EVFSM*> charging_order_pos; ///< List of cars sorted after their charging urgency for dispatching the charging process; This list contains all cars that can be charged
};

#endif
