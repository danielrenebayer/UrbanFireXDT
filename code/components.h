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
#include <variant>

#include "vehicles.h"


// State structures for component save/restore functionality
/**
 * @brief State structure for PhotoVoltaic (PV) component.
 * 
 * This structure stores all internal state variables of a PV component
 * that need to be preserved during save/restore operations.
 */
struct PVComponentState {
    float currentGeneration_kW;
    double generation_cumsum_total_kWh;
    double generation_cumsum_cweek_kWh;
};
/**
 * @brief State structure for Battery Storage (BS) component.
 * 
 * This structure stores all internal state variables of a battery storage component
 * that need to be preserved during save/restore operations.
 */
struct BSComponentState {
    double SOC;
    double currentE_kWh;
    double currentP_kW;
    double currentE_from_grid_kWh;
    double currentP_from_grid_kW;
    double currentE_from_surplus_kWh;
    double cweek_E_withdrawn_kWh;
    double total_E_withdrawn_kWh;
    double cweek_E_withdrawn_from_grid_kWh;
    double total_E_withdrawn_from_grid_kWh;
    unsigned long n_ts_SOC_empty;
    unsigned long n_ts_SOC_full;
};
/**
 * @brief State structure for Heat Pump (HP) component.
 * 
 * This structure stores all internal state variables of a heat pump component
 * that need to be preserved during save/restore operations.
 */
struct HPComponentState {
    double currentDemand_kW;
    double total_consumption_kWh;
    double cweek_consumption_kWh;
};
/**
 * @brief State structure for Electric Vehicle Finite State Machine (EVFSM) component.
 * 
 * This structure stores all internal state variables of an EV component
 * that need to be preserved during save/restore operations.
 */
struct EVFSMComponentState {
    EVState current_state;
    unsigned long carID;
    unsigned long current_ts;
    double energy_demand_per_tour_ts;
    double current_P_kW;
    double sum_of_driving_distance_km;
    double sum_of_E_used_for_driving_kWh;
    double sum_of_E_charged_home_kWh;
    double sum_of_E_discharged_home_kWh;
    unsigned long sum_of_ts_EV_is_connected;
    BSComponentState battery_state; // Nested battery state
};
/**
 * @brief State structure for Charging Station (CS) component.
 * 
 * This structure stores all internal state variables of a charging station component
 * that need to be preserved during save/restore operations.
 */
struct CSComponentState {
    double current_demand_kW;
    double total_consumption_kWh;
    double cweek_consumption_kWh;
    std::map<unsigned long, EVFSMComponentState> ev_states_by_id; ///< nested EV states mapped by carID
};
/**
 * @brief Variant type for all possible component states.
 * 
 * This variant can hold the state of any component type.
 */
using ComponentStateVariant = std::variant<
    PVComponentState,
    BSComponentState,
    HPComponentState,
    CSComponentState,
    EVFSMComponentState
>;

// The following classes are defined in this header file:
class BaseComponent;
class BaseComponentSemiFlexible;
class RoofSectionPV;
class ComponentPV;
class ComponentBS;
class ComponentHP;
class ComponentCS;
class EVFSM;

#include "global.h"

/*!
 * An abstract class acting as a super-class for all existing components.
 */
class BaseComponent {
    public:
        virtual ~BaseComponent() {}
        // modifieres (on structural level of the simulation)
        virtual void resetInternalState() = 0; ///< Resets the internal simulation state (like counters), but retains the structural state like attached sub-components
        // modifiers (in the course of simulation time)
        virtual void resetWeeklyCounter() = 0; ///< Resets the internal, weekly counters
        
        // State save/restore functionality - each component implements its specific state type
        
        /**
         * @brief Saves the current internal state of the component.
         * 
         * This pure virtual method is implemented by each component to save
         * its current internal state variables into a component-specific state structure.
         * The returned variant contains the appropriate state type for the component.
         * 
         * @return ComponentStateVariant containing the component's current state
         */
        virtual ComponentStateVariant saveInternalState() const = 0;
        
        /**
         * @brief Restores the internal state of the component from a saved state.
         * 
         * This pure virtual method is implemented by each component to restore
         * its internal state variables from a previously saved state structure.
         * 
         * @param state The state variant containing the component's saved state
         */
        virtual void restoreInternalState(const ComponentStateVariant& state) = 0;
};

/*!
 * An abstract class acting as super-class for all components that can be controlled, but still have to fulfil a task.
 * An example would be a heat pump or an individual EV. Both offer flexibility to some extend, but their main purpose
 * is to generate heat or to drive. Thus their control is semi-flexible, as the main purpose still has to be fulfilled.
 * 
 * Instead, a battery storage is a fully-flexible controlable element, as it has no other task to fulfil.
 */
class BaseComponentSemiFlexible : public BaseComponent {
    public:
        BaseComponentSemiFlexible();
        ~BaseComponentSemiFlexible();
        virtual double get_currentDemand_kW() const = 0;

        /**
         * Returns the maximum electricity consumption of a component for the next n time steps (given some flexibility).
         * This means that the value at position 0 returns the maximum cummulated consumption AT THE END of the current time step.
         * This method assumes to start in the current time step, i.e. the consumption up to the current time step is considered to be 0.
         * The length of the resulting vector is set using BaseComponentSemiFlexible::set_horizon_in_ts().
         * Attention: The object the returned pointer referes to is overwritten on subsequent calls!
         * @return: Returns a pointer to a vector storing min/max demand per future time step - READ THE ATTENTION NOTICE
         */
        const std::vector<double>* get_future_max_consumption_kWh() const { return &future_maxE_storage; }

        /**
         * Returns the minimum electricity consumption of a component for the next n time steps (given some flexibility).
         * This means that the value at position 0 returns the minimum cummulated consumption AT THE END of the current time step.
         * This method assumes to start in the current time step, i.e. the consumption up to the current time step is considered to be 0.
         * The length of the resulting vector is set using BaseComponentSemiFlexible::set_horizon_in_ts().
         * Attention: The object the returned pointer referes to is overwritten on subsequent calls!
         * @return: Returns a pointer to a vector storing min/max demand per future time step - READ THE ATTENTION NOTICE
         */
        const std::vector<double>* get_future_min_consumption_kWh() const { return &future_minE_storage; }

        /**
         * Sets another horizon for the number of time steps returned by BaseComponentSemiFlexible::get_future_min_max_demand()
         */
        void set_horizon_in_ts(unsigned int new_horizon);

        /**
         * Sets the demand value to the value as given in the data.
         * If this function is called, BaseComponentSemiFlexible::setDemandToGivenValues() must not be called anymore.
         */
        virtual void setDemandToProfileData(unsigned long ts) = 0;

        /**
         * Changes the (shiftable) demand of the component for the current time step.
         * It must be called **before** the current demand and total consumption values are set correctly.
         * If this function is called, BaseComponentSemiFlexible::setDemandToProfileData() must not be called anymore.
         * @param new_demand_kW: The new power for this time step in kW.
         * @return true if no error occured or false, if the upper or lower bound was violated
         */
        virtual bool setDemandToGivenValue(double new_demand_kW) = 0;

    protected:
        // cached member variables
        /**
         * Internal storage of an object that is outputted by BaseComponentSemiFlexible::get_future_max_consumption(),
         * only required that we do not need to reserve / allocate a new vector with every call of
         * BaseComponentSemiFlexible::get_future_max_consumption()
         */
        std::vector<double> future_maxE_storage;
        /**
         * Internal storage of an object that is outputted by BaseComponentSemiFlexible::get_future_min_consumption(),
         * only required that we do not need to reserve / allocate a new vector with every call of
         * BaseComponentSemiFlexible::get_future_min_consumption()
         */
        std::vector<double> future_minE_storage;
        
};


class RoofSectionPV {
    public:
        RoofSectionPV(size_t locationID, float this_section_kWp, const std::string& orientation);
        float get_generation_at_ts_kW(unsigned long ts) const; ///< Returns the current feed in at a given time step. This method does not change the object. The validity of the arguments is NOT checked anymore!
        float get_section_kWp()    const { return this_section_kWp; }
        const std::string& get_orientation()   const { return orientation;      }
        size_t get_profile_index() const { return profile_index;    }
        // structural modifiers
        void update_section_kWp(float new_kWp); ///< This structural modifier sets a new kWp for the roof section.
    private:
        // constant member variables
        const float* profile_data; ///< Reference to the array of size Global::get_n_timesteps(), where the profile is stored. Should be a part of global::pv_profiles_data
        float this_section_kWp; // not const, because it can be changed (e.g., by the optimization when integrating PV optimization)
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
        ComponentPV(float kWp_per_m2, float min_kWp_sec, float max_kWp_sec, float max_kWp_unit, unsigned long locationID); ///< Constructor in the case of dynamic kWp computation. If max_kWp is set to a value <= 0, it will be ignored
        // getter methods
        double get_kWp() const           { return total_kWp; } ///< Returns the total peak power, summed over all sections. See also ComponentPV::get_kWp_per_section().
        float get_currentGeneration_kW() const { return currentGeneration_kW; }
        double get_cweek_generation_kWh(){ return generation_cumsum_cweek_kWh; } ///< Returns the produced energy in kWh from the start of the current week until the current time step
        double get_total_generation_kWh(){ return generation_cumsum_total_kWh; } ///< Returns the total produced energy in kWh from the start of the simulation run until the current time step
        float get_generation_at_ts_kW(unsigned long ts) const; ///< Returns the generation at a given time step in kW. If there is no data available for this time step, 0.0 is returned.
        std::string* get_section_string(const std::string& prefix_per_line); ///< Returns a string listing information about all existing, simulated roof sections - one line per section
        /**
         * Returns the future electricity generation power for the future time steps over the complete (control) time horizon.
         * Attention: The object the returned pointer referes to is overwritten on subsequent calls! Do NOT delete the returned object.
         * @return: Returns a vector storing future generation in kW per future time step - READ THE ATTENTION NOTICE
         */
        const std::vector<double>* get_future_generation_kW() const { return &future_generation_kW; }
        /**
         * Computes and returns the total electricity generation for the future, starting at time step 0 and ending with the final time step of the simulation.
         * @returns: A list of orientations, containing the vector with the length of the full simulation horizon - with the same ordering as the internal variable ComponentPV::roof_sections.
         */
        std::list<std::vector<double>> get_total_generation_by_section_kW() const;
        /**
         * Computes and returns the installed power in kWp per section (not in total).
         * See also ComponentPV::get_kWp().
         * @returns: A list of kWp per orientation with the same ordering as the internal variable ComponentPV::roof_sections.
         */
        std::list<double> get_kWp_per_section() const;
        //
        // update / action methods during simulation run
        void  calculateCurrentFeedin(unsigned long ts);
        // methods for simulation reset
        void resetWeeklyCounter();
        void resetInternalState();
        // State save/restore implementation
        using BaseComponent::saveInternalState;
        /**
         * @brief Saves the current internal state of the PV component.
         * 
         * @return ComponentStateVariant containing PVComponentState 
         */
        ComponentStateVariant saveInternalState() const override;
        using BaseComponent::restoreInternalState;
        /**
         * @brief Restores the internal state of the PV component.
         * 
         * @param state ComponentStateVariant that must contain a PVComponentState
         */
        void restoreInternalState(const ComponentStateVariant& state) override;
        void set_horizon_in_ts(unsigned int new_horizon); ///< Sets another horizon for the number of time steps returned by ComponentPV::get_future_generation_kW()
        // modification methods for structural modifications
        /**
         * This method updates the PV size per section, according to the given list. The order must be the ssame as it is in the internal variable roof_sections.
         */
        void set_kWp_per_section(const std::list<double>& new_values);
    private:
        // constant members
        std::vector<RoofSectionPV> roof_sections;
        // semi-constant member variables, i.e. they might change for parameter variations
        float total_kWp;
        // member variables that can change over time
        float currentGeneration_kW;
        std::vector<double> future_generation_kW; ///< Internal storage of the object that is outputted by ComponentPV::get_future_generation_kW(), only required that we do not need to reserve / allocate a new vector with every call
        double generation_cumsum_total_kWh; ///< Total generation over the complete simulation run
        double generation_cumsum_cweek_kWh; ///< Generation of energy starting at the beginning of the current week
};

class ComponentBS : public BaseComponent {
    friend class EVFSM;
    public:
        ComponentBS(double maxE_kWh, double maxP_kW, double E_over_P_ratio, double discharge_rate_per_step, double efficiency_in, double efficiency_out, double initial_SoC); ///< Default constructor if the battery should represent a large-scale or household-style battery system
        ComponentBS(double maxE_kWh, double discharge_rate_per_step, double efficiency_in, double initial_SoC); ///< Constructor to use if it is the battery of an EV
        // getter methods
        double get_SOC() const               { return SOC; }          ///< Returns the current state of charge of the battery
        double get_SOE() const               { return currentE_kWh; } ///< Returns the current amount of energy stored in the battery
        double get_currentCharge_kWh() const { return currentE_kWh; } ///< Returns the charge of the battery (i.e. how much kWh are stored inside)
        double get_currentLoad_kW() const    { return currentP_kW;  } ///< Returns the load of the battery (from outside perspective)
        double get_gridOnly_discharge_kW() const { return currentP_from_grid_kW; } ///< Returns the amount of energy that is currently discharged, see ComponentBS::get_currentLoad_kW() if smaller 0.0, in kW that has been stored in the BS from the grid and not from surplus PV. With values from 0.0 to maxP.
        double get_maxE_kWh()       const    { return maxE_kWh;     }
        double get_maxP_kW()        const    { return maxP_kW;      }
        double get_current_EFC()   const    { return (maxE_kWh > 0) ? total_E_withdrawn_kWh / maxE_kWh : 0.0; } ///< Returns the equivalent full cycles (EFC)
        double get_cweek_EFC()     const    { return (maxE_kWh > 0) ? cweek_E_withdrawn_kWh / maxE_kWh : 0.0; } ///< Returns the equivalent full cycles (EFC) for the current week only
        unsigned long get_n_ts_empty() const { return n_ts_SOC_empty; } ///< Returns the number of time steps where the battery is empty
        unsigned long get_n_ts_full()  const { return n_ts_SOC_full;  } ///< Returns the number of time steps where the battery is fully charged
        double get_cweek_withdrawn_E_kWh() const { return cweek_E_withdrawn_kWh; } ///< Returns the total energy that is taken from the battery from the beginning of the current week until now
        double get_total_withdrawn_E_kWh() const { return total_E_withdrawn_kWh; } ///< Returns the total energy that is taken from the battery from the beginning of the simulation run until now
        double get_cweek_withdrawn_E_gridOnly_kWh() const { return cweek_E_withdrawn_from_grid_kWh; } ///< Returns the total energy that is taken from the battery and that was initially charged from the grid (not PV!) from the beginning of the current week until now
        double get_total_withdrawn_E_gridOnly_kWh() const { return total_E_withdrawn_from_grid_kWh; } ///< Returns the total energy that is taken from the battery and that was initially charged from the grid (not PV!) from the beginning of the simulation run until now
        double get_chargeRequest() const { return charge_request_kW; } ///< Returns the current charge request in kW
        double get_currentE_from_surplus() const { return currentE_from_surplus_kWh; } ///< Returns the amount of energy that was charged from surplus PV in kWh
        // setter methods
        void  set_chargeRequest(double requested_charge_kW) { charge_request_kW = requested_charge_kW; }
        void  set_grid_charged_amount(double grid_charged_kW); ///< Sets the amount that has been charged from the grid (and not from the PV) for the current step. Must be executed AFTER ComponentBS::calculateActions() has been called. This call is optional, but at maximum it must not be called more than once, before calculateActions() is called again!
        void set_surplus_charged_amount(double surplus_charged_kW); ///< Sets the amount that has been charged from global surplus controller for the current step. Must be executed AFTER ComponentBS::calculateActions() has been called. This call is optional, but at maximum it must not be called more than once, before calculateActions() is called again! This also consideres energy that was "discharged" by the surplus controller to directly cover load.
        void reset_surplus_charged_amount(); ///< Resets the amount of energy that was charged from global surplus controller to zero. This can be called from time to time fo prevent the amount drifting too much due to numerical inaccuracies, when the surplus controller does not knkow BESS states.
    protected:
        void  set_SOE_without_computations(double new_SOE_kWh); ///< Changes the current SOE without updating other internal variables! Should only be used by class EVFSM for the pre-computation!
    public:
        void  set_maxE_kWh(double value);
        void  set_maxP_kW (double value);
        void  set_maxP_by_EPRatio(double EP_ratio);
        void  set_efficiency_in(double value);
        void  set_efficiency_out(double value);
        void  set_self_discharge_rate(double value);
        // update / action methods
        void calculateActions();
        void resetWeeklyCounter();
        void resetInternalState();
        double const validateNoSurplusChargeRequest(double charge_request_kW); ///< Validates and constrains a charge request to be within battery charge/discharge capabilities and only allows discharging from non-surplus-charged energy
        // State save/restore implementation
        using BaseComponent::saveInternalState;
        /**
         * @brief Saves the current internal state of the battery storage component.
         * 
         * @return ComponentStateVariant containing BSComponentState 
         */
        ComponentStateVariant saveInternalState() const override;
        using BaseComponent::restoreInternalState;
        /**
         * @brief Restores the internal state of the battery storage component.
         * 
         * @param state ComponentStateVariant that must contain a BSComponentState
         */
        void restoreInternalState(const ComponentStateVariant& state) override;
    private:
        // semi-constant member variables, i.e. they might change for parameter variations
        double maxE_kWh;
        double maxP_kW;
        double E_over_P_ratio;
        double discharge_rate_per_step;
        double efficiency_in;
        double efficiency_out;
        const double initial_SoC;
        // member variables that can change over time
        double SOC;
        double currentE_kWh; ///< current energy inside the battery
        double currentP_kW;  ///< current power from the perspective of outside (internally efficiency has to be added)
        double currentE_from_grid_kWh;  ///< current energy stored inside the battery that was feed-in by the grid, not by the PV surplus - this value is always greater or equal 0.0, and smaller or equal to currentE_kWh
        double currentE_from_surplus_kWh; ///< current energy stored inside the battery that was feed-in by global surplus controller - this value is always greater or equal 0.0, and smaller or equal to currentE_kWh, and larger or equal to currentE_from_grid_kWh
        double currentP_from_grid_kW;   ///< current power that was sourced not by surplus PV, but by grid charging
        double charge_request_kW;
        double cweek_E_withdrawn_kWh; ///< summation variable for EFC computation for this current week
        double total_E_withdrawn_kWh; ///< summation variable for EFC computation from the beginning of the simulation run
        double cweek_E_withdrawn_from_grid_kWh; ///< summation variable for EFC computation for this current week - only considering energy that was charged from the grid (and not the PV)
        double total_E_withdrawn_from_grid_kWh; ///< summation variable for EFC computation from the beginning of the simulation run - only considering energy that was charged from the grid (and not the PV)
        unsigned long n_ts_SOC_empty; ///< Number of time steps where battery is empty (from the beginning of the simulation run)
        unsigned long n_ts_SOC_full;  ///< Number of time steps where battery is full (from the beginning of the simulation run)
};

class ComponentHP : public BaseComponentSemiFlexible {
    public:
        ComponentHP(const ControlUnit* connected_unit, float annual_econs_kWh);
        // getter methods
        using BaseComponentSemiFlexible::get_currentDemand_kW;
        double get_currentDemand_kW() const { return currentDemand_kW; }
        /**
         * Returns the total consumption from the beginning of the simulation up to the latest executed step in kWh.
         */
        double get_total_consumption_kWh()   const { return total_consumption_kWh; }
        /**
         * Returns the total consumption from the beginning of the current week up to the latest executed step in kWh.
         */
        double get_cweek_consumption_kWh()   const { return cweek_consumption_kWh; }
        /**
         * Returns the (simulated) rated power of the heat pump **without** the AUX heating.
         */
        float get_rated_power_without_AUX()  const { return rated_power_kW; }
        //float  get_demand_at_ts_kW(unsigned long ts) const; ///< Returns the demand of the heat pump at a given time step in kW. If there is no data available for this time step, 0.0 is returned.
        using BaseComponentSemiFlexible::get_future_max_consumption_kWh;
        using BaseComponentSemiFlexible::get_future_min_consumption_kWh;
        /**
         * Returns the maximum power of this component for the next n time steps (given some flexibility).
         * The length of the resulting vector is set using BaseComponentSemiFlexible::set_horizon_in_ts().
         * Attention: The object the returned pointer referes to is overwritten on subsequent calls!
         * @return: Returns a pointer to a vector storing min/max power per future time step - READ THE ATTENTION NOTICE
         */
        const std::vector<double>* get_future_max_power_kW() const { return &future_maxP_storage; }
        /**
         * Returns the minimum power of this component for the next n time steps (given some flexibility).
         * The length of the resulting vector is set using BaseComponentSemiFlexible::set_horizon_in_ts().
         * Attention: The object the returned pointer referes to is overwritten on subsequent calls!
         * @return: Returns a pointer to a vector storing min/max power per future time step - READ THE ATTENTION NOTICE
         */
        const std::vector<double>* get_future_min_power_kW() const { return &future_minP_storage; }

        //
        // update / action methods
        using BaseComponentSemiFlexible::set_horizon_in_ts;
        void set_horizon_in_ts(unsigned int new_horizon);
        /**
         * Computes the internal state for the first / next time step with a given timestep ID as parameter.
         */
        void computeNextInternalState(unsigned long ts);
        using BaseComponentSemiFlexible::setDemandToProfileData;
        void setDemandToProfileData(unsigned long ts);
        using BaseComponentSemiFlexible::setDemandToGivenValue;
        bool setDemandToGivenValue(double new_demand_kW);
        void resetWeeklyCounter();
        void resetInternalState();
        // State save/restore implementation
        using BaseComponent::saveInternalState;
        /**
         * @brief Saves the current internal state of the heat pump component.
         * 
         * @return ComponentStateVariant containing HPComponentState 
         */
        ComponentStateVariant saveInternalState() const override;
        using BaseComponent::restoreInternalState;
        /**
         * @brief Restores the internal state of the heat pump component.
         * 
         * @param state ComponentStateVariant that must contain a HPComponentState
         */
        void restoreInternalState(const ComponentStateVariant& state) override;
        //
        // static methods for initializing the random generators
        static void InitializeRandomGenerator();
        static void VacuumStaticVariables();
    private:
        // constant member variables
        const ControlUnit* connected_unit;
        const float yearly_electricity_consumption_kWh;
        const float scaling_factor; ///< Factor to scale the profile to fit the yearly_electricity_consumption_kWh
        const float* profile_data; ///< Reference to the profile of the demand in kW per time step, should be one of global::hp_profiles
        const double* profile_cumsum; ///< Reference to the cumsum of the profile_data
        float rated_power_kW; ///< Rated power of the heat pump in kW
        // member variables that can change over time
        double currentDemand_kW;
        double total_consumption_kWh; ///< Total consumption since the beginning of the simulation period
        double cweek_consumption_kWh; ///< Total consumption since the beginning of the the currently simulated week
        /**
         * Internal storage of the object that is outputted by BaseComponentSemiFlexible::get_future_unshiftable_demand_kW(),
         * only required that we do not need to reserve / allocate a new vector with every call
         */
        std::vector<double> future_maxP_storage; ///< Storage of maximum power over the considered time horizon
        std::vector<double> future_minP_storage; ///< Storage of minimum power over the considered time horizon
#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
        // Variables for access protection of non-const methods during the simulation run
        bool state_s1;
#endif
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
class ComponentCS : public BaseComponent {
    public:
        ComponentCS(ControlUnit* calling_control_unit, unsigned int number_of_flats); ///< Constructs a new Charging Station component. There will be one charging point per flat. @param calling_control_unit: Reference to the calling unit, @param number_of_flats: The number of flats connected to the charging station
        ~ComponentCS();
        // getters
        bool is_enabled() const { return enabled; }
        double get_max_P_kW() const; ///< Returns the maximum available charging power for this station in kW (regardless if this is currently available, as some EVs might already be fully charged or not available)
        double get_currentDemand_kW() const { return current_demand_kW; } ///< Returns the current power in kW for the current time step. Only valid after the call of set_charging_value(). This value might differ from the request set using set_charging_value().
        double get_total_consumption_kWh() const { return total_consumption_kWh;  } ///< Returns the total consumed energy in kWh after the current time step. Only valid after the call of set_charging_value()
        double get_cweek_consumption_kWh() const { return cweek_consumption_kWh;  } ///< Returns the total consumed energy in kWh after the current time step. Only valid after the call of set_charging_value()
        unsigned long get_n_EVs_pc()  const; ///< Returns the number of EVs that are currently at home AND connected with the station
        unsigned long get_n_EVs_pnc() const; ///< Returns the number of EVs that are currently at home AND NOT connected with the station
        unsigned long get_n_EVs()     const; ///< Returns the number of connected EVs if the component is enabled, otherwise 0 is returned.
        unsigned long get_possible_n_EVs() const; ///< Returns the number of possible connected EVs if the component would be enabled
        unsigned long get_control_unit_id() const; ///< Returns the control unit ID of the installation place
        std::vector<const EVFSM*> get_listOfEVs() const {
            std::vector<const EVFSM*> const_list;
            const_list.reserve(listOfEVs.size());
            for (EVFSM* ev : listOfEVs) const_list.push_back(ev);
            return const_list;
        } ///< Returns a reference to the internal list of EVs
        unsigned long get_n_chargers() const { return n_chargers; } ///< Returns the number of chargers available at this charging station (typically the number of flats of the building)
        /**
         * Returns the maximum electricity consumption of this component for the next n time steps (given some flexibility).
         * This means that the value at position 0 returns the maximum cummulated consumption AT THE END of the current time step.
         * This method assumes to start in the current time step, i.e. the consumption up to the current time step is considered to be 0.
         * The length of the resulting vector is set using ComponentCS::set_horizon_in_ts().
         * Attention: The object the returned pointer referes to is overwritten on subsequent calls!
         * @return: Returns a pointer to a vector storing min/max demand per car (index 0) and per future time step (index 0) - READ THE ATTENTION NOTICE
         */
        const std::vector<const std::vector<double>*>* get_future_max_consumption_kWh() const { return &future_maxE_storage; }
        /**
         * Returns the minimum electricity consumption of a component for the next n time steps (given some flexibility).
         * This means that the value at position 0 returns the minimum cummulated consumption AT THE END of the current time step.
         * This method assumes to start in the current time step, i.e. the consumption up to the current time step is considered to be 0.
         * The length of the resulting vector is set using ComponentCS::set_horizon_in_ts().
         * Attention: The object the returned pointer referes to is overwritten on subsequent calls!
         * @return: Returns a pointer to a vector storing min/max demand per car (index 0) and per future time step (index 0) - READ THE ATTENTION NOTICE
         */
        const std::vector<const std::vector<double>*>* get_future_min_consumption_kWh() const { return &future_minE_storage; }
        /**
         * Returns the maximum power per time step and EV for the future horizon.
         * Attention: The object the returned pointer referes to is overwritten on subsequent calls!
         * @return: Returns a pointer to a vector storing min/max demand per car (index 0) and per future time step (index 0) - READ THE ATTENTION NOTICE
         */
        const std::vector<const std::vector<double>*>* get_future_max_power_kW() const { return &future_maxP_storage; }
        // modifieres (on structural level of the simulation)
        void enable_station();
        void disable_station();
        void resetWeeklyCounter();
        void resetInternalState();
        // State save/restore implementation
        using BaseComponent::saveInternalState;
        /**
         * @brief Saves the current internal state of the charging station component.
         * 
         * @return ComponentStateVariant containing CSComponentState 
         */
        ComponentStateVariant saveInternalState() const override;
        using BaseComponent::restoreInternalState;
        /**
         * @brief Restores the internal state of the charging station component.
         * 
         * @param state ComponentStateVariant that must contain a CSComponentState
         */
        void restoreInternalState(const ComponentStateVariant& state) override;
        void add_ev(unsigned long carID);
        void set_horizon_in_ts(unsigned int new_horizon); ///< Sets another horizon for the number of time steps returned by ComponentCS::get_future_min_consumption_kWh() and ComponentCS::get_future_max_consumption_kWh()
        void preprocess_ev_data(); ///< Preprocesses EV data, calls EVFSM::preprocessTourInformation() for all attached EVs. Call this method only once right before the first simulation run.
        // modifiers (in the course of simulation time)
        void setCarStatesForTimeStep(unsigned long ts); ///< Sets the car states for all attached EVs for a new time step 'ts'. This method must be called with strictly consecutive values ​​of parameter 'ts'.
        void setDemandToProfileData(unsigned long ts); ///< Sets the charging demand to the (immediate charging) profile value. If this method is called once during a simulation run, ComponentCS::setDemandToGivenValues() MUST NOT be called afterwards.
        bool setDemandToGivenValues(std::vector<double>& charging_power_per_EV_kW); ///< Sets the charging demand (in kW) per EV as given by any contoller. If this method is called once during a simulation run, ComponentCS::setDemandToProfileData() MUST NOT be called afterwards. Returns true if no error occured for all connected EVs or false otherwise.
    private:
        // constant members
        const ControlUnit* installation_place;
        const unsigned long n_chargers; ///< The number of chargers available at this charging station
        double max_charging_power; // Not const, as the value might not be known on pre-initialization
        // variable members, constant during one simulation run
        bool enabled;
        std::vector<EVFSM*> listOfEVs;
        // variable members, variable during simulation run
        double current_demand_kW;
        double total_consumption_kWh;
        double cweek_consumption_kWh; ///< Total demand since the beginning of the the currently simulated week
        std::vector<const std::vector<double>*> future_maxE_storage; ///< Internal storage of future max energy consumption per EV and time step
        std::vector<const std::vector<double>*> future_minE_storage;
        std::vector<const std::vector<double>*> future_maxP_storage;
#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
        // Variables for access protection of non-const methods during the simulation run
        bool is_callable_setCarStatesForTimeStep;
        bool is_callable_set_charging_value;
#endif
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
class EVFSM : public BaseComponentSemiFlexible {

    public:
        EVFSM(unsigned long carID, ComponentCS* homeStation);
        ~EVFSM();
        // getters
        unsigned long get_carID() const { return carID; }
        EVState get_current_state() const { return current_state; }
        // EVStateIfConnAtHome get_current_state_icah() const { return current_state_icah; } ///< Returns the current state of the EV iff it is connected at home
        using BaseComponentSemiFlexible::get_future_max_consumption_kWh;
        using BaseComponentSemiFlexible::get_future_min_consumption_kWh;
        using BaseComponentSemiFlexible::get_currentDemand_kW;
        const std::vector<double>* get_future_max_power_kW() const { return &future_maxP_storage; } ///< Returns the maximum power of this EV per time step in the controller horizon. Attention: Returned object will be overwritten after calling EVFSM::setCarStateForTimeStep().
        double get_currentDemand_kW() const { return current_P_kW; } ///< Gets the current charging power in kW; Only valid after calling EVFSM::setDemandToProfileData() or EVFSM::setDemandToGivenValue()
        std::string* get_metrics_string_annual(); ///< Returns some metrics as string (useful for the output). Header see EVFSM::MetricsStringHeaderAnnual. Call this function only if simulation run is finished!
        const ComponentBS* get_battery() const { return battery; } ///< Returns a reference to the battery component of the EV
        // modifiers (on structural level of the simulation)
        using BaseComponentSemiFlexible::set_horizon_in_ts;
        void set_horizon_in_ts(unsigned int new_horizon);
        void add_weekly_tour(short weekday, unsigned long departure_ts_of_day, unsigned long ts_duration, double tour_length_km, bool with_work); ///< This method adds a home-centered car tour to the current car. All parameters that represent a time must have the same alignment as the global time information.
        void preprocessTourInformation(); ///< Transforms the list of weekly tours into a list of tours for the complete simulation time span and computes upper and lower cumlulative energy required (i.e., fill variables BaseComponentSemiFlexible::future_maxE_storage and BaseComponentSemiFlexible::future_minE_storage from upper class). THis method MUST be called before the main simulation run starts, but it MUST be called after EVFSM::add_weekly_tour() is called for the last time, i.e., all tours have been added.
        void resetInternalState(); ///< Resets the internal state
        // State save/restore implementation
        using BaseComponentSemiFlexible::saveInternalState;
        /**
         * @brief Saves the current internal state of the EV finite state machine component.
         * 
         * @return ComponentStateVariant containing EVFSMComponentState 
         */
        ComponentStateVariant saveInternalState() const override;
        using BaseComponentSemiFlexible::restoreInternalState;
        /**
         * @brief Restores the internal state of the EV finite state machine component.
         * 
         * @param state ComponentStateVariant that must contain an EVFSMComponentState
         */
        void restoreInternalState(const ComponentStateVariant& state) override;
        // modifiers (in the course of simulation time)
        void setCarStateForTimeStep(unsigned long ts); ///< Sets the car state for a new time step 'ts'. This method must be called with strictly consecutive values ​​of parameter 'ts'. It uses the precomputed vectors for internal processing.
        using BaseComponentSemiFlexible::setDemandToProfileData;
        void setDemandToProfileData(unsigned long ts);
        using BaseComponentSemiFlexible::setDemandToGivenValue;
        bool setDemandToGivenValue(double new_demand_kW);
        using BaseComponentSemiFlexible::resetWeeklyCounter;
        void resetWeeklyCounter() {} ///< Do nothing. The EV has no weekly metrics.
        // static methods
        static const std::map<unsigned long, EVFSM*>& GetArrayOfInstances() { return list_of_cars; } ///< Returns the map of all existing instances. The objects itself are mutable, but the map reference is const.
        static unsigned long GetNumberOfEVs() { return list_of_cars.size(); }
        static void AddWeeklyTour(unsigned long carID, short weekday, unsigned long departure_ts_of_day, unsigned long ts_duration, double tour_length_km, bool with_work); ///< This class method adds a home-centered car tour to the car with ID carID. All parameters that represent a time must have the same alignment as the global time information.
        static void VacuumStaticVariables();
        static void SetSeed(unsigned int seed); ///< Sets the seed for the EVFSM-class random number generator

    private:
        // constant members
        const unsigned long carID;       ///< ID of the car. The ID does not necessarily be successive but it must be unique.
        const float econs_kWh_per_km;    ///< Energy consumption of the current EV in kWh/km
        ComponentCS const* homeStation;  ///< Reference to the home station of the EV
        // variable members (constant during a simulation run)
        std::vector<std::vector<WeeklyVehicleTour>*> list_of_tours_pd; ///< Vector of vector of tours, one individual vector per week day (0 -> monday, 7 -> sunday)
        std::vector<WeeklyVehicleTour*> list_of_all_tours; ///< Vector of all weekly tours. Is the same as `list_of_tours_pd` in a flattened order, i.e., witout day information.
        std::vector<EVState> prec_vec_of_states; ///< precomputed list of states per time step
        std::vector<double> prec_vec_of_minE; ///< Complete precomputed vector of minimum charged electricity up to the end of time step ts
        std::vector<double> prec_vec_of_maxE; ///< Complete precomputed vector of maximum charged electricity up to the end of time step ts
        std::vector<double> prec_vec_of_driving_distance_km; ///< Complete precomputed vector of driven distance of the EV in km up to the end of time step ts
        std::vector<float>  prec_vec_of_maxP_kW; ///< Complete precomputed vector of maximum charging power per time step ts
        std::vector<double> prec_vec_of_curr_BS_E_cons_kWh; ///< Complete precomputed vector of electricity taken from the battery for driving (at time step ts, or 0, if the car is parked at home in any state)
        // cached member variables
        std::vector<double> future_maxP_storage;  ///< Internal storage of maximum charging power over the future time steps
        // attached members
        ComponentBS* battery;
        // variable members, variable during a simulation run
        EVState current_state;           ///< Internal current state of the EV
        unsigned long current_ts;        ///< Internal variable storing the current time step
        //EVStateIfConnAtHome current_state_icah; ///< Internal current state of the EV iff it is connected at home
        double energy_demand_per_tour_ts; ///< The mean energy demand per tour time step. This is the demand of the total tour divided by the number of time steps of the tour -> We assume a linear decay of the battery SOC, ignoring stops
        double current_P_kW;  ///< The current charging power in kW
        // variables for the final metrics calculation
        double sum_of_driving_distance_km;    ///< Sum of driven distance in km (only updated at the end of a tour, when the home place is reached again)
        double sum_of_E_used_for_driving_kWh; ///< Sum of electricity consumed by the EV required for driving in kWh
        double sum_of_E_charged_home_kWh;     ///< Sum of charged electricity in kWh
        double sum_of_E_discharged_home_kWh;  ///< Sum of discharged electricity in kWh from the EV
        ulong  sum_of_ts_EV_is_connected;     ///< Number of time steps the EV is connected as at the home charging point
#ifdef ADD_METHOD_ACCESS_PROTECTION_VARS
        // Variables for access protection of non-const methods during the simulation run
        bool state_s1; // directly after initialization, tours can be added
        bool state_s2; // ready for new/first time step, i.e., all tour information has been processed, no tours can be added anymore
        bool state_s3; // EV internal state set to new time step, current demand can be set
#endif
        //
        // class members
        static std::map<unsigned long, EVFSM*> list_of_cars;
        static std::default_random_engine*            random_generator; ///< Generator required for random sampling
        static std::uniform_real_distribution<float>* distribution    ; ///< Required for random sampling
    public:
        static const std::string MetricsStringHeaderAnnual; ///< The header for the output string produced by `EVFSM::get_metrics_string_annual()`
};

#endif
