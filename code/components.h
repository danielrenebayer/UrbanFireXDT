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
class BaseComponent;
class BaseComponentSemiFlexible;
class RoofSectionPV;
class ComponentPV;
class ComponentBS;
class ComponentHP;
class ComponentCS;
class EVFSM;

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
        virtual void resetWeeklyCounter() = 0; ///< Resets the internal, weekly counters
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
        virtual float get_currentDemand_kW() const = 0;

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
         */
        virtual void setDemandToGivenValue(float new_demand_kW) = 0;

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
        RoofSectionPV(float this_section_kWp, std::string& orientation);
        float get_generation_at_ts_kW(unsigned long ts) const; ///< Returns the current feed in at a given time step. This method does not change the object. The validity of the arguments is NOT checked anymore!
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
        double get_cweek_generation_kWh(){ return generation_cumsum_cweek_kWh; } ///< Returns the produced energy in kWh from the start of the current week until the current time step
        double get_total_generation_kWh(){ return generation_cumsum_total_kWh; } ///< Returns the total produced energy in kWh from the start of the simulation run until the current time step
        float get_generation_at_ts_kW(unsigned long ts) const; ///< Returns the generation at a given time step in kW. If there is no data available for this time step, 0.0 is returned.
        std::string* get_section_string(const std::string& prefix_per_line); ///< Returns a string listing information about all existing, simulated roof sections - one line per section
        /**
         * Returns the future electricity generation power for the future time steps over the complete time horizon.
         * Attention: The object the returned pointer referes to is overwritten on subsequent calls! Do NOT delete the returned object.
         * @return: Returns a vector storing future generation in kW per future time step - READ THE ATTENTION NOTICE
         */
        const std::vector<double>* get_future_generation_kW() const { return &future_generation_kW; }
        //
        // update / action methods during simulation run
        void  calculateCurrentFeedin(unsigned long ts);
        // methods for simulation reset
        void resetWeeklyCounter();
        void resetInternalState();
        void set_horizon_in_ts(unsigned int new_horizon); ///< Sets another horizon for the number of time steps returned by ComponentPV::get_future_generation_kW()
        // modification methods for structural modifications
        //void  set_kWp(float value);
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
        ComponentBS(float maxE_kWh, float maxP_kW, float E_over_P_ratio, float discharge_rate_per_step, float efficiency_in, float efficiency_out, float initial_SoC); ///< Default constructor if the battery should represent a large-scale or household-style battery system
        ComponentBS(float maxE_kWh, float discharge_rate_per_step, float efficiency_in, float initial_SoC); ///< Constructor to use if it is the battery of an EV
        // getter methods
        float get_SOC() const               { return SOC; }          ///< Returns the current state of charge of the battery
        float get_SOE() const               { return currentE_kWh; } ///< Returns the current amount of energy stored in the battery
        float get_currentCharge_kWh() const { return currentE_kWh; } ///< Returns the charge of the battery (i.e. how much kWh are stored inside)
        float get_currentLoad_kW() const    { return currentP_kW;  } ///< Returns the load of the battery (from outside perspective)
        float get_maxE_kWh()       const    { return maxE_kWh;     }
        float get_maxP_kW()        const    { return maxP_kW;      }
        double get_current_EFC()   const    { return (maxE_kWh > 0) ? total_E_withdrawn_kWh / maxE_kWh : 0.0; } ///< Returns the equivalent full cycles (EFC)
        double get_cweek_EFC()     const    { return (maxE_kWh > 0) ? cweek_E_withdrawn_kWh / maxE_kWh : 0.0; } ///< Returns the equivalent full cycles (EFC) for the current week only
        unsigned long get_n_ts_empty() const { return n_ts_SOC_empty; } ///< Returns the number of time steps where the battery is empty
        unsigned long get_n_ts_full()  const { return n_ts_SOC_full;  } ///< Returns the number of time steps where the battery is fully charged
        double get_cweek_withdrawn_E_kWh() const { return cweek_E_withdrawn_kWh; } ///< Returns the total energy that is taken from the battery from the beginning of the current week until now
        double get_total_withdrawn_E_kWh() const { return total_E_withdrawn_kWh; } ///< Returns the total energy that is taken from the battery from the beginning of the simulation run until now
        // setter methods
        void  set_chargeRequest(float requested_charge_kW) { charge_request_kW = requested_charge_kW; }
    protected:
        void  set_SOE_without_computations(float new_SOE_kWh); ///< Changes the current SOE without updating other internal variables! Should only be used by class EVFSM for the pre-computation!
    public:
        void  set_maxE_kWh(float value);
        void  set_maxP_kW (float value);
        void  set_maxP_by_EPRatio(float EP_ratio);
        // update / action methods
        void calculateActions();
        void resetWeeklyCounter();
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
        double cweek_E_withdrawn_kWh; ///< summation variable for EFC computation for this current week
        double total_E_withdrawn_kWh; ///< summation variable for EFC computation from the beginning of the simulation run
        unsigned long n_ts_SOC_empty; ///< Number of time steps where battery is empty (from the beginning of the simulation run)
        unsigned long n_ts_SOC_full;  ///< Number of time steps where battery is full (from the beginning of the simulation run)
};

class ComponentHP : public BaseComponentSemiFlexible {
    public:
        ComponentHP(float yearly_econs_kWh);
        // getter methods
        using BaseComponentSemiFlexible::get_currentDemand_kW;
        float  get_currentDemand_kW() const { return currentDemand_kW; }
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
         * Returns the future electricity demand per time step that is not shiftable for the heat pump.
         * Attention: The object the returned pointer referes to is overwritten on subsequent calls! Do NOT delete the returned object.
         * @return: Returns a vector storing min/max demand per future time step - READ THE ATTENTION NOTICE
         */
        const std::vector<double>* get_future_unshiftable_demand_kW() const { return &future_unshiftable_storage; }

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
        void setDemandToGivenValue(float new_demand_kW);
        void resetWeeklyCounter();
        void resetInternalState();
        //
        // static methods for initializing the random generators
        static void InitializeRandomGenerator();
        static void VacuumStaticVariables();
    private:
        // constant member variables
        const float yearly_electricity_consumption_kWh;
        const float scaling_factor; ///< Factor to scale the profile to fit the yearly_electricity_consumption_kWh
        const float* profile_data_shiftable; ///< Reference to the profile for the     shiftable part of the demand, should be one of global::hp_profiles_shiftable
        const float* profile_data_not_shift; ///< Reference to the profile for the not-shiftable part of the demand, should be one of global::hp_profiles_not_shift
        const double* profile_shiftable_cumsum; ///< Reference to the cumsum of the profile_data_shiftable
        float rated_power_kW; ///< Rated power of the heat pump in kW
        // member variables that can change over time
        float currentDemand_kW;
        double total_consumption_kWh; ///< Total consumption since the beginning of the simulation period
        double cweek_consumption_kWh; ///< Total consumption since the beginning of the the currently simulated week
        /**
         * Internal storage of the object that is outputted by BaseComponentSemiFlexible::get_future_unshiftable_demand_kW(),
         * only required that we do not need to reserve / allocate a new vector with every call
         */
        std::vector<double> future_unshiftable_storage;
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
        float  get_max_P_kW() const; ///< Returns the maximum available charging power for this station in kW (regardless if this is currently available, as some EVs might already be fully charged or not available)
        float  get_currentDemand_kW() const { return current_demand_kW; } ///< Returns the current power in kW for the current time step. Only valid after the call of set_charging_value(). This value might differ from the request set using set_charging_value().
        double get_total_consumption_kWh() const { return total_consumption_kWh;  } ///< Returns the total consumed energy in kWh after the current time step. Only valid after the call of set_charging_value()
        double get_cweek_consumption_kWh() const { return cweek_consumption_kWh;  } ///< Returns the total consumed energy in kWh after the current time step. Only valid after the call of set_charging_value()
        //double get_min_curr_charging_power_kW() const; ///< Returns the minimal charging power at the current time step. The charging station requires at least this portion to fulfil all charging demands.
        //float  get_max_curr_charging_power_kW() const; ///< The maximal charging power that could be charged into the currently parking cars at the station.
        unsigned long get_n_EVs_pc()  const; ///< Returns the number of EVs that are currently at home AND connected with the station
        unsigned long get_n_EVs_pnc() const; ///< Returns the number of EVs that are currently at home AND NOT connected with the station
        unsigned long get_n_EVs()     const; ///< Returns the number of connected EVs if the component is enabled, otherwise 0 is returned.
        unsigned long get_possible_n_EVs() const; ///< Returns the number of possible connected EVs if the component would be enabled
        unsigned long get_control_unit_id() const; ///< Returns the control unit ID of the installation place
        /**
         * Returns the maximum electricity consumption of this component for the next n time steps (given some flexibility).
         * This means that the value at position 0 returns the maximum cummulated consumption AT THE END of the current time step.
         * This method assumes to start in the current time step, i.e. the consumption up to the current time step is considered to be 0.
         * The length of the resulting vector is set using ComponentCS::set_horizon_in_ts().
         * Attention: The object the returned pointer referes to is overwritten on subsequent calls!
         * @return: Returns a pointer to a vector storing min/max demand per future time step (index 0) and per car (index 1) - READ THE ATTENTION NOTICE
         */
        const std::vector<std::vector<double>>* get_future_max_consumption_kWh() const { return &future_maxE_storage; }
        /**
         * Returns the minimum electricity consumption of a component for the next n time steps (given some flexibility).
         * This means that the value at position 0 returns the minimum cummulated consumption AT THE END of the current time step.
         * This method assumes to start in the current time step, i.e. the consumption up to the current time step is considered to be 0.
         * The length of the resulting vector is set using ComponentCS::set_horizon_in_ts().
         * Attention: The object the returned pointer referes to is overwritten on subsequent calls!
         * @return: Returns a pointer to a vector storing min/max demand per future time step (index 0) and per car (index 1) - READ THE ATTENTION NOTICE
         */
        const std::vector<std::vector<double>>* get_future_min_consumption_kWh() const { return &future_minE_storage; }
        // modifieres (on structural level of the simulation)
        void enable_station();
        void disable_station();
        void resetWeeklyCounter();
        void resetInternalState();
        void add_ev(unsigned long carID);
        void set_horizon_in_ts(unsigned int new_horizon); ///< Sets another horizon for the number of time steps returned by ComponentCS::get_future_min_consumption_kWh() and ComponentCS::get_future_max_consumption_kWh()
        void preprocess_ev_data(); ///< Preprocesses EV data, calls EVFSM::preprocessTourInformation() for all attached EVs. Call this method only once right before the first simulation run.
        // modifiers (in the course of simulation time)
        void setCarStatesForTimeStep(unsigned long ts); ///< Sets the car states for all attached EVs for a new time step 'ts'. This method must be called with strictly consecutive values ​​of parameter 'ts'.
        void setDemandToProfileData(unsigned long ts); ///< Sets the charging demand to the (immediate charging) profile value. If this method is called once during a simulation run, ComponentCS::setDemandToGivenValues() MUST NOT be called afterwards.
        void setDemandToGivenValues(std::vector<float>& charging_power_per_EV_kW); ///< Sets the charging demand (in kW) per EV as given by any contoller. If this method is called once during a simulation run, ComponentCS::setDemandToProfileData() MUST NOT be called afterwards.
    private:
        // constant members
        const ControlUnit* installation_place;
        float max_charging_power; // Not const, as the value might not be known on pre-initialization
        // variable members, constant during one simulation run
        bool enabled;
        std::vector<EVFSM*> listOfEVs;
        // variable members, variable during simulation run
        float  current_demand_kW;
        double total_consumption_kWh;
        double cweek_consumption_kWh; ///< Total demand since the beginning of the the currently simulated week
        //std::list<EVFSM*> charging_order_req; ///< List of cars sorted after their charging urgency for dispatching the charging process; This list contains all cars that require charging
        //std::list<EVFSM*> charging_order_pos; ///< List of cars sorted after their charging urgency for dispatching the charging process; This list contains all cars that can be charged
        std::vector<std::vector<double>> future_maxE_storage;
        std::vector<std::vector<double>> future_minE_storage;
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
        EVState get_current_state() const { return current_state; }
        // EVStateIfConnAtHome get_current_state_icah() const { return current_state_icah; } ///< Returns the current state of the EV iff it is connected at home
        using BaseComponentSemiFlexible::get_future_max_consumption_kWh;
        using BaseComponentSemiFlexible::get_future_min_consumption_kWh;
        using BaseComponentSemiFlexible::get_currentDemand_kW;
        float get_currentDemand_kW() const { return current_P_kW; } ///< Gets the current charging power in kW; Only valid after calling EVFSM::setDemandToProfileData() or EVFSM::setDemandToGivenValue()
        // TODO get_min_SOC_per_time_step ...
        std::string* get_metrics_string_annual(); ///< Returns some metrics as string (useful for the output). Header see EVFSM::MetricsStringHeaderAnnual. Call this function only if simulation run is finished!
        // modifiers (on structural level of the simulation)
        using BaseComponentSemiFlexible::set_horizon_in_ts;
        void add_weekly_tour(short weekday, unsigned long departure_ts_of_day, unsigned long ts_duration, double tour_length_km, bool with_work); ///< This method adds a home-centered car tour to the current car. All parameters that represent a time must have the same alignment as the global time information.
        void preprocessTourInformation(); ///< Transforms the list of weekly tours into a list of tours for the complete simulation time span and computes upper and lower cumlulative energy required (i.e., fill variables BaseComponentSemiFlexible::future_maxE_storage and BaseComponentSemiFlexible::future_minE_storage from upper class). THis method MUST be called before the main simulation run starts, but it MUST be called after EVFSM::add_weekly_tour() is called for the last time, i.e., all tours have been added.
        void resetInternalState(); ///< Resets the internal state
        // modifiers (in the course of simulation time)
        void setCarStateForTimeStep(unsigned long ts); ///< Sets the car state for a new time step 'ts'. This method must be called with strictly consecutive values ​​of parameter 'ts'. It uses the precomputed vectors for internal processing.
        //void setInternalStateToSimRunStart(); ///< Sets the state of the EV to the beginning of a simulation run
        //void set_current_charging_power(float power_kW); ///< Sets the current charging power in kW. If this method is called, EVFSM::setDemandToProfileData() cannot be called.
        //void setDemandToProfileData(unsigned long ts); ///< Sets the energy consumption to the value for immediate charging. If this method is called, EVFSM::set_current_charging_power() cannot be called.
        using BaseComponentSemiFlexible::setDemandToProfileData;
        void setDemandToProfileData(unsigned long ts);
        using BaseComponentSemiFlexible::setDemandToGivenValue;
        void setDemandToGivenValue(float new_demand_kW);
        using BaseComponentSemiFlexible::resetWeeklyCounter;
        void resetWeeklyCounter() {} ///< TODO IMPLEMENT !
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
        std::vector<float>  prec_vec_of_curr_BS_E_cons_kWh; ///< Complete precomputed vector of electricity taken from the battery for driving (at time step ts, or 0, if the car is parked at home in any state)
        // attached members
        ComponentBS* battery;
        // variable members, variable during a simulation run
        EVState current_state;           ///< Internal current state of the EV
        //EVStateIfConnAtHome current_state_icah; ///< Internal current state of the EV iff it is connected at home
        float energy_demand_per_tour_ts; ///< The mean energy demand per tour time step. This is the demand of the total tour divided by the number of time steps of the tour -> We assume a linear decay of the battery SOC, ignoring stops
        float current_P_kW;  ///< The current charging power in kW
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
