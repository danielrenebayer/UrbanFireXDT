/*
 * components.h
 *
 * It contains all individual components
 *
 */

#include "global.h"

#ifndef COMPONENTS_H
#define COMPONENTS_H

#include <random>
#include <string>
#include <vector>

class RoofSectionPV {
    public:
        RoofSectionPV(float this_section_kWp, std::string& orientation);
        float get_currentFeedin_kW(int ts);
        const float        get_section_kWp()   const { return this_section_kWp; }
        const std::string& get_orientation()   const { return orientation;      }
        const size_t       get_profile_index() const { return profile_index;    }
        //
        // static methods for initializing the random generators
        static void InitializeRandomGenerator();
        static void VacuumStaticVariables();
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
        static bool random_generator_init;
        static std::map<std::string, std::default_random_engine>            random_generators;
        static std::map<std::string, std::uniform_int_distribution<size_t>> distributions;
};

class ComponentPV {
    public:
        ComponentPV(float kWp_static, unsigned long locationID); ///< Constructor in the case of static kWp computation
        ComponentPV(float kWp_per_m2, float min_kWp, float max_kWp, unsigned long locationID); ///< Constructor in the case of dynamic kWp computation
        // getter methods
        float get_kWp() const            { return total_kWp; }
        float get_currentGeneration_kW() { return currentGeneration_kW; }
        std::string* get_section_string(const std::string& prefix_per_line); ///< Returns a string listing information about all existing, simulated roof sections - one line per section
        // update / action methods
        void  calculateCurrentFeedin(int ts);
        //void  set_kWp(float value);
    private:
        // constant members
        std::vector<RoofSectionPV> roof_sections;
        // semi-constant member variables, i.e. they might change for parameter variations
        float total_kWp;
        // member variables that can change over time
        float currentGeneration_kW;
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
    public:
        ComponentHP(float yearly_econs_kWh);
        // getter methods
        float get_currentDemand_kW() { return currentDemand_kW; }
        // update / action methods
        void calculateCurrentFeedin(int ts);
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
        //
        // static data for selecting the next time series for expansion
        static size_t next_hp_idx;
        static bool random_generator_init;
        static std::default_random_engine* random_generator;
        static std::uniform_int_distribution<size_t>* distribution;
};

class ComponentWB {
    // TODO: Implement Wallboxes
};

#endif
