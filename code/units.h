/*
 * units.h
 *
 * It contains the definition of all units. A unit is
 * a compositions of individual components, that might
 * has more complex control rules or the additionl data
 * input.
 *
 */

#include <string>
#include <list>
#include <map>

#include "global.h"
#include "components.h"
#include "output.h"

#ifndef UNITS_H
#define UNITS_H


class Substation;
class ControlUnit;
class MeasurementUnit;


class Substation {
    /*
        This class represents a substation.
    */
    public:
        static Substation* InstantiateNewSubstation(unsigned long id, std::string* name) {
            return new Substation(id, name);
        }
    private:
        Substation(unsigned long id, std::string* name);
    public:
        ~Substation();
        // getter methods
        unsigned long get_id() { return id; }
        const std::string * get_name() { return name; }
        // modifiers bevore simulation start
        void add_unit(ControlUnit* unit);
        // methods for simulation run
        float calc_load();
        const std::list<ControlUnit*>* get_connected_units() { return connected_units; } ///< Return the list of connected units (as read only list)
        //
        // static functions
        // 1. Initializers and destructors
        static void InitializeStaticVariables(unsigned long n_substations);
        static void VacuumInstancesAndStaticVariables();
        // 2. getter functions
        static inline Substation* GetInstance(unsigned long id);
        static        Substation*const * GetArrayOfInstances() {return st__substation_list;}
        static        size_t GetNumberOfInstances() {return st__n_substations;}
    private:
        // constant member variables (other languages might call this 'final')
        const unsigned long id;
        const std::string *const name;
        // member variables that can change over time
        std::list<ControlUnit*>* connected_units;
        //
        // static list of substations
        static bool st__substation_list_init;
        static unsigned long st__n_substations;
        static unsigned long st__new_substation_position;
        static Substation** st__substation_list;
};


class ControlUnit {
    /*
        This class represents a control unit
        (i.e. a private house or a small company)
    */
    public:
        static ControlUnit* InstantiateNewControlUnit(unsigned long unitID, unsigned long substation_id, unsigned long locationID) {
            return new ControlUnit(unitID, substation_id, locationID);
        }
    private:
        ControlUnit(unsigned long unitID, unsigned long substation_id, unsigned long locationID);
    public:
        ~ControlUnit();
        void add_unit(MeasurementUnit* unit);
        bool has_electricity_demand();
        bool has_pv();
        bool has_bs();
        bool has_hp();
        bool has_wb();
        bool has_chp();
        int  get_exp_combi_bit_repr();
        int  get_exp_combi_bit_repr_from_MUs();
        int  get_exp_combi_bit_repr_sim_added();
        float get_current_load_vSMeter_kW() { return current_load_vSM_kW; }
        size_t get_n_MUs()     { return connected_units->size(); } // returns the number of MUs, that are connected to the given control unit
        size_t get_unitID()    { return unitID; } // returns the unit ID of this control unit
        float get_sim_comp_pv_kWp(); // returns the kWp of the PV-component that is added for the simulation, returns 0 if there is no added PV component
        float get_sim_comp_bs_P_kW(); // returns the power in kW of the battery storage component that is added for the simulation, returns 0 if there is no added battery
        float get_sim_comp_bs_E_kWh(); // returns the capacity in kWh of the battery storage component that is added for the simulation, returns 0 if there is no added battery
        double  get_SSR(); ///< Returns the SSR of the CU from the start of the simulation run until the time of function call; most usefull at the end of a simulation run
        double  get_SCR(); ///< Returns the SCR of the CU from the start of the simulation run until the time of function call; most usefull at the end of a simulation run
        double  get_NPV(); ///< Returns the net present value (NPV) of the CU from the start of the simulation run until the time of function call; most usefull at the end of a simulation run
        string* get_metrics_string(); // call this function only if simulation run is finished! It will the compute sums of flows,SSC,SSR and output this as a string
        string* get_pv_section_string(); // This function returns a string containing information about the sections of the sim. added PV component. If no PV component is added, it returns an empty string.
        // modifiers
        void add_exp_pv();
        void add_exp_bs();
        void add_exp_hp();
        void add_exp_wb();
        void set_output_object(CUOutput* output_obj);
        void set_exp_pv_params_A(float kWp_static); ///< Set the kWp of expanded PV installations in the case of static kWp computation per section
        void set_exp_pv_params_B(float kWp_per_m2, float min_kWp, float max_kWp); ///< Set the kWp of expanded PV installations in the case of dynamic kWp computation per section
        void set_exp_bs_maxE_kWh(float value);
        void set_exp_bs_maxP_kW (float value);
        void set_exp_bs_E_P_ratio(float value); //< Set the E:P-ratio for simulatively added BS components to @param value
        void remove_sim_added_components(); ///< Remove all components that are added simulatively
        // for simulation runs
        bool compute_next_value(unsigned long ts);
        //
        // static functions
        // 1. Initializers and destructors
        static void InitializeStaticVariables(unsigned long n_CUs);
        static void VacuumInstancesAndStaticVariables();
        // 2. getter functions
        static inline ControlUnit* GetInstance(unsigned long unitID);
        static ControlUnit*const * GetArrayOfInstances() {return st__cu_list;}
        static size_t GetNumberOfInstances() {return st__n_CUs;}
        // 3. modifiers for all created objects
        static void ResetAllInternalStates();
        static void RemoveAllSimAddedComponents(); ///< Removes all simulatively added components from all control units
    private:
        // constant member variables (other languages might call this 'final')
        const unsigned long unitID;
        Substation *const higher_level_subst;
        const unsigned long locationID;
        // member variables that can change over time
        std::list<MeasurementUnit*>* connected_units;
        bool has_sim_pv; ///< boolean variable that states if a PV installation is simulatively added
        bool has_sim_bs; ///< boolean variable that states if a battery storage is simulatively added
        bool has_sim_hp; ///< boolean variable that states if a heap pump is simulatively added
        bool has_sim_wb; ///< boolean variable that states if a wallbox is simulatively added
        ComponentPV* sim_comp_pv; ///< Reference to the simulated PV-Component (if it exists)
        ComponentBS* sim_comp_bs; ///< Reference to the simulated battery storage component (if it exists)
        ComponentHP* sim_comp_hp; ///< Reference to the simulated Heat Pump Component (if it exists)
        ComponentWB* sim_comp_wb; ///< Reference to the simulated Wallbox Component (if it exists)
        CUOutput*    output_obj;
        double sum_of_consumption_kWh;    ///< The sum of consumed energy in kWh starting from the beginning of the current simulation run
        double sum_of_self_cons_kWh;      ///< The sum of self-consumed energy in kWh starting from the beginning of the current simulation run
        double sum_of_PV_generated_kWh;   ///< The sum of PV-generated energy in kWh starting from the beginning of the current simulation run
        double sum_of_feed_into_grid_kWh; ///< The sum of energy in kWh, that is fed into the grid, starting from the beginning of the current simulation run
        //float* history_self_prod_load_kW; ///< Array for later analysis, holding the historical values of the self-produced load, that is directly consumed
        //float* history_pv_generation_kW; ///< Array for later analysis, holding the historical values of the PV-generation
        //float* history_avg_consumption_load_kW; ///< Array for later analysis, holding the historical (average) consumption load for the corresponding time step
        //bool   create_history_output; ///< True, if a history output should be created for this control unit.
        //
        float current_load_vSM_kW; ///< Current load at the virtual smart meter
        float self_produced_load_kW; ///< Load [in kW] that is produced by the PV / taken from Battery / El. vehicle AND directly consumed by the measurement units
        //
        // static list of CUs
        static bool st__cu_list_init;
        static unsigned long st__n_CUs;
        static unsigned long st__new_CU_position;
        static ControlUnit** st__cu_list;
};


class MeasurementUnit {
    /*
        The measurement unit represents an existing unit, which is measured by an existing
        smart meter where we have real measured data for.
        It can have additional components, like PV, if it is known that these additional
        components exists in the real unit as well.
    */
    public:
        static MeasurementUnit* InstantiateNewMeasurementUnit(size_t meUID, size_t unitID, std::string * meterPointName, size_t locID, 
                bool has_demand, bool has_feedin, bool has_pv_resid, bool has_pv_opens,
                bool has_bess,   bool has_hp,     bool has_wb,       bool has_chp) {
                    return new MeasurementUnit(meUID, unitID, meterPointName, locID, has_demand, has_feedin, has_pv_resid, has_pv_opens, has_bess, has_hp, has_wb, has_chp);
                }
    private:
        // initialization and destruction
        MeasurementUnit(size_t meUID, size_t unitID, std::string * meterPointName, size_t locID, 
                        bool has_demand, bool has_feedin, bool has_pv_resid, bool has_pv_opens,
                        bool has_bess,   bool has_hp,     bool has_wb,       bool has_chp);
    public:
        ~MeasurementUnit();
        bool load_data(const char * filepath);
        bool has_demand() { return rsm_has_demand; }
        bool has_feedin() { return rsm_has_feedin; }
        bool has_pv()     { return rsm_with_pv_residential || rsm_with_pv_open_space; }
        bool has_pv_residential() { return rsm_with_pv_residential; }
        bool has_pv_open_space()  { return rsm_with_pv_open_space; }
        bool has_bs()     { return rsm_with_bess;}
        bool has_hp()     { return rsm_with_hp;  }
        bool has_wb()     { return rsm_with_wb;  }
        bool has_chp()    { return rsm_with_chp; }
        int  get_expansion_combination() { return expansion_combination; }
        inline const std::string * get_meterPointName() const;
        inline unsigned long get_meUID() const;
        inline unsigned long get_locationID() const;
        // for simulation runs
        bool compute_next_value(unsigned long ts);
        float get_current_ts_rsm_value() { return current_load_rsm_kW; }
        //
        // Class (i.e. static) functions
        // 1. Initializers and destructors
        static void InitializeStaticVariables(unsigned long n_MUs);
        static void VacuumInstancesAndStaticVariables();
        static size_t GetNumberOfInstances() {return st__n_MUs;}
        //
        // 2. Class getter methods
    private:
        // constant member variables (other languages might call this 'final')
        const unsigned long meUID;
        ControlUnit *const higher_level_cu;
        const std::string *const meterPointName;
        const unsigned long locationID;
        // member variables that can change over time
        float current_load_rsm_kW;
        bool rsm_has_demand; ///< RSM stands for Real Smart Meter
        bool rsm_has_feedin;
        bool rsm_with_pv_residential;
        bool rsm_with_pv_open_space;
        bool rsm_with_bess;
        bool rsm_with_hp;
        bool rsm_with_wb;
        bool rsm_with_chp;
        int expansion_combination;
        // member variables storing the data
        bool   data_loaded;
        int*   data_timestepID;
        float* data_value_demand;
        float* data_value_feedin;
        //char*  data_status_demand;
        //char*  data_status_feedin;
        //
        // static list of MUs
        static bool st__mu_list_init;
        static size_t st__n_MUs;
        static size_t st__new_MU_position;
        static MeasurementUnit** st__mu_list;
};



enum OpenSpacePVOrWindType {
    PV,
    Wind
};

class OpenSpacePVOrWind {
    public:
        OpenSpacePVOrWind(float kWp, OpenSpacePVOrWindType type);
        float get_current_feedin_kW() { return current_feedin_kW; }
        // for simulation runs
        bool compute_next_value(unsigned long ts);
    private:
        // constant member variables
        const float kWp;
        const float* profile_data; ///< Reference to the array of size Global::get_n_timesteps(), where the profile is stored. Should be global::wind_profile or global::pv_profile
        // members that might change during the time
        float current_feedin_kW;
};

#endif

