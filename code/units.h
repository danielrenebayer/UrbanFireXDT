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
        Substation(int id, std::string* name);
        ~Substation();
        // getter methods
        const int get_id() { return id; }
        const std::string *const get_name() { return name; }
        // modifiers bevore simulation start
        void add_unit(ControlUnit* unit);
        // methods for simulation run
        float calc_load();
        //
        // static functions
        // 1. Initializers and destructors
        static void InitializeStaticVariables(int n_substations);
        static void VacuumInstancesAndStaticVariables();
        // 2. getter functions
        static inline Substation* GetInstance(int id);
        static        Substation*const * GetArrayOfInstances() {return st__substation_list;}
        static        const int GetNumberOfInstances() {return st__n_substations;}
    private:
        // constant member variables (other languages might call this 'final')
        const int id;
        const std::string *const name;
        // member variables that can change over time
        std::list<ControlUnit*>* connected_units;
        //
        // static list of substations
        static bool st__substation_list_init;
        static int st__n_substations;
        static int st__new_substation_position;
        static Substation** st__substation_list;
};


class ControlUnit {
    /*
        This class represents a control unit
        (i.e. a private house or a small company)
    */
    public:
        ControlUnit(int unitID, int substation_id);
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
        const int get_unitID() { return unitID; } // returns the unit ID of this control unit
        float get_sim_comp_pv_kWp(); // returns the kWp of the PV-component that is added for the simulation, returns 0 if there is no added PV component
        float get_sim_comp_bs_P_kW(); // returns the power in kW of the battery storage component that is added for the simulation, returns 0 if there is no added battery
        float get_sim_comp_bs_E_kWh(); // returns the capacity in kWh of the battery storage component that is added for the simulation, returns 0 if there is no added battery
        // modifiers
        void add_exp_pv();
        void add_exp_bs();
        void add_exp_hp();
        void add_exp_wb();
        void set_output_object(CUOutput* output_obj);
        // for simulation runs
        bool compute_next_value(int ts);
        //
        // static functions
        // 1. Initializers and destructors
        static void InitializeStaticVariables(int n_CUs);
        static void VacuumInstancesAndStaticVariables();
        // 2. getter functions
        static inline ControlUnit* GetInstance(int unitID);
        static inline ControlUnit*const * GetArrayOfInstances() {return st__cu_list;}
        static inline const int GetNumberOfInstances() {return st__n_CUs;}
    private:
        // constant member variables (other languages might call this 'final')
        const int unitID;
        Substation *const higher_level_subst;
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
        //
        float current_load_vSM_kW; ///< Current load at the virtual smart meter
        float self_produced_load_kW; ///< Load [in kW] that is produced by the PV / taken from Battery / El. vehicle AND directly consumed by the measurement units
        //
        // static list of CUs
        static bool st__cu_list_init;
        static int st__n_CUs;
        static int st__new_CU_position;
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
        // initialization and destruction
        MeasurementUnit(int meloID, int unitID, std::string * melo, int locID, 
                        bool has_demand, bool has_feedin, bool has_pv_resid, bool has_pv_opens,
                        bool has_bess,   bool has_hp,     bool has_wb,       bool has_chp);
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
        inline const std::string * get_melo() const;
        inline const int get_meloID() const;
        inline const int get_locationID() const;
        // for simulation runs
        bool compute_next_value(int ts);
        float get_current_ts_rsm_value() { return current_load_rsm_kW; }
        //
        // Class (i.e. static) functions
        // 1. Initializers and destructors
        static void InitializeStaticVariables(int n_MUs);
        static void VacuumInstancesAndStaticVariables();
        //
        // 2. Class getter methods
    private:
        // constant member variables (other languages might call this 'final')
        const int meloID;
        ControlUnit *const higher_level_cu;
        const std::string *const melo;
        const int locationID;
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
        static int st__n_MUs;
        static int st__new_MU_position;
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
        bool compute_next_value(int ts);
    private:
        // constant member variables
        const float kWp;
        const float* profile_data; ///< Reference to the array of size Global::get_n_timesteps(), where the profile is stored. Should be global::wind_profile or global::pv_profile
        // members that might change during the time
        float current_feedin_kW;
};

#endif

