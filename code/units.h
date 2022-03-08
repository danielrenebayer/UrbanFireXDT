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
        void add_unit(ControlUnit* unit);
        //
        //
        static void PreInitializeStaticVariables();
        static void InitializeStaticVariables(int n_substations);
        static void VacuumStaticVariables();
    private:
        // constant member variables (other languages might call this 'final')
        const int id;
        const std::string* name;
        // member variables that can change over time
        std::list<ControlUnit*>* connected_units;
        //
        // static list of substations
        static bool st__substation_list_init;
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
        //
        //
        static void PreInitializeStaticVariables();
        static void InitializeStaticVariables(int n_CUs);
        static void VacuumStaticVariables();
    private:
        // constant member variables (other languages might call this 'final')
        const int unitID;
        const Substation* higher_level_subst;
        // member variables that can change over time
        std::list<MeasurementUnit*>* connected_units;
        //
        // static list of CUs
        static bool st__cu_list_init;
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
        MeasurementUnit(int meloID, std::string * melo, int locID);
        ~MeasurementUnit();
        bool load_data(const char * filepath);
        //
        //
        static void PreInitializeStaticVariables();
        static void InitializeStaticVariables(int n_MUs);
        static void VacuumStaticVariables();
        //
        // getter methods
        inline const std::string * get_melo() const;
        inline const int get_meloID() const;
        inline const int get_locationID() const;
        inline bool has_feedin();
        inline bool has_demand();
    private:
        // constant member variables (other languages might call this 'final')
        const int locationID;
        const int meloID;
        const std::string * melo;
        const ControlUnit* higher_level_cu;
        // member variables that can change over time
        float current_load_rsm_kW;
        bool rsm_has_demand;
        bool rsm_has_feedin;
        bool rsm_with_pv;
        bool rsm_with_bess;
        bool rsm_with_hp;
        bool rsm_with_wb;
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
        static MeasurementUnit** st__mu_list;
};


#endif