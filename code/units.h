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

#include "global.h"

#ifndef UNITS_H
#define UNITS_H


class ControlUnit {

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
};


#endif