
#include "simulation_logic.h"
using namespace simulation;

#include <ctime>
#include <iostream>
#include <sstream>

#include "global.h"
#include "output.h"
#include "units.h"


bool simulation::runSimulation() {
    //
    // This function loops over all time steps as they are defined in the data
    //
    std::cout << "Run simulation for a complete episode ..." << std::endl;

    int n_tsteps = Global::get_n_timesteps();
    struct tm* tm_start = Global::get_ts_start_tm();
    struct tm* tm_end   = Global::get_ts_end_tm();
    time_t t_start = mktime(tm_start);
    time_t t_end   = mktime(tm_end);
    struct tm* current_tm;
    bool sim_started = false; // gets true, if simulation range (as given by tm_start) has been reached
    // main loop
    for (int ts = 1; ts <= n_tsteps; ts++) {
        // get current time as struct tm
        current_tm = global::time_localtime_str->at(ts - 1);
        // jump time steps if they are not inside the simulation range
        if (sim_started) {
            if (difftime(t_end, mktime(current_tm)) <= 0) {
                std::cout << "End of the simulation range (as defined in the scenario) has been reached." << std::endl;
                break;
            }
        } else {
            if (difftime(mktime(current_tm), t_start) >= 0) {
                sim_started = true;
                std::cout << "Start of the simulation range (as defined in the scenario) is reached." << std::endl;
            } else {
                continue;
            }
        }
        //
        // execute one step
        if (!oneStep(ts)) return false;
        /*
        // flush output buffers every 64 steps, so that RAM consumption does not increase too much
        if ((ts % 64) == 0)
            output::flushBuffers();
        */
    }

    std::cout << " ... run finished." << std::endl;
    return true;

}

bool simulation::oneStep(int ts) {
    //
    // Run one time step of the simulation
    // Return false if an error occurs during execution.
    // Parameter ts is the current time step
    //

    //int tsID = ts - 1;

    // loop over all control units:
    //    set new values and execute next actions
    ControlUnit*const* cuList = ControlUnit::GetArrayOfInstances();
    const int nCUs = ControlUnit::GetNumberOfInstances();
	for (int i = 0; i < nCUs; i++) {
        if (!cuList[i]->compute_next_value(ts))
            return false;
	}
    //
    // set new global radiation values for PV and wind speed
    float total_load = 0.0;
    global::unit_open_space_pv->compute_next_value(ts);
    global::unit_open_space_wind->compute_next_value(ts);
    total_load -= global::unit_open_space_pv->get_current_feedin_kW();
    total_load -= global::unit_open_space_wind->get_current_feedin_kW();
    //
    // loop over all substations: compute new load values
    // and calculate total grid load
    Substation*const* subList = Substation::GetArrayOfInstances();
    const int nSubst = Substation::GetNumberOfInstances();
    *(output::substation_output) << ts << ","; // add timestep to output
	for (int i = 0; i < nSubst; i++) {
        float current_station_load = subList[i]->calc_load();
		total_load += current_station_load;
        // stuff for output
        *(output::substation_output) << current_station_load << ",";
	}
    *(output::substation_output) << global::unit_open_space_pv->get_current_feedin_kW()   << ",";
    *(output::substation_output) << global::unit_open_space_wind->get_current_feedin_kW() << ",";
    *(output::substation_output) << total_load << "\n"; // add total load to output

    std::cout << ".";

    return true;

}
