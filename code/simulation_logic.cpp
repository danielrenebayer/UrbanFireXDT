
#include "simulation_logic.h"
using namespace simulation;

#include <ctime>
#include <iostream>
#include <sstream>

#include "helper.h"
#include "global.h"
#include "output.h"
#include "units.h"


bool simulation::runSimulationForOneParamSetting() {
    //
    // This function loops over all time steps as they are defined in the data.
    // If multiple settings (e.g. a parameter variation) should be applied,
    // one has to call this function more often.
    //
    std::cout << "Run simulation for a complete episode ..." << std::endl;

    int n_tsteps = Global::get_n_timesteps();
    struct tm* tm_start = Global::get_ts_start_tm();
    struct tm* tm_end   = Global::get_ts_end_tm();
    //time_t t_start = mktime(tm_start);
    //time_t t_end   = mktime(tm_end);
    struct tm* current_tm;
    bool sim_started = false; // gets true, if simulation range (as given by tm_start) has been reached
    // main loop
    for (int ts = 1; ts <= n_tsteps; ts++) {
        // get current time as struct tm
        current_tm = global::time_localtime_str->at(ts - 1);
        // jump time steps if they are not inside the simulation range
        if (sim_started) {
            //if (difftime(t_end, mktime(current_tm)) <= 0) {
            if (compare_struct_tm(current_tm, tm_end) >= 0) {
                std::cout << "End of the simulation range (as defined in the scenario) has been reached." << std::endl;
                break;
            }
        } else {
            //if (difftime(mktime(current_tm), t_start) >= 0) {
            if (compare_struct_tm(current_tm, tm_start) >= 0) {
                sim_started = true;
                std::cout << "Start of the simulation range (as defined in the scenario) is reached." << std::endl;
            } else {
                continue;
            }
        }
        //
        // execute one step
        if (!oneStep(ts)) return false;
        // flush output buffers every configurable step, so that RAM consumption does not increase too much
        if ((ts % global::n_ts_between_flushs) == 0)
            output::flushBuffers();

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

bool simulation::runSimulationForAllVariations(int scenario_id) {
    //
    // This function runs the simulation for all given parameter
    // variations. If no variation is selected, it executes the
    // simulation once only.
    // Returns false, if an error occurs during the simulation,
    // else true.
    //

    if (Global::is_parameter_variation()) {
        // in case of a parameter variation:
        // loop over all parameter combinations
        unsigned int param_vari_combi_ind = 0;
        for (auto& list_of_variables : *(global::parameter_var_list)) {
            //
            // 0. set global variable
            global::curr_param_vari_combi_index = param_vari_combi_ind;
            // 0.1. initialize struct required for outputting the current parameter setting
            output::CurrentParamValues cParamVals;
            // 0.2. reset internal variables for control units
            ControlUnit::ResetAllInternalStates();
            //
            // 1. set current variable values
            //    i.e. iterate over all variable-name / value combinations and
            //    set values accordingly
            ControlUnit*const* cuList = ControlUnit::GetArrayOfInstances();
            for (auto& var_name_and_val : list_of_variables) {
                if        (var_name_and_val.first.compare("expansion PV kWp") == 0) {
                    for (int i = 0; i < ControlUnit::GetNumberOfInstances(); i++)
                        cuList[i]->set_exp_pv_kWp(var_name_and_val.second);
                    cParamVals.exp_pv_kWp = var_name_and_val.second;
                    cParamVals.exp_pv_kWp_set = true;

                } else if (var_name_and_val.first.compare("expansion BS P in kW")  == 0) {
                    for (int i = 0; i < ControlUnit::GetNumberOfInstances(); i++)
                        cuList[i]->set_exp_bs_maxP_kW(var_name_and_val.second);
                    cParamVals.exp_bs_maxP_kW = var_name_and_val.second;
                    cParamVals.exp_bs_maxP_kW_set = true;

                } else if (var_name_and_val.first.compare("expansion BS E in kWh") == 0) {
                    for (int i = 0; i < ControlUnit::GetNumberOfInstances(); i++)
                        cuList[i]->set_exp_bs_maxE_kWh(var_name_and_val.second);
                    cParamVals.exp_bs_maxE_kWh = var_name_and_val.second;
                    cParamVals.exp_bs_maxE_kWh_set = true;

                } else if (var_name_and_val.first.compare("expansion BS initial SOC") == 0) {
                    std::cerr << "This is not implemented!" << std::endl;

                } else {
                    std::cerr << "Unknown parameter variable to vary: " << var_name_and_val.first << std::endl;
                    return false;
                }
            }
            //
            // 2. open output files
            output::initializeDirectoriesPerPVar(scenario_id);
            output::initializeSubstationOutput(scenario_id);
            output::initializeCUOutput(scenario_id);
            // 2.b output the current parameter variation combination
            output::outputCurrentParamVariCombi(cParamVals);
            //
            // 3. run the simulation
            bool no_error = runSimulationForOneParamSetting();
            if (!no_error) return false;
            //
            // 4a. output metrics (if selected)
            output::outputMetrics();
            // 4b. close output files
            output::closeOutputs();
            //
            // 5. increment combination counter
            param_vari_combi_ind++;
        }
    } else {
        //
        // 1. open output files
        output::initializeDirectoriesPerPVar(scenario_id);
        output::initializeSubstationOutput(scenario_id);
        output::initializeCUOutput(scenario_id);
        // 1.b output the current parameter variation combination
        output::CurrentParamValues cParamVals;
        output::outputCurrentParamVariCombi(cParamVals);
        //
        // 2. Run the simulation
        bool no_error = runSimulationForOneParamSetting();
        if (!no_error) return false;
        //
        // 3a. output metrics (if selected)
        output::outputMetrics();
        // 3b. close output files
        output::closeOutputs();
    }

    return true;
}
