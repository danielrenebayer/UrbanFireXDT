
#include "simulation_logic.h"
using namespace simulation;

#include <algorithm>
#include <ctime>
#include <iostream>
#include <sstream>
#include <vector>

#include "helper.h"
#include "global.h"
#include "output.h"
#include "sac_planning.h"
#include "units.h"

//
// Internal variables required for output control
unsigned int output_counter = 0;
// #define OUTPUT_STATUS_OUTPUT_FREQ 24*3
#define OUTPUT_STATUS_OUTPUT_FREQ 24*24


bool simulation::runSimulationForOneParamSetting(std::vector<ControlUnit*>* subsection /* = NULL */, const char* output_prefix /* = "" */) {
    //
    // This function loops over all time steps as they are defined in the data.
    // If multiple settings (e.g. a parameter variation) should be applied,
    // one has to call this function more often.
    //
    if (output_prefix[0] == '\0')
        std::cout << output_prefix << "Run simulation for a complete episode ..." << std::endl;

    // get variable values that remain equal during this simulation run
    double totalBatteryCapacity_kWh = ControlUnit::GetAllSimCompBatteriesCapacity_kWh();

    // get time data
    unsigned long n_tsteps = Global::get_n_timesteps();
    struct tm* tm_start = Global::get_ts_start_tm();
    struct tm* tm_end   = Global::get_ts_end_tm();
    //time_t t_start = mktime(tm_start);
    //time_t t_end   = mktime(tm_end);
    struct tm* current_tm;
    unsigned int last_step_weekday = 0; // required to identify new weeks
    unsigned long week_number = 1;
    //
    bool sim_started = false; // gets true, if simulation range (as given by tm_start) has been reached
    // main loop
    for (unsigned long ts = 1; ts <= n_tsteps; ts++) {
        // get current time as struct tm
        current_tm = global::time_localtime_str->at(ts - 1);
        // jump time steps if they are not inside the simulation range
        if (sim_started) {
            //if (difftime(t_end, mktime(current_tm)) <= 0) {
            if (compare_struct_tm(current_tm, tm_end) >= 0) {
                if (output_prefix[0] == '\0')
                    std::cout << output_prefix << "\rEnd of the simulation range (as defined in the scenario) has been reached." << std::endl;
                break;
            }
        } else {
            //if (difftime(mktime(current_tm), t_start) >= 0) {
            if (compare_struct_tm(current_tm, tm_start) >= 0) {
                sim_started = true;
                if (output_prefix[0] == '\0')
                    std::cout << output_prefix << "Start of the simulation range (as defined in the scenario) is reached." << std::endl;
            } else {
                continue;
            }
        }
        //
        // day and time handling
        unsigned int dayOfWeek_r = (unsigned int) ((current_tm->tm_wday + 6) % 7); // get day of week in the format 0->Monday, 6->Sunday
        unsigned int hourOfDay_r = (unsigned int) (current_tm->tm_hour);
        // convert from rigth-alignment to left-alignment
        unsigned int dayOfWeek_l = dayOfWeek_r;
        unsigned int hourOfDay_l = hourOfDay_r - 1;
        if (hourOfDay_r <= 0) {
            hourOfDay_l = 23;
            if (dayOfWeek_r >= 1)
                dayOfWeek_l = dayOfWeek_r - 1;
            else
                dayOfWeek_l = 6;
        }
        // has a new week started?
        if (last_step_weekday > dayOfWeek_l) {
            // only execute this during the main run (checked by subsection == NULL)
            if (Global::get_compute_weekly_metrics() && subsection == NULL) {
                std::list<std::string*> *output_list = new std::list<std::string*>();
                // iterate over all control units, get the weekly metrics string (and reset weekly summation variables automatically)
                ControlUnit*const* cuList = ControlUnit::GetArrayOfInstances();
                const size_t nCUs = ControlUnit::GetNumberOfInstances();
                for (size_t i = 0; i < nCUs; i++) {
                    output_list->push_back(
                        cuList[i]->get_metrics_string_weekly_wr(week_number)
                    );
                }
                // output results
                output::outputWeeklyMetricsStrList(output_list, week_number); // attention: this function deletes output_list
            }
            week_number++;
        }
        last_step_weekday = dayOfWeek_l;
        //
        // execute one step
        if (!oneStep(ts, dayOfWeek_l, hourOfDay_l, totalBatteryCapacity_kWh, output_prefix, subsection)) return false;
        // flush output buffers every configurable step, so that RAM consumption does not increase too much
        if ((ts % global::n_ts_between_flushs) == 0)
            output::flushBuffers();

    }

    if (output_prefix[0] == '\0') {
        std::cout << output_prefix << "... run finished." << "\n";
        std::cout << global::output_section_delimiter << std::endl;
    }
    return true;

}

bool simulation::oneStep(unsigned long ts, unsigned int dayOfWeek_l, unsigned int hourOfDay_l,
                        double totalBatteryCapacity_kWh, const char* output_prefix /* = "" */,
                        std::vector<ControlUnit*>* subsection /* = NULL */) {
    //
    // Run one time step of the simulation
    // Return false if an error occurs during execution.
    // Parameter ts is the current time step
    //

    //int tsID = ts - 1;

    // TODO: parallelization of the unit calls
    // loop over all control units:
    //    set new values and execute next actions
    if (subsection == NULL) {
        ControlUnit*const* cuList = ControlUnit::GetArrayOfInstances();
        const size_t nCUs = ControlUnit::GetNumberOfInstances();
        for (size_t i = 0; i < nCUs; i++) {
            if (!cuList[i]->compute_next_value(ts, dayOfWeek_l, hourOfDay_l))
                return false;
        }
    } else {
        // execute simulation only for selected units
        for (ControlUnit* cu : *subsection) {
            if (!cu->compute_next_value(ts, dayOfWeek_l, hourOfDay_l))
                return false;
        }
    }
    //
    // set new global radiation values for PV and wind speed
    float total_load = 0.0;
    float total_residential_load   = 0.0;
    float total_residential_demand = 0.0;
    global::unit_open_space_pv->compute_next_value(ts);
    global::unit_open_space_wind->compute_next_value(ts);
    total_load -= global::unit_open_space_pv->get_current_feedin_kW();
    total_load -= global::unit_open_space_wind->get_current_feedin_kW();
    total_load += global::residual_gridload_kW[ts-1];
    //
    // loop over all substations: compute new load values
    // and calculate total grid load
    Substation*const* subList = Substation::GetArrayOfInstances();
    const size_t nSubst = Substation::GetNumberOfInstances();

    if (output::substation_output != NULL && output::substation_output_details != NULL) {
        //
        // compute SOC over all batteries
        // TODO: inlcude batteries that are directly attached to an substation
        double totalBatteryCharge_kWh = ControlUnit::GetAllSimCompBatteriesCharge_kWh();
        double totalBatterySOC = 0.0;
        if (totalBatteryCapacity_kWh > 0)
            totalBatterySOC = totalBatteryCharge_kWh / totalBatteryCapacity_kWh;
        //
        // generate output
        *(output::substation_output) << ts << ","; // add timestep to output
        *(output::substation_output_details) << ts << ",";
        for (size_t i = 0; i < nSubst; i++) {
            subList[i]->calc_load(); // TODO: This could be moved outside the test, if substation output exist -> Propably we need information on grid load later
            float current_station_load = subList[i]->get_station_load();
            float current_station_resident_load   = subList[i]->get_residential_load();
            float current_station_resident_demand = subList[i]->get_residential_demand();
            total_load += current_station_load;
            total_residential_load   += current_station_resident_load;
            total_residential_demand += current_station_resident_demand;
            // stuff for output
            *(output::substation_output) << current_station_load << ",";
            *(output::substation_output_details) << current_station_resident_load   << ",";
            *(output::substation_output_details) << current_station_resident_demand << ",";
        }
        *(output::substation_output) << global::unit_open_space_pv->get_current_feedin_kW()   << ",";
        *(output::substation_output) << global::unit_open_space_wind->get_current_feedin_kW() << ",";
        *(output::substation_output) << totalBatterySOC << ",";
        *(output::substation_output) << total_load << "\n"; // add total load to output
        *(output::substation_output_details) << total_residential_load << ",";
        *(output::substation_output_details) << total_residential_demand << "\n";
    }

    // Output current time step
    output_counter++;
    if (output_counter >= OUTPUT_STATUS_OUTPUT_FREQ) {
        output_counter = 0;
        if (output_prefix[0] != '\0')
            std::cout << "\r" << output_prefix;
        std::cout << "    Step " << ts << " of " << Global::get_n_timesteps() << " has been computed.";
        if (output_prefix[0] != '\0')
            std::cout << std::flush;
        else
            std::cout << std::endl;
    } /*else {
        std::cout << ".";
    }*/

    return true;

}

bool simulation::runSimulationForAllVariations(unsigned long scenario_id) {
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
        unsigned long n_param_vari_combis = global::parameter_var_list->size();
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

                if        (var_name_and_val.first.compare("expansion PV kWp static") == 0) {
                    for (size_t i = 0; i < ControlUnit::GetNumberOfInstances(); i++)
                        cuList[i]->set_exp_pv_params_A(var_name_and_val.second);
                    cParamVals.exp_pv_kWp_static     = var_name_and_val.second;
                    cParamVals.exp_pv_kWp_static_set = true;

                } else if (var_name_and_val.first.compare("expansion PV kWp min kWp for section usage") == 0) {
                    cParamVals.exp_pv_min_kWp_roof_sec     = var_name_and_val.second;
                    cParamVals.exp_pv_min_kWp_roof_sec_set = true;

                } else if (var_name_and_val.first.compare("expansion PV kWp max inst kWp per section") == 0) {
                    cParamVals.exp_pv_max_kWp_roof_sec     = var_name_and_val.second;
                    cParamVals.exp_pv_max_kWp_roof_sec_set = true;

                } else if (var_name_and_val.first.compare("expansion PV kWp per roof area in m2") == 0) {
                    cParamVals.exp_pv_kWp_per_m2     = var_name_and_val.second;
                    cParamVals.exp_pv_kWp_per_m2_set = true;

                } else if (var_name_and_val.first.compare("expansion BS P in kW")  == 0) {
                    for (size_t i = 0; i < ControlUnit::GetNumberOfInstances(); i++)
                        cuList[i]->set_exp_bs_maxP_kW(var_name_and_val.second);
                    cParamVals.exp_bs_maxP_kW = var_name_and_val.second;
                    cParamVals.exp_bs_maxP_kW_set = true;

                } else if (var_name_and_val.first.compare("expansion BS E in kWh") == 0) {
                    for (size_t i = 0; i < ControlUnit::GetNumberOfInstances(); i++)
                        cuList[i]->set_exp_bs_maxE_kWh(var_name_and_val.second);
                    cParamVals.exp_bs_maxE_kWh = var_name_and_val.second;
                    cParamVals.exp_bs_maxE_kWh_set = true;

                } else if (var_name_and_val.first.compare("expansion BS E:P ratio") == 0) {
                    for (size_t i = 0; i < ControlUnit::GetNumberOfInstances(); i++)
                        cuList[i]->set_exp_bs_E_P_ratio(var_name_and_val.second);
                    cParamVals.exp_bs_EP_ratio = var_name_and_val.second;
                    cParamVals.exp_bs_EP_ratio_set = true;

                } else if (var_name_and_val.first.compare("expansion BS initial SOC") == 0) {
                    std::cerr << "This is not implemented!" << std::endl;

                } else {
                    std::cerr << "Unknown parameter variable to vary: " << var_name_and_val.first << std::endl;
                    return false;
                }
            }
            // 1b. if expansion PV kWp min/max or kWp per m2 is set, apply it here
            //     it cannot be set directly at the place above (like the other parameters)
            //     as these values can only be set together
            if (cParamVals.exp_pv_min_kWp_roof_sec_set || cParamVals.exp_pv_max_kWp_roof_sec_set || cParamVals.exp_pv_kWp_per_m2_set) {
                float kWp_per_m2 = Global::get_exp_pv_kWp_per_m2();
                float min_kWp    = Global::get_exp_pv_min_kWp_roof_sec();
                float max_kWp    = Global::get_exp_pv_max_kWp_roof_sec();
                if (cParamVals.exp_pv_kWp_per_m2_set)       kWp_per_m2 = cParamVals.exp_pv_kWp_per_m2;
                if (cParamVals.exp_pv_min_kWp_roof_sec_set) min_kWp    = cParamVals.exp_pv_min_kWp_roof_sec;
                if (cParamVals.exp_pv_max_kWp_roof_sec_set) max_kWp    = cParamVals.exp_pv_max_kWp_roof_sec;
                for (size_t i = 0; i < ControlUnit::GetNumberOfInstances(); i++)
                    cuList[i]->set_exp_pv_params_B(kWp_per_m2, min_kWp, max_kWp);
            }

            //
            // 2. open output files
            output::initializeDirectoriesPerPVar();
            output::initializeSubstationOutput(scenario_id);
            output::initializeCUOutput(scenario_id);
            // 2.b output the current parameter variation combination
            output::outputCurrentParamVariCombi(cParamVals);
            output::outputCurrentCUSettings();
            //
            // 3. run the simulation
            cout << "Simulation run for parameter variation " << std::to_string(param_vari_combi_ind+1) << " of " << std::to_string(n_param_vari_combis) << "\n";
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
        output::initializeDirectoriesPerPVar();
        output::initializeSubstationOutput(scenario_id);
        output::initializeCUOutput(scenario_id);
        // 1.b output the current parameter variation combination
        output::CurrentParamValues cParamVals;
        output::outputCurrentParamVariCombi(cParamVals);
        output::outputCurrentCUSettings();
        //
        // 2.0 Reset internal states
        ControlUnit::ResetAllInternalStates();
        // 2.1 Run the simulation
        cout << "Main simulation run:\n";
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

bool simulation::runSimulationFAVsAndSAC(float expansion_matrix_rel_freq[16][16], unsigned long expansion_matrix_abs_freq[16][16], unsigned long scenario_id) {
    //
    // This function executes the expansion / adds sim. added components to selected CUs
    // and then calls runSimulationForAllVariations(int)
    //
    // 1) add expansion to units in the other cases
    //
    expansion::add_expansion_to_units(expansion_matrix_rel_freq, expansion_matrix_abs_freq);
    //
    // 2) run the simulation (for all parameter variations or a single run)
    return runSimulationForAllVariations(scenario_id);
}

bool simulation::runCompleteSimulation(float expansion_matrix_rel_freq[16][16], unsigned long expansion_matrix_abs_freq[16][16], unsigned long scenario_id) {
    if (Global::get_repetitions_selected()) {
        // if repetition is selected, this function will handle it
        for (unsigned int cRep = 1; cRep <= Global::get_n_repetitions(); cRep++) {
            // set current repetion number to a global variable, so that output paths and so on can be adjusted
            global::current_repetition_counter = cRep;
            // Initialize global output directories
            output::initializeDirectoriesBase(scenario_id);
            // run the main part
            if (!simulation::runSimulationFAVsAndSAC(expansion_matrix_rel_freq, expansion_matrix_abs_freq, scenario_id)) {
                return false;
            }
            // Remove all added components and reset all internal states
            ControlUnit::RemoveAllSimAddedComponents();
            ControlUnit::ResetAllInternalStates();
            // increment seed, if a seed is set, otherwise we would get the same result in the next repetition and thus a repetition would be useless
            if (Global::is_seed_set()) {
                Global::increment_seed();
            }
        }
        return true;
    } else {
        // Initialize global output directories
        output::initializeDirectoriesBase(scenario_id);
        // run the main part
        return simulation::runSimulationFAVsAndSAC(expansion_matrix_rel_freq, expansion_matrix_abs_freq, scenario_id);
    }
}
