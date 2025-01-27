
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
#include "worker_threads.hpp"

//
// Internal variables required for output control
unsigned int output_counter = 0;
// #define OUTPUT_STATUS_OUTPUT_FREQ 24*3
#define OUTPUT_STATUS_OUTPUT_FREQ 24*24


bool simulation::runSimulationForOneParamSetting(CUControllerThreadGroupManager* thread_manager /* = NULL*/, std::vector<ControlUnit*>* subsection /* = NULL */, const char* output_prefix /* = "" */) {
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
    unsigned long ts_start = Global::get_first_timestep();
    unsigned long ts_end   = Global::get_last_timestep();
    struct tm* current_tm_l; // left-aligned time stamp as struct tm
    unsigned int last_step_weekday = 0; // required to identify new weeks
    unsigned long week_number = 1;
    //
    // main loop
    for (unsigned long ts = ts_start; ts <= ts_end; ts++) {
        // get current time as struct tm
        current_tm_l = global::time_localtime_l->at(ts - 1);
        //
        // day and time handling
        unsigned int dayOfWeek_l = (unsigned int) (current_tm_l->tm_wday); // get day of week in the format 0->Monday, 6->Sunday
        // has a new week started?
        if (last_step_weekday > dayOfWeek_l) {
            // only execute this during the main run (checked by subsection == NULL)
            if (Global::get_compute_weekly_metrics() && subsection == NULL) {
                std::list<std::string*> *output_list = new std::list<std::string*>();
                // iterate over all control units, get the weekly metrics string (and reset weekly summation variables automatically)
                const std::vector<ControlUnit*>& units_list = ControlUnit::GetArrayOfInstances();
                for (ControlUnit* current_unit : units_list) {
                    output_list->push_back(
                        current_unit->get_metrics_string_weekly_wr(week_number)
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
        if (!oneStep(ts, totalBatteryCapacity_kWh, thread_manager, output_prefix, subsection)) return false;
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

bool simulation::oneStep(unsigned long ts,
                        double totalBatteryCapacity_kWh,
                        CUControllerThreadGroupManager* thread_manager /* = NULL */,
                        const char* output_prefix /* = "" */,
                        std::vector<ControlUnit*>* subsection /* = NULL */) {
    //
    // Run one time step of the simulation
    // Return false if an error occurs during execution.
    // Parameter ts is the current time step
    //

    //int tsID = ts - 1;

    // loop over all control units:
    //    set new values and execute next actions
    if (thread_manager == NULL)
    {
        // Case 1: No parallelization
        if (subsection == NULL) {
            // Case 1a: Execute simulation for all substations
            const std::vector<ControlUnit*>& units_list = ControlUnit::GetArrayOfInstances();
            for (ControlUnit* current_unit : units_list) {
                if (!current_unit->compute_next_value(ts))
                    return false;
            }
        } else {
            // Case 1b: Execute simulation only for selected units
            for (ControlUnit* cu : *subsection) {
                if (!cu->compute_next_value(ts))
                    return false;
            }
        }
    }
    else
    {
        // Case 2: Parallelization - subsection is ignored!
        // thread_manager MUST be initialized with subsection as argument!
        thread_manager->executeOneStep(ts);
        thread_manager->waitForWorkersToFinish();
        // TODO: If there occurs one false from cu->compute_next_values -> Break complete run![]
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
    const std::vector<Substation*>& substations_list = Substation::GetArrayOfInstances();

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
        for (Substation* s : substations_list) {
            s->calc_load(); // TODO: This could be moved outside the test, if substation output exist -> Propably we need information on grid load later
            float current_station_load = s->get_station_load();
            float current_station_resident_load   = s->get_residential_load();
            float current_station_resident_demand = s->get_residential_demand();
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

bool simulation::runSimulationForAllVariations(unsigned long scenario_id, CUControllerThreadGroupManager* thread_manager /*= NULL */) {
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
            const std::vector<ControlUnit*>& units_list = ControlUnit::GetArrayOfInstances();
            for (auto& var_name_and_val : list_of_variables) {

                if        (var_name_and_val.first.compare("expansion PV kWp static") == 0) {
                    for (ControlUnit* current_unit : units_list)
                        current_unit->set_exp_pv_params_A(var_name_and_val.second);
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
                    for (ControlUnit* current_unit : units_list)
                        current_unit->set_exp_bs_maxP_kW(var_name_and_val.second);
                    cParamVals.exp_bs_maxP_kW = var_name_and_val.second;
                    cParamVals.exp_bs_maxP_kW_set = true;

                } else if (var_name_and_val.first.compare("expansion BS E in kWh") == 0) {
                    for (ControlUnit* current_unit : units_list)
                        current_unit->set_exp_bs_maxE_kWh(var_name_and_val.second);
                    cParamVals.exp_bs_maxE_kWh = var_name_and_val.second;
                    cParamVals.exp_bs_maxE_kWh_set = true;

                } else if (var_name_and_val.first.compare("expansion BS E:P ratio") == 0) {
                    for (ControlUnit* current_unit : units_list)
                        current_unit->set_exp_bs_E_P_ratio(var_name_and_val.second);
                    cParamVals.exp_bs_EP_ratio = var_name_and_val.second;
                    cParamVals.exp_bs_EP_ratio_set = true;

                } else if (var_name_and_val.first.compare("expansion BS initial SOC") == 0) {
                    std::cerr << "This is not implemented!" << std::endl;

                } else if (var_name_and_val.first.compare("control horizont in ts") == 0) {
                    ControlUnit::ChangeControlHorizonInTS( (uint) (var_name_and_val.second) );
                    cParamVals.control_horizon_in_ts     = (uint) (var_name_and_val.second);
                    cParamVals.control_horizon_in_ts_set = true;

                } else if (var_name_and_val.first.compare("control update freq in ts") == 0) {
                    Global::UnlockAllVariables();
                    Global::set_control_update_freq_in_ts( (uint) (var_name_and_val.second) );
                    Global::LockAllVariables();
                    cParamVals.control_update_freq_in_ts     = (uint) (var_name_and_val.second);
                    cParamVals.control_update_freq_in_ts_set = true;

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
                for (ControlUnit* current_unit : units_list)
                    current_unit->set_exp_pv_params_B(kWp_per_m2, min_kWp, max_kWp);
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
            bool no_error = runSimulationForOneParamSetting(thread_manager);
            if (!no_error) return false;
            //
            // 4a. output metrics (if selected)
            output::outputMetrics();
            output::outputEVMetrics();
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
        bool no_error = runSimulationForOneParamSetting(thread_manager);
        if (!no_error) return false;
        //
        // 3a. output metrics (if selected)
        output::outputMetrics();
        output::outputEVMetrics();
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
    // 2a) Initialize and start the CUControllerThreadGroup if multi-threading is selected
    CUControllerThreadGroupManager* tgm = NULL;
    if (Global::get_n_threads() >= 1) {
        tgm = new CUControllerThreadGroupManager();
        tgm->startAllWorkerThreads();
    }
    // 2b) run the simulation (for all parameter variations or a single run)
    bool return_value = runSimulationForAllVariations(scenario_id, tgm);
    // 2c) delete the thread group manager
    if (tgm != NULL) {
        tgm->stopAllWorkerThreads();
        delete tgm;
    }
    //
    return return_value;
}

bool simulation::runCompleteSimulation(float expansion_matrix_rel_freq[16][16], unsigned long expansion_matrix_abs_freq[16][16], unsigned long scenario_id) {
    // Initialize and precompute EV states
    ControlUnit::PreprocessEVData();
    //
    bool return_value = true;
    if (Global::get_repetitions_selected()) {
        // if repetition is selected, this function will handle it
        for (unsigned int cRep = 1; cRep <= Global::get_n_repetitions(); cRep++) {
            // set current repetion number to a global variable, so that output paths and so on can be adjusted
            global::current_repetition_counter = cRep;
            // Initialize global output directories
            output::initializeDirectoriesBase(scenario_id);
            // run the main part
            if (!simulation::runSimulationFAVsAndSAC(expansion_matrix_rel_freq, expansion_matrix_abs_freq, scenario_id)) {
                return_value = false;
                break;
            }
            // Remove all added components and reset all internal states
            ControlUnit::RemoveAllSimAddedComponents();
            ControlUnit::ResetAllInternalStates();
            // increment seed, if a seed is set, otherwise we would get the same result in the next repetition and thus a repetition would be useless
            if (Global::is_seed_set()) {
                Global::increment_seed();
            }
        }
    } else {
        // Initialize global output directories
        output::initializeDirectoriesBase(scenario_id);
        // run the main part
        return_value = simulation::runSimulationFAVsAndSAC(expansion_matrix_rel_freq, expansion_matrix_abs_freq, scenario_id);
    }

    return return_value;
}
