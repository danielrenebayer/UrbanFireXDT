
#include "simulation_logic.h"
using namespace simulation;

#include <algorithm>
#include <ctime>
#include <iostream>
#include <sstream>

#include "helper.h"
#include "global.h"
#include "output.h"
#include "sac_planning.h"
#include "units.h"


bool simulation::runSimulationForOneParamSetting() {
    //
    // This function loops over all time steps as they are defined in the data.
    // If multiple settings (e.g. a parameter variation) should be applied,
    // one has to call this function more often.
    //
    std::cout << "Run simulation for a complete episode ..." << std::endl;

    unsigned long n_tsteps = Global::get_n_timesteps();
    struct tm* tm_start = Global::get_ts_start_tm();
    struct tm* tm_end   = Global::get_ts_end_tm();
    //time_t t_start = mktime(tm_start);
    //time_t t_end   = mktime(tm_end);
    struct tm* current_tm;
    bool sim_started = false; // gets true, if simulation range (as given by tm_start) has been reached
    // main loop
    for (unsigned long ts = 1; ts <= n_tsteps; ts++) {
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

bool simulation::oneStep(unsigned long ts) {
    //
    // Run one time step of the simulation
    // Return false if an error occurs during execution.
    // Parameter ts is the current time step
    //

    //int tsID = ts - 1;

    // loop over all control units:
    //    set new values and execute next actions
    ControlUnit*const* cuList = ControlUnit::GetArrayOfInstances();
    const size_t nCUs = ControlUnit::GetNumberOfInstances();
    for (size_t i = 0; i < nCUs; i++) {
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
    const size_t nSubst = Substation::GetNumberOfInstances();

    if (output::substation_output != NULL) {
        *(output::substation_output) << ts << ","; // add timestep to output
        for (size_t i = 0; i < nSubst; i++) {
            float current_station_load = subList[i]->calc_load();
            total_load += current_station_load;
            // stuff for output
            *(output::substation_output) << current_station_load << ",";
        }
        *(output::substation_output) << global::unit_open_space_pv->get_current_feedin_kW()   << ",";
        *(output::substation_output) << global::unit_open_space_wind->get_current_feedin_kW() << ",";
        *(output::substation_output) << total_load << "\n"; // add total load to output
    }

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

bool simulation::runSimulationFAVsAndSAC(float expansion_matrix_rel_freq[16][16], long expansion_matrix_abs_freq[16][16], int scenario_id) {
    //
    // This function executes the expansion / adds sim. added components to selected CUs
    // and then calls runSimulationForAllVariations(int)
    //
    if (Global::get_cu_selection_mode_fca() == global::CUSModeFCA::BestSSR) {
        //
        // 1) If selection mode is taking those with best SSR
        //
        // 1.0) initial checks
        /*if (!Global::get_comp_eval_metrics()) {
            cerr << "Error: Option 'metrics / m' not set even though CU selection for sim. add. components should be done according to best SSR.\nThis is impossible. Set 'metrics' as parameter or change CU selection mode." << endl;
            return false;
        }*/
        // 1.1) plan expansion as they would happen, but with random shuffling whatever is selected
        expansion::add_expansion_to_units(expansion_matrix_rel_freq, expansion_matrix_abs_freq, true);
        // 1.2a) add PV installations to all CUs, if they do not already have one
        ControlUnit*const* cuList = ControlUnit::GetArrayOfInstances();
        const size_t nCUs = ControlUnit::GetNumberOfInstances();
        for (size_t i = 0; i < nCUs; i++) {
            if (!cuList[i]->has_pv())
                cuList[i]->add_exp_pv();
        }
        // 1.2b) add Battery to all CUs (that have a sim. added PV-installation), if they do not already have a battery
        for (size_t i = 0; i < nCUs; i++) {
            if ( !cuList[i]->has_bs() && ( cuList[i]->get_exp_combi_bit_repr_sim_added() & expansion::MaskBS ) )
                cuList[i]->add_exp_bs();
        }
        // 1.3) execute the simulation once for the given scenario
        bool no_error = runSimulationForOneParamSetting();
        if (!no_error) { 
            cerr << "Error during selection of the CUs for adding simulated components." << endl;
            return false;
        }
        // 1.4) sort control units according to the SSR which has been computed in 1.2)
        // 1.4.a) Collect the SSR values
        vector<pair<double, ControlUnit*>> ssr_cu_pair_vector;
        ssr_cu_pair_vector.reserve( ControlUnit::GetNumberOfInstances() );
        for (size_t i = 0; i < nCUs; i++) {
            ssr_cu_pair_vector.emplace_back(cuList[i]->get_SSR(), cuList[i]); // SSR, pointer to the object
        }
        // 1.4.b) Sort the SSR / CU-Pointer vactor
        sort(ssr_cu_pair_vector.begin(),
             ssr_cu_pair_vector.end(),
             [](pair<double, ControlUnit*> a, pair<double, ControlUnit*> b) { return a.first > b.first; });
        //
        // 1.4.c) output metrics
        output::outputMetrics(true);
        //
        // 1.5) Reset internal variables
        ControlUnit::ResetAllInternalStates();
        //
        // 1.6) Remove sim. added components from all CUs
        ControlUnit::RemoveAllSimAddedComponents();
        //
        // 1.7) Plan expansion with the new order
        vector<ControlUnit*> ordered_cu_list;
        ordered_cu_list.reserve(ssr_cu_pair_vector.size());
        transform(ssr_cu_pair_vector.begin(), ssr_cu_pair_vector.end(),
                  back_inserter(ordered_cu_list),
                  [](auto& pair){ return pair.second; });
        expansion::add_expansion_to_units(expansion_matrix_rel_freq, expansion_matrix_abs_freq, false, &ordered_cu_list);
    } else {
        //
        // 2) add expansion to units in the other cases
        //
        expansion::add_expansion_to_units(expansion_matrix_rel_freq, expansion_matrix_abs_freq);
    }
    //
    // 3) run the simulation (for all parameter variations or a single run)
    return runSimulationForAllVariations(scenario_id);
}
