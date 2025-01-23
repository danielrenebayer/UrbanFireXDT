/*
 * simulation_logic.h
 *
 * This contains all functions required for running the
 * simulation, once it is created and the desired expansion
 * scenario is introduced.
 *
 */

#pragma once

#include <iostream>
#include <vector>

#include "units.h"

namespace simulation {

    /**
     * Runs one step of the simulation.
     * Should only be called by simulation::runSimulationForOneParamSetting().
     * 
     * Attention: If a thread_manager is given, the subsection-argument is ignored!
     * Thus, the thread manager must be initialized with the subsection!
     * 
     * @param ts current time step (starting from 1)
     * @param totalBatteryCapacity_kWh Total capacity over all batteries in kWh; required for computation of overall BS SOC
     * @param thread_manager: If multi-threading is active, an existing thread manager can be passed.
     * @param output_prefix (optional, default "") A prefix that should be added to the output file names
     * @param subsection (optional, default NULL) Execut the simulation for the current step only for a given sub-set of all known control units
     */
    bool oneStep(unsigned long ts,
                 double totalBatteryCapacity_kWh,
                 CUControllerThreadGroupManager* thread_manager = NULL,
                 const char* output_prefix = "",
                 std::vector<ControlUnit*>* subsection = NULL);

    /**
     * Runs the simulation for the complete defined time span for one parameter
     * setting and one repetition.
     * It can be either called by simulation::runSimulationForAllVariations()
     * or during the initialization phase for the SAC planning.
     * This function creates a a thread manager object, if multi-threading is
     * selected and no existing thread manager instance is passed.
     * 
     * @param thread_manager: If multi-threading is active, an existing thread manager can be passed.
     * @param subsection: If not NULL, the simulation is only executed for the given control units
     * @param output_prefix: The prefix for stdout information. If not empty, if changes the writing to stdout. 
     */
    bool runSimulationForOneParamSetting(CUControllerThreadGroupManager* thread_manager = NULL, std::vector<ControlUnit*>* subsection = NULL, const char* output_prefix = "");

    /**
     * Runs the simulation for all parameter variations.
     * If no variations as selected, it will only call
     * simulation::runSimulationForOneParamSetting() once.
     * This function also handels the correct initialization
     * of the output folders.
     */
    bool runSimulationForAllVariations(unsigned long scenario_id, CUControllerThreadGroupManager* thread_manager = NULL);

    /**
     * This function is called by simulation::runSimulationFAVsAndSAC().
     * It starts the SAC (simulatively added components) planning and then 
     * run the simulation for all parameter variations (FAVs).
     * It also initializes the thread group manager for the main run of the
     * simulation (i.e., NOT for the SAC planning, this is done by 
     * the SAC planning function itselfe) if multi-processing is selected.
     */
    bool runSimulationFAVsAndSAC(float expansion_matrix_rel_freq[16][16], unsigned long expansion_matrix_abs_freq[16][16], unsigned long scenario_id); ///< Runs the complete simulation for all parameter variations and add sim. added components to the CUs

    /**
     * Run the complete simulation for a given scenario.
     * If cmd-line argument 'repetitions' is set to a value higher 1, this function will handel this point.
     * It will call runSimulationFAVsAndSAC(...).
     * It also initializes the global output directory.
     * 
     * @param expansion_matrix_rel_freq reference to the exp. matrix with relative values
     * @param expansion_matrix_abs_freq reference to the exp. matrix with absolute values
     * @param scenario_id scenario id to use
     * 
     * @return false, if an error occurs during simulation run or expansion planning, otherwise true
     */
    bool runCompleteSimulation(float expansion_matrix_rel_freq[16][16], unsigned long expansion_matrix_abs_freq[16][16], unsigned long scenario_id);

}
