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

    bool runSimulationForOneParamSetting(std::vector<ControlUnit*>* subsection = NULL, const char* output_prefix = "");

    /**
     * Runs one step of the simulation.
     * Should only be called by simulation::runSimulationForOneParamSetting().
     * 
     * @param ts current time step (starting from 1)
     * @param dayOfWeek_l The day of week of the current time step (left aligned)
     * @param hourOfWeek_l The hour of the day of the current time step (left aligned)
     * @param totalBatteryCapacity_kWh Total capacity over all batteries in kWh; required for computation of overall BS SOC
     * @param output_prefix (optional, default "") A prefix that should be added to the output file names
     * @param subsection (optional, default NULL) Execut the simulation for the current step only for a given sub-set of all known control units
     */
    bool oneStep(unsigned long ts, unsigned int dayOfWeek_l, unsigned int hourOfDay_l, double totalBatteryCapacity_kWh, const char* output_prefix = "", std::vector<ControlUnit*>* subsection = NULL);

    bool runSimulationForAllVariations(int scenario_id);
    bool runSimulationFAVsAndSAC(float expansion_matrix_rel_freq[16][16], unsigned long expansion_matrix_abs_freq[16][16], int scenario_id); ///< Runs the complete simulation for all parameter variations and add sim. added components to the CUs

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
    bool runCompleteSimulation(float expansion_matrix_rel_freq[16][16], unsigned long expansion_matrix_abs_freq[16][16], int scenario_id);

}
