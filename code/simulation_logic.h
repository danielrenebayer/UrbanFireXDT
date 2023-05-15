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

    bool runSimulationForOneParamSetting(std::vector<ControlUnit*>* subsection = NULL);
    bool oneStep(unsigned long ts, double totalBatteryCapacity_kWh, std::vector<ControlUnit*>* subsection = NULL);
    bool runSimulationForAllVariations(int scenario_id);
    bool runSimulationFAVsAndSAC(float expansion_matrix_rel_freq[16][16], long expansion_matrix_abs_freq[16][16], int scenario_id); ///< Runs the complete simulation for all parameter variations and add sim. added components to the CUs

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
    bool runCompleteSimulation(float expansion_matrix_rel_freq[16][16], long expansion_matrix_abs_freq[16][16], int scenario_id);

}
