/*
 * simulation_logic.h
 *
 * This contains all functions required for running the
 * simulation, once it is created and the desired expansion
 * scenario is introduced.
 *
 */

#pragma once

namespace simulation {

    bool runSimulationForOneParamSetting();
    bool oneStep(unsigned long ts);
    bool runSimulationForAllVariations(int scenario_id);
    bool runSimulationFAVsAndSAC(float expansion_matrix_rel_freq[16][16], long expansion_matrix_abs_freq[16][16], int scenario_id); ///< Runs the complete simulation for all parameter variations and add sim. added components to the CUs

}
