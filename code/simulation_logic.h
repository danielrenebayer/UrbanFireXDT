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
    bool oneStep(int ts);
    bool runSimulationForAllVariations(int scenario_id);

}
