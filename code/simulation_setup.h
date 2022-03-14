/*
 * 
 * simulation_setup.h
 *
 * Contains all code required for initializing the
 * simulation and loading central data
 *
 * */

#ifndef SIMULATION_SETUP_H
#define SIMULATION_SETUP_H

#include "global.h"



namespace configld {
    /*
     * This namespace contains all functions required for loading central
     * infromation and the scenario file
     */

    bool parse_scenario_file(int scenario_id);
    bool load_data_from_central_database(const char* filepath);

}

namespace expansion {
    /*
     * This namespace contains all functions required for loading the expansion matrix,
     * storing information about which expansion possibility is possible and 
     * checking the expansion matrix itself.
     */

    //
    // General information about the expansion matrix is defined here
    //
    // There are two different ways in which a expansion is defined
    //  1. as the index it appears in the 
    const int MaskNothing = 0b0000;
    const int MaskPV      = 0b0001;
    const int MaskBS      = 0b0010;
    const int MaskHP      = 0b0100;
    const int MaskWB      = 0b1000;

    int expCombiMatrixOrderToBitRepr(int indexMatO);
    int expCombiBitReprToMatrixOrder(int bitRepr);
    int genExpCombiAsBitRepr(bool has_pv, bool has_bs, bool has_hp, bool bas_wb);
    bool isExpCombiPossible(int current_scenario_number, int b);

    bool load_expansion_matrix(float expansion_matrix[16][16]);
    bool verify_expansion_matrix(float expansion_matrix[16][16]);
    void add_expansion_to_units(float expansion_matrix_rel_freq[16][16], int expansion_matrix_abs_freq[16][16]);
}






#endif
