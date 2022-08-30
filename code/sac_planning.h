/*
 * 
 * sac_planning.h
 *
 * Planning of the Simulatively Added Components (to the control units)
 *
 * This file contains all code required for the planning of which
 * simulatively added components are added to which control unit.
 *
 * */

#ifndef SAC_PLANNING_H
#define SAC_PLANNING_H

#include <vector>

#include "units.h"


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
    void add_expansion_to_units(float expansion_matrix_rel_freq[16][16],
                                int   expansion_matrix_abs_freq[16][16],
                                int   scenario_id,
                                bool  random_anyway_no_output = false,
                                std::vector<ControlUnit*>* ordered_list = NULL);
}




#endif
