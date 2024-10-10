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

#include <string>
#include <vector>

#include "units.h"

/**
 * This namespace contains all functions required for loading the expansion matrix,
 * storing information about which expansion possibility is possible and 
 * checking the expansion matrix itself.
 **/
namespace expansion {

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

    inline double final_cumsum_of_added_pv_kWp; ///< The cummulative sum of added roof-top PV power in kWp after the SAC planning
    inline double final_cumsum_of_added_bs_kWh; ///< The cummulative sum of added battery capacity in kWh after the SAC planning
    inline ulong final_count_of_added_hps;      ///< Total number of added HPs
    inline ulong final_count_of_added_evchsts;  ///< Total number of added EV charging stations
    inline ulong final_count_of_added_evs;      ///< Total number of added EVs

    int expCombiMatrixOrderToBitRepr(int indexMatO);
    int expCombiBitReprToMatrixOrder(int bitRepr);
    int genExpCombiAsBitRepr(bool has_pv, bool has_bs, bool has_hp, bool bas_wb);
    bool isExpCombiPossible(int current_scenario_number, int b);
    std::string expCombiMatrixOrderToString(int indexMatO); ///< Converts a expansion combination in matrix order to a human-readable string

    bool load_expansion_matrix(float expansion_matrix[16][16]);
    bool verify_expansion_matrix(float expansion_matrix[16][16]);
    void add_expansion_to_units(float expansion_matrix_rel_freq[16][16],
                                unsigned long expansion_matrix_abs_freq[16][16]);
}




#endif
