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
    const unsigned char EXPMAT_POS_   = 0;
    const unsigned char EXPMAT_POS_PV = 1;
    const unsigned char EXPMAT_POS_BS = 2;
    const unsigned char EXPMAT_POS_HP = 3;
    const unsigned char EXPMAT_POS_WB = 4;
    const unsigned char EXPMAT_POS_PV_BS = 5;
    const unsigned char EXPMAT_POS_PV_HP = 6;
    const unsigned char EXPMAT_POS_PV_WB = 7;
    const unsigned char EXPMAT_POS_BS_HP = 8;
    const unsigned char EXPMAT_POS_BS_WB = 9;
    const unsigned char EXPMAT_POS_HP_WB = 10;
    const unsigned char EXPMAT_POS_PV_BS_HP = 11;
    const unsigned char EXPMAT_POS_PV_BS_WB = 12;
    const unsigned char EXPMAT_POS_PV_HP_WB = 13;
    const unsigned char EXPMAT_POS_BS_HP_WB = 14;
    const unsigned char EXPMAT_POS_PV_BS_HP_WB = 15;
    
    bool is_expansion_combination_possible(int current_scenario_number, int b);

    bool load_expansion_matrix(float expansion_matrix[16][16]);
    bool verify_expansion_matrix(float expansion_matrix[16][16]);
}






#endif
