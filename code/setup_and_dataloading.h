/*
 * 
 * setup_and_dataloading.h
 *
 * Contains all code required for initializing the
 * simulation and loading central data
 *
 * */

#ifndef SETUP_AND_DATALOADING_H
#define SETUP_AND_DATALOADING_H

#include <string>



namespace configld {
    /*
     * This namespace contains all functions required for loading central
     * infromation and the scenario file
     */

    bool load_config_file(int scenario_id, std::string& filepath);
    //bool parse_scenario_file(int scenario_id);
    bool load_data_from_central_database(const char* filepath);

}




#endif
