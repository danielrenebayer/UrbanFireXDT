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


/**
 * This namespace contains all functions required for loading central
 * infromation and the scenario file
 **/
namespace configld {

    bool load_config_file(int scenario_id, std::string& filepath); ///< Load the config file, that is passed as command line argument

    bool load_data_from_central_database(const char* filepath); ///< Load the complete simulation structure from the central database, and also load the recorded time series for the measurement units

    void output_variable_values(); ///< Outputs the configured variables to stdout

}




#endif
