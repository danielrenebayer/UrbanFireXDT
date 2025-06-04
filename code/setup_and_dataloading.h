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

    bool load_config_file(unsigned long scenario_id, std::string& filepath); ///< Load the config file, that is passed as command line argument

    bool load_data_from_central_database(const char* filepath); ///< Load the complete simulation structure from the central database, and also load the recorded time series for the measurement units

    /**
     * @brief Outputs the current simulation configuration and all parameter settings to the specified output stream.
     *
     * This method prints detailed information about the simulation setup, including compiler settings,
     * scenario configuration, data paths, control strategies, selection criteria, and output modes.
     * The output format is human-readable and organized by category to facilitate inspection and debugging.
     *
     * The output can be redirected to any valid std::ostream, such as:
     * - std::cout (to print to console),
     * - std::ofstream (to write to a file),
     * - std::ostringstream (to capture as a string).
     *
     * @param out Reference to an output stream where the configuration should be written.
     */
    void output_variable_values(std::ostream& current_outstream);

}




#endif
