#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <limits>
//#include <format> // not available for GCC !!!
#include <cstdio>

#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/program_options.hpp>

#include "global.h"

#include "components.h"
#include "units.h"
#include "output.h"
#include "simulation_setup.h"
#include "simulation_logic.h"


using namespace std;
namespace bpopts = boost::program_options;



/*
Return values of the program:
0	Normal execution, no errors occured
1	Wrong parameters
2	A required file was not found or error during database connections
3	Errors duing the simulation
4	Erroneous input files
5   No simulation executed / e.g. help displayed
*/

int main(int argc, char* argv[]) {

	//
	// parsing command line arguments
	//
	int scenario_id;
    string config_filepath;
    //
    bpopts::options_description opts_desc("Options");
    opts_desc.add_options()
        ("help,h",                            "Show help")
        ("config",   bpopts::value<string>(), "Path to the json configuration file")
        ("pvar",     bpopts::value<int>(),    "ID of parameter variation, that should be applied")
        ("scenario", bpopts::value<int>(),    "ID of the scenario that should be used, regardless of parameter variation is selected or not")
        ("cu-output,c", bpopts::value<string>(), "Modify output behavior for individual control units: 'no' switches off output completly, 'single' creates a single output instead of one per unit");
    bpopts::positional_options_description opts_desc_pos;
    opts_desc_pos.add("scenario", -1);
    bpopts::variables_map opts_vals;
    try {
        bpopts::command_line_parser parser{argc, argv};
        parser.options(opts_desc).positional(opts_desc_pos);
        bpopts::parsed_options parsed_options = parser.run();
        bpopts::store( parsed_options, opts_vals );
        bpopts::notify(opts_vals);
    } catch (const bpopts::error &err) {
        cerr << "Error when parsing command line arguments:" << "\n";
        cerr << err.what() << endl;
        return 1;
    }
    // now, command line arguments are parsed
    // we now set the internal variables accordingly
    if (opts_vals.count("help") > 0) {
        cerr << opts_desc << endl;
        cerr << "Usage: simulation [-h] [--config PATH] [--pvar IDvar] [[--scenario] IDscenario]" << endl;
        return 5;
    }
    if (opts_vals.count("config") > 0) {
        config_filepath = opts_vals["config"].as<string>();
    } else {
        config_filepath = "../config/simulation_config.json";
    }
    if (opts_vals.count("pvar") > 0) {
        Global::set_pvar_vals(true, opts_vals["pvar"].as<int>() );
    } else {
        Global::set_pvar_vals(false, 0);
    }
    if (opts_vals.count("scenario") > 0) {
		scenario_id = opts_vals["scenario"].as<int>();
    } else {
		scenario_id = 1;
    }
    if (opts_vals.count("cu-output") > 0) {
        string cu_output = opts_vals["cu-output"].as<string>();
        if (cu_output == "no") {
            Global::set_output_mode_per_cu(global::OutputModePerCU::NoOutput);
        } else if (cu_output == "single") {
            Global::set_output_mode_per_cu(global::OutputModePerCU::SingleFile);
        } else {
            // invalid argument
            cerr << "Error when parsing command line arguments: invalid option for --cu-output given!" << endl;
            return 1;
        }
    } else {
        Global::set_output_mode_per_cu(global::OutputModePerCU::IndividualFile);
    }

	cout << "Initializing the simulation for scenario ID " << scenario_id << endl;

	//
	// creating variables and initializing the
	// static attributes in the classes
	//
	float expansion_matrix_rel_freq[16][16] = {0};
	int   expansion_matrix_abs_freq[16][16] = {0};
	Global::InitializeStaticVariables();

	//
	// open and parse global settings file
	//
	if (!configld::load_config_file(scenario_id, config_filepath)) {
		return 2;
	}

	//
	// Load the expansion matrix
	//
	if (!expansion::load_expansion_matrix(expansion_matrix_rel_freq)) {
		return 2;
	}

	//
	// Check the expansion matrix
	//
	if (!expansion::verify_expansion_matrix(expansion_matrix_rel_freq)) {
		return 4;
	}

	//
	// Create the control units and load the data for the connected measurement units
	//
	if (!configld::load_data_from_central_database((Global::get_input_path() + "Merged_Information.db").c_str())) {
		return 2;
	}

	// this is not allowed anymore, as the object are deleted twice!!!
	//MeasurementUnit mu(1, NULL, 0/*meloID,string melo, locID*/);
	// do this:
	//MeasurementUnit* mu = new MeasurementUnit(1, NULL, 0/*meloID,string melo, locID*/);
	//mu->load_data("../data/input/SeparatedSmartMeterData/2.csv");
	// TODO: Make all consturctors for the units private, so this will not happen !!!

	//
	// open output files
	//
    output::initializeDirectory(scenario_id);
	output::initializeSubstationOutput(scenario_id);
	output::initializeCUOutput(scenario_id);

	//
	// bevore starting the simulation:
	// check if all global variables are set -> if not, error!
	//
	if (!global::all_variables_initialized() || !Global::AllVariablesInitialized()) {
		cout << "Some global variables are not initialized!" << endl;
		return 3;
	}

	//
	// Add expansion[s] to the control units
	//
	expansion::add_expansion_to_units(expansion_matrix_rel_freq, expansion_matrix_abs_freq, scenario_id);

	//
	// Run the simulation
	//
	if (!simulation::runSimulation()) {
		cerr << "Error during simulation run!" << endl;
		return 3;
	}
	
	//
	// clean up
	//
	MeasurementUnit::VacuumInstancesAndStaticVariables();
	ControlUnit::VacuumInstancesAndStaticVariables();
	Substation::VacuumInstancesAndStaticVariables();
	output::closeOutputs();
	Global::DeleteStaticVariables();
	global::vacuum();

	return 0;
}

