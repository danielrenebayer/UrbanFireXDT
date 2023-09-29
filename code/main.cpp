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
#include "setup_and_dataloading.h"
#include "simulation_logic.h"
#include "vehicles.h"


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
        ("repetitions,r", bpopts::value<uint>(),    "Number of times the simulation shoud be repeated - only usefull is random variables are used.")
      //("metrics,m",                         "SSC and SSR will be computed for every control unit. Therefore the complete time series has to be stored - this requires more RAM.")
        ("suof",     bpopts::value<unsigned long>()->default_value(1000), "Steps until output will be flushed, i.e. written to disk. Defaults to 1000.")
        ("cu-output,c", bpopts::value<string>(), "Modify output behavior for individual control units:  'off' or 'no' switches off output completly, 'single' creates a single output instead of one per unit, 'sl' on substation level (default)")
        ("st-output,t", bpopts::value<string>(), "Modify output behavior for substations: 'off' or 'no' switches off substation output completly, 'on' substation output (default)")
        ("seed,s",bpopts::value<unsigned int>(), "Sets the seed for the simulation run. By default, no seed is used.")
        ("weekly-metrics,w", bpopts::value<string>()->default_value("off"), "Controls the generation of the metrics computation on weekly level. Default is 'off'.");
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
        if (cu_output == "no" || cu_output == "off") {
            Global::set_output_mode_per_cu(global::OutputModePerCU::NoOutput);
        } else if (cu_output == "single") {
            Global::set_output_mode_per_cu(global::OutputModePerCU::SingleFile);
        } else if (cu_output == "sl") {
            Global::set_output_mode_per_cu(global::OutputModePerCU::IndividualFile);
        } else {
            // invalid argument
            cerr << "Error when parsing command line arguments: invalid option for --cu-output / -c given!" << endl;
            return 1;
        }
    } else {
        Global::set_output_mode_per_cu(global::OutputModePerCU::IndividualFile);
    }
    if (opts_vals.count("st-output") > 0) {
        string st_output = opts_vals["st-output"].as<string>();
        if (st_output == "no" || st_output == "off") {
            Global::set_create_substation_output(false);
        } else if (st_output == "on") {
            Global::set_create_substation_output(true);
        } else {
            // invalid argument
            cerr << "Error when parsing command line arguments: invalid option for --st-output / -t given!" << endl;
            return 1;
        }
    } else {
        Global::set_create_substation_output(true);
    }
    if (opts_vals.count("repetitions") > 0) {
        unsigned int n_repetitions = opts_vals["repetitions"].as<unsigned int>();
        if (n_repetitions > 1) {
            Global::set_repetitions_selected(true);
            Global::set_n_repetitions(n_repetitions);
        } else {
            Global::set_repetitions_selected(false);
        }
    } else {
        Global::set_repetitions_selected(false);
    }
    if (opts_vals.count("seed") > 0) {
        unsigned int seed = opts_vals["seed"].as<unsigned int>();
        Global::set_seed(seed);
        EVFSM::SetSeed(seed);
    }
    if (opts_vals.count("weekly-metrics") > 0) {
        string wm_mode = opts_vals["weekly-metrics"].as<string>();
        if (wm_mode == "off" || wm_mode == "no") {
            Global::set_compute_weekly_metrics(false);
        } else if (wm_mode == "on" || wm_mode == "yes") {
            Global::set_compute_weekly_metrics(true);
        } else {
            cerr << "Error when parsing command line arguments: invalid option for --weekly-metrics / -w given!" << endl;
            return 1;
        }
    }
    /*
    if (opts_vals.count("metrics")) {
        Global::set_comp_eval_metrics(true);
    } else {
        Global::set_comp_eval_metrics(false);
    }
    */
    // get time for time measurement
    auto t1_A = std::chrono::steady_clock::now();
    auto t1_B = std::chrono::system_clock::now();
    // set output flush interval
    global::n_ts_between_flushs = opts_vals["suof"].as<unsigned long>();

	cout << "Initializing the simulation for scenario ID " << scenario_id << endl;

	//
	// creating variables and initializing the
	// static attributes in the classes
	//
    float expansion_matrix_rel_freq[16][16] = {0};
    unsigned long expansion_matrix_abs_freq[16][16] = {0};
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
	if (!configld::load_data_from_central_database((Global::get_input_path() + Global::get_structure_database_name()).c_str())) {
		return 2;
	}

	//
	// bevore starting the simulation:
	// check if all global variables are set -> if not, error!
	//
	if (!global::all_variables_initialized() || !Global::AllVariablesInitialized()) {
		cout << "Some global variables are not initialized!" << endl;
        Global::PrintUninitializedVariables();
        global::print_uninitialized_variables();
		return 3;
	}

    //
    // Output all variable values
    //
    configld::output_variable_values();

    // get time for time measurement
    auto t2_A = std::chrono::steady_clock::now();
    auto t2_B = std::chrono::system_clock::now();

    //
    // Plan and add which sim. added components should be added to which control units
    // and Run the simulation
    // - once (if no parameter variation is selected) or
    // - multiple times, if param. vari. is selected
    //
    if (!simulation::runCompleteSimulation(expansion_matrix_rel_freq, expansion_matrix_abs_freq, scenario_id)) {
        cerr << "Error during simulation run!" << endl;
        return 3;
    }

    // get time for time measurement
    auto t3_A = std::chrono::steady_clock::now();
    auto t3_B = std::chrono::system_clock::now();
	
	//
	// clean up
	//
    EVFSM::VaccuumStaticVariables();
    ComponentHP::VacuumStaticVariables();
	MeasurementUnit::VacuumInstancesAndStaticVariables();
	ControlUnit::VacuumInstancesAndStaticVariables();
	Substation::VacuumInstancesAndStaticVariables();
	Global::DeleteStaticVariables();
	global::vacuum();

    // get time for time measurement
    auto t4_A = std::chrono::steady_clock::now();
    auto t4_B = std::chrono::system_clock::now();
    cout << "Run-time information:\n";
    cout << "  Setup and data loading: " << std::chrono::duration_cast<std::chrono::seconds>(t2_A-t1_A).count() << "s\n";
    cout << "  Main run:               " << std::chrono::duration_cast<std::chrono::seconds>(t3_A-t2_A).count() << "s\n";
    cout << "  Clean up:               " << std::chrono::duration_cast<std::chrono::seconds>(t4_A-t3_A).count() << "s\n";
    cout << "  Complete run time:      " << std::chrono::duration_cast<std::chrono::seconds>(t4_A-t1_A).count() << "s\n";
    cout << "Run-time information:\n";
    cout << "  Setup and data loading: " << std::chrono::duration_cast<std::chrono::seconds>(t2_B-t1_B).count() << "s\n";
    cout << "  Main run:               " << std::chrono::duration_cast<std::chrono::seconds>(t3_B-t2_B).count() << "s\n";
    cout << "  Clean up:               " << std::chrono::duration_cast<std::chrono::seconds>(t4_B-t3_B).count() << "s\n";
    cout << "  Complete run time:      " << std::chrono::duration_cast<std::chrono::seconds>(t4_B-t1_B).count() << "s" << std::endl;

	return 0;
}

