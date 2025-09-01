#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <limits>
//#include <format> // not available for GCC !!!
#include <cstdio>

#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#ifndef PYTHON_MODULE
#include <boost/program_options.hpp>
#else
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // for converting C++ containers like std::map
#endif

#include "global.h"

#include "cache_helper.hpp"
#include "components.h"
#include "units.h"
#include "output.h"
#include "setup_and_dataloading.h"
#include "simulation_logic.h"
#include "status_output.hpp"
#include "vehicles.h"
#include "worker_threads.hpp"

#ifdef PYTHON_MODULE
#include "python_module.hpp"
#endif

using namespace std;
#ifndef PYTHON_MODULE
namespace bpopts = boost::program_options;
#endif



/**
 * @brief Entry point of the simulation.
 *
 * This function initializes and runs the simulation based on the provided command line parameters.
 * Finally, the storage is cleaned up again.
 *
 * @param argc Number of command line arguments.
 * @param argv Array of command line argument strings.
 *
 * @return Return code indicating the execution result of the program:
 *   - **0**  Normal execution, no errors occurred
 *   - **1**  Wrong parameters
 *   - **2**  Required file not found / error during database connections
 *   - **3**  Errors during the simulation
 *   - **4**  Erroneous input files
 *   - **5**  No simulation executed (e.g., help displayed)
 */
int main(int argc, char* argv[]) {

#ifndef PYTHON_MODULE
	//
	// parsing command line arguments
	//
	unsigned long scenario_id;
    string config_filepath;
    //
    bpopts::options_description opts_desc("Options");
    opts_desc.add_options()
        ("help,h",                            "Show help")
        ("config",   bpopts::value<string>(), "Path to the JSON configuration file")
        ("pvar",     bpopts::value<int>(),    "ID of parameter variation, that should be applied")
        ("scenario", bpopts::value<unsigned long>(),"ID of the scenario that should be used, regardless of parameter variation is selected or not")
        ("repetitions,r", bpopts::value<uint>(),    "Number of times the simulation should be repeated - only useful if random variables are used. Modifies the behavior of the 'seed' option.")
        ("stop-on-cc-err,e",                        "Stop simulation execution if an computation error occurs inside an optimization-based controller in any control unit.")
        ("n_threads,n",   bpopts::value<uint>()->default_value(3),    "Number of working threads. Defaults to 3. If all tasks should be done by the main thread, set this value to 0. If there is only one working thread, the concept is useless, as the main thread will wait for the workers to finish and does not do anything by itself.")
        ("work-stealing",                           "Enables work stealing between the worker threads if multiple threads are used.")
        ("max-parallel-opti-vars", bpopts::value<unsigned long>(), "The maximum number of optimization variables that will be called in parallel. Only effective if multiple threads are used. Due to algorithmic constraints, the maximum can be exceeded if there is a building with multiple EVs. See implementation for details.")
        ("suof",     bpopts::value<unsigned long>()->default_value(1000), "Steps until output will be flushed, i.e. written to disk. Defaults to 1000.")
        ("cu-output,c", bpopts::value<string>(), "Modify output behavior for individual control units:  'off' or 'no' switches off output completely, 'single' creates a single output instead of one per unit, 'sl' on substation level (default)")
        ("st-output,t", bpopts::value<string>(), "Modify output behavior for substations: 'off' or 'no' switches off substation output completely, 'on' substation output (default)")
        ("ccmd-output", bpopts::value<string>(), "Modify output behavior for details of the Control Commands (inside each control unit): 'off' disables all output (default), 'on' or 'all' outputs detailed information about the control commands for all control units for every time step (that have an optimized controller).")
        ("ev-output",   bpopts::value<string>(), "Modify output behavior for EVs: 'off' disables all output (default), 'all' outputs information for all EVs for every time step.")
        ("seed,s",bpopts::value<unsigned int>(), "Sets the seed for the simulation run. By default, no seed is used. If the simulation is repeated more than once, the seed gets incremented by one per repetition.")
        ("weekly-metrics,w", bpopts::value<string>()->default_value("off"), "Controls the generation of the metrics computation on weekly level. Default is 'off'.")
        ("cu-output-selection,o", bpopts::value<string>(), "If not given, all control units will generate output for the options `--cu-output`, `--ccmd-output` and `--ev-output`. If a comma separated list of control unit IDs is given, the output will only be generated for these units.");
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
		scenario_id = opts_vals["scenario"].as<unsigned long>();
    } else {
		scenario_id = 1;
    }
    if (opts_vals.count("stop-on-cc-err") > 0) {
		Global::set_stop_on_cc_err(true);
    } else {
		Global::set_stop_on_cc_err(false);
    }
    if (opts_vals.count("n_threads") > 0) {
        unsigned int n_threads = opts_vals["n_threads"].as<uint>();
        Global::set_n_threads( n_threads );
        if (n_threads == 1) {
            cerr << "Warning: Defining only 1 working thread is useless, as the main thread will wait until all workers are finished.\nPlease increase the number of working threads or disable multi-threading by setting n_threads to 0." << std::endl;
        }
    }
    if (opts_vals.count("work-stealing") > 0) {
        Global::set_work_stealing(true);
    }
    if (opts_vals.count("max-parallel-opti-vars") > 0) {
        unsigned long n_vars = opts_vals["max-parallel-opti-vars"].as<unsigned long>();
        Global::set_max_parallel_opti_vars( n_vars );
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
    if (opts_vals.count("ccmd-output") > 0) {
        string st_output = opts_vals["ccmd-output"].as<string>();
        if (st_output == "no" || st_output == "off") {
            Global::set_create_control_cmd_output(false);
        } else if (st_output == "all" || st_output == "on") {
            Global::set_create_control_cmd_output(true);
        } else {
            // invalid argument
            cerr << "Error when parsing command line arguments: invalid option for --ccmd-output given!" << endl;
            return 1;
        }
    } else {
        Global::set_create_control_cmd_output(false);
    }
    if (opts_vals.count("ev-output") > 0) {
        string st_output = opts_vals["ev-output"].as<string>();
        if (st_output == "no" || st_output == "off") {
            Global::set_create_ev_detailed_output(false);
        } else if (st_output == "all") {
            Global::set_create_ev_detailed_output(true);
        } else {
            // invalid argument
            cerr << "Error when parsing command line arguments: invalid option for --ev-output given!" << endl;
            return 1;
        }
    } else {
        Global::set_create_ev_detailed_output(false);
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
    if (opts_vals.count("cu-output-selection") > 0) {
        std::string selected_CUs_str = opts_vals["cu-output-selection"].as<string>();
        // parse the string
        std::stringstream selected_CUs_strstr( selected_CUs_str );
        std::string oneStrSection;
        while (std::getline(selected_CUs_strstr, oneStrSection, ',')) {
            try {
                global::unitIDs_selected_for_output.push_back( std::stoul( oneStrSection ) );
            } catch (const std::exception& e) {
                std::cerr << "Error when parsing values of command line argument --cu-output-selection for token '" << oneStrSection << "'." << std::endl;
                return 1;
            }
        }
    }

    // get time for time measurement
    auto t1 = std::chrono::system_clock::now();
    global::time_of_simulation_start = t1;
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
    // Set cache file path
    //
    HPProfileIDCache::GetInstance().setCacheFilename(Global::get_input_path() + "cache_file_HP_profiles.json");
    PVProfileIDCache::GetInstance().setCacheFilename(Global::get_input_path() + "cache_file_PV_profiles.json");

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
    // Output all variable values (first time to stdout / cout)
    //
    configld::output_variable_values(std::cout);

    // get time for time measurement
    auto t2 = std::chrono::system_clock::now();

    //
    // Plan and add which sim. added components should be added to which control units
    // and Run the simulation
    // - once (if no parameter variation is selected) or
    // - multiple times, if param. vari. is selected
    //
    StatusOutput::start_status_updater_thread();
    if (!simulation::runCompleteSimulation(expansion_matrix_rel_freq, expansion_matrix_abs_freq, scenario_id)) {
        cerr << "Error during simulation run!" << endl;
        return 3;
    }
    StatusOutput::stop_status_updater_thread();

    // get time for time measurement and send first values to the file
    auto t3 = std::chrono::system_clock::now();
    long s_setup = std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count();
    long s_main  = std::chrono::duration_cast<std::chrono::seconds>(t3-t2).count();
    output::outputRuntimeInformation(s_setup, s_main);

    //
    // Output all variable values (second time to a file)
    //
    std::ofstream log_file(*global::current_global_output_dir / "parameter-settings-general.txt");
    configld::output_variable_values(log_file);
    log_file.close();
	
	//
	// clean up
	//
    HPProfileIDCache::GetInstance().saveCacheFile();
    PVProfileIDCache::GetInstance().saveCacheFile();
    EVFSM::VacuumStaticVariables();
    ComponentHP::VacuumStaticVariables();
	MeasurementUnit::VacuumInstancesAndStaticVariables();
	ControlUnit::VacuumInstancesAndStaticVariables();
	Substation::VacuumInstancesAndStaticVariables();
	Global::DeleteStaticVariables();
	global::vacuum();

    // get time for time measurement
    auto t4 = std::chrono::system_clock::now();
    cout << "Run-time information:\n";
    cout << "  Setup and data loading: " << s_setup << "s\n";
    cout << "  Main run:               " << s_main  << "s\n";
    cout << "  Clean up:               " << std::chrono::duration_cast<std::chrono::seconds>(t4-t3).count() << "s\n";
    cout << "  Complete run time:      " << std::chrono::duration_cast<std::chrono::seconds>(t4-t1).count() << "s" << std::endl;

#endif

	return 0;
}

#ifdef PYTHON_MODULE

using namespace pybind11::literals;

// Module definition
PYBIND11_MODULE(UrbanFireXDT, m) {
    m.doc() = "Python bindings for the UrbanFireDTX C++ simulation.";

    pybind11::class_<pyconn::SimulationControlUnitState>(m, "SimulationControlUnitState",
        "State of a control unit during a simulation step.")
        .def_readonly("controlUnitID", &pyconn::SimulationControlUnitState::controlUnitID,
                      "Public ID of the control unit.")
        .def_readonly("internalID", &pyconn::SimulationControlUnitState::internalID,
                      "Internal ID of the control unit.")
        .def_readonly("has_pv", &pyconn::SimulationControlUnitState::has_pv,
                      "True if a PV system (simulated or in data) is present.")
        .def_readonly("has_cntrl_bs", &pyconn::SimulationControlUnitState::has_cntrl_bs,
                      "True if a controllable (i.e., simulated) battery storage is present.")
        .def_readonly("has_cntrl_hp", &pyconn::SimulationControlUnitState::has_cntrl_hp,
                      "True if a controllable (i.e., simulated) heat pump is present.")
        .def_readonly("has_cntrl_evchst", &pyconn::SimulationControlUnitState::has_cntrl_evchst,
                      "True if a controllable (i.e., simulated) EV charging station is present.");

    m.def("initialize", &pyconn::initialize_simulation,
          "Initialize the simulation with a Python dict of options.");
    m.def("get_next_timestep_id", &pyconn::get_next_timestep_id,
         "Returns the timestep ID of the next time step.");
    m.def("simulation_finished", &pyconn::simulation_finished,
         "Returns true, if the simulation is finished, otherwise false. In the latter case, run_step() can be called.");
    m.def("get_state", &pyconn::getState, "controlUnitID"_a,
          "Get the state of a control unit.");
    m.def("send_commands", &pyconn::sendCommands, "controlUnitID"_a, "commands"_a,
          "Send a command to a control unit.");
    m.def("run_one_step", &pyconn::run_one_step,
          "Advance the simulation by one time step.");
    m.def("reset_simulation_to_start", &pyconn::reset_simulation_to_start,
          "Resets the simulation and all internal states to the first time step.");
    m.def("vacuum", &pyconn::vacuum,
          "Clean up and release resources.");
}

#endif
