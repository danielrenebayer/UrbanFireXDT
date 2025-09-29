/**
 * @file python_module.hpp
 * @brief Definition and implementation of functions that connect the C++ simulation with Python via pybind11.
 *
 * This is a header-only file providing Python bindings to initialize, control,
 * and finalize simulation runs from within Python.
 */

#ifndef PYTHON_MODULE_HPP
#define PYTHON_MODULE_HPP

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // for converting C++ containers like std::map

#include "cache_helper.hpp"
#include "components.h"
#include "global.h"
#include "setup_and_dataloading.h"
#include "simulation_logic.h"
#include "units.h"

namespace pyconn {

    /**
     * @struct SimulationControlUnitState
     * @brief Represents the state of a single control unit during a simulation time step.
     *
     * This object is returned to Python to provide information about the
     * simulated components connected to a specific control unit.
     */
    struct SimulationControlUnitState {
        unsigned long controlUnitID;
        unsigned long internalID;
        bool has_pv;
        bool has_cntrl_bs;
        bool has_cntrl_hp;
        bool has_cntrl_evchst;
        unsigned long timestepID;
        // BS
        double bs_soc;
        double bs_soe;
        // HP
        // double hp_cumulative_energy_kWh;
        // double temp_indoor;
        // double hp_cop;
        // // EV 1, 2, 3
        // // SOC, SOE, Whether EV is connected, time until departure, required energy until departure, time until next arrival
        // // environmental data
        // double temp_outdoor;
        // double pv_generation;
        // double household_demand;
        // double electricity_price;
    };

    // ---------------------------
    // Thread & shared state
    // ---------------------------
    inline std::thread g_sim_worker; ///< The worker thread for the simulation that is started by pyconn::initialize_simulation()
    inline std::atomic<bool> g_sim_started{false};   ///< States if the working thread for the simulation has been started
    inline std::atomic<bool> g_sim_finished{false};  ///< States if the simulation has reached its final state (and cannot proceed anymore, or must be reseted to proceed)
    
    inline std::atomic<bool> g_simulation_main_part_started; ///< Flag indicating whether the simulation has finished the SAC planning and addition.
    inline std::atomic<bool> g_simulation_idling;            ///< Flag indicating whether the simulation is currently idling or working
    inline std::mutex g_mtx_simulation;                      ///< Mutex protecting synchronization between the simulation thread and Python calls.
    inline std::condition_variable g_cv_sim_state_update;    ///< Condition variable to notify and to wait for changes in the simulation / python code state.
    inline std::atomic<bool> g_next_step_requested;          ///< Flag indicating that a new simulation step should be executed
    inline std::atomic<bool> g_worker_threads_shutdown_cmd{false}; ///< Flag indicating that all worker threads should shoutdown and the simulation should be finished. After that, no reset is possible anymore.
    inline std::atomic<unsigned long> g_atomic_next_tsID;    ///< Shared variable indicating the next time step ID to be executed when calling pyconn::run_one_step()
    inline std::atomic<bool> g_restart_simulation{false};          ///< Flag indicating to restart the simulation after a reset of all internal states.
    // keep these alive for the worker's entire lifetime
    inline float         expansion_matrix_rel_freq[16][16] = {0}; ///< Variable storing the expansion matrix with relative frequencies - if using as stand-alone simulation, this is placed inside main and no global variable
    inline unsigned long expansion_matrix_abs_freq[16][16] = {0}; ///< Variable storing the expansion matrix with absolute frequencies - if using as stand-alone simulation, this is placed inside main and no global variable
    inline unsigned long scenario_id; ///< Variable storing the selected scenario ID - if using as stand-alone simulation, this is placed inside main and no global variable

    /**
     * @brief Initializes and starts the simulation from Python.
     *
     * This function sets up the simulation environment based on a dictionary of options,
     * equivalent to command line arguments used in the standalone executable.
     * Supported options include:
     * - "config" (string): Path to JSON configuration file (default: "../config/simulation_config.json")
     * - "scenario" (unsigned long): Scenario ID (default: 1)
     * - "n_threads" (unsigned int): Number of worker threads
     * - "work-stealing" (bool): Enable/disable work stealing
     * - "max-parallel-opti-vars" (unsigned long): Max. number of optimization vars processed in parallel
     * - "stop-on-cc-err" (bool): Stop on optimization errors
     * - "seed" (unsigned int): Seed for random number generation
     * - "weekly-metrics" (string): "on"/"off" to enable weekly metrics
     * - "cu-output", "st-output", "ccmd-output", "ev-output" (string): Output modes
     * - "cu-output-selection" (string): Comma-separated list of control unit IDs for output selection
     * - "suof" (unsigned long): Steps until output is flushed (default: 1000)
     *
     * @param args A Python dictionary containing simulation options.
     * @throws std::invalid_argument If invalid options are provided.
     * @throws std::runtime_error If required files cannot be loaded or variables remain uninitialized.
     */
    inline void initialize_simulation(pybind11::dict args) {
        std::string config_filepath;

        if (g_sim_started.load()) {
            throw std::runtime_error("Simulation already running. Call vacuum() before re-initializing.");
        }

        // --------------------------
        // Equivalent of command line parsing
        // --------------------------

        if (args.contains("config")) {
            config_filepath = args["config"].cast<std::string>();
        } else {
            config_filepath = "../config/simulation_config.json";
        }

        if (args.contains("pvar")) {
            //Global::set_pvar_vals(true, args["pvar"].cast<int>());
            throw std::runtime_error("Parameter variations cannot be defined when using the python interface!");
        } else {
            Global::set_pvar_vals(false, 0);
        }

        if (args.contains("scenario")) {
            scenario_id = args["scenario"].cast<unsigned long>();
        } else {
            scenario_id = 1;
        }

        if (args.contains("stop-on-cc-err") && args["stop-on-cc-err"].cast<bool>()) {
            Global::set_stop_on_cc_err(true);
        } else {
            Global::set_stop_on_cc_err(false);
        }

        if (args.contains("n_threads")) {
            unsigned int n_threads = args["n_threads"].cast<uint>();
            Global::set_n_threads(n_threads);
            if (n_threads == 1) {
                std::cerr << "Warning: Defining only 1 working thread is useless.\n";
            }
        } else {
            Global::set_n_threads(3);
        }

        if (args.contains("work-stealing") && args["work-stealing"].cast<bool>()) {
            Global::set_work_stealing(true);
        }

        if (args.contains("max-parallel-opti-vars")) {
            unsigned long n_vars = args["max-parallel-opti-vars"].cast<unsigned long>();
            Global::set_max_parallel_opti_vars(n_vars);
        }

        if (args.contains("cu-output")) {
            std::string cu_output = args["cu-output"].cast<std::string>();
            if (cu_output == "no" || cu_output == "off") {
                Global::set_output_mode_per_cu(global::OutputModePerCU::NoOutput);
            } else if (cu_output == "single") {
                Global::set_output_mode_per_cu(global::OutputModePerCU::SingleFile);
            } else if (cu_output == "sl") {
                Global::set_output_mode_per_cu(global::OutputModePerCU::IndividualFile);
            } else {
                throw std::invalid_argument("Invalid option for cu-output");
            }
        } else {
            Global::set_output_mode_per_cu(global::OutputModePerCU::IndividualFile);
        }

        if (args.contains("st-output")) {
            std::string st_output = args["st-output"].cast<std::string>();
            if (st_output == "no" || st_output == "off") {
                Global::set_create_substation_output(false);
            } else if (st_output == "on") {
                Global::set_create_substation_output(true);
            } else {
                throw std::invalid_argument("Invalid option for st-output");
            }
        } else {
            Global::set_create_substation_output(true);
        }

        if (args.contains("ccmd-output")) {
            std::string val = args["ccmd-output"].cast<std::string>();
            if (val == "no" || val == "off") {
                Global::set_create_control_cmd_output(false);
            } else if (val == "all" || val == "on") {
                Global::set_create_control_cmd_output(true);
            } else {
                throw std::invalid_argument("Invalid option for ccmd-output");
            }
        } else {
            Global::set_create_control_cmd_output(false);
        }

        if (args.contains("ev-output")) {
            std::string val = args["ev-output"].cast<std::string>();
            if (val == "no" || val == "off") {
                Global::set_create_ev_detailed_output(false);
            } else if (val == "all") {
                Global::set_create_ev_detailed_output(true);
            } else {
                throw std::invalid_argument("Invalid option for ev-output");
            }
        } else {
            Global::set_create_ev_detailed_output(false);
        }

        if (args.contains("repetitions")) {
            throw std::runtime_error("Repetitions cannot be defined when using the python interface");
        } else {
            Global::set_repetitions_selected(false);
        }

        if (args.contains("seed")) {
            unsigned int seed = args["seed"].cast<unsigned int>();
            Global::set_seed(seed);
            EVFSM::SetSeed(seed);
        }

        if (args.contains("weekly-metrics")) {
            std::string wm_mode = args["weekly-metrics"].cast<std::string>();
            if (wm_mode == "off" || wm_mode == "no") {
                Global::set_compute_weekly_metrics(false);
            } else if (wm_mode == "on" || wm_mode == "yes") {
                Global::set_compute_weekly_metrics(true);
            } else {
                throw std::invalid_argument("Invalid option for weekly-metrics");
            }
        }

        if (args.contains("cu-output-selection")) {
            std::string selected_CUs_str = args["cu-output-selection"].cast<std::string>();
            std::stringstream strstr(selected_CUs_str);
            std::string token;
            while (std::getline(strstr, token, ',')) {
                try {
                    global::unitIDs_selected_for_output.push_back(std::stoul(token));
                } catch (...) {
                    throw std::invalid_argument("Invalid cu-output-selection token: " + token);
                }
            }
        }

        // Flush interval
        if (args.contains("suof")) {
            global::n_ts_between_flushs = args["suof"].cast<unsigned long>();
        } else {
            global::n_ts_between_flushs = 1000;
        }

        std::cout << "Initializing the simulation for scenario ID " << scenario_id << std::endl;

        // Initialize statics
        Global::InitializeStaticVariables();

        // Load configuration
        if (!configld::load_config_file(scenario_id, config_filepath)) {
            throw std::runtime_error("Error loading config file.");
        }
        if (!expansion::load_expansion_matrix(expansion_matrix_rel_freq)) {
            throw std::runtime_error("Error loading expansion matrix.");
        }
        if (!expansion::verify_expansion_matrix(expansion_matrix_rel_freq)) {
            throw std::runtime_error("Invalid expansion matrix.");
        }

        HPProfileIDCache::GetInstance().setCacheFilename(Global::get_input_path() + "cache_file_HP_profiles.json");
        PVProfileIDCache::GetInstance().setCacheFilename(Global::get_input_path() + "cache_file_PV_profiles.json");

        if (!configld::load_data_from_central_database((Global::get_input_path() + Global::get_structure_database_name()).c_str())) {
            throw std::runtime_error("Error loading database.");
        }

        if (!global::all_variables_initialized() || !Global::AllVariablesInitialized()) {
            std::cerr << "\nThese variables are NOT initialized:" << std::endl;
            Global::PrintUninitializedVariables();
            global::print_uninitialized_variables();
            std::cerr << "" << std::endl;
            throw std::runtime_error("Some global variables are not initialized!");
        }

        configld::output_variable_values(std::cout);

        // ---------------------------
        // Launch the simulation worker thread
        // ---------------------------
        g_sim_started.store(true);
        g_sim_finished.store(false);
        g_simulation_idling            = false;
        g_next_step_requested          = false; // stop before executing first step
        g_simulation_main_part_started = false; // for SAC planning, no interactive mode is given -> sim. should not wait there!
        g_restart_simulation           = false;

        g_sim_worker = std::thread([](){
            try {
                const bool retval = simulation::runCompleteSimulation(
                    expansion_matrix_rel_freq,
                    expansion_matrix_abs_freq,
                    scenario_id
                );
                if (!retval) {
                    // Throwing an error kills all other processes
                    throw std::runtime_error("Error during simulation run!");
                }
            } catch (const std::exception& e) {
                std::cerr << "Exception in simulation worker: " << e.what() << std::endl;
            } catch (...) {
                std::cerr << "Unknown exception in simulation worker." << std::endl;
            }
        });
    }

    /**
     * @brief Returns the timestep ID of the next time step.
     *
     * @return timestep ID of the next time step.
     */
    inline unsigned long get_next_timestep_id() {
        if (!g_sim_started) {
            throw std::runtime_error("Error: Simulation has not been started!");
        }

        return g_atomic_next_tsID;
    }

    /**
     * @brief Checks if the simulation reached its end.
     *
     * If the simulation can proceed, the function run_one_step() can be executed.
     *
     * @return True, if the simulation can proceed and false, if the simulation reached its final state.
     */
    inline bool simulation_finished() {
        // Fast path: if already finished, return immediately
        if (g_sim_finished.load()) return true;
        // Otherwise wait until the sim is no longer idling (or it finishes while waiting)
        std::unique_lock<std::mutex> lock_obj(g_mtx_simulation);
        g_cv_sim_state_update.wait(lock_obj, [] {
            return !g_next_step_requested.load() && ( g_simulation_idling.load() || g_sim_finished.load() );
        });
        return g_sim_finished.load();
    }

    /**
     * @brief Retrieves the state of a control unit.
     *
     * Returns information about the simulated components
     * associated with the given control unit ID.
     *
     * @param controlUnitID The external (public) ID of the control unit.
     * @return A SimulationControlUnitState struct containing the unit’s state.
     */
    inline SimulationControlUnitState getState(unsigned long controlUnitID) {
        // Wait until simulation is idling (or finished)
        std::unique_lock<std::mutex> lock_obj(g_mtx_simulation);
        g_cv_sim_state_update.wait(lock_obj, [] {
            return g_simulation_idling.load() || g_sim_finished.load();
        });

        SimulationControlUnitState s{};
        ControlUnit* cu = ControlUnit::GetInstancePublicIDWE( controlUnitID );
        // No check for validity of the pointer `cu` required, as otherwise an exception is thrown
        s.controlUnitID    = controlUnitID;
        s.internalID       = cu->get_internal_id();
        s.has_pv           = cu->has_pv();
        s.has_cntrl_bs     = cu->has_bs();
        s.has_cntrl_hp     = cu->has_hp();
        s.has_cntrl_evchst = cu->has_cs();
        s.timestepID       = get_next_timestep_id();
        // BS
        s.bs_soc           = cu->has_bs() ? cu->get_component_BS()->get_SOC() : -1.0;
        s.bs_soe           = cu->has_bs() ? cu->get_component_BS()->get_SOE() : -1.0;
        // HP
        //s.hp_cumulative_energy_kWh = cu->has_hp() ? cu->get_component_HP()->get_currentDemand_kW() : -1.0; //TODO: subtract actual demand
        //s.temp_indoor      = 
        //s.hp_cop           = 
        // EV 1, 2, 3
        // SOC, SOE, Whether EV is connected, time until departure, required energy until departure, time until next arrival
        // environmental data
        //s.temp_outdoor;
        //s.pv_generation    = cu->has_pv() ? cu->get_current_PV_generation_kW() : -1.0;
        //s.household_demand = cu->get_current_load_vSMeter_kW();
        //s.electricity_price = 
        return s;
    }

    /**
     * @brief Sends a control command to a given control unit.
     *
     * This function can be used from Python to override or inject
     * control signals into the simulation during runtime at the
     * level of the control units (see class ControlUnit for
     * details).
     *
     * @param controlUnitID The external (public) ID of the control unit.
     * @param command The control value/command to be applied.
     */
    inline void sendCommands(unsigned long controlUnitID, pybind11::dict commands) {
        if (!g_sim_started) {
            throw std::runtime_error("Error: Simulation has not been started!");
        }

        // Wait until simulation is idling (or finished)
        std::unique_lock<std::mutex> lock_obj(g_mtx_simulation);
        g_cv_sim_state_update.wait(lock_obj, [] {
            return !g_next_step_requested.load() && ( g_simulation_idling.load() || g_sim_finished.load() );
        });

        ControlUnit* cu = ControlUnit::GetInstancePublicIDWE( controlUnitID );
        // process the commands
        /* if (commands.contains("command_name")) {
            config_filepath = commands["command_name"].cast<unsigned long>();
        }
        */
    }

    /**
     * @brief Executes a single simulation time step.
     *
     * Advances the simulation by one time step and halts once
     * new control commands can be accepted from Python.
     *
     * This function returns control back to Python once the step
     * is complete and the system is ready for new interactions.
     */
    inline void run_one_step() {
        //std::cout << "    1) [ID = " << g_atomic_next_tsID << "] g_simulation_idling = " << g_simulation_idling << " -- g_sim_finished = " << g_sim_finished.load() << std::endl;
        if (!g_sim_started) {
            throw std::runtime_error("Error: Simulation has not been started!");
        }
        if (g_sim_finished) {
            throw std::runtime_error("Error: Simulation already finished!");
        }

        {
            std::unique_lock<std::mutex> lock_obj(pyconn::g_mtx_simulation);
            // Wait until the simulation consumed the flag (set it back to false)
            g_cv_sim_state_update.wait(lock_obj, [] {
                return g_simulation_idling.load() || g_sim_finished.load();
            });
            // update flag and notify all
            g_next_step_requested = true;
            lock_obj.unlock();
            g_cv_sim_state_update.notify_all();
            lock_obj.lock();
            // Wait (again) until the simulation consumed the flag (set it back to false)
            // only required by this function as it blocks until the simulation step is finished
            g_cv_sim_state_update.wait(lock_obj, [] {
                return !g_next_step_requested.load() && ( g_simulation_idling.load() || g_sim_finished.load() );
            });
        }
        //std::cout << "    2) [ID = " << g_atomic_next_tsID << "] g_simulation_idling = " << g_simulation_idling << " -- g_sim_finished = " << g_sim_finished.load() << std::endl;
    }

    /**
     * @brief Resets the simulation and all internal states to the first time step.
     *
     * All information and collected metrics inside the control units are deleted.
     */
    inline void reset_simulation_to_start() {
        if (!g_sim_started) {
            throw std::runtime_error("Error: Simulation has not been started!");
        }

        {
            // Wait until simulation is idling (or finished)
            std::unique_lock<std::mutex> lock_obj(g_mtx_simulation);
            g_cv_sim_state_update.wait(lock_obj, [] {
                return g_simulation_idling.load() || g_sim_finished.load();
            });

            ControlUnit::ResetAllInternalStates();
            g_sim_finished       = false;
            g_simulation_idling  = false;
            g_restart_simulation = true;
        }
        g_cv_sim_state_update.notify_all();
    }

    /**
     * @brief Cleans up the simulation, finished the working thread and releases resources.
     *
     * - Saves cache files for PV/HP profiles
     * - Releases all allocated static and dynamic simulation objects
     * - Resets global variables
     * - Joins the working thread
     *
     * This function should be called once a simulation is finished.
     */
    inline void vacuum() {
        if (!g_sim_started) {
            throw std::runtime_error("Error: Simulation has not been started!");
        }
        if (!g_sim_finished) {
            throw std::runtime_error("Error: Simulation has not been finished!");
        }

        // Send shutdown command
        g_worker_threads_shutdown_cmd = true;
        pyconn::g_cv_sim_state_update.notify_all();

        // Join worker thread if still running
        if (g_sim_worker.joinable()) {
            try {
                g_sim_worker.join();
            } catch (...) {
                // avoid throwing from vacuum; just ensure join attempt
            }
        }

        // cleanup like in main()
        HPProfileIDCache::GetInstance().saveCacheFile();
        PVProfileIDCache::GetInstance().saveCacheFile();
        EVFSM::VacuumStaticVariables();
        ComponentHP::VacuumStaticVariables();
        MeasurementUnit::VacuumInstancesAndStaticVariables();
        ControlUnit::VacuumInstancesAndStaticVariables();
        Substation::VacuumInstancesAndStaticVariables();
        Global::DeleteStaticVariables();
        global::vacuum();

        g_sim_started.store(false);
        g_sim_finished.store(false);
    }

}

#endif
