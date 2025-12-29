
#ifndef SURPLUS_CONTROLLER_HPP
#define SURPLUS_CONTROLLER_HPP

#include <unordered_map>
#include <vector>
#include <memory>
#include "worker_threads.hpp"

namespace surplus {

    /**
     * @brief Result structure for lookahead simulation
     */
    struct LookaheadResult {
        std::vector<double> total_load;  ///< Total load vector for each timestep in the horizon (in kWh)
        std::unordered_map<unsigned long, std::vector<double>> grid_demand_kWh;  ///< Grid demand per unit ID for all timesteps in horizon
        std::unordered_map<unsigned long, std::vector<double>> bs_stored_energy_kWh; ///< BESS stored energy per unit ID for all timesteps in horizon
        std::unordered_map<unsigned long, std::vector<double>> bs_power_kW; ///< BESS power per unit ID for all timesteps in horizon
    };

    /**
     * @brief Singleton controller for global surplus optimization across all control units
     * 
     * The SurplusController manages energy surplus optimization across all control units.
     * It implements the singleton pattern to ensure there's only one instance managing the global optimization process. 
     * The controller supplies charge requests for all units and performs periodic optimization to distribute surplus energy.
     */
    class SurplusController {
    private:
        static std::unique_ptr<SurplusController> instance; ///< Singleton instance pointer
        static bool initialized; ///< Flag indicating if singleton has been initialized
        
        // Data storage
        std::unordered_map<unsigned long, std::vector<double>> unit_charge_requests; ///< Charge requests per unit ID for all timesteps
        unsigned long last_optimization_ts;     ///< Last timestep when optimization was executed
        unsigned int optimization_frequency_ts; ///< Frequency of optimization in timesteps 
        unsigned int lookahead_horizon_ts;      ///< Lookahead horizon in timesteps for optimization
        bool enabled;                          ///< Whether surplus controller is enabled
        bool bess_knowledge;                    ///< Whether the surplus controller has knowledge of the state of charge and current power for all control units for all timesteps in the horizon

        std::vector<double> future_surplus_log; ///< Log of future surplus for analysis, TODO: remove later
        
        // Simulation parameters required for running the optimization
        CUControllerThreadGroupManager* thread_manager; 
        const char* output_prefix; 
        
        /**
         * @brief Private constructor for singleton
         */
        SurplusController();
        
    public:
        // Singleton access
        /**
         * @brief Get the singleton instance of SurplusController
         * @return Reference to the singleton SurplusController instance. If no instance exists, it is created.
         */
        static SurplusController& GetInstance();

        static double GetFutureSurplusLog(unsigned long ts); ///< Get logged future surplus for a specific timestep in the horizon, just for analysis, TODO: remove later
        
        /**
         * @brief Initialize the singleton instance
         * 
         * Creates and initializes the singleton instance if it doesn't exist.
         * This method should be called once at the beginning of the simulation.
         * @param thread_manager Optional thread manager for parallel execution
         * @param output_prefix Optional output prefix for simulation messages
         */
        static void Initialize(CUControllerThreadGroupManager* thread_manager = nullptr, const char* output_prefix = "");
        
        /**
         * @brief Clean up and destroy the singleton instance
         * 
         * Safely destroys the singleton instance and releases all resources.
         * This method should be called at the end of the simulation.
         */
        static void Cleanup();
        
        // Main optimization interface
        /**
         * @brief Check if optimization should be run at the current timestep
         * @param current_ts The current simulation timestep
         * @return true if optimization should be executed, false otherwise
         * 
         * Determines whether optimization should run based on the optimization
         * frequency and the last timestep when optimization was executed.
         */
        bool ShouldRunOptimization(unsigned long current_ts) const;
        
        /**
         * @brief Execute the surplus optimization algorithm
         * @param ts_horizon_start The starting timestep for the optimization horizon
         * @return true if optimization was successful, false otherwise
         * 
         * Runs the global surplus optimization algorithm across all control units
         * from ts_horizon_start to the length of the specified horizon. Updates internal data structures with
         * the charge requests that units can access later.
         */
        bool ExecuteOptimization(unsigned long ts_horizon_start);
        
        // Data access for ControlUnits
        /**
         * @brief Get the optimized charge request for a specific control unit
         * @param unit_id The unitID of the control unit
         * @return The charge request value for the unit, or 0.0 if not found
         * 
         * Retrieves the current optimized charge request for the specified control unit.
         * This value is set during the optimization process and represents the amount
         * of surplus energy allocated to this unit at the time step the surplus controller knows as current one.
         */
        double GetChargeRequest(unsigned long unit_id) const;

        // Static convenience methods 
        /**
         * @brief Get the total surplus energy scheduled to batteries
         * @return The total surplus energy scheduled to batteries across all units in the current time step
         */
        static double GetScheduledSurplusToBESS();

        /**
         * @brief Get the total surplus energy actually allocated to batteries
         * @return The total surplus energy allocated to batteries across all units in the current time step
         */
        static double GetActualSurplusToBESS();

        // TODO: just for analysis, very inefficient, remove later
        static double GetBESSChargeRequest();
        static double GetBESSLoad();
        
        /**
         * @brief Static convenience method to get charge request for a unit
         * @param unit_id The unitID of the control unit
         * @return The charge request value for the unit, or 0.0 if not found
         * 
         * Static wrapper around GetChargeRequest() that automatically uses the
         * singleton instance. Provides convenient access for ControlUnit objects.
         */
        static double GetChargeRequestForUnit(unsigned long unit_id);
        
        // Internal methods
        /**
         * @brief Shift time series data to the next timestep
         * 
         * Advances all time-dependent data structures by one timestep.
         * This method should be called from outside at the beginning of each new timestep, so that the surplus controller keeps track of the current time step.
         */
        void ShiftTimeSeriesData();
        void ResetAllData();
        
    private:
        /**
         * @brief Run the simulation lookahead over the optimization horizon and return required data
         * @param ts_horizon_start The starting timestep for the surplus calculation horizon
         * @param ts_horizon_end The ending timestep for the surplus calculation horizon
         * @return LookaheadResult (see struct definition)
         * 
         * Computes the expected surplus energy available at all timesteps in the horizon
         * by pre-running the simulation for the horizon length. Also returns the grid demand
         * for each control unit and the BESS stored energy for every timestep in the horizon.
         */
        LookaheadResult LookaheadSimulation(unsigned long ts_horizon_start, unsigned long ts_horizon_end);
    };

} // namespace surplus

#endif // SURPLUS_CONTROLLER_HPP