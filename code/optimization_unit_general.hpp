/**
 * optimization_unit_general.hpp
 *
 * This file contains all general classes / structs required
 * by the optimization implementations.
 */

#ifndef OPTIMIZATION_UNIT_GENERAL_HPP
#define OPTIMIZATION_UNIT_GENERAL_HPP

#include <vector>

/**
 * This class represents the base class for an optimized (i.e., not rule-based) controller.
 */
class BaseOptimizedController {
    protected:
        const unsigned long controlUnitID; ///< ID of the connected control unit
        const unsigned int n_cars; ///< Number of EVs with the given control unit as home place
        unsigned long time_horizon;
        std::vector<double> bs_power; ///< The battery storage power for every time step in the time horizon
        std::vector<double> hp_power; ///< The heat pump power for every time step in the time horizon
        std::vector<std::vector<double>> ev_power; ///< The charging power per EV (index 0) for every time step (index 1) in the time horizon

    public:
        /**
         * Constructs a new object with 0.0-initialized vectors with length of parameter time_horizon
         */
        BaseOptimizedController(unsigned long cuID, unsigned int time_horizon, unsigned int n_cars) :
            controlUnitID(cuID), n_cars(n_cars), bs_power(time_horizon, 0.0), hp_power(time_horizon, 0.0)
        {
            this->time_horizon = time_horizon;
            ev_power.assign(n_cars, std::vector<double>(time_horizon, 0.0));
        }

        virtual ~BaseOptimizedController() = default;

        /**
         * Returns the optimized control for the heat pump for the next time steps in the time_horizon
         * Attention: The returned object will be overwritten in future calls!
         */
        const std::vector<double>& get_future_hp_power_kW() { return hp_power; }

        /**
         * Returns the optimized control for the battery storage for the next time steps in the time_horizon
         * Attention: The returned object will be overwritten in future calls!
         */
        const std::vector<double>& get_future_bs_power_kW() { return bs_power; }

        /**
         * Returns the optimized control for every EV for the next time steps in the time_horizon
         * Attention: The returned object will be overwritten in future calls!
         */
        const std::vector<std::vector<double>>& get_future_ev_power_kW() { return ev_power; }

            /**
             * Resets the vectors with 0.0-initialized values with length of parameter (new) parameter time_horizon
             */
            void reset(unsigned int new_horizon) {
                this->time_horizon = new_horizon;
                bs_power.assign(new_horizon, 0.0);
                hp_power.assign(new_horizon, 0.0);
                ev_power.assign(this->n_cars, std::vector<double>(new_horizon, 0.0));
            }

            /**
             * Removes the first element, shifts the second one to the first place and so on and
             * adds a 0.0 at the last place.
             */
            void shiftVectorsByOnePlace() {
                std::move(bs_power.begin() + 1, bs_power.end(), bs_power.begin());
                bs_power.back() = 0.0;
                std::move(hp_power.begin() + 1, hp_power.end(), hp_power.begin());
                hp_power.back() = 0.0;
                for (unsigned int evIdx = 0; evIdx < n_cars; evIdx++) {
                    std::move(ev_power[evIdx].begin() + 1, ev_power[evIdx].end(), ev_power[evIdx].begin());
                    ev_power[evIdx].back() = 0.0;
                }
            }

        /**
         * Executes the controller with all (future) states in the current horizon and stores the results the member variables BaseOptimizedController::bs_power, BaseOptimizedController::hp_power and BaseOptimizedController::ev_power.
         *
         * @param ts: The current time step ID
         * @param max_p_bs_kW: Maximum battery power in kW
         * @param max_e_bs_kWh: Battery capacity in kWh
         * @param max_p_hp_kW: Maximum heat pump power in kW (only shiftable part!)
         * @param max_p_cs_kW: Maximum charging station power in kW
         * @param current_bs_charge_kWh: Current charge of the battery in kWh
         * @param future_resid_demand_kW: Vector with future residential energy demand in kW per time step in the horizon
         * @param future_pv_generation_kW: Vector with future PV generation power in kW per time step in the horizon
         * @param future_hp_unshiftable_kW: Vector with future unshiftable heat pump demand in kW per time step in the horizon
         * @param future_hp_shiftable_maxE: Vector with maximum accumulated energy consumption of the heat pump (without shiftable part) in kWh until the end of each time step in the horizon
         * @param future_hp_shiftable_minE: Vector with minimum accumulated energy consumption of the heat pump (without shiftable part) in kWh until the end of each time step in the horizon
         * @param future_ev_shiftable_maxE: Vector with maximum accumulated energy consumption of the charging station in kWh for every step over the considered horizon (first index) and every EV (second index)
         * @param future_ev_shiftable_minE: Vector with minimum accumulated energy consumption of the charging station in kWh for every step over the considered horizon (first index) and every EV (second index)
         * @param future_ev_maxP: Vector with maximum charging power per EV in kW for every step over the considered horizon (first index) and every EV (second index)
         *
         * @return: Returns a boolean value, with true indicating if the optimization succeded, or false indicating an error (e.g., an unfeasable problem)
         */
        virtual bool updateController(
            unsigned long ts,
            float max_p_bs_kW,
            float max_e_bs_kWh,
            float max_p_hp_kW,
            float max_p_cs_kW,
            float current_bs_charge_kWh,
            const std::vector<float>& future_resid_demand_kW,
            const std::vector<double>& future_pv_generation_kW,
            const std::vector<double>& future_hp_unshiftable_kW,
            const std::vector<double>& future_hp_shiftable_maxE,
            const std::vector<double>& future_hp_shiftable_minE,
            const std::vector<const std::vector<double>*>* future_ev_shiftable_maxE,
            const std::vector<const std::vector<double>*>* future_ev_shiftable_minE,
            const std::vector<const std::vector<double>*>* future_ev_maxP
        ) = 0;
};

#endif

