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
 * @brief State structure for Optimization Unit.
 * 
 * This structure stores all internal state variables of the optimization unit of a control unit
 * that need to be preserved during save/restore operations.
 */
struct OptimizationUnitState {
    std::vector<double> bs_power; 
    std::vector<double> hp_power;
    std::vector<std::vector<double>> ev_power;
};


/**
 * This class represents the base class for an optimized (i.e., not rule-based) controller.
 */
class BaseOptimizedController {
    protected:
        const unsigned long controlUnitID; ///< ID of the connected control unit
        const unsigned long n_cars; ///< Number of EVs with the given control unit as home place
        unsigned long time_horizon;
        std::vector<double> bs_power; ///< The battery storage power for every time step in the time horizon
        std::vector<double> hp_power; ///< The heat pump power for every time step in the time horizon
        std::vector<std::vector<double>> ev_power; ///< The charging power per EV (index 0) for every time step (index 1) in the time horizon
        std::list<double> optimal_pv_size_per_section_kW; ///< If the PV size should be optimized, this value contains the optimal PV size per section in kW
        double optimal_bs_size_kWh; ///< If the BS size should be optimized, this value contains the optimal BS size in kWh

    public:
        /**
         * Constructs a new object with 0.0-initialized vectors with length of parameter time_horizon
         */
        BaseOptimizedController(unsigned long cuID, unsigned int time_horizon, unsigned long n_cars) :
            controlUnitID(cuID), n_cars(n_cars), bs_power(time_horizon, 0.0), hp_power(time_horizon, 0.0)
        {
            this->time_horizon = time_horizon;
            ev_power.assign(n_cars, std::vector<double>(time_horizon, 0.0));
            optimal_bs_size_kWh = 0.0;
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
         * Returns the optimal BS size in kWh if the BS size was also part of the optimization
         */
        const std::list<double>& get_optimal_PV_size_per_section_kW() { return optimal_pv_size_per_section_kW; }

        /**
         * Returns the optimal BS size in kWh if the BS size was also part of the optimization
         */
        double get_optimal_BS_size_kWh() { return optimal_bs_size_kWh; }

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
         * Saves the internal state of the optimization unit.
         * @return: Returns an OptimizationUnitState struct containing all relevant internal state variables.
         */
        OptimizationUnitState saveInternalState() const {
            OptimizationUnitState state;
            state.bs_power = bs_power;
            state.hp_power = hp_power;
            state.ev_power = ev_power;
            return state;
        }

        /**
         * Restores the internal state of the optimization unit from a previously saved state.
         * @param state: The OptimizationUnitState struct containing all relevant internal state variables.
         */
        void restoreInternalState(const OptimizationUnitState& state) {
            bs_power = state.bs_power;
            hp_power = state.hp_power;
            ev_power = state.ev_power;
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
         * @brief Executes the controller with all (future) states in the current horizon and stores the results the member variables BaseOptimizedController::bs_power, BaseOptimizedController::hp_power and BaseOptimizedController::ev_power.
         * 
         * If parameter optimize_PV_size is set to true, the value of the parameter future_pv_generation_kW is ignored, and parameters total_PV_generation_per_section_kW and must be a valid reference.
         * Otherwise, the last two parameters are ignored.
         * If parameter optimize_BS_size is set to true, the value of max_e_bs_kWh and max_p_cs_kW are ignored. They are a result of the optimization. Moreover, in this case, current_bs_charge_kWh must be set to 0.0.
         *
         * @param ts: The current time step ID
         * @param max_p_bs_kW: Maximum battery power in kW
         * @param max_e_bs_kWh: Battery capacity in kWh
         * @param max_p_cs_kW: Maximum charging station power in kW
         * @param current_bs_charge_kWh: Current charge of the battery in kWh
         * @param future_resid_demand_kW: Vector with future residential energy demand in kW per time step in the horizon
         * @param future_pv_generation_kW: Vector with future PV generation power in kW per time step in the horizon
         * @param future_hp_shiftable_maxP: Vector with maximum power of the heat pump per time step
         * @param future_hp_shiftable_minP: Vector with minimum power of the heat pump per time step
         * @param future_hp_shiftable_maxE: Vector with maximum accumulated energy consumption of the heat pump in kWh until the end of each time step in the horizon
         * @param future_hp_shiftable_minE: Vector with minimum accumulated energy consumption of the heat pump in kWh until the end of each time step in the horizon
         * @param future_ev_shiftable_maxE: Vector with maximum accumulated energy consumption of the charging station in kWh for every step over the considered horizon (first index) and every EV (second index)
         * @param future_ev_shiftable_minE: Vector with minimum accumulated energy consumption of the charging station in kWh for every step over the considered horizon (first index) and every EV (second index)
         * @param future_ev_maxP: Vector with maximum charging power per EV in kW for every step over the considered horizon (first index) and every EV (second index)
         * @param optimize_PV_size: Flag indicating if the PV size (per section) should be part of the optimization
         * @param optimize_BS_size: Flag indicating if the BS size should be part of the optimization
         * @param total_PV_generation_per_section_kW: Total PV generation over the simulation horizon. Only required, if PV sizing (flag optimize_PV_size) is also a result of the optimization - otherwise, this parameter can be NULL
         * @param max_PV_power_per_section_kWp: The maximum installable power per roof section in kWp. Only required, if PV sizing (flag optimize_PV_size) is also a result of the optimization - otherwise, this parameter can be NULL
         *
         * @return: Returns a boolean value, with true indicating if the optimization succeded, or false indicating an error (e.g., an unfeasable problem)
         */
        virtual bool updateController(
            const unsigned long ts,
            double max_p_bs_kW,
            double max_e_bs_kWh,
            double max_p_cs_kW,
            double current_bs_charge_kWh,
            const std::vector<float>& future_resid_demand_kW,
            const std::vector<double>& future_pv_generation_kW,
            const std::vector<double>& future_hp_shiftable_maxP,
            const std::vector<double>& future_hp_shiftable_minP,
            const std::vector<double>& future_hp_shiftable_maxE,
            const std::vector<double>& future_hp_shiftable_minE,
            const std::vector<const std::vector<double>*>* future_ev_shiftable_maxE,
            const std::vector<const std::vector<double>*>* future_ev_shiftable_minE,
            const std::vector<const std::vector<double>*>* future_ev_maxP,
            const bool optimize_PV_size,
            const bool optimize_BS_size,
            const std::list<std::vector<double>>* total_PV_generation_per_section_kW,
            const std::list<double>* max_PV_power_per_section_kWp
        ) = 0;
};

#endif

