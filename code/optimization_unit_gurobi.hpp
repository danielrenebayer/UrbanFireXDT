/**
 * optimization_unit_gurobi.hpp
 *
 * This file contains all functions required to call gurobi as
 * optimizer inside the control units.
 */

#ifndef OPTIMIZATION_UNIT_GUROBI_HPP
#define OPTIMIZATION_UNIT_GUROBI_HPP

#include "global.h"
#include "optimization_unit_general.hpp"

#include "gurobi_c++.h"

class GurobiLPController : public BaseOptimizedController {

    private:
        static GRBEnv* env; ///< The global gurobi environment

    public:
        GurobiLPController(unsigned long cuID, unsigned int time_horizon, unsigned long n_cars) :
            BaseOptimizedController(cuID, time_horizon, n_cars)
        {}

        using BaseOptimizedController::updateController;
        bool updateController(
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
            );

        /**
         * Initializes the global environment.
         */
        static void InitializeGurobiEnvironment();

        /**
         * Deletes all global variables.
         */
        static void VaccumAllStaticVariables();

};

#endif
