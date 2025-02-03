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
        GurobiLPController(unsigned long cuID, unsigned int time_horizon, unsigned int n_cars) :
            BaseOptimizedController(cuID, time_horizon, n_cars)
        {}

        using BaseOptimizedController::updateController;
        bool updateController(
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
