/**
 * optimization_unit_coinor.hpp
 *
 * This file contains all functions required to call the CoinOR as
 * optimizer inside the control units.
 */

#ifndef OPTIMIZATION_UNIT_COINOR_HPP
#define OPTIMIZATION_UNIT_COINOR_HPP

//#include <span>
#include <vector>

#include "coin/ClpSimplex.hpp"

namespace CUOptimization {

    namespace {
        inline unsigned long n_opti_vars = 0; ///< Number of optimization variables
        inline std::vector<double> objective; ///< The objective function
        inline std::vector<double> A; ///< row-wise matrix A for inqeuality constraints (b1_lower <= Ax <= b1_upper) (as a vector of len = A_rows * n_opti_vars of A)
        inline std::vector<double> B; ///< row-wise matrix B for   qeuality constraints (Bx == b2) (as a vector of len = B_rows * n_opti_vars of B)
        inline unsigned long A_rows = 0; // number of rows in matrix A
        inline unsigned long B_rows = 0; // number of rows in matrix B
        inline unsigned long v_pos_p_hp_start = 0; ///< First position of p_hp in the vector of variables
        inline unsigned long v_pos_e_hp_start = 0; ///< First position of e_hp (cumsum) in the vector of variables
        inline unsigned long v_pos_p_cs_start = 0; ///< First position of p_cs in the vector of variables
        inline unsigned long v_pos_e_cs_start = 0; ///< First position of e_cs (cumsum) in the vector of variables
        inline unsigned long v_pos_p_bs_in_start  = 0; ///< First position of p_bs_in in the vector of variables
        inline unsigned long v_pos_p_bs_out_start = 0; ///< First position of p_bs_in in the vector of variables
        inline unsigned long v_pos_x_demand_start = 0; ///< First position of x_demand in the vector of variables
        inline unsigned long v_pos_x_feedin_start = 0; ///< First position of x_demand in the vector of variables
        inline unsigned long v_pos_e_bs_start = 0; ///< First position of e_bs in the vector of variables
        inline unsigned long v_pos_M = 0; // only valid, if PeakLoad is optimization target -> Position of M in variable vector
        inline unsigned long b2_pos_curr_e_bs = 0; ///< Position of current bs storage in kWh in b2
        inline unsigned long b2_pos_curr_p_local_start = 0; ///< First position of p_resid-p_pv+p_hp(unshift) in kW in b2
        inline unsigned long b1_pos_hp_minE_maxE_start = 0; ///< First position of hp_maxE / minE constraint in b1
        inline unsigned long b1_pos_cs_minE_maxE_start = 0; ///< First position of cs_maxE / minE constraint in b1
    }

    /**
     * Initializes and executes the optimization for one control unit using the CoinOR optimizer (open source).
     *
     * @param ts: The current time step ID
     * @param outputRef: A reference to the output object
     * @param max_p_hp_kW: Maximum heat pump power in kW (only shiftable part!)
     * @param max_p_cs_kW: Maximum charging station power in kW
     * @param max_p_bs_kW: Maximum battery power in kW
     * @param max_e_bs_kWh: Battery capacity in kWh
     * @param future_resid_demand_kW: Vector with future residential energy demand in kW per time step in the horizon
     * @param future_pv_generation_kW: Vector with future PV power in kW per time step in the horizon
     * @param future_hp_unshiftable_kW: Vector with future unshiftable heat pump demand in kW per time step in the horizon
     * @param future_hp_shiftable_maxE: Vector with maximum accumulated energy consumption of the heat pump (without shiftable part) in kWh until the end of each time step in the horizon
     * @param future_hp_shiftable_minE: Vector with minimum accumulated energy consumption of the heat pump (without shiftable part) in kWh until the end of each time step in the horizon
     * @param future_cs_shiftable_maxE: Vector with maximum accumulated energy consumption of the charging station in kWh until the end of each time step in the horizon
     * @param future_cs_shiftable_minE: Vector with minimum accumulated energy consumption of the charging station in kWh until the end of each time step in the horizon
     * @param current_bs_charge_kWh: Current charge of the battery in kWh
     *
     * @return: Returns a boolean value, with true indicating if the optimization succeded, or false indicating an error (e.g., an unfeasable problem)
     */
    bool runOptimizationForOneCU_CLP(
        unsigned long ts,
        CUOptimization::OptimalControlCommandsOverHorizon& outputRef,
        float max_p_hp_kW,
        float max_p_cs_kW,
        float max_p_bs_kW,
        float max_e_bs_kWh,
        const std::vector<float>& future_resid_demand_kW,
        const std::vector<double>& future_pv_generation_kW,
        const std::vector<double>& future_hp_unshiftable_kW,
        const std::vector<double>& future_hp_shiftable_maxE,
        const std::vector<double>& future_hp_shiftable_minE,
        const std::vector<double>& future_cs_shiftable_maxE,
        const std::vector<double>& future_cs_shiftable_minE,
        float current_bs_charge_kWh
    )
    {
        // Define the optimization target for the case of min. electr. prices or emissions
        // Hint: target set correctly for the case of peak demand reduction in prepareOrUpdateCLPModel()
        std::vector<double> local_objective_copy = objective; // copy constructor
        if (Global::get_controller_optimization_target() == global::ControllerOptimizationTarget::ElectricityCosts) {
            for (unsigned int t = 0; t < T; t++) {
                float current_demand_tariff = Global::get_demand_tariff();
                if (global::eprices_local_ts != NULL && ts + t < Global::get_n_timesteps())
                    current_demand_tariff = global::eprices_local_ts[ts - 1 + t];
                local_objective_copy[v_pos_x_demand_start + t] =   current_demand_tariff;
                local_objective_copy[v_pos_x_feedin_start + t] = - Global::get_feed_in_tariff();
            }
        } else if (Global::get_controller_optimization_target() == global::ControllerOptimizationTarget::Emissions) {
            for (unsigned int t = 0; t < T; t++) {
                float current_emissions = Global::get_emissions_g_CO2eq_per_kWh();
                    if (global::emission_ts != NULL && ts + t < Global::get_n_timesteps())
                        current_emissions = global::emission_ts[ts - 1 + t];
                local_objective_copy[v_pos_x_demand_start + t] =   current_emissions;
                // TODO: What about feedin emissions? Ignore? Substract?
            }
        }

        // Compute the right-side vectors b1 and b2
        // b1
        std::vector<double> b1_lower(A_rows, 0.0);
        std::vector<double> b1_upper(A_rows, 0.0);
        // lower / upper energy 
        for (unsigned int t = 0; t < T; t++) {
            b1_lower[b1_pos_hp_minE_maxE_start + t] = future_hp_shiftable_minE[t];
            b1_upper[b1_pos_hp_minE_maxE_start + t] = future_hp_shiftable_maxE[t];
            b1_lower[b1_pos_cs_minE_maxE_start + t] = future_cs_shiftable_minE[t];
            b1_upper[b1_pos_cs_minE_maxE_start + t] = future_cs_shiftable_maxE[t];
        }
        // b2
        std::vector<double> b2(B_rows, 0.0);
        // set current BS SOC
        b2[b2_pos_curr_e_bs] = current_bs_charge_kWh;
        // set current energy balance
        for (unsigned int t = 0; t < T; t++) {
            double resid_minus_pv_kW = future_resid_demand_kW[t] + future_hp_unshiftable_kW[t] - future_pv_generation_kW[t]; // i.e. local balance at time step t
            b2[b2_pos_curr_p_local_start + t] = resid_minus_pv_kW;
        }

        // Define the upper and lower bounds for the optimization variables
        std::vector<double> v_lower_bounds(n_opti_vars, 0.0);
        std::vector<double> v_upper_bounds(n_opti_vars, COIN_DBL_MAX);
        for (unsigned int t = 0; t < T; t++) {
            v_upper_bounds[t + v_pos_p_hp_start] = max_p_hp_kW;
            v_upper_bounds[t + v_pos_p_cs_start] = max_p_cs_kW;

            v_upper_bounds[t + v_pos_p_bs_in_start]  = max_p_bs_kW;
            v_upper_bounds[t + v_pos_p_bs_out_start] = max_p_bs_kW;
            /*
            v_upper_bounds[t + v_pos_e_hp_start] = future_hp_shiftable_maxE[t];
            v_lower_bounds[t + v_pos_e_hp_start] = future_hp_shiftable_minE[t];
            v_upper_bounds[t + v_pos_e_cs_start] = future_cs_shiftable_maxE[t];
            v_lower_bounds[t + v_pos_e_cs_start] = future_cs_shiftable_minE[t];
            */
        }
        for (unsigned int t = 0; t < T; t++) {
            v_upper_bounds[t + v_pos_e_bs_start] = max_e_bs_kWh;
        }
        if (v_pos_M > 0) {
            v_lower_bounds[v_pos_M] = 0.0;
            v_upper_bounds[v_pos_M] = COIN_DBL_MAX;
        }

        // Define the model
        ClpSimplex model;
        //CoinBigIndex rowStart[] = 
        model.loadProblem(n_opti_vars, A_rows + B_rows,
            columnStartIndices, // TODO !
            // see clp/examples/rowColumn.cpp
        );

        // Solve the model
        model.primal();

        // Check for optimal solution
        int model_status = model.status();
        if (model_status != 0) {
            std::cerr << "Optimization not resulting in optimal value.\n";
            std::cerr << "CoinOR CLP model status = " << model_status << std::endl;
            return false;
        }

        // Get the results and store them in the output reference
        const double* solution = model.primalColumnSolution();
        for (unsigned int t = 0; t < T; t++) {
            outputRef.hp_power[t] = solution[t];
            outputRef.cs_power[t] = p_cs_kW[t];
            outputRef.bs_power[t] = solution[t] - solution[t]; // TODO !!!
        }

        return true;
    }

    /**
     * Prepares the CoinOR model, i.e. defines the variables and the matricies.
     * If this function is called multiple times, it will reinitialize the variables and matricies
     * (required for example if the time horizon will change).
     */
    void prepareOrUpdateCLPModel() {
        // Parameters
        const unsigned int T  = Global::get_control_horizon_in_ts(); // number of time steps to consider
        const unsigned int Tp = T + 1; // required for the e_bs computation

        // Variable definition
        // - Define the first starting position per variable group
        v_pos_p_hp_start     = 0;
        v_pos_e_hp_start     = T * 1;
        v_pos_p_cs_start     = T * 2;
        v_pos_e_cs_start     = T * 3;
        v_pos_p_bs_in_start  = T * 4;
        v_pos_p_bs_out_start = T * 5;
        v_pos_x_demand_start = T * 6;
        v_pos_x_feedin_start = T * 7;
        v_pos_e_bs_start     = T * 8;
        v_pos_M = 0;
        // - Compute number of variables
        n_opti_vars = T * 8 + Tp;
        if (!Global::get_controller_allow_bs_grid_charging()) {
            n_opti_vars += T * 4;
        }
        if (Global::get_controller_optimization_target() == global::ControllerOptimizationTarget::PeakLoad) {
            v_pos_M = n_opti_vars;
            n_opti_vars += 1;
        }

        // Define the objective function for the case of peak load reduction
        objective.assign(n_opti_vars, 0.0); // set all values to 0.0 by default
        if (Global::get_controller_optimization_target() == global::ControllerOptimizationTarget::PeakLoad) {
            objective[v_pos_M] = 1.0;
        }

        // Define the system matrix
        A_rows = ;
        B_rows = T + 1 + T;
        A.assign(A_rows * n_opti_vars, 0.0);
        B.assign(B_rows * n_opti_vars, 0.0);
        // Battery balance equation -> Matrix B
        // -> n_opti_vars Equations = 0
        unsigned long cr = 0; // current row
        for (unsigned int t = 0; t < T; t++) {
            //std::span<double> this_row(A.data() + , );
            const unsigned long cri = cr * n_opti_vars; // current row ID
            B[cri + v_pos_e_bs_start     + t  ] =  1.0;
            B[cri + v_pos_e_bs_start     + t+1] = -1.0;
            B[cri + v_pos_p_bs_in_start  + t  ] =  Global::get_time_step_size_in_h() * Global::get_exp_bess_effi_in();
            B[cri + v_pos_p_bs_out_start + t  ] = -Global::get_time_step_size_in_h() * Global::get_exp_bess_effi_out();
            cr += 1;
        }
        // Battery initialization equation -> Matrix B
        // -> 1 euqation = current battery state
        B[cr * n_opti_vars + v_pos_e_bs_start] = 1.0;
        b2_pos_curr_e_bs = cr;
        cr += 1;
        // Power balance equation -> Matrix B
        // -> n_opti_vars Equations = ''local balance''
        b2_pos_curr_p_local_start = cr;
        if (Global::get_controller_allow_bs_grid_charging()) {
            // One balance equation
            for (unsigned int t = 0; t < T; t++) {
                const unsigned long cri = cr * n_opti_vars; // current row ID
                B[cri + v_pos_p_hp_start     + t] = -1.0;
                B[cri + v_pos_p_cs_start     + t] = -1.0;
                B[cri + v_pos_x_demand_start + t] =  1.0;
                B[cri + v_pos_x_feedin_start + t] = -1.0;
                B[cri + v_pos_p_bs_in_start  + t] = -1.0;
                B[cri + v_pos_p_bs_out_start + t] =  1.0;
                cr += 1;
            }
        } else {
            // Two balance euqations
            throw std::string("Not implemented!");
            // TODO: Mind also the right side!
        }
        // Ineqaulity constraints
        cr = 0;
        // Energy sum in min/max boundaries
        // -> n_opti_vars Equations: minE <= E cumsum <= maxE
        b1_pos_hp_minE_maxE_start = cr;
        for (unsigned int tm = 0; tm < T; tm++) {
            for (unsigned int t = 0; t <= tm; t++) {
                const unsigned long cri = cr * n_opti_vars; // current row ID
                A[cri + v_pos_p_hp_start + t] = Global::get_time_step_size_in_h();
            }
            cr += 1;
        }
        b1_pos_cs_minE_maxE_start = cr;
        for (unsigned int tm = 0; tm < T; tm++) {
            for (unsigned int t = 0; t <= tm; t++) {
                const unsigned long cri = cr * n_opti_vars; // current row ID
                A[cri + v_pos_p_cs_start + t] = Global::get_time_step_size_in_h();
            }
            cr += 1;
        }
    }

    /**
     * Deletes all global variables.
     */
    /*
    void vacuum() {
        if (objective != NULL) {
            delete objective;
            objective = NULL;
        }
    }*/

}

#endif
