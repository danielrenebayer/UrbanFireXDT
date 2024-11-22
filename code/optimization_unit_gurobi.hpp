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

namespace CUOptimization {

    /**
     * Initializes and executes the optimization for one control unit using gurobi optimizer (license required).
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
    bool runOptimizationForOneCU_Gurobi(
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
        try {
            GRBEnv env = GRBEnv(true);
            env.set(GRB_IntParam_OutputFlag, 0); // disable output
            env.start();
            GRBModel model = GRBModel(env);
            // create the variables
            const unsigned int T  = Global::get_control_horizon_in_ts(); // number of time steps to consider
            const unsigned int Tp = T + 1; // required for the e_bs
            std::vector<GRBVar> p_hp_kW(T);
            std::vector<GRBVar> e_hp_cum_kWh(T);
            std::vector<GRBVar> p_cs_kW(T);
            std::vector<GRBVar> e_cs_cum_kWh(T);
            std::vector<GRBVar> p_bs_in_kW(T);
            std::vector<GRBVar> p_bs_out_kW(T);
            std::vector<GRBVar> e_bs_kWh(T);
            std::vector<GRBVar> x_demand_kW(T); // power of grid demand in kW
            std::vector<GRBVar> x_feedin_kW(T); // power of grid feedin in kW
            for (unsigned int t = 0; t < T; t++) {
                const std::string tstr = to_string(t);
                p_hp_kW[t]      = model.addVar(0.0, max_p_hp_kW,  0.0, GRB_CONTINUOUS, "p_hp_kW" + tstr);
                e_hp_cum_kWh[t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "e_hp_cum_kWh" + tstr);
                p_cs_kW[t]      = model.addVar(0.0, max_p_cs_kW,  0.0, GRB_CONTINUOUS, "p_cs_kW" + tstr);
                e_cs_cum_kWh[t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "e_cs_cum_kWh" + tstr);
                p_bs_in_kW[t]   = model.addVar(0.0, max_p_bs_kW,  0.0, GRB_CONTINUOUS, "p_bs_in_kW" + tstr);
                p_bs_out_kW[t]  = model.addVar(0.0, max_p_bs_kW,  0.0, GRB_CONTINUOUS, "p_bs_out_kW" + tstr);
                x_demand_kW[t]  = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x_demand_kW" + tstr);
                x_feedin_kW[t]  = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x_feedin_kW" + tstr);
            }
            for (unsigned int t = 0; t < Tp; t++) {
                e_bs_kWh[t]     = model.addVar(0.0, max_e_bs_kWh, 0.0, GRB_CONTINUOUS, "e_bs_kWh"); // energy of the battery storage at the beginning of a time step
            }
            // init current battery state / initial SOC
            model.addConstr(e_bs_kWh[0] == current_bs_charge_kWh);
            // Battery balance equation
            for (unsigned int t = 0; t < T; t++) {
                model.addConstr(
                    e_bs_kWh[t + 1] == e_bs_kWh[t] * (1.0 - Global::get_exp_bess_self_ds_ts()) +
                                    p_bs_in_kW[t]  * Global::get_time_step_size_in_h() * Global::get_exp_bess_effi_in() -
                                    p_bs_out_kW[t] * Global::get_time_step_size_in_h() / Global::get_exp_bess_effi_out()
                );
            }
            // Power balance equation on control unit level
            if (Global::get_controller_allow_bs_grid_charging()) {
                // Case A: Battery charging from grid allowed
                for (unsigned int t = 0; t < T; t++) {
                    double resid_minus_pv_kW = future_resid_demand_kW[t] + future_hp_unshiftable_kW[t] - future_pv_generation_kW[t]; // i.e. local balance at time step t
                    model.addConstr(resid_minus_pv_kW
                                    + p_bs_in_kW[t]  - p_bs_out_kW[t]
                                    + p_hp_kW[t]     + p_cs_kW[t]
                                    + x_feedin_kW[t] - x_demand_kW[t] == 0);
                }
            } else {
                // Case B: No charging of the battery from grid allowed!
                /*
                // VARIANT WITH BINARY VARIABLES
                //
                std::vector<GRBVar> is_bs_charging(T);
                for (unsigned int t = 0; t < T; t++)
                    is_bs_charging[t] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "is_bs_charging");
                // compute Big M value
                float BigM = max_p_bs_kW + max_p_cs_kW + max_p_hp_kW + *std::max_element(future_resid_demand_kW.begin(), future_resid_demand_kW.end());
                BigM *= 2;
                if (BigM < 0)
                    BigM = 100;
                // add constraints linking the binary and the other variables
                for (unsigned int t = 0; t < T; t++) {
                    model.addConstr(p_bs_in_kW[t]  <= BigM *      is_bs_charging[t] );
                    model.addConstr(p_bs_out_kW[t] <= BigM * (1 - is_bs_charging[t]));
                    model.addConstr(x_feedin_kW[t] <= BigM *      is_bs_charging[t]) ;
                    model.addConstr(x_demand_kW[t] <= BigM * (1 - is_bs_charging[t]));
                }
                */
                //
                // VARIANT WITH TWO BALANCE EQUATIONS
                std::vector<GRBVar> p_hp_kW_1(T);
                std::vector<GRBVar> p_hp_kW_2(T);
                std::vector<GRBVar> p_cs_kW_1(T);
                std::vector<GRBVar> p_cs_kW_2(T);
                for (unsigned int t = 0; t < T; t++) {
                    p_hp_kW_1[t] = model.addVar(0.0, max_p_hp_kW,  0.0, GRB_CONTINUOUS, "p_hp_kW_1");
                    p_hp_kW_2[t] = model.addVar(0.0, max_p_hp_kW,  0.0, GRB_CONTINUOUS, "p_hp_kW_2");
                    p_cs_kW_1[t] = model.addVar(0.0, max_p_cs_kW,  0.0, GRB_CONTINUOUS, "p_cs_kW_1");
                    p_cs_kW_2[t] = model.addVar(0.0, max_p_cs_kW,  0.0, GRB_CONTINUOUS, "p_cs_kW_2");
                }
                for (unsigned int t = 0; t < T; t++) {
                    double resid_minus_pv_kW = future_resid_demand_kW[t] + future_hp_unshiftable_kW[t] - future_pv_generation_kW[t]; // i.e. local balance at time step t
                    double local_plus = 0, local_minus = 0;
                    if (resid_minus_pv_kW > 0) {
                        local_plus  =  resid_minus_pv_kW;
                    } else {
                        local_minus = -resid_minus_pv_kW;
                    }
                    // balance euqation 1
                    model.addConstr(local_plus + p_hp_kW_1[t] + p_cs_kW_1[t]
                                    - p_bs_out_kW[t] - x_demand_kW[t] == 0);
                    // balance equation 2
                    model.addConstr(p_hp_kW_2[t] + p_cs_kW_2[t] + x_feedin_kW[t] + p_bs_in_kW[t]
                                    - local_minus == 0);
                    // linkage of the power parts of HP and CS with actual power
                    model.addConstr(p_hp_kW[t] == p_hp_kW_1[t] + p_hp_kW_2[t]);
                    model.addConstr(p_cs_kW[t] == p_cs_kW_1[t] + p_cs_kW_2[t]);
                }
            }
            // Energy balance equation for heat pump and EV charging station
            for (unsigned int tm = 0; tm < T; tm++) {
                // t is the summation variable
                GRBLinExpr expr_hp_summation = 0.0;
                GRBLinExpr expr_cs_summation = 0.0;
                for (unsigned int t = 0; t <= tm; t++) {
                    // heat pump: summation variable
                    expr_hp_summation += p_hp_kW[t] * Global::get_time_step_size_in_h();
                    // heat pump balance
                    model.addConstr(e_hp_cum_kWh[tm] <= future_hp_shiftable_maxE[tm]);
                    model.addConstr(e_hp_cum_kWh[tm] >= future_hp_shiftable_minE[tm]);
                    // ev charging station: summation variable
                    expr_cs_summation += p_cs_kW[t] * Global::get_time_step_size_in_h();
                    // ev charging station balance
                    model.addConstr(e_cs_cum_kWh[tm] <= future_cs_shiftable_maxE[tm]);
                    model.addConstr(e_cs_cum_kWh[tm] >= future_cs_shiftable_minE[tm]);
                }
                model.addConstr(e_hp_cum_kWh[tm] == expr_hp_summation);
                model.addConstr(e_cs_cum_kWh[tm] == expr_cs_summation);
            }
            // Define the objective function
            GRBLinExpr obj = 0.0;
            if (Global::get_controller_optimization_target() == global::ControllerOptimizationTarget::ElectricityCosts) {
                for (unsigned int t = 0; t < T; t++) {
                    float current_demand_tariff = Global::get_demand_tariff();
                    if (global::eprices_local_ts != NULL && ts + t < Global::get_n_timesteps())
                        current_demand_tariff = global::eprices_local_ts[ts - 1 + t];
                    obj += x_demand_kW[t] * current_demand_tariff - x_feedin_kW[t] * Global::get_feed_in_tariff();
                }
            } else if (Global::get_controller_optimization_target() == global::ControllerOptimizationTarget::PeakLoad) {
                GRBVar max_demand_kW = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "max_demand_kW");
                obj = max_demand_kW;
                for (unsigned int t = 0; t < T; t++) {
                    model.addConstr(x_demand_kW[t] <= max_demand_kW);
                }
            } else /* if (Global::get_controller_optimization_target() == global::ControllerOptimizationTarget::Emissions) */ {
                for (unsigned int t = 0; t < T; t++) {
                    float current_emissions = Global::get_emissions_g_CO2eq_per_kWh();
                    if (global::emission_ts != NULL && ts + t < Global::get_n_timesteps())
                        current_emissions = global::emission_ts[ts - 1 + t];
                    obj += x_demand_kW[t] * current_emissions;
                    // TODO: What about feedin emissions? Ignore? Substract?
                }
            }
            model.setObjective(obj, GRB_MINIMIZE);
            //
            // Execute the optimization results and check results
            model.optimize();
            int model_status = model.get(GRB_IntAttr_Status);
            if (model_status != GRB_OPTIMAL) {
                std::cerr << "Optimization not resulting in optimal value.\n";
                std::cerr << "Gurobi model status = " << model_status << std::endl;
                return false;
            }
            //
            // Get the results and store them in the output reference
            for (unsigned int t = 0; t < T; t++) {
                outputRef.hp_power[t] = p_hp_kW[t].get(GRB_DoubleAttr_X);
                outputRef.cs_power[t] = p_cs_kW[t].get(GRB_DoubleAttr_X);
                outputRef.bs_power[t] = p_bs_in_kW[t].get(GRB_DoubleAttr_X) - p_bs_out_kW[t].get(GRB_DoubleAttr_X);
            }

        } catch (GRBException e) {
            std::cerr << "Error during optimization (code = " << e.getErrorCode() << ") with message:" << std::endl;
            std::cerr << e.getMessage() << std::endl;
            return false;
        } catch (...) {
            std::cerr << "Unknown exception during the optimization!" << std::endl;
            return false;
        }
        return true;
    }

}

#endif
