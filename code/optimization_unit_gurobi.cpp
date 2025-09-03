#include "optimization_unit_gurobi.hpp"

#include <optional>

#include "global.h"
#include "optimization_unit_general.hpp"

#include "gurobi_c++.h"

GRBEnv* GurobiLPController::env = NULL;


bool GurobiLPController::updateController(
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
)
{
    try {
        GRBModel model = GRBModel(*env);
        // only use the dual simplex (method=1), as barrier based methods lead to memory problems in a paralized setting
        model.set(GRB_IntParam_Method, 1);
        // create the variables
        const unsigned int T  = Global::get_control_horizon_in_ts(); // number of time steps to consider
        const unsigned int Tp = T + 1; // required for the e_bs
        const unsigned long n_pv_sections = optimize_PV_size ? total_PV_generation_per_section_kW->size() : 0;
        std::vector<GRBVar> p_resid_eq1(T);
        std::vector<GRBVar> p_pv_to_resid(T);
        std::vector<GRBVar> p_pv_eq2(T);
        std::vector<GRBVar> p_hp_kW(T);
        std::vector<GRBVar> e_hp_cum_kWh(T);
        std::vector<std::vector<GRBVar>> p_ev_kW(n_cars);      // first index -> EV, second index -> time step
        std::vector<std::vector<GRBVar>> e_ev_cum_kWh(n_cars);
        std::vector<GRBVar> p_bs_in_kW(T);
        std::vector<GRBVar> p_bs_out_kW(T);
        std::vector<GRBVar> e_bs_kWh(Tp);
        std::vector<GRBVar> x_demand_kW(T); // power of grid demand in kW
        std::vector<GRBVar> x_feedin_kW(T); // power of grid feedin in kW
        std::optional<GRBVar> dim_bs_E_kWh; // the optimal BS size
        std::optional<GRBVar> dim_bs_P_kW; // the (bound) actual battery power dependent on the (free) BS dimension
        std::vector<GRBVar> dim_pv_per_sec_kW(n_pv_sections); // optimal PV size per section
        // pre-initialization in case of PV or BS optimal sizing
        if (optimize_PV_size) {
            auto it  = max_PV_power_per_section_kWp->cbegin();
            auto it2 = total_PV_generation_per_section_kW->cbegin();
            for (unsigned long sectionIdx = 0; sectionIdx < n_pv_sections; sectionIdx++) {
                // safety checks
                if (it == max_PV_power_per_section_kWp->cend()) {
                    throw std::runtime_error("Error in the optimization unit: not enough values provided in argument 'max_PV_power_per_section_kWp' for method updateController().");
                }
                // another check for the content of the other list
                if ( (*it2).size() < T ) {
                    throw std::runtime_error("Error in the optimization unit: parameter 'max_PV_power_per_section_kWp' has a vector to small size for section with index " + to_string(sectionIdx) + ".");
                }
                // set upper limit for the PV size for this section
                double upper_limit_for_this_section_kWp = *it;
                dim_pv_per_sec_kW[sectionIdx] = model.addVar(0.0, upper_limit_for_this_section_kWp, 0.0, GRB_CONTINUOUS, "dim_pv_per_sec_kW_" + std::to_string(sectionIdx));
                // increment the iterators
                it++; it2++;
            }
        }
        if (optimize_BS_size) {
            max_e_bs_kWh = GRB_INFINITY; // set maximum BS size to infinity
            max_p_bs_kW  = GRB_INFINITY; // set maximum BS power to infinity
            dim_bs_E_kWh = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "dim_bs_E_kWh");
            dim_bs_P_kW  = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "dim_bs_P_kW");
            // bound the battery power to the free battery capacity
            model.addConstr( *dim_bs_P_kW * Global::get_exp_bess_E_P_ratio() == *dim_bs_E_kWh );
            // a final check
            if (current_bs_charge_kWh > 0.0) {
                std::cerr << "Warning: current_bs_charge_kWh > 0 --> Setting this value to 0.0" << std::endl;
                current_bs_charge_kWh = 0.0;
            }
        }
        // standard initialization
        for (unsigned int t = 0; t < T; t++) {
            const std::string tstr = to_string(t);
            p_resid_eq1[t]  = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "p_resid_eq1_"   + tstr);
            p_pv_to_resid[t]= model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "p_pv_to_resid_" + tstr);
            p_pv_eq2[t]     = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "p_pv_eq2_"      + tstr);
            p_hp_kW[t]      = model.addVar(future_hp_shiftable_minP[t], future_hp_shiftable_maxP[t], 0.0, GRB_CONTINUOUS, "p_hp_kW" + tstr);
            e_hp_cum_kWh[t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "e_hp_cum_kWh" + tstr);
            p_bs_in_kW[t]   = model.addVar(0.0, max_p_bs_kW,  0.0, GRB_CONTINUOUS, "p_bs_in_kW" + tstr);
            p_bs_out_kW[t]  = model.addVar(0.0, max_p_bs_kW,  0.0, GRB_CONTINUOUS, "p_bs_out_kW" + tstr);
            x_demand_kW[t]  = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x_demand_kW" + tstr);
            x_feedin_kW[t]  = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x_feedin_kW" + tstr);
            // add bounds for BS power if this BS size is an optimization variable
            if (optimize_BS_size) {
                model.addConstr( p_bs_in_kW[t] + p_bs_out_kW[t] <= *dim_bs_P_kW );
            }
        }
        // one additional time step of the battery storage
        for (unsigned int t = 0; t < Tp; t++) {
            const std::string tstr = to_string(t);
            e_bs_kWh[t]     = model.addVar(0.0, max_e_bs_kWh, 0.0, GRB_CONTINUOUS, "e_bs_kWh_" + tstr); // energy of the battery storage at the beginning of a time step
            if (optimize_BS_size) {
                model.addConstr( e_bs_kWh[t] <= *dim_bs_E_kWh );
            }
        }
        // separate initialization for the EVs
        for (unsigned long evIdx = 0; evIdx < n_cars; evIdx++) {
            p_ev_kW[evIdx].clear();
            e_ev_cum_kWh[evIdx].clear();
            for (unsigned int t = 0; t < T; t++) {
                const std::string tstr = to_string(t);
                const std::string estr = to_string(evIdx);
                //
                p_ev_kW[evIdx].push_back(
                    model.addVar(0.0, future_ev_maxP->at(evIdx)->at(t),  0.0, GRB_CONTINUOUS, "p_ev" + estr + "_kW" + tstr)
                );
                e_ev_cum_kWh[evIdx].push_back(
                    model.addVar(
                        future_ev_shiftable_minE->at(evIdx)->at(t), /* min */
                        future_ev_shiftable_maxE->at(evIdx)->at(t), /* max */
                        0.0, GRB_CONTINUOUS, "e_ev" + estr + "_cum_kWh" + tstr)
                );
                
            }
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
        // Maximum charging station power limits
        for (unsigned int t = 0; t < T; t++) {
            GRBLinExpr expr_cs_p_at_t = 0.0;
            for (unsigned long evIdx = 0; evIdx < n_cars; evIdx++) {
                expr_cs_p_at_t += p_ev_kW[evIdx][t];
            }
            model.addConstr( expr_cs_p_at_t <= max_p_cs_kW );
        }
        // Power balance for PV and residential demand
        for (unsigned int t = 0; t < T; t++) {
            model.addConstr(future_resid_demand_kW[t]  == p_resid_eq1[t] + p_pv_to_resid[t]);
          if (!optimize_PV_size) {
            model.addConstr(future_pv_generation_kW[t] == p_pv_eq2[t]    + p_pv_to_resid[t]);
          } else {
                GRBLinExpr expr_pv_gen_sum = 0.0;
                auto it = total_PV_generation_per_section_kW->cbegin();
                for (unsigned long sectionIdx = 0; sectionIdx < n_pv_sections; sectionIdx++) {
                    // safety checks not required, as n_pv_sections equals total_PV_generation_per_section_kW->size()
                    const std::vector<double>& time_series_for_this_section = *it++;
                    expr_pv_gen_sum += dim_pv_per_sec_kW[sectionIdx] * time_series_for_this_section[t];
                }
                model.addConstr( p_pv_eq2[t] + p_pv_to_resid[t] == expr_pv_gen_sum );
          }
        }
        // Power balance equation on control unit level
        if (Global::get_controller_bs_grid_charging_mode() == global::ControllerBSGridChargingMode::GridChargingAndDischarging) {
            // Case A: Battery charging from grid allowed
            for (unsigned int t = 0; t < T; t++) {
                GRBLinExpr expr_cs_p_at_t = 0.0;
                for (unsigned long evIdx = 0; evIdx < n_cars; evIdx++) {
                    expr_cs_p_at_t += p_ev_kW[evIdx][t];
                }
                model.addConstr(p_resid_eq1[t]   - p_pv_eq2[t]
                                + p_bs_in_kW[t]  - p_bs_out_kW[t]
                                + p_hp_kW[t]     + expr_cs_p_at_t
                                + x_feedin_kW[t] - x_demand_kW[t] == 0);
            }
        } else if (Global::get_controller_bs_grid_charging_mode() == global::ControllerBSGridChargingMode::OnlyGridCharging) {
            // Case B: Charging from the grid allowed, but no discharging into the grid possible
            std::vector<GRBVar> p_bs_in_kW_1(T);
            std::vector<GRBVar> p_bs_in_kW_2(T);
            std::vector<GRBVar> p_hp_kW_1(T);
            std::vector<GRBVar> p_hp_kW_2(T);
            std::vector<GRBVar> p_cs_kW_1(T);
            std::vector<GRBVar> p_cs_kW_2(T);
            for (unsigned int t = 0; t < T; t++) {
                const std::string tstr = to_string(t);
                p_bs_in_kW_1[t] = model.addVar(0.0, max_p_bs_kW,  0.0, GRB_CONTINUOUS, "p_bs_in_kW_BalEq1_"  + tstr);
                p_bs_in_kW_2[t] = model.addVar(0.0, max_p_bs_kW,  0.0, GRB_CONTINUOUS, "p_bs_in_kW_BalEq2_"  + tstr);
                p_hp_kW_1[t] = model.addVar(0.0, GRB_INFINITY,  0.0, GRB_CONTINUOUS, "p_hp_kW_BalEq1_" + tstr);
                p_hp_kW_2[t] = model.addVar(0.0, GRB_INFINITY,  0.0, GRB_CONTINUOUS, "p_hp_kW_BalEq2_" + tstr);
                p_cs_kW_1[t] = model.addVar(0.0, max_p_cs_kW,  0.0, GRB_CONTINUOUS, "p_cs_kW_BalEq1_"  + tstr);
                p_cs_kW_2[t] = model.addVar(0.0, max_p_cs_kW,  0.0, GRB_CONTINUOUS, "p_cs_kW_BalEq2_"  + tstr);
            }
            for (unsigned int t = 0; t < T; t++) {
                GRBLinExpr expr_cs_p_at_t = 0.0;
                for (unsigned long evIdx = 0; evIdx < n_cars; evIdx++) {
                    expr_cs_p_at_t += p_ev_kW[evIdx][t];
                }
                // balance euqation 1
                model.addConstr(p_resid_eq1[t] + p_hp_kW_1[t] + p_cs_kW_1[t]
                                + p_bs_in_kW_1[t]
                                - p_bs_out_kW[t] - x_demand_kW[t] == 0);
                // balance equation 2
                model.addConstr(p_hp_kW_2[t] + p_cs_kW_2[t] + x_feedin_kW[t] + p_bs_in_kW_2[t]
                                - p_pv_eq2[t] == 0);
                // linkage of the power parts of BS, HP and CS with actual power
                model.addConstr(p_bs_in_kW[t]  == p_bs_in_kW_1[t] + p_bs_in_kW_2[t]);
                model.addConstr(p_hp_kW[t]     == p_hp_kW_1[t]    + p_hp_kW_2[t]);
                model.addConstr(expr_cs_p_at_t == p_cs_kW_1[t]    + p_cs_kW_2[t]);
            }
        } else {
            // Case C: No charging of the battery from grid allowed!
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
                const std::string tstr = to_string(t);
                p_hp_kW_1[t] = model.addVar(0.0, GRB_INFINITY,  0.0, GRB_CONTINUOUS, "p_hp_kW_BalEq1_" + tstr);
                p_hp_kW_2[t] = model.addVar(0.0, GRB_INFINITY,  0.0, GRB_CONTINUOUS, "p_hp_kW_BalEq2_" + tstr);
                p_cs_kW_1[t] = model.addVar(0.0, max_p_cs_kW,  0.0, GRB_CONTINUOUS, "p_cs_kW_BalEq1_"  + tstr);
                p_cs_kW_2[t] = model.addVar(0.0, max_p_cs_kW,  0.0, GRB_CONTINUOUS, "p_cs_kW_BalEq2_"  + tstr);
            }
            for (unsigned int t = 0; t < T; t++) {
                GRBLinExpr expr_cs_p_at_t = 0.0;
                for (unsigned long evIdx = 0; evIdx < n_cars; evIdx++) {
                    expr_cs_p_at_t += p_ev_kW[evIdx][t];
                }
                // balance euqation 1
                model.addConstr(p_resid_eq1[t] + p_hp_kW_1[t] + p_cs_kW_1[t]
                                - p_bs_out_kW[t] - x_demand_kW[t] == 0);
                // balance equation 2
                model.addConstr(p_hp_kW_2[t] + p_cs_kW_2[t] + x_feedin_kW[t] + p_bs_in_kW[t]
                                - p_pv_eq2[t] == 0);
                // linkage of the power parts of HP and CS with actual power
                model.addConstr(p_hp_kW[t] == p_hp_kW_1[t] + p_hp_kW_2[t]);
                model.addConstr(expr_cs_p_at_t == p_cs_kW_1[t] + p_cs_kW_2[t]);
            }
        }
        // Energy balance equation for heat pump and EV charging station
        for (unsigned int tm = 0; tm < T; tm++) {
            // t is the summation variable
            GRBLinExpr expr_hp_summation = 0.0;
            for (unsigned int t = 0; t <= tm; t++) {
                // heat pump: summation variable
                expr_hp_summation += p_hp_kW[t] * Global::get_time_step_size_in_h();
            }
            model.addConstr(e_hp_cum_kWh[tm] == expr_hp_summation);
            // heat pump balance
            model.addConstr(e_hp_cum_kWh[tm] <= future_hp_shiftable_maxE[tm]);
            model.addConstr(e_hp_cum_kWh[tm] >= future_hp_shiftable_minE[tm]);
        }
        for (unsigned long evIdx = 0; evIdx < n_cars; evIdx++) {
            for (unsigned int tm = 0; tm < T; tm++) {
                GRBLinExpr expr_ev_E_cumsum = 0.0;
                for (unsigned int t = 0; t <= tm; t++) {
                    // Upper and lower limits for 'EV cumsum E' variable already set by definition
                    expr_ev_E_cumsum += p_ev_kW[evIdx][t] * Global::get_time_step_size_in_h();
                }
                model.addConstr(e_ev_cum_kWh[evIdx][tm] == expr_ev_E_cumsum);
            }
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
            if (optimize_PV_size) {
                for (unsigned long sectionIdx = 0; sectionIdx < n_pv_sections; sectionIdx++) {
                    obj += dim_pv_per_sec_kW[sectionIdx] * Global::get_annuity_factor() * Global::get_inst_cost_PV_per_kWp();
                }
            }
            if (optimize_BS_size) {
                obj += *dim_bs_E_kWh * Global::get_annuity_factor() * Global::get_inst_cost_BS_per_kWh();
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
            string filepath = "/tmp/gurobi_model_infeasible_";
            filepath += std::to_string(controlUnitID) + "_" + std::to_string(ts);
            filepath += ".lp";
            /*
            std::cerr << "    Computing irreducible inconsistent subsystem ...";
            model.write(filepath);
            filepath = "/tmp/gurobi_model_infeasible_";
            filepath += std::to_string(controlUnitID) + "_" + std::to_string(ts);
            filepath += "_iis.ilp";
            model.computeIIS();
            model.write(filepath);
            std::cerr << "    IIS written to " << filepath << std::endl;
            */
            return false;
        }
        //
        // Get the results and store them in the output reference
        for (unsigned int t = 0; t < T; t++) {
            bs_power[t] = p_bs_in_kW[t].get(GRB_DoubleAttr_X) - p_bs_out_kW[t].get(GRB_DoubleAttr_X);
            hp_power[t] = p_hp_kW[t].get(GRB_DoubleAttr_X);
            for (unsigned long evIdx = 0; evIdx < n_cars; evIdx++) {
                ev_power[evIdx][t] = p_ev_kW[evIdx][t].get(GRB_DoubleAttr_X);
            }
            if (optimize_PV_size) {
                optimal_pv_size_per_section_kW.clear();
                for (unsigned long sectionIdx = 0; sectionIdx < n_pv_sections; sectionIdx++) {
                    optimal_pv_size_per_section_kW.push_back( dim_pv_per_sec_kW[sectionIdx].get(GRB_DoubleAttr_X) );
                }
            }
            if (optimize_BS_size) {
                optimal_bs_size_kWh = dim_bs_E_kWh.value().get(GRB_DoubleAttr_X);
            }
        }
        model.terminate();

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

void GurobiLPController::InitializeGurobiEnvironment() {
    env = new GRBEnv(true);
    env->set(GRB_IntParam_OutputFlag, 0); // disable output
    env->start();
}

void GurobiLPController::VaccumAllStaticVariables() {
    if (env != NULL)
        delete env;
    env = NULL;
}

