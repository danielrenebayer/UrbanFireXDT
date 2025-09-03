/**
 * optimization_unit_or_tools.hpp
 *
 * This file contains all functions required to call gurobi as
 * optimizer inside the control units.
 */

#ifndef OPTIMIZATION_UNIT_OR_TOOLS_HPP
#define OPTIMIZATION_UNIT_OR_TOOLS_HPP

#include "global.h"
#include "optimization_unit_general.hpp"

#include "ortools/linear_solver/linear_solver.h"

using namespace operations_research;

class ORToolsLPController : public BaseOptimizedController {

    public:
    ORToolsLPController(unsigned long cuID, unsigned int time_horizon, unsigned long n_cars) :
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
    ) {
        // Initialize the solver
        MPSolver* model = MPSolver::CreateSolver("SCIP");
        //MPSolver* model = MPSolver::CreateSolver("GLOP");
        if (!model) {
                std::cout << "SCIP solver unavailable." << std::endl;
                //std::cout << "GLOP solver unavailable." << std::endl;
                return false;
        }
        const double infinity = model->infinity();
        //
        // Create the variables
        const unsigned int T  = Global::get_control_horizon_in_ts(); // number of time steps to consider
        const unsigned int Tp = T + 1; // required for the e_bs
        const unsigned long n_pv_sections = optimize_PV_size ? total_PV_generation_per_section_kW->size() : 0;
        std::vector<MPVariable*> p_resid_eq1(T);
        std::vector<MPVariable*> p_pv_to_resid(T);
        std::vector<MPVariable*> p_pv_eq2(T);
        std::vector<MPVariable*> p_hp_kW(T);
        std::vector<MPVariable*> e_hp_cum_kWh(T);
        std::vector<std::vector<MPVariable*>> p_ev_kW(n_cars);      // first index -> EV, second index -> time step
        std::vector<std::vector<MPVariable*>> e_ev_cum_kWh(n_cars);
        std::vector<MPVariable*> p_bs_in_kW(T);
        std::vector<MPVariable*> p_bs_out_kW(T);
        std::vector<MPVariable*> e_bs_kWh(Tp);
        std::vector<MPVariable*> x_demand_kW(T); // power of grid demand in kW
        std::vector<MPVariable*> x_feedin_kW(T); // power of grid feedin in kW
        MPVariable* dim_bs_E_kWh; // the optimal BS size
        MPVariable* dim_bs_P_kW; // the (bound) actual battery power dependent on the (free) BS dimension
        std::vector<MPVariable*> dim_pv_per_sec_kW(n_pv_sections); // optimal PV size per section
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
                dim_pv_per_sec_kW = model->MakeNumVar(0.0, upper_limit_for_this_section_kWp, "dim_pv_per_sec_kW_" + std::to_string(sectionIdx));
                // increment the iterators
                it++; it2++;
            }
        }
        if (optimize_BS_size) {
            max_e_bs_kWh = infinity; // set maximum BS size to infinity
            max_p_bs_kW  = infinity; // set maximum BS power to infinity
            dim_bs_E_kWh = model->MakeNumVar(0.0, infinity, "dim_bs_E_kWh");
            dim_bs_P_kW  = model->MakeNumVar(0.0, infinity, "dim_bs_P_kW");
            // bound the battery power to the free battery capacity
            MPConstraint* const c_bs_P = model->MakeRowConstraint(0.0, 0.0, "Bound BS power to capacity.");
            c_bs_P->SetCoefficient(dim_bs_P_kW,  -1.0);
            c_bs_P->SetCoefficient(dim_bs_E_kWh, 1.0 / Global::get_exp_bess_E_P_ratio());
            // a final check
            if (current_bs_charge_kWh > 0.0) {
                std::cerr << "Warning: current_bs_charge_kWh > 0 --> Setting this value to 0.0" << std::endl;
                current_bs_charge_kWh = 0.0;
            }
        }
        // standard initialization
        for (unsigned int t = 0; t < T; t++) {
            const std::string tstr = to_string(t);
            p_resid_eq1[t]  = model->MakeNumVar(0.0, infinity, "p_resid_eq1_"     + tstr);
            p_pv_to_resid[t]= model->MakeNumVar(0.0, infinity, "p_pv_to_resid_"   + tstr);
            p_pv_eq2[t]     = model->MakeNumVar(0.0, infinity, "p_pv_eq2_"        + tstr);
            p_hp_kW[t]      = model->MakeNumVar(future_hp_shiftable_minP[t], future_hp_shiftable_maxP[t],  "p_hp_kW_"     + tstr);
            // NOTICE: The next line is the only difference to the gurobi implementation: Bounds for cumsum E are set here directly
            e_hp_cum_kWh[t] = model->MakeNumVar(future_hp_shiftable_minE[t], future_hp_shiftable_maxE[t], "e_hp_cum_kWh_"+ tstr);
            p_bs_in_kW[t]   = model->MakeNumVar(0.0, max_p_bs_kW,  "p_bs_in_kW_"  + tstr);
            p_bs_out_kW[t]  = model->MakeNumVar(0.0, max_p_bs_kW,  "p_bs_out_kW_" + tstr);
            x_demand_kW[t]  = model->MakeNumVar(0.0, infinity,     "x_demand_kW_" + tstr);
            x_feedin_kW[t]  = model->MakeNumVar(0.0, infinity,     "x_feedin_kW_" + tstr);
            // add bounds for BS power if this BS size is an optimization variable
            if (optimize_BS_size) {
                MPConstraint* const c_bs_P_upper_limit = model->MakeRowConstraint(-infinity, 0.0, "Bound for BS power " + tstr);
                c_bs_P_upper_limit->SetCoefficient(p_bs_in_kW[t],  1.0);
                c_bs_P_upper_limit->SetCoefficient(p_bs_out_kW[t], 1.0);
                c_bs_P_upper_limit->SetCoefficient(dim_bs_P_kW,   -1.0);
            }
        }
        // one additional time step of the battery storage
        for (unsigned int t = 0; t < Tp; t++) {
            const std::string tstr = to_string(t);
            e_bs_kWh[t]     = model->MakeNumVar(0.0, max_e_bs_kWh, "e_bs_kWh_" + tstr); // energy of the battery storage at the beginning of a time step
            if (optimize_BS_size) {
                MPConstraint* const c_bs_E_upper_limit = model->MakeRowConstraint(-infinity, 0.0, "Bound for BS capacity " + tstr);
                c_bs_E_upper_limit->SetCoefficient(e_bs_kWh[t],   1.0);
                c_bs_E_upper_limit->SetCoefficient(dim_bs_E_kWh, -1.0);
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
                    model->MakeNumVar(0.0, future_ev_maxP->at(evIdx)->at(t), "p_ev" + estr + "_kW_" + tstr)
                );
                e_ev_cum_kWh[evIdx].push_back(
                    model->MakeNumVar(
                        future_ev_shiftable_minE->at(evIdx)->at(t), /* min */
                        future_ev_shiftable_maxE->at(evIdx)->at(t), /* max */
                        "e_ev" + estr + "_cum_kWh_" + tstr)
                );
                
            }
        }
        // init current battery state / initial SOC
        MPConstraint* const c0 = model->MakeRowConstraint(current_bs_charge_kWh, current_bs_charge_kWh, "Init BS state");
        c0->SetCoefficient(e_bs_kWh[0], 1.0);
        // Battery balance equation
        for (unsigned int t = 0; t < T; t++) {
            /*
            e_bs_kWh[t + 1] == e_bs_kWh[t] * (1.0 - Global::get_exp_bess_self_ds_ts()) +
                                p_bs_in_kW[t]  * Global::get_time_step_size_in_h() * Global::get_exp_bess_effi_in() -
                                p_bs_out_kW[t] * Global::get_time_step_size_in_h() / Global::get_exp_bess_effi_out()
            */
            const std::string tstr = to_string(t);
            MPConstraint* const c = model->MakeRowConstraint(0.0, 0.0, "BESS balance " + tstr);
            c->SetCoefficient(e_bs_kWh[t + 1], -1.0);
            c->SetCoefficient(e_bs_kWh[t],      1.0 - Global::get_exp_bess_self_ds_ts());
            c->SetCoefficient(p_bs_in_kW[t],    Global::get_time_step_size_in_h() * Global::get_exp_bess_effi_in());
            c->SetCoefficient(p_bs_out_kW[t],  -Global::get_time_step_size_in_h() / Global::get_exp_bess_effi_out());
            // Alternative:
            // auto c = model->MakeRowConstraint(-INF, INF, "BESS balance " + tstr);
            // c->SetOffset(0.0);
        }
        // Maximum charging station power limits
        for (unsigned int t = 0; t < T; t++) {
            const std::string tstr = to_string(t);
            MPConstraint* const c = model->MakeRowConstraint(-infinity, max_p_cs_kW, "Max CS power " + tstr);
            for (unsigned long evIdx = 0; evIdx < n_cars; evIdx++) {
                c->SetCoefficient( p_ev_kW[evIdx][t], 1.0 );
            }
        }
        // Power balance for PV and residential demand
        for (unsigned int t = 0; t < T; t++) {
            const std::string tstr = to_string(t);
            MPConstraint* const c1 = model->MakeRowConstraint(future_resid_demand_kW[t],  future_resid_demand_kW[t],  "Residential demand linkage " + tstr);
            c1->SetCoefficient(p_resid_eq1[t],    1.0);
            c1->SetCoefficient(p_pv_to_resid[t],  1.0);
          if (!optimize_PV_size) {
            MPConstraint* const c2 = model->MakeRowConstraint(future_pv_generation_kW[t], future_pv_generation_kW[t], "PV generation linkage " + tstr);
            c2->SetCoefficient(p_pv_eq2[t],       1.0);
            c2->SetCoefficient(p_pv_to_resid[t],  1.0);
          } else {
                MPConstraint* const c2_PVopti = model->MakeRowConstraint(0.0, 0.0, "PV generation linkage in case of PV size opti " + tstr);
                c2_PVopti->SetCoefficient(p_pv_eq2[t],       1.0);
                c2_PVopti->SetCoefficient(p_pv_to_resid[t],  1.0);
                auto it = total_PV_generation_per_section_kW->cbegin();
                for (unsigned long sectionIdx = 0; sectionIdx < n_pv_sections; sectionIdx++) {
                    // safety checks not required, as n_pv_sections equals total_PV_generation_per_section_kW->size()
                    const std::vector<double>& time_series_for_this_section = *it++;
                    c2_PVopti->SetCoefficient(dim_pv_per_sec_kW[sectionIdx], -time_series_for_this_section[t]);
                }
          }
        }
        // Power balance equation on control unit level
        if (Global::get_controller_bs_grid_charging_mode() == global::ControllerBSGridChargingMode::GridChargingAndDischarging) {
            // Case A: Battery charging from grid allowed
            for (unsigned int t = 0; t < T; t++) {
                const std::string tstr = to_string(t);
                double resid_minus_pv_kW = future_resid_demand_kW[t] - future_pv_generation_kW[t]; // i.e. local balance at time step t
                MPConstraint* const c = model->MakeRowConstraint(-resid_minus_pv_kW, -resid_minus_pv_kW, "CU Power Balance " + tstr);
                c->SetCoefficient(p_bs_in_kW[t],  1.0);
                c->SetCoefficient(p_bs_out_kW[t],-1.0);
                c->SetCoefficient(p_hp_kW[t],     1.0);
                c->SetCoefficient(x_feedin_kW[t], 1.0);
                c->SetCoefficient(x_demand_kW[t],-1.0);
                for (unsigned long evIdx = 0; evIdx < n_cars; evIdx++) {
                    c->SetCoefficient(p_ev_kW[evIdx][t], 1.0);
                }
            }
        } else if (Global::get_controller_bs_grid_charging_mode() == global::ControllerBSGridChargingMode::OnlyGridCharging) {
            // Case B: Charging from the grid allowed, but no discharging into the grid possible
            // VARIANT WITH TWO BALANCE EQUATIONS
            std::vector<MPVariable*> p_bs_in_kW_1(T);
            std::vector<MPVariable*> p_bs_in_kW_2(T);
            std::vector<MPVariable*> p_hp_kW_1(T);
            std::vector<MPVariable*> p_hp_kW_2(T);
            std::vector<MPVariable*> p_cs_kW_1(T);
            std::vector<MPVariable*> p_cs_kW_2(T);
            for (unsigned int t = 0; t < T; t++) {
                const std::string tstr = to_string(t);
                p_bs_in_kW_1[t] = model->MakeNumVar(0.0, max_p_bs_kW, "p_bs_kW_BalEq1_" + tstr);
                p_bs_in_kW_2[t] = model->MakeNumVar(0.0, max_p_bs_kW, "p_bs_kW_BalEq2_" + tstr);
                p_hp_kW_1[t] = model->MakeNumVar(0.0, infinity, "p_hp_kW_BalEq1_" + tstr);
                p_hp_kW_2[t] = model->MakeNumVar(0.0, infinity, "p_hp_kW_BalEq2_" + tstr);
                p_cs_kW_1[t] = model->MakeNumVar(0.0, max_p_cs_kW, "p_cs_kW_BalEq1_" + tstr);
                p_cs_kW_2[t] = model->MakeNumVar(0.0, max_p_cs_kW, "p_cs_kW_BalEq2_" + tstr);
            }
            for (unsigned int t = 0; t < T; t++) {
                const std::string tstr = to_string(t);
                // balance euqation 1
                MPConstraint* const c1 = model->MakeRowConstraint(0.0, 0.0, "CU Power Balance EQ1 " + tstr);
                c1->SetCoefficient(p_resid_eq1[t],  1.0);
                c1->SetCoefficient(p_hp_kW_1[t],    1.0);
                c1->SetCoefficient(p_cs_kW_1[t],    1.0);
                c1->SetCoefficient(p_bs_out_kW[t], -1.0);
                c1->SetCoefficient(x_demand_kW[t], -1.0);
                c1->SetCoefficient(p_bs_in_kW_1[t], 1.0);
                // balance equation 2
                MPConstraint* const c2 = model->MakeRowConstraint(0.0, 0.0, "CU Power Balance EQ2 " + tstr);
                c2->SetCoefficient(p_hp_kW_2[t],   1.0);
                c2->SetCoefficient(p_cs_kW_2[t],   1.0);
                c2->SetCoefficient(x_feedin_kW[t], 1.0);
                c2->SetCoefficient(p_bs_in_kW_2[t],1.0);
                c2->SetCoefficient(p_pv_eq2[t],   -1.0);
                // linkage of the power parts of BS, HP and CS with actual power
                MPConstraint* const cl0 = model->MakeRowConstraint(0.0, 0.0, "CU P Bal Linkage BS " + tstr);
                cl0->SetCoefficient(p_bs_in_kW[t],  -1.0);
                cl0->SetCoefficient(p_bs_in_kW_1[t], 1.0);
                cl0->SetCoefficient(p_bs_in_kW_2[t], 1.0);
                MPConstraint* const cl1 = model->MakeRowConstraint(0.0, 0.0, "CU P Bal Linkage HP " + tstr);
                cl1->SetCoefficient(p_hp_kW[t],  -1.0);
                cl1->SetCoefficient(p_hp_kW_1[t], 1.0);
                cl1->SetCoefficient(p_hp_kW_2[t], 1.0);
                MPConstraint* const cl2 = model->MakeRowConstraint(0.0, 0.0, "CU P Bal Linkage CS " + tstr);
                cl2->SetCoefficient(p_cs_kW_1[t], 1.0);
                cl2->SetCoefficient(p_cs_kW_2[t], 1.0);
                for (unsigned long evIdx = 0; evIdx < n_cars; evIdx++) {
                    cl2->SetCoefficient(p_ev_kW[evIdx][t], -1.0);
                }
            }
        } else {
            // Case C: No charging of the battery from grid allowed!
            // VARIANT WITH TWO BALANCE EQUATIONS
            std::vector<MPVariable*> p_hp_kW_1(T);
            std::vector<MPVariable*> p_hp_kW_2(T);
            std::vector<MPVariable*> p_cs_kW_1(T);
            std::vector<MPVariable*> p_cs_kW_2(T);
            for (unsigned int t = 0; t < T; t++) {
                const std::string tstr = to_string(t);
                p_hp_kW_1[t] = model->MakeNumVar(0.0, infinity, "p_hp_kW_BalEq1_" + tstr);
                p_hp_kW_2[t] = model->MakeNumVar(0.0, infinity, "p_hp_kW_BalEq2_" + tstr);
                p_cs_kW_1[t] = model->MakeNumVar(0.0, max_p_cs_kW, "p_cs_kW_BalEq1_" + tstr);
                p_cs_kW_2[t] = model->MakeNumVar(0.0, max_p_cs_kW, "p_cs_kW_BalEq2_" + tstr);
            }
            for (unsigned int t = 0; t < T; t++) {
                const std::string tstr = to_string(t);
                // balance euqation 1
                MPConstraint* const c1 = model->MakeRowConstraint(0.0, 0.0, "CU Power Balance EQ1 " + tstr);
                c1->SetCoefficient(p_resid_eq1[t],  1.0);
                c1->SetCoefficient(p_hp_kW_1[t],    1.0);
                c1->SetCoefficient(p_cs_kW_1[t],    1.0);
                c1->SetCoefficient(p_bs_out_kW[t], -1.0);
                c1->SetCoefficient(x_demand_kW[t], -1.0);
                // balance equation 2
                MPConstraint* const c2 = model->MakeRowConstraint(0.0, 0.0, "CU Power Balance EQ2 " + tstr);
                c2->SetCoefficient(p_hp_kW_2[t],   1.0);
                c2->SetCoefficient(p_cs_kW_2[t],   1.0);
                c2->SetCoefficient(x_feedin_kW[t], 1.0);
                c2->SetCoefficient(p_bs_in_kW[t],  1.0);
                c2->SetCoefficient(p_pv_eq2[t],   -1.0);
                // linkage of the power parts of HP and CS with actual power
                MPConstraint* const cl1 = model->MakeRowConstraint(0.0, 0.0, "CU P Bal Linkage 1 " + tstr);
                cl1->SetCoefficient(p_hp_kW[t],  -1.0);
                cl1->SetCoefficient(p_hp_kW_1[t], 1.0);
                cl1->SetCoefficient(p_hp_kW_2[t], 1.0);
                MPConstraint* const cl2 = model->MakeRowConstraint(0.0, 0.0, "CU P Bal Linkage 2 " + tstr);
                cl2->SetCoefficient(p_cs_kW_1[t], 1.0);
                cl2->SetCoefficient(p_cs_kW_2[t], 1.0);
                for (unsigned long evIdx = 0; evIdx < n_cars; evIdx++) {
                    cl2->SetCoefficient(p_ev_kW[evIdx][t], -1.0);
                }
            }
        }
        // Energy balance equation for heat pump and EV charging station
        for (unsigned int tm = 0; tm < T; tm++) {
            const std::string tmstr = to_string(tm);
            // NOTICE: Only difference to gurobi implementation:
            //         Upper and lower bounds for cumsum HP E are set in the definition already !
            MPConstraint* const c = model->MakeRowConstraint(0.0, 0.0, "HP Balance " + tmstr);
            // t is the summation variable for the previous time steps
            for (unsigned int t = 0; t <= tm; t++) {
                // heat pump: summation variable
                c->SetCoefficient(p_hp_kW[t], Global::get_time_step_size_in_h());
            }
            c->SetCoefficient(e_hp_cum_kWh[tm], -1.0);
        }
        for (unsigned long evIdx = 0; evIdx < n_cars; evIdx++) {
            for (unsigned int tm = 0; tm < T; tm++) {
                const std::string tmstr = to_string(tm);
                MPConstraint* const c = model->MakeRowConstraint(0.0, 0.0, "EV Balance " + tmstr);
                for (unsigned int t = 0; t <= tm; t++) {
                    // Upper and lower limits for 'EV cumsum E' variable already set by definition
                    c->SetCoefficient(p_ev_kW[evIdx][t], Global::get_time_step_size_in_h());
                }
                c->SetCoefficient(e_ev_cum_kWh[evIdx][tm], -1.0);
            }
        }
        // Define the objective function
        MPObjective* const objective = model->MutableObjective();
        if (Global::get_controller_optimization_target() == global::ControllerOptimizationTarget::ElectricityCosts) {
            for (unsigned int t = 0; t < T; t++) {
                float current_demand_tariff = Global::get_demand_tariff();
                if (global::eprices_local_ts != NULL && ts + t < Global::get_n_timesteps())
                    current_demand_tariff = global::eprices_local_ts[ts - 1 + t];
                objective->SetCoefficient(x_demand_kW[t],  current_demand_tariff);
                objective->SetCoefficient(x_feedin_kW[t], -Global::get_feed_in_tariff());
            }
            if (optimize_PV_size) {
                for (unsigned long sectionIdx = 0; sectionIdx < n_pv_sections; sectionIdx++) {
                    objective->SetCoefficient(dim_pv_per_sec_kW[sectionIdx], Global::get_annuity_factor() * Global::get_inst_cost_PV_per_kWp());
                }
            }
            if (optimize_BS_size) {
                objective->SetCoefficient(dim_bs_E_kWh, Global::get_annuity_factor() * Global::get_inst_cost_BS_per_kWh());
            }
        } else if (Global::get_controller_optimization_target() == global::ControllerOptimizationTarget::PeakLoad) {
            MPVariable *const max_demand_kW = model->MakeNumVar(0.0, infinity, "max_demand_kW");
            objective->SetCoefficient(max_demand_kW, 1.0);
            for (unsigned int t = 0; t < T; t++) {
                MPConstraint* const c = model->MakeRowConstraint(-infinity, 0.0);
                c->SetCoefficient(x_demand_kW[t], 1.0);
                c->SetCoefficient(max_demand_kW, -1.0);
            }
        } else /* if (Global::get_controller_optimization_target() == global::ControllerOptimizationTarget::Emissions) */ {
            for (unsigned int t = 0; t < T; t++) {
                float current_emissions = Global::get_emissions_g_CO2eq_per_kWh();
                if (global::emission_ts != NULL && ts + t < Global::get_n_timesteps())
                    current_emissions = global::emission_ts[ts - 1 + t];
                objective->SetCoefficient(x_demand_kW[t],  current_emissions);
                // TODO: What about feedin emissions? Ignore? Substract?
            }
        }
        objective->SetMinimization();
        //
        // Execute the optimization results and check results
        const MPSolver::ResultStatus result_status = model->Solve();
        if (result_status != MPSolver::OPTIMAL) {
            std::cerr << "Optimization not resulting in optimal value.\n";
            std::cerr << "Solver status = " << result_status << std::endl;
            return false;
        }
        //
        // Get the results and store them in the output reference
        for (unsigned int t = 0; t < T; t++) {
            bs_power[t] = p_bs_in_kW[t]->solution_value() - p_bs_out_kW[t]->solution_value();
            hp_power[t] = p_hp_kW[t]->solution_value();
            for (unsigned long evIdx = 0; evIdx < n_cars; evIdx++) {
                ev_power[evIdx][t] = p_ev_kW[evIdx][t]->solution_value();
            }
        }
        if (optimize_PV_size) {
            optimal_pv_size_per_section_kW.clear();
            for (unsigned long sectionIdx = 0; sectionIdx < n_pv_sections; sectionIdx++) {
                optimal_pv_size_per_section_kW.push_back( dim_pv_per_sec_kW->solution_value() );
            }
        }
        if (optimize_BS_size) {
            optimal_bs_size_kWh = dim_bs_E_kWh->solution_value();
        }
        //
        // Cleanup
        // done by the model itself: if (dim_bs_E_kWh != NULL) delete dim_bs_E_kWh;
        // done by the model itself: if (dim_bs_P_kW  != NULL) delete dim_bs_P_kW;
        delete model;
        return true;
    }

};

#endif
