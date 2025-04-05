#include "sac_planning.h"

using namespace expansion;


#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <queue>
#include <ranges>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>



using namespace std;



#include "global.h"
#include "units.h"
#include "simulation_logic.h"




int expansion::expCombiMatrixOrderToBitRepr(int indexMatO) {
    /*
     * This function maps an expansion combination number (as orderd in the
     * expansion matrix) to its bitwise representation.
     */
    switch (indexMatO) {
        case  0: return MaskNothing;
        case  1: return MaskPV;
        case  2: return        MaskBS;
        case  3: return               MaskHP;
        case  4: return                      MaskWB;
        case  5: return MaskPV|MaskBS;
        case  6: return MaskPV|       MaskHP;
        case  7: return MaskPV|              MaskWB;
        case  8: return        MaskBS|MaskHP;
        case  9: return        MaskBS|       MaskWB;
        case 10: return               MaskHP|MaskWB;
        case 11: return MaskPV|MaskBS|MaskHP;
        case 12: return MaskPV|MaskBS|       MaskWB;
        case 13: return MaskPV|       MaskHP|MaskWB;
        case 14: return        MaskBS|MaskHP|MaskWB;
        case 15: return MaskPV|MaskBS|MaskHP|MaskWB;
    }
    throw logic_error("Impossible index passed to function!");
    return 0;
}

int expansion::expCombiBitReprToMatrixOrder(int bitRepr) {
    /*
     * This function maps an bitwise representation of the expansion 
     * combination to the number (as orderd in the expansion matrix).
     * It is the inverse to expCombiMatrixOrderToBitRepr()
     */
    if      (bitRepr ==  MaskNothing)
        return  0;
    else if (bitRepr ==  MaskPV)
        return  1;
    else if (bitRepr ==         MaskBS)
        return  2;
    else if (bitRepr ==                MaskHP)
        return  3;
    else if (bitRepr ==                       MaskWB)
        return  4;
    else if (bitRepr == (MaskPV|MaskBS)               )
        return  5;
    else if (bitRepr == (MaskPV|       MaskHP)        )
        return  6;
    else if (bitRepr == (MaskPV|              MaskWB) )
        return  7;
    else if (bitRepr == (       MaskBS|MaskHP)        )
        return  8;
    else if (bitRepr == (       MaskBS|       MaskWB) )
        return  9;
    else if (bitRepr == (              MaskHP|MaskWB) )
        return 10;
    else if (bitRepr == (MaskPV|MaskBS|MaskHP)        )
        return 11;
    else if (bitRepr == (MaskPV|MaskBS|       MaskWB) )
        return 12;
    else if (bitRepr == (MaskPV|       MaskHP|MaskWB) )
        return 13;
    else if (bitRepr == (       MaskBS|MaskHP|MaskWB) )
        return 14;
    else if (bitRepr == (MaskPV|MaskBS|MaskHP|MaskWB) )
        return 15;

    throw logic_error("Invalid bit representation!");
}

int expansion::genExpCombiAsBitRepr(bool has_pv, bool has_bs, bool has_hp, bool has_wb) {
    /*
     * This function returns the binary representation of the expansion
     */
    int retval = 0;
    if (has_pv)
        retval = retval | MaskPV;
    if (has_bs)
        retval = retval | MaskBS;
    if (has_hp)
        retval = retval | MaskHP;
    if (has_wb)
        retval = retval | MaskWB;
    return retval;
}

bool expansion::isExpCombiPossible(int currExpNumber, int newExpNumber) {
    /*
    This function returns if a expansion from currExpNumber to newExpNumber (i.e. the future scenario number) is possible or not
    */
    if (currExpNumber == newExpNumber)
        return false; // this is not a expansion by definition!
    int currExpBitRepr = expCombiMatrixOrderToBitRepr(currExpNumber);
    int newExpBitRepr  = expCombiMatrixOrderToBitRepr(newExpNumber);
    // A combination is possible, if and only if
    //  - the bitwise and of curr and new state returns in the curr state
    // and
    //  - the bitwise or of curr and new state returns in the new state.
    // One can see, that it is enough to check one only.
    return (currExpBitRepr & newExpBitRepr) == currExpBitRepr;
}

std::string expansion::expCombiMatrixOrderToString(int indexMatO) {
    switch (indexMatO) {
        case  0: return std::string("Nothing");
        case  1: return std::string("PV only");
        case  2: return std::string("BS only");
        case  3: return std::string("HP only");
        case  4: return std::string("CS only");
        case  5: return std::string("PV+BS");
        case  6: return std::string("PV+HP");
        case  7: return std::string("PV+CS");
        case  8: return std::string("BS+HP");
        case  9: return std::string("BS+CS");
        case 10: return std::string("HP+CS");
        case 11: return std::string("PV+BS+HP");
        case 12: return std::string("PV+BS+CS");
        case 13: return std::string("PV+HP+CS");
        case 14: return std::string("BS+HP+CS");
        case 15: return std::string("BS+HP+CS");
    }
    throw logic_error("Impossible index passed to function!");
    return std::string("");
}


/*
Loads the expansion matrix into the first argument.

Returns false if an error occurs, else true.
*/
bool expansion::load_expansion_matrix(float expansion_matrix[16][16]) {
    ifstream expmat_input;
    stringstream expmat_input_path;
    expmat_input_path << "./expansion_scenarios/";
    expmat_input_path << setw(4) << setfill('0') << Global::get_expansion_scenario_id();
    expmat_input_path << ".csv";
    expmat_input.open(expmat_input_path.str());
    if (!expmat_input.good()) {
        cerr << "Error when loading the expansion matrix with path " << expmat_input_path.str() << endl;
        return false;
    } else {
        cout << "Opening expansion matrix file " << expmat_input_path.str() << endl;
        string currLineString;
        getline( expmat_input, currLineString ); // jump first line, as this is the header
        for (int i = 0; i < 16; i++) {
            getline( expmat_input, currLineString );
            stringstream currLineStream( currLineString );
            string currLineSplittedElement;
            float current_value;
            getline( currLineStream, currLineSplittedElement, ',' ); // jump first row, as this is the index
            for (int col = 0; col < 16; col++) {
                current_value = 0.0;
                // split this row on the ","
                if ( getline( currLineStream, currLineSplittedElement, ',' ) ) {
                    if (currLineSplittedElement.length() > 0) {
                        // remove \r if it occurs at the end
                        if (currLineSplittedElement.back() == '\r') {
                            currLineSplittedElement.erase(currLineSplittedElement.length()-1, 1);
                        }
                        // convert the value if there still is a element remaining
                        if (currLineSplittedElement.length() > 0)
                            current_value = stof( currLineSplittedElement );
                    }
                } else {
                    cout << "Warning: End of line reached before 16 elements are parsed in line " << i+2 << " in the expansion scenario file!" << endl;
                }
                expansion_matrix[i][col] = current_value;
            }
        }
        expmat_input.close();
    }
    #ifdef DEBUG_EXTRA_OUTPUT
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
            cout << setw(4) << setfill(' ') << expansion_matrix[i][j] << " ";
        }
        cout << endl;
    }
    #endif
    return true;
}




bool expansion::verify_expansion_matrix(float expansion_matrix[16][16]) {
    //
    // Check the expansion matrix
    // 1. 0 for impossible combinations (and diagonal, as this is computed below)
    //    and all values between 0 and 1
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
            if (expansion::isExpCombiPossible(i,j)) {
                // check if between 0 and 1
                if (expansion_matrix[i][j] < 0.0 || expansion_matrix[i][j] > 1.0) {
                    cerr << "Error in expansion matrix: value at position " << i << ", " << j << " is greater 1 or below 0!" << endl;
                    return false;
                }
            } else {
                // check if 0
                if (expansion_matrix[i][j] != 0.0) {
                    cerr << "Error in expansion matrix: value at position " << i << ", " << j << " is " << expansion_matrix[i][j] << " instead of 0!" << endl;
                    return false;
                }
            }
        }
    }
    // 2. calculate the diagonal values (i.e. percentage of unchanged values)
    //    and check if row sum is smaller or equal 1
    for (int i = 0; i < 16; i++) {
        float row_sum = 0.0;
        for (int j = i+1; j < 16; j++) {
            row_sum += expansion_matrix[i][j];
        }
        if (row_sum > 1) {
            cerr << "Error in expansion matrix: row sum > 1 for line " << i << endl;
            return false;
        }
        expansion_matrix[i][i] = 1 - row_sum;
    }

    return true;
}


bool add_expansion_to_units_random_or_data_order(
    unsigned long expansion_matrix_abs_freq[16][16],
    vector<vector<ControlUnit*>>& cuRefLstVectBitOrder
) {
    // cummulative sum of added kWp of residential PV nominal power in (kWp), bs capacity, and so on
    double cumsum_added_pv_kWp = 0.0;
    double cumsum_added_bs_kWh = 0.0;
    double cumsum_added_bs_kW  = 0.0;
    unsigned long cumsum_n_added_hps     = 0;
    unsigned long cumsum_n_added_evchsts = 0;
    unsigned long cumsum_n_added_evs     = 0;
    // are addition limits (if defined) reached?
    bool pv_addition_limit     = false;
    bool bs_addition_limit     = false;
    bool hp_addition_limit     = false;
    bool ev_addition_limit     = false;
    //
    // loop over all combinations (starting in the end to get rid of problems where PV limit is reached but other components should have to be added)
    for (long iMatOlong = 15; iMatOlong >= 0; iMatOlong--) {
        unsigned int iMatO = (unsigned int) iMatOlong; // this is required, as an unsigned int cannot be negative!
    //for (unsigned int iMatO = 0; iMatO < 16; iMatO++) {
        int iBitO = expCombiMatrixOrderToBitRepr( iMatO ); // get index in Bitwise Order (BitO)
        vector<ControlUnit*>* listOfCUs = &(cuRefLstVectBitOrder[ iBitO ]);
        if (Global::get_cu_selection_mode_fca() == global::CUSModeFCA::RandomSelection) {
            //
            // sort the list once according to the location ID, if a seed has been set later
            //     -> helpful, if a seed has been set so that the same control units are selected if 
            //        the simulation is executed for different years where the UnitIDs might change, but the LocationID remains stable
            if (Global::is_seed_set()) {
                std::ranges::sort(
                    *listOfCUs,
                    [](ControlUnit* a, ControlUnit* b) {
                        return a->get_location_id() < b->get_location_id();
                    }
                );
            }
            //
            // shuffle list if CU selection mode for comp. add. tells so (or random anyway is selected)
            random_device rndDevice;
            mt19937 rndGen( rndDevice() );
            if (Global::is_seed_set())
                rndGen.seed(Global::get_seed());
            shuffle(listOfCUs->begin(), listOfCUs->end(), rndGen);
        }
        //
        // convert the vector containint all list of CUs with the current expansion status to a set to be able to remove elements quickly later
        std::set<ControlUnit*> setOfCUs (listOfCUs->begin(), listOfCUs->end());
        // loop over all **target** expansion states
        // start at iMatO + 1, as all combinations bevore are impossible / or do not need any expansion
        // If one loops over them as well, the order of the ordered_list (in case it is set) would be useless!
        for (unsigned int jExpTargetMatO = 15; jExpTargetMatO >= iMatO + 1; jExpTargetMatO--) {
        //for (unsigned int jExpTargetMatO = iMatO + 1; jExpTargetMatO < 16; jExpTargetMatO++) {
            // 
            // get number of CUs that get the current expansion
            long numThisCombi_i_j = expansion_matrix_abs_freq[iMatO][jExpTargetMatO];
            // find out, which units we have to add for this i/j-combination
            int iBitRepr = expCombiMatrixOrderToBitRepr(iMatO);
            int jBitRepr = expCombiMatrixOrderToBitRepr(jExpTargetMatO);
            int ijXOR = iBitRepr ^ jBitRepr;
            bool expPV = false;
            bool expBS = false;
            bool expHP = false;
            bool expEV = false;
            if (ijXOR & MaskPV) expPV = true;
            if (ijXOR & MaskBS) expBS = true;
            if (ijXOR & MaskHP) expHP = true;
            if (ijXOR & MaskWB) expEV = true;
            // filter all elements from the list of control units with the current combination
            // that can be expanded to the current selected target combination
            std::queue<ControlUnit*> currExpandableSetOfCUs;
            for (ControlUnit* cu : setOfCUs) {
                bool include_this_unit = true;
                if (expPV && !cu->is_expandable_with_pv())
                    include_this_unit = false;
                if (expHP && !cu->is_expandable_with_hp())
                    include_this_unit = false;
                if (expBS && !cu->is_expandable_with_pv() && !cu->has_pv()) // if we add a BS, the building MUST contain a PV already or be able to add a PV
                    include_this_unit = false;
                if (expEV &&  cu->get_sim_comp_cs_possible_n_EVs() == 0)
                    include_this_unit = false;
                //
                if (include_this_unit)
                    currExpandableSetOfCUs.push( cu );
            }
            // loop over this number
            for (long n = 0; n < numThisCombi_i_j; n++) {
                if (currExpandableSetOfCUs.size() <= 0) {
                    cerr << "Warning: end of list for expansion reached before all expansion planing were fulfilled." << endl;
                    goto outer_loop_end;
                }
                // 0. check, if max global kWp addition is reached
                if (Global::get_exp_pv_max_kWp_total() >= 0.0 &&
                    expPV &&
                    cumsum_added_pv_kWp >= Global::get_exp_pv_max_kWp_total())
                {
                    pv_addition_limit = true;
                }
                if (Global::get_exp_bess_max_E_total() >= 0.0 &&
                    expBS &&
                    cumsum_added_bs_kWh >= Global::get_exp_bess_max_E_total())
                {
                    bs_addition_limit = true;
                }
                if (Global::get_exp_bess_max_P_total() >= 0.0 &&
                    expBS &&
                    cumsum_added_bs_kW  >= Global::get_exp_bess_max_P_total())
                {
                    bs_addition_limit = true;
                }
                if (Global::is_exp_hp_max_n_addition_set() &&
                    expHP &&
                    cumsum_n_added_hps >= Global::get_exp_hp_max_n_addition())
                {
                    hp_addition_limit = true;
                }
                if (Global::is_exp_ev_max_n_addition_set() &&
                    expEV &&
                    cumsum_n_added_evs >= Global::get_exp_ev_max_n_addition())
                {
                    ev_addition_limit = true;
                }
                // 0b. break inner loop if selected
                if (Global::get_break_sac_loop_if_limit_reached()) {
                    if (pv_addition_limit) {
                        cout << "Maximum of added roof-top PV power reached with " << cumsum_added_pv_kWp << " kWp" << endl;
                        break; // goto outer_loop_end;
                    }
                    if (bs_addition_limit) {
                        cout << "Max added battery storage capacity reached with " << cumsum_added_bs_kWh << " kWh and " << cumsum_added_bs_kW << " kW" << endl;
                        break;
                    }
                    if (hp_addition_limit) {
                        cout << "Max number of added heat pumps reached with " << cumsum_n_added_hps << " heat pumps" << endl;
                        break;
                    }
                    if (ev_addition_limit) {
                        cout << "Max number of added EVs reached with " << cumsum_n_added_evs << " EVs" << endl;
                        break;
                    }
                }
                // 0c. get the current unit
                ControlUnit* cu = currExpandableSetOfCUs.front();
                // 1. add components
                if (expPV && !pv_addition_limit) cu->add_exp_pv();
                if (expHP && !hp_addition_limit) {
                    cu->add_exp_hp();
                    cumsum_n_added_hps += 1;
                }
                if (expEV && !ev_addition_limit) {
                    cu->add_exp_cs();
                    cumsum_n_added_evchsts += 1;
                    cumsum_n_added_evs += cu->get_sim_comp_cs_n_EVs();
                }
                // bs is the last thing to add, as the sizing might depend (if config variable set) on the hp annual consumption size
                if (expBS) { // BS -> Do not install BS if BS addition limit is reached OR if PV addition is reached and this CU has no existing PV
                    if (!bs_addition_limit && (!pv_addition_limit || cu->has_pv()))
                        cu->add_exp_bs();
                }
                // 2. if Global::exp_pv_max_kWp_total_set is set, we have to stop if this value has been reached
                cumsum_added_pv_kWp += cu->get_sim_comp_pv_kWp();
                cumsum_added_bs_kWh += cu->get_sim_comp_bs_E_kWh();
                cumsum_added_bs_kW  += cu->get_sim_comp_bs_P_kW();
                // 3. remove element from both, the set and the list
                currExpandableSetOfCUs.pop();
                setOfCUs.erase( cu );
            }
        }
        outer_loop_end:;
    }
    //
    final_cumsum_of_added_pv_kWp = cumsum_added_pv_kWp;
    final_cumsum_of_added_bs_kWh = cumsum_added_bs_kWh;
    final_cumsum_of_added_bs_kW  = cumsum_added_bs_kW;
    final_count_of_added_hps     = cumsum_n_added_hps;
    final_count_of_added_evchsts = cumsum_n_added_evchsts;
    final_count_of_added_evs     = cumsum_n_added_evs;

    return true;
}

struct Combination {
    bool pv_added;
    bool bs_added;
    bool hp_added;
    bool ev_added;
};


/**
 * Internal helper function
 * for adding expansion to units according to a metric
 * -> Therfore, the metric of every combination has to be computed first
 * 
 * @param expansion_matrix_abs_freq: Expansion matrix with absolute values
 * @param cuRefLstVectBitOrder: List of all control units per expansion combination
 */
bool add_expansion_to_units_orderd_by_metric(
    unsigned long expansion_matrix_abs_freq[16][16],
    vector<vector<ControlUnit*>>& cuRefLstVectBitOrder
) {
    // cummulative sum of added kWp of residential PV nominal power in (kWp), bs capacity and so on
    double cumsum_added_pv_kWp = 0.0;
    double cumsum_added_bs_kWh = 0.0;
    double cumsum_added_bs_kW  = 0.0;
    unsigned long cumsum_n_added_hps     = 0;
    unsigned long cumsum_n_added_evchsts = 0;
    unsigned long cumsum_n_added_evs     = 0;
    bool pv_addition_limit     = false;
    bool bs_addition_limit     = false;
    bool hp_addition_limit     = false;
    bool ev_addition_limit     = false;
    //
    list<string*> output_str_collection; // only a string list required for the output
    auto sort_lambda = [](const pair<double, ControlUnit*> &a, const pair<double, ControlUnit*> &b) { return a.first >= b.first; };
    //
    // Loop over every current expansion / component combination
    // loop over all combinations (starting in the end to get rid of problems where PV limit is reached but other components should have to be added)
    for (long iMatOlong = 15; iMatOlong >= 0; iMatOlong--) {
        unsigned int iMatO = (unsigned int) iMatOlong; // this is required, as an unsigned int cannot be negative!
        string iStrO = expansion::expCombiMatrixOrderToString(iMatO); // used for later output
        int iBitO    = expCombiMatrixOrderToBitRepr( iMatO ); // get index in Bitwise Order (BitO)
        int iBitRepr = expCombiMatrixOrderToBitRepr(iMatO);
        vector<ControlUnit*>* listOfCUs = &(cuRefLstVectBitOrder[ iBitO ]);
        // skip if empty
        if (listOfCUs->size() == 0) {
            cout << "Skipping iMatO " << iMatO << endl;
            continue;
        }
        //
        // create a vector of lists 
        // This represents the following structure:
        // [ ( Target_Combination_Matrix_Order, [ (metric, corresponding CU reference) ] ) ]
        vector< pair< int, list<pair<double, ControlUnit*>>* > > sorted_target_vectors;
        //
        // Loop over every target combination
        for (unsigned int jExpTargetMatO = iMatO + 1; jExpTargetMatO < 16; jExpTargetMatO++) {
            // get number of CUs that get the current expansion
            long numThisCombi_i_j = expansion_matrix_abs_freq[iMatO][jExpTargetMatO];
            #ifdef DEBUG
            cout << "iMatO = " << iMatO << ", jExpTargetMatO = " << jExpTargetMatO << ": numThisCombi_i_j = " << numThisCombi_i_j << endl;
            #endif
            // jump loop, if no HP or EVCh St has to be added
            if (numThisCombi_i_j == 0)
                continue;
            // get current list of (metric, CU reference)
            list<pair<double, ControlUnit*>>* current_target_list = new list<pair<double, ControlUnit*>>();
            sorted_target_vectors.push_back(std::make_pair(jExpTargetMatO,  current_target_list));
            // find out, which units we have to add for this i/j-combination
            int jBitRepr = expCombiMatrixOrderToBitRepr(jExpTargetMatO);
            int ijXOR = iBitRepr ^ jBitRepr;
            bool expPV = false;
            bool expBS = false;
            bool expHP = false;
            bool expCS = false;
            if (ijXOR & MaskPV) expPV = true;
            if (ijXOR & MaskBS) expBS = true;
            if (ijXOR & MaskHP) expHP = true;
            if (ijXOR & MaskWB) expCS = true;
            string jStrO = expansion::expCombiMatrixOrderToString(jExpTargetMatO);
            // filter listOfCUs -> Only select those units where an addition of the selected components is possible
            auto filterLambda = [expPV,expBS,expHP,expCS](ControlUnit* cu){ 
                return 
                    ( (expPV) ? cu->is_expandable_with_pv() : true ) &&
                    ( (expBS) ? cu->is_expandable_with_pv() || cu->has_pv() : true ) &&
                    ( (expHP) ? cu->is_expandable_with_hp() : true ) &&
                    ( (expCS) ? cu->get_sim_comp_cs_possible_n_EVs() > 0 : true );
            };
            auto filteredListOfCUs = *listOfCUs | std::views::filter(filterLambda);
            // add the missing elements
            for (ControlUnit* cu : filteredListOfCUs) {
                if (expPV) cu->add_exp_pv();
                if (expHP) cu->add_exp_hp();
                if (expCS) cu->add_exp_cs();
                if (expBS) cu->add_exp_bs();
            }
            // run simulation for this combination
            cout << "\rSimulation pre-run for all CUs with current configuration " << iStrO << " -> target: " << jStrO << "\n";
            CUControllerThreadGroupManager* tgm = NULL;
            if (Global::get_n_threads() >= 1) {
                tgm = new CUControllerThreadGroupManager(listOfCUs);
                tgm->startAllWorkerThreads();
            }
            bool no_error = simulation::runSimulationForOneParamSetting(tgm, listOfCUs, "    ");
            if (tgm != NULL) {
                delete tgm;
                tgm = NULL;
            }
            if (!no_error) { 
                cerr << "Error during selection of the CUs for adding simulated components." << endl;
                return false;
            }
            // get metric for later sorting
            for (ControlUnit* cu : filteredListOfCUs) {
                double metric_result = (Global::get_cu_selection_mode_fca() == global::CUSModeFCA::BestSSR) ? cu->get_SSR() : cu->get_NPV();
                current_target_list->push_back(std::make_pair(metric_result, cu));
                // add to complete metrics output
                string * metrics_string = cu->get_metrics_string_annual();
                metrics_string->append(",");
                if (expPV) metrics_string->append("+PV");
                if (expBS) metrics_string->append("+BS");
                if (expHP) metrics_string->append("+HP");
                if (expCS) metrics_string->append("+CS");
                output_str_collection.push_back(metrics_string);
            }
            // sort the current_target_list
            current_target_list->sort(sort_lambda);
            // remove the added elements
            for (ControlUnit* cu : *listOfCUs) {
                cu->remove_sim_added_components();
                cu->reset_internal_state();
            }
        }
        //
        // Take always the uppermost element of every list and add the missing elements
        // Also respect addition limits
        size_t n_target_combinations = sorted_target_vectors.size();
        size_t n_already_expanded_CUs[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        // main loop, loop until as many CUs are expanded as it should be
        bool all_targets_reached = false;
        while (!all_targets_reached) {
            // for every target list, remove the uppermost elements that were already expanded
            for (size_t idx = 0; idx < n_target_combinations; idx++) {
                while (sorted_target_vectors[idx].second->size() > 0 &&
                       sorted_target_vectors[idx].second->front().second->is_sim_expanded)
                {
                    sorted_target_vectors[idx].second->pop_front();
                }
            }
            // get the control unit where the metric is maxmizied and expand this as it is required
            size_t max_metric_idx = 0; // the index of the target where the control unit with maximal metric is there
            bool   max_found      = false;
            double max_metric_val = -std::numeric_limits<double>::infinity(); // the maximal value of the metric up to now
            for (size_t idx = 0; idx < n_target_combinations; idx++) {
                if (sorted_target_vectors[idx].second->size() > 0) {
                    double curr_metric_val = sorted_target_vectors[idx].second->front().first;
                    if (max_metric_val < curr_metric_val) {
                        max_metric_val = curr_metric_val;
                        max_metric_idx = idx;
                        max_found      = true;
                    }
                }
            }
            // break loop if no min is found (as there are no elements left)
            if (!max_found)
                break;
            // take the uppermost element, if there are remaining elements
            if (sorted_target_vectors[max_metric_idx].second->size() > 0) {
                // get the current control unit which should be simulatively expanded and remove this from the list
                auto& first_pair = sorted_target_vectors[max_metric_idx].second->front();
                ControlUnit* cu_for_addition = first_pair.second;
                sorted_target_vectors[max_metric_idx].second->pop_front(); // this invalidates the reference first_pair
                // find out, which elements have to be added, Part 1
                int jExpTargetMatO = sorted_target_vectors[max_metric_idx].first;
                // find out, which elements have to be added, Part 2
                int jBitRepr = expCombiMatrixOrderToBitRepr(jExpTargetMatO);
                int ijXOR = iBitRepr ^ jBitRepr;
                bool expPV = false;
                bool expBS = false;
                bool expHP = false;
                bool expCS = false;
                if (ijXOR & MaskPV) expPV = true;
                if (ijXOR & MaskBS) expBS = true;
                if (ijXOR & MaskHP) expHP = true;
                if (ijXOR & MaskWB) expCS = true;
                // are already enough elements added?
                if ( n_already_expanded_CUs[jExpTargetMatO] >= expansion_matrix_abs_freq[iMatO][jExpTargetMatO] ) {
                    // empty this 'column'
                    sorted_target_vectors[max_metric_idx].second->clear();
                    continue;
                }
                // are global limits reached?
                if (Global::get_break_sac_loop_if_limit_reached()) {
                    if (expPV && pv_addition_limit) {
                        continue;
                    }
                    if (expBS && bs_addition_limit) {
                        continue;
                    }
                    if (expHP && hp_addition_limit) {
                        continue;
                    }
                    if (expCS && ev_addition_limit) {
                        continue;
                    }
                }
                // add the missing elements
                if (expPV && !pv_addition_limit) {
                    cu_for_addition->add_exp_pv();
                    cumsum_added_pv_kWp += cu_for_addition->get_sim_comp_pv_kWp(); // calculate new cumsum
                }
                if (expHP && !hp_addition_limit) {
                    cu_for_addition->add_exp_hp();
                    cumsum_n_added_hps += 1;
                }
                if (expCS && !ev_addition_limit) {
                    cu_for_addition->add_exp_cs();
                    cumsum_n_added_evchsts += 1;
                    cumsum_n_added_evs += cu_for_addition->get_sim_comp_cs_n_EVs();
                }
                if (expBS && !bs_addition_limit) {
                    cu_for_addition->add_exp_bs();
                    cumsum_added_bs_kWh += cu_for_addition->get_sim_comp_bs_E_kWh(); // calculate new cumsum
                    cumsum_added_bs_kW  += cu_for_addition->get_sim_comp_bs_P_kW();
                }
                cu_for_addition->is_sim_expanded = true;
                n_already_expanded_CUs[jExpTargetMatO]++;
            }
            // check, if addition limits are reached
            if (Global::get_exp_pv_max_kWp_total() >= 0.0 &&
                cumsum_added_pv_kWp >= Global::get_exp_pv_max_kWp_total())
            {
                pv_addition_limit = true;
            }
            if (Global::get_exp_bess_max_E_total() >= 0.0 &&
                cumsum_added_bs_kWh >= Global::get_exp_bess_max_E_total())
            {
                bs_addition_limit = true;
            }
            if (Global::get_exp_bess_max_P_total() >= 0.0 &&
                cumsum_added_bs_kW  >= Global::get_exp_bess_max_P_total())
            {
                bs_addition_limit = true;
            }
            if (Global::is_exp_hp_max_n_addition_set() &&
                cumsum_n_added_hps >= Global::get_exp_hp_max_n_addition())
            {
                hp_addition_limit = true;
            }
            if (Global::is_exp_ev_max_n_addition_set() &&
                cumsum_n_added_evs >= Global::get_exp_ev_max_n_addition())
            {
                ev_addition_limit = true;
            }
            // check, if all expansion targets are reached?
            all_targets_reached = true;
            for (size_t idx = 0; idx < n_target_combinations; idx++) {
                int jExpTargetMatO = sorted_target_vectors[idx].first;
                if (n_already_expanded_CUs[jExpTargetMatO] < expansion_matrix_abs_freq[iMatO][jExpTargetMatO])
                    all_targets_reached = false;
            }
        }
        // delete list
        for (auto &obj : sorted_target_vectors) {
            delete obj.second;
        }
    }
    //
    // output metrics
    output::outputMetricsStrListSACPlanning(output_str_collection);
    for (string* s : output_str_collection) delete s;
    // make nice output
    cout << "\r" << global::output_section_delimiter << std::endl;
    //
    final_cumsum_of_added_pv_kWp = cumsum_added_pv_kWp;
    final_cumsum_of_added_bs_kWh = cumsum_added_bs_kWh;
    final_cumsum_of_added_bs_kW  = cumsum_added_bs_kW;
    final_count_of_added_hps     = cumsum_n_added_hps;
    final_count_of_added_evchsts = cumsum_n_added_evchsts;
    final_count_of_added_evs     = cumsum_n_added_evs;

    return true;
}
double add_expansion_to_units_orderd_by_metric_OLD(
    unsigned long expansion_matrix_abs_freq[16][16],
    vector<vector<ControlUnit*>>& cuRefLstVectBitOrder
) {
    // cummulative sum of added kWp of residential PV nominal power in (kWp)
    double cumsum_added_pv_kWp = 0.0;
    double cumsum_added_bs_kWh = 0.0;
    list<string*> output_str_collection;
    //
    // Loop over every current expansion / component combination
    // loop over all combinations (starting in the end to get rid of problems where PV limit is reached but other components should have to be added)
    for (long iMatOlong = 15; iMatOlong >= 0; iMatOlong--) {
        unsigned int iMatO = (unsigned int) iMatOlong; // this is required, as an unsigned int cannot be negative!
        string iStrO = expansion::expCombiMatrixOrderToString(iMatO); // used for later output
    //for (unsigned int iMatO = 0; iMatO < 16; iMatO++) {
        int iBitO = expCombiMatrixOrderToBitRepr( iMatO ); // get index in Bitwise Order (BitO)
        int iBitRepr = expCombiMatrixOrderToBitRepr(iMatO);
        vector<ControlUnit*>* listOfCUs = &(cuRefLstVectBitOrder[ iBitO ]);
        //
        // shuffle list (for random addition of HP and EV Ch. St.)
        random_device rndDevice;
        mt19937 rndGen( rndDevice() );
        if (Global::is_seed_set())
            rndGen.seed(Global::get_seed());
        shuffle(listOfCUs->begin(), listOfCUs->end(), rndGen);
        //
        // 1) Add HP and EV Ch. St. to the units (if required)
        vector<ControlUnit*>::iterator iter = listOfCUs->begin(); // get the iterator
        list<ControlUnit*> list_of_CUs_added_HP_only;
        list<ControlUnit*> list_of_CUs_added_EV_only;
        list<ControlUnit*> list_of_CUs_added_HP_EV;
        list<ControlUnit*> list_of_CUs_nothing_added (listOfCUs->begin(), listOfCUs->end()) ; // copy of the existing list (which is a vector actually)
        //for (unsigned int jExpTargetMatO = 15; jExpTargetMatO >= iMatO + 1; jExpTargetMatO--) {
        for (unsigned int jExpTargetMatO = iMatO + 1; jExpTargetMatO < 16; jExpTargetMatO++) {
            // 
            // get number of CUs that get the current expansion
            long numThisCombi_i_j = expansion_matrix_abs_freq[iMatO][jExpTargetMatO];
            // find out, which units we have to add for this i/j-combination
            int jBitRepr = expCombiMatrixOrderToBitRepr(jExpTargetMatO);
            int ijXOR = iBitRepr ^ jBitRepr;
            //bool expPV = false;
            //bool expBS = false;
            bool expHP = false;
            bool expEV = false;
            //if (ijXOR & MaskPV) expPV = true;
            //if (ijXOR & MaskBS) expBS = true;
            if (ijXOR & MaskHP) expHP = true;
            if (ijXOR & MaskWB) expEV = true;
            //
            // jump loop, if no HP or EVCh St has to be added
            if (numThisCombi_i_j == 0 || !(expHP || expEV)) {
                continue;
            }
            //
            // loop over this number
            for (long n = 0; n < numThisCombi_i_j; n++) {
                if (iter == listOfCUs->end()) {
                    cerr << "Warning: end of list for expansion reached before all expansion planing were fulfilled." << endl;
                    goto outer_loop_end;
                }
                // 0b. if heat pump is added, check, if annual HP consumption is not exceeding addition clip level
                if (expHP) {
                    while ((*iter)->get_annual_hp_el_cons_kWh() <= 0) {
                        iter++;
                        if (iter == listOfCUs->end()) {
                            cerr << "Warning: end of list for expansion reached before all expansion planing were fulfilled (Pos. 2)." << endl;
                            goto outer_loop_end;
                        }
                    }
                }
                // 1. add components
                if (expHP) (*iter)->add_exp_hp();
                if (expEV) (*iter)->add_exp_cs();
                // 2. add element to correct lists (and remove from nothing_added-list)
                if      ( expHP || !expEV) list_of_CUs_added_HP_only.push_back(*iter);
                else if (!expHP ||  expEV) list_of_CUs_added_EV_only.push_back(*iter);
                else list_of_CUs_added_HP_EV.push_back(*iter);
                list_of_CUs_nothing_added.remove(*iter);
                // 3. increment iterator
                iter++;
            }
        }
        outer_loop_end:;
        //
        // 2) Add PV and PV+BS to every unit (with this combination) and compute required metrics
        //ControlUnit*const* cuList = ControlUnit::GetArrayOfInstances();
        //const size_t nCUs = ControlUnit::GetNumberOfInstances();
        // maps for storage
        map<ControlUnit*, pair<double, double>> metrics_no_HP_EV; // first argument: metric only with PV, second metric with PV and BS
        map<ControlUnit*, pair<double, double>> metrics_HP_only;
        map<ControlUnit*, pair<double, double>> metrics_EV_only;
        map<ControlUnit*, pair<double, double>> metrics_HP_EV;
        // get references to the maps for saving written code
        auto combinations = list {
            std::make_pair(&list_of_CUs_nothing_added, &metrics_no_HP_EV),
            std::make_pair(&list_of_CUs_added_HP_only, &metrics_HP_only),
            std::make_pair(&list_of_CUs_added_EV_only, &metrics_EV_only),
            std::make_pair(&list_of_CUs_added_HP_EV,   &metrics_HP_EV)
        };
        int combination_bitrepr[4] = { // same order as for 'combinations' defined above
            0,
            MaskHP,
            MaskWB,
            MaskHP | MaskWB
        };
        // initialize all values with (0,0)
        for ( auto m : combinations ) {
            auto m1 = *(m.first);  // list of control units with this combi
            auto m2 = *(m.second); // map of CUs to pair of metric results
            //
            for (ControlUnit* cu : m1) {
                m2[cu] = std::make_pair(0.0, 0.0);
            }
        }
        // 2.1) Case PV only
        // 2.1.1) Add PV only
        list<ControlUnit*> list_of_CUs_added_PV; // list of CUs with added PV (only for testing, thus it has to be removed later again)
        for (ControlUnit* cu : *listOfCUs) {
            if (!cu->has_pv()) {
                cu->add_exp_pv();
                list_of_CUs_added_PV.push_back(cu);
            }
        }
        // 2.1.2) Execute the simulation once
        cout << "\rSimulation pre-run for all CUs with current configuration " << iStrO << " -> target: +PV, no BS \n";
        CUControllerThreadGroupManager* tgm = NULL;
        if (Global::get_n_threads() >= 1) {
            tgm = new CUControllerThreadGroupManager(listOfCUs);
            tgm->startAllWorkerThreads();
        }
        bool no_error = simulation::runSimulationForOneParamSetting(tgm, listOfCUs, "    ");
        if (tgm != NULL) {
            delete tgm;
            tgm = NULL;
        }
        if (!no_error) { 
            cerr << "Error during selection of the CUs for adding simulated components." << endl;
            return false;
        }
        // 2.1.3) Collect the SSR (or NPV) values
        for ( auto m : combinations ) {
            for (ControlUnit* cu : *m.first) {
                double metric_result = (Global::get_cu_selection_mode_fca() == global::CUSModeFCA::BestSSR) ? cu->get_SSR() : cu->get_NPV();
                (*m.second)[cu].first = metric_result; // use first parameter for storing PV only metric value
            }
        }
        // 2.1.4) Output metric values
        for (ControlUnit* cu : *listOfCUs) {
            string * metrics_string = cu->get_metrics_string_annual();
            metrics_string->append(",PV only");
            output_str_collection.push_back(metrics_string);
        }
        // 2.1.5) Reset internal variables
        ControlUnit::ResetAllInternalStates();
        //
        // 2.2) Case PV and BS
        // 2.2.1) Add BS as well
        list<ControlUnit*> list_of_CUs_added_BS; // list of CUs with added BS (only for testing, thus it has to be removed later again)
        for (ControlUnit* cu : *listOfCUs) {
            if ( !cu->has_bs() /* && ( cu->get_exp_combi_bit_repr_sim_added() & expansion::MaskBS ) */ ) {
                cu->add_exp_bs();
                list_of_CUs_added_BS.push_back(cu);
            }
        }
        // 2.2.2) Execute the simulation once
        cout << "\rSimulation pre-run for all CUs with current configuration " << iStrO << " -> target: +PV+BS \n";
        if (Global::get_n_threads() >= 1) {
            tgm = new CUControllerThreadGroupManager(listOfCUs);
            tgm->startAllWorkerThreads();
        }
        /*bool*/ no_error = simulation::runSimulationForOneParamSetting(tgm, listOfCUs, "    ");
        if (tgm != NULL) {
            delete tgm;
            tgm = NULL;
        }
        if (!no_error) { 
            cerr << "Error during selection of the CUs for adding simulated components." << endl;
            return false;
        }
        // 2.2.3) Collect the SSR (or NPV) values
        for ( auto m : combinations ) {
            for (ControlUnit* cu : *m.first) {
                double metric_result = (Global::get_cu_selection_mode_fca() == global::CUSModeFCA::BestSSR) ? cu->get_SSR() : cu->get_NPV();
                (*m.second)[cu].second = metric_result; // use second parameter instead of first as some lines above
            }
        }
        // 2.2.4) Output metric values
        for (ControlUnit* cu : *listOfCUs) {
            string * metrics_string = cu->get_metrics_string_annual();
            metrics_string->append(",PV and BS");
            output_str_collection.push_back(metrics_string);
        }
        // 2.2.5) Reset internal variables
        ControlUnit::ResetAllInternalStates();
        // 2.3) Remove all PV and BS components
        for (ControlUnit* cu : list_of_CUs_added_PV) {
            cu->remove_sim_added_pv();
        }
        for (ControlUnit* cu : list_of_CUs_added_BS) {
            cu->remove_sim_added_bs();
        }
        //
        // 3) Select units with best SSR or NPV
        auto sort_lambda = [](const pair<double, ControlUnit*> &a, const pair<double, ControlUnit*> &b) { return a.first > b.first; };
        // Loop over all combinations (no HP and EV, HP only, EV only, HP + EV)
        unsigned int combination_idx = 0; // required for getting the number of elements to add
        for ( auto m : combinations ) {
            auto m1 = *m.first;  // list of control units with this combi
            auto m2 = *m.second; // map of CUs to pair of metric results
            // jump this combination, if emtpy
            if (m1.size() == 0)
                continue;
            // Decide what is better per CU (add PV or PV+BS)?
            //map<ControlUnit*, bool> internally_better_combi; // ture -> only PV, false -> PV + BS
            list<pair<double, ControlUnit*>> sorted_list_pv_only;   // elements, where metric with PV only is better than with PV and BS
            list<pair<double, ControlUnit*>> sorted_list_pv_and_bs; // vice versa
            for (auto cu_metric_pair : m2) {
                bool res = cu_metric_pair.second.first > cu_metric_pair.second.second;
                //internally_better_combi[cu_metric_pair.first] = res;
                if (res) {
                    sorted_list_pv_only.emplace_back(cu_metric_pair.second.first, cu_metric_pair.first);
                } else {
                    sorted_list_pv_and_bs.emplace_back(cu_metric_pair.second.second, cu_metric_pair.first);
                }
            }
            // get sorted lists for PV / PV+BS addition
            //sort(  sorted_list_pv_only.begin(),   sorted_list_pv_only.end(), sort_lambda);
            //sort(sorted_list_pv_and_bs.begin(), sorted_list_pv_and_bs.end(), sort_lambda);
            sorted_list_pv_only.sort(sort_lambda);
            sorted_list_pv_and_bs.sort(sort_lambda);
            auto iter_pv_only   =   sorted_list_pv_only.begin();
            auto iter_pv_and_bs = sorted_list_pv_and_bs.begin();
            // get number of PV / PV + BS to add
            unsigned int jExpTargetMatO1 = expCombiBitReprToMatrixOrder(combination_bitrepr[combination_idx] | MaskPV);
            unsigned int jExpTargetMatO2 = expCombiBitReprToMatrixOrder(combination_bitrepr[combination_idx] | MaskPV | MaskBS);
            long n1 = expansion_matrix_abs_freq[iMatO][jExpTargetMatO1]; // number of pv (no bs) to add
            long n2 = expansion_matrix_abs_freq[iMatO][jExpTargetMatO2]; // number of pv and bs to add
            long n1_done = 0;
            long n2_done = 0;
            int jBitRepr1 = expCombiMatrixOrderToBitRepr(jExpTargetMatO1);
            int jBitRepr2 = expCombiMatrixOrderToBitRepr(jExpTargetMatO2);
            bool expPV1 = (iBitRepr ^ jBitRepr1) & MaskPV;
            bool expPV2 = (iBitRepr ^ jBitRepr2) & MaskPV;
            bool expBS2 = (iBitRepr ^ jBitRepr2) & MaskBS;
            // limits for addition reached?
            bool pv_add_limit = false;
            bool bs_add_limit = false;
            // add PV+BS to those units, where PV+BS is better than PV only
            while (n2 > n2_done && iter_pv_and_bs != sorted_list_pv_and_bs.end()) {
                if (expPV2) {
                    iter_pv_and_bs->second->add_exp_pv();
                    cumsum_added_pv_kWp += iter_pv_and_bs->second->get_sim_comp_pv_kWp(); // calculate new cumsum
                }
                if (expBS2) {
                    iter_pv_and_bs->second->add_exp_bs();
                    cumsum_added_bs_kWh += iter_pv_and_bs->second->get_sim_comp_bs_E_kWh();
                }
                // increment iterators
                iter_pv_and_bs++;
                n2_done++;
                // test if addition limits are reached
                if (Global::get_exp_pv_max_kWp_total() >= 0.0 &&
                    cumsum_added_pv_kWp >= Global::get_exp_pv_max_kWp_total()) {
                        pv_add_limit = true;
                        break;
                }
                if (Global::get_exp_bess_max_E_total() >= 0.0 &&
                    cumsum_added_bs_kWh >= Global::get_exp_bess_max_E_total())
                {
                    bs_add_limit = true;
                    break;
                }
            }
            // add PV to those units, where PV is better than PV + BS
            while (expPV1 && n1 > n1_done && iter_pv_only != sorted_list_pv_only.end() && !pv_add_limit) {
                iter_pv_only->second->add_exp_pv();
                cumsum_added_pv_kWp += iter_pv_only->second->get_sim_comp_pv_kWp(); // calculate new cumsum
                // increment iterators
                iter_pv_only++;
                n1_done++;
                // test if addition limits are reached
                if (Global::get_exp_pv_max_kWp_total() >= 0.0 &&
                    cumsum_added_pv_kWp >= Global::get_exp_pv_max_kWp_total()) {
                        pv_add_limit = true;
                        break;
                }
            }
            // check, if we still have to add units to combinations, that are not locally optimal
            // we add them to the other list, as this situation can only occure, if one of the lists is empty
            // A) PV (using the list of units where PV+BS would be better, i.e. iter_pv_and_bs)
            while (expPV1 && n1 > n1_done && iter_pv_and_bs != sorted_list_pv_and_bs.end() && !pv_add_limit) {
                iter_pv_and_bs->second->add_exp_pv();
                cumsum_added_pv_kWp += iter_pv_and_bs->second->get_sim_comp_pv_kWp(); // calculate new cumsum
                // increment iterators
                iter_pv_and_bs++;
                n1_done++;
                // test if addition limits are reached
                if (Global::get_exp_pv_max_kWp_total() >= 0.0 &&
                    cumsum_added_pv_kWp >= Global::get_exp_pv_max_kWp_total()) {
                        pv_add_limit = true;
                        break;
                }
            }
            // B) PV + BS (using the list where PV would be better, i.e. iter_pv)
            while (n2 > n2_done && iter_pv_only != sorted_list_pv_only.end() && !pv_add_limit && !bs_add_limit) {
                if (expPV2) {
                    iter_pv_only->second->add_exp_pv();
                    cumsum_added_pv_kWp += iter_pv_only->second->get_sim_comp_pv_kWp(); // calculate new cumsum
                }
                if (expBS2) {
                    iter_pv_only->second->add_exp_bs();
                    cumsum_added_bs_kWh += iter_pv_and_bs->second->get_sim_comp_bs_E_kWh();
                }
                // increment iterators
                iter_pv_only++;
                n2_done++;
                // test if addition limits are reached
                if (Global::get_exp_pv_max_kWp_total() >= 0.0 &&
                    cumsum_added_pv_kWp >= Global::get_exp_pv_max_kWp_total()) {
                        pv_add_limit = true;
                        break;
                }
                if (Global::get_exp_bess_max_E_total() >= 0.0 &&
                    cumsum_added_bs_kWh >= Global::get_exp_bess_max_E_total())
                {
                    bs_add_limit = true;
                    break;
                }
            }
            // increment combination_idx
            combination_idx++;
        }
    }
    //
    // output metrics
    output::outputMetricsStrListSACPlanning(output_str_collection);
    for (string* s : output_str_collection) delete s;
    // make nice output
    cout << "\r" << global::output_section_delimiter << std::endl;

    return cumsum_added_pv_kWp;
}


void expansion::add_expansion_to_units(
    float expansion_matrix_rel_freq[16][16],
    unsigned long expansion_matrix_abs_freq[16][16]) {
    /*
     * This function adds the expansion given als relative counts in expansion_matrix_rel_freq
     * to the control units.
     * @param expansion_matrix_abs_freq contains the absolute counts afterwards
     * @param random_anyway_no_output: do random selection (whatever the Global::get_cu_selection_mode_fca() says) but switch off output, if set to true
     * @param ordered_list: if given, the order given in the list is taken (in this case, random selection is always neglected, regardeless of what random_anyway_no_output says)
     */

    long currExpCountsBitIndexed[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // order as given from bitwise representation (as this represents a sequential integer as well)
    long currExpCountsMatIndexed[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // order as in expansion matrix
    long newExpCountsMatIndexed[16]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    vector<vector<ControlUnit*>> cuRefLstVectBitOrder(16); // vector of 16 lists (actually also vectors), that contains the references to the CUs, where the index of the list in the vector corresponds to the expansion in bitwise order

    //
    // 1. count current expansion status (as it is in the given data)
    //    and store the references to the CUs for a given current expansion in a list
    const std::vector<ControlUnit*>& units_list = ControlUnit::GetArrayOfInstances();
    for (ControlUnit* current_unit : units_list) {
        /*
        // jump this unit, if it is not extensible (only if exp pv static mode is not available)
        if (!Global::get_exp_pv_static_mode() && !current_unit->is_expandable_with_pv_hp())
            continue;
        // exlude units with unknwon heat demand if requested
        if (Global::get_annual_heat_demand_limit_fsac() >= 1.0 &&
            current_unit->get_annual_heat_demand_th_kWh() > Global::get_annual_heat_demand_limit_fsac() ) {
                continue;
            }
        //
        // Only select units with known gas consumption
        if (Global::get_select_buildings_wg_heatd_only() && !current_unit->heat_demand_given_in_data())
            continue;
        */
        // Select Control Unit only if the number of added EVs would not exceed the individual limit (if it is set)
        if (Global::get_exp_cs_max_ev_per_cs() > 0) {
            if (current_unit->get_sim_comp_cs_possible_n_EVs() > Global::get_exp_cs_max_ev_per_cs())
                continue;
        }
        // select only residential buildings if selected
        if (Global::get_select_only_residential_buildings()) {
            if (!current_unit->is_residential())
                continue;
        }
        //
        int expCombi = current_unit->get_exp_combi_bit_repr();
        currExpCountsBitIndexed[ expCombi ]++;
        cuRefLstVectBitOrder[ expCombi ].push_back( current_unit ); // add unit to reference list
    }
    // change order
    for (int i = 0; i < 16; i++) {
        currExpCountsMatIndexed[i] = currExpCountsBitIndexed[ expCombiMatrixOrderToBitRepr(i) ];
    }
    // output delimiter for nice output
    std::cout << global::output_section_delimiter << std::endl;

    //
    // 2. calculate expansion matrix with absolute values based
    //    on the abolute counts from 1 by a row-wise multiplication
    //    Note: Impossible values will be jumped, and diagonal values
    //          will be computed afterwards as a difference because
    //          we the absolute numbers have to be integer (half units
    //          are impossible, obviously)
    for (int i = 0; i < 16; i++) {
        for (int j = i+1; j < 16; j++) {
            expansion_matrix_abs_freq[i][j] = (unsigned long) ( expansion_matrix_rel_freq[i][j] * (float) (currExpCountsMatIndexed[i]) );
        }
    }
    // calculate diagonal values
    for (int i = 0; i < 16; i++) {
        long sum_expanded_i = 0;
        for (int j = i+1; j < 16; j++){
            sum_expanded_i += expansion_matrix_abs_freq[i][j];
        }
        expansion_matrix_abs_freq[i][i] = currExpCountsMatIndexed[i] - sum_expanded_i;
    }

    //
    // 3. count new expansion status, that will be simulated now
    //    This is the column sum of the absolute expansion matrix
    for (int j = 0; j < 16; j++) { // sum over j (cols) instead of i (rows) first
        for (int i = 0; i < 16; i++) {
            newExpCountsMatIndexed[j] += expansion_matrix_abs_freq[i][j];
        }
    }

    //
    // 4. plan and execute expansion
    //
    if (Global::get_cu_selection_mode_fca() == global::CUSModeFCA::RandomSelection ||
        Global::get_cu_selection_mode_fca() == global::CUSModeFCA::OrderAsInData)
    {
        add_expansion_to_units_random_or_data_order(
            expansion_matrix_abs_freq,
            cuRefLstVectBitOrder
        );
    } else /* if (Global::get_cu_selection_mode_fca() == global::CUSModeFCA::BestSSR ||
        Global::get_cu_selection_mode_fca() == global::CUSModeFCA::BestNPV)*/ {
        // TODO process return value (false -> simulation error)
        add_expansion_to_units_orderd_by_metric(
            expansion_matrix_abs_freq,
            cuRefLstVectBitOrder
        );
    }

    /*
    
        if (ordered_list != NULL) {
            //
            // if ordered_list is given, this will be used
            // use the helper list that holds the filter elements of listOfCUs ...
            // ... filtering is acually required, as we have to filter the elements with a given expansion at that point
            helperList = new vector<ControlUnit*>();
            helperList->reserve(listOfCUs->size());
            for (ControlUnit* c : *ordered_list) {
                if ( find( listOfCUs->begin(), listOfCUs->end(), c) != listOfCUs->end() ) {
                    // listOfCUs contains c, so we add s to the new list
                    helperList->push_back(c);
                }
            }
            listOfCUs = helperList;
        }
    */
    

    //
    // additional output
    cout << "Overview of simulatively added components:\n";
    if (Global::get_exp_pv_max_kWp_total() >= 0.0)
    cout << "    Global::get_exp_pv_max_kWp_total()         = " << Global::get_exp_pv_max_kWp_total() << "\n";
    if (Global::get_exp_bess_max_E_total() >= 0.0)
    cout << "    Global::get_exp_bess_max_E_total()         = " << Global::get_exp_bess_max_E_total() << "\n";
    cout << "    Global::get_exp_bess_max_P_total()         = " << Global::get_exp_bess_max_P_total() << "\n";
    cout << "    Total cumsum of added PV kWp               = " << final_cumsum_of_added_pv_kWp << "\n";
    cout << "    Total cumsum of added BS capacity in kWh   = " << final_cumsum_of_added_bs_kWh << "\n";
    cout << "    Total cumsum of added BS power    in kW    = " << final_cumsum_of_added_bs_kW  << "\n";
    cout << "    Total number of added heat pumps           = " << final_count_of_added_hps << "\n";
    cout << "    Total number of added EV charging stations = " << final_count_of_added_evchsts << "\n";
    cout << "    Total number of added EVs                  = " << final_count_of_added_evs << "\n";
    cout << "    ControlUnit::GetNumberOfCUsWithSimCompPV() = " << ControlUnit::GetNumberOfCUsWithSimCompPV() << "\n";
    cout << "    ControlUnit::GetNumberOfCUsWithSimCompHP() = " << ControlUnit::GetNumberOfCUsWithSimCompHP() << "\n";
    cout << "    ControlUnit::GetNumberOfCUsWithSimCompEV() = " << ControlUnit::GetNumberOfCUsWithSimCompEV() << "\n";
    cout << global::output_section_delimiter << endl;

    //
    // finally: write expansion information to file
    // A. output expansion matrix with absolute numbers
    filesystem::path info_path_A (*(global::current_global_output_dir)); // for the expansion matrix it is acceptable to use the static output path (that does not change over time for different parameter variations, as expansion cannot change)
    info_path_A /= "expansion-matrix-abs-values.csv";
    ofstream output_exp_mat(info_path_A, std::ofstream::out);
    output_exp_mat << ",0. Nothing,1. PV,2. BS,3. HP,4. WB,5. PV+BS,6. PV+HP,7. PV+WB,8. BS+HP,9. BS+WB,10. HP+WB,11. PV+BS+HP,12. PV+BS+WB,13. PV+HP+WB,14. BS+HP+WB,15. PV+BS+HP+WB,Sum as in data" << endl;
    const char * first_column[16] = {"0. Nothing","1. PV","2. BS","3. HP","4. WB","5. PV+BS","6. PV+HP","7. PV+WB","8. BS+HP","9. BS+WB","10. HP+WB","11. PV+BS+HP","12. PV+BS+WB","13. PV+HP+WB","14. BS+HP+WB","15. PV+BS+HP+WB"};
    for (int i = 0; i < 16; i++) {
        output_exp_mat << first_column[i];
        for (int j = 0; j < 16; j++) {
            output_exp_mat << "," << expansion_matrix_abs_freq[i][j];
        }
        output_exp_mat << "," << currExpCountsMatIndexed[i] << endl;
    }
    output_exp_mat << "Sum as simulated";
    for (int i = 0; i < 16; i++)
        output_exp_mat << "," << newExpCountsMatIndexed[i];
    output_exp_mat << "," << endl;
    output_exp_mat.close();
    //
    // B. output information about added components per MeUID
    filesystem::path info_path_B (*(global::current_global_output_dir)); // same argument as 21 lines above
    info_path_B /= "expansion-per-cu.csv";
    ofstream output_per_cu(info_path_B, std::ofstream::out);
    output_per_cu << "UnitID,n_MUs,is_exp_with_pv,is_exp_with_hp,pv_orig,pv_added,bs_orig,bs_added,hp_orig,hp_added,cs_orig,cs_added,added_pv_kWp,added_bess_E_kWh,added_bess_P_kW,added_hp_AnnECons_kWh,added_n_EVs" << endl;
    // units_list defined above, at 1.
    for (ControlUnit* current_unit : units_list) {
        int expCombiAsInData    = current_unit->get_exp_combi_bit_repr_from_MUs();
        int expCombiAsSimulated = current_unit->get_exp_combi_bit_repr_sim_added();
        // output information
        output_per_cu <<        current_unit->get_unitID();
        output_per_cu << "," << current_unit->get_n_MUs();
        output_per_cu << "," << current_unit->is_expandable_with_pv();
        output_per_cu << "," << current_unit->is_expandable_with_hp();
        output_per_cu << "," << (0 < (expansion::MaskPV & expCombiAsInData));
        output_per_cu << "," << (0 < (expansion::MaskPV & expCombiAsSimulated));
        output_per_cu << "," << (0 < (expansion::MaskBS & expCombiAsInData));
        output_per_cu << "," << (0 < (expansion::MaskBS & expCombiAsSimulated));
        output_per_cu << "," << (0 < (expansion::MaskHP & expCombiAsInData));
        output_per_cu << "," << (0 < (expansion::MaskHP & expCombiAsSimulated));
        output_per_cu << "," << (0 < (expansion::MaskWB & expCombiAsInData));
        output_per_cu << "," << (0 < (expansion::MaskWB & expCombiAsSimulated));
        output_per_cu << "," << current_unit->get_sim_comp_pv_kWp();
        output_per_cu << "," << current_unit->get_sim_comp_bs_E_kWh();
        output_per_cu << "," << current_unit->get_sim_comp_bs_P_kW();
        output_per_cu << "," << ((0 < (expansion::MaskHP & expCombiAsSimulated)) ? current_unit->get_annual_hp_el_cons_kWh() : 0.0);
        output_per_cu << "," << current_unit->get_sim_comp_cs_n_EVs();
        output_per_cu << "\n";
    }
    output_per_cu.close();
}
