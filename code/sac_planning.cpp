#include "sac_planning.h"

using namespace expansion;


#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
/*#include <list>*/
#include <random>
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


double add_expansion_to_units_random_or_data_order(
    long  expansion_matrix_abs_freq[16][16],
    vector<vector<ControlUnit*>>& cuRefLstVectBitOrder
) {
    // cummulative sum of added kWp of residential PV nominal power in (kWp)
    double cumsum_added_pv_kWp = 0.0;
    //
    for (unsigned int iMatO = 0; iMatO < 16; iMatO++) {
        int iBitO = expCombiMatrixOrderToBitRepr( iMatO ); // get index in Bitwise Order (BitO)
        vector<ControlUnit*>* listOfCUs = &(cuRefLstVectBitOrder[ iBitO ]);
        if (Global::get_cu_selection_mode_fca() == global::CUSModeFCA::RandomSelection) {
            //
            // shuffle list if CU selection mode for comp. add. tells so (or random anyway is selected)
            random_device rndDevice;
            mt19937 rndGen( rndDevice() );
            shuffle(listOfCUs->begin(), listOfCUs->end(), rndGen);
        }
        // get the iterator
        vector<ControlUnit*>::iterator iter = listOfCUs->begin();
        // loop over all **target** expansion states
        // start at iMatO + 1, as all combinations bevore are impossible / or do not need any expansion
        // If one loops over them as well, the order of the ordered_list (in case it is set) would be useless!
        for (unsigned int jExpTargetMatO = iMatO + 1; jExpTargetMatO < 16; jExpTargetMatO++) {
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
            // loop over this number
            for (long n = 0; n < numThisCombi_i_j; n++) {
                if (iter == listOfCUs->end()) {
                    cerr << "Warning: end of list for expansion reached before all expansion planing were fulfilled." << endl;
                    goto outer_loop_end;
                }
                // 0. check, if max global kWp addition is reached
                if (Global::get_exp_pv_max_kWp_total() >= 0.0 &&
                    expPV &&
                    cumsum_added_pv_kWp >= Global::get_exp_pv_max_kWp_total()) {
                    cout << "Max added pv reached, with " << cumsum_added_pv_kWp << endl;
                    break; // goto outer_loop_end;
                }
                // 0b. if heat pump is added, check, if annual HP consumption is not exceeding addition clip level
                if (expHP) {
                    while ((*iter)->get_annual_hp_el_cons() <= 0 || (*iter)->get_annual_hp_el_cons() > 20000) {
                        // TODO: Make upper clip level configurable!
                        iter++;
                        if (iter == listOfCUs->end()) {
                            cerr << "Warning: end of list for expansion reached before all expansion planing were fulfilled (Pos. 2)." << endl;
                            goto outer_loop_end;
                        }
                    }
                }
                // 1. add components
                if (expPV) (*iter)->add_exp_pv();
                if (expBS) (*iter)->add_exp_bs();
                if (expHP) (*iter)->add_exp_hp();
                if (expEV) (*iter)->add_exp_evchst();
                // 2. if Global::exp_pv_max_kWp_total_set is set, we have to stop if this value has been reached
                if (Global::get_exp_pv_max_kWp_total() >= 0.0) {
                    cumsum_added_pv_kWp += (*iter)->get_sim_comp_pv_kWp();
                }
                // 3. remove from list (would be good, but not required - right now it does not happen)
                //    only the iterator is incremented
                iter++;
            }
        }
        outer_loop_end:;
    }

    return cumsum_added_pv_kWp;
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
double add_expansion_to_units_orderd_by_metric(
    long  expansion_matrix_abs_freq[16][16],
    vector<vector<ControlUnit*>>& cuRefLstVectBitOrder
) {
    // cummulative sum of added kWp of residential PV nominal power in (kWp)
    double cumsum_added_pv_kWp = 0.0;
    list<string*> output_str_collection;
    //
    // Loop over every current expansion / component combination
    for (unsigned int iMatO = 0; iMatO < 16; iMatO++) {
        int iBitO = expCombiMatrixOrderToBitRepr( iMatO ); // get index in Bitwise Order (BitO)
        int iBitRepr = expCombiMatrixOrderToBitRepr(iMatO);
        vector<ControlUnit*>* listOfCUs = &(cuRefLstVectBitOrder[ iBitO ]);
        //
        // shuffle list (for random addition of HP and EV Ch. St.)
        random_device rndDevice;
        mt19937 rndGen( rndDevice() );
        shuffle(listOfCUs->begin(), listOfCUs->end(), rndGen);
        //
        // 1) Add HP and EV Ch. St. to the units (if required)
        vector<ControlUnit*>::iterator iter = listOfCUs->begin(); // get the iterator
        list<ControlUnit*> list_of_CUs_added_HP_only;
        list<ControlUnit*> list_of_CUs_added_EV_only;
        list<ControlUnit*> list_of_CUs_added_HP_EV;
        list<ControlUnit*> list_of_CUs_nothing_added (listOfCUs->begin(), listOfCUs->end()) ; // copy of the existing list (which is a vector actually)
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
                    while ((*iter)->get_annual_hp_el_cons() <= 0 || (*iter)->get_annual_hp_el_cons() > 20000) {
                        // TODO: Make upper clip level configurable!
                        iter++;
                        if (iter == listOfCUs->end()) {
                            cerr << "Warning: end of list for expansion reached before all expansion planing were fulfilled (Pos. 2)." << endl;
                            goto outer_loop_end;
                        }
                    }
                }
                // 1. add components
                if (expHP) (*iter)->add_exp_hp();
                if (expEV) (*iter)->add_exp_evchst();
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
        bool no_error = simulation::runSimulationForOneParamSetting(listOfCUs);
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
            string * metrics_string = cu->get_metrics_string();
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
        /*bool*/ no_error = simulation::runSimulationForOneParamSetting(listOfCUs);
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
            string * metrics_string = cu->get_metrics_string();
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
            // add PV to those units, where PV is better than PV + BS
            while (expPV1 && n1 > n1_done && iter_pv_only != sorted_list_pv_only.end()) {
                iter_pv_only->second->add_exp_pv();
                iter_pv_only++;
                n1_done++;
            }
            // add PV+BS to those units, where PV is better than PV only
            while (n2 > n2_done && iter_pv_and_bs != sorted_list_pv_and_bs.end()) {
                if (expPV2)
                    iter_pv_and_bs->second->add_exp_pv();
                if (expBS2)
                    iter_pv_and_bs->second->add_exp_bs();
                iter_pv_and_bs++;
                n2_done++;
            }
            // check, if we still have to add units to combinations, that are not locally optimal
            // A) PV
            while (expPV1 && n1 > n1_done && iter_pv_and_bs != sorted_list_pv_and_bs.end()) {
                iter_pv_and_bs->second->add_exp_pv();
                iter_pv_and_bs++;
                n1_done++;
            }
            // B) PV + BS
            while (n2 > n2_done && iter_pv_only != sorted_list_pv_only.end()) {
                if (expPV2)
                    iter_pv_and_bs->second->add_exp_pv();
                if (expBS2)
                    iter_pv_and_bs->second->add_exp_bs();
                iter_pv_only++;
                n2_done++;
            }
            // TODO: Respect total addition limits for PV
            // increment combination_idx
            combination_idx++;
        }
    }
    //
    // output metrics
    output::outputMetricsStrList(output_str_collection);
    for (string* s : output_str_collection) delete s;

    return cumsum_added_pv_kWp;
}


void expansion::add_expansion_to_units(
    float expansion_matrix_rel_freq[16][16],
    long  expansion_matrix_abs_freq[16][16]) {
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
    ControlUnit *const * unit_list = ControlUnit::GetArrayOfInstances();
    const size_t n_CUs = ControlUnit::GetNumberOfInstances();
    for (size_t i = 0; i < n_CUs; i++) {
        ControlUnit* current_unit = unit_list[i];
        // jump this unit, if it is not extensible
        if (!current_unit->is_expandable_with_pv_hp())
            continue;
        //
        int expCombi = current_unit->get_exp_combi_bit_repr();
        currExpCountsBitIndexed[ expCombi ]++;
        cuRefLstVectBitOrder[ expCombi ].push_back( current_unit ); // add unit to reference list
    }
    // change order
    for (int i = 0; i < 16; i++) {
        currExpCountsMatIndexed[i] = currExpCountsBitIndexed[ expCombiMatrixOrderToBitRepr(i) ];
    }

    //
    // 2. calculate expansion matrix with absolute values based
    //    on the abolute counts from 1 by a row-wise multiplication
    //    Note: Impossible values will be jumped, and diagonal values
    //          will be computed afterwards as a difference because
    //          we the absolute numbers have to be integer (half units
    //          are impossible, obviously)
    for (int i = 0; i < 16; i++) {
        for (int j = i+1; j < 16; j++) {
            expansion_matrix_abs_freq[i][j] = (long) ( expansion_matrix_rel_freq[i][j] * (float) (currExpCountsMatIndexed[i]) );
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
    // cummulative sum of added kWp of residential PV nominal power in (kWp)
    double cumsum_added_pv_kWp = 0.0;
    if (Global::get_cu_selection_mode_fca() == global::CUSModeFCA::RandomSelection ||
        Global::get_cu_selection_mode_fca() == global::CUSModeFCA::OrderAsInData)
    {
        cumsum_added_pv_kWp = add_expansion_to_units_random_or_data_order(
            expansion_matrix_abs_freq,
            cuRefLstVectBitOrder
        );
    } else /* if (Global::get_cu_selection_mode_fca() == global::CUSModeFCA::BestSSR ||
        Global::get_cu_selection_mode_fca() == global::CUSModeFCA::BestNPV)*/ {
        cumsum_added_pv_kWp = add_expansion_to_units_orderd_by_metric(
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
    cout << "    Global::get_exp_pv_max_kWp_total()         = " << Global::get_exp_pv_max_kWp_total() << "\n";
    if (Global::get_exp_pv_max_kWp_total() >= 0.0)
    cout << "    Total cumsum of kWp                        = " << cumsum_added_pv_kWp << "\n";
    cout << "    ControlUnit::GetNumberOfCUsWithSimCompPV() = " << ControlUnit::GetNumberOfCUsWithSimCompPV() << "\n";
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
    output_per_cu << "UnitID,n_MUs,is_exp_with_pv_hp,pv_orig,pv_added,bs_orig,bs_added,hp_orig,hp_added,evchst_orig,evchst_added,added_pv_kWp,added_bess_E_kWh,added_bess_P_kW,added_hp_AnnECons_kWh" << endl;
    // n_CUs and unit_list defined above, at 1.
    for (unsigned long i = 0; i < n_CUs; i++) {
        ControlUnit* current_unit = unit_list[i];
        int expCombiAsInData    = current_unit->get_exp_combi_bit_repr_from_MUs();
        int expCombiAsSimulated = current_unit->get_exp_combi_bit_repr_sim_added();
        // output information
        output_per_cu <<        current_unit->get_unitID();
        output_per_cu << "," << current_unit->get_n_MUs();
        output_per_cu << "," << current_unit->is_expandable_with_pv_hp();
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
        output_per_cu << "," << ((0 < (expansion::MaskHP & expCombiAsSimulated)) ? current_unit->get_annual_hp_el_cons() : 0.0);
        output_per_cu << "\n";
    }
    output_per_cu.close();
}
