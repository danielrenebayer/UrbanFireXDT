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
    int currExpBitRepr = expCombiBitReprToMatrixOrder(currExpNumber);
    int newExpBitRepr  = expCombiBitReprToMatrixOrder(newExpNumber);
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
                    cout << "Warning: End of line reached before 16 elements are parsed in the expansion scenario file!" << endl;
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

void expansion::add_expansion_to_units(
    float expansion_matrix_rel_freq[16][16],
    long  expansion_matrix_abs_freq[16][16],
    bool  random_anyway_no_output /* = falsee */,
    vector<ControlUnit*>* ordered_list /* = NULL */) {
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
    for (unsigned int iMatO = 0; iMatO < 16; iMatO++) {
        int iBitO = expCombiMatrixOrderToBitRepr( iMatO ); // get index in Bitwise Order (BitO)
        vector<ControlUnit*>* listOfCUs = &(cuRefLstVectBitOrder[ iBitO ]);
        vector<ControlUnit*>* helperList = NULL;
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
        } else if (Global::get_cu_selection_mode_fca() == global::CUSModeFCA::RandomSelection || random_anyway_no_output) {
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
            bool expWB = false;
            if (ijXOR & MaskPV) expPV = true;
            if (ijXOR & MaskBS) expBS = true;
            if (ijXOR & MaskHP) expHP = true;
            if (ijXOR & MaskWB) expWB = true;
            // loop over this number
            for (long n = 0; n < numThisCombi_i_j; n++) {
                if (iter == listOfCUs->end()) {
                    cerr << "Warning: end of list for expansion reached before all expansion planing were fulfilled." << endl;
                    goto outer_loop_end;
                }
                // 1. add components
                if (expPV) (*iter)->add_exp_pv();
                if (expBS) (*iter)->add_exp_bs();
                if (expHP) (*iter)->add_exp_hp();
                if (expWB) (*iter)->add_exp_wb();
                // 2. remove from list (would be good, but not required)
                iter++;
            }
        }
        outer_loop_end:;
        if (helperList != NULL)
            delete helperList;
    }

    // exit, if no output is selected
    if (random_anyway_no_output)
        return;
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
    output_per_cu << "UnitID,n_MUs,pv_orig,pv_added,bs_orig,bs_added,hp_orig,hp_added,wb_orig,wb_added,added_pv_kWp,added_bess_E_kWh,added_bess_P_kW" << endl;
    // n_CUs and unit_list defined above, at 1.
    for (unsigned long i = 0; i < n_CUs; i++) {
        ControlUnit* current_unit = unit_list[i];
        int expCombiAsInData    = current_unit->get_exp_combi_bit_repr_from_MUs();
        int expCombiAsSimulated = current_unit->get_exp_combi_bit_repr_sim_added();
        // output information
        output_per_cu <<        current_unit->get_unitID();
        output_per_cu << "," << current_unit->get_n_MUs();
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
        output_per_cu << "\n";
    }
    output_per_cu.close();
}
