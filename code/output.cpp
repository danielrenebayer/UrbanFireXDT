#include "output.h"

#include <atomic>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ranges>
#include <sstream>
#include <string>
#include <vector>

#include "global.h"
#include "helper.h"
#include "vehicles.h"

using namespace std;
using namespace output;

//
// internal helper variables
//
bool base_directories_initialized = false;

//
// Helper function (not defined in header file).
// It removes a directory if this already exists
// and creates a new, empty directory
//
void create_dir_del_if_exists(filesystem::path dirpath) {
    // check if dir exists
    if (filesystem::is_directory(dirpath)) {
        // if yes, delete it
        filesystem::remove_all(dirpath);
    }
    // now, actually create output dir (maybe again)
    filesystem::create_directory(dirpath);
}

void output::initializeDirectoriesBase(unsigned long scenario_id) {
    //
    // create a subfolder for the current scenario, if it does not exist
    filesystem::path dirpath = Global::get_output_path();
    if (!filesystem::is_directory(dirpath)) {
        filesystem::create_directory(dirpath);
    }
    // create subfolder for selected scenario (if it does not exist)
    stringstream current_scenario_str;
    current_scenario_str << "S";
    current_scenario_str << setw(4) << setfill('0') << scenario_id;
    dirpath /= current_scenario_str.str();
    if (!filesystem::is_directory(dirpath)) {
        filesystem::create_directory(dirpath);
    }
    // if repetition is selected, create subfolder for repetitions
    if (Global::get_repetitions_selected()) {
        // create subfolder for all repetions as base dir
        dirpath /= "repetitions";
        if (!filesystem::is_directory(dirpath)) {
            filesystem::create_directory(dirpath);
        }
        // create subfolder for current repetition
        dirpath /= to_string(global::current_repetition_counter);
        if (!filesystem::is_directory(dirpath)) {
            filesystem::create_directory(dirpath);
        }
    }
    // copy to global variable
    //delete global::current_global_output_dir;
    global::current_global_output_dir = new filesystem::path(dirpath);
    //
    // create the output dir prefix
    filesystem::path param_vari_path = dirpath;
    stringstream current_vari_name;
    if (Global::is_parameter_variation()) {
        // if parameter variation is selected,
        // create a folder where subfolders are created for the output
        current_vari_name << "param vari ";
        current_vari_name << setw(4) << setfill('0') << Global::get_parameter_varID();
    } else {
        // if no parameter variuation is selected,
        // we just create one folder with "no param vari"
        current_vari_name << "no param vari";
    }
    param_vari_path /= current_vari_name.str();
    create_dir_del_if_exists(param_vari_path);
    // copy to global variable
    global::current_output_dir_prefix = new filesystem::path(param_vari_path);
    //
    // output build information and data information
    ofstream build_info_output( param_vari_path /= "build_and_run_info.txt", std::ofstream::out );
    build_info_output << "Simulation build at " << __DATE__ << " " << __TIME__ <<  "\n";
    #ifdef __GNUC__
    build_info_output << "GCC was used as compiler.\nGCC Version = " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << "\n";
    #endif
    #ifdef __VERSION__
    build_info_output << "Compiler version = " << __VERSION__ << "\n";
    #endif
    #ifdef __OPTIMIZE__
    build_info_output << "Optimization was enabled during compile time.\n";
    #endif
    build_info_output << "C++ standard = " << __cplusplus << "\n\n";
    //
    time_t current_time = time(nullptr);
    build_info_output << "Time of simulation start = " << put_time(localtime(&current_time), "%F %T") << "\n";
    // TODO: output information about when data was preprocessed
    build_info_output.close();
    base_directories_initialized = true;
}

void output::initializeDirectoriesPerPVar() {
    //
    // check, if initializeDirectoriesBase() is already called
    if (!base_directories_initialized) {
        throw logic_error("Error: output::initializeDirectoriesPerPVar() MUST be called after output::initializeDirectoriesBase() has been called.");
    }
    //
    // in case of a parameter variation: create subfolder for
    // the current variation (and delete old ones if the still exist)
    filesystem::path param_vari_path = *(global::current_output_dir_prefix);
    stringstream current_vari_name;
    if (Global::is_parameter_variation()) {
        // create folder for the current parameter variation combination
        current_vari_name << "variation index ";
        current_vari_name << setw(4) << setfill('0') << global::curr_param_vari_combi_index;
        param_vari_path /= current_vari_name.str();
        create_dir_del_if_exists(param_vari_path);
    } else {
        // do nothing, if no parameter variation is selected
    }
    //
    // copy output path to the global variable
    delete global::current_output_dir;
    global::current_output_dir = new filesystem::path(param_vari_path);
}

void output::initializeSubstationOutput(unsigned long scenario_id) {
    // initial check: should a output be created?
    if (!Global::get_create_substation_output())
        return;
    //
    // Part 1: The main file for the substation load time series
    //
    // initialize the output file
    stringstream output_path_subst;
    output_path_subst << setw(4) << setfill('0') << scenario_id;
    output_path_subst << "-substation-time-series.csv";
    filesystem::path output_path = *(global::current_output_dir);
    output_path /= output_path_subst.str();
    substation_output = new ofstream(output_path, std::ofstream::out);
    //
    // add header to output file
    *(substation_output) << "Timestep";
    const auto subList = Substation::GetArrayOfInstances();
    for (Substation* s : subList) {
        *(substation_output) << "," << s->get_name()->c_str();
    }
    *(substation_output) << ",pv_total_generation_kW,pv_feedin_expo_kW,bs_total_generation_kW,bs_feedin_expo_kW,chp_total_generation_kW,chp_feedin_expo_kW,wind_total_generation_kW,wind_feedin_expo_kW,unknown_total_generation_kW,unknown_feedin_expo_kW,total_demand_wo_BS_SC_kW,total_BS_charging_power_kW,OverallBatterySOC,total_load" << endl;
    //
    //
    // Part 2: The secondary file for additional information about the substations
    //
    // initialize the output file
    stringstream output_path_subst2;
    output_path_subst2 << "substation-detailed-time-series.csv";
    filesystem::path output_path2 = *(global::current_output_dir);
    output_path2 /= output_path_subst2.str();
    substation_output_details = new ofstream(output_path2, std::ofstream::out);
    //
    // add header to output file
    *(substation_output_details) << "Timestep";
    for (Substation* s : subList) {
        *(substation_output_details) << "," << s->get_name()->c_str() << "_resident_load_kW";
        *(substation_output_details) << "," << s->get_name()->c_str() << "_resident_demand_kW";
        *(substation_output_details) << "," << s->get_name()->c_str() << "_demand_wo_BESS_or_selfcons_kW";
    }
    *(substation_output_details) << ",total_residential_load,total_residential_demand" << endl;
}



void output::initializeSurplusOutput(unsigned long scenario_id) {
    // initial check: should a output be created?
    if (!Global::get_create_surplus_output())
        return;
    //
    // initialize the output file
    stringstream output_path_subst;
    output_path_subst << setw(4) << setfill('0') << scenario_id;
    output_path_subst << "-surplus-time-series.csv";
    filesystem::path output_path = *(global::current_output_dir);
    output_path /= output_path_subst.str();
    surplus_output = new ofstream(output_path, std::ofstream::out);
    //
    // add header to output file
    *(surplus_output) << "Timestep";
    *(surplus_output) << ",OverallBatterySOC,scheduled_surplus_discharge,scheduled_surplus_charge,actual_surplus_action,charge_request_BESS,load_BESS,BESS_surplus_energy,initial_total_load,total_load" << endl;
}


void output::initializeCUOutput(unsigned long scenario_id) {
    //
    // initializes the CU output, depending on the global setting
    if (Global::get_output_mode_per_cu() == global::OutputModePerCU::SingleFile) {
        // Case 1: One file for all CUs
        //
        cu_single_output = new CUOutputSingleFile(scenario_id);
        //
        // add reference to all CUs
        const std::vector<ControlUnit*>& cuList = ControlUnit::GetArrayOfInstances();
        for (ControlUnit* c : cuList) {
            c->set_output_object(cu_single_output);
        }
    } else if (Global::get_output_mode_per_cu() == global::OutputModePerCU::IndividualFile) {
        // Case 2: One file for each substation
        //
        // create output directory, and delete existing if present
        filesystem::path dirpath = *(global::current_output_dir);
        dirpath /= "ts-per-cu";
        if (filesystem::is_directory(dirpath)) {
            // clear existing directory
            filesystem::remove_all(dirpath);
        }
        // now, actually create output dir (maybe again)
        filesystem::create_directory(dirpath);
        //
        // create a output object for every (future) substation
        // and add reference to all connected CUs
        n_cu_multi_outputs = Global::get_n_substations();
        cu_multi_outputs = new CUOutputOneFilePerSubstation*[n_cu_multi_outputs];
        const auto substationList = Substation::GetArrayOfInstances();
        for (auto [subst_idx, currentS] : std::ranges::views::enumerate(substationList)) {
            cu_multi_outputs[subst_idx] = new CUOutputOneFilePerSubstation(currentS->get_name(), dirpath);
            // add output file to all connected CUs
            const list<ControlUnit*>* conn_units = currentS->get_connected_units();
            for (ControlUnit* cu : *conn_units) {
                cu->set_output_object(cu_multi_outputs[subst_idx]);
            }
        }
    } else {
        // Case 3: no output per CU selected
    }
    //
    // initialize additional outputs with details about the control units (if selected)
    if (Global::get_create_control_cmd_output()) {
        filesystem::path output_path = *(global::current_output_dir);
        output_path /= "CU-control-commands.json";
        cu_details_ccmd_output = new ofstream(output_path, std::ofstream::out);
        // activate buffers for speedup
        //buffer = new char[bufferSize];
        //cu_details_ccmd_output->rdbuf()->pubsetbuf(buffer, bufferSize);
        // add header to output file
        *(cu_details_ccmd_output) << "[\n";
    }
    if (Global::get_create_ev_detailed_output()) {
        filesystem::path output_path = *(global::current_output_dir);
        output_path /= "ev-details.csv";;
        cu_details_ev_output = new ofstream(output_path, std::ofstream::out);
        // activate buffers for speedup
        //buffer = new char[bufferSize];
        //cu_details_ev_output->rdbuf()->pubsetbuf(buffer, bufferSize);
        // add header to output file
        *(cu_details_ev_output) << "TimestepID,CarID,EVState,P_charging_kW,cumsum_E_charged_at_home_kWh,cumsum_E_min_kWh,cumsum_E_max_kWh,ev_bs_SOE_kWh" << endl;
    }
}

void output::closeOutputs() {
    // close output file for substations
    if (substation_output != NULL) {
        substation_output->close();
        delete substation_output;
        substation_output = NULL;
    }
    if (substation_output_details != NULL) {
        substation_output_details->close();
        delete substation_output_details;
        substation_output_details = NULL;
    }
    // close surplus output file
    if (surplus_output != NULL) {
        surplus_output->close();
        delete surplus_output;
        surplus_output = NULL;
    }
    //
    // close outputs for CUs if existing
    if (cu_single_output != NULL) {
        cu_single_output->close_buffer();
        delete cu_single_output;
        cu_single_output = NULL;
    }
    if (cu_multi_outputs != NULL) {
        for (size_t i = 0; i < n_cu_multi_outputs; i++) {
            cu_multi_outputs[i]->close_buffer();
            delete cu_multi_outputs[i];
        }
        delete[] cu_multi_outputs;
        cu_multi_outputs = NULL;
    }
    //
    // remove references to output objects from all CUs
    const std::vector<ControlUnit*>& cuList = ControlUnit::GetArrayOfInstances();
    for (ControlUnit* c : cuList) {
        c->set_output_object(NULL);
    }
    //
    if (cu_details_ccmd_output != NULL) {
        // add a closing bracket to the end
        (*cu_details_ccmd_output) << "]";
        // closing part
        cu_details_ccmd_output->close();
        delete cu_details_ccmd_output;
        cu_details_ccmd_output = NULL;
    }
    if (cu_details_ev_output != NULL) {
        cu_details_ev_output->flush();
        delete cu_details_ev_output;
        cu_details_ev_output = NULL;
    }
}

void output::flushBuffers() {
    if (substation_output != NULL) {
        substation_output->flush();
    }
    if (substation_output_details != NULL) {
        substation_output_details->flush();
    }
    if (surplus_output != NULL) {
        surplus_output->flush();
    }
    //
    // flush outputs of CUs if existing
    if (cu_single_output != NULL)
        cu_single_output->flush_buffer();
    if (cu_multi_outputs != NULL)
        for (size_t i = 0; i < n_cu_multi_outputs; i++)
            cu_multi_outputs[i]->flush_buffer();
    //
    if (cu_details_ccmd_output != NULL)
        cu_details_ccmd_output->flush();
    if (cu_details_ev_output != NULL)
        cu_details_ev_output->flush();
}

//
// This function outputs the current parameter variation combination
// to a file in the current directory.
// Changed variables have to be placed in the first argument.
//
void output::outputCurrentParamVariCombi(CurrentParamValues& cParamVals) {
    filesystem::path param_vari_output {*(global::current_output_dir)};
    param_vari_output /= "parameter-settings.csv";
    ofstream ofs(param_vari_output, std::ofstream::out);
    ofs << "Parameter Name,Value\n";
    //
    // If a variable is changed, the boolean indicator variable in the
    // argument cParamVals says so -> in this case one uses this value.
    // Otherwise (if variable is unchanged) use global value.
    ofs << "expansion PV min kWp for section usage,";
    if (cParamVals.exp_pv_min_kWp_roof_sec_set) ofs << cParamVals.exp_pv_min_kWp_roof_sec; else ofs << Global::get_exp_pv_min_kWp_roof_sec();
    ofs << "\n";
    ofs << "expansion PV max inst kWp per section,";
    if (cParamVals.exp_pv_max_kWp_roof_sec_set) ofs << cParamVals.exp_pv_max_kWp_roof_sec; else ofs << Global::get_exp_pv_max_kWp_roof_sec();
    ofs << "\n";
    ofs << "expansion PV max inst kWp per unit,";
    if (cParamVals.exp_pv_max_kWp_per_unit_set) ofs << cParamVals.exp_pv_max_kWp_per_unit; else ofs << Global::get_exp_pv_max_kWp_per_unit();
    ofs << "\n";
    ofs << "expansion PV kWp per roof area in m2,";
    if (cParamVals.exp_pv_kWp_per_m2_set)       ofs << cParamVals.exp_pv_kWp_per_m2;       else ofs << Global::get_exp_pv_kWp_per_m2();
    ofs << "\n";
    ofs << "expansion PV kWp static,";
    if (cParamVals.exp_pv_kWp_static_set)       ofs << cParamVals.exp_pv_kWp_static;       else ofs << Global::get_exp_pv_kWp_static();
    ofs << "\n";
    ofs << "expansion BS E in kWh,";
    if (cParamVals.exp_bs_maxE_kWh_set) ofs << cParamVals.exp_bs_maxE_kWh; else ofs << Global::get_exp_bess_kWh();
    ofs << "\n";
    ofs << "expansion BS E:P ratio,";
    if (cParamVals.exp_bs_EP_ratio_set) ofs << cParamVals.exp_bs_EP_ratio; else ofs << Global::get_exp_bess_E_P_ratio();
    ofs << "\n";
    ofs << "expansion BS initial SOC,";
    if (cParamVals.exp_bs_init_SOC_set) ofs << cParamVals.exp_bs_init_SOC; else ofs << Global::get_exp_bess_start_soc();
    ofs << "\n";
    ofs << "expansion BS efficiency in,";
    if (cParamVals.exp_bs_effi_in_and_out_set) ofs << cParamVals.exp_bs_effi_in_and_out; else ofs << Global::get_exp_bess_effi_in();
    ofs << "\n";
    ofs << "expansion BS efficiency out,";
    if (cParamVals.exp_bs_effi_in_and_out_set) ofs << cParamVals.exp_bs_effi_in_and_out; else ofs << Global::get_exp_bess_effi_out();
    ofs << "\n";
    ofs << "expansion BS self-discharge per ts,";
    if (cParamVals.exp_bs_self_ds_ts_set) ofs << cParamVals.exp_bs_self_ds_ts; else ofs << Global::get_exp_bess_self_ds_ts();
    ofs << "\n";
    ofs << "control horizont in ts,";
    if (cParamVals.control_horizon_in_ts_set) ofs << cParamVals.control_horizon_in_ts; else ofs << Global::get_control_horizon_in_ts();
    ofs << "\n";
    ofs << "control update freq in ts,";
    if (cParamVals.control_update_freq_in_ts_set) ofs << cParamVals.control_update_freq_in_ts; else ofs << Global::get_control_update_freq_in_ts();
    ofs << "\n";
    ofs << "surplus controller freq in ts,";
    if (cParamVals.surplus_controller_freq_in_ts_set) ofs << cParamVals.surplus_controller_freq_in_ts; else ofs << Global::get_surplus_controller_frequency_ts();
    ofs << "\n";
    ofs << "surplus controller lookahead horizon in ts,";
    if (cParamVals.surplus_controller_lookahead_horizon_in_ts_set) ofs << cParamVals.surplus_controller_lookahead_horizon_in_ts; else ofs << Global::get_surplus_controller_lookahead_horizon_ts();
    ofs << "\n";
    ofs << "seed_set," << Global::is_seed_set() << "\n";
    ofs << "seed," << Global::get_seed() << "\n";
    //
    ofs.close();
}

void output::outputCurrentCUSettings() {
    // 1)
    // output information per control unit about sim. added kWp of PV, BS ...
    // not required anymore -> this information is contained in the metrics file now
    const std::vector<ControlUnit*>& cuList = ControlUnit::GetArrayOfInstances();
    filesystem::path output_path {*(global::current_output_dir)};
    /*
    output_path /= "cu-parameters.csv";
    ofstream ofs(output_path, std::ofstream::out);
    ofs << "UnitID,PV kWp,BS P kW,BS E kWh,n_EVs,CS_max_P_kW\n";
    // now, iterate over all control units
    for (size_t i = 0; i < ControlUnit::GetNumberOfInstances(); i++) {
        ControlUnit* cu = cuList[i];
        ofs << cu->get_unitID()            << ",";
        ofs << cu->get_sim_comp_pv_kWp()   << ",";
        ofs << cu->get_sim_comp_bs_P_kW()  << ",";
        ofs << cu->get_sim_comp_bs_E_kWh() << ",";
        ofs << cu->get_sim_comp_cs_n_EVs() << ",";
        ofs << cu->get_sim_comp_cs_max_P_kW() << "\n";
    }
    ofs.close();
    */
    // 2)
    // output information about every sim. added roof section
    //output_path  = *(global::current_output_dir);
    output_path /= "sim-added-roof-sections-per-cu.csv";
    ofstream ofs2(output_path, std::ofstream::out);
    ofs2 << "UnitID,roof_section_number,section_kWp,orientation,profile_index\n";
    for (ControlUnit* c : cuList) {
        string* cu_string = c->get_pv_section_string();
        ofs2 << *cu_string;
        delete cu_string;
    }
    ofs2.close();
}

//
// This function outputs the computed metrics for all control units
// after the simulation has been finished.
// If no metric computetion is selected, the function does nothing
//
void output::outputMetrics(bool alt_fname /* = false */, string * fname_postfix /* = NULL */) {
        filesystem::path output_path;
        if (alt_fname) {
            if (fname_postfix == NULL)
                throw logic_error("Parameter fname_postfix of function output::outputMetrics is NULL!");
            output_path  = *(global::current_global_output_dir);
            output_path /= "metrics-" + *fname_postfix + ".csv";
        } else {
            output_path  = *(global::current_output_dir);
            output_path /= "metrics-per-cu.csv";
        }
        ofstream ofs(output_path, std::ofstream::out);
        ofs << ControlUnit::MetricsStringHeaderAnnual <<"\n";
        //
        // loop over all CUs and get metrics output string
        const std::vector<ControlUnit*>& cuList = ControlUnit::GetArrayOfInstances();
        for (ControlUnit* c : cuList) {
            string* output_str = c->get_metrics_string_annual();
            if (output_str != NULL)
                ofs << *output_str;
            ofs << "\n";
            delete output_str;
        }
        //
        ofs.close();
}


void output::outputMetricsStrListSACPlanning(list<string*> &output_list) {
    filesystem::path output_path;
    output_path  = *(global::current_global_output_dir);
    output_path /= "metrics-sac-planning-per-cu.csv";
    ofstream ofs(output_path, std::ofstream::out);
    ofs << ControlUnit::MetricsStringHeaderAnnual << ",Added components\n";
    //
    // output all elements from the list
    for (string* s : output_list) {
        ofs << *s << "\n";
    }
    //
    ofs.close();
}


void output::outputWeeklyMetricsStrList(list<string*> *output_list, unsigned long week_number) {
    filesystem::path output_path;
    output_path  = *(global::current_output_dir);
    output_path /= "weekly-metrics";
    // create subfolder if it does not exist already
    if (!filesystem::exists(output_path))
        filesystem::create_directory(output_path);
    output_path /= "week-" + std::to_string(week_number) + ".csv";
    ofstream ofs(output_path, std::ofstream::out);
    ofs << ControlUnit::MetricsStringHeaderWeekly << "\n";
    //
    // output all elements from the list
    for (string* s : *output_list) {
        ofs << *s << "\n";
    }
    //
    ofs.close();
    // delete the output list
    for (string* s : *output_list) delete s;
    delete output_list;
}

//
// This function outputs the computed metrics for all EVs
// after the simulation has been finished.
// If no EV has be simulated, the function does nothing
//
void output::outputEVMetrics(bool alt_fname /* = false */, string * fname_postfix /* = NULL */) {
    // return immediatley if there are not EVs
    if (EVFSM::GetNumberOfEVs() == 0)
        return;
    //
    filesystem::path output_path;
    if (alt_fname) {
        if (fname_postfix == NULL)
            throw logic_error("Parameter fname_postfix of function output::outputEVMetrics() is NULL!");
        output_path  = *(global::current_global_output_dir);
        output_path /= "metrics-per-ev-" + *fname_postfix + ".csv";
    } else {
        output_path  = *(global::current_output_dir);
        output_path /= "metrics-per-ev.csv";
    }
    ofstream ofs(output_path, std::ofstream::out);
    ofs << EVFSM::MetricsStringHeaderAnnual <<"\n";
    //
    // loop over all CUs and get metrics output string
    const std::map<unsigned long, EVFSM*>& evList = EVFSM::GetArrayOfInstances();
    for (auto const& [key, ev_ref] : evList) {
        string* output_str = ev_ref->get_metrics_string_annual();
        if (output_str != NULL)
            ofs << *output_str;
        ofs << "\n";
        delete output_str;
    }
    //
    ofs.close();
}

void output::outputRuntimeInformation(long seconds_setup, long seconds_main_run) {
    ofstream time_output( *global::current_output_dir_prefix /= "runtime-information.csv", std::ofstream::out );
    time_output << "key,value\n";
    time_output << "Setup time in s," << seconds_setup << "\n";
    time_output << "Main run time in s," << seconds_main_run << "\n";
    time_output.close();
}

// Template function to serialize a vector
template <typename T>
std::string serializeVector(const std::vector<T>& vec) {
    std::ostringstream voss;
    voss << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) voss << ", ";
        voss << vec[i];
    }
    voss << "]";
    return voss.str();
}

// Template function to serialize a vector of vectors
template <typename T>
std::string serializeVectorOfVectors(const std::vector<const std::vector<T>*>* vec) {
    if (!vec) return "[]";
    std::ostringstream voss;
    voss << "[";
    for (size_t i = 0; i < vec->size(); ++i) {
        if (i > 0) voss << ", ";
        voss << serializeVector(*((*vec)[i]));
    }
    voss << "]";
    return voss.str();
}

// Template function to serialize a vector of vectors
template <typename T>
std::string serializeVectorOfVectorsB(const std::vector<std::vector<T>>& vec) {
    std::ostringstream voss;
    voss << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) voss << ", ";
        voss << serializeVector(vec[i]);
    }
    voss << "]";
    return voss.str();
}

void output::outputControlCommandDetails(
    unsigned long ts, unsigned long cuID, bool optimization_state_ok,
    double inp_max_p_bs_kW,
    double inp_max_e_bs_kWh,
    double inp_max_p_cs_kW,
    double inp_current_bs_charge_kWh,
    const std::vector<float>&  inp_future_resid_demand_kW,
    const std::vector<double>& inp_future_pv_generation_kW,
    const std::vector<double>& inp_future_hp_shiftable_maxP,
    const std::vector<double>& inp_future_hp_shiftable_minP,
    const std::vector<double>& inp_future_hp_shiftable_maxE,
    const std::vector<double>& inp_future_hp_shiftable_minE,
    const std::vector<const std::vector<double>*>* inp_future_ev_shiftable_maxE,
    const std::vector<const std::vector<double>*>* inp_future_ev_shiftable_minE,
    const std::vector<const std::vector<double>*>* inp_future_ev_maxP,
    const std::vector<double>& out_future_bs_power_kW,
    const std::vector<double>& out_future_hp_power_kW,
    const std::vector<std::vector<double>>& out_future_ev_power_kW
) {
    std::unique_lock lock(mtx_cu_details_ccmd); // secure access by using a mutex to prevent multiple parallel function calls

    static std::atomic<bool> first_call = true; // persistent variable to track, if this is the first call or not
    //static bool first_call = true;

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6); // Ensure consistent floating-point precision
    oss << "{";
    oss << "\"ts\": " << ts << ", ";
    oss << "\"cuID\": " << cuID << ", ";
    oss << "\"optimization_state_ok\": " << (optimization_state_ok ? "true" : "false") << ", ";
    oss << "\"inp_max_p_bs_kW\": " << inp_max_p_bs_kW << ", ";
    oss << "\"inp_max_e_bs_kWh\": " << inp_max_e_bs_kWh << ", ";
    oss << "\"inp_max_p_cs_kW\": " << inp_max_p_cs_kW << ", ";
    oss << "\"inp_current_bs_charge_kWh\": " << inp_current_bs_charge_kWh << ", ";

    oss << "\"inp_future_resid_demand_kW\": " << serializeVector(inp_future_resid_demand_kW) << ", ";
    oss << "\"inp_future_pv_generation_kW\": " << serializeVector(inp_future_pv_generation_kW) << ", ";
    oss << "\"inp_future_hp_shiftable_maxP\": " << serializeVector(inp_future_hp_shiftable_maxP) << ", ";
    oss << "\"inp_future_hp_shiftable_minP\": " << serializeVector(inp_future_hp_shiftable_minP) << ", ";
    oss << "\"inp_future_hp_shiftable_maxE\": " << serializeVector(inp_future_hp_shiftable_maxE) << ", ";
    oss << "\"inp_future_hp_shiftable_minE\": " << serializeVector(inp_future_hp_shiftable_minE) << ", ";
    oss << "\"inp_future_ev_shiftable_maxE\": " << serializeVectorOfVectors(inp_future_ev_shiftable_maxE) << ", ";
    oss << "\"inp_future_ev_shiftable_minE\": " << serializeVectorOfVectors(inp_future_ev_shiftable_minE) << ", ";
    oss << "\"inp_future_ev_maxP\": " << serializeVectorOfVectors(inp_future_ev_maxP) << ", ";
    oss << "\"out_future_bs_power_kW\": " << serializeVector(out_future_bs_power_kW) << ", ";
    oss << "\"out_future_hp_power_kW\": " << serializeVector(out_future_hp_power_kW) << ", ";
    oss << "\"out_future_ev_power_kW\": " << serializeVectorOfVectorsB<double>(out_future_ev_power_kW);
    oss << "}";

    if (first_call) {
        first_call = false;
    } else {
        *(cu_details_ccmd_output) << ",";
    }
    *(cu_details_ccmd_output) << oss.str() << "\n";
}

void output::outputEVStateDetails(unsigned long ts, unsigned long carID, EVState ev_state, float p_charging_kW, float cumsum_E_ch_home, float cumsum_E_min, float cumsum_E_max, float ev_bs_SOE_kWh) {
    std::unique_lock lock(mtx_cu_details_ev); // secure access by using a mutex
    std::string ev_state_str;
    if (ev_state == EVState::ConnectedAtHome) {
        ev_state_str = "ConnectedAtHome";
    } else if (ev_state == EVState::DisconnectedAtHome) {
        ev_state_str = "DisconnectedAtHome";
    } else if (ev_state == EVState::Driving) {
        ev_state_str = "Driving";
    } else {
        ev_state_str = "unknown";
    }
    *(cu_details_ev_output) << ts << "," << carID << "," << ev_state_str << "," << p_charging_kW << "," << cumsum_E_ch_home << "," << cumsum_E_min << "," << cumsum_E_max << "," << ev_bs_SOE_kWh << "\n";
}


/////////////////////////////////
//    Implementation of all    //
//       sub-classes of        //
//          CUOutput           //
/////////////////////////////////

CUOutput::~CUOutput() {
    // close buffers, if still opend
    if (buffer_open)
        close_buffer();
    if (output_stream != NULL) {
        delete output_stream;
        output_stream = NULL;
    }
}

void CUOutput::close_buffer() {
    if (buffer_open)
        output_stream->close();
}

CUOutputSingleFile::CUOutputSingleFile(unsigned long scenario_id) {
    //
    // initialize the output file
    stringstream output_path_CUs;
    output_path_CUs << setw(4) << setfill('0') << scenario_id;
    output_path_CUs << "-CU-time-series.csv";
    filesystem::path output_path = *(global::current_output_dir);
    output_path /= output_path_CUs.str();
    output_stream = new ofstream(output_path, std::ofstream::out);
    //
    // activate buffers for speedup
    buffer = new char[bufferSize];
    output_stream->rdbuf()->pubsetbuf(buffer, bufferSize);
    //
    // add header to output file
    *(output_stream) << "Timestep,ControlUnitID,Load_vSmartMeter_kW,Load_rSmartMeters_kW,Load_self_produced_kW,PVFeedin_simulated_kW,BS_SOC,BS_load_kW,HP_load_kW,CS_load_kW,CS_n_EVs_conn,CS_n_EVs_not_conn" << endl;
    buffer_open = true;
}

CUOutputSingleFile::~CUOutputSingleFile() {
    delete[] buffer;
    buffer = NULL;
}

CUOutputOneFilePerCU::CUOutputOneFilePerCU(int cuID, filesystem::path& dirpath) {
    stringstream filename_i;
    filename_i << setw(5) << setfill('0') << cuID << "-CU-ts.csv";
    filesystem::path filepath_i;
    filepath_i /= dirpath;
    filepath_i /= filename_i.str();
    output_stream = new ofstream(filepath_i, std::ofstream::out);
    *(output_stream) << "Timestep,ControlUnitID,Load_vSmartMeter_kW,Load_rSmartMeters_kW,Load_self_produced_kW,PVFeedin_simulated_kW,BS_SOC,BS_load_kW,HP_load_kW,CS_load_kW,CS_n_EVs_conn,CS_n_EVs_not_conn\n";
    buffer_open = true;
}

CUOutputOneFilePerSubstation::CUOutputOneFilePerSubstation(const string* substName, filesystem::path& dirpath) {
    stringstream filename_i;
    filename_i << *(substName) << "-AllCUs-ts.csv";
    filesystem::path filepath_i;
    filepath_i /= dirpath;
    filepath_i /= filename_i.str();
    output_stream = new ofstream(filepath_i, std::ofstream::out);
    *(output_stream) << "Timestep,ControlUnitID,Load_vSmartMeter_kW,Load_rSmartMeters_kW,Load_self_produced_kW,PVFeedin_simulated_kW,BS_SOC,BS_load_kW,HP_load_kW,CS_load_kW,CS_n_EVs_conn,CS_n_EVs_not_conn\n";
    buffer_open = true;
}

void CUOutputSingleFile::output_for_one_cu(
        size_t cuID,     size_t ts,            double load_vsm,
        double load_rsm, double load_selfprod, double load_pv,
        double bs_SOC,   double load_bs,       double load_hp,
        double load_cs,  size_t n_cars_pc,     size_t n_cars_pnc)
{
    unique_lock lock(single_file_mutex); // secure access by using a mutex
    *(output_stream) << ts << "," << cuID << ","
         << round_float_5(load_vsm) << ","
         << round_float_5(load_rsm) << ","
         << round_float_5(load_selfprod) << ","
         << round_float_5(load_pv) << ","
         << round_float_5(bs_SOC) << ","
         << round_float_5(load_bs) << ","
         << round_float_5(load_hp) << ","
         << round_float_5(load_cs) << ","
         << n_cars_pc << "," << n_cars_pnc << "\n";
}

void CUOutputOneFilePerCU::output_for_one_cu(
        size_t cuID,     size_t ts,            double load_vsm,
        double load_rsm, double load_selfprod, double load_pv,
        double bs_SOC,   double load_bs,       double load_hp,
        double load_cs,  size_t n_cars_pc,     size_t n_cars_pnc)
{
    *(output_stream) << ts << "," << cuID << ","
        << round_float_5(load_vsm) << ","
        << round_float_5(load_rsm) << ","
        << round_float_5(load_selfprod) << ","
        << round_float_5(load_pv) << ","
        << round_float_5(bs_SOC) << ","
        << round_float_5(load_bs) << ","
        << round_float_5(load_hp) << ","
        << round_float_5(load_cs) << "," << n_cars_pc << "," << n_cars_pnc << "\n";
}

void CUOutputOneFilePerSubstation::output_for_one_cu(
        size_t cuID,     size_t ts,            double load_vsm,
        double load_rsm, double load_selfprod, double load_pv,
        double bs_SOC,   double load_bs,       double load_hp,
        double load_cs,  size_t n_cars_pc,     size_t n_cars_pnc)
{
    unique_lock lock(single_file_mutex); // secure access by using a mutex
    *(output_stream) << ts << "," << cuID << ","
        << round_float_5(load_vsm) << ","
        << round_float_5(load_rsm) << ","
        << round_float_5(load_selfprod) << ","
        << round_float_5(load_pv) << ","
        << round_float_5(bs_SOC) << ","
        << round_float_5(load_bs) << ","
        << round_float_5(load_hp) << ","
        << round_float_5(load_cs) << "," << n_cars_pc << "," << n_cars_pnc << "\n";
}

void CUOutputSingleFile::flush_buffer() {
    unique_lock lock(single_file_mutex); // secure access by using a mutex
    output_stream->flush();
}

void CUOutputOneFilePerCU::flush_buffer() {
    output_stream->flush();
}

void CUOutputOneFilePerSubstation::flush_buffer() {
    unique_lock lock(single_file_mutex); // secure access by using a mutex
    output_stream->flush();
}

