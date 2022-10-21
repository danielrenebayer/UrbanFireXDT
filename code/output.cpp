#include "output.h"

#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "global.h"

using namespace std;
using namespace output;

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

void output::initializeDirectoriesOnce(int scenario_id) {
    //
    // This function initializes the direcory (or directories)
    // for the output, for the current scenario.
    // Thus, it has to be callen for every parameter variation
    // setting individually again, if parameter variation
    // is selected at all.
    //
    // MUST be callen BEFORE initializeDirectoriesPerPVar()
    //
    //
    // create a subfolder for the current scenario, if it does not exist
    filesystem::path dirpath = Global::get_output_path();
    if (!filesystem::is_directory(dirpath)) {
        filesystem::create_directory(dirpath);
    }
    stringstream current_scenario_str;
    current_scenario_str << "S";
    current_scenario_str << setw(4) << setfill('0') << scenario_id;
    dirpath /= current_scenario_str.str();
    if (!filesystem::is_directory(dirpath)) {
        filesystem::create_directory(dirpath);
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
}

void output::initializeDirectoriesPerPVar() {
    //
    // This function initializes the direcory (or directories)
    // for the current parameter variation (if selected).
    // Thus, it has to be callen for every parameter variation
    // setting individually again.
    // Even in the case of no parameter variation, it has to be callen once.
    //
    // MUST be callen AFTER initializeDirectoriesOnce()
    //
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

void output::initializeSubstationOutput(int scenario_id) {
    //
    // This method initializes the substation output
    // file.
    //
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
    Substation*const* subList = Substation::GetArrayOfInstances();
    const size_t nSubst = Substation::GetNumberOfInstances();
    for (size_t i = 0; i < nSubst; i++) {
        *(substation_output) << "," << subList[i]->get_name()->c_str();
    }
    *(substation_output) << ",open_space_pv_feedin,wind_feedin,total_load" << endl;
    substation_output_init = true;
}

void output::initializeCUOutput(int scenario_id) {
    //
    // initializes the CU output, depending on the global setting
    if (Global::get_output_mode_per_cu() == global::OutputModePerCU::SingleFile) {
        // Case 1: One file for all CUs
        //
        cu_single_output = new CUOutputSingleFile(scenario_id);
        //
        // add reference to all CUs
        ControlUnit*const* cuList = ControlUnit::GetArrayOfInstances();
        for (size_t i = 0; i < ControlUnit::GetNumberOfInstances(); i++) {
            cuList[i]->set_output_object(cu_single_output);
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
        Substation*const* substationList = Substation::GetArrayOfInstances();
        for (size_t nSubst = 1; nSubst <= Global::get_n_substations(); nSubst++) {
            Substation* currentS = substationList[nSubst-1];
            cu_multi_outputs[nSubst-1] = new CUOutputOneFilePerSubstation(currentS->get_name(), dirpath);
            // add output file to all connected CUs
            const list<ControlUnit*>* conn_units = currentS->get_connected_units();
            for (auto cu : *conn_units) {
                cu->set_output_object(cu_multi_outputs[nSubst-1]);
            }
        }
        /*ControlUnit*const* cuList = ControlUnit::GetArrayOfInstances();
        for (int cuID = 1; cuID <= Global::get_n_CUs(); cuID++) {
            cu_multi_outputs[cuID-1] = new CUOutputOneFilePerCU(cuID, dirpath);
            cuList[cuID-1]->set_output_object(cu_multi_outputs[cuID-1]);
        }*/
    } else {
        // Case 3: no output per CU selected
    }
}

void output::closeOutputs() {
    // close output file for substations
    if (substation_output_init) {
        substation_output->close();
        delete substation_output;
        substation_output = NULL;
        substation_output_init = false;
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
    ControlUnit*const* cuList = ControlUnit::GetArrayOfInstances();
    for (size_t i = 0; i < ControlUnit::GetNumberOfInstances(); i++) {
        cuList[i]->set_output_object(NULL);
    }
}

void output::flushBuffers() {
    if (substation_output_init) {
        substation_output->flush();
    }
    //
    // flush outputs of CUs if existing
    if (cu_single_output != NULL)
        cu_single_output->flush_buffer();
    if (cu_multi_outputs != NULL)
        for (size_t i = 0; i < n_cu_multi_outputs; i++)
            cu_multi_outputs[i]->flush_buffer();
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
    ofs << "expansion PV kWp per roof area in m2,";
    if (cParamVals.exp_pv_kWp_per_m2_set)       ofs << cParamVals.exp_pv_kWp_per_m2;       else ofs << Global::get_exp_pv_kWp_per_m2();
    ofs << "\n";
    ofs << "expansion PV kWp static,";
    if (cParamVals.exp_pv_kWp_static_set)       ofs << cParamVals.exp_pv_kWp_static;       else ofs << Global::get_exp_pv_kWp_static();
    ofs << "\n";
    ofs << "expansion PV kWp static mode,";
    if (cParamVals.exp_pv_static_mode_set)      ofs << cParamVals.exp_pv_static_mode;      else ofs << Global::get_exp_pv_static_mode();
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
    //
    ofs.close();
}

void output::outputCurrentCUSettings() {
    // 1)
    // output information per control unit about sim. added kWp of PV, BS ...
    filesystem::path output_path {*(global::current_output_dir)};
    output_path /= "cu-parameters.csv";
    ofstream ofs(output_path, std::ofstream::out);
    ofs << "UnitID,PV kWp,BS P kW,BS E kWh\n";
    // now, iterate over all control units
    ControlUnit*const* cuList = ControlUnit::GetArrayOfInstances();
    for (size_t i = 0; i < ControlUnit::GetNumberOfInstances(); i++) {
        ControlUnit* cu = cuList[i];
        ofs << cu->get_unitID()            << ",";
        ofs << cu->get_sim_comp_pv_kWp()   << ",";
        ofs << cu->get_sim_comp_bs_P_kW()  << ",";
        ofs << cu->get_sim_comp_bs_E_kWh() << "\n";
    }
    ofs.close();
    // 2)
    // output information about every sim. added roof section
    output_path  = *(global::current_output_dir);
    output_path /= "sim-added-roof-sections-per-cu.csv";
    ofstream ofs2(output_path, std::ofstream::out);
    ofs2 << "UnitID,roof_section_number,section_kWp,orientation,profile_index\n";
    for (size_t i = 0; i < ControlUnit::GetNumberOfInstances(); i++) {
        string* cu_string = cuList[i]->get_pv_section_string();
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
    //if (Global::get_comp_eval_metrics()) { // option is disabled, check not required anymore
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
        ofs << "UnitID,SCR,SSR,NPV,Sum of demand [kWh],Sum of self-consumed e. [kWh],Sum of PV-generated e. [kWh],Sum of grid feed-in [kWh],BS EFC\n";
        //
        // loop over all CUs and get metrics output string
        ControlUnit*const* cuList = ControlUnit::GetArrayOfInstances();
        for (size_t i = 0; i < ControlUnit::GetNumberOfInstances(); i++) {
            string* output_str = cuList[i]->get_metrics_string();
            if (output_str != NULL)
                ofs << *output_str;
            ofs << "\n";
            delete output_str;
        }
        //
        ofs.close();
    //}
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

CUOutputSingleFile::CUOutputSingleFile(int scenario_id) {
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
    *(output_stream) << "Timestep,ControlUnitID,Load_vSmartMeter_kW,Load_rSmartMeters_kW,Load_self_produced_kW,PVFeedin_simulated_kW,BS_SOC,BS_load_kW" << endl;
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
    *(output_stream) << "Timestep,ControlUnitID,Load_vSmartMeter_kW,Load_rSmartMeters_kW,Load_self_produced_kW,PVFeedin_simulated_kW,BS_SOC,BS_load_kW,HP_load_kW,WB_load_kW\n";
    buffer_open = true;
}

CUOutputOneFilePerSubstation::CUOutputOneFilePerSubstation(const string* substName, filesystem::path& dirpath) {
    stringstream filename_i;
    filename_i << *(substName) << "-AllCUs-ts.csv";
    filesystem::path filepath_i;
    filepath_i /= dirpath;
    filepath_i /= filename_i.str();
    output_stream = new ofstream(filepath_i, std::ofstream::out);
    *(output_stream) << "Timestep,ControlUnitID,Load_vSmartMeter_kW,Load_rSmartMeters_kW,Load_self_produced_kW,PVFeedin_simulated_kW,BS_SOC,BS_load_kW,HP_load_kW,WB_load_kW\n";
    buffer_open = true;
}

void CUOutputSingleFile::output_for_one_cu(
        size_t cuID,    size_t ts,           float load_vsm,
        float load_rsm, float load_selfprod, float load_pv,
        float bs_SOC,   float load_bs,       float load_hp,
        float load_wb)
{
    unique_lock lock(single_file_mutex); // secure access by using a mutex
    *(output_stream) << ts << "," << cuID << "," << load_vsm << "," << load_rsm << "," << load_selfprod << "," << load_pv << "," << bs_SOC << "," << load_bs << "," << load_hp << "," << load_wb << "\n";
}

void CUOutputOneFilePerCU::output_for_one_cu(
        size_t cuID,    size_t ts,           float load_vsm,
        float load_rsm, float load_selfprod, float load_pv,
        float bs_SOC,   float load_bs,       float load_hp,
        float load_wb)
{
    *(output_stream) << ts << "," << cuID << "," << load_vsm << "," << load_rsm << "," << load_selfprod << "," << load_pv << "," << bs_SOC << "," << load_bs << "," << load_hp << "," << load_wb << "\n";
}

void CUOutputOneFilePerSubstation::output_for_one_cu(
        size_t cuID,    size_t ts,           float load_vsm,
        float load_rsm, float load_selfprod, float load_pv,
        float bs_SOC,   float load_bs,       float load_hp,
        float load_wb)
{
    *(output_stream) << ts << "," << cuID << "," << load_vsm << "," << load_rsm << "," << load_selfprod << "," << load_pv << "," << bs_SOC << "," << load_bs << "," << load_hp << "," << load_wb << "\n";
}

void CUOutputSingleFile::flush_buffer() {
    unique_lock lock(single_file_mutex); // secure access by using a mutex
    output_stream->flush();
}

void CUOutputOneFilePerCU::flush_buffer() {
    output_stream->flush();
}

void CUOutputOneFilePerSubstation::flush_buffer() {
    output_stream->flush();
}

