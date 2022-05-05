#include "output.h"

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
}

void output::initializeDirectoriesPerPVar(int scenario_id) {
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
    const int nSubst = Substation::GetNumberOfInstances();
    for (int i = 0; i < nSubst; i++) {
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
        for (int i = 0; i < ControlUnit::GetNumberOfInstances(); i++) {
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
        for (int nSubst = 1; nSubst <= Global::get_n_substations(); nSubst++) {
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
        for (int i = 0; i < Global::get_n_CUs(); i++)
            cu_multi_outputs[i]->flush_buffer();
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
    *(output_stream) << "Timestep,ControlUnitID,Load_vSmartMeter_kW,Load_rSmartMeters_kW,Load_self_produced_kW,PVFeedin_simulated_kW,BS_SOC,BS_load_kW\n";
    buffer_open = true;
}

CUOutputOneFilePerSubstation::CUOutputOneFilePerSubstation(const string* substName, filesystem::path& dirpath) {
    stringstream filename_i;
    filename_i << *(substName) << "-AllCUs-ts.csv";
    filesystem::path filepath_i;
    filepath_i /= dirpath;
    filepath_i /= filename_i.str();
    output_stream = new ofstream(filepath_i, std::ofstream::out);
    *(output_stream) << "Timestep,ControlUnitID,Load_vSmartMeter_kW,Load_rSmartMeters_kW,Load_self_produced_kW,PVFeedin_simulated_kW,BS_SOC,BS_load_kW\n";
    buffer_open = true;
}

void CUOutputSingleFile::output_for_one_cu(
        int   cuID,     int   ts,            float load_vsm,
        float load_rsm, float load_selfprod, float load_pv,
        float bs_SOC,   float load_bs)
{
    unique_lock lock(single_file_mutex); // secure access by using a mutex
    *(output_stream) << ts << "," << cuID << "," << load_vsm << "," << load_rsm << "," << load_selfprod << "," << load_pv << "," << bs_SOC << "," << load_bs << "\n";
}

void CUOutputOneFilePerCU::output_for_one_cu(
        int   cuID,     int   ts,            float load_vsm,
        float load_rsm, float load_selfprod, float load_pv,
        float bs_SOC,   float load_bs)
{
    *(output_stream) << ts << "," << cuID << "," << load_vsm << "," << load_rsm << "," << load_selfprod << "," << load_pv << "," << bs_SOC << "," << load_bs << "\n";
}

void CUOutputOneFilePerSubstation::output_for_one_cu(
        int   cuID,     int   ts,            float load_vsm,
        float load_rsm, float load_selfprod, float load_pv,
        float bs_SOC,   float load_bs)
{
    *(output_stream) << ts << "," << cuID << "," << load_vsm << "," << load_rsm << "," << load_selfprod << "," << load_pv << "," << bs_SOC << "," << load_bs << "\n";
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

