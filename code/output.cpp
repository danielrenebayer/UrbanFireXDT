#include "output.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "global.h"

using namespace std;
using namespace output;

void output::initializeSubstationOutput(int scenario_id) {
    //
    // This method initializes the substation output
    // file.
    //
    //
    // initialize the output file
    stringstream output_path_subst;
    output_path_subst << Global::get_output_path();
    output_path_subst << setw(4) << setfill('0') << scenario_id;
    output_path_subst << "-substation-time-series.csv";
    substation_output = new ofstream(output_path_subst.str().c_str(), std::ofstream::out);
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
        // Case 2: One file for each individual CU
        //
        // create output directory, and delete existing if present
        filesystem::path dirpath = Global::get_output_path();
        dirpath /= "ts-per-cu";
        if (filesystem::is_directory(dirpath)) {
            // clear existing directory
            filesystem::remove_all(dirpath);
        }
        // now, actually create output dir (maybe again)
        filesystem::create_directory(dirpath);
        //
        // create a output object for every (future) control unit
        // and add reference to all CUs
        cu_multi_outputs = new CUOutputOneFilePerCU*[Global::get_n_CUs()];
        ControlUnit*const* cuList = ControlUnit::GetArrayOfInstances();
        for (int cuID = 1; cuID <= Global::get_n_CUs(); cuID++) {
            cu_multi_outputs[cuID-1] = new CUOutputOneFilePerCU(cuID, dirpath);
            cuList[cuID-1]->set_output_object(cu_multi_outputs[cuID-1]);
        }
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
        for (int i = 0; i < Global::get_n_CUs(); i++) {
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
    output_stream->close();
}

CUOutputSingleFile::CUOutputSingleFile(int scenario_id) {
    //
    // initialize the output file
    stringstream output_path_CUs;
    output_path_CUs << Global::get_output_path();
    output_path_CUs << setw(4) << setfill('0') << scenario_id;
    output_path_CUs << "-CU-time-series.csv";
    output_stream = new ofstream(output_path_CUs.str().c_str(), std::ofstream::out);
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

void CUOutputSingleFile::flush_buffer() {
    unique_lock lock(single_file_mutex); // secure access by using a mutex
    output_stream->flush();
}

void CUOutputOneFilePerCU::flush_buffer() {
    output_stream->flush();
}

