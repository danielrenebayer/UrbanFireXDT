/*
 * output.h
 *
 * This contains all functions and classes for buffering outputs
 * and writing them to the disk finally.
 *
 */

#ifndef __OUTPUT_H
#define __OUTPUT_H

#include <filesystem>
#include <fstream>
#include <mutex>
#include <sstream>
#include <string>

using namespace std;

class CUOutput;
class CUOutputSingleFile;
class CUOutputOneFilePerCU;
class CUOutputOneFilePerSubstation;

namespace output {

    inline std::ofstream* substation_output;
    inline bool substation_output_init = false;
    inline CUOutputSingleFile* cu_single_output = NULL; ///< Reference to the single_output object, if one output for all CUs is selected
    inline CUOutputOneFilePerSubstation** cu_multi_outputs = NULL; ///< Reference to the array of CU ouputs, if one output per CU is selected
    inline size_t n_cu_multi_outputs = 0; ///< Number of elements in cu_multi_ouputs

    /**
     * This function initializes the base direcory (or directories) for the output,
     * for the current scenario but not for individual parameter variations.
     * It has to be callen once (or multiple times if repetition is selected as
     * cmd-line argument).
     * It must be callen before initializeDirectoriesPerPVar().
     * @param scenario_id The current scenario ID
     */
    void initializeDirectoriesBase(int scenario_id);
    /**
     * This function initializes the direcory (or directories)
     * for the current parameter variation (if selected).
     * Thus, it has to be callen for every parameter variation
     * setting individually again.
     * Even in the case of no parameter variation, it has to be callen once.
     * It MUST be callen AFTER initializeDirectoriesBase().
     */
    void initializeDirectoriesPerPVar();
    /**
     * This method initializes the substation output file.
     * @param scenario_id The current scenario ID
     */
    void initializeSubstationOutput(int scenario_id);
    /**
     * This method initializes the individual output files for the control
     * units. Depending on the globally selected mode, one file per CU is
     * created or some CUs share one output file.
     * @param scenario_id The current scenario ID
     */
    void initializeCUOutput(int scenario_id);

    void closeOutputs();

    void flushBuffers();

    // struct for the below defined funtion:
    // it holds the current parameter values
    struct CurrentParamValues {
        bool  exp_pv_static_mode    = false; bool exp_pv_static_mode_set      = false;
        float exp_pv_kWp_static       = 0.0; bool exp_pv_kWp_static_set       = false;
        float exp_pv_kWp_per_m2       = 0.0; bool exp_pv_kWp_per_m2_set       = false;
        float exp_pv_min_kWp_roof_sec = 0.0; bool exp_pv_min_kWp_roof_sec_set = false;
        float exp_pv_max_kWp_roof_sec = 0.0; bool exp_pv_max_kWp_roof_sec_set = false;
        float exp_bs_maxP_kW  = 0.0;  bool exp_bs_maxP_kW_set  = false;
        float exp_bs_maxE_kWh = 0.0;  bool exp_bs_maxE_kWh_set = false;
        float exp_bs_init_SOC = 0.0;  bool exp_bs_init_SOC_set = false;
        float exp_bs_EP_ratio = 0.0;  bool exp_bs_EP_ratio_set = false;
    };
    void outputCurrentParamVariCombi(CurrentParamValues&);

    void outputCurrentCUSettings(); ///< This function outputs current settings of the control units (like PV kWp, BS capacity and power, ...) and also which roof sections exist per sim. added PV component

    void outputMetrics(bool alt_fname = false, string * fname_postfix = NULL); ///< This function computed metrics for all control units after the simulation has been finished; if @param alt_fname is set to true, the output file will have the file name 'metrics-{fname_postfix}.csv' instead of 'metrics-per-cu.csv'

}

class CUOutput {
    /*
     * This (virtual) class represents the output for one (or more)
     * control units.
     * This is an abstract class, as most methods are 'implemented'
     * as pure virtual functions here.
     */
    public:
        virtual ~CUOutput();
        virtual void output_for_one_cu(
                size_t cuID,         size_t ts,
                float load_vsm,      float load_rsm,
                float load_selfprod, float load_pv,
                float bs_SOC,        float load_bs,
                float load_hp,       float load_wb) = 0;
        virtual void flush_buffer() = 0;
        void close_buffer();
    protected:
        bool            buffer_open;   ///< True, if buffer(s) is/are opened
        std::ofstream*  output_stream; ///< output stream
};

class CUOutputSingleFile : public CUOutput {
    /*
     * This class represents the output for the
     * control units, that is directed into one
     * output file.
     */
    public:
        CUOutputSingleFile(int scenario_id);
        ~CUOutputSingleFile();
        static const size_t bufferSize = 128*1024;
        //
        // definition of virtual methods from base class
        void output_for_one_cu(
                size_t cuID,         size_t ts,
                float load_vsm,      float load_rsm,
                float load_selfprod, float load_pv,
                float bs_SOC,        float load_bs,
                float load_hp,       float load_wb);
        void flush_buffer();
    private:
        mutex  single_file_mutex; ///< If a single file is selected as CU output, this variable holds the mutex to ensure correct concurrency behavior
        char*  buffer; // buffer for speedup
};

class CUOutputOneFilePerCU : public CUOutput {
    /*
     * This class represents the output for one
     * individual control unit.
     */
    public:
        CUOutputOneFilePerCU(int cuID, filesystem::path& dirpath);
        //
        // definition of virtual methods from base class
        void output_for_one_cu(
                size_t cuID,         size_t ts,
                float load_vsm,      float load_rsm,
                float load_selfprod, float load_pv,
                float bs_SOC,        float load_bs,
                float load_hp,       float load_wb);
        void flush_buffer();
};

class CUOutputOneFilePerSubstation : public CUOutput {
    /*
     * This class represents the output for all
     * control units, that are connected to one
     * substation.
     */
    public:
        CUOutputOneFilePerSubstation(const string* substName, filesystem::path& dirpath);
        //
        // definition of virtual methods from base class
        void output_for_one_cu(
                size_t cuID,         size_t ts,
                float load_vsm,      float load_rsm,
                float load_selfprod, float load_pv,
                float bs_SOC,        float load_bs,
                float load_hp,       float load_wb);
        void flush_buffer();
};

#endif

