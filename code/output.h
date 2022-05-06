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

    void initializeDirectoriesOnce(int scenario_id);
    void initializeDirectoriesPerPVar(int scenario_id);
    void initializeSubstationOutput(int scenario_id);
    void initializeCUOutput(int scenario_id);

    void closeOutputs();

    void flushBuffers();

    // struct for the below defined funtion:
    // it holds the current parameter values
    struct CurrentParamValues {
        float exp_pv_kWp      = 0.0;  bool exp_pv_kWp_set      = false;
        float exp_bs_maxP_kW  = 0.0;  bool exp_bs_maxP_kW_set  = false;
        float exp_bs_maxE_kWh = 0.0;  bool exp_bs_maxE_kWh_set = false;
        float exp_bs_init_SOC = 0.0;  bool exp_bs_init_SOC_set = false;
    };
    void outputCurrentParamVariCombi(CurrentParamValues&);

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
                int cuID,            int ts,
                float load_vsm,      float load_rsm,
                float load_selfprod, float load_pv,
                float bs_SOC,        float load_bs) = 0;
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
                int cuID,            int ts,
                float load_vsm,      float load_rsm,
                float load_selfprod, float load_pv,
                float bs_SOC,        float load_bs);
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
                int cuID,            int ts,
                float load_vsm,      float load_rsm,
                float load_selfprod, float load_pv,
                float bs_SOC,        float load_bs);
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
                int cuID,            int ts,
                float load_vsm,      float load_rsm,
                float load_selfprod, float load_pv,
                float bs_SOC,        float load_bs);
        void flush_buffer();
};

#endif

