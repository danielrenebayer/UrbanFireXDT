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
#include <list>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

// The following classes are defined in this header file:
class CUOutput;
class CUOutputSingleFile;
class CUOutputOneFilePerCU;
class CUOutputOneFilePerSubstation;

#include "vehicles.h"

using namespace std;

namespace output {

    inline std::ofstream* substation_output         = NULL; ///< The main file for the substation load time series
    inline std::ofstream* substation_output_details = NULL; ///< The secondary file for additional information about the substations
    inline std::ofstream* cu_details_ccmd_output    = NULL; ///< Output file for command details per time step and control unit (if selected by --ccmd-output option)
    inline std::ofstream* cu_details_ev_output      = NULL; ///< Output file for ev      details per time step and control unit (if selected by --ev-output option)
    inline CUOutputSingleFile* cu_single_output = NULL; ///< Reference to the single_output object, if one output for all CUs is selected
    inline CUOutputOneFilePerSubstation** cu_multi_outputs = NULL; ///< Reference to the array of CU ouputs, if one output per CU is selected
    inline size_t n_cu_multi_outputs = 0; ///< Number of elements in cu_multi_ouputs

    inline std::mutex mtx_cu_details_ccmd; ///< Mutex to ensure proper working in parallel processing for output::cu_details_ccmd_output
    inline std::mutex mtx_cu_details_ev;   ///< Mutex to ensure proper working in parallel processing for output::cu_details_ev_output

    /**
     * This function initializes the base direcory (or directories) for the output,
     * for the current scenario but not for individual parameter variations.
     * It has to be called once (or multiple times if repetition is selected as
     * cmd-line argument).
     * It must be called before initializeDirectoriesPerPVar().
     * @param scenario_id The current scenario ID
     */
    void initializeDirectoriesBase(unsigned long scenario_id);
    /**
     * This function initializes the direcory (or directories)
     * for the current parameter variation (if selected).
     * Thus, it has to be called for every parameter variation
     * setting individually again.
     * Even in the case of no parameter variation, it has to be called once.
     * It MUST be called AFTER initializeDirectoriesBase().
     */
    void initializeDirectoriesPerPVar();
    /**
     * This method initializes the substation output file.
     * @param scenario_id The current scenario ID
     */
    void initializeSubstationOutput(unsigned long scenario_id);
    /**
     * This method initializes the individual output files for the control
     * units. Depending on the globally selected mode, one file per CU is
     * created or some CUs share one output file.
     * Moreover, it initializes the output files for the inidivudal control
     * commands and the EV states per time step if selected by the command
     * line arguments --ccmd-output and --ev-output.
     * @param scenario_id The current scenario ID
     */
    void initializeCUOutput(unsigned long scenario_id);
    /**
     * This function closes all outputs that are opend by initializeCUOutput() or initializeSubstationOutput().
     * It further informs all instances of the class ControlUnit that the output object is deleted.
     */
    void closeOutputs();
    /**
     * This function flushes all buffers that are currently opend and were opened by initializeCUOutput() or initializeSubstationOutput().
     */
    void flushBuffers();

    /**
     * The struct CurrentParamValues is used for storing the current parameter values.
     * As they might change during a parameter variation, it is good to keep
     * track of the currently processe values.
     * 
     * This struct is important especially for outputting the current parameter setting.
     */ 
    struct CurrentParamValues {
        float exp_pv_kWp_static       = 0.0; bool exp_pv_kWp_static_set       = false;
        float exp_pv_kWp_per_m2       = 0.0; bool exp_pv_kWp_per_m2_set       = false;
        float exp_pv_min_kWp_roof_sec = 0.0; bool exp_pv_min_kWp_roof_sec_set = false;
        float exp_pv_max_kWp_roof_sec = 0.0; bool exp_pv_max_kWp_roof_sec_set = false;
        float exp_pv_max_kWp_per_unit = 0.0; bool exp_pv_max_kWp_per_unit_set = false;
        float exp_bs_maxP_kW  = 0.0;  bool exp_bs_maxP_kW_set  = false;
        float exp_bs_maxE_kWh = 0.0;  bool exp_bs_maxE_kWh_set = false;
        float exp_bs_init_SOC = 0.0;  bool exp_bs_init_SOC_set = false;
        float exp_bs_EP_ratio = 0.0;  bool exp_bs_EP_ratio_set = false;
        uint  control_horizon_in_ts     = 24; bool control_horizon_in_ts_set     = false;
        uint  control_update_freq_in_ts =  1; bool control_update_freq_in_ts_set = false;
    };
    void outputCurrentParamVariCombi(CurrentParamValues&);

    void outputCurrentCUSettings(); ///< This function outputs current settings of the control units (like PV kWp, BS capacity and power, ...) and also which roof sections exist per sim. added PV component

    void outputMetrics(bool alt_fname = false, string * fname_postfix = NULL); ///< This function computed metrics for all control units after the simulation has been finished; if @param alt_fname is set to true, the output file will have the file name 'metrics-{fname_postfix}.csv' instead of 'metrics-per-cu.csv'

    /**
     * Outputs every line that is collected in 'output_list' to 'metrics-sac-planning-per-cu.csv'.
     * The output holds an additional last column about the added components for computing the metrics.
     * 
     * @param output_list: Reference to the list containing each line as string (lines must not end with a new line)
     */
    void outputMetricsStrListSACPlanning(list<string*> &output_list);

    /**
     * Outputs every line that is collected in 'output_list' to 'weekly-metrics-per-cu.csv'.
     * Attention: This methods deletes the output list!
     * 
     * @param output_list: Reference to the list containing each line as string (lines must not end with a new line)
     * @param week_number: The week number of the output
     */
    void outputWeeklyMetricsStrList(list<string*> *output_list, unsigned long week_number);

    /**
     * Outputs the computed metrics for all simulated EVs after the simulation has been finished,
     * similar to output::outputMetrics().
     * Creates an output file only, if there is at least one EV.
     * 
     * @param alt_fname if true, the output file will have the file name 'metrics-per-ev-{fname_postfix}.csv'
     * instead of 'metrics-per-ev.csv'.
     * @param fname_postfix: Alternative ending of the file name if required.
     */
    void outputEVMetrics(bool alt_fname = false, string * fname_postfix = NULL);

    /**
     * Output information on the run time to a file.
     * Thus function must not be called before initializeDirectoriesBase().
     * 
     * @param seconds_setup: The duration of the setup and data loading in seconds
     * @param seconds_main_run: The duration of the main run in seconds
     */
    void outputRuntimeInformation(long seconds_setup, long seconds_main_run);

    /**
     * Output the optimization input and output for one time step and control unit.
     * Call this method only if Global::get_create_control_cmd_output() == true.
     * This function is thread safe.
     * All parameters starting with "inp_" are the input parameters for the optimization as defined in BaseOptimizedController::updateController().
     * All parameters starting with "out_" denote the output parameters returned by the controller / optimization.
     */
     void outputControlCommandDetails(
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
    );

    /**
     * Outputs information about one EV for a given time step.
     * Call this function only if Global::get_create_ev_detailed_output() == true.
     * This function is thread safe.
     */
    void outputEVStateDetails(unsigned long ts, unsigned long carID, EVState ev_state, float p_charging_kW, float cumsum_E_ch_home, float cumsum_E_min, float cumsum_E_max, float ev_bs_SOE_kWh);

}

/**
 * This (virtual) class represents the output for one (or more)
 * control units.
 * This is an abstract class, as most methods are 'implemented'
 * as pure virtual functions here.
 **/
class CUOutput {
    public:
        virtual ~CUOutput();
        virtual void output_for_one_cu(
                size_t cuID,          size_t ts,
                double load_vsm,      double load_rsm,
                double load_selfprod, double load_pv,
                double bs_SOC,        double load_bs,
                double load_hp,       double load_cs,
                size_t n_cars_pc,     size_t n_cars_pnc) = 0;
        virtual void flush_buffer() = 0;
        void close_buffer();
    protected:
        bool            buffer_open;   ///< True, if buffer(s) is/are opened
        std::ofstream*  output_stream; ///< output stream
};

/**
 * This class represents the output for the
 * control units, that is directed into one
 * output file.
 * 
 * This class is thread-safe.
 **/
class CUOutputSingleFile : public CUOutput {
    public:
        CUOutputSingleFile(unsigned long scenario_id);
        ~CUOutputSingleFile();
        static const size_t bufferSize = 128*1024;
        //
        // definition of virtual methods from base class
        void output_for_one_cu(
                size_t cuID,          size_t ts,
                double load_vsm,      double load_rsm,
                double load_selfprod, double load_pv,
                double bs_SOC,        double load_bs,
                double load_hp,       double load_cs,
                size_t n_cars_pc,     size_t n_cars_pnc);
        void flush_buffer();
    private:
        mutex  single_file_mutex; ///< If a single file is selected as CU output, this variable holds the mutex to ensure correct concurrency behavior
        char*  buffer; // buffer for speedup
};

/**
 * This class represents the output for one
 * individual control unit.
 **/
class CUOutputOneFilePerCU : public CUOutput {
    public:
        CUOutputOneFilePerCU(int cuID, filesystem::path& dirpath);
        //
        // definition of virtual methods from base class
        void output_for_one_cu(
                size_t cuID,          size_t ts,
                double load_vsm,      double load_rsm,
                double load_selfprod, double load_pv,
                double bs_SOC,        double load_bs,
                double load_hp,       double load_cs,
                size_t n_cars_pc,     size_t n_cars_pnc);
        void flush_buffer();
};

/**
 * This class represents the output for all
 * control units, that are connected to one
 * substation.
 * 
 * This class is thread-safe.
 **/
class CUOutputOneFilePerSubstation : public CUOutput {
    public:
        CUOutputOneFilePerSubstation(const string* substName, filesystem::path& dirpath);
        //
        // definition of virtual methods from base class
        void output_for_one_cu(
                size_t cuID,         size_t ts,
                double load_vsm,      double load_rsm,
                double load_selfprod, double load_pv,
                double bs_SOC,        double load_bs,
                double load_hp,       double load_cs,
                size_t n_cars_pc,    size_t n_cars_pnc);
        void flush_buffer();
    private:
        std::mutex single_file_mutex; ///< If a single file is selected as CU output, this variable holds the mutex to ensure correct concurrency behavior
};

#endif

