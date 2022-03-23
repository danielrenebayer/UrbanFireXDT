/*
 * output.h
 *
 * This contains all functions and classes for buffering outputs
 * and writing them to the disk finally.
 *
 */


#include <fstream>
#include <mutex>
#include <sstream>

using namespace std;


namespace output {

    inline std::ofstream* substation_output;
    inline std::ofstream* cu_output;
    inline mutex cu_output_mutex;

    void initializeSubstationOutput(int scenario_id);
    void initializeCUOutput(int scenario_id);

    void closeOutputs();

    void flushBuffers();

    void output_for_one_cu(int ts, int cuID, float load_vsm, float load_rsm, float load_selfprod, float load_pv, float bs_SOC, float load_bs);

    // TODO: these lines are a speedup test, only
    inline const size_t bufferSize = 64*1024;
    inline char buffer[bufferSize];

}

