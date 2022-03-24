#include "output.h"

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
	output_path_subst << "../data/output/";
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
	*(substation_output) << ",total_load" << endl;
	// TODO: make output path configurable and use absolute path
}

void output::initializeCUOutput(int scenario_id) {
	//
    // This method initializes the control unit
    // output file.
    //
    //
    // initialize the output file
    stringstream output_path_CUs;
	output_path_CUs << "../data/output/";
	output_path_CUs << setw(4) << setfill('0') << scenario_id;
	output_path_CUs << "-CU-time-series.csv";
	cu_output = new ofstream(output_path_CUs.str().c_str(), std::ofstream::out);
	// activate buffers for speedup
	cu_output->rdbuf()->pubsetbuf(buffer, bufferSize);
    //
	// add header to output file
	*(cu_output) << "Timestep,ControlUnitID,Load_vSmartMeter_kW,Load_rSmartMeters_kW,Load_self_produced_kW,PVFeedIn_Simulated_kW,BESS_SOC,BESS_load_kW" << endl;
	// TODO: make output path configurable and use absolute path
}

void output::closeOutputs() {
	// close output file for substations
    substation_output->close();
    delete substation_output;

    //close output file for control units
	cu_output->close();
	delete cu_output;
}

void output::flushBuffers() {
	substation_output->flush();
	unique_lock lock(cu_output_mutex); // secure access by using a mutex
	cu_output->flush();
}

void output::output_for_one_cu(int ts, int cuID, float load_vsm, float load_rsm, float load_selfprod, float load_pv, float bs_SOC, float load_bs) {
	unique_lock lock(cu_output_mutex); // secure access by using a mutex
	*(cu_output) << ts << "," << cuID << "," << load_vsm << "," << load_rsm << "," << load_selfprod << "," << load_pv << "," << bs_SOC << "," << load_bs << "\n";
}
