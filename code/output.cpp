#include "output.h"

#include <fstream>
#include <iomanip>
#include <sstream>

#include "global.h"

using namespace std;

void output::initializeSubstationOutput(int scenario_id) {
    //
    // This method initializes the substation output
    // file. It is destructed by the vacuum-function from namespace
    // global in global.h.
    //
    //
    // initialize the output file
    stringstream output_path_subst;
	output_path_subst << "../data/output/";
	output_path_subst << setw(4) << setfill('0') << scenario_id;
	output_path_subst << ".csv";
	global::substation_output = new ofstream(output_path_subst.str().c_str(), std::ofstream::out);
    //
	// add header to output file
	*(global::substation_output) << "Timestep";
	Substation*const* subList = Substation::GetArrayOfInstances();
    const int nSubst = Substation::GetNumberOfInstances();
	for (int i = 0; i < nSubst; i++) {
        *(global::substation_output) << "," << subList[i]->get_name()->c_str();
	}
	*(global::substation_output) << ",total_load" << endl;
	// TODO: make output path configurable and use absolute path
	global::substation_output_init = true;
}

