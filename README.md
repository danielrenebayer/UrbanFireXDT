# UrbanFireXDT: the Urban Future Integrated Renewable Energy systems eXploration Digital Twin

<!-- Project logo -->
<p align="center">
    <img src="documentation/logos/UrbanFireXDT_Logo_v1.png" alt="UrbanFireXDT Logo" />
</p>

<!-- Badges -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17302743.svg)](https://doi.org/10.5281/zenodo.17302743)
![C++23](https://img.shields.io/badge/C++-23-blue.svg?logo=c%2B%2B)

<!-- About the repository -->
## About the repository

This repository contains the code for the simulation component that is part of the digital twin for the analysis of future integrated renewable energy systems, as presented for example in [1].
For details, please read the scientific publications linked below.
This README contains all important information to compile the code and run the examples.
Parts of this code were developed as part of the *DigiSWM* research project founded by the Bavarian State Ministry of Economic Affairs, Regional Development and Energy under grant number DIK-2103-0017 / DIK0298/02.

<p align="center">
    <img src="documentation/logos/DigiSWM_Logo_v3_768x413px.png" alt="Logo of the research project DigiSWM" style="max-width: 200px; width: 50%; height: auto;"/>
</p>


<!-- Structure of the repository -->
## Structure of the repository

- The directory `code` contains all C++ code
- The directory `config-example` contains an example config file
- The directory `data` in an dummy directory for storing input data and (if not changed in the config) output data
- The directory `documentation` contains the complete input and output documentation as well as a documentation of all configuration variables
- The directory `test` contains a small example serving as an exemplary setup and for testing the code


<!-- Installation -->
## Installation and Requirements

This code is mainly written in C++ in standard **C++23** and developed for a **Debian or Ubuntu** environment.
For compilation of the code, **gcc in version 14 or above** is required.

### Dependencies

This project requires the following libraries to build and run:
* **libsqlite3-dev:** For SQLite3 database interaction
* **libboost-dev:** For parsing the JSON-based configuration file and the command line options
* **libboost-program-options-dev:** Specifically for command-line option parsing
* **libboost-json-dev:** Specifically for parsing JSON files (that are used as cache files)
* **pybind11-dev:** Exporting the simulation as python module (only required for compilation target `python_module`)
* **python3-pybind11:** The target python module (only required for compilation target `python_module`)

These dependencies can be installed on Debian/Ubuntu with:
```
apt-get install libsqlite3-dev libboost-dev libboost-program-options-dev
```

### Requirements for optimization-based strategies

If the program should include an optimization-based control strategy, an LP-solver must be installed on the host.
The simulation currently supports two optimization backends:
- Gurobi
- Google OR-Tools (using the internal GLOP solver)

### Installation of OR-Tools

- Download Google OR-Tools from [https://developers.google.com/optimization/install/cpp/binary_linux](https://developers.google.com/optimization/install/cpp/binary_linux)
- Unpack the OR-Tools and save the set global variable `OPTIMIZER_PATH`
- Install dependencies for the OR-Tools, i.e., `libprotobuf-dev`
```
apt-get install libprotobuf-dev
```

<!-- Usage options -->
## Usage options

This simulation can be used in two ways:

- **Standalone executable**: run directly from the command line, using configuration files and command-line options to define scenarios.  
  This mode is well suited for batch experiments, scenario studies, and reproducible evaluations.

- **Python module**: import the simulation into Python, e.g. for integration with data analysis pipelines or the training of reinforcement learning (RL)-based control strategies.  
  This mode provides step-by-step control and access to internal simulation states.

For both versions, examples can be found in the `test` folder.


<!-- First steps -->
## First steps

1. Install all dependencies including at least one of the optimizers
2. Modify `code/makefile` to set the correct optimizer and optimizer path
3. Add optimizer path to `LD_LIBRARY_PATH` using `export LD_LIBRARY_PATH=$OPTIMIZER_PATH/lib:LD_LIBRARY_PATH` (assuming that `OPTIMIZER_PATH` is manually set correctly)
4. Run `cd code; make all`
5. Run all tests `cd test; ./run-all-tests.sh`


<!-- Compilation options: Make targets -->
## Compilation options: Make targets

The makefile `code/makefile` contains different targets.

| Target name | Description |
| ---         | ---         |
| `debug`     | Compiles the code with all debugging symbols (for `gdb`) and without any optimizaation. Additionally, the access protection variables will be included that ensure the correct chronological sequence of method calls for the individual components. |
| `opti`      | Compiles the code for the final usage. Includes code optimization, and no access protection stuff. |
| `all`       | `debug` and `opti` together |
| `python_module` | Compile the simulation as python module. Not part of the target `all`. |
| `verbose_debug` | Like `debug`, with extra output |
| `clean`     | Removes compiled files in the output directory |
| `print-vars` | Prints out all finally (i.e., after reading both the makefile and the local makefile) effective configuration variables |

**Please note that all variables inside `code/makefile` can be overwritten with an local makefile `code/makefile_local`, if present.**
If using a local makefile, the local paths and settings will be ignored by git, as it is in the list of ignored files.


<!-- Scientific References -->
## Scientific References
[1] D. Bayer and M. Pruckner, "A digital twin of a local energy system based on real smart meter data," *Energy Informatics*, vol. 6, no. 1, Mar. 2023. doi: [10.1186/s42162-023-00263-6](https://dx.doi.org/10.1186/s42162-023-00263-6).


<!-- LICENSE -->
## License

This code is distributed under the MIT License. See `LICENSE` for more details.


<!-- Logo License -->
## Logo License

The logos included in `documentation/logos` are distributed under the Creative Commons Attribution-NoDerivatives (CC BY-ND) license.
See [https://creativecommons.org/licenses/by-nd/4.0/](https://creativecommons.org/licenses/by-nd/4.0/) for more details.


