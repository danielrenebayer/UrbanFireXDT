FLAGS = -std=c++23 -pedantic -Wall
FLAGS_EXTRA = -Wconversion -Wextra
LIBS  = -pthread -l sqlite3 -l boost_program_options
FILES = main.cpp components.cpp units.cpp global.cpp sac_planning.cpp setup_and_dataloading.cpp simulation_logic.cpp output.cpp helper.cpp worker_threads.cpp
OUTPUT_DIR = ../bin

# Optimizataion configuration: Select the requested optimization backend
MACROS = -D USE_GUROBI
#MACROS = -D USE_OR_TOOLS
#MACROS = -D USE_GLPK
# Gurobi settings
ifeq ($(findstring USE_GUROBI,$(MACROS)),USE_GUROBI)
	# Load additional file with paths
ifneq ("$(wildcard makefile_paths.sh)","")
include makefile_paths.sh
else
$(error makefile_paths.sh file not found. Please create this file holding the path to your Gurobi installation inside GUROBI_PATH.)
endif
ifeq ($(GUROBI_PATH),)
$(error Variable GUROBI_PATH has not been set in makefile_paths.sh. Please add this variable there.)
endif
	# Set variables
	FILES += optimization_unit_gurobi.cpp
	GUROBI_LIBS = -lgurobi_c++ -lgurobi110 -lm
	GUROBI_INCLUDE_PATH = $(GUROBI_PATH)/include/
	GUROBI_LIB_PATH = $(GUROBI_PATH)/lib/
	OPTIMIZATION_OPTIONS = -I$(GUROBI_INCLUDE_PATH) -L$(GUROBI_LIB_PATH) $(GUROBI_LIBS)
endif
# GLPK settings
ifeq ($(findstring USE_GLPK,$(MACROS)),USE_GLPK)
	OPTIMIZATION_OPTIONS = -lglpk
endif
# OR-Tools settings
ifeq ($(findstring USE_OR_TOOLS,$(MACROS)),USE_OR_TOOLS)
# Load additional file with paths
ifneq ("$(wildcard makefile_paths.sh)","")
include makefile_paths.sh
else
$(error makefile_paths.sh file not found. Please create this file holding the path to your OR-Tools installation inside OR_TOOLS_PATH.)
endif
ifeq ($(OR_TOOLS_PATH),)
$(error Variable OR_TOOLS_PATH has not been set in makefile_paths.sh. Please add this variable there.)
endif
	# Set variables
	OR_TOOLS_INCLUDE_PATH = $(OR_TOOLS_PATH)/include
	OR_TOOLS_LIB_PATH = $(OR_TOOLS_PATH)/lib
	OR_TOOLS_LIBS = -lortools -labsl_base -labsl_raw_logging_internal -labsl_spinlock_wait
	OPTIMIZATION_OPTIONS = -I$(OR_TOOLS_INCLUDE_PATH) -L$(OR_TOOLS_LIB_PATH) $(OR_TOOLS_LIBS)
endif


debug: directories
	g++ $(FILES) -o $(OUTPUT_DIR)/simulation-dbg $(FLAGS) $(FLAGS_EXTRA) $(MACROS) -g -D DEBUG -D ADD_METHOD_ACCESS_PROTECTION_VARS $(LIBS) $(OPTIMIZATION_OPTIONS)

directories:
	mkdir -p $(OUTPUT_DIR)

all: debug opti

opti: directories
	g++ $(FILES) -o $(OUTPUT_DIR)/simulation-opti $(FLAGS) -O3 $(MACROS) $(LIBS) $(OPTIMIZATION_OPTIONS)

verbose_debug: directories
	g++ $(FILES) -o $(OUTPUT_DIR)/simulation-dbg $(FLAGS) $(FLAGS_EXTRA) $(MACROS) -g -D DEBUG -D ADD_METHOD_ACCESS_PROTECTION_VARS -D DEBUG_EXTRA_OUTPUT $(LIBS) $(OPTIMIZATION_OPTIONS)

clean:
	rm -f $(OUTPUT_DIR)/simulation-dbg $(OUTPUT_DIR)/simulation-opti

