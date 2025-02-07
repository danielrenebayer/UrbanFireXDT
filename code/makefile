CC = g++
FLAGS = -std=c++23 -pedantic -Wall
FLAGS_EXTRA = -Wconversion -Wextra
LIBS  = -pthread -l sqlite3 -l boost_program_options
FILES = main.cpp components.cpp units.cpp global.cpp sac_planning.cpp setup_and_dataloading.cpp simulation_logic.cpp output.cpp helper.cpp worker_threads.cpp
OUTPUT_DIR = ../bin

# Optimizataion configuration: Select the requested optimization backend
#MACROS = -D USE_GUROBI
MACROS = -D USE_OR_TOOLS

# Reading additional input file
ifneq ("$(wildcard makefile_local)","")
include makefile_local
$(info Reading additional local settings from file makefile_local.)
endif

# Optimization settings
# Gurobi settings
ifeq ($(findstring USE_GUROBI,$(MACROS)),USE_GUROBI)
ifeq ($(GUROBI_PATH),)
$(error Variable GUROBI_PATH has not been set. Please define this variable in a file called makefile_local.)
endif
	# Set variables
	FILES += optimization_unit_gurobi.cpp
	GUROBI_LIBS = -lgurobi_c++ -lgurobi110 -lm
	GUROBI_INCLUDE_PATH = $(GUROBI_PATH)/include/
	GUROBI_LIB_PATH = $(GUROBI_PATH)/lib/
	OPTIMIZATION_OPTIONS = -I$(GUROBI_INCLUDE_PATH) -L$(GUROBI_LIB_PATH) $(GUROBI_LIBS)
endif
# OR-Tools settings
ifeq ($(findstring USE_OR_TOOLS,$(MACROS)),USE_OR_TOOLS)
ifeq ($(OR_TOOLS_PATH),)
$(error Variable OR_TOOLS_PATH has not been set. Please define this variable in a file called makefile_local.)
endif
	# Set variables
	OR_TOOLS_INCLUDE_PATH = $(OR_TOOLS_PATH)/include
	OR_TOOLS_LIB_PATH = $(OR_TOOLS_PATH)/lib
	OR_TOOLS_LIBS = -lortools -labsl_base -labsl_raw_logging_internal -labsl_spinlock_wait
	OPTIMIZATION_OPTIONS = -I$(OR_TOOLS_INCLUDE_PATH) -L$(OR_TOOLS_LIB_PATH) $(OR_TOOLS_LIBS)
endif


debug: directories
	$(CC) $(FILES) -o $(OUTPUT_DIR)/simulation-dbg $(FLAGS) $(FLAGS_EXTRA) $(MACROS) -g -D DEBUG -D ADD_METHOD_ACCESS_PROTECTION_VARS $(LIBS) $(OPTIMIZATION_OPTIONS)

directories:
	mkdir -p $(OUTPUT_DIR)

all: debug opti

opti: directories
	$(CC) $(FILES) -o $(OUTPUT_DIR)/simulation-opti $(FLAGS) -O3 $(MACROS) $(LIBS) $(OPTIMIZATION_OPTIONS)

verbose_debug: directories
	$(CC) $(FILES) -o $(OUTPUT_DIR)/simulation-dbg $(FLAGS) $(FLAGS_EXTRA) $(MACROS) -g -D DEBUG -D ADD_METHOD_ACCESS_PROTECTION_VARS -D DEBUG_EXTRA_OUTPUT $(LIBS) $(OPTIMIZATION_OPTIONS)

clean:
	rm -f $(OUTPUT_DIR)/simulation-dbg $(OUTPUT_DIR)/simulation-opti

print-vars:
	@$(foreach V,$(sort $(.VARIABLES)), $(if $(filter-out environment% default automatic, $(origin $V)), $(info $V=$($V))))

