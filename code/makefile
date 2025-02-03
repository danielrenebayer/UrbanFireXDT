FLAGS = -std=c++23 -pedantic -Wall
FLAGS_EXTRA = -Wconversion -Wextra
LIBS  = -pthread -l sqlite3 -l boost_program_options
FILES = main.cpp components.cpp units.cpp global.cpp sac_planning.cpp setup_and_dataloading.cpp simulation_logic.cpp output.cpp helper.cpp worker_threads.cpp
OUTPUT_DIR = ../bin

# Stuff for optimization
MACROS = -D USE_GUROBI
#MACROS = -D USE_GLPK
#MACROS = -D USE_COINOR
# Gurobi settings
ifeq ($(findstring USE_GUROBI,$(MACROS)),USE_GUROBI)
	FILES += optimization_unit_gurobi.cpp
	GUROBI_LIBS = -lgurobi_c++ -lgurobi110 -lm
	GUROBI_INCLUDE_PATH = /home/daniel/User-Software/gurobi1103/linux64/include/
	GUROBI_LIB_PATH = /home/daniel/User-Software/gurobi1103/linux64/lib/
	OPTIMIZATION_OPTIONS = -I$(GUROBI_INCLUDE_PATH) -L$(GUROBI_LIB_PATH) $(GUROBI_LIBS)
endif
# GLPK settings
ifeq ($(findstring USE_GLPK,$(MACROS)),USE_GLPK)
	OPTIMIZATION_OPTIONS = -lglpk
endif
# COIN-OR settings
ifeq ($(findstring USE_COINOR,$(MACROS)),USE_COINOR)
	OPTIMIZATION_OPTIONS = -lClp -lCoinUtils -lOsiClp -lOsi
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

