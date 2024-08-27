FLAGS = -std=c++23 -pedantic -Wall
FLAGS_EXTRA = -Wconversion -Wextra
LIBS  = -pthread -l sqlite3 -l boost_program_options
FILES = main.cpp components.cpp units.cpp global.cpp sac_planning.cpp setup_and_dataloading.cpp simulation_logic.cpp output.cpp helper.cpp vehicles.cpp worker_threads.cpp
OUTPUT_DIR = ../bin

directories:
	mkdir -p $(OUTPUT_DIR)

default: directories
	g++ $(FILES) -o $(OUTPUT_DIR)/simulation $(FLAGS) $(LIBS)

all: debug default opti

opti: directories
	g++ $(FILES) -o $(OUTPUT_DIR)/simulation-opti $(FLAGS) -O3 $(LIBS)

debug: directories
	g++ $(FILES) -o $(OUTPUT_DIR)/simulation-dbg $(FLAGS) $(FLAGS_EXTRA) -g -D DEBUG $(LIBS)

verbose_debug: directories
	g++ $(FILES) -o $(OUTPUT_DIR)/simulation-dbg $(FLAGS) -g -D DEBUG -D DEBUG_EXTRA_OUTPUT $(LIBS)

clean:
	rm -f $(OUTPUT_DIR)/simulation $(OUTPUT_DIR)/simulation-dbg $(OUTPUT_DIR)/simulation-opti

