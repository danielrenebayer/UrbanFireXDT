FLAGS = -std=c++23 -pedantic -Wall
FLAGS_EXTRA = -Wconversion -Wextra
LIBS  = -pthread -l sqlite3 -l boost_program_options
FILES = main.cpp components.cpp units.cpp global.cpp sac_planning.cpp setup_and_dataloading.cpp simulation_logic.cpp output.cpp helper.cpp vehicles.cpp worker_threads.cpp

default:
	g++ $(FILES) -o main $(FLAGS) $(LIBS)

all: debug default opti

opti:
	g++ $(FILES) -o main-opti $(FLAGS) -O3 $(LIBS)

debug:
	g++ $(FILES) -o main-dbg $(FLAGS) $(FLAGS_EXTRA) -g -D DEBUG $(LIBS)

verbose_debug:
	g++ $(FILES) -o main-dbg $(FLAGS) -g -D DEBUG -D DEBUG_EXTRA_OUTPUT $(LIBS)

clean:
	rm -f main main-dbg main-opti

