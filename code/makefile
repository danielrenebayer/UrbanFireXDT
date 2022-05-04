FLAGS = -std=c++20 -pedantic -Wall
LIBS  = -l sqlite3 -l boost_program_options
FILES = main.cpp components.cpp units.cpp global.cpp simulation_setup.cpp simulation_logic.cpp output.cpp helper.cpp

default:
	g++ $(FILES) -o main $(FLAGS) $(LIBS)

all: debug default opti

opti:
	g++ $(FILES) -o main-opti $(FLAGS) -O3 $(LIBS)

debug:
	g++ $(FILES) -o main-dbg $(FLAGS) -g -D DEBUG $(LIBS)

verbose_debug:
	g++ $(FILES) -o main-dbg $(FLAGS) -g -D DEBUG -D DEBUG_EXTRA_OUTPUT $(LIBS)

clean:
	rm -f main main-dbg main-opti

