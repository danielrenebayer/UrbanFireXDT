FLAGS = -std=c++20 -l sqlite3 -pedantic -Wall
FILES = main.cpp components.cpp units.cpp global.cpp

all:
	g++ $(FILES) -o main $(FLAGS)

debug:
	g++ $(FILES) -o main-dbg $(FLAGS) -g -D DEBUG

verbose_debug:
	g++ $(FILES) -o main-dbg $(FLAGS) -g -D DEBUG -D DEBUG_EXTRA_OUTPUT

clean:
	rm main main-dbg

