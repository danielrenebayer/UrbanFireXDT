FLAGS = -std=c++20 -l sqlite3 -pedantic -Wall
FILES = main.cpp components.cpp units.cpp

all:
	g++ $(FILES) -o main $(FLAGS)

debug:
	g++ $(FILES) -o main-dbg $(FLAGS) -g

clean:
	rm main main-dbg

