.PHONY: all

FLAGS=-ggdb -DTOOLCHAIN_CLANG -W -Wall -fvisibility=hidden -I/Users/kpenev/projects/git/AstroWISP/src/third_party_libs/Eigen
LIBS=-L$(HOME)/lib/ -lastrowisp

all: noiseless_psffit

noiseless_psffit: noiseless_psffit.cpp
	g++ $(FLAGS) noiseless_psffit.cpp -o noiseless_psffit  $(LIBS)

apphot: apphot.cpp
	g++ $(FLAGS) apphot.cpp -o apphot $(LIBS)

io_tree_updates: test_io_tree_updates.cpp
	g++ $(FLAGS) test_io_tree_updates.cpp -o test_io_tree_updates $(LIBS)
