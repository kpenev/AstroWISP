.PHONY: all

FLAGS=-ggdb -DTOOLCHAIN_CLANG -W -Wall -fvisibility=hidden -I/Users/kpenev/projects/git/SuperPhot/src/third_party_libs/Eigen
LIBS=-L$(HOME)/lib/ -lsuperphot

all: noiseless_psffit

noiseless_psffit: noiseless_psffit.cpp
	g++ $(FLAGS) noiseless_psffit.cpp -o noiseless_psffit  $(LIBS)

apphot: apphot.cpp
	g++ $(FLAGS) apphot.cpp -o apphot $(LIBS)
