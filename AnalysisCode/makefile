CXX      := c++
CXXSTD   := -std=c++17 -O2

FASTJET_CONFIG := /data/ALEPH/MC/mcgen/fastjet-install/bin/fastjet-config

ROOT_FLAGS    := $(shell root-config --cflags --libs)
FASTJET_FLAGS := $(shell $(FASTJET_CONFIG) --cxxflags --libs)

.PHONY: all clean

all: unbias_weights

unbias_weights: unbias_weights.cc subjet_basis.h basis_functions.h
	$(CXX) $(CXXSTD) unbias_weights.cc $(ROOT_FLAGS) $(FASTJET_FLAGS) -o unbias_weights

clean:
	rm -f unbias_weights