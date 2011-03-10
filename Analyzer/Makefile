C := g++ -Wall
L := g++

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

CXX_FILES := $(wildcard src/*.cxx)
OBJ_FILES := $(patsubst src/%.cxx,obj/%.o,$(CXX_FILES))

Analysis: $(OBJ_FILES)
	$(C) $(ROOTLIBS) -o $@ $(OBJ_FILES)

obj/%.o: src/%.cxx
	$(L) $(ROOTCFLAGS) -Iinc/ -c -o $@ $<

