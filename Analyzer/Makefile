C := g++ -Wall
L := g++

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

CXX_FILES := $(wildcard src/*.cxx)
CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(patsubst src/%.cxx,obj/%.o,$(CXX_FILES))
OBJCPP_FILES := $(patsubst src/%.cpp,obj/%.o,$(CPP_FILES))

TopMass: $(OBJ_FILES) $(OBJCPP_FILES)
	$(C) $(ROOTLIBS) -o $@ $(OBJ_FILES) $(OBJCPP_FILES)

obj/%.o: src/%.cxx
	$(L) $(ROOTCFLAGS) -Iinc/ -c -o $@ $<

obj/%.o: src/%.cpp
	$(L) $(ROOTCFLAGS) -Iinc/ -c -o $@ $<
