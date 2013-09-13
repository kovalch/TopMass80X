C := g++ -O3 -Wall -fPIC -std=c++0x
L := g++ -O3

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLIBS    := $(shell root-config --libs) -lRooFit -lRooFitCore -lTreePlayer
ROOTLDFLAGS := $(shell root-config --ldflags)

# files for TopMassAnalyzer
CXX_FILES := $(wildcard src/*.cxx)
OBJ_FILES := $(patsubst src/%.cxx,obj/1/%.o,$(filter-out src/main.cxx, $(CXX_FILES)))
# files from CMSSW for classes in the tree
CXX_FILES_CMSSW := $(wildcard ../TopEventTree/src/*.cc)
OBJ_FILES_CMSSW := $(patsubst ../TopEventTree/src/%.cc,obj/2/%.o,$(filter-out ../TopEventTree/src/TreeRegistryService.cc, $(CXX_FILES_CMSSW)))
OBJ_FILES_CMSSW_DICT := $(patsubst ../TopEventTree/src/%.cc,obj/2/%_dict.o,$(filter-out ../TopEventTree/src/TreeRegistryService.cc, $(CXX_FILES_CMSSW)))
# files needed to run calibration program
CXX_FILES_CALIBRATION := $(wildcard root/*.cc)
OBJ_FILES_CALIBRATION := $(patsubst root/%.cc,obj/3/%.o,$(filter-out root/calibration_main.cc, $(CXX_FILES_CALIBRATION))) obj/1/ProgramOptionsReader.o

.PHONY: clean all top cal

all: top cal

top: TopMass

cal: Calibration

clean:
	@rm -rf obj
	@rm -f TopMass
obj:
	@mkdir -p obj/1
	@mkdir -p obj/2
	@mkdir -p obj/3

TopMass: obj $(OBJ_FILES) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT) obj/1/main.o
	$(L) $(ROOTLDFLAGS) $(ROOTLIBS) -Wl,-rpath -Wl,lib -L/afs/naf.desy.de/group/cms/sw/slc5_amd64_gcc462/external/boost/1.50.0-cms2/lib -lboost_program_options -o $@ obj/1/main.o $(OBJ_FILES) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT)
	@echo -e "\033[32m\033[1mSuccessfully compiled TopMass\033[0m"

Calibration: obj $(OBJ_FILES_CALIBRATION) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT) obj/3/calibration_main.o
	$(L) $(ROOTLDFLAGS) $(ROOTLIBS) -Wl,-rpath -Wl,lib -L/afs/naf.desy.de/group/cms/sw/slc5_amd64_gcc462/external/boost/1.50.0-cms2/lib -lboost_program_options -o $@ obj/3/calibration_main.o $(OBJ_FILES_CALIBRATION) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT)
	@echo -e "\033[32m\033[1mSuccessfully compiled Calibration\033[0m"

#rules
obj/1/%.o: src/%.cxx
	$(C) $(ROOTCFLAGS) -I/afs/naf.desy.de/group/cms/sw/slc5_amd64_gcc462/external/boost/1.50.0-cms2/include -I../.. -Iinc/ -MMD -c -o $@ $<
#@sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
#     -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.d

obj/2/%_dict.cc: ../TopEventTree/interface/%.h ../TopEventTree/interface/%_LinkDef.h
	rootcint -l -f $@ -c -DSTANDALONE $^

obj/2/%_dict.o: obj/2/%_dict.cc
	$(C) $(ROOTCFLAGS) -I../.. -I. -MMD -c $< -o $@

obj/2/%.o: ../TopEventTree/src/%.cc ../TopEventTree/interface/%.h
	$(C) $(ROOTCFLAGS) -I../.. -MMD -c -o $@ $<

obj/3/%.o: root/%.cc
	$(C) $(ROOTCFLAGS) -I../.. -Iinc/ -MMD -c -o $@ $<
