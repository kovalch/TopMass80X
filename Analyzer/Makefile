C := g++ -g -Wall -fPIC
L := g++ -g

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLIBS    := $(shell root-config --libs) -lRooFit
ROOTLDFLAGS := $(shell root-config --ldflags)

# files for TopMassAnalyzer
CXX_FILES := $(wildcard src/*.cxx)
OBJ_FILES := $(patsubst src/%.cxx,obj/%.o,$(filter-out src/main.cxx, $(CXX_FILES)))
# files from CMSSW for classes in the tree
CXX_FILES_CMSSW := $(wildcard ../TopEventTree/src/*.cc)
OBJ_FILES_CMSSW := $(patsubst ../TopEventTree/src/%.cc,obj2/%.o,$(filter-out ../TopEventTree/src/TreeRegistryService.cc, $(CXX_FILES_CMSSW)))
OBJ_FILES_CMSSW_DICT := $(patsubst ../TopEventTree/src/%.cc,obj2/%_dict.o,$(filter-out ../TopEventTree/src/TreeRegistryService.cc, $(CXX_FILES_CMSSW)))
# files needed to run calibration program
CXX_FILES_CALIBRATION := $(wildcard root/*.cc)
OBJ_FILES_CALIBRATION := $(patsubst root/%.cc,obj3/%.o,$(filter-out root/calibration_main.cc, $(CXX_FILES_CALIBRATION)))

.PHONY: clean all top cal

all: top cal

top: TopMass

cal: Calibration

clean:
	@rm -rf obj
	@rm -rf obj2
	@rm -rf obj3
	@rm -f TopMass
obj:
	@mkdir -p obj
	@mkdir -p obj2
	@mkdir -p obj3

TopMass: obj $(OBJ_FILES) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT) obj/main.o
	$(L) $(ROOTLDFLAGS) $(ROOTLIBS) -Wl,-rpath -Wl,lib -lboost_program_options -L/afs/naf.desy.de/user/e/eschliec/wd/LHAPDF/lib -lLHAPDF -o $@ obj/main.o $(OBJ_FILES) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT)
	@echo -e "\033[32m\033[1mSuccessfully compiled TopMass\033[0m"

Calibration: obj $(OBJ_FILES_CALIBRATION) $(OBJ_FILES) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT) obj3/calibration_main.o
	$(L) $(ROOTLDFLAGS) $(ROOTLIBS) -Wl,-rpath -Wl,lib -lboost_program_options -L/afs/naf.desy.de/user/e/eschliec/wd/LHAPDF/lib -lLHAPDF -o $@ obj3/calibration_main.o $(OBJ_FILES_CALIBRATION) $(OBJ_FILES) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT)
	@echo -e "\033[32m\033[1mSuccessfully compiled Calibration\033[0m"

#rules
obj/%.o: src/%.cxx
	$(C) $(ROOTCFLAGS) -I../.. -Iinc/ -I/afs/naf.desy.de/user/e/eschliec/wd/LHAPDF/include -MMD -c -o $@ $<
#@sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
#     -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.d

obj2/%_dict.cc: ../TopEventTree/interface/%.h ../TopEventTree/interface/%_LinkDef.h
	rootcint -l -f $@ -c -DSTANDALONE $^

obj2/%_dict.o: obj2/%_dict.cc
	$(C) $(ROOTCFLAGS) -I../.. -I. -MMD -c $< -o $@

obj2/%.o: ../TopEventTree/src/%.cc
	$(C) $(ROOTCFLAGS) -I../.. -MMD -c -o $@ $<

obj3/%.o: root/%.cc
	$(C) $(ROOTCFLAGS) -I../.. -Iinc/ -MMD -c -o $@ $<
	