C := g++ -O3 -Wall -fPIC -std=c++0x
L := g++ -O3 

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLIBS    := $(shell root-config --libs) -lTreePlayer
ROOFITLIBS  := -lRooFit -lRooFitCore
ROOTLDFLAGS := $(shell root-config --ldflags)

# files for TopMassAnalyzer
CXX_FILES := $(wildcard src/*.cxx)
OBJ_FILES := $(patsubst src/%.cxx,obj/1/%.o,$(filter-out src/main.cxx, $(CXX_FILES)))
# files from CMSSW for classes in the tree
CXX_FILES_CMSSW := $(wildcard ../TopEventTree/src/*.cc)
OBJ_FILES_CMSSW := $(patsubst ../TopEventTree/src/%.cc,obj/2/%.o,$(filter-out ../TopEventTree/src/TreeRegistryService.cc, $(CXX_FILES_CMSSW)))
OBJ_FILES_CMSSW_DICT := $(patsubst ../TopEventTree/src/%.cc,obj/2/%_dict.o,$(filter-out ../TopEventTree/src/TreeRegistryService.cc, $(CXX_FILES_CMSSW)))
# files needed to run template derivation program
CXX_FILES_CALIBRATION := $(wildcard root/calibration*.cc)
OBJ_FILES_CALIBRATION := $(patsubst root/%.cc,obj/3/%.o,$(filter-out $(wildcard root/*_main.cc), $(CXX_FILES_CALIBRATION))) obj/1/ProgramOptionsReader.o
# files needed to run control plot program
CXX_FILES_CONTROLPLOT := $(wildcard root/controlPlots*.cc)
OBJ_FILES_CONTROLPLOT := $(patsubst root/%.cc,obj/3/%.o,$(filter-out $(wildcard root/*_main.cc), $(CXX_FILES_CONTROLPLOT))) obj/1/ProgramOptionsReader.o obj/1/Helper.o
# files needed to run duplicate removal program
CXX_FILES_DUPREMOVAL := $(wildcard root/duplicateRemoval*.cc)
OBJ_FILES_DUPREMOVAL := $(patsubst root/%.cc,obj/3/%.o,$(filter-out $(wildcard root/*_main.cc), $(CXX_FILES_DUPREMOVAL))) obj/1/ProgramOptionsReader.o
# files needed to run infinity removal program
CXX_FILES_INFREMOVAL := $(wildcard root/infRemoval*.cc)
OBJ_FILES_INFREMOVAL := $(patsubst root/%.cc,obj/3/%.o,$(filter-out $(wildcard root/*_main.cc), $(CXX_FILES_INFREMOVAL))) obj/1/ProgramOptionsReader.o
# files needed to run skimmer program
CXX_FILES_SKIMMER := $(wildcard root/skimmer*.cc)
OBJ_FILES_SKIMMER := $(patsubst root/%.cc,obj/3/%.o,$(filter-out $(wildcard root/*_main.cc), $(CXX_FILES_SKIMMER))) obj/1/ProgramOptionsReader.o
# files needed to run add PDF weights program
CXX_FILES_ADDPDFWEIGHTS := $(wildcard root/addPDFweights*.cc)
OBJ_FILES_ADDPDFWEIGHTS := $(patsubst root/%.cc,obj/4/%.o,$(filter-out $(wildcard root/*_main.cc), $(CXX_FILES_ADDPDFWEIGHTS))) obj/1/ProgramOptionsReader.o obj/1/Helper.o
# files needed to run SystematicUncertainties program
CXX_FILES_SYSTEMATICUNCERTAINTIES := $(wildcard root/SystematicUncertainties*.cc)
OBJ_FILES_SYSTEMATICUNCERTAINTIES := $(patsubst root/%.cc,obj/3/%.o,$(filter-out $(wildcard root/*_main.cc), $(CXX_FILES_SYSTEMATICUNCERTAINTIES)))

.PHONY: clean all top tem con dup ski pdf sys inf

all: top tem con dup ski pdf sys inf

top: TopMass

tem: TemplateDerivation

con: ControlPlots

dup: DuplicateRemoval

inf: InfRemoval

ski: Skimmer

pdf: AddPDFweights

sys: SystematicUncertainties

clean:
	@rm -rf obj
	@rm -f TopMass
	@rm -f TemplateDerivation
	@rm -f ControlPlots
	@rm -f DuplicateRemoval
	@rm -f InfRemoval
	@rm -f Skimmer
	@rm -f AddPDFweights
	@rm -f SystematicUncertainties
obj:
	@mkdir -p obj/1
	@mkdir -p obj/2
	@mkdir -p obj/3
	@mkdir -p obj/4

TopMass: obj $(OBJ_FILES) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT) obj/1/main.o
	$(L) $(ROOTLDFLAGS) $(ROOTLIBS) $(ROOFITLIBS) -Wl,-rpath -Wl,lib -L/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/boost/1.51.0/lib -lboost_program_options -o $@ obj/1/main.o $(OBJ_FILES) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT)
	@echo -e "\033[32m\033[1mSuccessfully compiled TopMass\033[0m"

TemplateDerivation: obj $(OBJ_FILES) $(OBJ_FILES_CALIBRATION) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT) obj/3/calibration_main.o
	$(L) $(ROOTLDFLAGS) $(ROOTLIBS) $(ROOFITLIBS) -Wl,-rpath -Wl,lib -L/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/boost/1.51.0/lib -lboost_program_options -o $@ obj/3/calibration_main.o $(OBJ_FILES_CALIBRATION) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT)
	@echo -e "\033[32m\033[1mSuccessfully compiled TemplateDerivation\033[0m"

ControlPlots: obj $(OBJ_FILES_CONTROLPLOT) obj/3/controlPlots_main.o
	$(L) $(ROOTLDFLAGS) $(ROOTLIBS) -Wl,-rpath -Wl,lib -L/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/boost/1.51.0/lib -lboost_program_options -o $@ obj/3/controlPlots_main.o $(OBJ_FILES_CONTROLPLOT)
	@echo -e "\033[32m\033[1mSuccessfully compiled ControlPlots\033[0m"

DuplicateRemoval: obj $(OBJ_FILES) $(OBJ_FILES_DUPREMOVAL) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT) obj/3/duplicateRemoval_main.o
	$(L) $(ROOTLDFLAGS) $(ROOTLIBS) -Wl,-rpath -Wl,lib -L/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/boost/1.51.0/lib -lboost_program_options -o $@ obj/3/duplicateRemoval_main.o $(OBJ_FILES_DUPREMOVAL) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT)
	@echo -e "\033[32m\033[1mSuccessfully compiled DuplicateRemoval\033[0m"

InfRemoval: obj $(OBJ_FILES) $(OBJ_FILES_INFREMOVAL) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT) obj/3/infRemoval_main.o
	$(L) $(ROOTLDFLAGS) $(ROOTLIBS) -Wl,-rpath -Wl,lib -L/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/boost/1.51.0/lib -lboost_program_options -o $@ obj/3/infRemoval_main.o $(OBJ_FILES_INFREMOVAL) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT)
	@echo -e "\033[32m\033[1mSuccessfully compiled InfRemoval\033[0m"

Skimmer: obj $(OBJ_FILES) $(OBJ_FILES_SKIMMER) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT) obj/3/skimmer_main.o
	$(L) $(ROOTLDFLAGS) $(ROOTLIBS) -Wl,-rpath -Wl,lib -L/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/boost/1.51.0/lib -lboost_program_options -o $@ obj/3/skimmer_main.o $(OBJ_FILES_SKIMMER) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT)
	@echo -e "\033[32m\033[1mSuccessfully compiled Skimmer\033[0m"

AddPDFweights: obj $(OBJ_FILES) $(OBJ_FILES_ADDPDFWEIGHTS) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT) obj/4/addPDFweights_main.o
	$(L) $(ROOTLDFLAGS) $(ROOTLIBS) -Wl,-rpath -Wl,lib -L/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/boost/1.51.0/lib -lboost_program_options -L/afs/desy.de/user/e/eschliec/xxl-af-cms/wd/LHAPDF-605/lib -lLHAPDF -o $@ obj/4/addPDFweights_main.o $(OBJ_FILES_ADDPDFWEIGHTS) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT)
	@echo -e "\033[32m\033[1mSuccessfully compiled AddPDFweights\033[0m"

SystematicUncertainties: obj $(OBJ_FILES_SYSTEMATICUNCERTAINTIES) obj/3/SystematicUncertainties_main.o
	$(L) $(ROOTLDFLAGS) $(ROOTLIBS) -Wl,-rpath -Wl,lib -o $@ obj/3/SystematicUncertainties_main.o $(OBJ_FILES_SYSTEMATICUNCERTAINTIES)
	@echo -e "\033[32m\033[1mSuccessfully compiled SystematicUncertainties\033[0m"

#rules
obj/1/%.o: src/%.cxx
	$(C) $(ROOTCFLAGS) -I/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/boost/1.51.0/include -I../.. -Iinc/ -MMD -c -o $@ $<
#@sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
#     -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.d

obj/2/%_dict.cc: ../TopEventTree/interface/%.h ../TopEventTree/interface/%_LinkDef.h
	rootcint -l -f $@ -c -DSTANDALONE $^

obj/2/%_dict.o: obj/2/%_dict.cc
	$(C) $(ROOTCFLAGS) -I../.. -I. -MMD -c $< -o $@

obj/2/%.o: ../TopEventTree/src/%.cc ../TopEventTree/interface/%.h
	$(C) $(ROOTCFLAGS) -I../.. -MMD -c -o $@ $<

obj/3/%.o: root/%.cc
	$(C) $(ROOTCFLAGS) -I/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/boost/1.51.0/include -I../.. -Iinc/ -MMD -c -o $@ $<

obj/4/%.o: root/%.cc
	$(C) $(ROOTCFLAGS) -I/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/boost/1.51.0/include -I/afs/desy.de/user/e/eschliec/xxl-af-cms/wd/LHAPDF-605/include -I../.. -Iinc/ -MMD -c -o $@ $<
