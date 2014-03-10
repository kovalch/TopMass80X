C := g++ -O3 -Wall -fPIC -std=c++0x
L := g++ -O3

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLIBS    := $(shell root-config --libs)
ROOTLDFLAGS := $(shell root-config --ldflags)

# files needed to run add PDF weights program
CXX_FILES_ADDPDFWEIGHTS := $(wildcard addPDFweights*.cc)
OBJ_FILES_ADDPDFWEIGHTS := $(patsubst %.cc,obj/1/%.o,$(filter-out $(wildcard *_main.cc), $(CXX_FILES_ADDPDFWEIGHTS)))

.PHONY: clean all

all: AddPDFweights

clean:
	@rm -rf obj
	@rm -f AddPDFweights
obj:
	@mkdir -p obj/1

AddPDFweights: obj $(OBJ_FILES) $(OBJ_FILES_ADDPDFWEIGHTS) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT) obj/1/addPDFweights_main.o
	$(L) $(ROOTLDFLAGS) $(ROOTLIBS) $(ROOFITLIBS) -Wl,-rpath -Wl,lib -L/nfs/dust/cms/group/tophh/LHAPDF-605/lib -lLHAPDF -o $@ obj/1/addPDFweights_main.o $(OBJ_FILES_ADDPDFWEIGHTS) $(OBJ_FILES_CMSSW) $(OBJ_FILES_CMSSW_DICT)
	@echo -e "\033[32m\033[1mSuccessfully compiled AddPDFweights\033[0m"

#rules
obj/1/%.o: %.cc
	$(C) $(ROOTCFLAGS) -I/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/boost/1.51.0/include -I/nfs/dust/cms/group/tophh/LHAPDF-605/include -I../.. -Iinc/ -MMD -c -o $@ $<