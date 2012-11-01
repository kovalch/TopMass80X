C := g++ -O3 -Wall -fPIC
L := g++ 

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs) -lRooFit
ROOTLDFLAGS := $(shell root-config --ldflags)

CXX_FILES := $(wildcard src/*.cxx)
OBJ_FILES := $(patsubst src/%.cxx,obj/%.o,$(filter-out src/main.cxx, $(CXX_FILES)))


.PHONY: clean all

all: TopMass

clean:
	@rm -rf obj
	@rm -f lib/libTopMass.so
	@rm -f TopMass
obj:
	@mkdir -p obj


libTopMass.so:	obj $(OBJ_FILES) 
	$(L) $(ROOTLDFLAGS) -shared $(OBJ_FILES) -o lib/$@

# export LHAPATH=/afs/naf.desy.de/user/e/eschliec/wd/LHAPDF/share/lhapdf/PDFsets
# export LD_LIBRARY_PATH=/afs/naf.desy.de/user/e/eschliec/wd/LHAPDF/lib:$LD_LIBRARY_PATH
TopMass: obj  $(OBJ_FILES) obj/main.o
	$(L) $(ROOTLDFLAGS) $(ROOTLIBS) -Wl,-rpath -Wl,lib -lboost_program_options -L/afs/naf.desy.de/user/e/eschliec/wd/LHAPDF/lib -lLHAPDF -o $@ obj/main.o $(OBJ_FILES)
	@echo "done"

#obj/%.o: src/%.cxx
#	$(C) $(ROOTCFLAGS) -Iinc/ -c -o $@ $<

#rules
obj/%.o: src/%.cxx 
	$(C) $(ROOTCFLAGS) -Iinc/ -I/afs/naf.desy.de/user/e/eschliec/wd/LHAPDF/include  -MMD -c -o $@ $<
#@sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
#     -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.d

-include $(patsubst src/%.cxx,obj/%.d,$(CXX_FILES))
