#include <vector>
#include <cmath>
#include <fstream>
#include <time.h>
#include <boost/program_options.hpp>

#include "Analysis.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TH1.h"
#include "TH3F.h"
#include "TVector.h"
#include "TLegend.h"

#include "tinyxml.h"

#include "Helper.h"

namespace po = boost::program_options;


bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}

class TopMass {
  private:
    TString fMethod;
    double fLumi;
    
    Analysis* aSim;
      
  public:
    TopMass(po::variables_map vm);
    
    void WriteEnsembleTest(po::variables_map vm, std::vector<float> vBinning);
};
