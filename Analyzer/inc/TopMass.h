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

struct massPoint {
  double genMass;
  double genJES;
  double genLumi;
  TString identifier;
  TString fileName;
  Analysis* analysis;
  
  TH1D* hMass;
  TH1D* hMassError;
  TH1D* hMassPull;
  TH2F* h2Mass;
  TH2F* h2MassError;
  TH3F* h3Mass;
  TH3F* h3MassError;
  TH3F* h3MassPull;
  TH3F* h3JES;
  TH3F* h3JESError;
  TH3F* h3JESPull;
};

class TopMass {
  private:
    TString fMethod;
    int fBins;
    double fLumi;
    
    Analysis* aSim;
      
  public:
    TopMass(po::variables_map vm);
    
    void WriteEnsembleTest(po::variables_map vm, std::vector<float> vBinning);
};
