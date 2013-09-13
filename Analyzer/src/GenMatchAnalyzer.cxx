#include "GenMatchAnalyzer.h"

#include "TCanvas.h"

#include <boost/lexical_cast.hpp>

void GenMatchAnalyzer::Analyze(const std::string& cuts, int i, int j) {
  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();

  TF1* voigt = new TF1("voigt", "[0]*TMath::Voigt(x-[1], 12*1.1, 2)");
  voigt->SetParameters(1000, 170);

  fTree_->Fit("voigt", "hadTopMass", cuts.c_str(), "LEM");
  
  std::string path("plot/GenMatch/"); path+= fIdentifier_; path += std::string("_"); path += boost::lexical_cast<std::string>(i); path += std::string("_"); path += boost::lexical_cast<std::string>(j); path += std::string(".png");
  ctemp->Print(path.c_str());
        
  SetValue("mass", voigt->GetParameter(1), voigt->GetParError (1));
}
