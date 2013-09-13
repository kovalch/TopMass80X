#include "MVAAnalyzer.h"

#include "TCanvas.h"

#include <boost/lexical_cast.hpp>

void MVAAnalyzer::Analyze(const std::string& cuts, int i, int j) {
  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();

  gaus = new TF1("gaus", "gaus");

  fTree_->Draw("hadTopMass", cuts.c_str());

  fTree_->Fit("gaus", "hadTopMass", cuts.c_str());
  
  std::string path("plot/MVA/"); path+= fIdentifier_; path += std::string("_"); path += boost::lexical_cast<std::string>(i); path += std::string("_"); path += boost::lexical_cast<std::string>(j); path += std::string(".png");
  ctemp->Print(path.c_str());
        
  SetValue("mass", gaus->GetParameter(1), gaus->GetParError (1));
}
