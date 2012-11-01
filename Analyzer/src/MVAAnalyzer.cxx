#include "MVAAnalyzer.h"

#include "TCanvas.h"

void MVAAnalyzer::Analyze(const TString& cuts, int i, int j) {
  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();

  gaus = new TF1("gaus", "gaus");

  fTree_->Draw("hadTopMass", cuts);

  fTree_->Fit("gaus", "hadTopMass", cuts);
  
  TString path("plot/MVA/"); path+= fIdentifier_; path += "_"; path += i; path += "_"; path += j; path += ".png";
  ctemp->Print(path);
        
  SetValue("mass", gaus->GetParameter(1), gaus->GetParError (1));
}
