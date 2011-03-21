#include "MVAAnalyzer.h"

double MVAAnalyzer::GetMass() {
  return fMass;
}

void MVAAnalyzer::Analyze(TString cuts, int i, int j) {
  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();

  gaus = new TF1("gaus", "gaus");

  fTree->Draw("hadTopMass", cuts);

  fTree->Fit("gaus", "hadTopMass", cuts);
  
  TString path("plot/MVA/"); path+= fIdentifier; path += "_"; path += i; path += "_"; path += j; path += ".png";
  ctemp->Print(path);
        
  fMass      = gaus->GetParameter(1);
  fMassError = gaus->GetParError(1);
  fMassSigma = gaus->GetParameter(2);
}
