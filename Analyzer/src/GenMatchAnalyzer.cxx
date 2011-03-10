#include "GenMatchAnalyzer.h"

void GenMatchAnalyzer::Analyze(TString cuts, int i, int j) {
  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();

  gaus = new TF1("gaus", "gaus");

  fChain->Draw("hadTopMass", cuts);

  fChain->Fit("gaus", "hadTopMass", cuts);
  
  TString path("plot/GenMatch/"); path += i; path += "_"; path += j; path += ".png";
  ctemp->Print(path);
        
  fMass      = gaus->GetParameter(1);
  fMassError = gaus->GetParError(1);
  fMassSigma = gaus->GetParameter(2);
}
