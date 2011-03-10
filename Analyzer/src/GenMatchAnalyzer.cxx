#include "GenMatchAnalyzer.h"

void GenMatchAnalyzer::Analyze(TString cuts) {
  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();

  gaus = new TF1("gaus", "gaus");

  fChain->Draw("hadTopMass", cuts);

  fChain->Fit("gaus", "hadTopMass", cuts);

  ctemp->Print("plot/temp.ps");
        
  fMass      = gaus->GetParameter(1);
  fMassError = gaus->GetParError(1);
  fMassSigma = gaus->GetParameter(2);
}
