#include "GenMatchAnalyzer.h"

double GenMatchAnalyzer::GetMass() {
  return fMass;
}

void GenMatchAnalyzer::Analyze(TString cuts, int i, int j) {
  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();

  gaus = new TF1("gaus", "gaus");

  fTree->Draw("hadTopMass", cuts);

  fTree->Fit("gaus", "hadTopMass", cuts);
  
  if (gaus->GetChisquare()/gaus->GetNDF() > 3) {
    gaus = new TF1("gaus", "gaus(0)+gaus(3)");
    gaus->SetParameters(100, 170, 10, 10, 250, 50);
    fTree->Fit("gaus", "hadTopMass", cuts);
  }
  
  TString path("plot/GenMatch/"); path+= fIdentifier; path += "_"; path += i; path += "_"; path += j; path += ".png";
  ctemp->Print(path);
        
  fMass      = gaus->GetParameter(1);
  fMassError = gaus->GetParError(1);
  fMassSigma = gaus->GetParameter(2);
}
