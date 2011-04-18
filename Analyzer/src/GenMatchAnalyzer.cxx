#include "GenMatchAnalyzer.h"

double GenMatchAnalyzer::GetMass() {
  return fMass;
}

void GenMatchAnalyzer::Analyze(TString cuts, int i, int j) {
  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();

  TF1* voigt = new TF1("voigt", "[0]*TMath::Voigt(x-[1], 12, 2)");
  voigt->SetParameters(1000, 170);

  fTree->Fit("voigt", "hadTopMass", cuts, "LEM");
  
  TString path("plot/GenMatch/"); path+= fIdentifier; path += "_"; path += i; path += "_"; path += j; path += ".png";
  ctemp->Print(path);
        
  fMass      = voigt->GetParameter(1);
  fMassError = voigt->GetParError(1);
  fMassSigma = voigt->GetParameter(2);
}
