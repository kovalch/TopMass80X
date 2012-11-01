#include "GenMatchAnalyzer.h"

#include "TCanvas.h"

void GenMatchAnalyzer::Analyze(const TString& cuts, int i, int j) {
  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();

  TF1* voigt = new TF1("voigt", "[0]*TMath::Voigt(x-[1], 12*1.1, 2)");
  voigt->SetParameters(1000, 170);

  fTree_->Fit("voigt", "hadTopMass", cuts, "LEM");
  
  TString path("plot/GenMatch/"); path+= fIdentifier_; path += "_"; path += i; path += "_"; path += j; path += ".png";
  ctemp->Print(path);
        
  SetValue("mass", voigt->GetParameter(1), voigt->GetParError (1));
}
