#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "THStack.h"

void expand(TString sFile)
{
  // ---
  //    open input file
  // ---
  TFile* file = new TFile(sFile, "UPDATE");
  
  // ---
  //    Get hist
  // ---
  pileup = (TH1F*) file->Get("pileup");
  TH1F* pileup71 = new TH1F("pileup71", "data-derived pileup distrubution", 71, -0.5, 70.5);
  
  for (int i=0; i<71; i++) {
    pileup71->SetBinContent(i, pileup->GetBinContent(i));
  }
  
  pileup71->Draw();
  
  file->Write();
}
