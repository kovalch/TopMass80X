#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"

#include "tdrstyle.C"

void massByJESBinned()
{
  // Open input files
  TFile* fTTJets097 = new TFile("/scratch/hh/current/cms/user/mseidel/TTJets1725-S4_1.00_0.97/analyzeTop.root");
  TFile* fTTJets100 = new TFile("/scratch/hh/current/cms/user/mseidel/TTJets1725-S4_1.00_1.00/analyzeTop.root");
  TFile* fTTJets103 = new TFile("/scratch/hh/current/cms/user/mseidel/TTJets1725-S4_1.00_1.03/analyzeTop.root");
  
  // Get trees
  TTree* t097 = (TTree*) fTTJets097->Get("analyzeHitFit/eventTree");
  TTree* t100 = (TTree*) fTTJets100->Get("analyzeHitFit/eventTree");
  TTree* t103 = (TTree*) fTTJets103->Get("analyzeHitFit/eventTree");
  
  // Observable setup  
  TString binObservable("genHadBE");
  double firstBinLowEdge = 30;
  double lastBinHighEdge = 230;
  double binWidth        = 20;
  int    nBins = (lastBinHighEdge-firstBinLowEdge)/binWidth;
  
  // Graph preps
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("bJES influence;genHadBE;mass");
  
  TVectorD mass(nBins);
  TVectorD massError(nBins);
  TVectorD bin(nBins);
  TVectorD binError(nBins);
  
  for (int iJES=0; iJES<3; iJES++) {
    for (int i=0; i<nBins; i++) {
      fit = new TF1("fit", "gaus");
      TString cut("target==1&"); cut+=binObservable; cut+=">"; cut+=firstBinLowEdge+binWidth*i; cut +="&"; cut+=binObservable; cut+="<"; cut+=firstBinLowEdge+binWidth*(i+1);
      
      switch (iJES) {
        case 0: {
          t097->Fit("fit", "hadTopMass", cut);
          break;
        }
        case 1: {
          t100->Fit("fit", "hadTopMass", cut);
          break;
        }
        case 2: {
          t103->Fit("fit", "hadTopMass", cut);
          break;
        }
      }
      
      mass[i]      = fit->GetParameter(1);
      massError[i] = fit->GetParError(1);
      bin[i]       = firstBinLowEdge+binWidth*(i+0.5);
      binError[i]  = sqrt(binWidth/2.);
    }
    
    TGraphErrors* gr = new TGraphErrors(bin, mass, binError, massError);
    gr->SetLineColor(kRed+iJES);
    
    mg->Add(gr);
  }
  
  c1->Clear();
  mg->Draw("AC");
}
