#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include <Math/Boost.h>

void testL7AllJets() {
  TFile* file = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_1.00_muon/analyzeTop.root", "READ");
  TTree* eventTree = (TTree*) file->Get("analyzeGenMatch/eventTree");
  
  int event; // event number
  eventTree->SetBranchAddress("event", &event);
  
  int combi; // permutation number inside event
  eventTree->SetBranchAddress("combi", &combi);
  
  int target; // 1 = correct permutation, 0 = wrong, -10 = unmatched
              // algorithm unambiguousOnly, maxDist = 0.3
  eventTree->SetBranchAddress("target", &target);
  
  // b-discriminator, medium working point = 0.679
  // preselection: maxBDiscLightJets = 0.85
  //               minBDiscBJets     = 0.50
  double hadQBCSV, hadQBarBCSV, hadBBCSV, lepBBCSV;
  eventTree->SetBranchAddress("hadQBCSV", &hadQBCSV);
  eventTree->SetBranchAddress("hadQBarBCSV", &hadQBarBCSV);
  eventTree->SetBranchAddress("hadBBCSV", &hadBBCSV);
  eventTree->SetBranchAddress("lepBBCSV", &lepBBCSV);
  
  // goodness-of-fit from HitFit
  double hitFitProb;
  eventTree->SetBranchAddress("hitFitProb", &hitFitProb);
  
  /* Objects after kinematic fit (used 4 leading jets for quarks)
     hadTop, hadB, hadW, hadQ, hadQBar, lepTop, lepB, lepW, lepton, nu
  */
  
  /* Objects before kinematic fit (used 4 leading jets for quarks)
     hadTopRaw, hadBRaw, hadWRaw, hadQRaw, hadQBarRaw, lepTopRaw, lepBRaw, lepWRaw, leptonRaw, nuRaw
  */
  TLorentzVector* hadWRaw = new TLorentzVector(); 
  eventTree->SetBranchAddress("hadWRaw.", &hadWRaw);
  TLorentzVector* hadQRaw = new TLorentzVector(); 
  eventTree->SetBranchAddress("hadQRaw.", &hadQRaw);
  TLorentzVector* hadQBarRaw = new TLorentzVector(); 
  eventTree->SetBranchAddress("hadQBarRaw.", &hadQBarRaw);
  
  /* MC truth
     hadTopGen, hadBGen, hadWGen, hadQGen, hadQBarGen, lepTopGen, lepBGen, lepWGen, leptonGen, nuGen
  */
  
  // Initialize histogram
  TH1F* hTest = new TH1F("hTest", "hTest", 100, 0, 200);
  TH1F* hCorr = new TH1F("hCorr", "hCorr", 100, 0, 200);
  hCorr->SetLineColor(kRed+1);
  
  double mW = 80.4;
  double p0 = 90.9183;
  double p1 = 4.21458e-2;
  double p2 = -1.62001;
  double xi = 1.01988005;
  double jecFactor1, jecFactor2;
  
  // Loop over all permutations
  for (int i = 0; i < eventTree->GetEntries(); i++) {
    eventTree->GetEntry(i);
    
    // Use only permutations passing the CSVM working point
    if (hadQBCSV > 0.679 && hadQBarBCSV > 0.679 && hadBBCSV < 0.679 && lepBBCSV < 0.679) continue;
    
    // Print
    
    jecFactor1 = ((mW*xi)/(p0 + p1*hadQRaw->Pt() + p2*log(hadQRaw->Pt())));
    jecFactor2 = ((mW*xi)/(p0 + p1*hadQBarRaw->Pt() + p2*log(hadQBarRaw->Pt())));
    
    // Fill histogram
    hTest->Fill(sqrt(hadWRaw->M2()));
    hCorr->Fill(sqrt(hadWRaw->M2()*jecFactor1*jecFactor2));
  }
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 500, 500);
  canvas->cd();
  hCorr->Draw("");
  hTest->Draw("SAME");
}
