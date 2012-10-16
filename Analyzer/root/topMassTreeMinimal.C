#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"
#include <Math/Boost.h>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > PxPyPzEVector;

void topMassTreeMinimal() {
  TFile* file = new TFile("analyzeTop.root", "READ");
  TTree* eventTree = (TTree*) file->Get("analyzeHitFit/eventTree");
  
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
  PxPyPzEVector* hadTop    = new PxPyPzEVector(); 
  eventTree->SetBranchAddress("hadTop.", &hadTop);
  
  /* Objects before kinematic fit (used 4 leading jets for quarks)
     hadTopRaw, hadBRaw, hadWRaw, hadQRaw, hadQBarRaw, lepTopRaw, lepBRaw, lepWRaw, leptonRaw, nuRaw
  */
  PxPyPzEVector* hadTopRaw = new PxPyPzEVector(); 
  eventTree->SetBranchAddress("hadTopRaw.", &hadTopRaw);
  
  /* MC truth
     hadTopGen, hadBGen, hadWGen, hadQGen, hadQBarGen, lepTopGen, lepBGen, lepWGen, leptonGen, nuGen
  */
  PxPyPzEVector* hadTopGen = new PxPyPzEVector(); 
  eventTree->SetBranchAddress("hadTopGen.", &hadTopGen);
  
  // Initialize histogram
  TH1F* hTest = new TH1F("hTest", "hTest", 100, 0, 500);
  
  // Loop over all permutations
  for (int i = 0; i < eventTree->GetEntries(); i++) {
    eventTree->GetEntry(i);
    
    // Use only permutations passing the CSVM working point
    if (hadQBCSV < 0.679 && hadQBarBCSV < 0.679 && hadBBCSV > 0.679 && lepBBCSV > 0.679) continue;
    
    // Print
    std::cout << hadTop->M() << std::endl;
    
    // Fill histogram
    hTest->Fill(hadTop->M());
  }
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 500, 500);
  canvas->cd();
  hTest->Draw("");
}
