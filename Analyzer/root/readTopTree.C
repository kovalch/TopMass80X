#include <iostream>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "TClassTable.h"
#include "TSystem.h"

#include <vector>

// Header file for the classes stored in the TTree if any.
#include <TObject.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TVector2.h>

#include <TopMass/TopEventTree/interface/JetEvent.h>
#include <TopMass/TopEventTree/interface/BRegJetEvent.h>
#include <TopMass/TopEventTree/interface/TopEvent.h>
#include <TopMass/TopEventTree/interface/WeightEvent.h>

// old C++ standard with ACliC so far... ?!

void readTopTree(){

//  TFile *f = new TFile("/nfs/dust/test/cms/user/kirschen/BRegression_MCS3250MinEvtsBReg/Summer12_TTJets1725_1.00_muon/merged.root");
//  TTree *chain = (TTree*)f->Get("analyzeHitFit/eventTree");


// In principle, this should also work for TChains... 
// solution:   Long64_t nentries = chain->GetEntries();   chain->GetEntry(5); before chain->GetBranch(X);
  TChain* chain = new TChain("analyzeHitFit/eventTree");
  //  std::cout << "added " << chain->Add("/nfs/dust/test/cms/user/kirschen/BRegression_MCS3250MinEvtsBReg/Summer12_TTJets1725_1.00_muon/job_*_analyzeTop.root") << " files" << std::endl;
  std::cout << "added " << chain->Add("/nfs/dust/cms/user/mseidel/trees/Summer12_TTJets1725_1.00_muon/job_*_analyzeTop.root") << " files" << std::endl;

  Long64_t nentries = chain->GetEntries();

  // create a pointer to an event object. This will be used
  // to read the branch values.
  JetEvent        *jet            = new JetEvent();
  BRegJetEvent    *BRegJet        = new BRegJetEvent(); //not available in all trees
  TopEvent        *top            = new TopEvent();
  TopEvent        *bRegTop        = new TopEvent(); //not available in all trees
  WeightEvent     *weight         = new WeightEvent();
  
  chain->SetBranchStatus("*",1);

  chain->GetEntry(5); //needed to have full access to TChain

  //  TBranch *branchBRegJet = chain->GetBranch("BRegJet.");
  //  branchBRegJet->SetAddress(&BRegJet);
  TBranch *branchjet = chain->GetBranch("jet.");
  branchjet->SetAddress(&jet);
  TBranch *branchtop = chain->GetBranch("top.");
  branchtop->SetAddress(&top);
  //  TBranch *branchbregtop = chain->GetBranch("bRegTop.");
  //  branchbregtop->SetAddress(&bRegTop);
  TBranch *branchweight = chain->GetBranch("weight.");
  branchweight->SetAddress(&weight);


  //  comment out to run over a small number of events
  nentries = 10;

  std::cout <<"nentries: " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    chain->GetEntry(jentry);
    if(jentry%20000==0)std::cout << "event " << jentry << std::endl;

    //access to all members of  jet, BRegJet, top bRegTop, weight
    //example: std::vector <TLorentzVector> fitTop1 -->
    for(unsigned int idx = 0; idx < top->fitTop1.size(); idx++){
      std::cout << "fitTop1.at(idx).M() at idx "<< idx << ": " << top->fitTop1.at(idx).M() << std::endl;
    }

  }

}
