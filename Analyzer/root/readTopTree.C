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

  TChain* chain = new TChain("analyzeHitFit/eventTree");
  std::cout << "added " << chain->Add("/nfs/dust/cms/user/mseidel/trees/Summer12_TTJets1725_1.00_muon/job_*_analyzeTop.root") << " files" << std::endl;

  Long64_t nentries = chain->GetEntries();

  // create a pointer to an event object. This will be used
  // to read the branch values.
  JetEvent        *jet            = new JetEvent();
  TopEvent        *top            = new TopEvent();
  WeightEvent     *weight         = new WeightEvent();
  
  chain->SetBranchStatus("*",1);

  chain->GetEntry(5); //needed to have full access to TChain; otherwise chain->GetBranch does not work

  TBranch *branchjet = chain->GetBranch("jet.");
  branchjet->SetAddress(&jet);
  TBranch *branchtop = chain->GetBranch("top.");
  branchtop->SetAddress(&top);
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
