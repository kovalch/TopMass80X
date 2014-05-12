#include "duplicateRemoval.h"

#include "TopMass/TopEventTree/interface/JetEvent.h"
#include "TopMass/TopEventTree/interface/TopEvent.h"
#include "TopMass/TopEventTree/interface/WeightEvent.h"

#include "ProgramOptionsReader.h"

#include "iostream"
#include "fstream"
#include "sstream"

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

typedef ProgramOptionsReader po;

DuplicateRemover::DuplicateRemover()
{
  removeDuplicates();
}

void DuplicateRemover::removeDuplicates()
{
  TChain* fromTree = new TChain("analyzeHitFit/eventTree");
  fromTree->Add(po::GetOption<std::string>("input"));
  
  fromTree->SetBranchStatus("*", 0);
  fromTree->SetBranchStatus("jet.*"   , 1);
  fromTree->SetBranchStatus("top.*"   , 1);
  fromTree->SetBranchStatus("weight.*", 1);

  JetEvent    *jetEvent    = new JetEvent();
  TopEvent    *topEvent    = new TopEvent();
  WeightEvent *weightEvent = new WeightEvent();

  fromTree->SetBranchAddress("jet."   , &jetEvent);
  fromTree->SetBranchAddress("top."   , &topEvent);
  fromTree->SetBranchAddress("weight.", &weightEvent);
  
  std::set<std::pair<int, int> > runevt;
  int nEntries = fromTree->GetEntries();
  int nDuplicates = 0;
  
  TFile* toFile = new TFile(po::GetOption<std::string>("outPath"), "RECREATE");
  toFile->mkdir("analyzeHitFit");
  toFile->cd("analyzeHitFit");
  
  TTree* toTree = fromTree->CloneTree(0);
  toTree->SetBranchAddress("jet."   , &jetEvent);
  toTree->SetBranchAddress("top."   , &topEvent);
  toTree->SetBranchAddress("weight.", &weightEvent);
  
  for (int i = 0; i < nEntries; i++) {
    fromTree->GetEntry(i);
    //std::cout << i << " " << topEvent->run << " " << topEvent->event << std::endl;

    // if we have not seen this event yet, add it to the set
    // and to the entry list
    if (runevt.count(std::pair<int, int>(topEvent->run, topEvent->event)) == 0)
    {
        runevt.insert(std::pair<int, int>(topEvent->run, topEvent->event));
        toTree->Fill();
    }
    else {
      ++nDuplicates;
      //std::cout << "Duplicate found: " << run << " : " << event << std::endl;
    }
  }
  
  std::cout << "Total number of events: " << nEntries << "   Duplicates: " << nDuplicates << std::endl;
  
  toFile->Write();
  toFile->Close();
}


