#include "infRemoval.h"

#include "TopMass/TopEventTree/interface/JetEvent.h"
#include "TopMass/TopEventTree/interface/TopEvent.h"
#include "TopMass/TopEventTree/interface/WeightEvent.h"

#include "ProgramOptionsReader.h"

#include "TTree.h"
#include "TChain.h"
#include "TH2F.h"
#include "TFile.h"

#include "iostream"
#include "fstream"
#include "sstream"
#include <sys/stat.h>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

typedef ProgramOptionsReader po;

InfRemover::InfRemover()
{
  removeInfinities();
  
  //mkdir("backup", 0755);
  //rename("test2", "backup/test");
}

void InfRemover::removeInfinities()
{
  TFile* fromFile = new TFile(po::GetOption<std::string>("input").c_str(), "READ");
  TTree* fromTree = (TTree*) fromFile->Get("analyzeHitFit/eventTree");
  
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
  
  int nEntries   = fromTree->GetEntries("");
  int infEntries = fromTree->GetEntries("weight.combinedWeight>9999.");
  int nanEntries = fromTree->GetEntries("weight.combinedWeight!=weight.combinedWeight");
  std::cout << nEntries << " entries, " << infEntries << " infinities, " << nanEntries << " nan" << std::endl;
  if (infEntries+nanEntries == 0) return;
  
  int lastindex       = po::GetOption<std::string>("input").find_last_of("."); 
  int lastdir         = po::GetOption<std::string>("input").find_last_of("/");
  
  std::string rawname = po::GetOption<std::string>("input").substr(lastdir+1, lastindex-lastdir-1); 
  std::string dirname = "";
  if (lastdir > -1) dirname = po::GetOption<std::string>("input").substr(0, lastdir+1);
  std::string bakname = dirname + std::string("backup/") + rawname + std::string(".root");
  
  TFile* toFile = new TFile((dirname+rawname+std::string("_noInf.root")).c_str(), "RECREATE");
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
    if (weightEvent->combinedWeight > 9999.
     || weightEvent->combinedWeight != weightEvent->combinedWeight) {
        weightEvent->combinedWeight = 0.;
    }
    
    toTree->Fill();
  }
    
  toFile->Write();
  toFile->Close();
  
  mkdir((dirname+std::string("backup")).c_str(), 0755);
  rename(po::GetOption<std::string>("input").c_str(), bakname.c_str());
}


