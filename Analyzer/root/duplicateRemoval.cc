#include "duplicateRemoval.h"

//#include "TArrow.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
//#include "TF2.h"
#include "TFile.h"
//#include "TGraphErrors.h"
//#include "TGraph2DErrors.h"
//#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
//#include "TPaveText.h"
//#include "TMath.h"
//#include "TRandom3.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
//#include "TSystem.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TEventList.h"
#include "TTreeFormula.h"

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
  TFile* fromFile = new TFile("/scratch/hh/dust/naf/cms/user/mseidel/trees/Run2012_electron/job_Run2012D_analyzeTop.root", "READ");
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
  
  std::set<std::pair<int, int> > runevt;
  int nEntries = fromTree->GetEntries();
  int nDuplicates = 0;
  
  TFile* toFile = new TFile("/scratch/hh/dust/naf/cms/user/mseidel/trees/Run2012_electron/job_Run2012Dv2_analyzeTop.root", "RECREATE");
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
  fromFile->Close();
}


