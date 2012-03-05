#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "THStack.h"

#include "tdrstyle.C"

void subTreeEvent2b() {
  /*
  makeSubTreeEvent2b("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_1.00/");
  makeSubTreeEvent2b("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_scaleup/");
  makeSubTreeEvent2b("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_scaledown/");
  makeSubTreeEvent2b("/scratch/hh/current/cms/user/mseidel/Summer11_QCD/");
  makeSubTreeEvent2b("/scratch/hh/current/cms/user/mseidel/Summer11_WJets/");
  makeSubTreeEvent2b("/scratch/hh/current/cms/user/mseidel/Summer11_ZJets/");
  makeSubTreeEvent2b("/scratch/hh/current/cms/user/mseidel/Summer11_Tbar_s-channel/");
  makeSubTreeEvent2b("/scratch/hh/current/cms/user/mseidel/Summer11_Tbar_t-channel/");
  makeSubTreeEvent2b("/scratch/hh/current/cms/user/mseidel/Summer11_Tbar_tW-channel/");
  makeSubTreeEvent2b("/scratch/hh/current/cms/user/mseidel/Summer11_T_s-channel/");
  makeSubTreeEvent2b("/scratch/hh/current/cms/user/mseidel/Summer11_T_t-channel/");
  makeSubTreeEvent2b("/scratch/hh/current/cms/user/mseidel/Summer11_T_tW-channel/");
  //*/
  makeSubTreeEvent2b("/scratch/hh/current/cms/user/mseidel/Run2011_MuHad/");
  makeSubTreeEvent2b("/scratch/hh/current/cms/user/mseidel/Run2011/");
}

void makeSubTreeEvent2b(TString sDir) {
  std::cout << sDir << std::endl;
  
  TString sFromFile = sDir; sFromFile += "analyzeTop.root";
  TString sToFile = sDir; sToFile += "analyzeTop_event2bW.root";
  
  TFile* fromFile = new TFile(sFromFile, "READ");
  fromTree = (TTree*) fromFile->Get("analyzeMVADisc/eventTree");
  fromTree->SetBranchStatus("*", 0);
  fromTree->SetBranchStatus("run", 1);
  fromTree->SetBranchStatus("luminosityBlock", 1);
  fromTree->SetBranchStatus("event", 1);
  fromTree->SetBranchStatus("combi", 1);
  fromTree->SetBranchStatus("leptonPt", 1);
  fromTree->SetBranchStatus("hadQBSSV", 1);
  fromTree->SetBranchStatus("hadQBarBSSV", 1);
  fromTree->SetBranchStatus("hadBBSSV", 1);
  fromTree->SetBranchStatus("lepBBSSV", 1);
  fromTree->SetBranchStatus("bottomSSVJetMultiplicity", 1);
  fromTree->SetBranchStatus("jetMultiplicity", 1);
  fromTree->SetBranchStatus("nVertex", 1);
  fromTree->SetBranchStatus("leptonEta", 1);
  fromTree->SetBranchStatus("nuRawPt", 1);
  fromTree->SetBranchStatus("hadBRawPt", 1);
  fromTree->SetBranchStatus("hadQRawPt", 1);
  fromTree->SetBranchStatus("hadQBarRawPt", 1);
  fromTree->SetBranchStatus("hadBEta", 1);
  fromTree->SetBranchStatus("hadQEta", 1);
  fromTree->SetBranchStatus("hadWRawMass", 1);
  fromTree->SetBranchStatus("lepWRawMass", 1);
  fromTree->SetBranchStatus("hadTopRawMass", 1);
  fromTree->SetBranchStatus("lepTopRawMass", 1);
  fromTree->SetBranchStatus("hadTopMass", 1);
  fromTree->SetBranchStatus("hitFitProb", 1);
  fromTree->SetBranchStatus("hitFitChi2", 1);
  fromTree->SetBranchStatus("MCWeight", 1);
  fromTree->SetBranchStatus("target", 1);
  
  TFile* toFile = new TFile(sToFile, "RECREATE");
  TTree* toTree = fromTree->CloneTree(0);
  
  int event;
  fromTree->SetBranchAddress("event", &event);
  toTree->SetBranchAddress("event", &event);
  
  int combi;
  fromTree->SetBranchAddress("combi", &combi);
  toTree->SetBranchAddress("combi", &combi);
  
  double leptonPt;
  fromTree->SetBranchAddress("leptonPt", &leptonPt);
  toTree->SetBranchAddress("leptonPt", &leptonPt);
  
  double hadQBSSV;
  fromTree->SetBranchAddress("hadQBSSV", &hadQBSSV);
  toTree->SetBranchAddress("hadQBSSV", &hadQBSSV);
  
  double hadQBarBSSV;
  fromTree->SetBranchAddress("hadQBarBSSV", &hadQBarBSSV);
  toTree->SetBranchAddress("hadQBarBSSV", &hadQBarBSSV);
  
  double hadBBSSV;
  fromTree->SetBranchAddress("hadBBSSV", &hadBBSSV);
  toTree->SetBranchAddress("hadBBSSV", &hadBBSSV);
  
  double lepBBSSV;
  fromTree->SetBranchAddress("lepBBSSV", &lepBBSSV);
  toTree->SetBranchAddress("lepBBSSV", &lepBBSSV);
  
  double hitFitProb;
  fromTree->SetBranchAddress("hitFitProb", &hitFitProb);
  toTree->SetBranchAddress("hitFitProb", &hitFitProb);
  
  int entries = fromTree->GetEntries();  
  int currentEvent = -1;
  
  for (int i = 0; i < entries; i++) {
    if (i%100000 == 0) std::cout << entries - i << "... " << std::endl;
    
    fromTree->GetEntry(i);
    
    if (event!=currentEvent && leptonPt > 30 && hadQBSSV < 1.74 && hadQBarBSSV < 1.74 && hadBBSSV > 1.74 && lepBBSSV > 1.74 && hitFitProb > 0.2) {
      currentEvent = event;
      toTree->Fill();
    }
  }
  
  toFile->Write();
  toFile->Close();
  fromFile->Close();
}
