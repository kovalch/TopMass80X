#include "skimmer.h"

//#include "TArrow.h"
//#include "TAxis.h"
//#include "TCanvas.h"
#include "TChain.h"
//#include "TF1.h"
//#include "TF2.h"
#include "TFile.h"
//#include "TGraphErrors.h"
//#include "TGraph2DErrors.h"
//#include "TH1F.h"
//#include "TLatex.h"
//#include "TLegend.h"
//#include "TPaveText.h"
//#include "TMath.h"
//#include "TRandom3.h"
//#include "TROOT.h"
//#include "TString.h"
//#include "TStyle.h"
//#include "TSystem.h"
#include "TTreeFormula.h"

#include "TopMass/TopEventTree/interface/JetEvent.h"
#include "TopMass/TopEventTree/interface/TopEvent.h"
#include "TopMass/TopEventTree/interface/WeightEvent.h"

#include "ProgramOptionsReader.h"

#include "string"
#include "iostream"

//#include <boost/algorithm/string/split.hpp>
//#include <boost/algorithm/string/classification.hpp>

typedef ProgramOptionsReader po;

Skimmer::Skimmer()
{
  skim();
}

void Skimmer::skim()
{
  //std::string samplePath(po::GetOption<std::string>("analysisConfig.samplePath"));
  //std::string samplePath("/scratch/hh/dust/naf/cms/user/eschliec/GRID-CONTROL_JOBS/TopMassTreeWriter_02_Data06/");
  std::string samplePath("dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/user/eschliec/TopMassTreeWriter_02_Data06/");
  std::string selection(po::GetOption<std::string>("analysisConfig.selection"));
  //int maxPermutations(po::GetOption<int>("analysisConfig.maxPermutations"));
  int maxPermutations(1);
  std::string input(po::GetOption<std::string>("input"));

  TChain* fromChain = new TChain("analyzeKinFit/eventTree");
  int nFiles = 0;
  nFiles += fromChain->Add((samplePath+std::string("QCDMixing_MJPS12B_v1_data/*" )+input+std::string(".root")).c_str());
  nFiles += fromChain->Add((samplePath+std::string("QCDMixing_MJPS12C1_v1_data/*")+input+std::string(".root")).c_str());
  nFiles += fromChain->Add((samplePath+std::string("QCDMixing_MJPS12C2_v1_data/*")+input+std::string(".root")).c_str());
  nFiles += fromChain->Add((samplePath+std::string("QCDMixing_MJPS12D1_v1_data/*")+input+std::string(".root")).c_str());

  std::cout << "Number of files to be skimmed: " << nFiles << std::endl;

  fromChain->SetBranchStatus("*", 0);
  fromChain->SetBranchStatus("jet.*"   , 1);
  fromChain->SetBranchStatus("top.*"   , 1);
  fromChain->SetBranchStatus("weight.*", 1);

  JetEvent    *jetEvent    = new JetEvent();
  TopEvent    *topEvent    = new TopEvent();
  WeightEvent *weightEvent = new WeightEvent();

  fromChain->SetBranchAddress("jet."   , &jetEvent);
  fromChain->SetBranchAddress("top."   , &topEvent);
  fromChain->SetBranchAddress("weight.", &weightEvent);
  
  TFile* toFile = new TFile((std::string("QCDMixing_MJPS12_v1_data_dCache_")+input+std::string(".root")).c_str(), "RECREATE");
  toFile->mkdir("analyzeKinFit")->cd();
  
  TTree* toTree = fromChain->CloneTree(0);
  toTree->SetBranchAddress("jet."   , &jetEvent);
  toTree->SetBranchAddress("top."   , &topEvent);
  toTree->SetBranchAddress("weight.", &weightEvent);
  
  std::cout << "Starting to shrink tree ..." << std::endl;

  TTreeFormula *sel = new TTreeFormula("sel", selection.c_str(), fromChain);

  //int nEntries = fromChain->GetEntries();
  int numberOfSkimmedFiles = -1;

  for(int i = 0; ; ++i){
    //if(i%1000000==0) std::cout << i << " / " << nEntries << " (" << i*100/nEntries << "%)" << std::endl;
    long entry = fromChain->LoadTree(i);
    if(entry  < 0) break;
    if(entry == 0){
      sel->UpdateFormulaLeaves();
      std::cout << ++numberOfSkimmedFiles << " / " << nFiles << " skimmed ..." << std::endl;
    }
    fromChain->GetTree()->GetEntry(entry, 1);
    if(!sel->GetNdata()) continue;
    int shrink = -1;
    for(int j = 0, l = std::min(maxPermutations, sel->GetNdata()); j < l; ++j){
      if(sel->EvalInstance(j)) shrink = j;
    }
    if(shrink>=0){
      topEvent->shrink(shrink+1);
      toTree->Fill();
    }
  }

  std::cout << "Finished shrinking tree!" << std::endl;
  
  toFile->Write();
  toFile->Close();
}


