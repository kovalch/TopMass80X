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

#include "iostream"

#include <boost/algorithm/string/replace.hpp>
//#include <boost/algorithm/string/split.hpp>
//#include <boost/algorithm/string/classification.hpp>

typedef ProgramOptionsReader po;

Skimmer::Skimmer()
{
  std::string samplePath(po::GetOption<std::string>("analysisConfig.samplePath"));
  //std::string samplePath("/nfs/dust/cms/user/eschliec/TopMass/2012/04/");
  //std::string samplePath("dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/user/eschliec/TopMassTreeWriter_04_DataMix01/");
  std::string outputPath("/nfs/dust/cms/user/eschliec/TopMass/2012/Skim_04/");

  std::string inputFolder = po::GetOption<std::string>("input");
  std::string inputSample = po::GetOption<std::string>("outPath");

  samplePath += inputFolder + std::string("/");

  std::vector<std::string> samples {
    // test
    //"*.root"
    // test sample
    //"Z2_S12_MadSpin_triggerTest_ABS_JES_100_172_5_sig.root"
    // data samples
    //"MJP12*_v1_data.root"
    //// background samples
    //,"QCDMixing_MJPS12*_v1_data.root"
    //// mass samples (including central sample)
    //"Z2_S12_ABS_JES_096_166_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_096_169_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_096_171_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_096_172_5_MadSpin_sig*.root"
    //,"Z2_S12_ABS_JES_096_173_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_096_175_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_096_178_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_098_166_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_098_169_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_098_171_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_098_172_5_MadSpin_sig*.root"
    //,"Z2_S12_ABS_JES_098_173_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_098_175_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_098_178_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_100_166_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_100_169_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_100_171_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_100_172_5_MadSpin_sig*.root"
    //,"Z2_S12_ABS_JES_100_173_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_100_175_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_100_178_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_102_166_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_102_169_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_102_171_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_102_172_5_MadSpin_sig*.root"
    //,"Z2_S12_ABS_JES_102_173_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_102_175_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_102_178_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_104_166_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_104_169_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_104_171_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_104_172_5_MadSpin_sig*.root"
    //,"Z2_S12_ABS_JES_104_173_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_104_175_5_MadSpin_sig.root"
    //,"Z2_S12_ABS_JES_104_178_5_MadSpin_sig.root"
    //// systematic samples from central sample
    //,"Z2_S12_JER_Down_MadSpin_sig.root"
    //,"Z2_S12_JER_Up_MadSpin_sig.root"
    //,"Z2_S12_BJES_Down_MadSpin_sig.root"
    //,"Z2_S12_BJES_Up_MadSpin_sig.root"
    //,"Z2_S12_CORFLAVQUARKJES_Down_MadSpin_sig.root"
    //,"Z2_S12_CORFLAVQUARKJES_Up_MadSpin_sig.root"
    //,"Z2_S12_CORTOTNOFLAVJES_Down_MadSpin_sig.root"
    //,"Z2_S12_CORTOTNOFLAVJES_Up_MadSpin_sig.root"
    //,"Z2_S12_COR1JES_Down_MadSpin_sig.root"
    //,"Z2_S12_COR1JES_Up_MadSpin_sig.root"
    //,"Z2_S12_COR2JES_Down_MadSpin_sig.root"
    //,"Z2_S12_COR2JES_Up_MadSpin_sig.root"
    //,"Z2_S12_COR3JES_Down_MadSpin_sig.root"
    //,"Z2_S12_COR3JES_Up_MadSpin_sig.root"
    //,"Z2_S12_COR4JES_Down_MadSpin_sig.root"
    //,"Z2_S12_COR4JES_Up_MadSpin_sig.root"
    //// systematic samples individual
    //,"Z2_S12_Matching_Down_MadSpin_sig.root"
    //,"Z2_S12_Matching_Up_MadSpin_sig.root"
    //,"Z2_S12_Scale_Down_MadSpin_sig.root"
    //,"Z2_S12_Scale_Up_MadSpin_sig.root"
    //// samples for CR and UE
    //,"Z2_S12*_ABS_JES_100_172_5_sig.root"
    //,"Z2_S12*_P11NoCR_sig.root"
    //,"Z2_S12*_P11TeV_sig.root"
    //,"Z2_S12*_P11_sig.root"
    //,"Z2_S12*_P11mpiHi_sig.root"
    //// generator samples
    //,"Z2_S12_MCNLO*_sig.root"
    //,"Z2_S12_POWHEG_sig.root"
    //,"Z2_S12_POWHER_sig.root"
    //,"Z2_S12_ABS_JES_100_172_5_MassiveBinDecay_sig.root"
    //// QCD samples
    //,"QCD_Alp.root"
    //,"QCD_Mad_Pt0100To0250.root"
    //,"QCD_Mad_Pt0250To0500.root"
    //,"QCD_Mad_Pt0500To1000.root"
    //,"QCD_Mad_Pt1000ToInfi.root"
  };

  samples.push_back(inputSample);

  for(auto& mySample : samples)
    skim(samplePath, outputPath, mySample);
}

void Skimmer::skim(std::string inputPath, std::string outputPath, std::string sample)
{
  std::string selection(po::GetOption<std::string>("analysisConfig.selection"));
  int maxPermutations(po::GetOption<int>("analysisConfig.maxPermutations"));
  //int maxPermutations(1);

  std::cout << "Selection: " << selection << std::endl;
  std::cout << "Max Permutations: " << maxPermutations << std::endl;

  TChain* fromChain = new TChain("analyzeKinFit/eventTree");
  int nFiles = 0;
  std::string fullInPath = inputPath+sample;
  std::cout << "Sample: " << fullInPath << std::endl;

  std::string prunedSampleName = po::GetOption<std::string>("input")+sample;
  boost::replace_all(prunedSampleName,"*","");
  std::string fullOutPath = outputPath+prunedSampleName;
  std::cout << "Output: " << fullOutPath << std::endl;

  nFiles += fromChain->Add(fullInPath.c_str());
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
  
  TFile* toFile = new TFile(fullOutPath.c_str(), "CREATE");
  if(toFile->IsZombie()) return;
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
    int iter = -1;
    for(int j = 0, l = sel->GetNdata(); j < l; ++j){
      if(sel->EvalInstance(j)) {
        ++iter;

        topEvent->recoTTBar  .at(iter) = topEvent->recoTTBar  .at(j);
        topEvent->recoTop1   .at(iter) = topEvent->recoTop1   .at(j);
        topEvent->recoTop2   .at(iter) = topEvent->recoTop2   .at(j);
        topEvent->recoW1     .at(iter) = topEvent->recoW1     .at(j);
        topEvent->recoW2     .at(iter) = topEvent->recoW2     .at(j);
        topEvent->recoB1     .at(iter) = topEvent->recoB1     .at(j);
        topEvent->recoW1Prod1.at(iter) = topEvent->recoW1Prod1.at(j);
        topEvent->recoW1Prod2.at(iter) = topEvent->recoW1Prod2.at(j);
        topEvent->recoB2     .at(iter) = topEvent->recoB2     .at(j);
        topEvent->recoW2Prod1.at(iter) = topEvent->recoW2Prod1.at(j);
        topEvent->recoW2Prod2.at(iter) = topEvent->recoW2Prod2.at(j);

        topEvent->recoJetIdxB1     .at(iter) = topEvent->recoJetIdxB1     .at(j);
        topEvent->recoJetIdxW1Prod1.at(iter) = topEvent->recoJetIdxW1Prod1.at(j);
        topEvent->recoJetIdxW1Prod2.at(iter) = topEvent->recoJetIdxW1Prod2.at(j);
        topEvent->recoJetIdxB2     .at(iter) = topEvent->recoJetIdxB2     .at(j);
        topEvent->recoJetIdxW2Prod1.at(iter) = topEvent->recoJetIdxW2Prod1.at(j);
        topEvent->recoJetIdxW2Prod2.at(iter) = topEvent->recoJetIdxW2Prod2.at(j);

        topEvent->fitTTBar  .at(iter) = topEvent->fitTTBar  .at(j);
        topEvent->fitTop1   .at(iter) = topEvent->fitTop1   .at(j);
        topEvent->fitTop2   .at(iter) = topEvent->fitTop2   .at(j);
        topEvent->fitW1     .at(iter) = topEvent->fitW1     .at(j);
        topEvent->fitW2     .at(iter) = topEvent->fitW2     .at(j);
        topEvent->fitB1     .at(iter) = topEvent->fitB1     .at(j);
        topEvent->fitW1Prod1.at(iter) = topEvent->fitW1Prod1.at(j);
        topEvent->fitW1Prod2.at(iter) = topEvent->fitW1Prod2.at(j);
        topEvent->fitB2     .at(iter) = topEvent->fitB2     .at(j);
        topEvent->fitW2Prod1.at(iter) = topEvent->fitW2Prod1.at(j);
        topEvent->fitW2Prod2.at(iter) = topEvent->fitW2Prod2.at(j);

        topEvent->combinationType.at(iter) = topEvent->combinationType.at(j);

        topEvent->fitProb.at(iter) = topEvent->fitProb.at(j);
        topEvent->fitChi2.at(iter) = topEvent->fitChi2.at(j);
        // treat sigMT specially as it is not used in all-jets channel up to now
        if(topEvent->fitSigMT.size() ==  topEvent->fitProb.size()) topEvent->fitSigMT.at(iter) = topEvent->fitSigMT.at(j);
      }
    }
    if(iter>=0){
      topEvent->shrink(maxPermutations);
      toTree->Fill();
    }
  }

  toFile->Write();
  toFile->Close();

  std::cout << "Finished shrinking tree!" << std::endl;
}


