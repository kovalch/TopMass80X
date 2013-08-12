#include "RandomSubsetCreatorNewInterface.h"

#include "Analysis.h"
#include "Helper.h"
#include "ProgramOptionsReader.h"

#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/progress.hpp>

#include "TChain.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TSystem.h"

#include "LHAPDF/LHAPDF.h"

typedef ProgramOptionsReader po;

RandomSubsetCreatorNewInterface::RandomSubsetCreatorNewInterface() :
    // filled from program options
    selection_  (po::GetOption<std::string>("analysisConfig.selection")),
    samplePath_ (po::GetOption<std::string>("analysisConfig.samplePath")),
    fIdentifier_(po::GetOption<std::string>("input")),
    //fChannel_   (po::GetOption<std::string>("channel")),
    fWeight_    (po::GetOption<std::string>("weight")),
    fLumi_  (po::GetOption<double>("lumi")),
    fSig_   (po::GetOption<double>("fsig")),
    fBDisc_ (po::GetOption<double>("bdisc")),
    topEvent_(new TopEvent()),
    weightEvent_(new WeightEvent()),
    random_(0)
{
  channelID_ = Helper::channelID();

  if (channelID_ == Helper::kAllJets) {
    fTreesSig_.push_back(PrepareTree(samplePath_+fIdentifier_+TString(".root")));
    if(fLumi_>0) fTreesBkg_.push_back(PrepareTree(samplePath_+"QCDEstimationMix.root"));
  }
  if (channelID_ == Helper::kMuonJets || channelID_ == Helper::kLeptonJets) {
    fTreesSig_.push_back(PrepareTree(samplePath_+fIdentifier_+TString("_muon/analyzeTop.root")));
    if(fLumi_>0) {
      fTreesBkg_.push_back(PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_Wbb_muon/analyzeTop.root"));
      fTreesBkg_.push_back(PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_T_muon/analyzeTop.root"));
    }
    //fTreeTTmu = PrepareTree(samplePath_+fIdentifier_+TString("_muon/analyzeTop.root"));
    //fTreeWmu  = PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_Wbb_muon/analyzeTop.root");
    //fTreeSTmu = PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_T_muon/analyzeTop.root");
  }
  if (channelID_ == Helper::kElectronJets || channelID_ == Helper::kLeptonJets) {
    fTreesSig_.push_back(PrepareTree(samplePath_+fIdentifier_+TString("_electron/analyzeTop.root")));
    if(fLumi_>0) {
      fTreesBkg_.push_back(PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_Wbb_electron/analyzeTop.root"));
      fTreesBkg_.push_back(PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_T_electron/analyzeTop.root"));
    }
    //fTreeTTe  = PrepareTree(samplePath_+fIdentifier_+TString("_electron/analyzeTop.root"));
    //fTreeWe   = PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_Wbb_electron/analyzeTop.root");
    //fTreeSTe  = PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_T_electron/analyzeTop.root");
  }
  random_ = new TRandom3(0);
  std::cout << "Random seed: " << random_->GetSeed() << std::endl;
}

RandomSubsetCreatorNewInterface::~RandomSubsetCreatorNewInterface() {
}

TTree* RandomSubsetCreatorNewInterface::CreateRandomSubset() {
  if (fLumi_>0) {
    std::cout << "Create random subset..." << std::endl;

    time_t start, end;
    time(&start);
    time(&end);

    fTree_ = new TTree("fTree", "fTree");
    fTree_->Branch("top.", topEvent_);
    fTree_->Branch("weight.", weightEvent_);

    // DATA
    double nEventsDataAllJets  = 2767.;
    double nEventsDataMuon     = 2906.;
    double nEventsDataElectron = 2268.;

    int eventsPEAllJets  = random_->Poisson(nEventsDataAllJets /3544.844*fLumi_);
    int eventsPEMuon     = random_->Poisson(nEventsDataMuon    /5000.000*fLumi_);
    int eventsPEElectron = random_->Poisson(nEventsDataElectron/5000.000*fLumi_);

    if (channelID_ == Helper::kAllJets) {
      DrawEvents(fTreesSig_[0], eventsPEAllJets*    fSig_ );
      DrawEvents(fTreesBkg_[0], eventsPEAllJets*(1.-fSig_));
    }
    if (channelID_ == Helper::kMuonJets || channelID_ == Helper::kLeptonJets) {
      DrawEvents(fTreesSig_[0], eventsPEMuon*    fSig_       );
      DrawEvents(fTreesBkg_[0], eventsPEMuon*(1.-fSig_)*1./4.);
      DrawEvents(fTreesBkg_[1], eventsPEMuon*(1.-fSig_)*3./4.);
    }
    if (channelID_ == Helper::kElectronJets || channelID_ == Helper::kLeptonJets) {
      short offset = 0;
      if(channelID_ == Helper::kLeptonJets) offset = 1;
      DrawEvents(fTreesSig_[  offset+0], eventsPEElectron*    fSig_       );
      DrawEvents(fTreesBkg_[2*offset+0], eventsPEElectron*(1.-fSig_)*1./4.);
      DrawEvents(fTreesBkg_[2*offset+1], eventsPEElectron*(1.-fSig_)*3./4.);
    }

    time(&end);
    std::cout << "Created random subset in " << difftime(end, start) << " seconds." << std::endl;
  }
  else {
    TList treeList;

    treeList.Add(fTreesSig_[0]);
    if (channelID_ == Helper::kLeptonJets) treeList.Add(fTreesSig_[1]);

    fTree_ = TTree::MergeTrees(&treeList);
  }

  return fTree_;
}

void RandomSubsetCreatorNewInterface::DrawEvents(TTree* tempTree, double nEventsPE) {
  std::cout << "nEventsPE: " << nEventsPE << std::endl;

  fTree_->CopyAddresses(tempTree);

  int permsMC = tempTree->GetEntries("");

  std::vector<std::string> vWeight;
  boost::split( vWeight, fWeight_, boost::is_any_of("-*"));

  double maxMCWeight = 1.; // calculate upper bound for combined MCWeight

  for (unsigned int i = 0; i < vWeight.size(); ++i) {
    std::cout << vWeight[i] << std::endl;
    if (strncmp((TString) vWeight[i], "pdfWeights", 10)) maxMCWeight *= tempTree->GetMaximum((TString) vWeight[i]);
    else maxMCWeight *= tempTree->GetMaximum("pdfWeights");
  }

  std::cout << "maxMCWeight(" << fWeight_ << "):" << maxMCWeight  << std::endl;

  if (maxMCWeight ==  0) { std::cout << "Weight not active?" << std::endl; }
  if (maxMCWeight == -1) { std::cout << "Running over data?" << std::endl; }

  std::cout << "while (eventsDrawn < eventsPE)..." << std::endl;

  int eventsDrawn = 0;
  int nAttempts = 0;

  boost::progress_display progress((int)nEventsPE, std::cout);

  while (eventsDrawn < (int)nEventsPE) {
    int drawn = random_->Integer(permsMC);
    tempTree->GetEntry(drawn);
    ++nAttempts;
    double eventWeight = 1.;
    for (unsigned int i = 0; i < vWeight.size(); ++i) {
      if (!strcmp((TString) vWeight[i], "mcWeight")) eventWeight *= weightEvent_->mcWeight;
      else if (!strcmp((TString) vWeight[i], "muWeight")) eventWeight *= weightEvent_->muWeight;
      else if (!strcmp((TString) vWeight[i], "puWeight")) eventWeight *= weightEvent_->puWeight;
      else if (!strcmp((TString) vWeight[i], "puWeightUp")) eventWeight *= weightEvent_->puWeightUp;
      else if (!strcmp((TString) vWeight[i], "puWeightDown")) eventWeight *= weightEvent_->puWeightDown;
      //TODO still needs to be fixed in WeightEvent
      //else if (!strcmp((TString) vWeight[i], "bWeight")) eventWeight *= bWeight;
      //else if (!strcmp((TString) vWeight[i], "bWeight_bTagSFUp")) eventWeight *= bWeight_bTagSFUp;
      //else if (!strcmp((TString) vWeight[i], "bWeight_bTagSFDown")) eventWeight *= bWeight_bTagSFDown;
      //else if (!strcmp((TString) vWeight[i], "bWeight_misTagSFUp")) eventWeight *= bWeight_misTagSFUp;
      //else if (!strcmp((TString) vWeight[i], "bWeight_misTagSFDown")) eventWeight *= bWeight_misTagSFDown;
      //else if (!strncmp((TString) vWeight[i], "pdfWeights", 10)) {
      std::string sub = vWeight[i].substr(11);
      int pdfWeightN = atol(sub.c_str());
      //std::cout << "pdfWeightN: " << pdfWeightN << std::endl;
      eventWeight *= weightEvent_->pdfWeight[pdfWeightN];
    }

    if (eventWeight > random_->Uniform(0, maxMCWeight)) {
      fTree_->Fill();

      if (weightEvent_->mcWeight < 0) {
        eventsDrawn += -1;
        //progress    += -1;
      }
      else {
        ++eventsDrawn;
        ++progress;
      }
    }
  }
  std::cout << eventsDrawn << " events drawn in " << nAttempts << " attempts." << std::endl;
}

TTree* RandomSubsetCreatorNewInterface::PrepareTree(TString file) {
  TChain* chain = new TChain("analyzeKinFit/eventTree");
  chain->Add(file);

  //chain->SetBranchStatus("*",0);
  //chain->SetBranchStatus("top", 1);
  //chain->SetBranchStatus("weight", 1);

  TString selection(selection_);
  std::cout << file << ": " << chain->GetEntries(selection) << " events" << std::endl;

  return chain->CopyTree(selection);
}
