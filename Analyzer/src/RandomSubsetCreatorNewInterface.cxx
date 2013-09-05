#include "RandomSubsetCreatorNewInterface.h"

#include "Analysis.h"
#include "Helper.h"
#include "ProgramOptionsReader.h"

#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/progress.hpp>

//#include "TCanvas.h"
#include "TChain.h"
//#include "TEntryList.h"
//#include "TLorentzVector.h"
#include "TROOT.h"
//#include "TSystem.h"
#include "TTreeFormula.h"

//#include "LHAPDF/LHAPDF.h"

typedef ProgramOptionsReader po;

RandomSubsetCreatorNewInterface::RandomSubsetCreatorNewInterface() :
    // filled from program options
    selection_  (po::GetOption<std::string>("analysisConfig.selection")),
    samplePath_ (po::GetOption<std::string>("analysisConfig.samplePath")),
    fIdentifier_(po::GetOption<std::string>("input")),
    //fChannel_   (po::GetOption<std::string>("channel")),
    fVar1_      (po::GetOption<std::string>("analysisConfig.var1")),
    fVar2_      (po::GetOption<std::string>("analysisConfig.var2")),
    fVar3_      (po::GetOption<std::string>("analysisConfig.var3")),
    fWeight_    (po::GetOption<std::string>("weight")),
    fLumi_  (po::GetOption<double>("lumi")),
    fSig_   (po::GetOption<double>("fsig")),
    fBDisc_ (po::GetOption<double>("bdisc")),
    //topEvent_(new TopEvent()),
    //weightEvent_(new WeightEvent()),
    //fChain_(0),
    random_(0)
{
  channelID_ = Helper::channelID();

  //if (channelID_ == Helper::kAllJets) {
  //  fChainsSig_.push_back(PrepareEvents(samplePath_+fIdentifier_+TString(".root")));
  //  //if(fLumi_>0) fChainsBkg_.push_back(PrepareChain(samplePath_+"QCDEstimationMix.root"));
  //  if(fLumi_>0) fChainsBkg_.push_back(PrepareEvents(samplePath_+"QCDMixing_MJPS12*_data.root"));
  //}
  //if (channelID_ == Helper::kMuonJets || channelID_ == Helper::kLeptonJets) {
  //  fChainsSig_.push_back(PrepareEvents(samplePath_+fIdentifier_+TString("_muon/analyzeTop.root")));
  //  if(fLumi_>0) {
  //    fChainsBkg_.push_back(PrepareEvents("/scratch/hh/lustre/cms/user/mseidel/Fall11_Wbb_muon/analyzeTop.root"));
  //    fChainsBkg_.push_back(PrepareEvents("/scratch/hh/lustre/cms/user/mseidel/Fall11_T_muon/analyzeTop.root"));
  //  }
  //  //fTreeTTmu = PrepareTree(samplePath_+fIdentifier_+TString("_muon/analyzeTop.root"));
  //  //fTreeWmu  = PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_Wbb_muon/analyzeTop.root");
  //  //fTreeSTmu = PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_T_muon/analyzeTop.root");
  //}
  //if (channelID_ == Helper::kElectronJets || channelID_ == Helper::kLeptonJets) {
  //  fChainsSig_.push_back(PrepareEvents(samplePath_+fIdentifier_+TString("_electron/analyzeTop.root")));
  //  if(fLumi_>0) {
  //    fChainsBkg_.push_back(PrepareEvents("/scratch/hh/lustre/cms/user/mseidel/Fall11_Wbb_electron/analyzeTop.root"));
  //    fChainsBkg_.push_back(PrepareEvents("/scratch/hh/lustre/cms/user/mseidel/Fall11_T_electron/analyzeTop.root"));
  //  }
  //  //fTreeTTe  = PrepareTree(samplePath_+fIdentifier_+TString("_electron/analyzeTop.root"));
  //  //fTreeWe   = PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_Wbb_electron/analyzeTop.root");
  //  //fTreeSTe  = PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_T_electron/analyzeTop.root");
  //}
  if (channelID_ == Helper::kAllJets) {
    PrepareEvents(samplePath_+fIdentifier_+TString(".root"));
    if(fLumi_>0) PrepareEvents(samplePath_+"QCDMixing_MJPS12*_data.root");
  }
  if (channelID_ == Helper::kMuonJets || channelID_ == Helper::kLeptonJets) {
    PrepareEvents(samplePath_+fIdentifier_+TString("_muon/analyzeTop.root"));
    if(fLumi_>0) {
      PrepareEvents("/scratch/hh/lustre/cms/user/mseidel/Fall11_Wbb_muon/analyzeTop.root");
      PrepareEvents("/scratch/hh/lustre/cms/user/mseidel/Fall11_T_muon/analyzeTop.root");
    }
  }
  if (channelID_ == Helper::kElectronJets || channelID_ == Helper::kLeptonJets) {
    PrepareEvents(samplePath_+fIdentifier_+TString("_electron/analyzeTop.root"));
    if(fLumi_>0) {
      PrepareEvents("/scratch/hh/lustre/cms/user/mseidel/Fall11_Wbb_electron/analyzeTop.root");
      PrepareEvents("/scratch/hh/lustre/cms/user/mseidel/Fall11_T_electron/analyzeTop.root");
    }
  }
  random_ = new TRandom3(0);
  std::cout << "Random seed: " << random_->GetSeed() << std::endl;
}

RandomSubsetCreatorNewInterface::~RandomSubsetCreatorNewInterface()
{
  delete random_;
}

TTree* RandomSubsetCreatorNewInterface::CreateRandomSubset() {
  subset_.Clear();
  if (fLumi_>0) {
    std::cout << "Create random subset..." << std::endl;

    time_t start, end;
    time(&start);
    time(&end);

    //fChain_ = new TChain("analyzeKinFit/eventTree");

    // DATA
    double nEventsDataAllJets  = 11428.;
    double nEventsDataMuon     = 2906.;
    double nEventsDataElectron = 2268.;

    int eventsPEAllJets  = random_->Poisson(nEventsDataAllJets /18352.0 *fLumi_);
    int eventsPEMuon     = random_->Poisson(nEventsDataMuon    /5000.000*fLumi_);
    int eventsPEElectron = random_->Poisson(nEventsDataElectron/5000.000*fLumi_);

    if (channelID_ == Helper::kAllJets) {
      DrawEvents(events_.at(0), eventsPEAllJets*    fSig_ );
      DrawEvents(events_.at(1), eventsPEAllJets*(1.-fSig_));
    }
    if (channelID_ == Helper::kMuonJets || channelID_ == Helper::kLeptonJets) {
      DrawEvents(events_.at(0), eventsPEMuon*    fSig_       );
      DrawEvents(events_.at(1), eventsPEMuon*(1.-fSig_)*1./4.);
      DrawEvents(events_.at(2), eventsPEMuon*(1.-fSig_)*3./4.);
    }
    if (channelID_ == Helper::kElectronJets || channelID_ == Helper::kLeptonJets) {
      short offset = 0;
      if(channelID_ == Helper::kLeptonJets) offset = 1;
      DrawEvents(events_.at(3*offset+0), eventsPEElectron*    fSig_       );
      DrawEvents(events_.at(3*offset+1), eventsPEElectron*(1.-fSig_)*1./4.);
      DrawEvents(events_.at(3*offset+2), eventsPEElectron*(1.-fSig_)*3./4.);
    }

    time(&end);
    std::cout << "Created random subset in " << difftime(end, start) << " seconds." << std::endl;
  }
  else {
    subset_ = events_.at(0);
    if (channelID_ == Helper::kLeptonJets) {
      subset_ += events_.at(1);
    }
  }
  return 0;
}

void RandomSubsetCreatorNewInterface::DrawEvents(const DataSample& sample, double nEventsPE) {
  std::cout << "nEventsPE: " << nEventsPE << std::endl;

  //currentChain->SetBranchStatus("*",0);
  //currentChain->SetBranchStatus("top.*",1);
  //currentChain->SetBranchStatus("weight.*",1);
  //
  //currentChain->SetBranchAddress("top."   , &topEvent_);
  //currentChain->SetBranchAddress("weight.", &weightEvent_);
  //
  //TEntryList *selectedEvents   = currentChain->GetEntryList();
  //TEntryList *selectedEventsPE = new TEntryList(currentChain);

  int perms = sample.nEvents;

  //std::vector<std::string> vWeight;
  ////FIXME still needs to be implemented correctly, seems to work somehow
  //boost::split( vWeight, fWeight_, boost::is_any_of("-*"));
  //
  //double maxMCWeight = 1.; // calculate upper bound for combined MCWeight
  //
  //for (unsigned int i = 0; i < vWeight.size(); ++i) {
  //  std::cout << vWeight[i] << std::endl;
  //  if (strncmp((TString) vWeight[i], "weight.pdfWeights", 17)) maxMCWeight *= currentChain->GetMaximum((TString) vWeight[i]);
  //  else maxMCWeight *= currentChain->GetMaximum("weight.pdfWeights");
  //}

  double maxMCWeight = sample.maxWeight;

  std::cout << "maxMCWeight(" << fWeight_ << "):" << maxMCWeight  << std::endl;

  if (maxMCWeight ==  0) { std::cout << "Weight not active?" << std::endl; }
  if (maxMCWeight == -1) { std::cout << "Running over data?" << std::endl; }

  std::cout << "while (eventsDrawn < nEventsPE)..." << std::endl;

  int eventsDrawn = 0;
  int nAttempts = 0;

  boost::progress_display progress((int)nEventsPE, std::cout);

  while (eventsDrawn < (int)nEventsPE) {
    int drawn = random_->Integer(perms);

    //topEvent_->init();
    //weightEvent_->init();
    //
    //int treeNumber = 0;
    //int chainEntry = selectedEvents->GetEntryAndTree(drawn,treeNumber);
    //chainEntry += currentChain->GetTreeOffset()[treeNumber];
    //currentChain->GetEntry(chainEntry);
    //
    //++nAttempts;
    //double eventWeight = 1.;
    //for (unsigned int i = 0; i < vWeight.size(); ++i) {
    //  if      (!strcmp ((TString) vWeight[i], "weight.combinedWeight"         )) eventWeight *= weightEvent_->combinedWeight;
    //  else if (!strcmp ((TString) vWeight[i], "weight.mcWeight"               )) eventWeight *= weightEvent_->mcWeight;
    //  else if (!strcmp ((TString) vWeight[i], "weight.muWeight"               )) eventWeight *= weightEvent_->muWeight;
    //  else if (!strcmp ((TString) vWeight[i], "weight.puWeight"               )) eventWeight *= weightEvent_->puWeight;
    //  else if (!strcmp ((TString) vWeight[i], "weight.puWeightUp"             )) eventWeight *= weightEvent_->puWeightUp;
    //  else if (!strcmp ((TString) vWeight[i], "weight.puWeightDown"           )) eventWeight *= weightEvent_->puWeightDown;
    //  else if (!strcmp ((TString) vWeight[i], "weight.bTagWeight"             )) eventWeight *= weightEvent_->bTagWeight;
    //  else if (!strcmp ((TString) vWeight[i], "weight.bTagWeight_bTagSFUp"    )) eventWeight *= weightEvent_->bTagWeight_bTagSFUp;
    //  else if (!strcmp ((TString) vWeight[i], "weight.bTagWeight_bTagSFDown"  )) eventWeight *= weightEvent_->bTagWeight_bTagSFDown;
    //  else if (!strcmp ((TString) vWeight[i], "weight.bTagWeight_misTagSFUp"  )) eventWeight *= weightEvent_->bTagWeight_misTagSFUp;
    //  else if (!strcmp ((TString) vWeight[i], "weight.bTagWeight_misTagSFDown")) eventWeight *= weightEvent_->bTagWeight_misTagSFDown;
    //  else if (!strncmp((TString) vWeight[i], "weight.pdfWeights", 17)) {
    //    std::string sub = vWeight[i].substr(11);
    //    int pdfWeightN = atol(sub.c_str());
    //    //std::cout << "pdfWeightN: " << pdfWeightN << std::endl;
    //    eventWeight *= weightEvent_->pdfWeight[pdfWeightN];
    //  }
    //}
    //
    //if (eventWeight > random_->Uniform(0, maxMCWeight)) {
    //  if(selectedEventsPE->Enter(chainEntry, currentChain)){
    //    if (weightEvent_->mcWeight < 0) {
    //      eventsDrawn += -1;
    //      //progress    += -1;
    //    }
    //    else {
    //      ++eventsDrawn;
    //      ++progress;
    //    }
    //  }
    //}

    if (std::abs(sample.weights.at(drawn)) > random_->Uniform(0, maxMCWeight)) {
        if (sample.weights.at(drawn) < 0) {
          eventsDrawn += -1;
          //progress    += -1;
        }
        else {
          ++eventsDrawn;
          ++progress;
      }
    }
  }

  //fChain_->Add(currentChain);
  //TEntryList* entryList = fChain_->GetEntryList();
  //if(entryList){
  //  entryList->Add(selectedEventsPE);
  //  fChain_->SetEntryList(entryList);
  //}
  //else{
  //  fChain_->SetEntryList(selectedEventsPE);
  //}

  std::cout << eventsDrawn << " events drawn in " << nAttempts << " attempts." << std::endl;
}

void RandomSubsetCreatorNewInterface::PrepareEvents(TString file) {
  TChain* chain = new TChain("analyzeKinFit/eventTree");
  int nFiles = chain->Add(file);
  std::cout << "Adding " << nFiles << " files for " << file << std::flush;

  TTreeFormula *f1     = new TTreeFormula("f1"    , fVar1_  .c_str(), chain);
  TTreeFormula *f2     = new TTreeFormula("f2"    , fVar2_  .c_str(), chain);
  TTreeFormula *f3     = new TTreeFormula("f3"    , fVar3_  .c_str(), chain);
  TTreeFormula *weight = new TTreeFormula("weight", fWeight_.c_str(), chain);
  TTreeFormula *index  = new TTreeFormula("index" , "Iteration$"    , chain);
  TTreeFormula *sel    = new TTreeFormula("sel"   , selection_      , chain);

  DataSample sample;

  int selected = 0;
  for(int i = 0; ; ++i){
    if(chain->LoadTree(i) < 0) break;
    if(!sel->GetNdata()) continue;
    //std::cout << f1->GetNdata() << " " << f2->GetNdata() << " " << f3->GetNdata() << " " << weight->GetNdata() << " " << index->GetNdata() << " " << sel->GetNdata() << std::endl;
    //std::cout << f1->GetExpFormula() << " " << f2->GetExpFormula() << " " << f3->GetExpFormula() << " " << weight->GetExpFormula() << " " << index->GetExpFormula() << " " << sel->GetExpFormula() << std::endl;
    //std::cout << f1->EvalInstance(0) << " " << f2->EvalInstance(0) << " " << f3->EvalInstance(0) << " " << weight->EvalInstance(0) << " " << index->EvalInstance(0) << " " << sel->EvalInstance(0) << std::endl;
    for(unsigned char j = 0, l = index->GetNdata(); j < l; ++j)
      sample.Fill(f1->EvalInstance(j), f2->EvalInstance(j), f3->EvalInstance(j), weight->EvalInstance(j), (unsigned char)(index->EvalInstance(j)));
    ++selected;
  }

  sample.nEvents = selected;
  std::cout << ": " << selected << " events \n(" << selection_ << ")" << std::endl;

  //TString selectionName = file;
  //selectionName.ReplaceAll("*","");
  //selectionName.ReplaceAll("/","");
  //selectionName.ReplaceAll(".","");
  //
  //chain->Draw(">>selectedEvents_"+selectionName,selection_,"entrylist");
  //TEntryList *selectedEvents = (TEntryList*)gDirectory->Get("selectedEvents_"+selectionName);
  //std::cout << ": " << selectedEvents->GetN() << " events \n(" << selection_ << ")" << std::endl;
  //chain->SetEntryList(selectedEvents);

  events_.push_back(sample);

  //return chain;
}
