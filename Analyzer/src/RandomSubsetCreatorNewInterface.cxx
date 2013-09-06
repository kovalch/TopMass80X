#include "RandomSubsetCreatorNewInterface.h"

#include "Analysis.h"
#include "Helper.h"
#include "ProgramOptionsReader.h"

#include <iostream>
#include <boost/progress.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "TChain.h"
#include "TTreeFormula.h"

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
    activeBranches_(po::GetOption<std::string>("analysisConfig.activeBranches")),
    fLumi_  (po::GetOption<double>("lumi")),
    fSig_   (po::GetOption<double>("fsig")),
    fBDisc_ (po::GetOption<double>("bdisc")),
    random_(0)
{
  channelID_ = Helper::channelID();

  std::cout << "Reading data from disk ..." << std::endl;
  std::cout << "Event selection: " << selection_ << std::endl;

  time_t start, end;
  time(&start);
  time(&end);

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
  time(&end);
  std::cout << "Read data from disk in " << difftime(end, start) << " seconds." << std::endl;

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

  int perms = sample.nEvents;

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
    ++nAttempts;

    if (std::abs(sample.events.at(drawn).weight) > random_->Uniform(0, maxMCWeight)) {
      subset_.AddEvent(sample.events.at(drawn));
      if (sample.events.at(drawn).weight < 0) {
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

void RandomSubsetCreatorNewInterface::PrepareEvents(TString file) {

  TChain* chain = new TChain("analyzeKinFit/eventTree");
  int nFiles = chain->Add(file);
  std::cout << "Adding " << nFiles << " files for " << file << std::flush;

  chain->SetBranchStatus("*", 0);
  std::vector<std::string> vActiveBanches;
  boost::split( vActiveBanches, activeBranches_, boost::is_any_of("|"));
  for(auto branch : vActiveBanches){
    chain->SetBranchStatus(branch.c_str(), 1);
  }

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
    if(chain->LoadTree(i) == 0){
      f1    ->UpdateFormulaLeaves();
      f2    ->UpdateFormulaLeaves();
      f3    ->UpdateFormulaLeaves();
      weight->UpdateFormulaLeaves();
      index ->UpdateFormulaLeaves();
      sel   ->UpdateFormulaLeaves();
    }
    if(!sel->GetNdata()) continue;
    if(!sel->EvalInstance(0)) continue;
    for(int j = 0, l = index->GetNdata(); j < l; ++j){
      if(!sel->EvalInstance(j)) continue;
      sample.Fill(f1->EvalInstance(j), f2->EvalInstance(j), f3->EvalInstance(j), weight->EvalInstance(j), (int)(index->EvalInstance(j)));
      //if(i < 10) std::cout << f1->EvalInstance(j) << " " << f2->EvalInstance(j) << " " << f3->EvalInstance(j) << std::endl;
    }
    ++selected;
  }

  sample.nEvents = selected;
  std::cout << ": " << selected << " events" << std::endl;

  events_.push_back(sample);
  delete chain;
}
