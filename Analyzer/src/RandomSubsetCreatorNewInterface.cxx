#include "RandomSubsetCreatorNewInterface.h"

//#include "Analysis.h"
#include "Helper.h"
#include "ProgramOptionsReader.h"

#include <iostream>
//#include <boost/progress.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "TChain.h"
#include "TRandom3.h"
#include "TTreeFormula.h"

typedef ProgramOptionsReader po;

RandomSubsetCreatorNewInterface::RandomSubsetCreatorNewInterface(const std::vector<float>& v) :
    // filled from program options
    selection_  (po::GetOption<std::string>("analysisConfig.selection")),
    samplePath_ (po::GetOption<std::string>("analysisConfig.samplePath")),
    fIdentifier_(po::GetOption<std::string>("input")),
    fVar1_      (po::GetOption<std::string>("analysisConfig.var1")),
    fVar2_      (po::GetOption<std::string>("analysisConfig.var2")),
    fVar3_      (po::GetOption<std::string>("analysisConfig.var3")),
    fVar4_      (po::GetOption<std::string>("analysisConfig.var4")),
    fWeight_    (po::GetOption<std::string>("weight")),
    activeBranches_(po::GetOption<std::string>("analysisConfig.activeBranches")),
    fBinning_   (po::GetOption<std::string>("binning")),
    vBinning_   (v),
    fLumi_  (po::GetOption<double>("lumi")),
    fSig_   (po::GetOption<double>("fsig")),
    //fBDisc_ (po::GetOption<double>("bdisc")),
    maxPermutations_(po::GetOption<int>("analysisConfig.maxPermutations")),
    random_(0)
{
  channelID_ = Helper::channelID();

  std::cout << "Reading data from disk ..." << std::endl;
  std::cout << "Event selection: " << selection_ << std::endl;
  std::cout << "Variable 1: " << fVar1_ << std::endl;
  std::cout << "Variable 2: " << fVar2_ << std::endl;
  std::cout << "Variable 3: " << fVar3_ << std::endl;
  std::cout << "Variable 4: " << fVar4_ << std::endl;
  std::cout << "Binning: " << fBinning_ << std::endl;
  std::cout << "Weight: " << fWeight_ << std::endl;
  std::cout << "Lumi: " << fLumi_ << std::endl;
  std::cout << "Signal Fraction: " << fSig_ << std::endl;

  time_t start, end;
  time(&start);
  time(&end);

  if (channelID_ == Helper::kAllJets) {
    PrepareEvents(samplePath_+fIdentifier_+std::string(".root"));
    if(fLumi_>0) PrepareEvents(samplePath_+"QCDMixing_MJPS12_v1_data.root");
  }
  if (channelID_ == Helper::kMuonJets || channelID_ == Helper::kLeptonJets) {
    PrepareEvents(samplePath_+fIdentifier_+std::string("_muon/job_*.root"));
    if(fLumi_>0 && fSig_<1.) {
      PrepareEvents(""+samplePath_+"Summer12_WJets_muon/job_*.root");
      PrepareEvents(""+samplePath_+"Summer12_singleTop_muon/job_*.root");
    }
  }
  if (channelID_ == Helper::kElectronJets || channelID_ == Helper::kLeptonJets) {
    PrepareEvents(samplePath_+fIdentifier_+std::string("_electron/job_*.root"));
    if(fLumi_>0 && fSig_<1.) {
      PrepareEvents(""+samplePath_+"Summer12_WJets_electron/job_*.root");
      PrepareEvents(""+samplePath_+"Summer12_singleTop_electron/job_*.root");
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
    double nEventsDataAllJets  =  6792.;
    double nEventsDataMuon     = 15172.;
    double nEventsDataElectron = 13937.;

    int eventsPEAllJets  = random_->Poisson(nEventsDataAllJets /18291.000*fLumi_);
    int eventsPEMuon     = random_->Poisson(nEventsDataMuon    /19712.000*fLumi_);
    int eventsPEElectron = random_->Poisson(nEventsDataElectron/19712.000*fLumi_);

    if (channelID_ == Helper::kAllJets) {
      DrawEvents(events_.at(0), eventsPEAllJets*    fSig_ );
      DrawEvents(events_.at(1), eventsPEAllJets*(1.-fSig_));
    }
    if (channelID_ == Helper::kMuonJets || channelID_ == Helper::kLeptonJets) {
      DrawEvents(events_.at(0), eventsPEMuon*    fSig_       );
      if (fSig_<1.) {
        DrawEvents(events_.at(1), eventsPEMuon*(1.-fSig_)*2./5.);
        DrawEvents(events_.at(2), eventsPEMuon*(1.-fSig_)*3./5.);
      }
    }
    if (channelID_ == Helper::kElectronJets || channelID_ == Helper::kLeptonJets) {
      short offset = 0;
      if(channelID_ == Helper::kLeptonJets) offset = 1;
      if (fSig_<1.) {
        DrawEvents(events_.at(3*offset+0), eventsPEElectron*    fSig_       );
        DrawEvents(events_.at(3*offset+1), eventsPEElectron*(1.-fSig_)*1./3.);
        DrawEvents(events_.at(3*offset+2), eventsPEElectron*(1.-fSig_)*2./3.);
      }
      else {
        DrawEvents(events_.at(1), eventsPEElectron);
      }
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
  //std::cout << "nEventsPE: " << nEventsPE << std::endl;

  int perms = sample.nEvents;

  double maxMCWeight = sample.maxWeight;

  std::cout << "maxMCWeight(" << fWeight_ << "): " << maxMCWeight  << std::endl;

  if (maxMCWeight ==  0) { std::cout << "Weight not active?" << std::endl; }
  if (maxMCWeight == -1) { std::cout << "Running over data?" << std::endl; }

  //std::cout << "while (eventsDrawn < nEventsPE)..." << std::endl;

  int eventsDrawn = 0;
  int nAttempts = 0;

  //boost::progress_display progress((int)nEventsPE, std::cout);

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
        //++progress;
      }
    }
  }

  std::cout << eventsDrawn << " events drawn in " << nAttempts << " attempts." << std::endl;
}

void RandomSubsetCreatorNewInterface::PrepareEvents(const std::string& file) {

  TChain* chain;
  if (channelID_ == Helper::kAllJets) {
    chain = new TChain("analyzeKinFit/eventTree");
  }
  else {
    chain = new TChain("analyzeHitFit/eventTree");
  }
  int nFiles = chain->Add(file.c_str());
  std::cout << "Adding " << nFiles << " files for " << file << std::flush;

  chain->SetBranchStatus("*", 0);
  std::vector<std::string> vActiveBanches;
  boost::split(vActiveBanches, activeBranches_, boost::is_any_of("|"));
  for(const auto& branch : vActiveBanches){
    chain->SetBranchStatus(branch.c_str(), 1);
  }

  TTreeFormula *f1     = new TTreeFormula("f1"    , fVar1_    .c_str(), chain);
  TTreeFormula *f2     = new TTreeFormula("f2"    , fVar2_    .c_str(), chain);
  TTreeFormula *f3     = new TTreeFormula("f3"    , fVar3_    .c_str(), chain);
  TTreeFormula *f4     = new TTreeFormula("f4"    , fVar4_    .c_str(), chain);
  TTreeFormula *binning= new TTreeFormula("binning", fBinning_.c_str(), chain);
  TTreeFormula *weight = new TTreeFormula("weight", fWeight_  .c_str(), chain);
  TTreeFormula *sel    = new TTreeFormula("sel"   , selection_.c_str(), chain);

  DataSample sample;

  int selected = 0;
  for(int i = 0; ; ++i){
    long entry = chain->LoadTree(i);
    if(entry  < 0) break;
    if(entry == 0){
      f1    ->UpdateFormulaLeaves();
      f2    ->UpdateFormulaLeaves();
      f3    ->UpdateFormulaLeaves();
      f4    ->UpdateFormulaLeaves();
      binning->UpdateFormulaLeaves();
      weight->UpdateFormulaLeaves();
      sel   ->UpdateFormulaLeaves();
    }
    if(!f1    ->GetNdata()) continue;
    if(!f2    ->GetNdata()) continue;
    if(!f3    ->GetNdata()) continue;
    if(!f4    ->GetNdata()) continue;
    if(!binning->GetNdata()) continue;
    if(!weight->GetNdata()) continue;
    if(!sel   ->GetNdata()) continue;
    int filledPermutations = 0;
    for(int j = 0, l = std::min(maxPermutations_, sel->GetNdata()); j < l; ++j){
      if(!sel->EvalInstance(j)) continue;
      
      int bin = 0;
      for (const auto& boundary : vBinning_) {
        if (binning->EvalInstance(j) > boundary) ++bin;
      }

      sample.Fill(f1->EvalInstance(j), f2->EvalInstance(j), f3->EvalInstance(j), f4->EvalInstance(j), weight->EvalInstance(j), filledPermutations++, bin);
    }
    if(filledPermutations) ++selected;
  }

  sample.nEvents = selected;
  std::cout << ": " << selected << " events" << std::endl;

  events_.push_back(sample);
  delete chain;
  delete f1;
  delete f2;
  delete f3;
  delete f4;
  delete binning;
  delete weight;
  delete sel;
}
