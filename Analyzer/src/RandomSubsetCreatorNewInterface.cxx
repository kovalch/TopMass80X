#include "RandomSubsetCreatorNewInterface.h"

#include "ProgramOptionsReader.h"

#include <iostream>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "TChain.h"
#include "TRandom3.h"
#include "TTreeFormula.h"

typedef ProgramOptionsReader po;

RandomSubsetCreatorNewInterface::RandomSubsetCreatorNewInterface(const std::vector<float>& v) :
    // filled from program options
    selection_     (po::GetOption<std::string>("analysisConfig.selection")),
    selectionLept_ (po::GetOption<std::string>("analysisConfig.selectionLept")),
    selectionJets_ (po::GetOption<std::string>("analysisConfig.selectionJets")),
    samplePath_    (po::GetOption<std::string>("analysisConfig.samplePath")),
    samplePathLept_(po::GetOption<std::string>("analysisConfig.samplePathLept")),
    samplePathJets_(po::GetOption<std::string>("analysisConfig.samplePathJets")),
    fIdentifier_   (po::GetOption<std::string>("input")),
    fVar1_         (po::GetOption<std::string>("analysisConfig.var1")),
    fVar2_         (po::GetOption<std::string>("analysisConfig.var2")),
    fVar3_         (po::GetOption<std::string>("analysisConfig.var3")),
    fVar4_         (po::GetOption<std::string>("analysisConfig.var4")),
    fWeight_       (po::GetOption<std::string>("weight")),
    activeBranches_(po::GetOption<std::string>("analysisConfig.activeBranches")),
    fBinning_      (po::GetOption<std::string>("binning")),
    vBinning_   (v),
    fLumi_    (po::GetOption<double>("lumi")),
    fSig_     (po::GetOption<double>("fsig")),
    fLumiLept_(po::GetOption<double>("lumiLept")),
    fLumiJets_(po::GetOption<double>("lumiJets")),
    fSigLept_ (po::GetOption<double>("fsigLept")),
    fSigJets_ (po::GetOption<double>("fsigJets")),
    //fBDisc_ (po::GetOption<double>("bdisc")),
    maxPermutations_(po::GetOption<int>("analysisConfig.maxPermutations")),
    random_(0)
{
  channelID_ = Helper::channelID();
  mergedsample_.nEvents = 0;
  mergedsample_.maxWeight = 0.;

  std::cout << "Reading data from disk ..." << std::endl;
  if(channelID_ == Helper::kHamburg) {
    std::cout << "Event selection LepJets: " << selectionLept_ << std::endl;
    std::cout << "Event selection AllJets: " << selectionJets_ << std::endl;
  }
  else{
    std::cout << "Event selection: " << selection_ << std::endl;
  }
  std::cout << "Variable 1: " << fVar1_ << std::endl;
  std::cout << "Variable 2: " << fVar2_ << std::endl;
  std::cout << "Variable 3: " << fVar3_ << std::endl;
  std::cout << "Variable 4: " << fVar4_ << std::endl;
  std::cout << "Binning: " << fBinning_ << std::endl;
  std::cout << "Weight: " << fWeight_ << std::endl;
  if(channelID_ == Helper::kHamburg) {
    std::cout << "Lumi LepJets: " << fLumiLept_ << std::endl;
    std::cout << "Lumi AllJets: " << fLumiJets_ << std::endl;
  }
  else{
    std::cout << "Lumi: " << fLumi_ << std::endl;
    if(channelID_ == Helper::kMuonJets || channelID_ == Helper::kElectronJets || channelID_ == Helper::kLeptonJets) fLumiLept_ = fLumi_;
    if(channelID_ == Helper::kAllJets) fLumiJets_ = fLumi_;
  }
  if(channelID_ == Helper::kHamburg) {
    std::cout << "Signal Fraction LepJets: " << fSigLept_ << std::endl;
    std::cout << "Signal Fraction AllJets: " << fSigJets_ << std::endl;
  }
  else{
    std::cout << "Signal Fraction: " << fSig_ << std::endl;
    if(channelID_ == Helper::kMuonJets || channelID_ == Helper::kElectronJets || channelID_ == Helper::kLeptonJets) fSigLept_ = fSig_;
    if(channelID_ == Helper::kAllJets) fSigJets_ = fSig_;
  }

  time_t start, end;
  time(&start);
  time(&end);
  
  Helper* helper = new Helper();
  std::vector<std::string> backgroundSamplesLept = helper->readParametersString("analysisConfig.backgroundSamplesLept");
  std::vector<double>      backgroundFactorsLept = helper->readParameters("analysisConfig.backgroundFactorsLept");
  if(backgroundSamplesLept.size()!=backgroundFactorsLept.size()) std::cerr << "Error: Background samples and factors do not match. You have to check this!" << std::endl;

  // TODO: Background samples and normalization from config in kHamburg and kAllJets
  if (channelID_ == Helper::kHamburg) {
    PrepareEvents(samplePathJets_+fIdentifier_+std::string("_alljets.root"), Helper::kAllJets);
    PrepareEvents(samplePathLept_+fIdentifier_+std::string("_muon/job_*.root"), Helper::kLeptonJets);
    PrepareEvents(samplePathLept_+fIdentifier_+std::string("_electron/job_*.root"), Helper::kLeptonJets);
    if(fLumiLept_>0 || fLumiJets_>0){
      //PrepareEvents(samplePathJets_+"QCDMixing_MJPS12_v1_data.root", Helper::kAllJets);
      //PrepareEvents(samplePathJets_+"Run2012_Background_alljets.root", Helper::kAllJets);
      PrepareEvents(samplePathJets_+"Run2012_Mixing8_fix_alljets.root", Helper::kAllJets);
      if(fSigLept_<1.) {
        PrepareEvents(""+samplePathLept_+"Summer12_WJets_muon/job_*.root", Helper::kLeptonJets);
        PrepareEvents(""+samplePathLept_+"Summer12_singleTop_muon/job_*.root", Helper::kLeptonJets);
        PrepareEvents(""+samplePathLept_+"Summer12_WJets_electron/job_*.root", Helper::kLeptonJets);
        PrepareEvents(""+samplePathLept_+"Summer12_singleTop_electron/job_*.root", Helper::kLeptonJets);
      }
    }
  }
  else{
    if (channelID_ == Helper::kAllJets) {
      PrepareEvents(samplePath_+fIdentifier_+std::string(".root"));
      if(fLumiJets_>0){
        if(fIdentifier_.find("BackgroundSystematic")!=std::string::npos)
          PrepareEvents(samplePath_+"QCDMixing_Z2_S12_Madspin_sig.root");
        else
          //PrepareEvents(samplePath_+"QCDMixing_MJPS12_v1_data.root");
          //PrepareEvents(samplePath_+"Run2012_Background_alljets.root");
          //PrepareEvents(samplePath_+"Run2012_Mixing8_alljets.root");
          PrepareEvents(samplePath_+"Run2012_Mixing8_fix_alljets.root");
      }
    }
    if (channelID_ == Helper::kMuonJets || channelID_ == Helper::kLeptonJets) {
      PrepareEvents(samplePath_+fIdentifier_+std::string("_muon/job_*.root"), Helper::kLeptonJets, po::GetOption<double>("analysisConfig.signalFactor"));
      // Get background samples and normalization from config
      if(fLumiLept_>0 && fSigLept_<1.) {
        int iBkg = 0;
        for (std::string backgroundSampleLept : backgroundSamplesLept) {
          PrepareEvents(""+samplePath_+backgroundSampleLept+"_muon/job_*.root", Helper::kLeptonJets, backgroundFactorsLept[iBkg]);
          ++iBkg;
        }
      }
    }
    if (channelID_ == Helper::kElectronJets || channelID_ == Helper::kLeptonJets) {
      PrepareEvents(samplePath_+fIdentifier_+std::string("_electron/job_*.root"), Helper::kLeptonJets, po::GetOption<double>("analysisConfig.signalFactor"));
       if(fLumiLept_>0 && fSigLept_<1.) {
        int iBkg = 0;
        for (std::string backgroundSampleLept : backgroundSamplesLept) {
          PrepareEvents(""+samplePath_+backgroundSampleLept+"_electron/job_*.root", Helper::kLeptonJets, backgroundFactorsLept[iBkg]);
          ++iBkg;
        }
      }
    }
  }
  time(&end);
  std::cout << "Read data from disk in " << difftime(end, start) << " seconds." << std::endl;

  random_ = new TRandom3(po::GetOption<int>("seed")+1);
  std::cout << "Random seed: " << random_->GetSeed() << std::endl;
}

RandomSubsetCreatorNewInterface::~RandomSubsetCreatorNewInterface()
{
  delete random_;
}

TTree* RandomSubsetCreatorNewInterface::CreateRandomSubset() {
  subset_.Clear();
  if (fLumiLept_>0 || fLumiJets_>0) {
    std::cout << "Create random subset..." << std::endl;

    time_t start, end;
    time(&start);
    time(&end);

    // DATA
    // TODO: update AllJets event yields with latest JEC
    double nEventsDataAllJets  =  7049.;
    double nEventsDataMuon     = 14685.;
    double nEventsDataElectron = 13514.;

    int eventsPEAllJets  = random_->Poisson(nEventsDataAllJets /18192.000*fLumiJets_);
    int eventsPEMuon     = random_->Poisson(nEventsDataMuon    /19712.000*fLumiLept_);
    int eventsPEElectron = random_->Poisson(nEventsDataElectron/19712.000*fLumiLept_);

    // TODO: Use mergedsample_ in kHamburg and kAllJets
    if (channelID_ == Helper::kHamburg) {
      DrawEvents(events_.at(0), eventsPEAllJets *fSigJets_);
      DrawEvents(events_.at(1), eventsPEMuon    *fSigLept_);
      DrawEvents(events_.at(2), eventsPEElectron*fSigLept_);

      DrawEvents(events_.at(3), eventsPEAllJets*(1.-fSigJets_));
      if (fSigLept_<1.) {
        DrawEvents(events_.at(4), eventsPEMuon*(1.-fSigLept_)*2./5.);
        DrawEvents(events_.at(5), eventsPEMuon*(1.-fSigLept_)*3./5.);
        DrawEvents(events_.at(6), eventsPEElectron*(1.-fSigLept_)*1./3.);
        DrawEvents(events_.at(7), eventsPEElectron*(1.-fSigLept_)*2./3.);
      }
    }
    if (channelID_ == Helper::kAllJets) {
      DrawEvents(events_.at(0), eventsPEAllJets*    fSigJets_ );
      DrawEvents(events_.at(1), eventsPEAllJets*(1.-fSigJets_));
    }
    if (channelID_ == Helper::kMuonJets) {
      DrawEvents(mergedsample_, eventsPEMuon);
    }
    if (channelID_ == Helper::kElectronJets) {
      DrawEvents(mergedsample_, eventsPEElectron);
    }
    if (channelID_ == Helper::kLeptonJets) {
      DrawEvents(mergedsample_, eventsPEMuon+eventsPEElectron);
    }

    time(&end);
    std::cout << "Created random subset in " << difftime(end, start) << " seconds." << std::endl;
  }
  else {
    subset_ = events_.at(0);
    if (channelID_ == Helper::kLeptonJets) {
      subset_ += events_.at(1);
    }
    if (channelID_ == Helper::kHamburg) {
      subset_ += events_.at(1);
      subset_ += events_.at(2);
    }
  }
  return 0;
}

void RandomSubsetCreatorNewInterface::DrawEvents(const DataSample& sample, double nEventsPE) {
  int perms = sample.nEvents;
  std::cout << "Total number of permuations: " << perms << std::endl;
  
  if (perms == 0) return;

  double maxMCWeight = sample.maxWeight;

  std::cout << "maxMCWeight(" << fWeight_ << "): " << maxMCWeight  << std::endl;

  if (maxMCWeight ==  0) { std::cout << "Weight not active?" << std::endl; }
  if (maxMCWeight == -1) { std::cout << "Running over data?" << std::endl; }

  int eventsDrawn = 0;
  int nAttempts = 0;

  while (eventsDrawn < (int)nEventsPE) {
    int drawn = random_->Integer(perms);
    ++nAttempts;

    if (std::abs(sample.events.at(drawn).weight) > random_->Uniform(0, maxMCWeight)) {
      subset_.AddEvent(sample.events.at(drawn));
      if (sample.events.at(drawn).weight < 0) {
        eventsDrawn += -1;
      }
      else {
        ++eventsDrawn;
      }
    }
  }

  std::cout << eventsDrawn << " events drawn in " << nAttempts << " attempts." << std::endl;
}

void RandomSubsetCreatorNewInterface::PrepareEvents(const std::string& file, const Helper::ChannelID currentID, double sampleFactor) {

  TChain* chain;
  if (channelID_ == Helper::kAllJets || currentID == Helper::kAllJets) {
    chain = new TChain("analyzeKinFit/eventTree");
    if (Helper::getCMSEnergy() == 7) {
      chain = new TChain("FullHadTreeWriter/tree");
    }
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
  TTreeFormula *sel    = 0;
  if(channelID_ == Helper::kHamburg) {
    if(currentID == Helper::kAllJets) {
      sel = new TTreeFormula("sel", selectionJets_.c_str(), chain);
    }
    else if(currentID == Helper::kLeptonJets) {
      sel = new TTreeFormula("sel", selectionLept_.c_str(), chain);
    }
  }
  else{
    sel = new TTreeFormula("sel"   , selection_.c_str(), chain);
  }

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
    for(int j = 0, l = std::min(((currentID == Helper::kAllJets) ? maxPermutations_ : 4), sel->GetNdata()); j < l; ++j){
      if(!sel->EvalInstance(j)) continue;
      if(!weight->EvalInstance(j)) continue;
      
      int bin = 0;
      for (const auto& boundary : vBinning_) {
        if (binning->EvalInstance(j) > boundary) ++bin;
      }
      sample.Fill(f1->EvalInstance(j), f2->EvalInstance(j), f3->EvalInstance(j), f4->EvalInstance(j), weight->EvalInstance(j)*sampleFactor, filledPermutations++, bin);
    }
    if(filledPermutations) ++selected;
  }

  sample.nEvents = selected;
  std::cout << ": " << selected << " events" << std::endl;

  events_.push_back(sample);
  mergedsample_ += sample;
  delete chain;
  delete f1;
  delete f2;
  delete f3;
  delete f4;
  delete binning;
  delete weight;
  delete sel;
}
