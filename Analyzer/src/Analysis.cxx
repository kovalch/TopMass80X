#include "Analysis.h"

#include <iostream>
#include <sstream>
#include <csignal>

#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"

#include "GenMatchAnalyzer.h"
#include "Helper.h"
#include "IdeogramAnalyzer.h"
#include "IdeogramAnalyzerNewInterface.h"
#include "IdeogramAnalyzerMinimizer.h"
#include "MVAAnalyzer.h"
#include "ProgramOptionsReader.h"
#include "RandomSubsetCreatorLeptonJets.h"
#include "RandomSubsetCreatorAllJets.h"
#include "RandomSubsetCreatorNewInterface.h"
#include "RooFitTemplateAnalyzer.h"

typedef ProgramOptionsReader po;

Analysis::Analysis(const std::vector<float>& v):
  fIdentifier_(po::GetOption<std::string>("input" )),
  fMethod_    (po::GetOption<std::string>("method")),
  fBinning_   (po::GetOption<std::string>("binning")),
  vBinning_   (v),
  fChannelID_ (Helper::channelID()),
  fMethodID_  (Helper::methodID()),
  fAnalyzer_(0),
  fCreator_(0),
  fTree_(0)
{
}

Analysis::~Analysis()
{
  for(std::map<std::string, TH1F*>::iterator hist = histograms_.begin(); hist != histograms_.end(); ++hist){
    delete hist->second;
  }
  histograms_.clear();
  //delete fTree_; // deletion is taken care of by MassAnalyzer
  if (fMethodID_ == Helper::kIdeogramNew || fMethodID_ == Helper::kIdeogramMin) delete fAnalyzer_;
  delete fCreator_;
}

void Analysis::Analyze() {

  std::cout << "Analyze " << fIdentifier_ << " with method " << fMethod_ << std::endl;

  // random subset creation
  if(!fCreator_){
    if (fMethodID_ == Helper::kIdeogramNew || fMethodID_ == Helper::kIdeogramMin) {
      fCreator_ = new RandomSubsetCreatorNewInterface(vBinning_);
    }
    else{
      if (fChannelID_ == Helper::kElectronJets || fChannelID_ == Helper::kMuonJets || fChannelID_ == Helper::kLeptonJets) {
        fCreator_ = new RandomSubsetCreatorLeptonJets();
      }
      else if (fChannelID_ == Helper::kAllJets) {
        fCreator_ = new RandomSubsetCreatorAllJets();
      }
      else {
        std::cerr << "Stopping analysis! Specified decay channel *" << po::GetOption<std::string>("channel") << "* not known!" << std::endl;
        return;
      }
    }
  }
  fTree_ = fCreator_->CreateRandomSubset();

  if      (fMethodID_ == Helper::kGenMatch   ) fAnalyzer_ = new GenMatchAnalyzer            (fIdentifier_, fTree_);
  else if (fMethodID_ == Helper::kMVA        ) fAnalyzer_ = new MVAAnalyzer                 (fIdentifier_, fTree_);
  else if (fMethodID_ == Helper::kIdeogram   ) fAnalyzer_ = new IdeogramAnalyzer            (fIdentifier_, fTree_);
  else if (fMethodID_ == Helper::kIdeogramNew) { if(!fAnalyzer_) { fAnalyzer_ = new IdeogramAnalyzerNewInterface(fIdentifier_, fTree_); } }
  else if (fMethodID_ == Helper::kIdeogramMin) { if(!fAnalyzer_) { fAnalyzer_ = new IdeogramAnalyzerMinimizer   (fIdentifier_, fTree_); } }
  else if (fMethodID_ == Helper::kRooFit     ) fAnalyzer_ = new RooFitTemplateAnalyzer      (fIdentifier_, fTree_);
  else {
    std::cerr << "Stopping analysis! Specified analysis method *" << fMethod_ << "* not known!" << std::endl;
    return;
  }
   
  Helper* helper = new Helper(fBinning_, vBinning_);
  helper->SetTDRStyle();
  delete helper;

  TCanvas* canvas = new TCanvas("canvas", "Top mass", 900, 600);
  canvas->cd();
  
  for(unsigned int i = 0; i < vBinning_.size()-1; ++i) {
    int entries = ((RandomSubsetCreatorNewInterface*)fCreator_)->GetDataSample().nEvents;

    if (entries > 25) {
      if      (fMethodID_ == Helper::kIdeogramMin) ((IdeogramAnalyzerMinimizer   *)fAnalyzer_)->SetDataSample(((RandomSubsetCreatorNewInterface*)fCreator_)->GetDataSample());
      else if (fMethodID_ == Helper::kIdeogramNew) ((IdeogramAnalyzerNewInterface*)fAnalyzer_)->SetDataSample(((RandomSubsetCreatorNewInterface*)fCreator_)->GetDataSample());
      fAnalyzer_->Analyze("", i+1, 0);
      const std::map<std::string, std::pair<double, double>> values = fAnalyzer_->GetValues();

      for(const auto& value: values){
        CreateHisto(value.first);
        CreateHisto(value.first+std::string("_Error"));
        CreateHisto(value.first+std::string("_Pull"));
      }

      double genMass = po::GetOption<double>("mass");
      double genJES  = po::GetOption<double>("jes" );
      double genfSig = po::GetOption<double>("fsig");
      for(const auto& value: values){
        double val      = value.second.first;
        double valError = value.second.second;
        double gen = 0;
        if     (!strncmp(value.first.c_str(),"mass",4)) gen = genMass;
        else if(!strncmp(value.first.c_str(),"JES" ,3)) gen = genJES;
        else if(!strncmp(value.first.c_str(),"fSig",4)) gen = genfSig;
        GetH1(value.first                 ) ->SetBinContent(i+1, val);
        GetH1(value.first                 ) ->SetBinError  (i+1, valError);
        GetH1(value.first+std::string("_Error"))->SetBinContent(i+1, valError);
        GetH1(value.first+std::string("_Pull" ))->SetBinContent(i+1, (val - gen)/valError);
        std::cout << "Measured " << value.first << ": " << val << " +/- " << valError << std::endl;
      }
      std::cout << std::endl;
    }
  }

  canvas->Clear();
  canvas->Divide(2,2);

  canvas->cd(1);
  GetH1("Entries")->Draw("E1");

  canvas->cd(2);
  GetH1("mass_mTop_JES")->Draw("E1");
  GetH1("mass_mTop_JES")->SetAxisRange(GetH1("mass_mTop_JES")->GetMinimum(0.05)-1., GetH1("mass_mTop_JES")->GetMaximum()+1, "Z");
  GetH1("mass_mTop_JES")->Fit("pol0");

  canvas->cd(3);
  GetH1("JES_mTop_JES")->Draw("E1");
  GetH1("JES_mTop_JES")->SetAxisRange(GetH1("JES_mTop_JES")->GetMinimum(0.05)-0.01, GetH1("JES_mTop_JES")->GetMaximum()+0.01, "Z");
  GetH1("JES_mTop_JES")->Fit("pol0");
  
  canvas->cd(4);
  GetH1("mass_mTop")->Draw("E1");
  GetH1("mass_mTop")->SetAxisRange(GetH1("mass_mTop")->GetMinimum(0.05)-1., GetH1("mass_mTop")->GetMaximum()+1, "Z");
  GetH1("mass_mTop")->Fit("pol0");
  
  // clean binning to give a proper path
  TString binningForPath = fBinning_;
  binningForPath.ReplaceAll("[","_"); binningForPath.ReplaceAll("]","_"); binningForPath.ReplaceAll("(","_"); binningForPath.ReplaceAll(")","_");
  binningForPath.ReplaceAll(".","_"); binningForPath.ReplaceAll("/","_"); binningForPath.ReplaceAll("@","_");
  binningForPath.ReplaceAll(",","_"); binningForPath.ReplaceAll(" ","_"); binningForPath.ReplaceAll(":","_");
  TString localIdentifier = fIdentifier_; localIdentifier.ReplaceAll("*","_"); localIdentifier.ReplaceAll("/","_");
  std::string path("plot/"); path += fMethod_; path += "_"; path += localIdentifier; path += "_"; path += binningForPath; path += ".eps";
  canvas->Print(path.c_str());
  
  std::string pathr("plot/"); pathr += fMethod_; pathr += "_"; pathr += localIdentifier; pathr += "_"; pathr += binningForPath; pathr += ".root";
  canvas->Print(pathr.c_str());
  
  delete canvas;
  if (fMethodID_ != Helper::kIdeogramNew && fMethodID_ != Helper::kIdeogramMin) delete fAnalyzer_;
}

void Analysis::CreateHisto(std::string name) {
  std::map<std::string, TH1F*>::iterator hist = histograms_.find(name);
  if(hist == histograms_.end()){
    gROOT->cd();
    Helper* helper = new Helper(fBinning_, vBinning_);
    SetH1(name, helper->GetH1(name));
    delete helper;
  }
}

TH1F*
Analysis::GetH1(std::string histName){
  std::map<std::string, TH1F*>::const_iterator hist_iterator = histograms_.find(histName);
  if(hist_iterator != histograms_.end()){
    return hist_iterator->second;
  }
  else{
    std::cerr << "The searched histogram *" << histName << "* does not exist!" << std::endl;
    raise(SIGINT);
  }
  return 0;
}

const std::map<std::string, TH1F*>
Analysis::GetH1s() const{
  return histograms_;
}

void Analysis::SetH1(std::string histName, TH1F* hist){
  histograms_[histName] = hist;
}

std::string Analysis::GetIdentifier() {
  return fIdentifier_;
}
