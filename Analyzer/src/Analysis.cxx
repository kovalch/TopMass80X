#include "Analysis.h"

#include <string>

#include "TCanvas.h"
#include "TChain.h"
#include "TEntryList.h"
#include "TROOT.h"

#include "GenMatchAnalyzer.h"
#include "Helper.h"
#include "IdeogramAnalyzer.h"
#include "IdeogramAnalyzerNewInterface.h"
#include "MVAAnalyzer.h"
#include "ProgramOptionsReader.h"
#include "RandomSubsetCreatorLeptonJets.h"
#include "RandomSubsetCreatorAllJets.h"
#include "RandomSubsetCreatorNewInterface.h"
#include "RooFitTemplateAnalyzer.h"

#include "LHAPDF/LHAPDF.h"

typedef ProgramOptionsReader po;

Analysis::Analysis(std::vector<float> v):
  fIdentifier_(po::GetOption<std::string>("input" )),
  fMethod_    (po::GetOption<std::string>("method")),
  fBinning_   (po::GetOption<std::string>("binning")),
  vBinning_   (v),
  fChannelID_ (Helper::channelID()),
  fMethodID_  (Helper::methodID()),
  fCreator_(0),
  fTree_(0)
{
}

Analysis::~Analysis()
{
  for(std::map<TString, TH1F*>::iterator hist = histograms_.begin(); hist != histograms_.end(); ++hist){
    delete hist->second;
  }
  histograms_.clear();
  //delete fTree_; // deletion is taken care of by MassAnalyzer
  delete fCreator_;
}

void Analysis::Analyze() {

  std::cout << "Analyze " << fIdentifier_ << " with method " << fMethod_ << std::endl;

  // random subset creation
  if(!fCreator_){
    if (fChannelID_ == Helper::kElectronJets || fChannelID_ == Helper::kMuonJets || fChannelID_ == Helper::kLeptonJets) {
      fCreator_ = new RandomSubsetCreatorLeptonJets();
    }
    else if (fChannelID_ == Helper::kAllJets) {
      if (fMethodID_ == Helper::kIdeogramNew) {
        fCreator_ = new RandomSubsetCreatorNewInterface();
      }
      else{
        fCreator_ = new RandomSubsetCreatorAllJets();
      }
    }
    else {
      std::cerr << "Stopping analysis! Specified decay channel *" << po::GetOption<std::string>("channel") << "* not known!" << std::endl;
      return;
    }
  }
  fTree_ = fCreator_->CreateRandomSubset();

  MassAnalyzer* fAnalyzer = 0;

  if      (fMethodID_ == Helper::kGenMatch   ) fAnalyzer = new GenMatchAnalyzer            (fIdentifier_, fTree_);
  else if (fMethodID_ == Helper::kMVA        ) fAnalyzer = new MVAAnalyzer                 (fIdentifier_, fTree_);
  else if (fMethodID_ == Helper::kIdeogram   ) fAnalyzer = new IdeogramAnalyzer            (fIdentifier_, fTree_);
  else if (fMethodID_ == Helper::kIdeogramNew) fAnalyzer = new IdeogramAnalyzerNewInterface(fIdentifier_, fTree_);
  else if (fMethodID_ == Helper::kRooFit     ) fAnalyzer = new RooFitTemplateAnalyzer      (fIdentifier_, fTree_);
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
    // calculate cuts
    std::stringstream stream;
    stream << vBinning_[i] << " < " << fBinning_ << " & "
           << fBinning_ << " < " << vBinning_[i+1];
    TString cuts = stream.str();
    
    if (fMethodID_ == Helper::kGenMatch) {
      cuts += " & target == 1";
    }
    else if (fMethodID_ == Helper::kMVA) {
      cuts += " & mvaDisc > 0";
    }
    else if (fMethodID_ == Helper::kIdeogram) {
    }
    else if (fMethodID_ == Helper::kIdeogramNew) {
    }
    else if (fMethodID_ == Helper::kRooFit) {
    }

    std::cout << fTree_->GetEntryList()->GetN() << std::endl;

    TCanvas* canvDummy = new TCanvas();
    fTree_->Draw("top.fitTop1[0].M()");
    canvDummy->Print("canvDummy.eps");

    fTree_->Draw(">>selectedEventsBin",cuts,"entrylist");
    TEntryList *selectedEventsBin = (TEntryList*)gDirectory->Get("selectedEventsBin");
    fTree_->SetEntryList(selectedEventsBin);

    int entries = selectedEventsBin->GetN();
    std::cout << cuts << std::endl;
    std::cout << entries << std::endl;

    CreateHisto("Entries");
    GetH1("Entries")->SetBinContent(i+1, entries);

    if (entries > 25) {
      fAnalyzer->Analyze(cuts, i, 0);
      const std::map<TString, std::pair<double, double>> values = fAnalyzer->GetValues();

      for(std::map<TString, std::pair<double, double>>::const_iterator value = values.begin(); value != values.end(); ++value){
        CreateHisto(value->first);
        CreateHisto(value->first+TString("_Error"));
        CreateHisto(value->first+TString("_Pull"));
      }

      double genMass = po::GetOption<double>("mass");
      double genJES  = po::GetOption<double>("jes" );
      double genfSig = po::GetOption<double>("fsig");
      for(std::map<TString, std::pair<double, double>>::const_iterator value = values.begin(); value != values.end(); ++value){
        double val      = value->second.first;
        double valError = value->second.second;
        double gen = 0;
        if     (value->first.BeginsWith("mass")) gen = genMass;
        else if(value->first.BeginsWith("JES" )) gen = genJES;
        else if(value->first.BeginsWith("fSig")) gen = genfSig;
        GetH1(value->first                 ) ->SetBinContent(i+1, val);
        GetH1(value->first+TString("_Error"))->SetBinContent(i+1, valError);
        GetH1(value->first+TString("_Pull" ))->SetBinContent(i+1, (val - gen)/valError);
        std::cout << "Measured " << value->first << ": " << val << " +/- " << valError << std::endl;
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
  GetH1("mass_mTop_JES")->SetAxisRange(GetH1("mass_mTop_JES")->GetMinimum(0.05), GetH1("mass_mTop_JES")->GetMaximum(), "Z");
  GetH1("mass_mTop_JES")->Fit("pol0");

  canvas->cd(3);
  GetH1("mass_mTop_JES_Error")->Draw("E1");
  GetH1("mass_mTop_JES_Error")->SetAxisRange(0.05, 5, "Z");
  GetH1("mass_mTop_JES_Error")->Fit("pol0");
  
  /*
  canvas->cd(4);
  hJES->Draw("E1");
  hJES->Fit("pol0");
  */
  // clean binning to give a proper path
  TString binningForPath = fBinning_;
  binningForPath.ReplaceAll("[","_"); binningForPath.ReplaceAll("]","_"); binningForPath.ReplaceAll("(","_"); binningForPath.ReplaceAll(")","_");
  binningForPath.ReplaceAll(".","_"); binningForPath.ReplaceAll("/","_");
  TString localIdentifier = fIdentifier_; localIdentifier.ReplaceAll("*","_"); localIdentifier.ReplaceAll("/","_");
  TString path("plot/"); path += fMethod_; path += "_"; path += localIdentifier; path += "_"; path += binningForPath; path += ".eps";
  canvas->Print(path);
  
  TString pathr("plot/"); pathr += fMethod_; pathr += "_"; pathr += localIdentifier; pathr += "_"; pathr += binningForPath; pathr += ".root";
  canvas->Print(pathr);
  
  delete canvas;
  delete fAnalyzer;
}

void Analysis::CreateHisto(TString name) {
  std::map<TString, TH1F*>::iterator hist = histograms_.find(name);
  if(hist != histograms_.end()){
    hist->second->Reset();
  }
  else{
    gROOT->cd();
    Helper* helper = new Helper(fBinning_, vBinning_);
    SetH1(name, helper->GetH1(name));
    delete helper;
  }
}

TH1F*
Analysis::GetH1(TString histName){
  std::map<TString, TH1F*>::const_iterator hist_iterator = histograms_.find(histName);
  if(hist_iterator != histograms_.end()){
    return hist_iterator->second;
  }
  else{
    std::cerr << "The searched histogram *" << histName << "* does not exist!" << std::endl;
    raise(SIGINT);
  }
  return 0;
}

const std::map<TString, TH1F*>
Analysis::GetH1s() const{
  return histograms_;
}

void Analysis::SetH1(TString histName, TH1F* hist){
  histograms_[histName] = hist;
}

TString Analysis::GetIdentifier() {
  return fIdentifier_;
}
