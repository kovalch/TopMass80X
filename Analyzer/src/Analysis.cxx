#include "Analysis.h"

#include <iostream>
#include <sstream>
#include <csignal>

#include "TCanvas.h"
#include "TH1F.h"
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

//debug
//#include "TH1D.h"
//#include "TFile.h"

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
//debug
/*TFile* fDebug=new TFile("fDebug.root","UPDATE");
TH1D* hVal = (TH1D*)fDebug->Get("hVal");
TH1D* hGen =(TH1D*)fDebug->Get("hGen");
TH1D* hValErr = (TH1D*)fDebug->Get("hValErr");*/

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
  fTree_ = fCreator_->CreateRandomSubset(); //fuer RandomSubsetCreatorNewInterface == 0 !!!



  if      (fMethodID_ == Helper::kGenMatch   ) fAnalyzer_ = new GenMatchAnalyzer            (fIdentifier_, fTree_);
  else if (fMethodID_ == Helper::kMVA        ) fAnalyzer_ = new MVAAnalyzer                 (fIdentifier_, fTree_);
  else if (fMethodID_ == Helper::kIdeogram   ) fAnalyzer_ = new IdeogramAnalyzer            (fIdentifier_, fTree_);
  else if (fMethodID_ == Helper::kIdeogramNew) { if(!fAnalyzer_) { fAnalyzer_ = new IdeogramAnalyzerNewInterface(fIdentifier_, fTree_); } }
  else if (fMethodID_ == Helper::kIdeogramMin) { if(!fAnalyzer_) { fAnalyzer_ = new IdeogramAnalyzerMinimizer   (fIdentifier_,  0 /*fTree_*/); } } //fTree_ wird nicht verwendet (anders wenn nicht RandomSubsetCreatorNewInterface verwendet wird)
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
    // FIXME binning not yet implemented, still needs some implementation ...
    
    // FIXME Entries per bin     C: wie ist das Binning formatiert?? for-Schleife hat fuer fitTop1[0].M grad ja eh nur einen Durchlauf
    int entries = ((RandomSubsetCreatorNewInterface*)fCreator_)->GetDataSample().nEvents;
    CreateHisto("Entries");
    GetH1("Entries")->SetBinContent(i+1, 10+i);

//debugging
 //std::cout << "DEBUG Analysis DEBUG: her i am with ((RandomSubsetCreatorNewInterface*)fCreator_)->GetDataSample().nEvents = "<<entries << std::endl;

    if (entries > 25) {
      if      (fMethodID_ == Helper::kIdeogramMin) ((IdeogramAnalyzerMinimizer   *)fAnalyzer_)->SetDataSample(((RandomSubsetCreatorNewInterface*)fCreator_)->GetDataSample()); //hier bekommt er die Pseudodaten
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
//debug
	/*if(value.first=="mass_mTop_JES"){
	std::cout<<"DEBUG check the debug Fill "<<hVal->Fill(val)<<std::endl;
		hGen->Fill(gen);
		hValErr->Fill(valError);
	}*/
      }
      std::cout << std::endl;
    }

//debug
		/*std::cout<<"DEBUG is dDebug Open? "<<fDebug->IsOpen()<<fDebug->cd()<<std::endl;
		std::cout<<"DEBUG check the size of hval "<<hVal->Write()<<std::endl;
		hGen->Write();
		hValErr->Write();
fDebug->Close("R");*/
  }

  canvas->Clear();
  canvas->Divide(2,3);

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
  
  if (po::GetOption<bool>("constrainJSF")) {
    canvas->cd(5);
    GetH1("JES_mTop_JES_jsfc")->Draw("E1");
    GetH1("JES_mTop_JES_jsfc")->SetAxisRange(GetH1("JES_mTop_JES_jsfc")->GetMinimum(0.05)-0.01, GetH1("JES_mTop_JES_jsfc")->GetMaximum()+0.01, "Z");
    GetH1("JES_mTop_JES_jsfc")->Fit("pol0");
    
    canvas->cd(6);
    GetH1("mass_mTop_JES_jsfc")->Draw("E1");
    GetH1("mass_mTop_JES_jsfc")->SetAxisRange(GetH1("mass_mTop_JES_jsfc")->GetMinimum(0.05)-1., GetH1("mass_mTop_JES_jsfc")->GetMaximum()+1, "Z");
    GetH1("mass_mTop_JES_jsfc")->Fit("pol0");
  }

  std::string path("plot/"); path += fMethod_; path += "_"; path += HelperFunctions::cleanedName(fIdentifier_); path += "_"; path += HelperFunctions::cleanedName(fBinning_);
  path += ".eps";
  canvas->Print(path.c_str());
  boost::replace_all(path,".eps",".root");
  canvas->Print(path.c_str());
  
  delete canvas;
  if (fMethodID_ != Helper::kIdeogramNew && fMethodID_ != Helper::kIdeogramMin) delete fAnalyzer_;  //wird also atm nicht aufgerufen?? (wuerde eh nur "fTree_=0" deleten...)
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
