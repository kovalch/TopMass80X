#include "Analysis.h"

#include <string>

#include "TCanvas.h"
#include "TROOT.h"

#include "GenMatchAnalyzer.h"
#include "Helper.h"
#include "IdeogramAnalyzer.h"
#include "MVAAnalyzer.h"
#include "RandomSubsetCreatorAllJets.h"
#include "RooFitTemplateAnalyzer.h"
#include "XMLConfigReader.h"

#include "LHAPDF/LHAPDF.h"

typedef XMLConfigReader xml;

Analysis::Analysis(po::variables_map vm) :
  fIdentifier_(vm["input" ].as<std::string>()),
  fMethod_    (vm["method"].as<std::string>()),
  fBins_      (vm["bins"  ].as<int>()),
  fChannel_(xml::GetParameter("decayChannel")),
  fTree_(0)
{
}

Analysis::Analysis(TString identifier, TString file, TString method, int bins, double lumi) :
  fIdentifier_(identifier),
  fMethod_(method), fBins_(bins),
  fChannel_(xml::GetParameter("decayChannel")),
  fTree_(0)
{
}

Analysis::~Analysis()
{
  for(std::map<TString, TH2F*>::iterator hist = histograms_.begin(); hist != histograms_.end(); ++hist){
    delete hist->second;
  }
  histograms_.clear();
  //delete fTree_; // deletion is taken care of by MassAnalyzer
}

void Analysis::Analyze(po::variables_map vm) {

  std::cout << "Analyze " << fIdentifier_ << " with method " << fMethod_ << std::endl;

  // random subset creation
  RandomSubsetCreator* fCreator = 0;
  if (!strcmp(fChannel_, "AllJets")) {
    fCreator = new RandomSubsetCreatorAllJets(vm);
  }
  else {
    std::cerr << "Stopping analysis! Specified decay channel *" << fChannel_ << "* not known!" << std::endl;
    return;
  }
  fTree_ = fCreator->CreateRandomSubset();

  MassAnalyzer* fAnalyzer = 0;

  if (!strcmp(fMethod_, "GenMatch")) {
    fAnalyzer = new GenMatchAnalyzer(fIdentifier_, fTree_);
  }
  else if (!strcmp(fMethod_, "MVA")) {
    fAnalyzer = new MVAAnalyzer(fIdentifier_, fTree_);
  }
  else if (!strcmp(fMethod_, "Ideogram")) {
    fAnalyzer = new IdeogramAnalyzer(fIdentifier_, fTree_);
  }
  else if (!strcmp(fMethod_, "RooFit")) {
    fAnalyzer = new RooFitTemplateAnalyzer(fIdentifier_, fTree_);
  }
  else {
    std::cerr << "Stopping analysis! Specified analysis method *" << fMethod_ << "* not known!" << std::endl;
    return;
  }

  Helper* helper = new Helper(fBins_);
  helper->SetTDRStyle();
  delete helper;

  TCanvas* canvas = new TCanvas("canvas", "Top mass", 900, 600);
  canvas->cd();

  double smearBins = 1.;
  double rangeX = 10.;
  double rangeY = 10.;

  double smearX = smearBins/fBins_*rangeX;
  double smearY = smearBins/fBins_*rangeY;

  double minEntries = 25;

  TString observableX = "dRbb";
  TString observableY = "dRbb";

  for(int i = 0; i < fBins_; i++) {
    for(int j = 0; j < fBins_; j++) {
      // calculate cuts
      std::stringstream stream;
      stream << (rangeX/fBins_)*(i)-smearX << "<" << observableX << "&"
             << observableX << "<" << (rangeX/fBins_)*(i+1)+smearX << " & "
             << rangeY/fBins_*(j)-smearY << "<" << observableY << "&"
             << observableY << "<" << rangeY/fBins_*(j+1)+smearY;
      TString cuts(stream.str());

      if (!strcmp(fMethod_, "GenMatch")) {
        cuts += " & comboType == 1";
      }
      else if (!strcmp(fMethod_, "MVA")) {
        cuts += " & mvaDisc > 0";
      }
      else if (!strcmp(fMethod_, "Ideogram")) {
//        //cuts += " & target == 1";
//        //cuts += " & (target == 0 | target == 1)";
//        //cuts += " & Njet == 6";
//        cuts += " & probs[0] > 0.09";
//        cuts += " & dRbb > 1.5";
//        //cuts += " & nVertex >= 5";
//        //cuts += " & nVertex <= 5";
//        cuts += " & topMasses[0] > 100. & topMasses[0] < 550.";
      }
      else if (!strcmp(fMethod_, "RooFit")) {
      }

      int entries = fTree_->GetEntries(cuts);
      std::cout << cuts << std::endl;
      std::cout << entries << std::endl;

      CreateHisto("Entries");
      GetH2("Entries")->SetCellContent(i+1, j+1, entries);

      if (entries > minEntries) {
        fAnalyzer->Analyze(cuts, i, j, vm);
        const std::map<TString, std::pair<double, double> > values = fAnalyzer->GetValues();

        for(std::map<TString, std::pair<double, double> >::const_iterator value = values.begin(); value != values.end(); ++value){
          CreateHisto(value->first);
          CreateHisto(value->first+TString("_Error"));
          CreateHisto(value->first+TString("_Pull")); //FIXME this dummy histogram has to be ALPHABETICALLY after the two others (like it is now)
        }
        std::cout << std::endl;
        for(std::map<TString, std::pair<double, double> >::const_iterator value = values.begin(); value != values.end(); ++value){
          double val      = value->second.first;
          double valError = value->second.second;
          GetH2(value->first                 ) ->SetCellContent(i+1, j+1, val);
          GetH2(value->first+TString("_Error"))->SetCellContent(i+1, j+1, valError);
          GetH2(value->first+TString("_Pull" ))->SetCellContent(i+1, j+1, -1.); //FIXME dummy histogram needed in top mass
          std::cout << "Measured " << value->first << ": " << val << " +/- " << valError << std::endl;
        }
        std::cout << std::endl;
      }
    }
  }

  canvas->Clear();
  canvas->Divide(2,2);

  canvas->cd(1);
  GetH2("Entries")->Draw("COLZ");

  canvas->cd(2);
  GetH2("mass_mTop_JES")->Draw("COLZ,TEXT");
  GetH2("mass_mTop_JES")->SetAxisRange(GetH2("mass_mTop_JES")->GetMinimum(0.05), GetH2("mass_mTop_JES")->GetMaximum(), "Z");

  canvas->cd(3);
  GetH2("mass_mTop_JES_Error")->Draw("COLZ,TEXT");
  GetH2("mass_mTop_JES_Error")->SetAxisRange(0.05, 5, "Z");

  TString path("plot/"); path += fMethod_; path += "_"; path += fIdentifier_; path += ".eps";
  canvas->Print(path);

  delete canvas;
  delete fAnalyzer;
  delete fCreator;
}

void Analysis::CreateHisto(TString name) {
  std::map<TString, TH2F*>::iterator hist = histograms_.find(name);
  if(hist != histograms_.end()){
    hist->second->Reset();
  }
  else{
    gROOT->cd();
    Helper* helper = new Helper(fBins_);
    SetH2(name, helper->GetH2(name));
    delete helper;
  }
}

TH2F*
Analysis::GetH2(TString histName){
  std::map<TString, TH2F*>::const_iterator hist_iterator = histograms_.find(histName);
  if(hist_iterator != histograms_.end()){
    return hist_iterator->second;
  }
  else{
    std::cerr << "The searched histogram *" << histName << "* does not exist!" << std::endl;
    assert(0);
  }
  return 0;
}

const std::map<TString, TH2F*>
Analysis::GetH2s() const{
  return histograms_;
}

void Analysis::SetH2(TString histName, TH2F* hist){
  histograms_[histName] = hist;
}

TString Analysis::GetIdentifier() {
  return fIdentifier_;
}
