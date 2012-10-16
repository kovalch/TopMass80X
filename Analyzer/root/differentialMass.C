#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TChain.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"

#include "tdrstyle.C"

TH1F* hData;
TH1F* hTTJets096;
TH1F* hTTJets100;
TH1F* hTTJets104;
TH1F* hTTJetsFlavorUp;
TH1F* hTTJetsFlavorDown;
TH1F* hTTJetsP11;
TH1F* hTTJetsP11noCR;

enum lepton           { kElectron, kMuon, kAll};
TString lepton_ [3] = { "electron", "muon", "all"};

TString sBinning[8] = {"hadTopPt", "hadTopEta", "hadBPt", "hadBEta", "TTBarMass", "TTBarPt", "deltaRHadQHadQBar", "deltaRHadBLepB"};
TString sBinNice[8] = {"p_{T,t,had} [GeV]", "#eta_{t,had}", "p_{T,b,had} [GeV]", "#eta_{b,had}", "m_{t#bar{t}} [GeV]", "p_{T,t#bar{t}} [GeV]", "#DeltaR_{q#bar{q}}", "#DeltaR_{b#bar{b}}"};
int jBinning = 3;

TString sObservable[4] = {"Entries", "Mass", "MassAlt", "JES"};
TString sObsCanvas[4] = {"canvas_1", "canvas_2", "canvas_3", "canvas_4"};
TString sObsNice[4] = {"Fraction of entries / bin", "m_{t}^{2D} - <m_{t}^{2D}> [GeV]", "m_{t}^{1D} - <m_{t}^{1D}> [GeV]", "JES - <JES>"};
int jObs = 0;

int channel = 2;

int nbins = 10;
double xfirst = 0;
double xlast  = 1;

int color_ [] =  {kRed+1, kRed-7, kRed-10, kGray , kGreen+1, kRed+1 , kAzure-2, kGreen-3, 
                 kYellow , kMagenta, 10      , kBlack  , 
                 kYellow , kYellow , kYellow , kYellow , kYellow , kYellow ,
                 10      , 10      , 10      , 
                 kMagenta, kMagenta, kMagenta, kMagenta, kMagenta, kMagenta };

int marker_[] = {20, 22, 22, 20, 22, 23, 29, 23, 
                 21, 27, 28, 20, 
                 21, 21, 21, 21, 21, 21,
                 28, 28, 28,
                 27, 27, 27, 27, 27, 27};

struct sample {
  TChain* chain;
  TFile* file;
  TCanvas* canvas;
  TPad* pad;
  TH1F* profile;
  TH1F* error;
  TF1* fit;
  TGraphErrors* graph;
  const char* id;
  const char* label;
  int color;
  int line;
  sample(TH1F* p, const char* i, const char* l, int c = kBlack, int li = 1)
  : profile(p), id(i), label(l), color(c), line(li) {}
};

void differentialMass(int iBinning = jBinning, int iObs = jObs, bool batch = false)
{
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);
  
  std::vector<sample> samples;
  std::vector<sample>::iterator it;
  
  // Open input files
  TString sfBase("/afs/naf.desy.de/user/m/mseidel/scratch/CMSSW_4_2_4/src/TopMass/Analyzer/plot/Ideogram_");
  TString sfData; sfData += sfBase; sfData += "Run2011_"; sfData += sBinning[iBinning]; sfData += ".root";
  
  samples.push_back(sample(hData, "Run2011", "Data", kBlack));
  samples.push_back(sample(hTTJets100, "Fall11_TTJets1725_1.00", "Tune Z2", kRed+1));
  samples.push_back(sample(hTTJetsP11, "Fall11_TTJets1725_P11", "Tune P11", kMagenta+1, 9));
  samples.push_back(sample(hTTJetsP11noCR, "Fall11_TTJets1725_P11noCR", "Tune P11noCR", kCyan+1, 2));
    
  //"[0] + exp([1]+[2]*x)"
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    
    TString sf; sf += sfBase; sf += it->id; sf += "_"; sf += sBinning[iBinning]; sf += ".root";
    it->file    = new TFile(sf);
    it->canvas  = (TCanvas*) it->file  ->Get("canvas");
    it->pad     = (TPad*)    it->canvas->GetPrimitive(sObsCanvas[iObs]);
    it->profile = (TH1F*)    it->pad   ->GetPrimitive(sObservable[iObs]);
    
    if (iObs != 0) {
      it->fit     = (TF1*)     it->pad   ->GetPrimitive("pol0");
      it->fit->SetLineColor(it->color);
      it->fit->SetLineStyle(1);
      it->fit->SetLineWidth(0);
    }
    
    it->profile->SetLineColor(it->color);
    it->profile->SetLineWidth(2);
    it->profile->SetLineStyle(it->line);
    if (it->color == kBlack) it->profile->SetMarkerStyle(20);
    
    double integral = it->profile->Integral();
    
    for (int i = 0; i < it->profile->GetNbinsX()+2; i++) {
      if (iObs != 0) it->profile->SetBinContent(i, it->profile->GetBinContent(i) - it->fit->GetParameter(0));
      else {
        it->profile->SetBinContent(i, it->profile->GetBinContent(i) / integral);
        it->profile->SetBinError(i, it->profile->GetBinError(i) / integral);
      }
    }
    
    it->error = (TH1F*)it->profile->Clone();
    it->error->SetMarkerStyle(0);
    it->error->SetFillColor(it->color-11);
    it->error->SetFillStyle(3254);
  }
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
  
  TH1F* hNull = new TH1F("null", "", samples[0].profile->GetNbinsX(), samples[0].profile->GetBinLowEdge(1), samples[0].profile->GetBinLowEdge(samples[0].profile->GetNbinsX()+1));
  hNull->GetYaxis()->SetRangeUser(samples[0].profile->GetMinimum() - 2*samples[0].profile->GetBinError(samples[0].profile->GetMinimumBin()),
                                  samples[0].profile->GetMaximum() + 3*samples[0].profile->GetBinError(samples[0].profile->GetMaximumBin()));
  hNull->GetXaxis()->SetTitle(sBinNice[iBinning]);
  hNull->GetYaxis()->SetTitle(sObsNice[iObs]);
  hNull->Draw();
  
  TLegend *leg1 = new TLegend(0.6, 0.75, 0.9, 0.925);
  leg1->SetTextSize(0.03);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    if (it->color != kBlack) {
      it->error->Draw("SAME,E2");
    }
  }
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    if (it->color == kBlack) {
      it->profile->Draw("SAME,E");
      leg1->AddEntry(it->profile, it->label, "P");
    }
    else {
      it->profile->Draw("SAME,HIST");
      leg1->AddEntry(it->profile, it->label, "L");
    }
  }
  
  leg1->Draw();
  
  TString path("differential/diff"); path += sObservable[iObs]; path += "_"; path += sBinning[iBinning]; path += ".eps";
  canvas->Print(path);
  
  if (batch == true) {
    // cleanup
    canvas->Clear();
    leg1->Clear();
    delete canvas;
    delete leg1;
    
    for (it = samples.begin(); it != samples.end(); ++it) {
      it->file->Close();
    }
  }
}

void differentialBatch() {
  for (int b = 0; b < 8; ++b) {
    for (int o = 0; o < 4; ++o) {
      differentialMass(b, o, true);
    }
  }
}
