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

TChain* tTTJets096;
TChain* tTTJets100;
TChain* tTTJets104;
TChain* tTTJetsFlavorUp;
TChain* tTTJetsFlavorDown;
TChain* tTTJetsP11;
TChain* tTTJetsP11noCR;
TChain* tData;

enum lepton           { kElectron, kMuon, kAll};
TString lepton_ [3] = { "electron", "muon", "all"};

int channel = 0;

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
  TH1F* profile;
  const char* label;
  int color;
  sample(TChain* ch, const char* l, int c = kBlack)
  : chain(ch), label(l), color(c) {}
};


std::vector<sample> samples;
std::vector<sample>::iterator it;

TH1F* responseProfile(TString sDraw, TChain* chain, TString label, int color) {
  sDraw += " >> h2"; sDraw += label; sDraw += sDraw(0, 7); sDraw += "(15, 0, 150, 500, 0, 5)";
  chain->Draw(sDraw, "(0.5*MCWeight)*(hitFitProb > 0.2 & hadQBSSV<1.74 & hadQBarBSSV<1.74 & hadBBSSV>1.74 & lepBBSSV>1.74 & leptonPt > 30)");
  TString sH2("h2"); sH2 += label; sH2 += sDraw(0, 7);
  TH2F* h2 = (TH2F*) gDirectory->Get(sH2);
  
  h2->FitSlicesY(0, 0, -1, 30);
  /*
  hTTJets_1->GetYaxis()->SetRangeUser(ylow, yup);
  hTTJets_1->SetLineColor(color_[kSigCP]);
  hTTJets_1->SetMarkerStyle(marker_[kSigCP]);
  hTTJets_1->SetMarkerColor(color_[kSigCP]);
  hTTJets_1->Fit("fitTTJets", "EM");
  */
  TString sH2_pfx("h2"); sH2_pfx += label; sH2_pfx += sDraw(0, 7); sH2_pfx += "_1";
  TH1F* h2_pfx = (TH1F*) gDirectory->Get(sH2_pfx);
  
  TF1* constFit = new TF1("constFit", "[0]+(x-7)*[1]");
  constFit->SetParNames("offset", "slope");
  constFit->SetLineColor(color);
  constFit->SetLineWidth(1);
  constFit->SetLineStyle(7);
  h2_pfx->Fit("constFit");
  
  return h2_pfx;
}

void massResponse2D()
{
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  
  // Define chains
  TString sTree("analyzeHitFit/eventTree");
  tTTJets096        = new TChain(sTree);
  tTTJets100        = new TChain(sTree);
  tTTJets104        = new TChain(sTree);
  tTTJetsFlavorUp   = new TChain(sTree);
  tTTJetsFlavorDown = new TChain(sTree);
  tTTJetsP11        = new TChain(sTree);
  tTTJetsP11noCR    = new TChain(sTree);
  tData             = new TChain(sTree);
  
  
  // ---
  //    open input files
  // ---
  if (channel == kMuon || channel == kAll) {
    tTTJets096        ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_0.96_muon/analyzeTop.root");
    tTTJets100        ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_1.00_muon/analyzeTop.root");
    tTTJets104        ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_1.04_muon/analyzeTop.root");
    tTTJetsFlavorUp   ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_flavor:up_muon/analyzeTop.root");
    tTTJetsFlavorDown ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_flavor:down_muon/analyzeTop.root");
    tTTJetsP11        ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_P11_muon/analyzeTop.root");
    tTTJetsP11noCR    ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_P11noCR_muon/analyzeTop.root");
    tData             ->Add("/scratch/hh/current/cms/user/mseidel/Run2011_muon/analyzeTop.root");
  }
  if (channel == kElectron || channel == kAll) {
    tTTJets096        ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_0.96_electron/analyzeTop.root");
    tTTJets100        ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_1.00_electron/analyzeTop.root");
    tTTJets104        ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_1.04_electron/analyzeTop.root");
    tTTJetsFlavorUp   ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_flavor:up_electron/analyzeTop.root");
    tTTJetsFlavorDown ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_flavor:down_electron/analyzeTop.root");
    tTTJetsP11        ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_P11_electron/analyzeTop.root");
    tTTJetsP11noCR    ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_P11noCR_electron/analyzeTop.root");
    tData             ->Add("/scratch/hh/current/cms/user/mseidel/Run2011_electron/analyzeTop.root");
  }
  
  //samples.push_back(sample(tData, "Data", kBlack));
  samples.push_back(sample(tTTJets100, "Z2", kBlue+1));
  //samples.push_back(sample(tTTJets096, "Z2, JES -4%", kRed+1));
  //samples.push_back(sample(tTTJets104, "Z2, JES +4%", kGreen+1));
  //samples.push_back(sample(tTTJetsFlavorDown, "Z2_bDown", kRed));
  //samples.push_back(sample(tTTJetsFlavorUp, "Z2_bUp", kGreen));
  samples.push_back(sample(tTTJetsP11, "P11", kMagenta+1));
  samples.push_back(sample(tTTJetsP11noCR, "P11noCR", kYellow+1));
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    it->profile  = responseProfile("deltaRHadQHadQBar:hadQBarGen.Pt()", it->chain, it->label, it->color);
  }
  
  TLegend *leg1 = new TLegend(0.6, 0.75, 0.9, 0.925);
  leg1->SetTextSize(0.03);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    it->profile->GetYaxis()->SetRangeUser(150, 200);
    it->profile->SetLineColor(it->color);
    it->profile->SetMarkerStyle(1);
    leg1->AddEntry(it->profile, it->label, "PL");
  }
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    (it == samples.begin()) ? it->profile->Draw() : it->profile->Draw("SAME");
  }

  leg1->Draw();
}
