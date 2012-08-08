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

enum lepton           { kElectron, kMuon, kAll};
TString lepton_ [3] = { "electron", "muon", "all"};

int channel = 2;

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
  TGraphErrors* graph;
  const char* label;
  int color;
  sample(TChain* ch, const char* l, int c = kBlack)
  : chain(ch), label(l), color(c) {}
};


std::vector<sample> samples;
std::vector<sample>::iterator it;

TH1F* responseProfile(TString sDraw, TChain* chain, TString label) {
  sDraw += " >> h2"; sDraw += label; sDraw += sDraw(1, 7); sDraw += "(25, 0, 5)";
  chain->Draw(sDraw);//, "(target==1 & leptonPt > 30 & hadQBCSV<0.679 & hadQBarBCSV<0.679 & hadBBCSV>0.679 & lepBBCSV>0.679 & hitFitProb>0.2)*(MCWeight*hitFitProb)");
  TString sH2("h2"); sH2 += label; sH2 += sDraw(1, 7);
  TH1F* h2 = (TH1F*) gDirectory->Get(sH2);
  h2->Fit("gaus");
  return h2;
}

TGraphErrors* ratioGraph(TH1F* p1, TH1F* p2, int color) {
  int n = 25;
  double position[n];
  double posError[n];
  double ratio[n];
  double ratioError[n];
  
  for (int i = 0; i < n; ++i) {
    position[i] = 0.5 + 1*i;
    posError[i] = 0.5;
    
    if (p2->GetBinContent(i+1)>0) {
      ratio[i]      = p1->GetBinContent(i+1)/p2->GetBinContent(i+1) * p2->Integral()/p1->Integral();
      ratioError[i] = sqrt(pow(p1->GetBinError(i+1)/p2->GetBinContent(i+1), 2) + pow(p1->GetBinContent(i+1)/pow(p2->GetBinContent(i+1), 2) * p2->GetBinError(i+1), 2))  * p2->Integral()/p1->Integral();
    }
    /*
    else {
      ratio[i] = 1;
      ratioError[i] = 0.00001;
    }
    */
  }
  
  TGraphErrors* graph = new TGraphErrors(n, position, ratio, posError, ratioError);
  
  TF1* constFit = new TF1("constFit", "[0]+x*[1]");
  constFit->SetParNames("offset");
  constFit->SetLineColor(color);
  constFit->SetLineWidth(1);
  constFit->SetLineStyle(7);
  graph->Fit("constFit");
  
  return graph;
}

void massResponse()
{
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  
  // Define chains
  TString sTree("analyzeGenMatch/eventTree");
  tTTJets096        = new TChain(sTree);
  tTTJets100        = new TChain(sTree);
  tTTJets104        = new TChain(sTree);
  tTTJetsFlavorUp   = new TChain(sTree);
  tTTJetsFlavorDown = new TChain(sTree);
  tTTJetsP11        = new TChain(sTree);
  tTTJetsP11noCR    = new TChain(sTree);
  
  
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
  }
  if (channel == kElectron || channel == kAll) {
    tTTJets096        ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_0.96_electron/analyzeTop.root");
    tTTJets100        ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_1.00_electron/analyzeTop.root");
    tTTJets104        ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_1.04_electron/analyzeTop.root");
    tTTJetsFlavorUp   ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_flavor:up_electron/analyzeTop.root");
    tTTJetsFlavorDown ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_flavor:down_electron/analyzeTop.root");
    tTTJetsP11        ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_P11_electron/analyzeTop.root");
    tTTJetsP11noCR    ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_P11noCR_electron/analyzeTop.root");
  }
  
  samples.push_back(sample(tTTJets100, "Z2", kBlue+1));
  samples.push_back(sample(tTTJets096, "Z2_jesDown", kRed+1));
  samples.push_back(sample(tTTJets104, "Z2_jesUp", kGreen+1));
  //samples.push_back(sample(tTTJetsFlavorDown, "Z2_bDown", kRed));
  //samples.push_back(sample(tTTJetsFlavorUp, "Z2_bUp", kGreen));
  samples.push_back(sample(tTTJetsP11, "P11", kMagenta+1));
  samples.push_back(sample(tTTJetsP11noCR, "P11noCR", kYellow+1));
  
  //"[0] + exp([1]+[2]*x)"
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    
    //it->profile  = responseProfile("hadTopMass", it->chain, it->label);
    it->profile  = responseProfile("deltaThetaHadQHadQBar", it->chain, it->label);
    //it->profile  = responseProfile("(hadQRaw.Pt()+hadQBarRaw.Pt())/(hadQGen.Pt()+hadQBarGen.Pt()):(hadQGen.Pt()+hadQBarGen.Pt())/2", it->chain, it->label);
    //it->profileB = responseProfile("(lepBRaw.Pt()/lepBGen.Pt() + hadBRaw.Pt()/hadBGen.Pt())/2:(hadBGen.Pt()+lepBGen.Pt())/2", it->chain, it->label);
    
    //it->profile  = responseProfile("hadQRaw.Pt()/hadQGen.Pt():hadQGen.Pt()", it->chain, it->label);
    //it->profile  = responseProfile("hadWRaw.E()/hadWGen.E():hadWGen.Pt()", it->chain, it->label);
    //it->profileB = responseProfile("(hadBRaw.Pt()/hadBGen.Pt()):hadBGen.Pt()", it->chain, it->label);
    //it->profileB = responseProfile("(lepBRaw.Pt()/lepBGen.Pt()):lepBGen.Pt()", it->chain, it->label);
  }
  
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle(";[bin];mass response rel. to Z2");
  
  TLegend *leg1 = new TLegend(0.6, 0.75, 0.9, 0.925);
  leg1->SetTextSize(0.03);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    it->graph = ratioGraph(it->profile, samples.begin()->profile, it->color);
    it->graph->SetLineColor(it->color);
    it->graph->SetMarkerStyle(1);
    mg->Add(it->graph);
    leg1->AddEntry(it->graph, it->label, "PL");
  }
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
  
  mg->Draw("AP");
  //*
  mg->SetMinimum(0.);
  mg->SetMaximum(2);
  //*/
  leg1->Draw();
}
