#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "THStack.h"

void drawcutline(double cutval, double maximum)
{
  TLine *cut = new TLine();
  cut->SetLineWidth(2);
  cut->SetLineStyle(7);
  cut->SetLineColor(1);
  cut->DrawLine(cutval, 0.,cutval, maximum);
}

void ideogramWeightCut()
{
  TCanvas* canvas = new TCanvas("canvas", "canvas");
  canvas->SetLogy(1);
  canvas->cd();
  
  gStyle->SetOptStat(0);
  
  enum styles             {kSig,   kBkg,   kWjets,   kData,  };
  int color_      [ 4 ] = {kRed+1, kRed-7, kGreen-3, kBlack, };
  int markerStyle_[ 4 ] = {20,     22,     23,       20,     };

  //TFile* testFile = new TFile("./analyzeTop_1725.root");
  //TTree* testTree = testFile->Get("analyzeHitFit/eventTree"));
  
  // ---
  //    open input files
  // ---
  TFile* fTTJets = new TFile("./analyzeTop_1725.root");
  TFile* fWJets  = new TFile("./WJets.root");
  TFile* fData   = new TFile("./analyzeTop_Run2010B.root");
  
  // ---
  //    Get trees
  // ---
  TTree* tT = (TTree*) fTTJets->Get("analyzeHitFit/eventTree");
  TTree* tW = (TTree*) fWJets ->Get("analyzeHitFit/eventTree");
  TTree* tD = (TTree*) fData  ->Get("analyzeHitFit/eventTree");

  // ---
  // define weights concerning luminosity
  // ---
  double luminosity = 23;

  double normalization = tT->GetEntries("combi==0");
  
  std::vector<double> lumiweight_;
  // add scaling factors here!

  // 7 TeV Monte Carlo fall 10 samples
  // -----------------------------------
  for(unsigned int idx=0; idx<3; ++idx) {
    // for current ttbar(lept.mu on gen level and other) Madgraph 
    if (idx==kSig || idx==kBkg) lumiweight_.push_back(188./35.9*luminosity);
    // for current W+jets MADGRAPH sample
    if (idx==kWjets) lumiweight_.push_back(21./35.9*luminosity);
  }

  for (unsigned int idx=0; idx<3; ++idx) {
    //trees_[idx]->Draw("bProbSSV*hitFitProb >> prob(150, -100, 50)");
    
    //TH1F *prob = (TH1F*) gDirectory->Get("prob");
    
    //prob->Draw();
  }
  
  tT->Draw("log10(bProbSSV*hitFitProb) >> hTSig(40, -10, 0)", "target==1");
  tT->Draw("log10(bProbSSV*hitFitProb) >> hTBkg(40, -10, 0)", "target!=1");
  tW->Draw("log10(bProbSSV*hitFitProb) >> hW   (40, -10, 0)");
  tD->Draw("log10(bProbSSV*hitFitProb) >> hD   (40, -10, 0)");
  
  TH1F* hNull = new TH1F("null", "", 40, -10, 0);
  hNull->GetYaxis()->SetRangeUser(0.1, 1000);
  hNull->GetXaxis()->SetTitle("log(w_{i})");
  
  TH1F* hTSig = (TH1F*) gDirectory->Get("hTSig");
  TH1F* hTBkg = (TH1F*) gDirectory->Get("hTBkg");
  TH1F* hW    = (TH1F*) gDirectory->Get("hW");
  TH1F* hD    = (TH1F*) gDirectory->Get("hD");
  
  hTBkg->Scale(1/normalization * 188/35.9*luminosity);
  hTSig->Scale(1/normalization * 188/35.9*luminosity);
  hW   ->Scale(1/hW->Integral() *  21/35.9*luminosity);
  
  hTSig->SetFillColor(color_[kSig]);
  hTBkg->SetFillColor(color_[kBkg]);
  hW   ->SetFillColor(color_[kWjets]);
  hD   ->SetMarkerStyle(markerStyle_[kData]);
  
  THStack* stack = new THStack("stack", "");
  stack->Add(hW);
  stack->Add(hTSig);
  stack->Add(hTBkg);
  
  // ---
  //    create legend
  // ---
  // samples: separate canvas
  TLegend *leg0 = new TLegend(0.6, 0.7, 0.85, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  leg0->AddEntry( hD , "Data (23 pb ^{-1})"       , "PL");
  leg0->AddEntry( hW , "W#rightarrowl#nu"              , "F" );
  leg0->AddEntry( hTBkg , "t#bar{t} wp"                , "F" );
  leg0->AddEntry( hTSig , "t#bar{t} cp"  , "F" );

  hNull->Draw();
  stack->Draw("SAME");
  hD   ->Draw("E,SAME");
  leg0->Draw("");
  
  drawcutline(-1.3, 100);
}
