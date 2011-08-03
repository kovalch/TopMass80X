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

void permutationWeights()
{
  TCanvas* cPermutationWeights = new TCanvas("cPermutationWeights", "cPermutationWeights", 600, 600);
  cPermutationWeights->SetLogy(1);
  cPermutationWeights->cd();
  
  gStyle->SetOptStat(0);
  
  enum styles             {kSig,   kBkg,   kBkg2,   kWjets,   kData,  };
  int color_      [ 5 ] = {kRed+1, kRed-7, kRed-10, kGreen-3, kBlack, };
  int markerStyle_[ 5 ] = {20,     22,     22,      23,       20,     };

  //TFile* testFile = new TFile("./analyzeTop_1725.root");
  //TTree* testTree = testFile->Get("analyzeHitFit/eventTree"));
  
  // ---
  //    open input files
  // ---
  TFile* fTTJets = new TFile("/scratch/hh/current/cms/user/mseidel/TTJets1725-S4_1.00_1.00/analyzeTop.root");
  TFile* fWJets  = new TFile("/scratch/hh/current/cms/user/mseidel/WJets_abs/WJets.root");
  TFile* fData   = new TFile("/scratch/hh/current/cms/user/mseidel/Run2011A.root");
  
  // ---
  //    Get trees
  // ---
  TTree* tT = (TTree*) fTTJets->Get("analyzeHitFit/eventTree");
  TTree* tW = (TTree*) fWJets ->Get("analyzeHitFit/eventTree");
  TTree* tD = (TTree*) fData  ->Get("analyzeHitFit/eventTree");

  // ---
  // define weights concerning luminosity
  // ---
  double luminosity = 1100;

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

  tT->Draw("log10(bProbSSV*hitFitProb) >> hTCP(40, -5, 0)", "(PUWeight)*(target==1)");
  tT->Draw("log10(bProbSSV*hitFitProb) >> hTWP(40, -5, 0)", "(PUWeight)*(target==0)");
  tT->Draw("log10(bProbSSV*hitFitProb) >> hTUN(40, -5, 0)", "(PUWeight)*(target==-10)");
  tW->Draw("log10(bProbSSV*hitFitProb) >> hW   (40, -5, 0)", "");
  tD->Draw("log10(bProbSSV*hitFitProb) >> hD   (40, -5, 0)", "");
  
  TH1F* hNull = new TH1F("null", "", 40, -5, 0);
  hNull->GetYaxis()->SetRangeUser(1, 50000);
  hNull->GetXaxis()->SetTitle("log(w_{i})");
  
  TH1F* hTCP = (TH1F*) gDirectory->Get("hTCP");
  TH1F* hTWP = (TH1F*) gDirectory->Get("hTWP");
  TH1F* hTUN = (TH1F*) gDirectory->Get("hTUN");
  TH1F* hW   = (TH1F*) gDirectory->Get("hW");
  TH1F* hD   = (TH1F*) gDirectory->Get("hD");
  
  hTCP->Scale(1/normalization * 188/35.9*luminosity);
  hTWP->Scale(1/normalization * 188/35.9*luminosity);
  hTUN->Scale(1/normalization * 188/35.9*luminosity);
  hW  ->Scale(1/hW->Integral() *  21/35.9*luminosity);
  
  hTCP->SetFillColor(color_[kSig]);
  hTWP->SetFillColor(color_[kBkg]);
  hTUN->SetFillColor(color_[kBkg2]);
  hW  ->SetFillColor(color_[kWjets]);
  hD  ->SetMarkerStyle(markerStyle_[kData]);
  
  THStack* stack = new THStack("stack", "");
  stack->Add(hW);
  stack->Add(hTCP);
  stack->Add(hTWP);
  stack->Add(hTUN);
  
  // ---
  //    create legend
  // ---
  TLegend *leg0 = new TLegend(0.6, 0.7, 0.85, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  leg0->AddEntry( hTCP , "t#bar{t} correct", "F" );
  leg0->AddEntry( hTWP , "t#bar{t} wrong", "F" );
  leg0->AddEntry( hTUN , "t#bar{t} unmatched", "F" );
  leg0->AddEntry( hW , "W#rightarrowl#nu", "F" );
  leg0->AddEntry( hD , "Data (1.0 fb ^{-1})", "PL");


  hNull->Draw();
  stack->Draw("SAME");
  hD   ->Draw("E,SAME");
  leg0->Draw("");
  
  drawcutline(-1.3, 3000);
}
