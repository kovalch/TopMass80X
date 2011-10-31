#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "THStack.h"

TTree* tTTJets;
TTree* tWJets;
TTree* tData;

std::vector<double> lumiWeight;

enum styles             {kSigCP, kSigWP, kSigUN,  kWJets,   kData,  };
int color_      [ 5 ] = {kRed+1, kRed-7, kRed-10, kGreen-3, kBlack, };
int markerStyle_[ 5 ] = {20,     22,     22,      23,       20,     };

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
  gStyle->SetOptStat(0);

  //TFile* testFile = new TFile("./analyzeTop_1725.root");
  //TTree* testTree = testFile->Get("analyzeHitFit/eventTree"));
  
  // ---
  //    open input files
  // ---
  TFile* fTTJets = new TFile("analyzeTop.root");
  TFile* fWJets  = new TFile("WJets.root");
  TFile* fData   = new TFile("Run2011A.root");
  
  // ---
  //    Get trees
  // ---
  tTTJets = (TTree*) fTTJets->Get("analyzeHitFit/eventTree");
  tWJets  = (TTree*) fWJets ->Get("analyzeHitFit/eventTree");
  tData   = (TTree*) fData  ->Get("analyzeHitFit/eventTree");

  // ---
  // define weights concerning luminosity
  // ---
  double luminosity   = 1090;
  double crossSection = 0;
  double sampleSize   = 0;

  // 7 TeV Monte Carlo samples
  for(unsigned int idx=0; idx<4; ++idx) {
    switch(idx) {
      case kSigCP:
      case kSigWP:
      case kSigUN: {
        crossSection = 158.; // combined 2010 CMS XSec
        sampleSize   = 3701947.; // Summer11
        break;
      }
      case kWJets: {
        crossSection = 31314.; // combined 2010 CMS XSec
        sampleSize   = 14805546.; // Fall10
        break;
      }
    }
    lumiWeight.push_back(luminosity*crossSection/sampleSize);
  }
  
  makeControlPlot("hadTopMass", "m_{t}", 50, 50, 400);
  //makeControlPlot("hadBPt", "pT", 50, 0, 200);
  //makeControlPlot("hadBEta", "pT", 50, -2.5, 2.5);
}

void makeControlPlot(TString sObservable, TString sObservableShort, int nbinsx, double xlow, double xup, bool logY = false) {
  TCanvas* cPermutationWeights = new TCanvas("cPermutationWeights", "cPermutationWeights", 600, 600);
  cPermutationWeights->SetLogy(logY);
  cPermutationWeights->cd();
  
  TString sTTJetsCP = sObservable; sTTJetsCP += " >> hTTJetsCP("; sTTJetsCP += nbinsx;
    sTTJetsCP += ","; sTTJetsCP += xlow; sTTJetsCP += ","; sTTJetsCP += xup; sTTJetsCP += ")";
  TString sTTJetsWP = sObservable; sTTJetsWP += " >> hTTJetsWP("; sTTJetsWP += nbinsx;
    sTTJetsWP += ","; sTTJetsWP += xlow; sTTJetsWP += ","; sTTJetsWP += xup; sTTJetsWP += ")";
  TString sTTJetsUN = sObservable; sTTJetsUN += " >> hTTJetsUN("; sTTJetsUN += nbinsx;
    sTTJetsUN += ","; sTTJetsUN += xlow; sTTJetsUN += ","; sTTJetsUN += xup; sTTJetsUN += ")";
  TString sWJets = sObservable; sWJets += " >> hWJets("; sWJets += nbinsx;
    sWJets += ","; sWJets += xlow; sWJets += ","; sWJets += xup; sWJets += ")";
  TString sData = sObservable; sData += " >> hData("; sData += nbinsx;
    sData += ","; sData += xlow; sData += ","; sData += xup; sData += ")";

  tTTJets->Draw(sTTJetsCP, "(PUWeight)*(target==1)");
  tTTJets->Draw(sTTJetsWP, "(PUWeight)*(target==0)");
  tTTJets->Draw(sTTJetsUN, "(PUWeight)*(target==-10)");
  tWJets ->Draw(sWJets, "");
  tData  ->Draw(sData, "");
  
  TH1F* hNull = new TH1F("null", "", nbinsx, xlow, xup);
  hNull->GetYaxis()->SetRangeUser(1, hData->GetMaximum()*1.2);
  hNull->GetXaxis()->SetTitle(sObservableShort);
  
  hTTJetsCP->Scale(lumiWeight[kSigCP]);
  hTTJetsWP->Scale(lumiWeight[kSigWP]);
  hTTJetsUN->Scale(lumiWeight[kSigUN]);
  hWJets   ->Scale(lumiWeight[kWJets]);
  
  hTTJetsCP->SetFillColor(color_[kSigCP]);
  hTTJetsWP->SetFillColor(color_[kSigWP]);
  hTTJetsUN->SetFillColor(color_[kSigUN]);
  hWJets   ->SetFillColor(color_[kWJets]);
  hData    ->SetMarkerStyle(markerStyle_[kData]);
  
  THStack* stack = new THStack("stack", "");
  stack->Add(hWJets);
  stack->Add(hTTJetsCP);
  stack->Add(hTTJetsWP);
  stack->Add(hTTJetsUN);
  
  // ---
  //    create legend
  // ---
  TLegend *leg0 = new TLegend(0.6, 0.7, 0.85, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  leg0->AddEntry( hTTJetsCP, "t#bar{t} correct", "F" );
  leg0->AddEntry( hTTJetsWP, "t#bar{t} wrong", "F" );
  leg0->AddEntry( hTTJetsUN, "t#bar{t} unmatched", "F" );
  leg0->AddEntry( hWJets, "W#rightarrowl#nu", "F" );
  leg0->AddEntry( hData, "Data (1.0 fb ^{-1})", "PL");


  hNull->Draw();
  stack->Draw("SAME");
  hData->Draw("E,SAME");
  leg0->Draw("");
  
  //drawcutline(-1.3, 3000);
}
