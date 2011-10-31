#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "THStack.h"

#include "tdrstyle.C"

TTree* tTTJets;
TTree* tWJets;
TTree* tSAntiTops;
TTree* tSAntiTopt;
TTree* tSAntiToptW;
TTree* tSTops;
TTree* tSTopt;
TTree* tSToptW;
TTree* tData;

std::vector<double> lumiWeight;

enum samples    {kSigCP, kSigWP, kSigUN,  kZjets  , kWjets  , kQCD   , kSTop   , kDiBos, kData , kWW, kWZ, kZZ, kSAntiTops, kSAntiTopt, kSAntiToptW, kSTops  , kSTopt  , kSToptW };
int color_ [] = {kRed+1, kRed-7, kRed-10, kAzure-2, kGreen-3, kYellow, kMagenta, 10    , kBlack, 10 , 10 , 10 , kMagenta  , kMagenta  , kMagenta   , kMagenta, kMagenta, kMagenta};
int marker_[] = {20,     22,     22,      29      , 23      , 21     , 27      , 28    , 20    , 28 , 28 , 28 , 27        , 27        , 27         , 27      , 27      , 27      };


void drawcutline(double cutval, double maximum)
{
  TLine *cut = new TLine();
  cut->SetLineWidth(2);
  cut->SetLineStyle(7);
  cut->SetLineColor(1);
  cut->DrawLine(cutval, 0.,cutval, maximum);
}

void controlPlots()
{
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
	gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);

  //TFile* testFile = new TFile("./analyzeTop_1725.root");
  //TTree* testTree = testFile->Get("analyzeHitFit/eventTree"));
  
  // ---
  //    open input files
  // ---
  TFile* fTTJets = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_1.00_2b/analyzeTop.root");
  TFile* fWJets  = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_WJets_1.00_1.00_2b/analyzeTop.root");
  TFile* fSAntiTops  = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_Tbar_s-channel_1.00_1.00_2b/analyzeTop.root");
  TFile* fSAntiTopt  = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_Tbar_t-channel_1.00_1.00_2b/analyzeTop.root");
  TFile* fSAntiToptW = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_Tbar_tW-channel_1.00_1.00_2b/analyzeTop.root");
  TFile* fSTops  = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_T_s-channel_1.00_1.00_2b/analyzeTop.root");
  TFile* fSTopt  = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_T_t-channel_1.00_1.00_2b/analyzeTop.root");
  TFile* fSToptW = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_T_tW-channel_1.00_1.00_2b/analyzeTop.root");
  TFile* fData   = new TFile("/scratch/hh/current/cms/user/mseidel/Run2011A_2b.root");
  
  // ---
  //    Get trees
  // ---
  tTTJets = (TTree*) fTTJets->Get("analyzeHitFit/eventTree");
  tWJets  = (TTree*) fWJets ->Get("analyzeHitFit/eventTree");
  tSAntiTops  = (TTree*) fSAntiTops ->Get("analyzeHitFit/eventTree");
  tSAntiTopt  = (TTree*) fSAntiTopt ->Get("analyzeHitFit/eventTree");
  tSAntiToptW = (TTree*) fSAntiToptW->Get("analyzeHitFit/eventTree");
  tSTops  = (TTree*) fSTops ->Get("analyzeHitFit/eventTree");
  tSTopt  = (TTree*) fSTopt ->Get("analyzeHitFit/eventTree");
  tSToptW = (TTree*) fSToptW->Get("analyzeHitFit/eventTree");
  
  tData   = (TTree*) fData  ->Get("analyzeHitFit/eventTree");

  // ---
  // define weights concerning luminosity
  // ---
  double luminosity   = 1132;
  double crossSection = 0;
  double sampleSize   = 0;

  // 7 TeV Monte Carlo samples
  for(unsigned int idx=0; idx<18; ++idx) {
    switch(idx) {
      case kSigCP:
      case kSigWP:
      case kSigUN: {
        crossSection = 158.; // combined 2010 CMS XSec
        sampleSize   = 3701947.; // Summer11
        break;
      }
      case kWjets: {
        crossSection = 31314.;
        sampleSize   = 81352581.; // Summer11
        break;
      }
      case kSAntiTops: {
        crossSection = 1.44; //
        sampleSize   = 137980.; // Summer11
        break;
      }
      case kSAntiTopt: {
        crossSection = 22.65; //
        sampleSize   = 1944826.; // Summer11
        break;
      }
      case kSAntiToptW: {
        crossSection = 7.87; //
        sampleSize   = 809984.; // Summer11
        break;
      }
      case kSTops: {
        crossSection = 3.19; //
        sampleSize   = 259971.; // Summer11
        break;
      }
      case kSTopt: {
        crossSection = 41.92; //
        sampleSize   = 3900171.; // Summer11
        break;
      }
      case kSToptW: {
        crossSection = 7.87; //
        sampleSize   = 814390.; // Summer11
        break;
      }
      default: {
        crossSection = 0.; // 
        sampleSize   = 1.; // Dummy
      }
    }
    lumiWeight.push_back(luminosity*crossSection/sampleSize);
  }
  
  //* Raw
  makeControlPlot("hadWRawMass", "m_{W,had}", "hadWRawMass", 50, 0, 300);
  makeControlPlot("lepWRawMass", "m_{W,lep}", "lepWRawMass", 50, 0, 300);
  makeControlPlot("hadBBSSV", "b-disc (SSVHE)", "hadBBSSV", 30, 1, 7, true);
  makeControlPlot("hadTopRawMass", "m_{t,had}", "hadTopRawMass", 50, 50, 400);
  makeControlPlot("lepTopRawMass", "m_{t,lep}", "lepTopRawMass", 50, 50, 400);
  makeControlPlot("jetMultiplicity", "Number of jets", "jetMultiplicity", 5, 4, 9);
  makeControlPlot("nVertex", "Number of vertices", "nVertex", 16, 0, 16);
  //*/
  
  /*
  makeControlPlot("hadTopMass", "m_{t}", "hadTopMass", 50, 50, 400);
  makeControlPlot("hadBBSSV", "B(SSVHE)", "hadBBSSV", 50, -1, 5, true);
  makeControlPlot("log10(bProbSSV*hitFitProb)", "log_{10}(w_{i})", "weightFactor", 50, -5, 0, true, -1.3);
  */
  
  /* Fitted
  makeControlPlot("hadBPt", "p_{T,b}", "BPt", 20, 0, 200);
  makeControlPlot("hadBEta", "#eta_{b}", "BEta", 25, -2.5, 2.5);
  makeControlPlot("hadQPt", "p_{T,q}", "QPt", 20, 0, 200);
  makeControlPlot("hadQEta", "#eta_{q}", "QEta", 25, -2.5, 2.5);
  
  makeControlPlot("hadTopMass", "m_{t}", "hadTopMass", 50, 50, 400);
  makeControlPlot("hadTopPt", "p_{T,t,had}", "hadTopPt", 20, 0, 200);
  makeControlPlot("hadTopEta", "#eta_{t,had}", "hadTopEta", 25, -2.5, 2.5);
  makeControlPlot("TTBarMass", "m_{tt}", "TTBarMass", 50, 100, 1300);
  makeControlPlot("hitFitProb", "P_{fit}", "hitFitProb", 20, 0, 1, true, 0.2);
  makeControlPlot("hitFitChi2", "#chi^{2}_{fit}", "hitFitChi2", 20, 0, 10, false, 3.218875825);
  //*/
  
  /* cut&weight
  makeControlPlot("hadBPt", "p_{T,b}", "BPt_weighted", 20, 0, 200, false, -999, true);
  makeControlPlot("hadBEta", "#eta_{b}", "BEta_weighted", 25, -2.5, 2.5, false, -999, true);
  makeControlPlot("hadQPt", "p_{T,q}", "QPt_weighted", 20, 0, 200, false, -999, true);
  makeControlPlot("hadQEta", "#eta_{q}", "QEta_weighted", 25, -2.5, 2.5, false, -999, true);
  
  makeControlPlot("hadTopMass", "m_{t}", "hadTopMass_weighted", 50, 50, 400, false, -999, true);
  makeControlPlot("hadWRawMass", "m_{W}", "hadWRawMass_weighted", 50, 0, 300, false, -999, true);
  makeControlPlot("hadTopPt", "p_{T,t,had}", "hadTopPt_weighted", 20, 0, 200, false, -999, true);
  makeControlPlot("hadTopEta", "#eta_{t,had}", "hadTopEta_weighted", 25, -2.5, 2.5, false, -999, true);
  makeControlPlot("TTBarMass", "m_{tt}", "TTBarMass_weighted", 50, 100, 1300, false, -999, true);
  //*/

}

void makeControlPlot(TString sObservable, TString sObservableShort, TString sFileName, int nbinsx, double xlow, double xup, bool logY = false, double cut = -999, bool weighted = false) {
  TCanvas* cControlPlots = new TCanvas("cControlPlots", "cControlPlots", 600, 600);
  cControlPlots->SetLogy(logY);
  cControlPlots->cd();
  
  TString sTTJetsCP = sObservable; sTTJetsCP += " >> hTTJetsCP("; sTTJetsCP += nbinsx;
    sTTJetsCP += ","; sTTJetsCP += xlow; sTTJetsCP += ","; sTTJetsCP += xup; sTTJetsCP += ")";
  TString sTTJetsWP = sObservable; sTTJetsWP += " >> hTTJetsWP("; sTTJetsWP += nbinsx;
    sTTJetsWP += ","; sTTJetsWP += xlow; sTTJetsWP += ","; sTTJetsWP += xup; sTTJetsWP += ")";
  TString sTTJetsUN = sObservable; sTTJetsUN += " >> hTTJetsUN("; sTTJetsUN += nbinsx;
    sTTJetsUN += ","; sTTJetsUN += xlow; sTTJetsUN += ","; sTTJetsUN += xup; sTTJetsUN += ")";
  
  TString sWJets = sObservable; sWJets += " >> hWJets("; sWJets += nbinsx;
    sWJets += ","; sWJets += xlow; sWJets += ","; sWJets += xup; sWJets += ")";
  
  TString sSAntiTops = sObservable; sSAntiTops += " >> hSAntiTops("; sSAntiTops += nbinsx;
    sSAntiTops += ","; sSAntiTops += xlow; sSAntiTops += ","; sSAntiTops += xup; sSAntiTops += ")";
  TString sSAntiTopt = sObservable; sSAntiTopt += " >> hSAntiTopt("; sSAntiTopt += nbinsx;
    sSAntiTopt += ","; sSAntiTopt += xlow; sSAntiTopt += ","; sSAntiTopt += xup; sSAntiTopt += ")";
  TString sSAntiToptW = sObservable; sSAntiToptW += " >> hSAntiToptW("; sSAntiToptW += nbinsx;
    sSAntiToptW += ","; sSAntiToptW += xlow; sSAntiToptW += ","; sSAntiToptW += xup; sSAntiToptW += ")";
  TString sSTops = sObservable; sSTops += " >> hSTops("; sSTops += nbinsx;
    sSTops += ","; sSTops += xlow; sSTops += ","; sSTops += xup; sSTops += ")";
  TString sSTopt = sObservable; sSTopt += " >> hSTopt("; sSTopt += nbinsx;
    sSTopt += ","; sSTopt += xlow; sSTopt += ","; sSTopt += xup; sSTopt += ")";
  TString sSToptW = sObservable; sSToptW += " >> hSToptW("; sSToptW += nbinsx;
    sSToptW += ","; sSToptW += xlow; sSToptW += ","; sSToptW += xup; sSToptW += ")";
    
  TString sData = sObservable; sData += " >> hData("; sData += nbinsx;
    sData += ","; sData += xlow; sData += ","; sData += xup; sData += ")";
  
  if (!weighted) {
    tTTJets->Draw(sTTJetsCP, "(MCWeight)*(target==1 & combi==0)");
    tTTJets->Draw(sTTJetsWP, "(MCWeight)*(target==0 & combi==0)");
    tTTJets->Draw(sTTJetsUN, "(MCWeight)*(target==-10 & combi==0)");
    tWJets ->Draw(sWJets, "(MCWeight)*(combi==0)");
    tSAntiTops ->Draw(sSAntiTops,  "(MCWeight)*(combi==0)");
    tSAntiTopt ->Draw(sSAntiTopt,  "(MCWeight)*(combi==0)");
    tSAntiToptW->Draw(sSAntiToptW, "(MCWeight)*(combi==0)");
    tSTops ->Draw(sSTops,  "(MCWeight)*(combi==0)");
    tSTopt ->Draw(sSTopt,  "(MCWeight)*(combi==0)");
    tSToptW->Draw(sSToptW, "(MCWeight)*(combi==0)");
    tData  ->Draw(sData, "hadQBSSV<1.74 & hadQBarBSSV<1.74 & hadBBSSV>1.74 & lepBBSSV>1.74 & run<168000 & combi==0");
  }
  else {
    tTTJets->Draw(sTTJetsCP, "(MCWeight*hitFitProb)*(target==1 & hitFitProb>0.2)");
    tTTJets->Draw(sTTJetsWP, "(MCWeight*hitFitProb)*(target==0 & hitFitProb>0.2)");
    tTTJets->Draw(sTTJetsUN, "(MCWeight*hitFitProb)*(target==-10 & hitFitProb>0.2)");
    tWJets ->Draw(sWJets, "(MCWeight*hitFitProb)*(hitFitProb>0.2)");
    tSAntiTops ->Draw(sSAntiTops,  "(MCWeight*hitFitProb)*(hitFitProb>0.2)");
    tSAntiTopt ->Draw(sSAntiTopt,  "(MCWeight*hitFitProb)*(hitFitProb>0.2)");
    tSAntiToptW->Draw(sSAntiToptW, "(MCWeight*hitFitProb)*(hitFitProb>0.2)");
    tSTops ->Draw(sSTops,  "(MCWeight*hitFitProb)*(hitFitProb>0.2)");
    tSTopt ->Draw(sSTopt,  "(MCWeight*hitFitProb)*(hitFitProb>0.2)");
    tSToptW->Draw(sSToptW, "(MCWeight*hitFitProb)*(hitFitProb>0.2)");
    tData  ->Draw(sData, "(hitFitProb)*(hadQBSSV<1.74 & hadQBarBSSV<1.74 & hadBBSSV>1.74 & lepBBSSV>1.74 & hitFitProb>0.2 & run<168000)");
  }
  
  TH1F* hNull = new TH1F("null", "", nbinsx, xlow, xup);
  
  if (logY) hNull->GetYaxis()->SetRangeUser(1, hData->GetMaximum()*10);
  else hNull->GetYaxis()->SetRangeUser(1, hData->GetMaximum()*1.5);
  
  hNull->GetXaxis()->SetTitle(sObservableShort);
  
	if (weighted) hNull->GetYaxis()->SetTitle("Sum of weights");
	else hNull->GetYaxis()->SetTitle("Number of permutations");
  
  hTTJetsCP->Scale(lumiWeight[kSigCP]);
  hTTJetsWP->Scale(lumiWeight[kSigWP]);
  hTTJetsUN->Scale(lumiWeight[kSigUN]);
  hWJets   ->Scale(lumiWeight[kWjets]);
  hSAntiTops   ->Scale(lumiWeight[kSAntiTops]);
  hSAntiTopt   ->Scale(lumiWeight[kSAntiTopt]);
  hSAntiToptW  ->Scale(lumiWeight[kSAntiToptW]);
  hSTops   ->Scale(lumiWeight[kSTops]);
  hSTopt   ->Scale(lumiWeight[kSTopt]);
  hSToptW  ->Scale(lumiWeight[kSToptW]);
    
  TH1F* hSTop = new TH1F("hSTop", "", nbinsx, xlow, xup);
  for (int i = 0; i < nbinsx+2; i++) {
    hSTop->SetBinContent(i, hSAntiTops ->GetBinContent(i)
                          + hSAntiTopt ->GetBinContent(i)
                          + hSAntiToptW->GetBinContent(i)
                          + hSTops ->GetBinContent(i)
                          + hSTopt ->GetBinContent(i)
                          + hSToptW->GetBinContent(i)
                        );
  }
  
  hTTJetsCP->SetFillColor(color_[kSigCP]);
  hTTJetsWP->SetFillColor(color_[kSigWP]);
  hTTJetsUN->SetFillColor(color_[kSigUN]);
  hSTop    ->SetFillColor(color_[kSTop]);
  hWJets   ->SetFillColor(color_[kWjets]);
  hData    ->SetMarkerStyle(marker_[kData]);
  
  THStack* stack = new THStack("stack", "");
  stack->Add(hSTop);
  stack->Add(hWJets);
  stack->Add(hTTJetsCP);
  stack->Add(hTTJetsWP);
  stack->Add(hTTJetsUN);
  
  // ---
  //    create legend
  // ---
  TLegend *leg0 = new TLegend(0.25, 0.75, 0.55, 0.925);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  leg0->AddEntry( hTTJetsCP, "t#bar{t} correct", "F" );
  leg0->AddEntry( hTTJetsWP, "t#bar{t} wrong", "F" );
  leg0->AddEntry( hTTJetsUN, "t#bar{t} unmatched", "F" );
  
  TLegend *leg1 = new TLegend(0.6, 0.75, 0.9, 0.925);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->AddEntry( hWJets, "W#rightarrowl#nu", "F" );
  leg1->AddEntry( hSTop, "single top", "F" );
  leg1->AddEntry( hData, "Data (1.1 fb ^{-1})", "PL");


  hNull->Draw();
  stack->Draw("SAME");
  hData->Draw("E,SAME");
  leg0->Draw("");
  leg1->Draw();
  
  if (cut != -999) drawcutline(cut, hData->GetMaximum());
  
  TString path("controlPlots/"); path+= sFileName; path += ".eps";
  cControlPlots->Print(path);
}
