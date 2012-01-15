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
TTree* tTTJetsCP;
TTree* tTTJetsWP;
TTree* tTTJetsUN;
TTree* tQCD;
TTree* tWJets;
TTree* tZJets;
TTree* tSAntiTops;
TTree* tSAntiTopt;
TTree* tSAntiToptW;
TTree* tSTops;
TTree* tSTopt;
TTree* tSToptW;
TTree* tData;

std::vector<double> lumiWeight;

enum plotType   {kEvent, kEvent2b, kPerm2b, kPerm2bW};

enum samples    {kSigCP, kSigWP, kSigUN,  kSig  , kZjets  , kWjets  , kQCD   , kSTop   , kDiBos, kData , kWW, kWZ, kZZ, kSAntiTops, kSAntiTopt, kSAntiToptW, kSTops  , kSTopt  , kSToptW };
int color_ [] = {kRed+1, kRed-7, kRed-10, kGray , kAzure-2, kGreen-3, kYellow, kMagenta, 10    , kBlack, 10 , 10 , 10 , kMagenta  , kMagenta  , kMagenta   , kMagenta, kMagenta, kMagenta};
int marker_[] = {20,     22,     22,      20,     29      , 23      , 21     , 27      , 28    , 20    , 28 , 28 , 28 , 27        , 27        , 27         , 27      , 27      , 27      };


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
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetNdivisions(505, "X");
	gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  
  //TH1::SetDefaultSumw2(true);

  //TFile* testFile = new TFile("./analyzeTop_1725.root");
  //TTree* testTree = testFile->Get("analyzeHitFit/eventTree"));
  
  // ---
  //    open input files
  // ---
  TFile* fTTJets = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_1.00/analyzeTop.root");
  TFile* fQCD    = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_QCD/analyzeTop.root");
  TFile* fWJets  = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_WJets/analyzeTop.root");
  TFile* fZJets  = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_ZJets/analyzeTop.root");
  TFile* fSAntiTops  = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_Tbar_s-channel/analyzeTop.root");
  TFile* fSAntiTopt  = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_Tbar_t-channel/analyzeTop.root");
  TFile* fSAntiToptW = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_Tbar_tW-channel/analyzeTop.root");
  TFile* fSTops  = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_T_s-channel/analyzeTop.root");
  TFile* fSTopt  = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_T_t-channel/analyzeTop.root");
  TFile* fSToptW = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_T_tW-channel/analyzeTop.root");
  //TFile* fData   = new TFile("/scratch/hh/current/cms/user/mseidel/Run2011_2b.root");
  TFile* fData   = new TFile("/scratch/hh/current/cms/user/mseidel/Run2011.root");
  
  // ---
  //    Get trees
  // ---
  tTTJets = (TTree*) fTTJets->Get("analyzeHitFit/eventTree");
  tQCD    = (TTree*) fQCD   ->Get("analyzeHitFit/eventTree");
  tWJets  = (TTree*) fWJets ->Get("analyzeHitFit/eventTree");
  tZJets  = (TTree*) fZJets ->Get("analyzeHitFit/eventTree");
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
  //double luminosity   = 1132;
  double luminosity   = 4700;
  //double luminosity   = 47.4;
  double crossSection = 0;
  double sampleSize   = 0;

  // 7 TeV Monte Carlo samples
  for(unsigned int idx=0; idx<19; ++idx) {
    switch(idx) {
      case kSigCP:
      case kSigWP:
      case kSigUN: {
        crossSection = 158.; // combined 2010 CMS XSec
        sampleSize   = 3701947.; // Summer11
        break;
      }
      case kSig: {
        crossSection = 158.; // combined 2010 CMS XSec + uncertainty
        sampleSize   = 3701947.; // Summer11
        break;
      }
      case kQCD: {
        crossSection = 296600000.*0.0002855; // generator crossSection * prefilter efficiency
        sampleSize   = 20416038.; // Summer11
        break;
      }
      case kWjets: {
        crossSection = 31314.;
        sampleSize   = 81352581.; // Summer11
        break;
      }
      case kZjets: {
        crossSection = 3048.;
        sampleSize   = 35032553.; // Summer11
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
  
  
  //* Test
  //makeControlPlot("event", "jetMultiplicity", "Number of jets", "jetMultiplicity", 5, 4, 9, kEvent2b);
  //makeControlPlot("event", "bottomSSVJetMultiplicity", "Number of b-jets", "bottomSSVJetMultiplicity", 5, 0, 5, kEvent);
  //makeControlPlot("event", "nuRawPt", "MET [GeV]", "nuRawPt_event", 50, 0, 500, kEvent2b, true);
  makeControlPlot("event", "nVertex", "Number of vertices", "nVertex_event", 25, 0, 25, kEvent2b);
  //makeControlPlot("test", "hadTopMass", "m_{t}^{fit} [GeV]", "hadTopMass", 35, 50, 400, kEvent2b);
  //makeControlPlot("permutation", "hitFitProb", "P_{fit}", "hitFitProb", 20, 0, 1, kPerm2b, true, 0.2);
  //makeControlPlot("permutation", "hitFitChi2", "#chi^{2}_{fit}", "hitFitChi2", 20, 0, 10, kPerm2b, false, 3.218875825);
  //makeControlPlot("permutation", "hadTopRawMass", "m_{t,had}^{raw} [GeV]", "hadTopRawMass", 35, 50, 400, kPerm2b);
  //*/
  
  /* Event
  makeControlPlot("event", "jetMultiplicity", "Number of jets", "jetMultiplicity_event", 5, 4, 9, kEvent2b);
  makeControlPlot("event", "nVertex", "Number of vertices", "nVertex_event", 25, 0, 25, kEvent2b);
  
  makeControlPlot("muon", "leptonPt", "p_{T,#mu} [GeV]", "leptonPt_event", 50, 0, 200, kEvent2b);
  makeControlPlot("muon", "leptonEta", "#eta_{#mu}", "leptonEta_event", 25, -2.5, 2.5, kEvent2b);
  makeControlPlot("event", "nuRawPt", "MET [GeV]", "nuRawPt_event", 50, 0, 200, kEvent2b);
  
  makeControlPlot("jet", "hadBRawPt", "p_{T,b} [GeV]", "BRawPt_event", 20, 0, 200, kEvent2b);
  makeControlPlot("jet", "hadQRawPt", "p_{T,q} [GeV]", "QRawPt_event", 20, 0, 200, kEvent2b);
  makeControlPlot("jet", "hadQBarRawPt", "p_{T,#bar{q}} [GeV]", "QBarRawPt_event", 20, 0, 200, kEvent2b);
  makeControlPlot("jet", "hadBEta", "#eta_{b}", "BEta_event", 25, -2.5, 2.5, kEvent2b);
  makeControlPlot("jet", "hadQEta", "#eta_{q}", "QEta_event", 25, -2.5, 2.5, kEvent2b);
  makeControlPlot("jet", "hadBBSSV", "b-disc (SSVHE)", "hadBBSSV_event", 30, 1, 7, kEvent2b);
  //*/
  
  /* Permutation
  makeControlPlot("permutation", "hadWRawMass", "m_{W,had}^{raw} [GeV]", "hadWRawMass", 30, 0, 300, kPerm2b);
  makeControlPlot("permutation", "lepWRawMass", "m_{W,lep}^{raw} [GeV]", "lepWRawMass", 30, 0, 300, kPerm2b);
  makeControlPlot("permutation", "hadTopRawMass", "m_{t,had}^{raw} [GeV]", "hadTopRawMass", 35, 50, 400, kPerm2b);
  makeControlPlot("permutation", "lepTopRawMass", "m_{t,lep}^{raw} [GeV]", "lepTopRawMass", 35, 50, 400, kPerm2b);
  makeControlPlot("permutation", "hadTopMass", "m_{t}^{fit} [GeV]", "hadTopMass", 35, 50, 400, kPerm2b);
  makeControlPlot("permutation", "hitFitProb", "P_{fit}", "hitFitProb", 20, 0, 1, kPerm2b, true, 0.2);
  makeControlPlot("permutation", "hitFitChi2", "#chi^{2}_{fit}", "hitFitChi2", 20, 0, 10, kPerm2b, false, 3.218875825);
  //*/
  
  /* Permutation, weighted
  makeControlPlot("permutation", "hadWRawMass", "m_{W,had}^{raw} [GeV]", "hadWRawMass_weighted", 30, 0, 300, kPerm2bW);
  makeControlPlot("permutation", "lepWRawMass", "m_{W,lep}^{raw} [GeV]", "lepWRawMass_weighted", 30, 0, 300, kPerm2bW);
  makeControlPlot("permutation", "hadTopRawMass", "m_{t,had}^{raw} [GeV]", "hadTopRawMass_weighted", 35, 50, 400, kPerm2bW);
  makeControlPlot("permutation", "lepTopRawMass", "m_{t,lep}^{raw} [GeV]", "lepTopRawMass_weighted", 35, 50, 400, kPerm2bW);
  makeControlPlot("permutation", "hadTopMass", "m_{t}^{fit} [GeV]", "hadTopMass_weighted", 35, 50, 400, kPerm2bW);
  //*/
}

void makeControlPlot(TString typeForTitle, TString sObservable, TString sObservableShort, TString sFileName,
                     int nbinsx, double xlow, double xup, int plot, bool logY = false, double cut = -999) {
  TCanvas* cControlPlots = new TCanvas("cControlPlots", "cControlPlots", 600, 600);
  cControlPlots->SetLogy(logY);
  cControlPlots->cd();
  
  TString sAddObservable = sObservable;
  
  TString sTTJetsCP = sObservable; sTTJetsCP += " >> hTTJetsCP("; sTTJetsCP += nbinsx;
    sTTJetsCP += ","; sTTJetsCP += xlow; sTTJetsCP += ","; sTTJetsCP += xup; sTTJetsCP += ")";
  TString sTTJetsWP = sObservable; sTTJetsWP += " >> hTTJetsWP("; sTTJetsWP += nbinsx;
    sTTJetsWP += ","; sTTJetsWP += xlow; sTTJetsWP += ","; sTTJetsWP += xup; sTTJetsWP += ")";
  TString sTTJetsUN = sObservable; sTTJetsUN += " >> hTTJetsUN("; sTTJetsUN += nbinsx;
    sTTJetsUN += ","; sTTJetsUN += xlow; sTTJetsUN += ","; sTTJetsUN += xup; sTTJetsUN += ")";
  TString sTTJets = sObservable; sTTJets += " >> hTTJets("; sTTJets += nbinsx;
    sTTJets += ","; sTTJets += xlow; sTTJets += ","; sTTJets += xup; sTTJets += ")";
  
  TString sQCD = sObservable; sQCD += " >> hQCD("; sQCD += nbinsx;
    sQCD += ","; sQCD += xlow; sQCD += ","; sQCD += xup; sQCD += ")";
  TString sWJets = sObservable; sWJets += " >> hWJets("; sWJets += nbinsx;
    sWJets += ","; sWJets += xlow; sWJets += ","; sWJets += xup; sWJets += ")";
  TString sZJets = sObservable; sZJets += " >> hZJets("; sZJets += nbinsx;
    sZJets += ","; sZJets += xlow; sZJets += ","; sZJets += xup; sZJets += ")";
  
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
  
  TString sCut   = "";
  
  switch(plot) {
    case kEvent: {
      sCut += "(MCWeight)*(leptonPt > 30 & combi == 0)";
      break;
    }
    case kEvent2b: {
      sCut += "(MCWeight)*(leptonPt > 30 & combi == 0 & hadBBSSV>1.74 & lepBBSSV>1.74 & hadQBSSV<1.74 & hadQBarBSSV<1.74)";
      break;
    }
    case kPerm2b: {
      sCut += "(MCWeight)*(leptonPt > 30 & hadBBSSV>1.74 & lepBBSSV>1.74 & hadQBSSV<1.74 & hadQBarBSSV<1.74)";
      break;
    }
    case kPerm2bW: {
      sCut += "(MCWeight*hitFitProb)*(leptonPt > 30 & hadBBSSV>1.74 & lepBBSSV>1.74 & hadQBSSV<1.74 & hadQBarBSSV<1.74 & hitFitProb>0.2)";
      break;
    }
  }
  
  TString sCutCP = sCut; sCutCP += "*(target==1)";
  TString sCutWP = sCut; sCutWP += "*(target==0)";
  TString sCutUN = sCut; sCutUN += "*(target==-10)";
  
  // DRAW
  tTTJets     ->Draw(sTTJetsCP,   sCutCP);
  tTTJets     ->Draw(sTTJetsWP,   sCutWP);
  tTTJets     ->Draw(sTTJetsUN,   sCutUN);
  tTTJets     ->Draw(sTTJets,     sCut);
  tQCD        ->Draw(sQCD,        sCut);
  tWJets      ->Draw(sWJets,      sCut);
  tZJets      ->Draw(sZJets,      sCut);
  tSAntiTops  ->Draw(sSAntiTops,  sCut);
  tSAntiTopt  ->Draw(sSAntiTopt,  sCut);
  tSAntiToptW ->Draw(sSAntiToptW, sCut);
  tSTops      ->Draw(sSTops,      sCut);
  tSTopt      ->Draw(sSTopt,      sCut);
  tSToptW     ->Draw(sSToptW,     sCut);
  tData       ->Draw(sData,       sCut);

  
  TH1F* hNull = new TH1F("null", "", nbinsx, xlow, xup);
  
  if (logY) hNull->GetYaxis()->SetRangeUser(0.1, hData->GetMaximum()*10);
  else hNull->GetYaxis()->SetRangeUser(1, hData->GetMaximum()*1.5);
  
  hNull->GetXaxis()->SetTitle(sObservableShort);
  
  /*
	if (eventwise)     hNull->GetYaxis()->SetTitle("Number of events");
	else if (weighted) hNull->GetYaxis()->SetTitle("Sum of weights");
	else               hNull->GetYaxis()->SetTitle("Number of permutations");
	*/
	
	if (plot == kPerm2bW) {
	  TString sTitle = "Sum of "; sTitle += typeForTitle; sTitle += " weights";
	}
	else {
	  TString sTitle = "Number of "; sTitle += typeForTitle; sTitle += "s";
	}
	hNull->GetYaxis()->SetTitle(sTitle);
  
  hTTJetsCP->Scale(lumiWeight[kSigCP]);
  hTTJetsWP->Scale(lumiWeight[kSigWP]);
  hTTJetsUN->Scale(lumiWeight[kSigUN]);
  hTTJets  ->Scale(lumiWeight[kSig]);
  hQCD     ->Scale(lumiWeight[kQCD]);
  hWJets   ->Scale(lumiWeight[kWjets]);
  hZJets   ->Scale(lumiWeight[kZjets]);
  hSAntiTops   ->Scale(lumiWeight[kSAntiTops]);
  hSAntiTopt   ->Scale(lumiWeight[kSAntiTopt]);
  hSAntiToptW  ->Scale(lumiWeight[kSAntiToptW]);
  hSTops   ->Scale(lumiWeight[kSTops]);
  hSTopt   ->Scale(lumiWeight[kSTopt]);
  hSToptW  ->Scale(lumiWeight[kSToptW]);
    
  TH1F* hSTop = new TH1F("hSTop", "", nbinsx, xlow, xup);
  for (int i = 0; i < nbinsx+2; i++) {
    hSTop->SetBinContent(i, (hSAntiTops ->GetBinContent(i)
                          + hSAntiTopt ->GetBinContent(i)
                          + hSAntiToptW->GetBinContent(i)
                          + hSTops ->GetBinContent(i)
                          + hSTopt ->GetBinContent(i)
                          + hSToptW->GetBinContent(i))
                        );
  }
  
  TH1F* hMC = new TH1F("hMC", "", nbinsx, xlow, xup);
  for (int i = 0; i < nbinsx+2; i++) {
    hMC  ->SetBinContent(i, (hTTJets->GetBinContent(i)
                          + hWJets ->GetBinContent(i)
                          + hSTop  ->GetBinContent(i)) * 177./158.
                        );
  }
  hMC->SetFillColor(kBlack);
  hMC->SetFillStyle(3004);
  hMC->SetLineColor(0);
  
  hTTJetsCP->SetFillColor(color_[kSigCP]);
  hTTJetsWP->SetFillColor(color_[kSigWP]);
  hTTJetsUN->SetFillColor(color_[kSigUN]);
  hTTJets  ->SetFillColor(color_[kSigCP]);
  hSTop    ->SetFillColor(color_[kSTop]);
  hQCD     ->SetFillColor(color_[kQCD]);
  hWJets   ->SetFillColor(color_[kWjets]);
  hZJets   ->SetFillColor(color_[kZjets]);
  hData    ->SetMarkerStyle(marker_[kData]);
  
  THStack* stack = new THStack("stack", "");
  stack->Add(hSTop);
  stack->Add(hQCD);
  stack->Add(hWJets);
  stack->Add(hZJets);
  
  if (plot == kPerm2b || plot == kPerm2bW) {
    stack->Add(hTTJetsCP);
    stack->Add(hTTJetsWP);
    stack->Add(hTTJetsUN);
  }
  else stack->Add(hTTJets);
  
  // ---
  //    create legend
  // ---
  TLegend *leg0 = new TLegend(0.25, 0.75, 0.55, 0.925);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  if (plot == kPerm2b || plot == kPerm2bW) {
    leg0->AddEntry( hTTJetsCP, "t#bar{t} correct", "F" );
    leg0->AddEntry( hTTJetsWP, "t#bar{t} wrong", "F" );
    leg0->AddEntry( hTTJetsUN, "t#bar{t} unmatched", "F" );
  }
  else {
    leg0->SetY1(0.8);
    leg0->AddEntry( hTTJets, "t#bar{t}", "F" );
  }
  leg0->AddEntry( hMC,       "t#bar{t} + 1#sigma", "F" );
  
  TLegend *leg1 = new TLegend(0.6, 0.75, 0.9, 0.925);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->AddEntry( hQCD, "QCD", "F" );
  leg1->AddEntry( hWJets, "W#rightarrowl#nu", "F" );
  leg1->AddEntry( hZJets, "Z+jets", "F" );
  leg1->AddEntry( hSTop, "single top", "F" );
  leg1->AddEntry( hData, "Data (4.7 fb ^{-1})", "PL");

  TPaveText *pt = new TPaveText(0.6, 0.7, 0.9, 0.75, "NDC");
  pt->SetBorderSize(0);
  pt->SetFillColor(kBlack);
  pt->SetTextColor(kWhite);
  pt->AddText("Work in progress");
  
  hNull->Draw();
  hMC  ->Draw("SAME");
  stack->Draw("SAME");
  hData->Draw("E,SAME");
  leg0->Draw("");
  leg1->Draw();
  //pt->Draw();
    
  if (cut != -999) drawcutline(cut, hData->GetMaximum());
  
  TString path("controlPlots/"); path+= sFileName; path += ".eps";
  cControlPlots->Print(path);
}
