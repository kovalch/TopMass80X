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
TTree* tTTUp;
TTree* tTTDown;
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

bool upDown = true;
bool qcd = false;
//TString sTree("analyzeHitFit/eventTree");
TString sTree("analyzeHitFit/eventTree");

std::vector<double> lumiWeight;

enum plotType   {kEvent, kEvent2b, kEvent2bW, kPerm2b, kPerm2bW};

enum samples    {kSigCP, kSigWP, kSigUN,  kSig  , kUp     , kDown  , kZjets  , kWjets  , kQCD   , kSTop   , kDiBos, kData , kWW, kWZ, kZZ, kSAntiTops, kSAntiTopt, kSAntiToptW, kSTops  , kSTopt  , kSToptW };
int color_ [] = {kRed+1, kRed-7, kRed-10, kGray , kGreen+1, kRed+1 , kAzure-2, kGreen-3, kYellow, kMagenta, 10    , kBlack, 10 , 10 , 10 , kMagenta  , kMagenta  , kMagenta   , kMagenta, kMagenta, kMagenta};
int marker_[] = {20,     22,     22,      20,     22      , 23     , 29      , 23      , 21     , 27      , 28    , 20    , 28 , 28 , 28 , 27        , 27        , 27         , 27      , 27      , 27      };


void drawcutline(double cutval, double maximum)
{
  TLine *cut = new TLine();
  cut->SetLineWidth(2);
  cut->SetLineStyle(7);
  cut->SetLineColor(1);
  cut->DrawLine(cutval, 0.,cutval, maximum);
}

void controlPlots_TriggerComparison()
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
  if (upDown) {
    TFile* fTTUp   = new TFile("/scratch/hh/current/cms/user/mseidel/PAS/Run2011_MuHad/analyzeTop.root");
    TFile* fTTDown = new TFile("/scratch/hh/current/cms/user/mseidel/PAS/Run2011_ElectronHad/analyzeTop.root");
  }
  TFile* fData   = new TFile("/scratch/hh/current/cms/user/mseidel/PAS/Run2011/analyzeTop.root");
  
  // ---
  //    Get trees
  // ---
  if (upDown) {
    tTTUp   = (TTree*) fTTUp  ->Get(sTree);
    tTTDown = (TTree*) fTTDown->Get(sTree);
  }
  tData   = (TTree*) fData  ->Get(sTree);

  /* Test
  //makeControlPlot("jet", "hadBBSSV", "b-disc (SSVHE)", "", "hadBBSSV_event", 30, 1, 7, kEvent2b, true);
  makeControlPlot("event", "jetMultiplicity", "Number of jets", "", "jetMultiplicity", 5, 4, 9, kEvent2b);
  //makeControlPlot("event", "bottomSSVJetMultiplicity", "Number of b-jets", "", "bottomSSVJetMultiplicity", 5, 0, 5, kEvent);
  //makeControlPlot("event", "nVertex", "Number of vertices", "", "nVertex_event", 25, 0, 25, kPerm2bW);
  //makeControlPlot("jet", "hadBBCSV", "b-disc (CSV)", "", "hadBBCSV_event", 50, 0.5, 1, kEvent2b, true);
  //makeControlPlot("permutation", "hadTopMass", "m_{t}^{fit}","GeV", "hadTopMass_weighted", 70, 50, 400, kPerm2bW);
  //*/
  
  //* Event
  makeControlPlot("event", "bottomSSVJetMultiplicity", "Number of b-jets (SSV)", "", "bottomSSVJetMultiplicity_event", 5, 0, 5, kEvent2b);
  makeControlPlot("event", "bottomCSVJetMultiplicity", "Number of b-jets (CSV)", "", "bottomCSVJetMultiplicity_event", 5, 0, 5, kEvent2b);
  makeControlPlot("event", "jetMultiplicity", "Number of jets", "", "jetMultiplicity_event", 5, 4, 9, kEvent2b);
  makeControlPlot("event", "nVertex", "Number of vertices", "", "nVertex_event", 25, 0, 25, kEvent2b);
  
  makeControlPlot("lepton", "leptonPt", "p_{T,e/#mu}", "GeV", "leptonPt_event", 40, 0, 200, kEvent2b);
  makeControlPlot("lepton", "leptonEta", "#eta_{e/#mu}", "", "leptonEta_event", 25, -2.5, 2.5, kEvent2b);
  makeControlPlot("event", "nuRawPt", "MET", "GeV", "nuRawPt_event", 40, 0, 200, kEvent2b);
  
  makeControlPlot("jet", "hadBRawPt", "p_{T,b}", "GeV", "BRawPt_event", 40, 0, 200, kEvent2b);
  makeControlPlot("jet", "hadQRawPt", "p_{T,q}", "GeV", "QRawPt_event", 40, 0, 200, kEvent2b);
  makeControlPlot("jet", "hadQBarRawPt", "p_{T,#bar{q}}", "GeV", "QBarRawPt_event", 40, 0, 200, kEvent2b);
  makeControlPlot("jet", "hadBEta", "#eta_{b}", "", "BEta_event", 25, -2.5, 2.5, kEvent2b);
  makeControlPlot("jet", "hadQEta", "#eta_{q}", "", "QEta_event", 25, -2.5, 2.5, kEvent2b);
  makeControlPlot("jet", "hadBBSSV", "b-disc (SSVHE)", "", "hadBBSSV_event", 30, 1, 7, kEvent2b, true);
  makeControlPlot("jet", "hadBBCSV", "b-disc (CSV)", "", "hadBBCSV_event", 50, 0.5, 1, kEvent2b, true);
  //*/
  
  /* Permutation
  makeControlPlot("permutation", "hadWRawMass", "m_{W,had}^{reco}", "GeV", "hadWRawMass", 60, 0, 300, kPerm2b);
  makeControlPlot("permutation", "lepWRawMass", "m_{W,lep}^{reco}", "GeV", "lepWRawMass", 60, 0, 300, kPerm2b);
  makeControlPlot("permutation", "hadTopRawMass", "m_{t,had}^{reco}", "GeV", "hadTopRawMass", 70, 50, 400, kPerm2b);
  makeControlPlot("permutation", "lepTopRawMass", "m_{t,lep}^{reco}", "GeV", "lepTopRawMass", 70, 50, 400, kPerm2b);
  makeControlPlot("permutation", "hadTopMass", "m_{t}^{fit}", "GeV", "hadTopMass", 70, 50, 400, kPerm2b);
  makeControlPlot("permutation", "hitFitProb", "P_{fit}", "", "hitFitProb", 20, 0, 1, kPerm2b, true, 0.2);
  makeControlPlot("permutation", "hitFitChi2", "#chi^{2}_{fit}", "", "hitFitChi2", 20, 0, 10, kPerm2b, false, 3.218875825);
  //*/
  
  /* Permutation, weighted
  makeControlPlot("permutation", "hadWRawMass", "m_{W,had}^{reco}", "GeV", "hadWRawMass_weighted", 60, 0, 300, kPerm2bW);
  makeControlPlot("permutation", "lepWRawMass", "m_{W,lep}^{reco}", "GeV", "lepWRawMass_weighted", 60, 0, 300, kPerm2bW);
  makeControlPlot("permutation", "hadTopRawMass", "m_{t,had}^{reco}","GeV", "hadTopRawMass_weighted", 70, 50, 400, kPerm2bW);
  makeControlPlot("permutation", "lepTopRawMass", "m_{t,lep}^{reco}","GeV", "lepTopRawMass_weighted", 70, 50, 400, kPerm2bW);
  makeControlPlot("permutation", "hadTopMass", "m_{t}^{fit}","GeV", "hadTopMass_weighted", 70, 50, 400, kPerm2bW);
  //*/
  
  /* Systematics up/down
  makeControlPlot("event", "jetMultiplicity", "Number of jets", "", "jetMultiplicity_scale", 5, 4, 9, kEvent2b);
  makeControlPlot("event", "hadWRawMass", "m_{W,had}^{reco}", "GeV", "hadWRawMass_scale", 30, 0, 300, kEvent2b);
  makeControlPlot("event", "hadTopMass", "m_{t}^{fit}", "GeV", "hadTopMass_scale", 35, 50, 400, kEvent2b);
  //*/
}

void makeControlPlot(TString typeForTitle, TString sObservable, TString sObservableShort, TString sUnit, TString sFileName,
                     int nbinsx, double xlow, double xup, int plot, bool logY = false, double cut = -999) {
  TCanvas* cControlPlots = new TCanvas("cControlPlots", "cControlPlots", 600, 600);
  cControlPlots->SetLogy(logY);
  cControlPlots->cd();
  
  TString sAddObservable = sObservable;
  
  TString sTTUp = sObservable; sTTUp += " >> hTTUpL("; sTTUp += nbinsx;
    sTTUp += ","; sTTUp += xlow; sTTUp += ","; sTTUp += xup; sTTUp += ")";
  TString sTTDown = sObservable; sTTDown += " >> hTTDownL("; sTTDown += nbinsx;
    sTTDown += ","; sTTDown += xlow; sTTDown += ","; sTTDown += xup; sTTDown += ")";
    
  TString sData = sObservable; sData += " >> hDataL("; sData += nbinsx;
    sData += ","; sData += xlow; sData += ","; sData += xup; sData += ")";
  
  TString sCut   = "";
  
  switch(plot) {
    case kEvent: {
      sCut += "(MCWeight)*(leptonPt > 30 & combi == 0)";
      break;
    }
    case kEvent2b: {
      sCut += "(MCWeight)*(leptonPt > 30 & combi == 0 & bottomCSVJetMultiplicity == 2)";
      break;
    }
    case kEvent2bW: {
      sCut += "(MCWeight)*(leptonPt > 30 & combi == 0 & hadQBCSV<0.679 & hadQBarBCSV<0.679 & hadBBCSV>0.679 & lepBBCSV>0.679 & hitFitProb>0.2)";
      break;
    }
    case kPerm2b: {
      sCut += "(MCWeight)*(leptonPt > 30 & hadQBCSV<0.679 & hadQBarBCSV<0.679 & hadBBCSV>0.679 & lepBBCSV>0.679)";
      break;
    }
    case kPerm2bW: {
      sCut += "(MCWeight*hitFitProb)*(leptonPt > 30 & hadQBCSV<0.679 & hadQBarBCSV<0.679 & hadBBCSV>0.679 & lepBBCSV>0.679 & hitFitProb>0.2)";
      break;
    }
  }
  
  // DRAW
  if (upDown) {
    tTTUp       ->Draw(sTTUp,       sCut);
    tTTDown     ->Draw(sTTDown,     sCut);
  }
  tData       ->Draw(sData,       sCut);
  
  std::cout << "=============== Event yields ===============" << std::endl;
  std::cout << "IsoMu24: "   << hDataL->Integral()    << std::endl;
  std::cout << "IsoMu17 + TriJet30: "   << hTTUpL->Integral()    << std::endl;
  std::cout << "Ele25 + TriJet30: "   << hTTDownL->Integral()    << std::endl;
  
  TH1F* hTTUp   = new TH1F("hTTUp", "hTTUp", nbinsx, xlow, xup);
  TH1F* hTTDown = new TH1F("hTTDown", "hTTDown", nbinsx, xlow, xup);
  TH1F* hData   = new TH1F("hData", "hData", nbinsx, xlow, xup);
  
  for (int i = 0; i < nbinsx+2; i++) {
      hTTUp   ->SetBinContent(i, (hTTUpL  ->GetBinContent(i)));
      hTTDown ->SetBinContent(i, (hTTDownL->GetBinContent(i)));
      hData   ->SetBinContent(i, (hDataL  ->GetBinContent(i)));
      hTTUp   ->SetBinError(i, (sqrt(hTTUpL->GetBinContent(i))));
      hTTDown ->SetBinError(i, (sqrt(hTTDownL->GetBinContent(i))));
      hData   ->SetBinError(i, (sqrt(hDataL->GetBinContent(i))));
  }
  
  hTTDown->Scale(1./hTTDown->Integral());
  hTTUp  ->Scale(1./hTTUp->Integral());
  hData  ->Scale(1./hData->Integral());
  
  TH1F* hNull = new TH1F("null", "", nbinsx, xlow, xup);
  
  if (logY) hNull->GetYaxis()->SetRangeUser(1, hData->GetMaximum()*10);
  else hNull->GetYaxis()->SetRangeUser(0, hData->GetMaximum()*1.6);
  
  if (sUnit.Length()>0) hNull->GetXaxis()->SetTitle(sObservableShort + " [" + sUnit + "]");
  else hNull->GetXaxis()->SetTitle(sObservableShort);
  
  /*
	if (eventwise)     hNull->GetYaxis()->SetTitle("Number of events");
	else if (weighted) hNull->GetYaxis()->SetTitle("Sum of weights");
	else               hNull->GetYaxis()->SetTitle("Number of permutations");
	*/
	
	TString sBinning(""); sBinning += (xup-xlow)/nbinsx;
	if (sBinning.Length() > 6) sBinning.Resize(4);
	if (plot == kPerm2bW) {
	  TString sTitle = "Fraction of "; sTitle += typeForTitle; sTitle += " weights";
	}
	else {
	  TString sTitle = "Fraction of "; sTitle += typeForTitle; sTitle += "s";
	}
	sTitle += " / "; sTitle += sBinning; sTitle+= " "; sTitle += sUnit;
	hNull->GetYaxis()->SetTitle(sTitle);
  
    
  if (upDown) {
    hTTUp    ->SetLineColor(color_[kUp]);
    hTTUp    ->SetMarkerColor(color_[kUp]);
    hTTUp    ->SetMarkerStyle(marker_[kUp]);
    hTTDown  ->SetLineColor(color_[kDown]);
    hTTDown  ->SetMarkerColor(color_[kDown]);
    hTTDown  ->SetMarkerStyle(marker_[kDown]);
  }

  hData    ->SetMarkerStyle(marker_[kData]);
  
  // ---
  //    create legend
  // ---
  TLegend *leg0 = new TLegend(0.25, 0.75, 0.55, 0.925);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  if (upDown) {
    leg0->AddEntry( hData,    "IsoMu24 (5.0 fb ^{-1})", "PL" );
    leg0->AddEntry( hTTUp,    "IsoMu17 + TriJet30", "PL" );
    leg0->AddEntry( hTTDown,  "Ele25 + TriJet30", "PL" );
  }
  
  hNull->Draw();
  hData->Draw("E,SAME");
  hTTDown->Draw("E,SAME");
  hTTUp->Draw("E,SAME");
  tdrStyle->SetErrorX(0.5);
  leg0->Draw("");
    
  if (cut != -999) drawcutline(cut, hData->GetMaximum());
  
  gPad->RedrawAxis();
  
  DrawCMSPrel();
  
  TString path("controlPlots_TriggerComparison/"); path+= sFileName; path += ".eps";
  cControlPlots->Print(path);
  
}
