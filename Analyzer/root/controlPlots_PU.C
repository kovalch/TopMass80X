#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "THStack.h"
#include "TPaveStats.h"

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
int marker_[] = {24,     22,     22,      20,     29      , 23      , 21     , 27      , 28    , 20    , 28 , 28 , 28 , 27        , 27        , 27         , 27      , 27      , 27      };


void drawcutline(double cutval, double maximum)
{
  TLine *cut = new TLine();
  cut->SetLineWidth(2);
  cut->SetLineStyle(7);
  cut->SetLineColor(1);
  cut->DrawLine(cutval, 0.,cutval, maximum);
}

void controlPlots_PU()
{
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetNdivisions(505, "X");
	gStyle->SetTitleYOffset(1.25);
  gStyle->SetOptStat(0);
  
  //TH1::SetDefaultSumw2(true);

  //TFile* testFile = new TFile("./analyzeTop_1725.root");
  //TTree* testTree = testFile->Get("analyzeHitFit/eventTree"));
  
  // ---
  //    open input files
  // ---
  TFile* fTTJets = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_1.00/analyzeTop.root");
  /*
  TFile* fQCD    = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_QCD/analyzeTop.root");
  TFile* fWJets  = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_WJets/analyzeTop.root");
  TFile* fZJets  = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_ZJets/analyzeTop.root");
  TFile* fSAntiTops  = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_Tbar_s-channel/analyzeTop.root");
  TFile* fSAntiTopt  = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_Tbar_t-channel/analyzeTop.root");
  TFile* fSAntiToptW = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_Tbar_tW-channel/analyzeTop.root");
  TFile* fSTops  = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_T_s-channel/analyzeTop.root");
  TFile* fSTopt  = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_T_t-channel/analyzeTop.root");
  TFile* fSToptW = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_T_tW-channel/analyzeTop.root");
  */
  TFile* fData   = new TFile("/scratch/hh/current/cms/user/mseidel/Run2011/analyzeTop.root");
  
  // ---
  //    Get trees
  // ---
  tTTJets = (TTree*) fTTJets->Get("analyzeHitFit/eventTree");
  /*
  tQCD    = (TTree*) fQCD   ->Get("analyzeHitFit/eventTree");
  tWJets  = (TTree*) fWJets ->Get("analyzeHitFit/eventTree");
  tZJets  = (TTree*) fZJets ->Get("analyzeHitFit/eventTree");
  tSAntiTops  = (TTree*) fSAntiTops ->Get("analyzeHitFit/eventTree");
  tSAntiTopt  = (TTree*) fSAntiTopt ->Get("analyzeHitFit/eventTree");
  tSAntiToptW = (TTree*) fSAntiToptW->Get("analyzeHitFit/eventTree");
  tSTops  = (TTree*) fSTops ->Get("analyzeHitFit/eventTree");
  tSTopt  = (TTree*) fSTopt ->Get("analyzeHitFit/eventTree");
  tSToptW = (TTree*) fSToptW->Get("analyzeHitFit/eventTree");
  */
  tData   = (TTree*) fData  ->Get("analyzeHitFit/eventTree");
  
  makeControlPlot("hadWRawMass0", "m_{W}^{raw} [GeV]", "hadWRawMass0", 70, 95);
  //makeControlPlot("hadTopMass", "m_{t}^{fit} [GeV]", "hadTopMass", 165, 190);
}

void makeControlPlot(TString sObservable, TString sObservableShort, TString sFileName,
                     double ylow, double yup, double uncertaintyScale = 1.) {
  
  TCanvas* cControlPlots = new TCanvas("cControlPlots", "cControlPlots", 600, 600);
  
  TF1* fitTTJets = new TF1("fitTTJets", "[0]+(x-7)*[1]");
  fitTTJets->SetLineWidth(2);
  fitTTJets->SetLineColor(color_[kSigCP]);
  
  TF1* fitData = new TF1("fitData", "[0]+(x-7)*[1]", 3, 15);
  fitData->SetLineWidth(2);
  fitData->SetLineColor(color_[kData]);
  fitData->SetParameter(1, 0);
  fitData->SetParLimits(1, -1, 1);
  
  TF1* fitDataUp = new TF1("fitDataUp", "[0]+(x-7)*[1]+sqrt([2]^2+((x-7)*[3])^2)", 0, 15);
  TF1* fitDataDown = new TF1("fitDataDown", "[0]+(x-7)*[1]-sqrt([2]^2+((x-7)*[3])^2)", 0, 15);
  TF1* fitTTJetsUp = new TF1("fitTTJetsUp", "[0]+(x-7)*[1]+sqrt([2]^2+((x-7)*[3])^2)", 0, 15);
  TF1* fitTTJetsDown = new TF1("fitTTJetsDown", "[0]+(x-7)*[1]-sqrt([2]^2+((x-7)*[3])^2)", 0, 15);
  
  TString selection("(0.5*MCWeight)*(hitFitProb > 0.2 & hadQBSSV<1.74 & hadQBarBSSV<1.74 & hadBBSSV>1.74 & lepBBSSV>1.74 & leptonPt > 30)");
  
  TString sTTJets = sObservable; sTTJets += ":nVertex >> hTTJets(15, 0, 15, 500, 0, 500)";
  TString sData   = sObservable; sData   += ":nVertex >> hData  (15, 0, 15, 500, 0, 500)";
    
  tTTJets->Draw(sTTJets, selection);
  hTTJets->FitSlicesY(0, 0, -1, 30);
  //hTTJets->ProfileX("pTTJetsW");
  hTTJets_1->GetYaxis()->SetRangeUser(ylow, yup);
  hTTJets_1->SetLineColor(color_[kSigCP]);
  hTTJets_1->SetMarkerStyle(marker_[kSigCP]);
  hTTJets_1->SetMarkerColor(color_[kSigCP]);
  hTTJets_1->Fit("fitTTJets", "EM");

  gPad->Update();
  
  TPaveStats* statsTTJets = (TPaveStats*) hTTJets_1->GetListOfFunctions()->FindObject("stats");
  statsTTJets->SetTextColor(color_[kSigCP]);
  statsTTJets->SetX1NDC(0.66);
  statsTTJets->SetY1NDC(0.83);
  statsTTJets->SetX2NDC(0.93);
  statsTTJets->SetY2NDC(0.93);
  
  tData->Draw(sData, selection);
  hData->FitSlicesY(0, 0, -1, 30);
  //hData->ProfileX("pDataW");
  hData_1->GetYaxis()->SetRangeUser(ylow, yup);
  hData_1->SetMarkerStyle(marker_[kData]);
  hData_1->Fit("fitData", "EMR");

  gPad->Update();
  
  TPaveStats* statsData = (TPaveStats*) hData_1->GetListOfFunctions()->FindObject("stats");
  statsData->SetTextColor(color_[kData]);
  statsData->SetX1NDC(0.66);
  statsData->SetY1NDC(0.71);
  statsData->SetX2NDC(0.93);
  statsData->SetY2NDC(0.81);
  
  hData_1->GetXaxis()->SetTitle("Vertex multiplicity");
  hData_1->GetYaxis()->SetTitle(sObservableShort);
  
  /*
  hData_1->SetStats(0);
  hTTJets_1->SetStats(1);
  */
  
  hData_1->Draw();
  hTTJets_1->Draw("same");
  
  fitData->SetRange(0, 15);
  fitData->Draw("same");
  
  fitDataUp->SetParameters(fitData->GetParameter(0), fitData->GetParameter(1), fitData->GetParError(0), fitData->GetParError(1));
  fitDataUp->SetLineColor(color_[kData]);
  fitDataUp->SetLineStyle(2);
  fitDataUp->Draw("same");
  
  fitDataDown->SetParameters(fitData->GetParameter(0), fitData->GetParameter(1), fitData->GetParError(0), fitData->GetParError(1));
  fitDataDown->SetLineColor(color_[kData]);
  fitDataDown->SetLineStyle(2);
  fitDataDown->Draw("same");
  
  fitTTJetsUp->SetParameters(fitTTJets->GetParameter(0), fitTTJets->GetParameter(1), fitTTJets->GetParError(0), fitTTJets->GetParError(1));
  fitTTJetsUp->SetLineColor(color_[kSigCP]);
  fitTTJetsUp->SetLineStyle(2);
  fitTTJetsUp->Draw("same");
  
  fitTTJetsDown->SetParameters(fitTTJets->GetParameter(0), fitTTJets->GetParameter(1), fitTTJets->GetParError(0), fitTTJets->GetParError(1));
  fitTTJetsDown->SetLineColor(color_[kSigCP]);
  fitTTJetsDown->SetLineStyle(2);
  fitTTJetsDown->Draw("same");
  
  
  // ---
  //    create legend
  // ---
  TLegend *leg0 = new TLegend(0.25, 0.75, 0.55, 0.925);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  leg0->AddEntry(hTTJets_1, "t#bar{t} MC", "PL");
  leg0->AddEntry(hData_1, "Data (4.7 fb ^{-1})", "PL");
  leg0->Draw();
  
  TString path("controlPlots/"); path+= sFileName; path += "-nVertex.eps";
  cControlPlots->Print(path);
}

