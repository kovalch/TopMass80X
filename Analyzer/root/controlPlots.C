#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "THStack.h"
#include "THStack.h"

#include "tdrstyle.C"

TChain* tTTJets;
TChain* tTTJetsCP;
TChain* tTTJetsWP;
TChain* tTTJetsUN;
TChain* tQCD;
TChain* tQCDEM1;
TChain* tQCDEM2;
TChain* tQCDEM3;
TChain* tQCDBCE1;
TChain* tQCDBCE2;
TChain* tQCDBCE3;
TChain* tWJets;
TChain* tZJets;
TChain* tSAntiTops;
TChain* tSAntiTopt;
TChain* tSAntiToptW;
TChain* tSTops;
TChain* tSTopt;
TChain* tSToptW;
TChain* tData;

bool qcd = false;

enum lepton           { kElectron, kMuon, kAll};
TString lepton_ [3] = { "electron", "muon", "all"};

int channel = 2;

//TString sTree("analyzeHitFit/eventTree");
TString sTree("analyzeHitFit/eventTree");

std::vector<double> lumiWeight;

enum plotType   {kEvent, kEvent2b, kEvent2bW, kPerm2b, kPerm2bW};
enum enumForPUWeights {kSummer11_v4, kFall10, kFall11Old_v6, kFall11_v6, kFall11Plus08_v6, kFall11Minus08_v6}; 

/*
enum samples    {kSigCP, kSigWP, kSigUN,  kSig  , kUp     , kDown  , kZjets  , kWjets  , kQCD   , kQCDe1 , kQCDe2 , kQCDe3 , kQCDe4 , kQCDe5 , kQCDe6 , kSTop   , kDiBos, kData , kWW, kWZ, kZZ, kSAntiTops, kSAntiTopt, kSAntiToptW, kSTops  , kSTopt  , kSToptW };
int color_ [] = {kRed+1, kRed-7, kRed-10, kGray , kGreen+1, kRed+1 , kAzure-2, kGreen-3, kYellow, kYellow, kYellow, kYellow, kYellow, kYellow, kYellow, kMagenta, 10    , kBlack, 10 , 10 , 10 , kMagenta  , kMagenta  , kMagenta   , kMagenta, kMagenta, kMagenta};
int marker_[] = {20,     22,     22,      20,     22      , 23     , 29      , 23      , 21     , 21     , 21     , 21     , 21     , 21     , 21     , 27      , 28    , 20    , 28 , 28 , 28 , 27        , 27        , 27         , 27      , 27      , 27      };
*/
                   /*0:*/    /*1:*/    /*2:*/    /*3:*/    
  enum samples    {kSigCP, kSigWP, kSigUN,  kSig  , kUp     , kDown  , kZjets  , kWjets  , 
                   /*4:*/    /*5:*/    /*6:*/    /*7:*/  
                   kQCD    , kSTop   , kDiBos  , kData,   
                   /*8*/     /*9*/     /*10*/    /*11*/    /*12*/    /*13*/
                   kQCDEM1 , kQCDEM2 , kQCDEM3 , kQCDBCE1, kQCDBCE2, kQCDBCE3,  
                   /*14*/    /*15*/    /*16*/
                   kWW     , kWZ     , kZZ     , 
                   /*17*/    /*18*/    /*19*/    /*20*/    /*21*/    /*22*/
                   kSTops  , kSATops , kSTopt  , kSATopt , kSToptW , kSAToptW,
                   /*23*/
                   ENDOFSAMPLEENUM};

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
  gStyle->SetHatchesLineWidth(1.5);
  
  //TH1::SetDefaultSumw2(true);

  //TFile* testFile = new TFile("./analyzeTop_1725.root");
  //TTree* testTree = testFile->Get("analyzeHitFit/eventTree"));
  
  // Define chains
  tTTJets     = new TChain("analyzeHitFit/eventTree");
  tTTJetsCP   = new TChain("analyzeHitFit/eventTree");
  tTTJetsWP   = new TChain("analyzeHitFit/eventTree");
  tTTJetsUN   = new TChain("analyzeHitFit/eventTree");
  tQCD        = new TChain("analyzeHitFit/eventTree");
  tQCDEM1     = new TChain("analyzeHitFit/eventTree");
  tQCDEM2     = new TChain("analyzeHitFit/eventTree");
  tQCDEM3     = new TChain("analyzeHitFit/eventTree");
  tQCDBCE1    = new TChain("analyzeHitFit/eventTree");
  tQCDBCE2    = new TChain("analyzeHitFit/eventTree");
  tQCDBCE3    = new TChain("analyzeHitFit/eventTree");
  tWJets      = new TChain("analyzeHitFit/eventTree");
  tZJets      = new TChain("analyzeHitFit/eventTree");
  tSAntiTops  = new TChain("analyzeHitFit/eventTree");
  tSAntiTopt  = new TChain("analyzeHitFit/eventTree");
  tSAntiToptW = new TChain("analyzeHitFit/eventTree");
  tSTops      = new TChain("analyzeHitFit/eventTree");
  tSTopt      = new TChain("analyzeHitFit/eventTree");
  tSToptW     = new TChain("analyzeHitFit/eventTree");
  tData       = new TChain("analyzeHitFit/eventTree");
  
  
  // ---
  //    open input files
  // ---
  if (channel == kMuon || channel == kAll) {
    tTTJets    ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_1.00_muon/analyzeTop.root");
    if (qcd) tQCD->Add("/scratch/hh/current/cms/user/mseidel/Fall11_QCD_muon/analyzeTop.root");
    tWJets     ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_WJets_muon/analyzeTop.root");
    tZJets     ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_ZJets_muon/analyzeTop.root");
    tSAntiTops ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_Tbar_s_1.00_muon/analyzeTop.root");
    tSAntiTopt ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_Tbar_t_1.00_muon/analyzeTop.root");
    tSAntiToptW->Add("/scratch/hh/current/cms/user/mseidel/Fall11_Tbar_tW_1.00_muon/analyzeTop.root");
    tSTops     ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_T_s_1.00_muon/analyzeTop.root");
    tSTopt     ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_T_t_1.00_muon/analyzeTop.root");
    tSToptW    ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_T_tW_1.00_muon/analyzeTop.root");
    tData      ->Add("/scratch/hh/current/cms/user/mseidel/Run2011_muon/analyzeTop.root");
  }
  if (channel == kElectron || channel == kAll) {
    tTTJets    ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_1.00_electron/analyzeTop.root");
    if (qcd) {
      tQCDBCE1->Add("/scratch/hh/current/cms/user/mseidel/Fall11_QCD_Pt-20to30_BCtoE_electron/analyzeTop.root");
      tQCDBCE2->Add("/scratch/hh/current/cms/user/mseidel/Fall11_QCD_Pt-30to80_BCtoE_electron/analyzeTop.root");
      tQCDBCE3->Add("/scratch/hh/current/cms/user/mseidel/Fall11_QCD_Pt-80to170_BCtoE_electron/analyzeTop.root");
      tQCDEM1 ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_QCD_Pt-20to30_EMEnriched_electron/analyzeTop.root");
      tQCDEM2 ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_QCD_Pt-30to80_EMEnriched_electron/analyzeTop.root");
      tQCDEM3 ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_QCD_Pt-80to170_EMEnriched_electron/analyzeTop.root");
    }
    tWJets     ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_WJets_electron/analyzeTop.root");
    tZJets     ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_ZJets_electron/analyzeTop.root");
    tSAntiTops ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_Tbar_s_1.00_electron/analyzeTop.root");
    tSAntiTopt ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_Tbar_t_1.00_electron/analyzeTop.root");
    tSAntiToptW->Add("/scratch/hh/current/cms/user/mseidel/Fall11_Tbar_tW_1.00_electron/analyzeTop.root");
    tSTops     ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_T_s_1.00_electron/analyzeTop.root");
    tSTopt     ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_T_t_1.00_electron/analyzeTop.root");
    tSToptW    ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_T_tW_1.00_electron/analyzeTop.root");
    tData      ->Add("/scratch/hh/current/cms/user/mseidel/Run2011_electron/analyzeTop.root");
  }

  // ---
  // define weights concerning luminosity
  // ---
  //double luminosity   = 1132;
  double luminosity   = 5000;
  //double luminosity   = 47.4;
  double crossSection = 0;
  double sampleSize   = 0;

  // 7 TeV Monte Carlo samples
  for(unsigned int idx=0; idx<27; ++idx) {
    switch(idx) {
      case kSigCP:
      case kSigWP:
      case kSigUN: 
      case kSig: {
        crossSection = 163.; // 
        sampleSize   = 59613991.; // Fall11
        //sampleSize   = 3700000.; // Summer11
        //sampleSize = 8621453.; // P11
        break;
      }
      case kUp: {
        crossSection = 163.; // 
        //sampleSize   = 1062792.; // matching
        sampleSize   = 930483.; // scale 
        break;
      }
      case kDown: {
        crossSection = 163.; // 
        //sampleSize   = 1065323.; // matching
        sampleSize   = 967055.; // scale 
        break;
      }
      case kQCD: {
        crossSection = 296600000.*0.0002855; // generator crossSection * prefilter efficiency
        sampleSize   = 20416038.; // Fall11
        break;
      }
      case kQCDEM1: {
        crossSection = 0.0106*236100000; // generator crossSection * prefilter efficiency
        sampleSize   = 35721833.; // Fall11
        break;
      }
      case kQCDEM2: {
        crossSection = 0.061*59440000; // generator crossSection * prefilter efficiency
        sampleSize   = 70392060.; // Fall11
        break;
      }
      case kQCDEM3: {
        crossSection = 0.159*898200; // generator crossSection * prefilter efficiency
        sampleSize   = 8150672.; // Fall11
        break;
      }
      case kQCDBCE1: {
        crossSection = 0.00059*236100000; // generator crossSection * prefilter efficiency
        sampleSize   = 2071133.; // Fall11
        break;
      }
      case kQCDBCE2: {
        crossSection = 0.00242*59440000; // generator crossSection * prefilter efficiency
        sampleSize   = 2030033.; // Fall11
        break;
      }
      case kQCDBCE3: {
        crossSection = 0.0105*898200; // generator crossSection * prefilter efficiency
        sampleSize   = 1082691.; // Fall11
        break;
      }
      case kWjets: {
        crossSection = 31314.;
        sampleSize   = 81352581.; // Fall11
        break;
      }
      case kZjets: {
        crossSection = 3048.;
        sampleSize   = 35032553.; // Fall11
        break;
      }
      case kSATops: {
        crossSection = 1.44; //
        sampleSize   = 137980.; // Fall11
        break;
      }
      case kSATopt: {
        crossSection = 22.65; //
        sampleSize   = 1944826.; // Fall11
        break;
      }
      case kSAToptW: {
        crossSection = 7.87; //
        sampleSize   = 785246.; // Fall11
        break;
      }
      case kSTops: {
        crossSection = 3.19; //
        sampleSize   = 259971.; // Fall11
        break;
      }
      case kSTopt: {
        crossSection = 41.92; //
        sampleSize   = 3900171.; // Fall11
        break;
      }
      case kSToptW: {
        crossSection = 7.87; //
        sampleSize   = 814390.; // Fall11
        break;
      }
      default: {
        crossSection = 0.; // 
        sampleSize   = 1.; // Dummy
      }
    }
    lumiWeight.push_back(luminosity*crossSection/sampleSize);
  }
  
  /* Test
  //makeControlPlot("jet", "hadBBSSV", "b-disc (SSVHE)", "", "hadBBSSV_event", 30, 1, 7, kEvent2b, true);
  //makeControlPlot("event", "jetMultiplicity", "Number of jets", "", "jetMultiplicity_test", 15, 0, 15, kEvent2b);
  //makeControlPlot("event", "bottomSSVJetMultiplicity", "Number of b-jets", "", "bottomSSVJetMultiplicity", 5, 0, 5, kEvent);
  //makeControlPlot("event", "nVertex", "Number of vertices", "", "nVertex_test", 25, 0, 25, kEvent2b);
  //makeControlPlot("jet", "hadBBCSV", "b-disc (CSV)", "", "hadBBCSV_event", 50, 0.5, 1, kEvent2b, true);
  makeControlPlot("permutation", "hadTopMass", "m_{t}^{fit}","GeV", "hadTopMass_weighted", 70, 50, 400, kPerm2bW);
  //makeControlPlot("permutation", "hadWRawMass", "m_{W}^{reco}", "GeV", "hadWRawMass", 60, 0, 300, kPerm2b);
  //makeControlPlot("permutation", "hadWRawMass", "m_{W}^{reco}", "GeV", "hadWRawMass_weighted", 60, 0, 300, kPerm2bW);
  //*/
  
  /* Event
  makeControlPlot("event", "bottomSSVJetMultiplicity", "Number of b-jets (SSV)", "", "bottomSSVJetMultiplicity_event", 5, 0, 5, kEvent2b);
  makeControlPlot("event", "bottomCSVJetMultiplicity", "Number of b-jets (CSV)", "", "bottomCSVJetMultiplicity_event", 5, 0, 5, kEvent2b);
  makeControlPlot("event", "jetMultiplicity", "Number of jets", "", "jetMultiplicity_event", 5, 4, 9, kEvent2b);
  makeControlPlot("event", "nVertex", "Number of vertices", "", "nVertex_event", 25, 0, 25, kEvent2b);
  
  makeControlPlot("lepton", "leptonPt", "p_{T,l}", "GeV", "leptonPt_event", 40, 0, 200, kEvent2b);
  makeControlPlot("lepton", "leptonEta", "#eta_{l}", "", "leptonEta_event", 25, -2.5, 2.5, kEvent2b);
  makeControlPlot("event", "nuRawPt", "MET", "GeV", "nuRawPt_event", 40, 0, 200, kEvent2b);
  
  makeControlPlot("jet", "hadBRawPt", "p_{T,b}", "GeV", "BRawPt_event", 40, 0, 200, kEvent2b);
  makeControlPlot("jet", "hadQRawPt", "p_{T,q}", "GeV", "QRawPt_event", 40, 0, 200, kEvent2b);
  makeControlPlot("jet", "hadQBarRawPt", "p_{T,#bar{q}}", "GeV", "QBarRawPt_event", 40, 0, 200, kEvent2b);
  makeControlPlot("jet", "hadBEta", "#eta_{b}", "", "BEta_event", 25, -2.5, 2.5, kEvent2b);
  makeControlPlot("jet", "hadQEta", "#eta_{q}", "", "QEta_event", 25, -2.5, 2.5, kEvent2b);
  makeControlPlot("jet", "hadBBSSV", "b-disc (SSVHE)", "", "hadBBSSV_event", 30, 1, 7, kEvent2b, true);
  makeControlPlot("jet", "hadBBCSV", "b-disc (CSV)", "", "hadBBCSV_event", 50, 0.5, 1, kEvent2b, true);
  //*/
  
  //* Permutation
  makeControlPlot("permutation", "hadWRawMass", "m_{W}^{reco}", "GeV", "hadWRawMass", 60, 0, 300, kPerm2b);
  makeControlPlot("permutation", "lepWRawMass", "m_{W,lep}^{reco}", "GeV", "lepWRawMass", 60, 0, 300, kPerm2b);
  makeControlPlot("permutation", "hadTopRawMass", "m_{t,had}^{reco}", "GeV", "hadTopRawMass", 70, 50, 400, kPerm2b);
  makeControlPlot("permutation", "lepTopRawMass", "m_{t,lep}^{reco}", "GeV", "lepTopRawMass", 70, 50, 400, kPerm2b);
  makeControlPlot("permutation", "hadTopMass", "m_{t}^{fit}", "GeV", "hadTopMass", 70, 50, 400, kPerm2b);
  makeControlPlot("permutation", "hitFitProb", "P_{gof}", "", "hitFitProb", 20, 0, 1, kPerm2b, true, 0.2);
  makeControlPlot("permutation", "hitFitChi2", "#chi^{2}_{fit}", "", "hitFitChi2", 20, 0, 10, kPerm2b, false, 3.218875825);
  //*/
  
  //* Permutation, weighted
  makeControlPlot("permutation", "hadWRawMass", "m_{W}^{reco}", "GeV", "hadWRawMass_weighted", 60, 0, 300, kPerm2bW);
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

void makeControlPlot(TString typeForTitle, TString sObservable, TString sObservableShort, TString sUnit, TString sFileName, int nbinsx, double xlow, double xup, int plot, bool logY = false, double cut = -999) {
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
  TString sTTUp = sObservable; sTTUp += " >> hTTUp("; sTTUp += nbinsx;
    sTTUp += ","; sTTUp += xlow; sTTUp += ","; sTTUp += xup; sTTUp += ")";
  TString sTTDown = sObservable; sTTDown += " >> hTTDown("; sTTDown += nbinsx;
    sTTDown += ","; sTTDown += xlow; sTTDown += ","; sTTDown += xup; sTTDown += ")";
  
  TString sQCD = sObservable; sQCD += " >> hQCD("; sQCD += nbinsx;
    sQCD += ","; sQCD += xlow; sQCD += ","; sQCD += xup; sQCD += ")";
  TString sQCDEM1 = sObservable; sQCDEM1 += " >> hQCDEM1("; sQCDEM1 += nbinsx;
    sQCDEM1 += ","; sQCDEM1 += xlow; sQCDEM1 += ","; sQCDEM1 += xup; sQCDEM1 += ")";
  TString sQCDEM2 = sObservable; sQCDEM2 += " >> hQCDEM2("; sQCDEM2 += nbinsx;
    sQCDEM2 += ","; sQCDEM2 += xlow; sQCDEM2 += ","; sQCDEM2 += xup; sQCDEM2 += ")";
  TString sQCDEM3 = sObservable; sQCDEM3 += " >> hQCDEM3("; sQCDEM3 += nbinsx;
    sQCDEM3 += ","; sQCDEM3 += xlow; sQCDEM3 += ","; sQCDEM3 += xup; sQCDEM3 += ")";
  TString sQCDBCE1 = sObservable; sQCDBCE1 += " >> hQCDBCE1("; sQCDBCE1 += nbinsx;
    sQCDBCE1 += ","; sQCDBCE1 += xlow; sQCDBCE1 += ","; sQCDBCE1 += xup; sQCDBCE1 += ")";
  TString sQCDBCE2 = sObservable; sQCDBCE2 += " >> hQCDBCE2("; sQCDBCE2 += nbinsx;
    sQCDBCE2 += ","; sQCDBCE2 += xlow; sQCDBCE2 += ","; sQCDBCE2 += xup; sQCDBCE2 += ")";
  TString sQCDBCE3 = sObservable; sQCDBCE3 += " >> hQCDBCE3("; sQCDBCE3 += nbinsx;
    sQCDBCE3 += ","; sQCDBCE3 += xlow; sQCDBCE3 += ","; sQCDBCE3 += xup; sQCDBCE3 += ")";
  
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
      sCut += "(MCWeight)*(leptonPt > 30 & combi == 0 & hadQBCSV<0.679 & hadQBarBCSV<0.679 & hadBBCSV>0.679 & lepBBCSV>0.679)";
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
  
  TString sCutCP = sCut; sCutCP += "*(target==1)";
  TString sCutWP = sCut; sCutWP += "*(target==0)";
  TString sCutUN = sCut; sCutUN += "*(target==-10)";
  
  // DRAW
  tTTJets     ->Draw(sTTJetsCP,   sCutCP);
  tTTJets     ->Draw(sTTJetsWP,   sCutWP);
  tTTJets     ->Draw(sTTJetsUN,   sCutUN);
  tTTJets     ->Draw(sTTJets,     sCut);
  if (qcd)  {
    if (channel == kMuon) tQCD        ->Draw(sQCD,        sCut);
    else if (channel == kElectron) {
      std::cout << "Drawing electron QCD..." << std::endl;
      /*
      tQCDBCE1  ->Draw(sQCDBCE1, sCut);
      tQCDBCE2  ->Draw(sQCDBCE2, sCut);
      tQCDBCE3  ->Draw(sQCDBCE3, sCut);
      tQCDEM1   ->Draw(sQCDEM1, sCut);
      tQCDEM2   ->Draw(sQCDEM2, sCut);
      tQCDEM3   ->Draw(sQCDEM3, sCut);
      */
      tQCDEM2   ->Draw(sQCDEM2, sCut);
    }
  }
  else TH1F* hQCD = new TH1F("hQCD", "hQCD", nbinsx, xlow, xup);
  tWJets      ->Draw(sWJets,      sCut);
  tZJets      ->Draw(sZJets,      sCut);
  tSAntiTops  ->Draw(sSAntiTops,  sCut);
  tSAntiTopt  ->Draw(sSAntiTopt,  sCut);
  tSAntiToptW ->Draw(sSAntiToptW, sCut);
  tSTops      ->Draw(sSTops,      sCut);
  tSTopt      ->Draw(sSTopt,      sCut);
  tSToptW     ->Draw(sSToptW,     sCut);
  tData       ->Draw(sData,       sCut);
  //hData->Scale(-1./100.);
  
  // Double bin error due to correlations
  if (plot == kPerm2b || plot == kPerm2bW) {
    for (int i = 0; i < hData->GetNbinsX(); i++) {
      int binError = hData->GetBinError(i);
      hData->SetBinError(i, binError * 2);
    }
  }
  
  TH1F* hNull = new TH1F("null", "", nbinsx, xlow, xup);
  
  if (logY) hNull->GetYaxis()->SetRangeUser(1, hData->GetMaximum()*10);
  else hNull->GetYaxis()->SetRangeUser(1, hData->GetMaximum()*1.6);
  
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
	  TString sTitle = "Sum of "; sTitle += typeForTitle; sTitle += " weights";
	}
	else {
	  TString sTitle = "Number of "; sTitle += typeForTitle; sTitle += "s";
	}
	sTitle += " / "; sTitle += sBinning; sTitle+= " "; sTitle += sUnit;
	hNull->GetYaxis()->SetTitle(sTitle);
  
  hTTJetsCP->Scale(lumiWeight[kSigCP]);
  hTTJetsWP->Scale(lumiWeight[kSigWP]);
  hTTJetsUN->Scale(lumiWeight[kSigUN]);
  hTTJets  ->Scale(lumiWeight[kSig]);
  if (qcd)  {
    if (channel == kMuon) hQCD->Scale(lumiWeight[kQCD]);
    else if (channel == kElectron) {
      /*
      hQCDBCE1->Scale(lumiWeight[kQCDBCE1]);
      hQCDBCE2->Scale(lumiWeight[kQCDBCE2]);
      hQCDBCE3->Scale(lumiWeight[kQCDBCE3]);
      hQCDEM1->Scale(lumiWeight[kQCDEM1]);
      hQCDEM2->Scale(lumiWeight[kQCDEM2]);
      hQCDEM3->Scale(lumiWeight[kQCDEM3]);
      */
      hQCDEM2->Scale(lumiWeight[kQCDEM2]);
    }
  }
  hWJets   ->Scale(lumiWeight[kWjets]);
  hZJets   ->Scale(lumiWeight[kZjets]);
  hSAntiTops   ->Scale(lumiWeight[kSATops]);
  hSAntiTopt   ->Scale(lumiWeight[kSATopt]);
  hSAntiToptW  ->Scale(lumiWeight[kSAToptW]);
  hSTops   ->Scale(lumiWeight[kSTops]);
  hSTopt   ->Scale(lumiWeight[kSTopt]);
  hSToptW  ->Scale(lumiWeight[kSToptW]);
  
  double nQCD = 1.;
  if (channel == kMuon) nQCD = hQCD->GetEntries();
  else if (channel == kElectron) hQCDEM2->GetEntries();
  double nSTop = hSAntiTops->GetEntries() + hSAntiTopt->GetEntries() + hSAntiToptW->GetEntries() + hSTops->GetEntries() + hSTopt->GetEntries() + hSToptW->GetEntries();
    
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
  
  if (channel == kElectron && qcd) {
    TH1F* hQCD = new TH1F("hQCD", "", nbinsx, xlow, xup);
    for (int i = 0; i < nbinsx+2; i++) {
      hQCD->SetBinContent(i, (/*hQCDBCE1->GetBinContent(i)
                            + hQCDBCE2->GetBinContent(i)
                            + hQCDBCE3->GetBinContent(i)
                            + hQCDEM1 ->GetBinContent(i)
                            + hQCDEM2 ->GetBinContent(i)
                            + hQCDEM3 ->GetBinContent(i)*/
                            hQCDEM2 ->GetBinContent(i))
                          );
    }
  }
  
  TH1F* hMC = new TH1F("hMC", "", nbinsx, xlow, xup); 
  for (int i = 0; i < nbinsx+2; i++) {
    hMC  ->SetBinContent(i, (hTTJets->GetBinContent(i)
                          + hWJets ->GetBinContent(i)
                          + hZJets ->GetBinContent(i)
                          + hSTop  ->GetBinContent(i))
                        );
    hMC  ->SetBinError(i, (hTTJets->GetBinContent(i)
                          + hWJets ->GetBinContent(i)
                          + hZJets ->GetBinContent(i)
                          + hSTop  ->GetBinContent(i)) * 11./163.
                        );
  }
  hMC->SetFillColor(kBlack);
  //hMC->SetLineColor(kWhite);
  hMC->SetFillStyle(3254);
  hMC->SetMarkerStyle(0);
  
  hTTJetsCP->SetFillColor(color_[kSigCP]);
  hTTJetsWP->SetFillColor(color_[kSigWP]);
  hTTJetsUN->SetFillColor(color_[kSigUN]);
  hTTJets  ->SetFillColor(color_[kSigCP]);
  hSTop    ->SetFillColor(color_[kSTop]);
  if (qcd) hQCD     ->SetFillColor(color_[kQCD]);
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
  
  double chi2 = 0;
  for (int i = 0; i < nbinsx; ++i) {
    if (hData->GetBinError(i)) {
      chi2 += abs(hData->GetBinContent(i) - hMC->GetBinContent(i))/hData->GetBinError(i);
      //std::cout << i << " " << abs(hData->GetBinContent(i) - hMC->GetBinContent(i))/hData->GetBinError(i) << std::endl;
    }
  }
  
  char cChi2[6]; sprintf(cChi2, "%3.1f", chi2);
  TString sConstFit("#chi^{2}/ndf="); sConstFit+=cChi2; sConstFit+="/"; sConstFit+=nbinsx;
  
  // ---
  //    create legend
  // ---
  TLegend *leg0 = new TLegend(0.25, 0.75, 0.55, 0.925);
  leg0->SetTextSize(0.03);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  if (plot == kPerm2b || plot == kPerm2bW) {
    leg0->AddEntry( hTTJetsUN, "t#bar{t} unmatched", "F" );
    leg0->AddEntry( hTTJetsWP, "t#bar{t} wrong", "F" );
    leg0->AddEntry( hTTJetsCP, "t#bar{t} correct", "F" );
    leg0->AddEntry( hMC,   "t#bar{t} uncertainty", "F" );
  }
  else {
    //leg0->SetY1(0.8);
    leg0->AddEntry( hTTJets, "t#bar{t}", "F" );
    leg0->AddEntry( hMC,   "t#bar{t} uncertainty", "F" );
    leg0->AddEntry( hData, "Data (5.0 fb ^{-1})", "PL");
  }
  
  TLegend *leg1 = new TLegend(0.6, 0.75, 0.9, 0.925);
  leg1->SetTextSize(0.03);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  if (qcd) leg1->AddEntry( hQCD, "QCD", "F" );
  leg1->AddEntry( hWJets, "W+jets", "F" );
  leg1->AddEntry( hZJets, "Z+jets", "F" );
  leg1->AddEntry( hSTop, "single top", "F" );
  if (plot == kPerm2b || plot == kPerm2bW) leg1->AddEntry( hData, "Data (5.0 fb ^{-1})", "PL");

  TPaveText *pt = new TPaveText(0.6, 0.7, 0.9, 0.75, "NDC");
  pt->SetBorderSize(0);
  pt->SetFillColor(kBlack);
  pt->SetTextColor(kWhite);
  pt->AddText("Work in progress");
  
  hNull->Draw();
  stack->Draw("SAME");

  hData->Draw("E,SAME");
  tdrStyle->SetErrorX(0.5);
  hMC  ->Draw("SAME, E2");
  leg0->Draw("");
  leg1->Draw();
    
  if (cut != -999) drawcutline(cut, hData->GetMaximum());
  
  gPad->RedrawAxis();
  
  DrawCMS(channel);
  
  TString path("controlPlots/"); path+= sFileName; path += "_"; path += lepton_[channel]; path += ".eps";
  cControlPlots->Print(path);
  
  // Uncertainties
  double uTTJets  = sqrt(1./hTTJets->GetEntries() + pow(11./163., 2));
  double uQCD     = sqrt(1./nQCD + pow(0./1., 2));
  double uZJets   = sqrt(1./hZJets->GetEntries() + pow(132./3048., 2));
  double uWJets   = sqrt(1./hWJets->GetEntries() + pow(1558./31314., 2));
  double uSTop    = sqrt(1./nSTop + pow(3.50/79.61, 2));
  
  double iMC = hTTJets->Integral() + hQCD->Integral() + hZJets->Integral() + hWJets->Integral() + hSTop->Integral();
  double eMC = sqrt(pow(hTTJets->Integral()*uTTJets, 2) + pow(hQCD->Integral()*uQCD, 2) + pow(hZJets->Integral()*uZJets, 2) + pow(hWJets->Integral()*uWJets, 2) + pow(hSTop->Integral()*uSTop, 2));
  
  std::cout << "=============== Event yields ===============" << std::endl;
  std::cout << "tt: "     << hTTJets->Integral()  << " +/- " << hTTJets->Integral()*uTTJets << std::endl;
  std::cout << "QCD: "    << hQCD->Integral()     << " +/- " << hQCD->Integral()*uQCD << std::endl;
  std::cout << "Z+jets: " << hZJets->Integral()   << " +/- " << hZJets->Integral()*uZJets << std::endl;
  std::cout << "W+jets: " << hWJets->Integral()   << " +/- " << hWJets->Integral()*uWJets << std::endl;
  std::cout << "st: "     << hSTop->Integral()    << " +/- " << hSTop->Integral()*uSTop << std::endl;
  std::cout << "mc: "     << iMC                  << " +/- " << eMC << std::endl;
  std::cout << "data: "   << hData->Integral()    << std::endl;
  
  std::cout << "Global chi2: " << chi2 << std::endl;
  std::cout << "KS prob: " << hData->KolmogorovTest(hMC) << std::endl;
  std::cout << "KS prob (N): " << hData->KolmogorovTest(hMC, "N") << std::endl;
  
  for (int i = 0; i < nbinsx+2; i++) {
    hMC  ->SetBinError(i, 1e-3);
  }
  
  std::cout << "KS prob (2): " << hData->KolmogorovTest(hMC) << std::endl;
  std::cout << "KS prob (N2): " << hData->KolmogorovTest(hMC, "N") << std::endl;
  
}
