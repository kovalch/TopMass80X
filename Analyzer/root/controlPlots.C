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

TTree* tTTJets;
TTree* tTTJetsCP;
TTree* tTTJetsWP;
TTree* tTTJetsUN;
TTree* tTTUp;
TTree* tTTDown;
TTree* tQCD;
TTree* tQCDEM1;
TTree* tQCDEM2;
TTree* tQCDEM3;
TTree* tQCDBCE1;
TTree* tQCDBCE2;
TTree* tQCDBCE3;
TTree* tWJets;
TTree* tZJets;
TTree* tSAntiTops;
TTree* tSAntiTopt;
TTree* tSAntiToptW;
TTree* tSTops;
TTree* tSTopt;
TTree* tSToptW;
TTree* tData;

bool upDown = false;
bool qcd = true;
bool muon = true;
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

// return the PU weights for the different samples
double calcPUWeight_(double weight)
{
  // ----
  // default weight to be used for fall11 samples with S6 PU scenario
  // ----
  // New prescription
  double weightsPUFall11_v6[] = {0.171907, 0.0210348, 0.0518887, 0.52666, 1.33025, 1.57157, 1.48441, 1.378, 1.31076, 1.33109, 1.36731, 1.4117, 1.46509, 1.46766, 1.35429, 1.11039, 0.787917, 0.479174, 0.247366, 0.108575, 0.040871, 0.013655, 0.0041918, 0.00124305, 0.000367625, 0.000111436, 3.44048e-05, 1.08332e-05, 3.70591e-06, 1.8736e-06, 2.1086e-06, 3.90444e-06, 8.16201e-06, 1.70748e-05, 3.50119e-05, 7.07827e-05, 0.000138927, 0.000274857, 0.000519194, 0.000960464, 0.00179002, 0.00310324, 0.00552995, 0.00888662, 0.0157832, 0.0236279, 0.0442009, 0.0578902, 0.100113, 0};
  
  // ----
  // default weight to be used for fall11 samples with S6 PU scenario
  // ----
  // Old prescription
  double weightsPUFall11Old_v6[] = {0.00148423, 0.00868081, 0.0324432, 0.379216, 1.12028, 1.51638, 1.57128, 1.52891, 1.41227, 1.36669, 1.38674, 1.44307, 1.5097, 1.5106, 1.36793, 1.06773, 0.692711, 0.366921, 0.156193, 0.0534647, 0.014939, 0.00359882, 0.000806089, 0.00018234, 4.16128e-05, 9.09789e-06, 1.79256e-06, 3.10132e-07, 4.63742e-08, 5.95318e-09, 6.60696e-10, 6.29379e-11, 5.18397e-12, 3.64969e-13, 2.1941e-14, 1.13633e-15, 4.99361e-17, 1.93311e-18, 6.24369e-20, 1.72576e-21, 4.19912e-23, 8.30451e-25, 1.47498e-26, 2.06395e-28, 2.78826e-30, 2.7731e-32, 3.00977e-34, 1.99696e-36, 1.52742e-38, 37324.4};

  // ----
  // fall11 weight with S6 PU scenario shifted up by 8%
  // ----
  // 3.36fb-1 !!! precale weighted !!!
  double weightsPUFall11Plus08_v6[] = {0.212137, 0.481649, 0.777523, 1.07205, 1.32273, 1.48335, 1.61874, 1.62923, 1.61937, 1.5503, 1.46811, 1.36091, 1.20766, 1.05833, 0.899226, 0.759718, 0.619675, 0.500285, 0.396308, 0.310774, 0.242083, 0.183319, 0.138713, 0.102823, 0.0763742, 0.0557689, 0.0399365, 0.0285313, 0.0205054, 0.0143279, 0.00993792, 0.00685962, 0.00475568, 0.00326323, 0.0022337, 0.00148687, 0.00100628, 0.000653758, 0.000451267, 0.000298719, 0.000203806, 0.000130117, 8.29142e-05, 5.36985e-05, 3.74398e-05, 2.55464e-05, 1.86526e-05, 1.10057e-05, 6.52171e-06, 5.03728e-06, 2.82379e-06, 1.83314e-06, 1.28827e-06, 1.29117e-06, 4.67501e-07, 3.67301e-07, 1.44861e-06, 0, 2.71084e-07, 5.78719e-08, 4.8986e-08, 2.05502e-08, 8.54523e-09, 0, 0, 5.82806e-10};

  // ----
  // fall11 weight with S6 PU scenario shifted down by 8%
  // ----
  // 3.36fb-1 !!! precale weighted !!!
  double weightsPUFall11Minus08_v6[] = {0.460591, 0.93254, 1.34303, 1.6592, 1.84558, 1.87888, 1.87462, 1.73659, 1.59777, 1.42193, 1.25474, 1.08425, 0.895585, 0.728338, 0.571947, 0.44451, 0.331913, 0.244145, 0.17543, 0.12428, 0.0871451, 0.0592156, 0.0400974, 0.0265369, 0.0175643, 0.0114108, 0.00726104, 0.00460527, 0.00293658, 0.00181994, 0.00111955, 0.000685499, 0.000421757, 0.000256985, 0.000156325, 9.25565e-05, 5.57688e-05, 3.22874e-05, 1.98781e-05, 1.17452e-05, 7.15671e-06, 4.08204e-06, 2.3241e-06, 1.34456e-06, 8.37012e-07, 5.09544e-07, 3.31601e-07, 1.74184e-07, 9.17676e-08, 6.2928e-08, 3.12724e-08, 1.79704e-08, 1.11627e-08, 9.87504e-09, 3.15173e-09, 2.18e-09, 7.5605e-09, 0, 1.09066e-09, 2.04152e-10, 1.51395e-10, 5.56032e-11, 2.02289e-11, 0, 0, 9.2097e-13};

  // ----
  // default weight to be used for summer11 samples with S4 PU scenario
  // ----
  // 3.36 fb-1 !!! precale weighted !!!
  double weightsPUSummer11_v4[] = {0.0234688, 0.20157, 0.461962, 0.819595, 1.18053, 1.439, 1.60761, 1.62866, 1.57094, 1.46643, 1.38996, 1.28788, 1.18575, 1.11281, 1.02391, 0.954518, 0.888418, 0.835735, 0.767096, 0.719189, 0.664107, 0.625866, 0.567383, 0.534894, 0.497094, 0.458134, 0.397003, 0.372562, 0.337936, 0.346835, 0.263635, 0.227641, 0.224675, 0.173008, 0.111392, 0.143859, 0.070867, 0.19547, 0.0564552, 0.0268947, 0.0634434, 0.0296601};
  
  for (int i = 0; i < 50; ++i) {
    if (abs(weight-weightsPUFall11Old_v6[i]) < 1e-6) return weightsPUFall11_v6[i];
  }
}

void AddWeights(TTree* tempTree) {
  std::cout << "Adding alternative weight to tree ..." << std::endl;

  double AltWeight = -1.;
  double AltPUWeight = -1.;
  double PUWeight;
  double muWeight;
  double bWeight;
  int nPU[3];
  TBranch * br = tempTree->Branch("AltWeight", &AltWeight, "AltWeight/D");
  tempTree->SetBranchAddress("nPU", &nPU);
  tempTree->SetBranchAddress("PUWeight", &PUWeight);
  tempTree->SetBranchAddress("muWeight", &muWeight);
  tempTree->SetBranchAddress("bWeight", &bWeight);
  
  for(int i = 0; i < tempTree->GetEntries(); ++i){
    tempTree->GetEntry(i);
    AltPUWeight = calcPUWeight_(PUWeight);
    AltWeight = AltPUWeight * muWeight * bWeight * 9754./7485.8;
    br->Fill();
  }

  std::cout << "Added alternative weight to tree" << std::endl;
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
  if (muon) {
    TFile* fTTJets = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_1.00_muon/analyzeTop.root");
    if (upDown) {
      TFile* fTTUp   = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_scaleup_muon/analyzeTop.root");
      TFile* fTTDown = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_scaledown_muon/analyzeTop.root");
    }
    if (qcd) TFile* fQCD    = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_QCD_muon/analyzeTop.root");
    TFile* fWJets  = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_WJets_muon/analyzeTop.root");
    TFile* fZJets  = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_ZJets_muon/analyzeTop.root");
    TFile* fSAntiTops  = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_Tbar_s_1.00_muon/analyzeTop.root");
    TFile* fSAntiTopt  = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_Tbar_t_1.00_muon/analyzeTop.root");
    TFile* fSAntiToptW = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_Tbar_tW_1.00_muon/analyzeTop.root");
    TFile* fSTops  = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_T_s_1.00_muon/analyzeTop.root");
    TFile* fSTopt  = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_T_t_1.00_muon/analyzeTop.root");
    TFile* fSToptW = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_T_tW_1.00_muon/analyzeTop.root");
    TFile* fData   = new TFile("/scratch/hh/current/cms/user/mseidel/Run2011_muon/analyzeTop.root");
  }
  else {
    TFile* fTTJets = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_1.00_electron/analyzeTop.root");
    if (upDown) {
      TFile* fTTUp   = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_scaleup_electron/analyzeTop.root");
      TFile* fTTDown = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_scaledown_electron/analyzeTop.root");
    }
    if (qcd) {
      TFile* fQCDBCE1       = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_QCD_Pt-20to30_BCtoE_electron/analyzeTop.root");
      TFile* fQCDBCE2       = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_QCD_Pt-30to80_BCtoE_electron/analyzeTop.root");
      TFile* fQCDBCE3      = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_QCD_Pt-80to170_BCtoE_electron/analyzeTop.root");
      TFile* fQCDEM1  = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_QCD_Pt-20to30_EMEnriched_electron/analyzeTop.root");
      TFile* fQCDEM2  = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_QCD_Pt-30to80_EMEnriched_electron/analyzeTop.root");
      TFile* fQCDEM3 = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_QCD_Pt-80to170_EMEnriched_electron/analyzeTop.root");
    }
    TFile* fWJets  = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_WJets_electron/analyzeTop.root");
    TFile* fZJets  = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_ZJets_electron/analyzeTop.root");
    TFile* fSAntiTops  = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_Tbar_s_1.00_electron/analyzeTop.root");
    TFile* fSAntiTopt  = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_Tbar_t_1.00_electron/analyzeTop.root");
    TFile* fSAntiToptW = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_Tbar_tW_1.00_electron/analyzeTop.root");
    TFile* fSTops  = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_T_s_1.00_electron/analyzeTop.root");
    TFile* fSTopt  = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_T_t_1.00_electron/analyzeTop.root");
    TFile* fSToptW = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_T_tW_1.00_electron/analyzeTop.root");
    TFile* fData   = new TFile("/scratch/hh/current/cms/user/mseidel/Run2011_electron/analyzeTop.root");
  }
  
  
  // ---
  //    Get trees
  // ---
  tTTJets = (TTree*) fTTJets->Get(sTree);
  if (upDown) {
    tTTUp   = (TTree*) fTTUp  ->Get(sTree);
    tTTDown = (TTree*) fTTDown->Get(sTree);
  }
  if (qcd) {
    if (muon) tQCD    = (TTree*) fQCD   ->Get(sTree);
    else {
      tQCDBCE1  = (TTree*) fQCDBCE1->Get(sTree);
      tQCDBCE2  = (TTree*) fQCDBCE2->Get(sTree);
      tQCDBCE3  = (TTree*) fQCDBCE3->Get(sTree);
      tQCDEM1   = (TTree*) fQCDEM1->Get(sTree);
      tQCDEM2   = (TTree*) fQCDEM2->Get(sTree);
      tQCDEM3   = (TTree*) fQCDEM3->Get(sTree);
    }
  }
  tWJets  = (TTree*) fWJets ->Get(sTree);
  tZJets  = (TTree*) fZJets ->Get(sTree);
  tSAntiTops  = (TTree*) fSAntiTops ->Get(sTree);
  tSAntiTopt  = (TTree*) fSAntiTopt ->Get(sTree);
  tSAntiToptW = (TTree*) fSAntiToptW->Get(sTree);
  tSTops  = (TTree*) fSTops ->Get(sTree);
  tSTopt  = (TTree*) fSTopt ->Get(sTree);
  tSToptW = (TTree*) fSToptW->Get(sTree);
  tData   = (TTree*) fData  ->Get(sTree);
  
  //AddWeights(tTTJets);

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
  
  //* Test
  //makeControlPlot("jet", "hadBBSSV", "b-disc (SSVHE)", "", "hadBBSSV_event", 30, 1, 7, kEvent2b, true);
  makeControlPlot("event", "jetMultiplicity", "Number of jets", "", "jetMultiplicity_test", 15, 0, 15, kEvent2b);
  //makeControlPlot("event", "bottomSSVJetMultiplicity", "Number of b-jets", "", "bottomSSVJetMultiplicity", 5, 0, 5, kEvent);
  //makeControlPlot("event", "nVertex", "Number of vertices", "", "nVertex_test", 25, 0, 25, kEvent2b);
  //makeControlPlot("jet", "hadBBCSV", "b-disc (CSV)", "", "hadBBCSV_event", 50, 0.5, 1, kEvent2b, true);
  //makeControlPlot("permutation", "hadTopMass", "m_{t}^{fit}","GeV", "hadTopMass_weighted", 70, 50, 400, kPerm2bW);
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
  if (upDown) {
    tTTUp       ->Draw(sTTUp,       sCut);
    tTTDown     ->Draw(sTTDown,     sCut);
  }
  if (qcd)  {
    if (muon) tQCD        ->Draw(sQCD,        sCut);
    else {
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
  if (upDown) {
    hTTUp    ->Scale(lumiWeight[kUp]);
    hTTDown  ->Scale(lumiWeight[kDown]);
  }
  if (qcd)  {
    if (muon) hQCD->Scale(lumiWeight[kQCD]);
    else {
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
  
  double nQCD; muon ? nQCD = hQCD->GetEntries() : nQCD = hQCDEM2->GetEntries();
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
  
  if (!muon && qcd) {
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
                          + hQCD   ->GetBinContent(i)
                          + hWJets ->GetBinContent(i)
                          + hZJets ->GetBinContent(i)
                          + hSTop  ->GetBinContent(i))
                        );
    hMC  ->SetBinError(i, (hTTJets->GetBinContent(i)
                          + hQCD   ->GetBinContent(i)
                          + hWJets ->GetBinContent(i)
                          + hZJets ->GetBinContent(i)
                          + hSTop  ->GetBinContent(i)) * 11./163.
                        );
  }
  if (!upDown) {
    hMC->SetFillColor(kBlack);
    //hMC->SetLineColor(kWhite);
    hMC->SetFillStyle(3004);
    hMC->SetMarkerStyle(0);
  }
  
  if (upDown) {
    TH1F* hMCUp = new TH1F("hMCUp", "", nbinsx, xlow, xup);
    TH1F* hMCDown = new TH1F("hMCDown", "", nbinsx, xlow, xup);
    for (int i = 0; i < nbinsx+2; i++) {
      hMCUp->SetBinContent(i, (hTTUp->GetBinContent(i)
                            + hQCD   ->GetBinContent(i)
                            + hWJets ->GetBinContent(i)
                            + hZJets ->GetBinContent(i)
                            + hSTop  ->GetBinContent(i))
                          );
      hMCDown->SetBinContent(i, (hTTDown->GetBinContent(i)
                            + hQCD   ->GetBinContent(i)
                            + hWJets ->GetBinContent(i)
                            + hZJets ->GetBinContent(i)
                            + hSTop  ->GetBinContent(i))
                          );
    }
    hMCUp    ->SetLineColor(color_[kUp]);
    hMCUp    ->SetMarkerColor(color_[kUp]);
    hMCUp    ->SetMarkerStyle(marker_[kUp]);
    hMCDown  ->SetLineColor(color_[kDown]);
    hMCDown  ->SetMarkerColor(color_[kDown]);
    hMCDown  ->SetMarkerStyle(marker_[kDown]);
  }
  
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
  if (upDown) {
    leg0->AddEntry( hMC,     "nominal", "F" );
    leg0->AddEntry( hMCUp,   "scale up", "PL" );
    leg0->AddEntry( hMCDown, "scale down", "PL" );
  }
  else {
    if (plot == kPerm2b || plot == kPerm2bW) {
    leg0->AddEntry( hTTJetsCP, "t#bar{t} correct", "F" );
    leg0->AddEntry( hTTJetsWP, "t#bar{t} wrong", "F" );
    leg0->AddEntry( hTTJetsUN, "t#bar{t} unmatched", "F" );
    leg0->AddEntry( hMC,   "t#bar{t} uncertainty", "F" );
    }
    else {
      //leg0->SetY1(0.8);
      leg0->AddEntry( hTTJets, "t#bar{t}", "F" );
      leg0->AddEntry( hMC,   "t#bar{t} uncertainty", "F" );
      leg0->AddEntry( hData, "Data (5.0 fb ^{-1})", "PL");
    }
    //leg0->AddEntry((TObject*)0, sConstFit, "");
  }
  
  TLegend *leg1 = new TLegend(0.6, 0.75, 0.9, 0.925);
  leg1->SetTextSize(0.03);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  if (!upDown) {
    if (qcd) leg1->AddEntry( hQCD, "QCD", "F" );
    leg1->AddEntry( hWJets, "W+jets", "F" );
    leg1->AddEntry( hZJets, "Z+jets", "F" );
    leg1->AddEntry( hSTop, "single top", "F" );
    if (plot == kPerm2b || plot == kPerm2bW) leg1->AddEntry( hData, "Data (5.0 fb ^{-1})", "PL");
  }

  TPaveText *pt = new TPaveText(0.6, 0.7, 0.9, 0.75, "NDC");
  pt->SetBorderSize(0);
  pt->SetFillColor(kBlack);
  pt->SetTextColor(kWhite);
  pt->AddText("Work in progress");
  
  hNull->Draw();
  if (!upDown) stack->Draw("SAME");
  else {
    hMCDown->Draw("E,SAME");
    hMCUp->Draw("E,SAME");
  }
  hData->Draw("E,SAME");
  tdrStyle->SetErrorX(0.5);
  hMC  ->Draw("SAME, E2");
  leg0->Draw("");
  leg1->Draw();
    
  if (cut != -999) drawcutline(cut, hData->GetMaximum());
  
  gPad->RedrawAxis();
  
  if (muon) DrawCMSPrelMuon();
  else DrawCMSPrelElectron();
  
  TString path("controlPlots/"); path+= sFileName; if (!muon) path += "_electron"; path += ".eps";
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
  
}
