#include <vector>
#include <iostream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TChain.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TMath.h"

#include "tdrstyle_new.C"
#include "CMS_lumi.C"

enum lepton          { kElectron, kMuon, kAll};
TString lepton_ [] = { "electron", "muon", "all"};
enum cat             { kCR, kRad, kBoth};

const int nKin = 17;

TString sBinning[nKin] = {"TopPt", "TopEta", "B1Pt", "B1Eta", "TTBarMass", "TTBarPt", "DeltaRqq", "DeltaRbb", "HT", "NJet", "Q1Pt", "Q1Eta", "WPt", "WEta", "NVertex", "LeptonCharge", "Alpha"};
bool paper[nKin] = {true, false, true, true, true, true, true, true, false, true, false, false, false, false, false, false, false};
TString sData[nKin] = {"top_fitTop1_Pt__", "top_fitTop1_Eta__", "top_fitB1_Pt__", "top_fitB1_Eta__", "top_fitTTBar_M__", "top_fitTTBar_Pt__", "sqrt_pow_top_fitW1Prod1_Eta__-top_fitW1Prod2_Eta__2____pow_TVector2Phi_mpi_pi_top_fitW1Prod1_Phi__-top_fitW1Prod2_Phi___2__", "sqrt_pow_top_fitB1_Eta__-top_fitB2_Eta__2____pow_TVector2Phi_mpi_pi_top_fitB1_Phi__-top_fitB2_Phi___2__", "top_fitB1_Pt___top_fitB2_Pt___top_fitW1Prod1_Pt___top_fitW1Prod2_Pt__", "Max__Iteration___jet_jet_Pt___gt_30_____1", "top_fitW1Prod1_Pt__", "top_fitW1Prod1_Eta__", "top_fitW1_Pt__", "top_fitW1_Eta__", "weight_nVertex", "-top_leptonFlavour_abs_top_leptonFlavour_", "jet_jet_4__Pt____top_recoB1_0__Pt___top_recoB2_0__Pt____2"};
TString sBinNice[nKin] = {"p_{T}^{t,had} [GeV]", "|#eta^{t,had}|", "p_{T}^{b,had} [GeV]", "|#eta^{b,had}|", "m_{t#bar{t}} [GeV]", "p_{T}^{t#bar{t}} [GeV]", "#DeltaR_{q#bar{q}}", "#DeltaR_{b#bar{b}}", "H_{T}^{4} [GeV]", "Number of jets", "p_{T}^{q} [GeV]", "|#eta^{q}|", "p_{T}^{W} [GeV]", "|#eta^{W}|", "N(PV)", "lepton charge", "#alpha"};
TString sBinLatex[nKin] = {"$p_{\\rm T,t,had}$", "$\\left|\\eta_{\\rm t,had}\\right|$", "$p_{\\rm T,b,had}$", "$\\left|\\eta_{\\rm b,had}\\right|$", "$m_{\\ttbar}$", "$p_{\\rm T,\\ttbar}$", "$\\Delta R_{\\qqbar}$", "$\\Delta R_{\\bbbar}$", "$\\HT^{4}$", "Number of jets", "$p_{\\rm T,q,had}^{1}$", "$\\left|\\eta_{\\rm q,had}^{1}\\right|$", "$p_{\\rm T,W,had}$", "$\\left|\\eta_{\\rm W,had}\\right|$", "N(PV)", "lepton charge", "$\\alpha$"};

const int nObs = 6;

TString sObservable[nObs] = {"Entries", "mass_mTop_JES", "mass_mTop", "JES_mTop_JES", "mass_mTop_JES_jsfc", "JES_mTop_JES_jsfc"};
TString sObsLowCase[nObs] = {"Entries", "mass_mTop_JES", "mass_mTop", "JES_mTop_JES", "mass_mTop_JES_jsfc", "JES_mTop_JES_jsfc"};
TString sObsFile[nObs]    = {"Entries", "MT2D", "MT1D", "JES", "MTH", "JESH"};
TString sObsCanvas[nObs] = {"canvas_1", "canvas_2", "canvas_4", "canvas_3", "canvas_6", "canvas_5"};
TString sObsNice[nObs] = {"Number of permutations / bin width", "m_{t}^{2D} #font[122]{\55} <m_{t}^{2D}> [GeV]", "m_{t}^{1D} #font[122]{\55} <m_{t}^{1D}> [GeV]", "JSF #font[122]{\55} <JSF>", "m_{t}^{hyb} #font[122]{\55} <m_{t}^{hyb}> [GeV]", "JSF^{hyb} #font[122]{\55} <JSF^{hyb}>"};
TString sObsNiceCal[nObs] = {"Number of permutations / bin width", "m_{t,cal}^{2D} #font[122]{\55} <m_{t}^{2D}> [GeV]", "m_{t,cal}^{1D} #font[122]{\55} <m_{t}^{1D}> [GeV]", "JSF_{cal} #font[122]{\55} <JSF>", "m_{t,cal}^{hyb} #font[122]{\55} <m_{t}^{hyb}> [GeV]", "JSF_{cal}^{hyb} #font[122]{\55} <JSF^{hyb}>"};
TString sObsNice2[nObs] = {"Number of permutations / bin width", "(m_{t}^{2D} #font[122]{\55} <m_{t}^{2D}>)_{MC} - (m_{t}^{2D} #font[122]{\55} <m_{t}^{2D}>)_{data} [GeV]", "(m_{t}^{1D} #font[122]{\55} <m_{t}^{1D}>)_{MC} - (m_{t}^{1D} #font[122]{\55} <m_{t}^{1D}>)_{data} [GeV]", "(JSF #font[122]{\55} <JSF>)_{MC} - (JSF #font[122]{\55} <JSF>)_{data}", "fixme1", "fixme2"};

std::vector<TString> sampleNames;

double chi2matrix[nKin][nObs][30];

double crossSection   = 230;
double peLumi         = 20000.;

int channel = 2;
bool doCalibration = true;

struct resultContainer {
  TPad* pad;
  TH1F* profile;
  TH1F* profileDiff;
  TH1F* profile2;
  TH1F* error;
  double incl;
};

struct sample {
  bool data;
  bool systematic;
  bool cr;
  bool rad;
  bool drawerror;
  bool calibration;
  TChain* chain;
  TFile* file;
  TCanvas* canvas;
  std::vector<resultContainer> results;
  TF1* fit;
  const char* id;
  const char* label;
  int color;
  int line;
  double size; // EFFECTIVE sample size
  int marker;
  double genMass, genJES;
  sample(bool d, bool sys, bool c_, bool r, bool de, bool cal, const char* i, const char* l, int c = kBlack, int li = 1, double s = 0, int ma = 1)
  : data(d), systematic(sys), cr(c_), rad(r), drawerror(de), calibration(cal), id(i), label(l), color(c), line(li), size(s), marker(ma) {}
};

std::string itos(int number)
{
   std::stringstream ss;
   ss << number;
   return ss.str();
}

TPaveLabel* DrawLabel(TString text, const double x1, const double y1, const double x2, Color_t color = kBlack)
{
  // function to directly draw a label into the active canvas
  double y2 = y1 + 0.05;
  double yOffset = 0.02;
  TPaveLabel *label = new TPaveLabel(x1, y1+yOffset, x2, y2+yOffset, text, "br NDC");
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  label->SetTextSize(0.75);
  label->SetTextAlign(12);
  label->SetTextColor(color);
  label->Draw("same");
  
  return label;
}

void differentialMass(int iBinning = 4)
{
  bool fixme = false;
  
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetPadTopMargin(0.08);
  //gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);
  gStyle->SetEndErrorSize(4);
  
  std::vector<sample> samples;
  std::vector<sample>::iterator it;
  std::vector<sample> calibsamples;
  std::vector<sample>::iterator cit;
  //std::vector<std::vector<double> > values;
  
  samples.push_back(sample(true,  false,  true,  true, false, false, "Run2012", "Data", kBlack, 0, 0, 20));
  
  samples.push_back(sample(false, false,  true,  true,  true, false, "Summer12_TTJetsMS1725_1.00", "MG+MS, Pythia Z2*", kRed+1, 0, 62131965./1.2, 24));
  samples.push_back(sample(false,  true,  true,  true,  true, false, "Summer12_TTJets1725_MGDecays_Z2", "MG, Pythia Z2*", kRed+1, 0, 56000000./1.2));
  
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_flavor:up_FlavorPureBottom", "MG, b-JES up", kBlue+1, 2, 62131965./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_flavor:down_FlavorPureBottom", "MG, b-JES down", kGreen+1, 2, 62131965./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_flavor:up_FlavorPureGluon", "MG, Gluon-JES up", kBlue+1, 2, 62131965./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_flavor:down_FlavorPureGluon", "MG, Gluon-JES down", kGreen+1, 2, 62131965./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_flavor:up_FlavorPureCharm", "MG, Charm-JES up", kBlue+1, 2, 62131965./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_flavor:down_FlavorPureCharm", "MG, Charm-JES down", kGreen+1, 2, 62131965./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_flavor:up_FlavorPureQuark", "MG, Quark-JES up", kBlue+1, 2, 62131965./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_flavor:down_FlavorPureQuark", "MG, Quark-JES down", kGreen+1, 2, 62131965./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_source:up_Total", "MG, JES up", kBlue+1, 2, 62131965./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_source:down_Total", "MG, JES down", kGreen+1, 2, 62131965./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_jer:1", "MG, JER up", kBlue+1, 2, 62131965./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_jer:-1", "MG, JER down", kGreen+1, 2, 62131965./1.2));
  //samples.push_back(sample(false, true,  false, false, false, "Summer12_TTJets1725_source:up_PileUpPtBB", "MG, PUJES up", kBlue+1, 2, 62131965./1.2));
  //samples.push_back(sample(false, true,  false, false, false, "Summer12_TTJets1725_source:down_PileUpPtBB", "MG, PUJES down", kGreen+1, 2, 62131965./1.2));
  
  //*
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_scaleup", "MG, scale up", kRed-1, 5, 3696269./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_scaledown", "MG, scale down", kRed-7, 6, 4004587./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_matchingup", "MG, matching up", kRed-1, 7, 4029823./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_matchingdown", "MG, matching down", kRed-7, 7, 1545688./1.2));
  //*/
  
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_1.00_puUp", "MG, PU up", kBlue+1, 2, 62131965./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_1.00_puDown", "MG, PU down", kGreen+1, 2, 62131965./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_1.00_fNuUp", "MG, fNu up", kBlue+1, 2, 62131965./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_1.00_fNuDown", "MG, fNu down", kGreen+1, 2, 62131965./1.2));
  samples.push_back(sample(false,  true, false, false, false, false, "Summer12_TTJetsMS1725_1.00_frag", "MG, rbLEP", kGreen+1, 2, 62131965./1.2));
  
  //*
  samples.push_back(sample(false, false,  true,  true,  true, false, "Summer12_TTJets1725_MGDecays_P11", "MG, Pythia P11", kMagenta+1, 5, 27000000./1.2, 25));
  samples.push_back(sample(false, false,  true,  true,  true, false, "Summer12_TTJets1725_MGDecays_P11noCR", "MG, Pythia P11noCR", kCyan+1, 6, 27000000./1.2, 26));
  samples.push_back(sample(false, false,  true,  true,  true, false, "Summer12_TTJets1725_powheg", "Powheg, Pythia Z2*", kGreen+1, 9, 21675970./1.2, 27));
  samples.push_back(sample(false, false,  true,  true,  true, false, "Summer12_TTJets1725_powheg_herwig", "Powheg, Herwig 6", kOrange+2, 7, 27684235./1.2, 28));
  samples.push_back(sample(false, false,  true,  true,  true, false, "Summer12_TTJets1725_mcatnlo_herwig", "MC@NLO, Herwig 6", kBlue+1, 2, 32852589.*0.7718/1.2, 29));
  samples.push_back(sample(false, false,  true,  true,  true, false, "Summer12_TTJets1725_sherpa", "Sherpa", kYellow+1, 3, 44000000./1.2, 30));
  
  //*/
  
  /*
  // Debug samples
  samples.push_back(sample(false, false,  true,  true,  true, false, "Summer12_TTJetsMS1715_1.00", "MG+MS, mt = 171.5", kRed, 2, 24439341./1.2));
  samples.push_back(sample(false, false,  true,  true,  true, false, "Summer12_TTJetsMS1735_1.00", "MG+MS, mt = 173.5", kGreen, 2, 26489020./1.2));
  samples.push_back(sample(false, false,  true,  true,  true, false, "Summer12_TTJetsMS1725_1.02", "JES up", kBlack, 1, 62131965./1.2));
  samples.push_back(sample(false, false,  true,  true,  true, false, "Summer12_TTJetsMS1755_1.00", "Mass up", kGray, 1, 40244328./1.2));
  //*/
  
  /*
  // FastSim samples
  samples.push_back(sample(false, false,  true,  true,  true, false, "Summer12_TTJets1725_FSIM_parj81_nominal", "FSIM MG+MS", kRed, 2, 24439341./1.2));
  samples.push_back(sample(false, false,  true,  true,  true, false, "Summer12_TT1725_FSIM_pythia8", "FSIM Pythia8", kGreen, 2, 26489020./1.2));
  samples.push_back(sample(false, false,  true,  true,  true, false, "Summer12_TT1725_amcatnlo_herwigpp:pythia8", "FSIM aMC+H++:P8", kBlack, 1, 62131965./1.2));
  samples.push_back(sample(false, false,  true,  true,  true, false, "Summer12_TT1725_amcatnlo_herwigpp", "FSIM aMC+H++", kGray, 1, 40244328./1.2));
  samples.push_back(sample(false, false,  true,  true,  true, false, "Summer12_TT1725_amcatnlo_pythia8", "FSIM aMC+P8", kGray, 1, 40244328./1.2));
  //*/
  
  // calibration samples
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1665_0.98", "a", kBlack, 1, 27078777./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1665_1.00", "b", kBlack, 1, 27078777./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1665_1.02", "c", kBlack, 1, 27078777./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1695_0.98", "d", kBlack, 1, 39518234./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1695_1.00", "e", kBlack, 1, 39518234./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1695_1.02", "f", kBlack, 1, 39518234./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1715_0.98", "g", kBlack, 1, 24439341./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1715_1.00", "h", kBlack, 1, 24439341./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1715_1.02", "i", kBlack, 1, 24439341./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1725_0.98", "j", kBlack, 1, 62131965./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1725_1.00", "k", kBlack, 1, 62131965./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1725_1.02", "l", kBlack, 1, 62131965./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1735_0.98", "m", kBlack, 1, 26489020./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1735_1.00", "n", kBlack, 1, 26489020./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1735_1.02", "o", kBlack, 1, 26489020./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1755_0.98", "p", kBlack, 1, 40244328./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1755_1.00", "q", kBlack, 1, 40244328./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1755_1.02", "r", kBlack, 1, 40244328./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1785_0.98", "s", kBlack, 1, 24359161./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1785_1.00", "t", kBlack, 1, 24359161./1.2));
  samples.push_back(sample(false, false, false, false, false,  true, "Summer12_TTJetsMS1785_1.02", "u", kBlack, 1, 24359161./1.2));
  
  double minVal[nObs] = {0., 0., 0., 0., 0., 0.};
  double maxVal[nObs] = {0., 0., 0., 0., 0., 0.};
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    sampleNames.push_back(it->label);
    std::cout << "Reading " << it->id << std::endl;
    for (int iObs = 0; iObs < nObs; ++iObs) it->results.push_back(resultContainer());
    
    if (it->data) {
      TFile* inclFile     = new TFile("/afs/desy.de/user/m/mseidel/xxl/CMSSW_5_3_11/src/TopMass/Analyzer/plot/IdeogramMin_Run2012_JEC_Winter14_V8_incl.root");
      TCanvas* inclCanvas = (TCanvas*) inclFile  ->Get("canvas");
      
      for (int iObs = 0; iObs < nObs; ++iObs) {
        TPad* inclPad       = (TPad*)    inclCanvas->GetPrimitive(sObsCanvas[iObs]);
        TH1F* inclProfile   = (TH1F*)    inclPad   ->GetPrimitive(sObservable[iObs]);
        it->results[iObs].incl = inclProfile->GetBinContent(1);
      }
      inclFile->Close();
      
      it->file    = new TFile("/afs/desy.de/user/m/mseidel/xxl/CMSSW_5_3_11/src/TopMass/Analyzer/plot/IdeogramMin_Run2012_JEC_Winter14_V8_" + sData[iBinning] + ".root");
      it->canvas  = (TCanvas*) it->file  ->Get("canvas");
      
      for (int iObs = 0; iObs < nObs; ++iObs) {
        it->results[iObs].pad     = (TPad*)    it->canvas->GetPrimitive(sObsCanvas[iObs]);
        it->results[iObs].profile = (TH1F*)    it->results[iObs].pad   ->GetPrimitive(sObservable[iObs]);
      }
    }
    else {
      if (!doCalibration && it->calibration) continue;
      TString sfBase("/nfs/dust/cms/user/mseidel/pseudoexperiments/topmass_paper_diff_jsfconstraint/");
      
      TChain* inclChain = new TChain("tree");
      inclChain->Add(sfBase + "lepton/" + it->id + "/TopMass/job_*_ensemble.root");
      
      it->chain   = new TChain("tree");
      it->chain->Add(sfBase + "lepton/" + it->id + "/" + sBinning[iBinning] + "/job_*_ensemble.root");
      
      for (int iObs = 0; iObs < nObs; ++iObs) {
        int minFit = 0;
        int maxFit = 100000;
        if (iObs == 1 || iObs == 2 || iObs == 4) {
          minFit = 100;
          maxFit = 250;
        }
        if (iObs == 3 || iObs == 5) maxFit = 2;
        
        // inclusive values -> doubles
        
        TF1* inclGaus = new TF1("inclGaus", "gaus");
        inclChain->Fit("inclGaus", sObsLowCase[iObs], sObsLowCase[iObs] + " < " + itos(maxFit) + " & " + sObsLowCase[iObs] + " > " + itos(minFit), "LEMQ");
        it->results[iObs].incl = inclGaus->GetParameter(1);
        std::cout << "\t" << sObsFile[iObs] << " " << it->results[iObs].incl << std::endl;
        
        // binned values -> hists
        
        it->results[iObs].profile = (TH1F*) samples[0].results[iObs].profile->Clone(it->id); // TODO: unique hist name
        
        for (int i = 1; i < it->results[iObs].profile->GetNbinsX()+1; i++) {
          TF1* gaus = new TF1("gaus", "gaus");
          int entries = it->chain->GetEntries(sObsLowCase[iObs] + " < " + itos(maxFit) + " & " + sObsLowCase[iObs] + " > " + itos(minFit) + " & bin == " + itos(i));
          double pullWidth = 1.;
          if (entries < 100) {
              std::cout << "WARNING: Too less entries: " << itos(entries) << " in " << it->id << " obs " << itos(iObs) << ", bin " << itos(i) << std::endl;
              fixme = true;
          }
          if (entries > 10) {
            // pull
            it->chain->Fit("gaus", sObsLowCase[iObs] + "_Pull", sObsLowCase[iObs] + " < " + itos(maxFit) + " & " + sObsLowCase[iObs] + " > " + itos(minFit) + " & bin == " + itos(i), "LEMQ");
            pullWidth = gaus->GetParameter(2);
            
            // result
            it->chain->Fit("gaus", sObsLowCase[iObs], sObsLowCase[iObs] + " < " + itos(maxFit) + " & " + sObsLowCase[iObs] + " > " + itos(minFit) + " & bin == " + itos(i), "LEMQ");
            it->results[iObs].profile->SetBinContent(i, gaus->GetParameter(1));
            it->results[iObs].profile->SetBinError(i, gaus->GetParameter(2) / sqrt(it->size/(crossSection*peLumi)) * pullWidth);
          }
          else {
            it->results[iObs].profile->SetBinContent(i, -1);
            it->results[iObs].profile->SetBinError(i, 1e6);
          }
          
          // PULL correction for data
          if (strcmp(it->id, "Summer12_TTJetsMS1725_1.00") == 0 && !it->calibration && iObs != 0) {
            samples[0].results[iObs].profile->SetBinError(i, samples[0].results[iObs].profile->GetBinError(i) * pullWidth);
            std::cout << "Pull correction:" << pullWidth << std::endl;
          }
        }
        
        double genMass, genJES;
        it->chain->SetBranchAddress("genMass", &genMass);
        it->chain->SetBranchAddress("genJES",  &genJES);
        it->chain->GetEntry(0);
        it->genMass = genMass;
        it->genJES  = genJES;
      }
      
      inclChain->Delete();
    }
  }
  
  const int iterations = 5;
  double calibration[iterations][nObs][10][4]; // iteration/obs/bin/param
  
  // TODO: on-the-fly calibration
  //*
  if (doCalibration) {
    for (int iter = 0; iter < iterations; ++iter) {
      for (int iObs = 1; iObs < nObs; ++iObs) {
        for (int i = 1; i < samples[0].results[iObs].profile->GetNbinsX()+1; i++) {
          // mt2d, jsf
          if (iObs == 1 || iObs == 3) {
            TH2D* h2Bias = new TH2D("h2Bias", "h2Bias", 1000, 150, 200, 2000, 0.8, 1.2);
            
            for (it = samples.begin(); it != samples.end(); ++it) {
              if (!it->calibration) continue;
              
              double genValue = it->genMass;
              if (iObs == 3) genValue = it->genJES;
              //std::cout << it->id << " " << iObs << ", bin " << i << ", fitted " << it->results[iObs].profile->GetBinContent(i) << ", generated " << genValue << std::endl;
              double m = it->results[1].profile->GetBinContent(i);
              double j = it->results[3].profile->GetBinContent(i);
              if (iter > 0) {
                for (int iterprev = 0; iterprev < iter; ++iterprev) {
                  double mnext = m - calibration[iterprev][1][i][0] - calibration[iterprev][1][i][1]*(m-172.5) - calibration[iterprev][1][i][2]*(j-1.) - calibration[iterprev][1][i][3]*(m-172.5)*(j-1.);
                  double jnext = j - calibration[iterprev][3][i][0] - calibration[iterprev][3][i][1]*(m-172.5) - calibration[iterprev][3][i][2]*(j-1.) - calibration[iterprev][3][i][3]*(m-172.5)*(j-1.);
                  m = mnext;
                  j = jnext;
                }
              }
              
              h2Bias->Fill(m, j, (iObs == 1 ? m : j) - genValue);
            }
            
            TF2* fit2D = new TF2("fit2D", "[0] + [1]*(x-172.5) + [2]*(y-1) + [3]*(x-172.5)*(y-1)");
            fit2D->SetParNames("offset", "slopeMass", "slopeJES");
            h2Bias->Fit("fit2D");
            
            for (int iPar = 0; iPar < 4; ++iPar) calibration[iter][iObs][i][iPar] = fit2D->GetParameter(iPar);
            h2Bias->Delete();
          }
          // mt1d
          else if (iObs == 2) {
            TH1D* hBias = new TH1D("hBias", "hBias", 1000, 150, 200);
            
            for (it = samples.begin(); it != samples.end(); ++it) {
              if (!it->calibration || it->genJES != 1.) continue;
              
              double genValue = it->genMass;
              double m = it->results[2].profile->GetBinContent(i);
              if (iter > 0) {
                for (int iterprev = 0; iterprev < iter; ++iterprev) {
                  double mnext = m - calibration[iterprev][2][i][0] - calibration[iterprev][2][i][1]*(m-172.5);
                  m = mnext;
                }
              }
              hBias->Fill(m, m - genValue);
            }
            
            TF1* fit = new TF1("fit", "[0] + [1]*(x-172.5)");
            fit->SetParNames("offset", "slopeMass");
            hBias->Fit("fit");
            
            for (int iPar = 0; iPar < 2; ++iPar) calibration[iter][iObs][i][iPar] = fit->GetParameter(iPar);
            for (int iPar = 2; iPar < 4; ++iPar) calibration[iter][iObs][i][iPar] = 0.;
            hBias->Delete();
          }
          // hybrid
          //*
          else if (iObs == 4 || iObs == 5) {
            TH2D* hhBias = new TH2D("hhBias", "hhBias", 1000, 150, 200, 2000, 0.8, 1.2);
            
            for (it = samples.begin(); it != samples.end(); ++it) {
              if (!it->calibration) continue;
              
              double genValue = it->results[iObs].incl; //it->genMass;
              //if (iObs == 5) genValue = it->genJES; //(it->genJES+1.)/2.;
              double m = it->results[4].profile->GetBinContent(i);
              double j = it->results[5].profile->GetBinContent(i); //(it->results[5].profile->GetBinContent(i)-1.)*2. + 1.;
              if (iter > 0) {
                for (int iterprev = 0; iterprev < iter; ++iterprev) {
                  double mnext = m - calibration[iterprev][4][i][0] - calibration[iterprev][4][i][1]*(m-172.5) - calibration[iterprev][4][i][2]*(j-1.) - calibration[iterprev][4][i][3]*(m-172.5)*(j-1.);
                  double jnext = j - calibration[iterprev][5][i][0] - calibration[iterprev][5][i][1]*(m-172.5) - calibration[iterprev][5][i][2]*(j-1.) - calibration[iterprev][5][i][3]*(m-172.5)*(j-1.);
                  m = mnext;
                  j = jnext;
                }
              }
              std::cout << it->id << " " << iObs << ", bin " << i << ", fitted " << m << " " << j << ", generated " << genValue << std::endl;
              
              hhBias->Fill(m, j, (iObs == 4 ? m : j) - genValue);
            }
            
            TF2* fit2D = new TF2("fit2D", "[0] + [1]*(x-172.5) + [2]*(y-1) + [3]*(x-172.5)*(y-1)");
            fit2D->SetParNames("offset", "slopeMass", "slopeJES");
            hhBias->Fit("fit2D");
            
            for (int iPar = 0; iPar < 4; ++iPar) calibration[iter][iObs][i][iPar] = fit2D->GetParameter(iPar);
            hhBias->Delete();
          }
          //*/
        }
      }
    }
  }
  //*/
  
  // Do the (calibrated) difference
  for (it = samples.begin(); it != samples.end(); ++it) {
    if (it->calibration) continue;
    for (int iObs = 0; iObs < nObs; ++iObs) {
      //std::cout << it->genMass << it->genJES << std::endl;
      TString fitid("constFit"); fitid += it->id;
      it->fit = new TF1(fitid, "[0]");
      it->fit->SetLineColor(kBlack);
      //std::cout << it->id << ", incl: " << it->incl << std::endl;
      it->results[iObs].profile->Fit(fitid, "Q0");
      
      it->results[iObs].profile->SetLineColor(it->color);
      it->results[iObs].profile->SetLineWidth(2);
      //it->results[iObs].profile->SetLineStyle(it->line);
      it->results[iObs].profile->SetMarkerStyle(it->marker);
      it->results[iObs].profile->SetMarkerColor(it->color);
      
      it->results[iObs].profileDiff = (TH1F*) it->results[iObs].profile->Clone();
      
      for (int i = 1; i < it->results[iObs].profile->GetNbinsX()+1; i++) {
        if (iObs != 0) { // mass, JES
          double value = it->results[iObs].profile->GetBinContent(i);
          double vaerr = it->results[iObs].profile->GetBinError(i);
          if (doCalibration) {
            int mStart = 1;
            int jStart = 3;
            int mIndex = 1;
            int jIndex = 3;
            int vIndex = iObs;
            if (iObs == 2) {
              mStart = 2;
              mIndex = 2;
            }
            else if (iObs > 3) { // hybrid
              mStart = 4;
              jStart = 5;
              mIndex = 4;
              jIndex = 5;
            }
            double m = it->results[mStart].profile->GetBinContent(i);
            double j = it->results[jStart].profile->GetBinContent(i);
            for (int iter = 0; iter < iterations; ++iter) {
              value = value - calibration[iter][vIndex][i][0] - calibration[iter][vIndex][i][1]*(m-172.5) - calibration[iter][vIndex][i][2]*(j-1.) - calibration[iter][vIndex][i][3]*(m-172.5)*(j-1.); // /2 brings jsfh in ~agreement
              double mnext = m - calibration[iter][mIndex][i][0]/2. - calibration[iter][mIndex][i][1]*(m-172.5) - calibration[iter][mIndex][i][2]*(j-1.) - calibration[iter][mIndex][i][3]*(m-172.5)*(j-1.);
              double jnext = j - calibration[iter][5][i][0]/2. - calibration[iter][jIndex][i][1]*(m-172.5) - calibration[iter][jIndex][i][2]*(j-1.) - calibration[iter][jIndex][i][3]*(m-172.5)*(j-1.);
              m = mnext;
              j = jnext;
              
              if (iObs != 3 && iObs !=5) {
                vaerr = sqrt(pow(vaerr,2) + pow(vaerr * abs(calibration[iter][vIndex][i][1]),2)); // mass slope
              }
              else vaerr = sqrt(pow(vaerr,2) + pow(vaerr * abs(calibration[iter][vIndex][i][2]),2)); // jsf slope
            }
          }
          it->results[iObs].profileDiff->SetBinContent(i, value - it->results[iObs].incl);
          it->results[iObs].profileDiff->SetBinError(i, sqrt(pow(vaerr, 2) - pow(it->fit->GetParError(0), 2)));
        }
        else {           // entries, normalize, double uncertainty due to in-event correlations
          double integral = it->results[iObs].profile->Integral() / 29109.;
          it->results[iObs].profileDiff->SetBinContent(i, it->results[iObs].profile->GetBinContent(i) / (integral * it->results[iObs].profile->GetBinWidth(i)));
          it->results[iObs].profileDiff->SetBinError(i, 2*it->results[iObs].profile->GetBinError(i) / (integral * it->results[iObs].profile->GetBinWidth(i)));
        }
        
        // Prepare Y axis scaling
        //if (it->systematic) continue;
        //if ((iObs != 0 && fabs(it->results[iObs].profile->GetBinContent(i)) > 100) || it->systematic) continue;
        //std::cout << it->id << " " << iObs << " " << it->results[iObs].profile->GetBinContent(i) << " " << it->results[iObs].profile->GetBinError(i) << std::endl;
        double newMin = it->results[iObs].profileDiff->GetBinContent(i) - it->results[iObs].profileDiff->GetBinError(i);
        if (newMin < minVal[iObs]) minVal[iObs] = newMin;
        double newMax = it->results[iObs].profileDiff->GetBinContent(i) + it->results[iObs].profileDiff->GetBinError(i);
        if (newMax > maxVal[iObs]) maxVal[iObs] = newMax;
      }
      
      // error histos
      it->results[iObs].error = (TH1F*)it->results[iObs].profileDiff->Clone();
      //it->results[iObs].error->SetMarkerStyle(0);
      it->results[iObs].error->SetFillColor(it->color-11);
      it->results[iObs].error->SetFillStyle(3254);
      
      // shift MC points
      //TF1 *f1 = new TF1("f1","1+0.1*x",-2,2);
      //it->results[iObs].profileDiff->Multiply(f1,1);
      //Double_t shift = 20.; //0.5*it->results[iObs].profileDiff->GetBinWidth(1);
      //it->results[iObs].profileDiff->GetXaxis()->SetLimits(
      //  it->results[iObs].profileDiff->GetXaxis()->GetXmin()+shift,
      //  it->results[iObs].profileDiff->GetXaxis()->GetXmax()+shift
      //);
      
      //*
      const int nbins = samples[0].results[iObs].profile->GetNbinsX()+1;
      std::cout << "nbins " << nbins << std::endl;
      double bins[nbins];
      double shift = 0.;
      if (it->marker > 20) shift = (it->results[iObs].profileDiff->GetXaxis()->GetXmax() - it->results[iObs].profileDiff->GetXaxis()->GetXmin()) / 100. / nbins * 4. * (it->marker-27.);
      if (it->marker >= 27) shift = (it->results[iObs].profileDiff->GetXaxis()->GetXmax() - it->results[iObs].profileDiff->GetXaxis()->GetXmin()) / 100. / nbins * 4. * (it->marker-27.+1.);
      for (int i = 1; i < nbins+1; i++) {
        bins[i-1] = it->results[iObs].profileDiff->GetXaxis()->GetBinLowEdge(i)+shift;
        std::cout << it->results[iObs].profileDiff->GetXaxis()->GetBinLowEdge(i) << std::endl;
      }
      it->results[iObs].profileDiff->GetXaxis()->Set(nbins-1, bins);
      //*/
    }
  }
  
  // SYSTEMATICS
  //*
  std::vector<TH1F*> hNull;
  
  for (int iObs = 0; iObs < nObs; ++iObs) {
    hNull.push_back((TH1F*) samples[0].results[iObs].profileDiff->Clone("hNull"));
    hNull[iObs]->SetLineColor(kBlack);
    hNull[iObs]->SetLineStyle(1);
    //samples[0].profile->SetLineWidth(4);
    
    int sysRef = 1; // MadSpin sample
    
    for (int i = 1; i < samples[0].results[iObs].profileDiff->GetNbinsX()+1; i++) {
      double syst2 = 0;
      syst2 += pow(max(fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+2].results[iObs].profileDiff->GetBinContent(i)),
                       fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+3].results[iObs].profileDiff->GetBinContent(i))), 2);
      syst2 += pow(max(fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+4].results[iObs].profileDiff->GetBinContent(i)),
                       fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+5].results[iObs].profileDiff->GetBinContent(i))), 2);
      syst2 += pow(max(fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+6].results[iObs].profileDiff->GetBinContent(i)),
                       fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+7].results[iObs].profileDiff->GetBinContent(i))), 2);
      syst2 += pow(max(fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+8].results[iObs].profileDiff->GetBinContent(i)),
                       fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+9].results[iObs].profileDiff->GetBinContent(i))), 2);
      syst2 += pow(max(fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+10].results[iObs].profileDiff->GetBinContent(i)),
                       fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+11].results[iObs].profileDiff->GetBinContent(i))), 2);
      syst2 += pow(max(fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+12].results[iObs].profileDiff->GetBinContent(i)),
                       fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+13].results[iObs].profileDiff->GetBinContent(i))), 2);
      // TODO: use max(up, down, stat)?
      syst2 += pow(max(
                max(fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+14].results[iObs].profileDiff->GetBinContent(i)),
                    fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+15].results[iObs].profileDiff->GetBinContent(i))),
                //max(sqrt(pow(samples[sysRef].results[iObs].profileDiff->GetBinError(i),2) + pow(samples[sysRef+8].results[iObs].profileDiff->GetBinError(i), 2)),
                //    sqrt(pow(samples[sysRef].results[iObs].profileDiff->GetBinError(i),2) + pow(samples[sysRef+9].results[iObs].profileDiff->GetBinError(i), 2)))
                0.), 2);
      syst2 += pow(max(
                max(fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+16].results[iObs].profileDiff->GetBinContent(i)),
                    fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+17].results[iObs].profileDiff->GetBinContent(i))),
                //max(sqrt(pow(samples[sysRef].results[iObs].profileDiff->GetBinError(i),2) + pow(samples[sysRef+10].results[iObs].profileDiff->GetBinError(i), 2)),
                //    sqrt(pow(samples[sysRef].results[iObs].profileDiff->GetBinError(i),2) + pow(samples[sysRef+11].results[iObs].profileDiff->GetBinError(i), 2)))
                 0.), 2);
      syst2 += pow(max(fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+18].results[iObs].profileDiff->GetBinContent(i)),
                       fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+19].results[iObs].profileDiff->GetBinContent(i))), 2);
      syst2 += pow(max(fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+20].results[iObs].profileDiff->GetBinContent(i)),
                       fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+21].results[iObs].profileDiff->GetBinContent(i))), 2);
      syst2 += pow(fabs(samples[sysRef].results[iObs].profileDiff->GetBinContent(i) - samples[sysRef+22].results[iObs].profileDiff->GetBinContent(i)), 2);
      
      hNull[iObs]->SetBinError(i, sqrt(pow(samples[0].results[iObs].profileDiff->GetBinError(i), 2) + syst2));
      if (iObs == 0) hNull[iObs]->SetBinError(i, hNull[iObs]->GetBinError(i) / hNull[iObs]->GetBinWidth(i));
    }
  }
  //*/
  
  // CALCULATE CHI2
  // TODO: for all generator setups
  //*
  //if (batch) {
  for (int iSample = 0; iSample < 30; ++iSample) {
    std::cout << "Comparing data to " << samples[iSample].label << std::endl;
    for (int iObs = 0; iObs < nObs; ++iObs) {
      double chi2 = 0;
      int ndf = samples[0].results[iObs].profileDiff->GetNbinsX()-1;
      for (int i = 1; i < samples[0].results[iObs].profileDiff->GetNbinsX()+1; i++) {
        chi2 += pow((hNull[iObs]->GetBinContent(i) - samples[iSample].results[iObs].profileDiff->GetBinContent(i)), 2) / (pow(hNull[iObs]->GetBinError(i), 2) + pow(samples[iSample].results[iObs].profileDiff->GetBinError(i), 2));
      }
      std::cout << sBinning[iBinning] << " " << sObservable[iObs] << ": "
                << "chi2 = " << chi2 
                << ", ndof = " << ndf 
                << ", prob = " << TMath::Prob(chi2, ndf)
                << std::endl;
      if (iObs == 0) chi2matrix[iBinning][iObs][iSample] = (double) ndf;
      else chi2matrix[iBinning][iObs][iSample] = chi2;
    }
  }
  //}
  //*/
  
  // FINISH
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
  
  for (int iObs = 0; iObs < nObs; ++iObs) {
    double pad2size = 0.3;
    if (iObs != 0 && !doCalibration) {
      TPad *pad1 = new TPad("pad1","pad1",0,pad2size,1,1);
      pad1->SetTopMargin(1./(1.-pad2size)*pad1->GetTopMargin());
      pad1->SetBottomMargin(0);
      pad1->Draw();
      pad1->cd();
    }
    
    //int ref = 0; // DATA
    double range = maxVal[iObs] - minVal[iObs];
    double minRange = minVal[iObs] - range * 0.1;
    double maxRange = maxVal[iObs] + range * 0.5;
    //*
    if (doCalibration && iBinning == 16) {
      if (iObs == 1 || iObs == 2 || iObs == 4) {
        minRange = -2.0;
        maxRange =  4.5;
      }
      if (iObs == 3) {
        minRange = -0.05;
        maxRange =  0.05;
      }
    }
    //*/
    //std::cout << "min: " << minRange << ", max: " << maxRange << std::endl;
    if (iObs != 0) hNull[iObs]->GetYaxis()->SetRangeUser(minRange, maxRange);
    else           hNull[iObs]->GetYaxis()->SetRangeUser(0, maxRange);
    hNull[iObs]->GetXaxis()->SetTitleOffset(1.2);
    hNull[iObs]->GetXaxis()->SetTitle(sBinNice[iBinning]);
    if (doCalibration == false) hNull[iObs]->GetYaxis()->SetTitle(sObsNice[iObs]);
    else                        hNull[iObs]->GetYaxis()->SetTitle(sObsNiceCal[iObs]);
    if (iObs == 0) hNull[iObs]->GetYaxis()->SetTitleOffset(1.75);
    hNull[iObs]->Draw("E");
    
    //TF1* f0 = new TF1("f0","0", 0., 2000.);
    //f0->SetLineColor(kGray+2);
    //f0->Draw("same");
    
    TLegend* leg0 = new TLegend(0.25, 0.73, 0.55, 0.9);
    leg0->SetTextSize(0.03);
    leg0->SetFillStyle(0);
    leg0->SetBorderSize(0);
    
    TLegend *leg1 = new TLegend(0.6, 0.73, 0.9, 0.9);
    leg1->SetTextSize(0.03);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    
    //if (!doCalibration) {
      for (it = samples.begin(); it != samples.end(); ++it) {
        if (!it->data && !it->systematic && !it->calibration && it->drawerror) {
          //it->results[iObs].error->Draw("SAME,E2");
        }
      }
    //}
    
    int legEntry = 0;
    for (it = samples.begin(); it != samples.end(); ++it) {
      if (it->data) {
        it->results[iObs].profileDiff->Draw("SAME,E1,X0");
        ++legEntry;
        leg0->AddEntry(it->results[iObs].profileDiff, it->label, "P");
      }
      else if (!it->systematic && !it->calibration //&& !doCalibration
              ) {
        it->results[iObs].profileDiff->Draw("E,X0,SAME");
        ++legEntry;
        // 2 columns
        if (legEntry > 4) leg1->AddEntry(it->results[iObs].error, it->label, "P");
        else              leg0->AddEntry(it->results[iObs].error, it->label, "P");
      }
    }
    hNull[iObs]->Draw("SAME,E,X0");
    gPad->RedrawAxis();
    
    if (iObs != 0 && !doCalibration) {
      canvas->cd();
      TPad *pad0 = new TPad("pad0","pad0",0,0,1,1);
      pad0->Draw();
      pad0->cd();
      pad0->SetFillStyle(0);
    }
    
    leg0->Draw();
    leg1->Draw();
    
    if (!doCalibration) {
      if (iObs == 1 || iObs == 2) {
        TPaveLabel* gevLabel = DrawLabel("[GeV]", 0.06, 0.26, 0.2);
        gevLabel->SetTextFont(42);
        gevLabel->SetTextAngle(90);
        gevLabel->SetTextSize(0.85);
      }
    }
    CMS_lumi(canvas, 21, 0.);
    
    if (fixme) DrawLabel("WARNING: Potential errors detected", 0.2, 0.00, 0.9);
    
    // MC-DATA
    if (iObs != 0 && !doCalibration) {
      canvas->cd();
      TPad *pad2 = new TPad("pad2","pad2",0,0,1,pad2size);
      pad2->SetTopMargin(0);
      pad2->SetBottomMargin(1./pad2size*pad2->GetBottomMargin());
      pad2->Draw();
      pad2->cd();
      
      TF1* null = new TF1("null", "0", 0, 1000);
      null->SetLineColor(kBlack); null->SetLineStyle(1);
      
      TH1F* hNull2 = (TH1F*) hNull[iObs]->Clone();
      
      hNull[iObs]->GetYaxis()->SetLabelSize(1./(1.-pad2size)*hNull[iObs]->GetYaxis()->GetLabelSize());
      //hNull[iObs]->GetYaxis()->SetTickLength(1./(1.-pad2size)*hNull[iObs]->GetYaxis()->GetTickLength());
      hNull[iObs]->GetYaxis()->SetTitleSize(1./(1.-pad2size)*hNull[iObs]->GetYaxis()->GetTitleSize());
      hNull[iObs]->GetYaxis()->SetTitleOffset(1.0);
      hNull[iObs]->GetXaxis()->SetLabelSize(1./(1.-pad2size)*hNull[iObs]->GetXaxis()->GetLabelSize());
      hNull[iObs]->GetXaxis()->SetTickLength(1./(1.-pad2size)*hNull[iObs]->GetXaxis()->GetTickLength());
      hNull[iObs]->GetXaxis()->SetTitleSize(5./4.*1./(1.-pad2size)*hNull[iObs]->GetXaxis()->GetTitleSize());
      
      double minVal2 = 0;
      double maxVal2 = 0;
      for (int i = 1; i < samples[0].results[iObs].profileDiff->GetNbinsX()+1; i++) {
        double val = hNull2->GetBinContent(i) - samples[1].results[iObs].profileDiff->GetBinContent(i);
        double error = hNull2->GetBinError(i);
        hNull2->SetBinContent(i, val);
        
        if (val-error < minVal2) minVal2 = val-error;
        if (val+error > maxVal2) maxVal2 = val+error;
      }
      
      double range2 = 0.015;
      if (max(abs(minVal2), maxVal2) > 0.014) range2 = 0.025;
      if (max(abs(minVal2), maxVal2) > 0.024) range2 = 0.03;
      if (max(abs(minVal2), maxVal2) > 0.029) range2 = 0.06;
      if (max(abs(minVal2), maxVal2) > 0.059) range2 = 0.12;
      if (max(abs(minVal2), maxVal2) > 0.119) range2 = 1.5;
      if (max(abs(minVal2), maxVal2) > 1.4) range2 = 2.5;
      if (max(abs(minVal2), maxVal2) > 2.4) range2 = 6;
      if (max(abs(minVal2), maxVal2) > 5.9) range2 = 12;
      
      hNull2->GetYaxis()->SetRangeUser(-range2, range2);
      hNull2->GetYaxis()->SetTitle("data - MG Z2* ");
      hNull2->GetYaxis()->SetLabelSize(1./pad2size*hNull2->GetYaxis()->GetLabelSize());
      //hNull2->GetYaxis()->SetTickLength(1./pad2size*hNull2->GetYaxis()->GetTickLength());
      hNull2->GetYaxis()->SetTitleSize(5./4.*1./pad2size*hNull2->GetYaxis()->GetTitleSize()*0.5);
      hNull2->GetYaxis()->SetTitleOffset(1.75*pad2size*hNull2->GetYaxis()->GetTitleOffset());
      hNull2->GetYaxis()->SetNdivisions(503);
      hNull2->GetXaxis()->SetLabelSize(1./pad2size*hNull2->GetXaxis()->GetLabelSize());
      hNull2->GetXaxis()->SetTickLength(1./pad2size*hNull2->GetXaxis()->GetTickLength());
      hNull2->GetXaxis()->SetTitleSize(5./4.*1./pad2size*hNull2->GetXaxis()->GetTitleSize()*0.8);
      
      //hNull2->SetFillColor(kYellow);
      //hNull2->SetMarkerColor(kYellow);
      hNull2->Draw("E,X0");
      
      for (it = samples.begin(); it != samples.begin()+2; ++it) {
        it->results[iObs].profile2 = (TH1F*) it->results[iObs].error->Clone();
        for (int i = 1; i < samples[0].results[iObs].profileDiff->GetNbinsX()+1; i++) {
          it->results[iObs].profile2->SetBinContent(i, it->results[iObs].profileDiff->GetBinContent(i) - samples[1].results[iObs].profileDiff->GetBinContent(i));
        }
        if (it->data) {
          it->results[iObs].profile2->SetMarkerStyle(20);
          it->results[iObs].profile2->Draw("SAME,E1,X0,A");
        }
        else {
          it->results[iObs].profile2->Draw("SAME,E2");
          it->results[iObs].profile2->Draw("SAME,E1");
        }
      }
      //null->Draw("SAME");
      hNull2->Draw("SAME,E,X0");
    }
    
    canvas->Update();
    
    TString path("../plot/differential/diff_" + sBinning[iBinning] + "_" + sObsFile[iObs] + (doCalibration ? "_cal" : "") + ".eps");
    canvas->Print(path);
    TString pathr("../plot/differential/root/diff_" + sBinning[iBinning] + "_" + sObsFile[iObs] + (doCalibration ? "_cal" : "") + ".root");
    canvas->Print(pathr);
    canvas->Clear();
  }
}

void differentialBatch() {
  for (unsigned int b = 0; b < sizeof(sBinning) / sizeof(sBinning[0]); ++b) {    
    // 1D mass, JES
    differentialMass(b);
  }
  
  for (int iSample = 0; iSample < 30; ++iSample) {
    double chi2_m1d = 0.;
    double chi2_m2d = 0.;
    double chi2_mh  = 0.;
    double chi2_jes = 0.;
    double sum_ndf  = 0.;
    for (unsigned int b = 0; b < 14; ++b) {
      if (!paper[b]) continue;
      printf("%s & %2.2f & %2.2f & %2.2f & %2.2f & %1.0f \\\\ \n", sBinLatex[b].Data(), chi2matrix[b][2][iSample], chi2matrix[b][3][iSample], chi2matrix[b][1][iSample], chi2matrix[b][4][iSample], chi2matrix[b][0][iSample]);
      
      chi2_m1d += chi2matrix[b][2][iSample];
      chi2_jes += chi2matrix[b][3][iSample];
      chi2_m2d += chi2matrix[b][1][iSample];
      chi2_mh  += chi2matrix[b][4][iSample];
      sum_ndf  += chi2matrix[b][0][iSample];
    }
    
    printf("%s & %2.2f & %2.2f & %2.2f & %2.2f & %1.0f & %2.2f & %2.2f & %2.2f \\\\ \n\n", sampleNames[iSample].Data(), chi2_m1d, chi2_jes, chi2_m2d, chi2_mh, sum_ndf, sqrt(2)*TMath::ErfInverse(1-TMath::Prob(chi2_m1d+chi2_jes, 2*sum_ndf)), sqrt(2)*TMath::ErfInverse(1-TMath::Prob(chi2_m2d, sum_ndf)), sqrt(2)*TMath::ErfInverse(1-TMath::Prob(chi2_mh, sum_ndf)));
  }
}

