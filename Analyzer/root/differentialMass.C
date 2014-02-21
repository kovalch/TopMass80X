#include <vector>
#include <iostream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TF1.h"
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

#include "tdrstyle.C"

TH1F* hData;
TH1F* hTTJets096;
TH1F* hTTJets100;
TH1F* hTTJets104;
TH1F* hTTJetsFlavorUp;
TH1F* hTTJetsFlavorDown;
TH1F* hTTJetsJESUp;
TH1F* hTTJetsJESDown;
TH1F* hTTJetsJERUp;
TH1F* hTTJetsJERDown;
TH1F* hTTJetsScaleUp;
TH1F* hTTJetsScaleDown;
TH1F* hTTJetsMatchingUp;
TH1F* hTTJetsMatchingDown;
TH1F* hTTJetsP11;
TH1F* hTTJetsP11noCR;
TH1F* hTTJetsMCatNLO;
TH1F* hTTJetsPowheg;

enum lepton          { kElectron, kMuon, kAll};
TString lepton_ [] = { "electron", "muon", "all"};
enum cat             { kCR, kRad, kBoth};

TString sBinning[] = {"TopPt", "TopEta", "B1Pt", "B1Eta", "TTBarMass", "TTBarPt", "DeltaRqq", "DeltaRbb", "HT", "NJet", "Q1Pt", "Q1Eta", "WPt", "WEta"};
TString sData[] = {"top_fitTop1_Pt__", "top_fitTop1_Eta__", "top_fitB1_Pt__", "top_fitB1_Eta__", "top_fitTTBar_M__", "top_fitTTBar_Pt__", "sqrt_pow_top_fitW1Prod1_Eta__-top_fitW1Prod2_Eta___2__+_pow_TVector2__Phi_mpi_pi_top_fitW1Prod1_Phi__-top_fitW1Prod2_Phi____2__", "sqrt_pow_top_fitB1_Eta__-top_fitB2_Eta___2__+_pow_TVector2__Phi_mpi_pi_top_fitB1_Phi__-top_fitB2_Phi____2__", "top_fitB1_Pt__+top_fitB1_Pt__+top_fitW1Prod1_Pt__+top_fitW1Prod2_Pt__", "jet__jet_size__", "top_fitW1Prod1_Pt__", "top_fitW1Prod1_Eta__", "top_fitW1_Pt__", "top_fitW1_Eta__"};
TString sBinNice[] = {"p_{T,t,had} [GeV]", "|#eta_{t,had}|", "p_{T,b,had} [GeV]", "|#eta_{b,had}|", "m_{t#bar{t}} [GeV]", "p_{T,t#bar{t}} [GeV]", "#DeltaR_{q#bar{q}}", "#DeltaR_{b#bar{b}}", "H_{T}^{4} [GeV]", "Number of jets", "p_{T,q} [GeV]", "|#eta_{q}|", "p_{T,W} [GeV]", "|#eta_{W}|"};
int binCat[] = {/*"p_{T,t,had} [GeV]"*/ kCR, /*"|#eta_{t,had}|"*/ kCR, /*"p_{T,b,had} [GeV]"*/ kBoth, /*"|#eta_{b,had}|"*/ kBoth, /*"m_{t#bar{t}} [GeV]"*/ kRad, /*"p_{T,t#bar{t}} [GeV]"*/ kRad, /*"#DeltaR_{q#bar{q}}"*/ kCR, /*"#DeltaR_{b#bar{b}}"*/ kBoth, /*"p_{T,b,lep} [GeV]"*/ kBoth, /*"|#eta_{b,lep}|"*/ kBoth, /*"#Delta#phi_{q#bar{q}}"*/ kCR, /*"#Delta#phi_{b#bar{b}}"*/ kBoth, /*"p_{T,lep} [GeV]"*/ kRad, /*"|#eta_{lep}|"*/ kRad, /*"MET [GeV]"*/ kRad, /*"H_{T} [GeV]"*/ kRad, /*"Number of jets"*/ kRad, /*"Number of b-tagged jets"*/ kRad, kBoth, kBoth};
bool pas[] = {/*"p_{T,t,had} [GeV]"*/ true, /*"|#eta_{t,had}|"*/ true, /*"p_{T,b,had} [GeV]"*/ true, /*"|#eta_{b,had}|"*/ true, /*"m_{t#bar{t}} [GeV]"*/ true, /*"p_{T,t#bar{t}} [GeV]"*/ true, /*"#DeltaR_{q#bar{q}}"*/ true, /*"#DeltaR_{b#bar{b}}"*/ true, /*"p_{T,b,lep} [GeV]"*/ false, /*"|#eta_{b,lep}|"*/ false, /*"#Delta#phi_{q#bar{q}}"*/ true, /*"#Delta#phi_{b#bar{b}}"*/ true, /*"p_{T,lep} [GeV]"*/ false, /*"|#eta_{lep}|"*/ false, /*"MET [GeV]"*/ false, /*"H_{T} [GeV]"*/ true, /*"Number of jets"*/ true, /*"Number of b-tagged jets"*/ false, false, false};

TString sObservable[4] = {"Entries", "mass_mTop_JES", "mass_mTop", "JES_mTop_JES"};
TString sObsLowCase[4] = {"Entries", "mass_mTop_JES", "mass_mTop", "JES_mTop_JES"};
TString sObsFile[4]    = {"Entries", "MT2D", "MT1D", "JES"};
TString sObsCanvas[4] = {"canvas_1", "canvas_2", "canvas_4", "canvas_3"};
TString sObsNice[4] = {"Number of permutations / bin width", "m_{t}^{2D} - <m_{t}^{2D}> [GeV]", "m_{t}^{1D} - <m_{t}^{1D}> [GeV]", "JES - <JES>"};
TString sObsNice2[4] = {"Number of permutations / bin width", "(m_{t}^{2D} - <m_{t}^{2D}>)_{MC} - (m_{t}^{2D} - <m_{t}^{2D}>)_{data} [GeV]", "(m_{t}^{1D} - <m_{t}^{1D}>)_{MC} - (m_{t}^{1D} - <m_{t}^{1D}>)_{data} [GeV]", "(JES - <JES>)_{MC} - (JES - <JES>)_{data}"};

double chi2matrix[14][4];

double crossSection   = 164.4;
double peLumi         = 5000.;

int channel = 2;

int nbins = 10;
double xfirst = 0;
double xlast  = 1;

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

struct sample {
  bool data;
  bool systematic;
  bool cr;
  bool rad;
  bool drawerror;
  TChain* chain;
  TFile* file;
  TCanvas* canvas;
  TPad* pad;
  TH1F* profile;
  TH1F* profile2;
  double incl;
  //std::vector<double> vector;
  TH1F* error;
  TF1* fit;
  TGraphErrors* graph;
  const char* id;
  const char* label;
  int color;
  int line;
  double size; // EFFECTIVE sample size
  sample(bool d, bool sys, bool c_, bool r, bool de, const char* i, const char* l, int c = kBlack, int li = 1, double s = 0)
  : data(d), systematic(sys), cr(c_), rad(r), drawerror(de), id(i), label(l), color(c), line(li), size(s) {}
};

string itos(int number)
{
   stringstream ss;
   ss << number;
   return ss.str();
}

void differentialMass(int iBinning = 0, int iObs = 2, bool batch = false)
{
  bool fixme = false;
  
  TStyle* tdrStyle = setTDRStyle();
  tdrStyle->SetPadLeftMargin(0.2);
  tdrStyle->SetPadRightMargin(0.05);
  tdrStyle->SetNdivisions(505, "X");
  tdrStyle->SetTitleYOffset(1.75);
  tdrStyle->SetOptStat(0);
  tdrStyle->SetErrorX(0.5);
  tdrStyle->SetEndErrorSize(4);
  
  std::vector<sample> samples;
  std::vector<sample>::iterator it;
  //std::vector<std::vector<double> > values;
  
  samples.push_back(sample(true,  false, true, true, false, "Run2012", "Data", kBlack, 0));  
  samples.push_back(sample(false, false, true, true, true,  "1.00", "MG, Pythia Z2*", kRed+1, 0, 7000000./1.7));
  
  samples.push_back(sample(false, true,  false, false, false, "flavor:up_FlavorQCD", "MG, b-JES up", kBlue+1, 2, 59613991./1.7));
  samples.push_back(sample(false, true,  false, false, false, "flavor:down_FlavorQCD", "MG, b-JES down", kGreen+1, 2, 59613991./1.7));
  samples.push_back(sample(false, true,  false, false, false, "source:up_Total", "MG, JES up", kBlue+1, 2, 59613991./1.7));
  samples.push_back(sample(false, true,  false, false, false, "source:down_Total", "MG, JES down", kGreen+1, 2, 59613991./1.7));
  samples.push_back(sample(false, true,  false, false, false, "jer:up", "MG, JER up", kBlue+1, 2, 59613991./1.7));
  samples.push_back(sample(false, true,  false, false, false, "jer:down", "MG, JER down", kGreen+1, 2, 59613991./1.7));
  samples.push_back(sample(false, true,  false, false, false, "source:up_PileUpPtBB", "MG, PUJES up", kBlue+1, 2, 59613991./1.7));
  samples.push_back(sample(false, true,  false, false, false, "source:down_PileUpPtBB", "MG, PUJES down", kGreen+1, 2, 59613991./1.7));
  
  /*
  samples.push_back(sample(false, false, false, false, false, "scaleup", "MG, scale up", kRed-1, 5, 3696269./1.7));
  samples.push_back(sample(false, false, false, false, false, "scaledown", "MG, scale down", kRed-7, 6, 4004587./1.7));
  samples.push_back(sample(false, false, false, false, false, "matchingup", "MG, matching up", kRed-1, 7, 4029823./1.7));
  samples.push_back(sample(false, false, false, false, false, "matchingdown", "MG, matching down", kRed-7, 7, 1545688./1.7));
  //*/
  
  //*
  samples.push_back(sample(false, false, true, true, true, "MGDecays_P11", "MG, Pythia P11", kMagenta+1, 5, 27000000./1.7));
  samples.push_back(sample(false, false, true, true, true, "MGDecays_P11noCR", "MG, Pythia P11noCR", kCyan+1, 6, 27000000./1.7));
  samples.push_back(sample(false, false, true, true, true, "powheg", "Powheg, Pythia Z2*", kGreen+1, 9, 21675970./1.7));
  samples.push_back(sample(false, false, true, true, true, "powheg_herwig", "Powheg, Herwig 6", kOrange+2, 7, 27684235./1.7));
  samples.push_back(sample(false, false, true, true, true, "mcatnlo_herwig", "MC@NLO, Herwig 6", kBlue+1, 2, 32852589./1.7));
  samples.push_back(sample(false, false, true, true, true, "sherpa", "Sherpa", kYellow+1, 3, 44000000./1.7));
  //*/
  //
  
  double minVal = 0;
  double maxVal = 0;
  
  TH1F* reference;
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    if (it->data) {
      TFile* inclFile     = new TFile("/afs/desy.de/user/m/mseidel/xxl/CMSSW_5_3_11/src/TopMass/Analyzer/plot/IdeogramMin_Run2012_top_fitTop1_0__M__.root");
      TCanvas* inclCanvas = (TCanvas*) inclFile  ->Get("canvas");
      TPad* inclPad       = (TPad*)    inclCanvas->GetPrimitive(sObsCanvas[iObs]);
      TH1F* inclProfile   = (TH1F*)    inclPad   ->GetPrimitive(sObservable[iObs]);
      it->incl            = inclProfile->GetBinContent(1);
      inclFile->Close();
      
      it->file    = new TFile("/afs/desy.de/user/m/mseidel/xxl/CMSSW_5_3_11/src/TopMass/Analyzer/plot/IdeogramMin_Run2012_" + sData[iBinning] + ".root");
      it->canvas  = (TCanvas*) it->file  ->Get("canvas");
      it->pad     = (TPad*)    it->canvas->GetPrimitive(sObsCanvas[iObs]);
      it->profile = (TH1F*)    it->pad   ->GetPrimitive(sObservable[iObs]);
      if (iBinning == 16 || iBinning ==17) it->profile->SetNdivisions(203);
      reference = (TH1F*) it->profile->Clone();
      /* TODO remove?
      for (int i = 0; i < it->profile->GetNbinsX()+2; i++) {
        if (it->profile->GetBinContent(i) == -1 || ((iObs == 1 || iObs ==2) && it->profile->GetBinContent(i) < 100)) {
          //it->profile->SetBinContent(i, 400);
          it->profile->SetBinError(i, 1e6);
        }
      }
      */
    }
    else {
      TString sfBase("/nfs/dust/cms/user/mseidel/pseudoexperiments/topmass_131012_diff/");
      
      int minFit = 0;
      int maxFit = 100000;
      if (iObs == 1 || iObs == 2) {
        minFit = 100;
        maxFit = 250;
      }
      if (iObs == 3) maxFit = 2;
      
      TChain* inclChain = new TChain("tree");
      inclChain->Add(sfBase + "lepton/Summer12_TTJets1725_" + it->id + "/TopMass/job_*_ensemble.root");
      TF1* inclGaus = new TF1("inclGaus", "gaus");
      inclChain->Fit("inclGaus", sObsLowCase[iObs], sObsLowCase[iObs] + " < " + itos(maxFit) + " & " + sObsLowCase[iObs] + " > " + itos(minFit), "LEMQ");
      it->incl = inclGaus->GetParameter(1);
      inclChain->Delete();
      
      it->chain   = new TChain("tree");
      it->chain->Add(sfBase + "lepton/Summer12_TTJets1725_" + it->id + "/" + sBinning[iBinning] + "/job_*_ensemble.root");
      
      it->profile = (TH1F*) samples[0].profile->Clone(it->id);
      
      for (int i = 1; i < it->profile->GetNbinsX()+1; i++) {
        TF1* gaus = new TF1("gaus", "gaus");
        int entries = it->chain->GetEntries(sObsLowCase[iObs] + " < " + itos(maxFit) + " & " + sObsLowCase[iObs] + " > " + itos(minFit) + " & bin == " + itos(i));
        if (entries < 100) {
            std::cout << "WARNING: Too less entries: " << itos(entries) << " in " << it->id << " bin " << itos(i) << std::endl;
            fixme = true;
        }
        if (entries > 10) {
          it->chain->Fit("gaus", sObsLowCase[iObs], sObsLowCase[iObs] + " < " + itos(maxFit) + " & " + sObsLowCase[iObs] + " > " + itos(minFit) + " & bin == " + itos(i), "LEMQ");
          it->profile->SetBinContent(i, gaus->GetParameter(1));
          it->profile->SetBinError(i, gaus->GetParameter(2) / sqrt(it->size/(crossSection*peLumi)));
        }
        else {
          it->profile->SetBinContent(i, -1);
          it->profile->SetBinError(i, 1e6);
        }
      }
    }
    
    TString fitid("constFit"); fitid += it->id;
    it->fit = new TF1(fitid, "[0]");
    it->fit->SetLineColor(kBlack);
    //std::cout << it->id << ", incl: " << it->incl << std::endl;
    it->profile->Fit(fitid, "Q0");
    
    it->profile->SetLineColor(it->color);
    it->profile->SetLineWidth(2);
    it->profile->SetLineStyle(it->line);
    if (it->data) it->profile->SetMarkerStyle(20);
    else          it->profile->SetMarkerStyle(1);
    
    double integral = it->profile->Integral() / 29109.;
    
    for (int i = 1; i < it->profile->GetNbinsX()+1; i++) {
      if (iObs != 0) { // mass, JES
        it->profile->SetBinContent(i, it->profile->GetBinContent(i) - it->incl);
        it->profile->SetBinError(i, sqrt(pow(it->profile->GetBinError(i), 2) - pow(it->fit->GetParError(0), 2)));
        /* TODO ratio?
        it->profile->SetBinContent(i, it->profile->GetBinContent(i) / reference->GetBinContent(i));
        it->profile->SetBinError(i, it->profile->GetBinError(i) / reference->GetBinContent(i));
        //*/
        //it->vector.push_back(it->profile->GetBinContent(i));
      }
      else {           // entries, double uncertainty due to in-event correlations
        //it->vector.push_back(it->profile->GetBinContent(i)/(integral*5174.));
        it->profile->SetBinContent(i, it->profile->GetBinContent(i) / (integral * it->profile->GetBinWidth(i)));
        it->profile->SetBinError(i, 2*it->profile->GetBinError(i) / (integral * it->profile->GetBinWidth(i)));
      }
      
      // Prepare scaling
      if ((iObs != 0 && fabs(it->profile->GetBinContent(i)) > 100) || it->systematic) continue;
      double newMin = it->profile->GetBinContent(i) - it->profile->GetBinError(i);
      if (newMin < minVal) minVal = newMin;
      double newMax = it->profile->GetBinContent(i) + it->profile->GetBinError(i);
      if (newMax > maxVal) maxVal = newMax;
    }
    
    it->error = (TH1F*)it->profile->Clone();
    it->error->SetMarkerStyle(0);
    it->error->SetFillColor(it->color-11);
    it->error->SetFillStyle(3254);
        
    //values.push_back(it->vector);
  }
  
  // SYSTEMATICS
  TH1F* hNull = (TH1F*) samples[0].profile->Clone("hNull");
  hNull->SetLineColor(kBlack);
  hNull->SetLineStyle(1);
  //samples[0].profile->SetLineWidth(4);
  
  //*
  for (int i = 1; i < samples[0].profile->GetNbinsX()+1; i++) {
    double syst2 = 0;
    syst2 += pow(max(fabs(samples[1].profile->GetBinContent(i) - samples[2].profile->GetBinContent(i)),
                     fabs(samples[1].profile->GetBinContent(i) - samples[3].profile->GetBinContent(i))), 2);
    syst2 += pow(max(fabs(samples[1].profile->GetBinContent(i) - samples[4].profile->GetBinContent(i)),
                      fabs(samples[1].profile->GetBinContent(i) - samples[5].profile->GetBinContent(i))), 2);
    syst2 += pow(max(fabs(samples[1].profile->GetBinContent(i) - samples[6].profile->GetBinContent(i)),
                     fabs(samples[1].profile->GetBinContent(i) - samples[7].profile->GetBinContent(i))), 2);
    //* PU (wrong NJet)
    if (iBinning != 9) {
    syst2 += pow(max(fabs(samples[1].profile->GetBinContent(i) - samples[8].profile->GetBinContent(i)),
                     fabs(samples[1].profile->GetBinContent(i) - samples[9].profile->GetBinContent(i))), 2);
    }
    //*/
    /*
    syst2 += pow(fabs(samples[1].profile->GetBinContent(i) - samples[8].profile->GetBinContent(i)), 2);
    syst2 += pow(max(fabs(samples[1].profile->GetBinContent(i) - samples[9].profile->GetBinContent(i)),
                     fabs(samples[1].profile->GetBinContent(i) - samples[10].profile->GetBinContent(i))), 2);
    syst2 += pow(max(fabs(samples[1].profile->GetBinContent(i) - samples[11].profile->GetBinContent(i)),
                     fabs(samples[1].profile->GetBinContent(i) - samples[12].profile->GetBinContent(i))), 2);
    */
    hNull->SetBinError(i, sqrt(pow(samples[0].profile->GetBinError(i), 2) + syst2));
    if (iObs == 0) hNull->SetBinError(i, hNull->GetBinError(i) / hNull->GetBinWidth(i));
  }
  //*/
  
  // CALCULATE CHI2
  if (batch) {
    double chi2 = 0;
    int ndf = samples[0].profile->GetNbinsX()-1;
    for (int i = 1; i < samples[0].profile->GetNbinsX()+1; i++) {
      chi2 += pow((hNull->GetBinContent(i) - samples[1].profile->GetBinContent(i)), 2) / (pow(hNull->GetBinError(i), 2) + pow(samples[1].profile->GetBinError(i), 2));
    }
    std::cout << sBinning[iBinning] << " " << sObservable[iObs] << ": "
              << "chi2 = " << chi2 
              << ", ndof = " << ndf 
              << ", prob = " << TMath::Prob(chi2, ndf)
              << std::endl;
    if (iObs == 0) chi2matrix[iBinning][iObs] = (double) ndf;
    else chi2matrix[iBinning][iObs] = chi2;
  }
  
  // FINISH
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
  double pad2size = 0.3;
  if (iObs != 0) {
    TPad *pad1 = new TPad("pad1","pad1",0,pad2size,1,1);
    pad1->SetTopMargin(1./(1.-pad2size)*pad1->GetTopMargin());
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
  }
  
  //int ref = 0; // DATA
  double range = maxVal - minVal;
  double minRange = minVal - range * 0.1;
  double maxRange = maxVal + range * 0.55;
  if (iObs != 0) hNull->GetYaxis()->SetRangeUser(minRange, maxRange);
  else hNull->GetYaxis()->SetRangeUser(0, maxRange);
  hNull->GetXaxis()->SetTitle(sBinNice[iBinning]);
  hNull->GetYaxis()->SetTitle(sObsNice[iObs]);
  if (iObs == 0) hNull->GetYaxis()->SetTitleOffset(1.75);
  hNull->Draw("E,X0");
  
  TLegend* leg0 = new TLegend(0.25, 0.75, 0.55, 0.925);
  leg0->SetTextSize(0.03);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  
  TLegend *leg1 = new TLegend(0.6, 0.75, 0.9, 0.925);
  leg1->SetTextSize(0.03);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    if (!it->data && !it->systematic && it->drawerror
        && ( (binCat[iBinning] == kCR && it->cr)
            || (binCat[iBinning] == kRad && it->rad)
            || binCat[iBinning] == kBoth ) ) {
      it->error->Draw("SAME,E2");
    }
  }
  
  int legEntry = 0;
  for (it = samples.begin(); it != samples.end(); ++it) {
    if (it->data) {
      it->profile->Draw("SAME,E1,X0");
      ++legEntry;
      leg0->AddEntry(it->profile, it->label, "P");
    }
    else if (!it->systematic
        && ( (binCat[iBinning] == kCR && it->cr)
            || (binCat[iBinning] == kRad && it->rad)
            || binCat[iBinning] == kBoth ) ) {
      it->profile->Draw("SAME");
      ++legEntry;
      if (legEntry > 4) leg1->AddEntry(it->error, it->label, "LF");
      else              leg0->AddEntry(it->error, it->label, "LF");
    }
  }
  hNull->Draw("SAME,E,X0");
  gPad->RedrawAxis();
  
  if (iObs != 0) {
    canvas->cd();
    TPad *pad0 = new TPad("pad0","pad0",0,0,1,1);
    pad0->Draw();
    pad0->cd();
    pad0->SetFillStyle(0);
  }
  
  leg0->Draw();
  leg1->Draw();
  DrawCMSPrel8LeptonJets();
  if (iObs == 1 || iObs == 2) {
    TPaveLabel* gevLabel = DrawLabel("[GeV]", 0.06, 0.26, 0.2);
    gevLabel->SetTextFont(42);
    gevLabel->SetTextAngle(90);
    gevLabel->SetTextSize(0.85);
  }
  
  if (fixme) DrawLabel("WARNING: Potential errors detected", 0.2, 0.00, 0.9);
  
  // MC-DATA
  if (iObs != 0) {
    canvas->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,pad2size);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(1./pad2size*pad2->GetBottomMargin());
    pad2->Draw();
    pad2->cd();
    
    TF1* null = new TF1("null", "0", 0, 1000);
    null->SetLineColor(kBlack); null->SetLineStyle(1);
    
    TH1F* hNull2 = (TH1F*) hNull->Clone();
    
    hNull->GetYaxis()->SetLabelSize(1./(1.-pad2size)*hNull->GetYaxis()->GetLabelSize());
    //hNull->GetYaxis()->SetTickLength(1./(1.-pad2size)*hNull->GetYaxis()->GetTickLength());
    hNull->GetYaxis()->SetTitleSize(1./(1.-pad2size)*hNull->GetYaxis()->GetTitleSize());
    hNull->GetYaxis()->SetTitleOffset(1.0);
    hNull->GetXaxis()->SetLabelSize(1./(1.-pad2size)*hNull->GetXaxis()->GetLabelSize());
    hNull->GetXaxis()->SetTickLength(1./(1.-pad2size)*hNull->GetXaxis()->GetTickLength());
    hNull->GetXaxis()->SetTitleSize(5./4.*1./(1.-pad2size)*hNull->GetXaxis()->GetTitleSize());
    
    double minVal2 = 0;
    double maxVal2 = 0;
    for (int i = 1; i < samples[0].profile->GetNbinsX()+1; i++) {
      double val = hNull2->GetBinContent(i) - samples[1].profile->GetBinContent(i);
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
      it->profile2 = (TH1F*) it->error->Clone();
      for (int i = 1; i < samples[0].profile->GetNbinsX()+1; i++) {
        it->profile2->SetBinContent(i, it->profile->GetBinContent(i) - samples[1].profile->GetBinContent(i));
      }
      if (it->data) {
        it->profile2->SetMarkerStyle(20);
        it->profile2->Draw("SAME,E1,X0,A");
      }
      if (!it->data && !it->systematic
          && ( (binCat[iBinning] == kCR && it->cr)
              || (binCat[iBinning] == kRad && it->rad)
              || binCat[iBinning] == kBoth ) ) {
        it->profile2->Draw("SAME,E2");
        it->profile2->Draw("SAME,E1");
      }
    }
    //null->Draw("SAME");
    hNull2->Draw("SAME,E,X0");
  }
  
  canvas->Update();
  
  TString path("../plot/differential/diff_" + sBinning[iBinning] + "_" + sObsFile[iObs] + ".eps");
  canvas->Print(path);
  
  // CLEANUP
  if (batch == true) {
    canvas->Clear();
    leg1->Clear();
    delete canvas;
    delete leg1;
    
    for (it = samples.begin(); it != samples.end(); ++it) {
      if (it->data) {
        it->file->Close();
      }
      else {
        it->chain->Delete();
      }
    }
  }
  
  
  //return values;
}

void differentialBatch() {
  for (unsigned int b = 0; b < sizeof(sBinning) / sizeof(sBinning[0]); ++b) {    
    // 1D mass, JES
    differentialMass(b, 0, true);
    differentialMass(b, 1, true);
    differentialMass(b, 2, true);
    differentialMass(b, 3, true);
  }
  
  double chi2_m1d = 0.;
  double chi2_m2d = 0.;
  double chi2_jes = 0.;
  double sum_ndf  = 0.;
  for (unsigned int b = 0; b < sizeof(sBinning) / sizeof(sBinning[0]); ++b) {
    printf("$%s$ & %2.2f & %2.2f & %2.2f & %1.0f \\\\ \n", sBinNice[b].Data(), chi2matrix[b][2], chi2matrix[b][3], chi2matrix[b][1], chi2matrix[b][0]);
    if (true) {
      chi2_m1d += chi2matrix[b][2];
      chi2_jes += chi2matrix[b][3];
      chi2_m2d += chi2matrix[b][1];
      sum_ndf  += chi2matrix[b][0];
    }
  }
  
  printf("PAS chi2 & %2.2f & %2.2f & %2.2f & %1.0f \\\\ \n", chi2_m1d, chi2_jes, chi2_m2d, sum_ndf);
}

/*
void differentialBatch() {
  for (unsigned int b = 0; b < sizeof(sBinning) / sizeof(sBinning[0]); ++b) {    
    // Entries, 2d mass: Calculate systematic uncertainty by reweighting
    std::vector<std::vector<double> > entries = differentialMass(b, 0, true);
    std::vector<std::vector<double> > mass    = differentialMass(b, 1, true);
    
    double massShift = 0;
    for (unsigned int i = 0; i < entries[0].size(); ++i) {
      //std::cout << "Data -- entries: " << entries[0][i] << " dmt: " << mass[0][i] << std::endl;
      //std::cout << "MGZ2 -- entries: " << entries[1][i] << " dmt: " << mass[1][i] << std::endl;
      massShift += (entries[0][i]/entries[1][i] - 1.) * mass[1][i] * entries[1][i];
    }
    std::cout << "Mass shift for reweighting to data: " << massShift << std::endl;
    
    // 1D mass, JES
    differentialMass(b, 2, true);
    differentialMass(b, 3, true);
  }
  
  double chi2_m1d = 0.;
  double chi2_m2d = 0.;
  double chi2_jes = 0.;
  double sum_ndf  = 0.;
  for (unsigned int b = 0; b < sizeof(sBinning) / sizeof(sBinning[0]); ++b) {
    printf("$%s$ & %2.2f & %2.2f & %2.2f & %1.0f \\\\ \n", sBinNice[b].Data(), chi2matrix[b][2], chi2matrix[b][3], chi2matrix[b][1], chi2matrix[b][0]);
    if (pas[b]) {
      chi2_m1d += chi2matrix[b][2];
      chi2_jes += chi2matrix[b][3];
      chi2_m2d += chi2matrix[b][1];
      sum_ndf  += chi2matrix[b][0];
    }
  }
  
  printf("PAS chi2 & %2.2f & %2.2f & %2.2f & %1.0f \\\\ \n", chi2_m1d, chi2_jes, chi2_m2d, sum_ndf);
}

void differentialSys(int iBinning = 0) {
  std::vector<std::vector<double> > entries = differentialMass(iBinning, 0, true);
  std::vector<std::vector<double> > mass    = differentialMass(iBinning, 1, true);
  
  double massShift = 0;
  
  for (unsigned int i = 0; i < entries[0].size(); ++i) {
    std::cout << "Data -- entries: " << entries[0][i] << " dmt: " << mass[0][i] << std::endl;
    std::cout << "MGZ2 -- entries: " << entries[1][i] << " dmt: " << mass[1][i] << std::endl;
    
    massShift += (entries[0][i]/entries[1][i] - 1.) * mass[1][i] * entries[1][i];
  }
  
  std::cout << "Mass shift for reweighting to data: " << massShift << std::endl;
}
*/
