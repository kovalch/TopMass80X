#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLatex.h"
#include "TBox.h"

#include "tdrstyle.C"

enum styles             { kStat,  kOverall, kTev};
int color_      [ 4 ] = { kRed+1, kBlack  , kGray};

TH1D* hError;
TH1D* hPull;
TTree* tree;

void DrawLabel(TString text, const double x1, const double y1, const double x2)
{
  // function to directly draw a label into the active canvas
  double y2 = y1 + 0.05;
  double yOffset = 0.02;
  TPaveLabel *label = new TPaveLabel(x1, y1+yOffset, x2, y2+yOffset, text, "br NDC");
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  label->SetTextSize(0.75);
  label->SetTextAlign(12);
  label->SetTextColor(kBlack);
  label->Draw("same");
}

void drawcutline(double cutval, double maximum)
{
  TLine *cut = new TLine();
  cut->SetLineWidth(2);
  cut->SetLineStyle(7);
  cut->SetLineColor(1);
  cut->DrawLine(cutval, 0.,cutval, maximum);
}

void drawArrow(double cutval, double maximum)
{
  TArrow *arrow = new TArrow();
  arrow->SetLineWidth(3);
  arrow->SetLineColor(kBlack);
  arrow->DrawArrow(cutval, 0., cutval, maximum, 0.05, "<");
}

void ensembleMassesJSF()
{
  setTDRStyle();
  tdrStyle->SetNdivisions(505, "X");
  tdrStyle->SetPadLeftMargin(0.2);
  tdrStyle->SetTitleYOffset(1.8);
  
  //// Get histos
  
  TChain* tree = new TChain("tree");
  tree->Add("/nfs/dust/cms/user/mseidel/pseudoexperiments/topmass_paper_uncalibrated/lepton/Summer12_TTJetsMS1725_0.98/job_*_ensemble.root");
  tree->Add("/nfs/dust/cms/user/mseidel/pseudoexperiments/topmass_paper_uncalibrated/lepton/Summer12_TTJetsMS1725_1.00/job_*_ensemble.root");
  tree->Add("/nfs/dust/cms/user/mseidel/pseudoexperiments/topmass_paper_uncalibrated/lepton/Summer12_TTJetsMS1725_1.02/job_*_ensemble.root");
  
  hError = new TH1D("hError", "hError", 400, 0.95, 1.05);
  hPull  = new TH1D("hPull", "hPull", 100, -5, 5);
  
  double mass;
  tree->SetBranchAddress("JES_mTop_JES", &mass);
  
  double genMass;
  tree->SetBranchAddress("genMass", &genMass);
  
  double genJES;
  tree->SetBranchAddress("genJES", &genJES);
  
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (mass>0) hError->Fill(mass);
  }
  
  hError->GetYaxis()->SetTitle("Pseudo-experiments / 0.00025");
  hError->GetYaxis()->SetTitleSize(0.05);
  hError->GetXaxis()->SetTitle("JSF_{extr}");
  hPull->GetYaxis()->SetTitle("Pseudo-experiments");
  hPull->GetXaxis()->SetTitle("Mass pull");
  
  //// Do something :)
  
  TCanvas* cError = new TCanvas("cMassesJSF", "cMassesJSF", 600, 600);
  
  //hError->Rebin(4);
  //hError->GetXaxis()->SetRangeUser(0.19, 0.20);
  //hError->GetYaxis()->SetRangeUser(0, 1150);
  hError->SetFillColor(kRed+1);
  hError->Draw();
  
  Draw8LeptonJets();
  
  /*
  TCanvas* cPull = new TCanvas("cPull", "cPull", 600, 600);
  
  //hPull->Rebin(4);
  //hPull->GetXaxis()->SetRangeUser(1., 1.7);
  //hPull->GetYaxis()->SetRangeUser(0, 400);
  hPull->Fit("gaus");
  //*/
  
}
