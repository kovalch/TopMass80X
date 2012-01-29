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

TH3F* h3MassError_1725;
TH3F* h3MassPull_1725;
TH1D* hError;
TH1D* hPull;

void DrawLabel(TString text, const double x1, const double y1, const double x2, Color_t color = kBlack)
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
  arrow->SetLineWidth(2);
  arrow->SetLineColor(kRed+1);
  arrow->DrawArrow(cutval, 0., cutval, maximum, 0.05, "<");
}

void ensemble()
{
  setTDRStyle();
  tdrStyle->SetNdivisions(505, "X");
  
  //// Get histos
  
  TFile* fEnsemble = new TFile("/scratch/hh/current/cms/user/mseidel/topmass_120120_corrected/ensemble.root");
  h3MassError_1725 = (TH3F*) fEnsemble->Get("h3MassError_1725");
  h3MassPull_1725 = (TH3F*) fEnsemble->Get("h3MassPull_1725");
  
  hError = new TH1D("hError", "hError", 200, 0, 2);
  hPull  = new TH1D("hPull", "hPull", 100, -5, 5);
  
  h3MassError_1725->ProjectionZ("hError");
  h3MassPull_1725->ProjectionZ("hPull");
  
  hError->GetYaxis()->SetTitle("Number of pseudo-experiments");
  hError->GetXaxis()->SetTitle("#sigma(m_{t}) [GeV]");
  hPull->GetYaxis()->SetTitle("Number of pseudo-experiments");
  hPull->GetXaxis()->SetTitle("Mass pull");
  
  //// Do something :)
  
  TCanvas* cError = new TCanvas("cError", "cError", 600, 600);
  
  //hError->Rebin(4);
  hError->GetXaxis()->SetRangeUser(0.6, 0.8);
  hError->GetYaxis()->SetRangeUser(0, 1500);
  hError->Draw();
  drawArrow(0.70, 1300);
  DrawLabel("4.7 fb^{-1} collision data", 0.36, 0.85, 0.46, kRed+1);
  
  /*
  TCanvas* cPull = new TCanvas("cPull", "cPull", 600, 600);
  
  //hPull->Rebin(4);
  //hPull->GetXaxis()->SetRangeUser(1., 1.7);
  //hPull->GetYaxis()->SetRangeUser(0, 400);
  hPull->Fit("gaus");
  //*/
  
}
