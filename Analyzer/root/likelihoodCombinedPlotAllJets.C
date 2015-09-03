#include <string>
#include <vector>
#include <iostream>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

#include "tdrstyle_new.C"
#include "CMS_lumi.C"


void likelihoodCombinedPlotAllJets() {
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadTopMargin(0.08);
  //gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  
  int lineColor2 = kAzure+1;
  int lineColorh = kGreen+1;
  int lineColor1 = kRed+1;
  
  int lineWidthThick = 4;
  int lineWidthThin  = 2;
  
  TFile f2("../plot/Ideogram/Run2012_alljets_mass_vs_JES_old.root","read");
  TCanvas* c2 = (TCanvas*) f2.Get("canv");
  TGraph* c2_gr2 = (TGraph*) c2->GetPrimitive("gr2");
          c2_gr2->SetLineColor(lineColor2);
          c2_gr2->SetLineWidth(lineWidthThin);
  TGraph* c2_gr1 = (TGraph*) c2->GetPrimitive("gr1");
          c2_gr1->SetLineColor(lineColor2);
          c2_gr1->SetLineWidth(lineWidthThick);
  TGraph* c2_gr0 = (TGraph*) c2->GetPrimitive("gr0");
          c2_gr0->SetLineColor(lineColor2);
          c2_gr0->SetLineWidth(lineWidthThick);
          c2_gr0->SetMarkerColor(lineColor2);
          c2_gr0->SetMarkerStyle(34);
          
  double mt2d[]         = {171.637};
  double mt2d_stat[]    = {  0.327759};
  double mt2d_stat2[]   = {2*0.327759};
  double mt2d_jsf[]     = {1.01074};
  double mt2d_jsferr[]  = {  0.00302304};
  double mt2d_jsferr2[] = {2*0.00302304};
  TGraphErrors* c2_c2 = new TGraphErrors(1, mt2d, mt2d_jsf, mt2d_stat2, mt2d_jsferr2);
                c2_c2->SetLineColor(lineColor2);
                c2_c2->SetLineWidth(lineWidthThin);
                c2_c2->SetMarkerColor(lineColor2);
                c2_c2->SetMarkerSize(2);
                c2_c2->SetMarkerStyle(34);
  TGraphErrors* c2_c1 = new TGraphErrors(1, mt2d, mt2d_jsf, mt2d_stat, mt2d_jsferr);
                c2_c1->SetLineColor(lineColor2);
                c2_c1->SetLineWidth(lineWidthThick);
                c2_c1->SetMarkerColor(lineColor2);
                c2_c1->SetMarkerSize(2);
                c2_c1->SetMarkerStyle(34);
  
  TFile fh("../plot/Ideogram/Run2012_alljets_mass_vs_JESc_old.root","read");
  TCanvas* ch = (TCanvas*) fh.Get("canv");
  TGraph* ch_gr2 = (TGraph*) ch->GetPrimitive("gr2");
          ch_gr2->SetLineColor(lineColorh);
          ch_gr2->SetLineWidth(lineWidthThin);
  TGraph* ch_gr1 = (TGraph*) ch->GetPrimitive("gr1");
          ch_gr1->SetLineColor(lineColorh);
          ch_gr1->SetLineWidth(lineWidthThick);
  TGraph* ch_gr0 = (TGraph*) ch->GetPrimitive("gr0");
          ch_gr0->SetLineColor(lineColorh);
          ch_gr0->SetLineWidth(lineWidthThick);
          ch_gr0->SetMarkerColor(lineColorh);
          ch_gr0->SetMarkerStyle(29);
  
  double mth[]         = {172.318};
  double mth_stat[]    = {  0.252721};
  double mth_stat2[]   = {2*0.252721};
  double mth_jsf[]     = {1.0019};
  double mth_jsferr[]  = {  0.00127085};
  double mth_jsferr2[] = {2*0.00127085};
  TGraphErrors* ch_c2 = new TGraphErrors(1, mth, mth_jsf, mth_stat2, mth_jsferr2);
                ch_c2->SetLineColor(lineColorh);
                ch_c2->SetLineWidth(lineWidthThin);
                ch_c2->SetMarkerColor(lineColorh);
                ch_c2->SetMarkerSize(2);
                ch_c2->SetMarkerStyle(29);
  TGraphErrors* ch_c1 = new TGraphErrors(1, mth, mth_jsf, mth_stat, mth_jsferr);
                ch_c1->SetLineColor(lineColorh);
                ch_c1->SetLineWidth(lineWidthThick);
                ch_c1->SetMarkerColor(lineColorh);
                ch_c1->SetMarkerSize(2);
                ch_c1->SetMarkerStyle(29);
  
  double mt1d[]         = {172.463};
  double mt1d_stat[]    = {  0.233541};
  double mt1d_stat2[]   = {2*0.233541};
  double mt1d_jsf[]     = {1.};
  double mt1d_jsferr[]  = {0.00001};
  TGraphErrors* c1_c2 = new TGraphErrors(1, mt1d, mt1d_jsf, mt1d_stat2, mt1d_jsferr);
                c1_c2->SetLineColor(lineColor1);
                c1_c2->SetLineWidth(lineWidthThin);
                c1_c2->SetMarkerColor(lineColor1);
                c1_c2->SetMarkerSize(2);
                c1_c2->SetMarkerStyle(20);
  TGraphErrors* c1_c1 = new TGraphErrors(1, mt1d, mt1d_jsf, mt1d_stat, mt1d_jsferr);
                c1_c1->SetLineColor(lineColor1);
                c1_c1->SetLineWidth(lineWidthThick);
                c1_c1->SetMarkerColor(lineColor1);
                c1_c1->SetMarkerSize(2);
                c1_c1->SetMarkerStyle(20);
  
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle(";m_{t} [GeV];JSF");
  mg->Add(c2_gr2);
  mg->Add(ch_gr2);
  mg->Add(c1_c2);
  
  TLegend *leg = new TLegend(0.75, 0.73, 0.93, 0.90);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(1001);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(c2_gr0, "2D", "P");
  leg->AddEntry(ch_gr0, "hybrid", "P");
  leg->AddEntry(c1_c1, "1D", "P");
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 500, 500);
  canvas->cd();
  
  mg->Draw("AC");
  
  c2_gr2->Draw("C,SAME");
  c2_gr1->Draw("C,SAME");
  c2_gr0->Draw("P,SAME");
  //c2_c2 ->Draw("P,SAME");
  //c2_c1 ->Draw("P,SAME");
  
  ch_gr2->Draw("C,SAME");
  ch_gr1->Draw("C,SAME");
  ch_gr0->Draw("P,SAME");
  //ch_c2 ->Draw("P,SAME");
  //ch_c1 ->Draw("P,SAME");
  
  c1_c2->Draw("P,SAME");
  c1_c1->Draw("P,SAME");
  
  leg->Draw();
  
  //DrawCMS8AllJets();
  CMS_lumi(canvas, 22, 0.);
  
  canvas->Print("likelihoodCombinedPlotAllJets.eps");
}
