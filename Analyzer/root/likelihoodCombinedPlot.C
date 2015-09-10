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


void likelihoodCombinedPlot() {
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  
  int lineColor2 = kAzure+1;
  int lineColorh = kGreen+1;
  int lineColor1 = kRed+1;
  
  int lineWidthThick = 4;
  int lineWidthThin  = 2;
  
  TFile f2("../plot/Ideogram/Run2012_JEC_Winter14_V8_mass_vs_JES.root","read");
  TCanvas* c2 = (TCanvas*) f2.Get("canv");
  TGraph* c2_gr2 = (TGraph*) c2->GetPrimitive("gr2");
          c2_gr2->SetLineColor(lineColor2);
          c2_gr2->SetLineWidth(lineWidthThin);
          //c2_gr2->SetLineStyle(7);
  TGraph* c2_gr1 = (TGraph*) c2->GetPrimitive("gr1");
          c2_gr1->SetLineColor(lineColor2);
          c2_gr1->SetLineWidth(lineWidthThick);
          //c2_gr1->SetLineStyle(7);
  TGraph* c2_gr0 = (TGraph*) c2->GetPrimitive("gr0");
          c2_gr0->SetLineColor(lineColor2);
          c2_gr0->SetLineWidth(lineWidthThick);
          c2_gr0->SetMarkerColor(lineColor2);
          c2_gr0->SetMarkerStyle(34);
  
  double mt2d[]         = {172.141};
  double mt2d_stat[]    = {  0.193914};
  double mt2d_stat2[]   = {2*0.193914};
  double mt2d_jsf[]     = {1.00485};
  double mt2d_jsferr[]  = {  0.00180311};
  double mt2d_jsferr2[] = {2*0.00180311};
  TGraphErrors* c2_c2 = new TGraphErrors(1, mt2d, mt2d_jsf, mt2d_stat2, mt2d_jsferr2);
                c2_c2->SetLineColor(kBlack);
                c2_c2->SetLineWidth(lineWidthThin);
                c2_c2->SetMarkerColor(lineColor2);
                c2_c2->SetMarkerSize(2);
                c2_c2->SetMarkerStyle(34);
  TGraphErrors* c2_c1 = new TGraphErrors(1, mt2d, mt2d_jsf, mt2d_stat, mt2d_jsferr);
                c2_c1->SetLineColor(kBlack);
                c2_c1->SetLineWidth(lineWidthThick);
                c2_c1->SetMarkerColor(lineColor2);
                c2_c1->SetMarkerSize(2);
                c2_c1->SetMarkerStyle(34);
  
  TFile fh("../plot/Ideogram/Run2012_JEC_Winter14_V8_mass_vs_JES_hyb.root","read");
  TCanvas* ch = (TCanvas*) fh.Get("canv");
  TGraph* ch_gr2 = (TGraph*) ch->GetPrimitive("gr2");
          ch_gr2->SetLineColor(lineColorh);
          ch_gr2->SetLineWidth(lineWidthThin);
          //ch_gr2->SetLineStyle(7);
  TGraph* ch_gr1 = (TGraph*) ch->GetPrimitive("gr1");
          ch_gr1->SetLineColor(lineColorh);
          ch_gr1->SetLineWidth(lineWidthThick);
          //ch_gr1->SetLineStyle(7);
  TGraph* ch_gr0 = (TGraph*) ch->GetPrimitive("gr0");
          ch_gr0->SetLineColor(lineColorh);
          ch_gr0->SetLineWidth(lineWidthThick);
          ch_gr0->SetMarkerColor(lineColorh);
          ch_gr0->SetMarkerStyle(29);
  
  double mth[]         = {172.349};
  double mth_stat[]    = {  0.160533};
  double mth_stat2[]   = {2*0.160533};
  double mth_jsf[]     = {1.00243};
  double mth_jsferr[]  = {  0.00127309};
  double mth_jsferr2[] = {2*0.00127309};
  TGraphErrors* ch_c2 = new TGraphErrors(1, mth, mth_jsf, mth_stat2, mth_jsferr2);
                ch_c2->SetLineColor(kBlack);
                ch_c2->SetLineWidth(lineWidthThin);
                ch_c2->SetMarkerColor(lineColorh);
                ch_c2->SetMarkerSize(2);
                ch_c2->SetMarkerStyle(29);
  TGraphErrors* ch_c1 = new TGraphErrors(1, mth, mth_jsf, mth_stat, mth_jsferr);
                ch_c1->SetLineColor(kBlack);
                ch_c1->SetLineWidth(lineWidthThick);
                ch_c1->SetMarkerColor(lineColorh);
                ch_c1->SetMarkerSize(2);
                ch_c1->SetMarkerStyle(29);
  
  double mt1d[]         = {172.558};
  double mt1d_stat[]    = {  0.117746};
  double mt1d_stat2[]   = {2*0.117746};
  double mt1d_jsf[]     = {1.};
  double mt1d_jsferr[]  = {0.00001};
  TGraphErrors* c1_c2 = new TGraphErrors(1, mt1d, mt1d_jsf, mt1d_stat2, mt1d_jsferr);
                c1_c2->SetLineColor(lineColor1);
                c1_c2->SetLineWidth(lineWidthThin);
                c1_c2->SetMarkerColor(lineColor1);
                c1_c2->SetMarkerSize(2);
                ch_gr2->SetMarkerStyle(20);
  TGraphErrors* c1_c1 = new TGraphErrors(1, mt1d, mt1d_jsf, mt1d_stat, mt1d_jsferr);
                c1_c1->SetLineColor(lineColor1);
                c1_c1->SetLineWidth(lineWidthThick);
                c1_c1->SetMarkerColor(lineColor1);
                c1_c1->SetMarkerSize(2);
                ch_gr1->SetMarkerStyle(20);
  
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle(";m_{t} [GeV];JSF");
  mg->Add(c2_gr2);
  mg->Add(ch_gr2);
  mg->Add( c1_c2);
  
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
  
  mg->Draw("A");
  
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
  
  //Draw8LeptonJets();
  CMS_lumi(canvas, 21, 0.);
  
  canvas->Print("likelihoodCombinedPlot.eps");
}
