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
#include "TPad.h"
#include "TPaveLabel.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TLegend.h"

#include "tdrstyle.C"

enum styles             { kStat,  kJEC,   kTev};
int color_      [ 4 ] = { kRed+1, kBlack, kGray};


void DrawLabel(TString text, double size, const double x1, const double y1, const double x2, Color_t color = kBlack)
{
  // function to directly draw a label into the active canvas
  double y2 = y1 + 0.05;
  TPaveLabel *label = new TPaveLabel(x1, y1, x2, y2, text, "NDC");
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  label->SetTextSize(size);
  label->SetTextAlign(12);
  label->SetTextColor(color);
  label->Draw("same");
}


void unc_correlation_jsf()
{
  TStyle* tdrStyle = setTDRStyle();
  tdrStyle->SetPadLeftMargin(0.2);
  tdrStyle->SetTitleYOffset(1.4);
  tdrStyle->SetNdivisions(505, "XYZ");
  tdrStyle->SetOptStat(0);
  tdrStyle->SetOptFit(0);
  
  TCanvas* corr1 = new TCanvas("corr1", "corr1", 600, 600);
  corr1->cd();
  
  TMultiGraph *mg = new TMultiGraph();
  
  double dummyv[] = {0.}; double dummyu[] = {0.01}; 
  
  TGraphErrors* gDummy = new TGraphErrors(1, dummyv, dummyv, dummyu, dummyu);
                gDummy->SetMarkerColor(kBlack);
                gDummy->SetLineColor  (kBlack);
                gDummy->SetMarkerStyle(1);
                gDummy->SetLineStyle(7);
  
  mg->Add(gDummy);
  
  TString uncs[] = {"bJES", "bres", "iJES", "LES", "MCGen", "PDF", "Q", "JPS", "JER", "BTag", "MET", "UE", "BKGDT", "BKGMC", "FitCal", "PU", "CR", "Trig", "Toppt"};
  
  //*
  mg->SetTitle(";#deltaJSF^{l+jets};#deltaJSF^{all-jets}");
  double lj2v[] = {+0.001, -0.001, 0.000, 0.000, -0.001, 0.001, -0.004, -0.002, +0.002, 0.000, 0.000, +0.002,  0.00, 0.000, 0.001, +0.002, -0.002, 0.00, -0.003};
  double lj2u[] = {0.00, 0.00,  0.002, 0.00,  0.001, 0.00,  0.001,  0.001,  0.00, 0.00, 0.00,  0.001,  0.00, 0.00, 0.00,  0.00,  0.001, 0.00,  0.00};
  double aj1v[] = {+0.001, -0.001, 0.00, 0.00, +0.002, 0.000, -0.005, -0.002, +0.001, 0.000, 0.00, +0.002, 0.007, 0.00, 0.001, +0.001, -0.003, 0.000, +0.001};
  double aj1u[] = {0.00, 0.00,  0.003, 0.00,  0.002, 0.00,  0.001,  0.001,  0.00, 0.00, 0.00,  0.002,  0.00, 0.00, 0.00,  0.00,  0.002, 0.00,  0.00};
  TString filename("unc_corr_lj2_aj2_jsf.pdf");
  //*/
  
  std::vector<double> lj2v_anti, lj2u_anti, aj1v_anti, aj1u_anti, lj2v_corr, lj2u_corr, aj1v_corr, aj1u_corr, lj2v_ambi, lj2u_ambi, aj1v_ambi, aj1u_ambi;
  
  for (int i = 0; i < 19; ++i) {
    std::cout << uncs[i] << "\t";
    std::cout << lj2v[i] << "\t";
    std::cout << aj1v[i] << "\t";
    // uncertainty larger than shift -> ambiguous
    if (fabs(lj2v[i])-lj2u[i]<0. || fabs(aj1v[i])-aj1u[i]<0.) {
      std::cout << "ambiguous (+/-" << lj2u[i] << "," << aj1u[i] << ")" << std::endl;
      lj2v_ambi.push_back(lj2v[i]);
      lj2u_ambi.push_back(lj2u[i]);
      aj1v_ambi.push_back(aj1v[i]);
      aj1u_ambi.push_back(aj1u[i]);
    }
    // correlated
    else if (lj2v[i]*aj1v[i]>=0.) {
      std::cout << "correlated" << std::endl;
      lj2v_corr.push_back(lj2v[i]);
      lj2u_corr.push_back(lj2u[i]);
      aj1v_corr.push_back(aj1v[i]);
      aj1u_corr.push_back(aj1u[i]);
    }
    // anti-correlated
    else if (lj2v[i]*aj1v[i]<0.) {
      std::cout << "anti-correlated" << std::endl;
      lj2v_anti.push_back(lj2v[i]);
      lj2u_anti.push_back(lj2u[i]);
      aj1v_anti.push_back(aj1v[i]);
      aj1u_anti.push_back(aj1u[i]);
    }
  }
  std::cout << std::endl;
  
  TGraphErrors* gAmbi = new TGraphErrors(lj2v_ambi.size(), &(lj2v_ambi[0]), &(aj1v_ambi[0]), &(lj2u_ambi[0]), &(aj1u_ambi[0]));
                gAmbi->SetMarkerColor(kOrange+1);
                gAmbi->SetLineColor  (kOrange+1);
                gAmbi->SetMarkerStyle(21);
  
  TGraphErrors* gCorr = new TGraphErrors(lj2v_corr.size(), &(lj2v_corr[0]), &(aj1v_corr[0]), &(lj2u_corr[0]), &(aj1u_corr[0]));
                gCorr->SetMarkerColor(kRed+1);
                gCorr->SetLineColor  (kRed+1);
                gCorr->SetMarkerStyle(20);
  
  TGraphErrors* gAnti = new TGraphErrors(lj2v_anti.size(), &(lj2v_anti[0]), &(aj1v_anti[0]), &(lj2u_anti[0]), &(aj1u_anti[0]));
                gAnti->SetMarkerColor(kBlue+1);
                gAnti->SetLineColor  (kBlue+1);
                gAnti->SetMarkerStyle(24);
  
  mg->Add(gAmbi);
  mg->Add(gCorr);
  mg->Add(gAnti);
  
  mg->Draw("AP");
  
  //gJEC->Fit("abslin");
  
  //mg->SetMinimum(0.0);
  //mg->SetMaximum(1.0);
  
  TLegend *leg1 = new TLegend(0.2, 0.775, 0.6, 0.92);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->AddEntry( gCorr, "#rho = +100%", "P" );
  leg1->AddEntry( gAmbi, "ambiguous", "P" );
  leg1->AddEntry( gAnti, "#rho = -100%", "P" );
  
  leg1->Draw();
  
  //DrawCMSSim();
  
  corr1->Print(filename);
}
