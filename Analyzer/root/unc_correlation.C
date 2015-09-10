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

#include "tdrstyle_new.C"
#include "CMS_lumi.C"

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


void unc_correlation()
{
  TStyle* tdrStyle = setTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadTopMargin(0.08);
  tdrStyle->SetTitleYOffset(1.4);
  tdrStyle->SetNdivisions(505, "XYZ");
  tdrStyle->SetOptStat(0);
  tdrStyle->SetOptFit(0);
  
  TCanvas* corr1 = new TCanvas("corr1", "corr1", 600, 600);
  corr1->cd();
  
  TMultiGraph *mg = new TMultiGraph();
  
  double dummyv[] = {0.}; double dummyu[] = {0.5}; 
  
  TGraphErrors* gDummy = new TGraphErrors(1, dummyv, dummyv, dummyu, dummyu);
                gDummy->SetMarkerColor(kBlack);
                gDummy->SetLineColor  (kBlack);
                gDummy->SetMarkerStyle(1);
                gDummy->SetLineStyle(7);
  
  mg->Add(gDummy);
  
  TString uncs[] = {"bJES", "bres", "dintJES", "dMPFJES", "dunJES", "iJES", "LES", "MCGen", "PDF", "Q", "JPS", "JER", "BTag", "MET", "UE", "BKGDT", "BKGMC", "FitCal", "PU", "CR", "Trig", "Toppt"};
  
  /*
  mg->SetTitle(";#deltaJSF^{l+jets};#deltaJSF^{all-jets}");
  double lj2v[] = {0.001, 0.001, 0.000, 0.003, 0.004, 0.000, 0.000, 0.001, 0.001, 0.004, 0.002, 0.002, 0.000, 0.000, 0.002,  0.00, 0.000, 0.001, 0.003, 0.002, 0.00, 0.003};
  double lj2u[] = {0.00, 0.00,  0.00,  0.00,  0.00, 0.002, 0.00,  0.001, 0.00,  0.001,  0.001,  0.00, 0.00, 0.00,  0.001,  0.00, 0.00, 0.00,  0.00,  0.001, 0.00,  0.00};
  double aj1v[] = {0.001, 0.001, 0.000, 0.000, 0.001, 0.00, 0.00, 0.002, 0.000, 0.005, 0.002, 0.001, 0.000, 0.00, 0.002, 0.007, 0.00, 0.001, 0.002, 0.003, 0.000, 0.001};
  double aj1u[] = {0.00, 0.00,  0.00,  0.00,  0.00, 0.003, 0.00,  0.002, 0.00,  0.001,  0.001,  0.00, 0.00, 0.00,  0.002,  0.00, 0.00, 0.00,  0.00,  0.002, 0.00,  0.00};
  TString filename("unc_corr_lj2_aj2_jsf.pdf");
  //*/
  
  /* paper!
  mg->SetTitle(";#deltam_{t (2D)}^{l+jets} [GeV];#deltam_{t (1D)}^{all-jets} [GeV]");
  double lj2v[] = {-0.40, -0.17, +0.00, -0.01, +0.09, 0.00, 0.01, -0.07, 0.09, +0.17, +0.11, -0.11, 0.06, 0.04, +0.15,  0.00, 0.05, 0.04, -0.13, +0.11, 0.00, +0.16};
  double lj2u[] = {0.00, 0.00,  0.00,  0.00,  0.00, 0.15, 0.00,  0.11, 0.00,  0.08,  0.09,  0.00, 0.00, 0.00,  0.15,  0.00, 0.00, 0.00,  0.00,  0.13, 0.00,  0.00};
  double aj1v[] = {-0.30, -0.12, -0.02, +0.22, -0.19, 0.00, 0.00, -0.18, 0.01, -0.19, +0.12, +0.03, 0.01, 0.00, +0.13, -0.14, 0.00, 0.06, +0.11, +0.14, 0.01, +0.08};
  double aj1u[] = {0.00, 0.00,  0.00,  0.00,  0.00, 0.00, 0.00,  0.14, 0.00,  0.11,  0.11,  0.00, 0.00, 0.00,  0.18,  0.00, 0.00, 0.00,  0.00,  0.16, 0.00,  0.00};
  TString filename("unc_corr_lj2_aj1.pdf");
  //*/
  
  /*
  mg->SetTitle(";#deltam_{t (2D)}^{l+jets} [GeV];#deltam_{t (1D)}^{l+jets} [GeV]");
  double lj2v[] = {0.40, 0.17, +0.00, -0.01, +0.09, 0.00, 0.01, -0.07, 0.09, +0.17, +0.11, -0.11, 0.06, 0.04, +0.15,  0.00, 0.05, 0.04, -0.13, +0.11, 0.00, +0.16};
  double lj2u[] = {0.00, 0.00,  0.00,  0.00,  0.00, 0.15, 0.00,  0.11, 0.00,  0.08,  0.09,  0.00, 0.00, 0.00,  0.15,  0.00, 0.00, 0.00,  0.00,  0.13, 0.00,  0.00};
  double aj1v[] = {0.30, 0.17, -0.02, +0.24, -0.26, 0.00, 0.01, -0.16, 0.06, -0.24, -0.07, +0.05, 0.04, 0.03, +0.07,  0.00, 0.01, 0.04, +0.12, -0.09, 0.00, -0.11};
  double aj1u[] = {0.00, 0.00,  0.00,  0.00,  0.00, 0.00, 0.00,  0.07, 0.00,  0.06,  0.06,  0.00, 0.00, 0.00,  0.09,  0.00, 0.00, 0.00,  0.00,  0.08, 0.00,  0.00};
  TString filename("unc_corr_lj2_lj1.pdf");
  //*/
  
  //* paper!
  mg->SetTitle(";#deltam_{t (hybrid)}^{l+jets} [GeV];#deltam_{t (1D)}^{all-jets} [GeV]");
  double lj2v[] = {-0.35, -0.16, +0.01, +0.12, -0.10, 0.00, 0.01, -0.12, 0.04, -0.09, +0.02, -0.03, 0.06, 0.04, +0.08, 00.00, 0.03, 0.04,  0.08, +0.01, 0.00, +0.02};
  double lj2u[] = {0.00, 0.00,  0.00,  0.00,  0.00, 0.00, 0.00,  0.08, 0.00,  0.08,  0.09,  0.00, 0.00, 0.00,  0.10,  0.00, 0.00, 0.00,  0.00,  0.09, 0.00,  0.00};
  double aj1v[] = {-0.30, -0.12, -0.02, +0.22, -0.19, 0.00, 0.00, -0.18, 0.01, -0.19, +0.12, +0.03, 0.01, 0.00, +0.13, -0.14, 0.00, 0.06, +0.11, +0.14, 0.01, +0.08};
  double aj1u[] = {0.00, 0.00,  0.00,  0.00,  0.00, 0.00, 0.00,  0.14, 0.00,  0.11,  0.11,  0.00, 0.00, 0.00,  0.18,  0.00, 0.00, 0.00,  0.00,  0.16, 0.00,  0.00};
  TString filename("unc_corr_ljh_aj1.pdf");
  //*/
  
  /*
  mg->SetTitle(";#deltam_{t (hybrid)}^{l+jets} [GeV];#deltam_{t (1D)}^{l+jets} [GeV]");
  double lj2v[] = {0.35, 0.16, +0.01, +0.12, -0.10, 0.00, 0.01, -0.12, 0.04, -0.09, +0.02, -0.03, 0.06, 0.04, +0.08, 00.00, 0.03, 0.04,  0.08, +0.01, 0.00, +0.02};
  double lj2u[] = {0.00, 0.00,  0.00,  0.00,  0.00, 0.00, 0.00,  0.08, 0.00,  0.08,  0.09,  0.00, 0.00, 0.00,  0.10,  0.00, 0.00, 0.00,  0.00,  0.09, 0.00,  0.00};
  double aj1v[] = {0.30, 0.17, -0.02, +0.24, -0.26, 0.00, 0.01, -0.16, 0.06, -0.24, -0.07, +0.05, 0.04, 0.03, +0.07,  0.00, 0.01, 0.04, +0.12, -0.09, 0.00, -0.11};
  double aj1u[] = {0.00, 0.00,  0.00,  0.00,  0.00, 0.00, 0.00,  0.07, 0.00,  0.06,  0.06,  0.00, 0.00, 0.00,  0.09,  0.00, 0.00, 0.00,  0.00,  0.08, 0.00,  0.00};
  TString filename("unc_corr_ljh_lj1.pdf");
  //*/
  
  std::vector<double> lj2v_anti, lj2u_anti, aj1v_anti, aj1u_anti, lj2v_corr, lj2u_corr, aj1v_corr, aj1u_corr, lj2v_ambi, lj2u_ambi, aj1v_ambi, aj1u_ambi;
  
  for (int i = 0; i < 22; ++i) {
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
  
  TLegend *leg1 = new TLegend(0.2, 0.755, 0.6, 0.90);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->AddEntry( gCorr, "#rho > 0", "P" );
  leg1->AddEntry( gAmbi, "ambiguous", "P" );
  leg1->AddEntry( gAnti, "#rho < 0", "P" );
  
  leg1->Draw();
  
  CMS_lumi(corr1, 12, 0.);
  
  corr1->Print(filename);
}
