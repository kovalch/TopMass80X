#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TBox.h"
#include "TArrow.h"
#include "TLine.h"

#include "tdrstyle.C"

//enum styles       { kStat,  kOverall, kTev};
//int color_[ 3 ] = { kRed+1, kBlack  , kGray};

//double genMasses[] = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5, 184.5};
//double genJESes [] = {0.96, 1.00, 1.04};

double genMasses[] = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5, 184.5};
double genJESes [] = {0.96, 0.98, 1.00, 1.02, 1.04};

std::vector<std::vector<TH1D*> > hErrors;
std::vector<std::vector<TH1D*> > hErrorsConstJES;
std::vector<std::vector<TH1D*> > hErrorsFSig;
std::vector<std::vector<TH1D*> > hPulls;
std::vector<std::vector<TH1D*> > hPullsConstJES;
std::vector<std::vector<TH1D*> > hPullsFSig;
std::vector<std::vector<TH1D*> > hMasses;
std::vector<std::vector<TH1D*> > hMassesConstJES;
std::vector<std::vector<TH1D*> > hMassesFSig;
TTree* tree;

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

void ensemblePlots()
{
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  //// Get histos
  
  TString sFile("/scratch/hh/dust/naf/cms/user/kirschen/BRegression_PE_NoCalibrationApplied/muon_BReg/");
  //  TString sFile("/scratch/hh/current/cms/user/eschliec/TopMass/");
  //TString sFile("/scratch/hh/dust/naf/cms/user/eschliec/TopMass/");
  //sFile += "ensemble_F11_ROOFIT_Recalibrated_NEWER.root";
  //sFile += "ensemble_F11_ROOFIT_Calibrated_BUGFIXED_5JES_NEW.root";
  //  sFile += "ensemble_F11_TemplateRooFit_Uncalibrated.root";
  //sFile += "ensemble_F11_TemplateRooFit_Calibrated.root";
  sFile += "ensemble.root";
  TFile* fEnsemble = new TFile(sFile);
  tree = (TTree*) fEnsemble->Get("tree");
  
  TH1D * hError = new TH1D("hError", "hError", 200, 0, 2);
  hError->GetYaxis()->SetTitle("Number of pseudo-experiments / 10 MeV");
  hError->GetYaxis()->SetTitleSize(0.05);
  hError->GetYaxis()->SetTitleOffset(1.5);
  hError->GetXaxis()->SetTitle("#sigma(m_{t}) [GeV]");  
  //hErrorConstJES = (TH1D*)hError->Clone("hErrorConstJES");

  TH1D * hPull = new TH1D("hPull", "hPull", 200, -10, 10);
  hPull->GetYaxis()->SetTitle("Number of pseudo-experiments");
  hPull->GetYaxis()->SetTitleSize(0.05);
  hPull->GetYaxis()->SetTitleOffset(1.5);
  hPull->GetXaxis()->SetTitle("pull m_{t}");

  TH1D * hMass = new TH1D("hMass", "hMass", 300, 140, 200);
  hMass->GetYaxis()->SetTitle("Number of pseudo-experiments / 1 GeV");
  hMass->GetYaxis()->SetTitleSize(0.05);
  hMass->GetYaxis()->SetTitleOffset(1.5);
  hMass->GetXaxis()->SetTitle("m_{t} [GeV]");

  for(unsigned int i = 0; i < sizeof(genMasses)/sizeof(double); ++i){
    std::vector<TH1D*> helper_hErrors        ;
    std::vector<TH1D*> helper_hErrorsConstJES;
    std::vector<TH1D*> helper_hErrorsFSig    ;
    std::vector<TH1D*> helper_hPulls         ;
    std::vector<TH1D*> helper_hPullsConstJES ;
    std::vector<TH1D*> helper_hPullsFSig     ;
    std::vector<TH1D*> helper_hMasses        ;
    std::vector<TH1D*> helper_hMassesConstJES;
    std::vector<TH1D*> helper_hMassesFSig    ;
    for(unsigned int j = 0; j < sizeof(genJESes)/sizeof(double); ++j){
      helper_hErrors        .push_back((TH1D*)hError->Clone());
      helper_hErrorsConstJES.push_back((TH1D*)hError->Clone());
      helper_hErrorsFSig    .push_back((TH1D*)hError->Clone());
      helper_hPulls         .push_back((TH1D*)hPull ->Clone());
      helper_hPullsConstJES .push_back((TH1D*)hPull ->Clone());
      helper_hPullsFSig     .push_back((TH1D*)hPull ->Clone());
      helper_hMasses        .push_back((TH1D*)hMass ->Clone());
      helper_hMassesConstJES.push_back((TH1D*)hMass ->Clone());  
      helper_hMassesFSig    .push_back((TH1D*)hMass ->Clone());
    }    
    hErrors        .push_back(helper_hErrors        );
    hErrorsConstJES.push_back(helper_hErrorsConstJES);
    hErrorsFSig    .push_back(helper_hErrorsFSig    );
    hPulls         .push_back(helper_hPulls         );
    hPullsConstJES .push_back(helper_hPullsConstJES );
    hPullsFSig     .push_back(helper_hPullsFSig     );
    hMasses        .push_back(helper_hMasses        );
    hMassesConstJES.push_back(helper_hMassesConstJES);
    hMassesFSig    .push_back(helper_hMassesFSig    );
  }
  
  double massError;
  tree->SetBranchAddress("massError", &massError);
  
  double massErrorConstJES;
  tree->SetBranchAddress("massConstJESError", &massErrorConstJES);
  
  double massErrorFSig;
  tree->SetBranchAddress("massFSigError", &massErrorFSig);
  
  double massPull;
  tree->SetBranchAddress("massPull", &massPull);
  
  double massPullConstJES;
  tree->SetBranchAddress("massConstJESPull", &massPullConstJES);
  
  double massPullFSig;
  tree->SetBranchAddress("massFSigPull", &massPullFSig);
  
  double mass;
  tree->SetBranchAddress("mass", &mass);
  
  double massConstJES;
  tree->SetBranchAddress("massConstJES", &massConstJES);
  
  double massFSig;
  tree->SetBranchAddress("massFSig", &massFSig);
  
  double genMass;
  tree->SetBranchAddress("genMass", &genMass);
  
  double genJES;
  tree->SetBranchAddress("genJES", &genJES);
  
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    int iMass = -1;
    int iJES  = -1;
    for(unsigned int idx = 0; idx < sizeof(genMasses)/sizeof(double); ++idx)
      if(genMass == genMasses[idx]) iMass = idx;
    for(unsigned int jdx = 0; jdx < sizeof(genJESes)/sizeof(double); ++jdx)
      if(genJES == genJESes[jdx]) iJES = jdx;

    if(iMass == -1 || iJES == -1) continue;

    if(mass < 0 || massError < 0 || massConstJES < 0 || massErrorConstJES < 0 || massFSig < 0 || massErrorFSig < 0) continue;

    if(mass != 172.5 || genMass == 172.5){
      hErrors        [iMass][iJES]->Fill(massError);
      hPulls         [iMass][iJES]->Fill(massPull);
      hMasses        [iMass][iJES]->Fill(mass);
    }
    if(massConstJES != 172.5 || genMass == 172.5){
      hErrorsConstJES[iMass][iJES]->Fill(massErrorConstJES);
      hPullsConstJES [iMass][iJES]->Fill(massPullConstJES);
      hMassesConstJES[iMass][iJES]->Fill(massConstJES);
    }
    if(massFSig != 172.5 || genMass == 172.5){
      hErrorsFSig    [iMass][iJES]->Fill(massErrorFSig);
      hPullsFSig     [iMass][iJES]->Fill(massPullFSig);
      hMassesFSig    [iMass][iJES]->Fill(massFSig);
    }
  }
  
  //// Do something :)
  
  for(unsigned int i = 0; i < sizeof(genMasses)/sizeof(double); ++i){
    for(unsigned int j = 0; j < sizeof(genJESes)/sizeof(double); ++j){

      TCanvas* cError = new TCanvas("cError", "cError", 600, 600);

      double max = hErrors[i][j]->GetMaximum();
      //hErrors[i][j]->GetXaxis()->SetRangeUser(0.6, 1.4);
      hErrors[i][j]->GetYaxis()->SetRangeUser(0, 1.2*max);
      hErrors[i][j]->Draw();
      hErrors[i][j]->Fit("gaus");//,"LIME");
      drawArrow(1.00, 1.05*max);
      DrawLabel("3.54 fb^{-1} collision data", 0.36, 0.85, 0.46, kRed+1);

      DrawCMSPrel();

      TString plotName = "plots/massError_"; plotName += 10*genMasses[i]; plotName += "_"; plotName += genJESes[j]==0.96 ? "096" : genJESes[j]==0.98 ? "098" : genJESes[j]==1. ? "100" : genJESes[j]==1.02 ? "102" : genJESes[j]==1.04 ? "104" : "999"; plotName += ".eps";
      cError->Print(plotName);

      TCanvas* cErrorConstJES = new TCanvas("cErrorConstJES", "cErrorConstJES", 600, 600);

      max = hErrorsConstJES[i][j]->GetMaximum();
      //hErrorsConstJES[i][j]->GetXaxis()->SetRangeUser(0.4, 1.);
      hErrorsConstJES[i][j]->GetYaxis()->SetRangeUser(0, 1.2*max);
      hErrorsConstJES[i][j]->Draw();
      hErrorsConstJES[i][j]->Fit("gaus");//,"LIME");
      drawArrow(0.69, 1.05*max);
      DrawLabel("3.54 fb^{-1} collision data", 0.36, 0.85, 0.46, kRed+1);

      DrawCMSPrel();

      plotName = "plots/massConstJESError_"; plotName += 10*genMasses[i]; plotName += "_"; plotName += genJESes[j]==0.96 ? "096" : genJESes[j]==0.98 ? "098" : genJESes[j]==1. ? "100" : genJESes[j]==1.02 ? "102" : genJESes[j]==1.04 ? "104" : "999"; plotName += ".eps";
      cErrorConstJES->Print(plotName);
 
      TCanvas* cErrorFSig = new TCanvas("cErrorFSig", "cErrorFSig", 600, 600);

      max = hErrorsFSig[i][j]->GetMaximum();
      //hErrorsFSig[i][j]->GetXaxis()->SetRangeUser(0.6, 1.4);
      hErrorsFSig[i][j]->GetYaxis()->SetRangeUser(0, 1.2*max);
      hErrorsFSig[i][j]->Draw();
      hErrorsFSig[i][j]->Fit("gaus");//,"LIME");
      drawArrow(1.00, 1.05*max);
      DrawLabel("3.54 fb^{-1} collision data", 0.36, 0.85, 0.46, kRed+1);

      DrawCMSPrel();

      plotName = "plots/massErrorFSig_"; plotName += 10*genMasses[i]; plotName += "_"; plotName += genJESes[j]==0.96 ? "096" : genJESes[j]==0.98 ? "098" : genJESes[j]==1. ? "100" : genJESes[j]==1.02 ? "102" : genJESes[j]==1.04 ? "104" : "999"; plotName += ".eps";
      cErrorFSig->Print(plotName);

      TCanvas* cPull = new TCanvas("cPull", "cPull", 600, 600);

      max = hPulls[i][j]->GetMaximum();
      hPulls[i][j]->GetXaxis()->SetRangeUser(-6., 6.);
      hPulls[i][j]->GetYaxis()->SetRangeUser(0, 1.2*max);
      hPulls[i][j]->Draw();
      hPulls[i][j]->Fit("gaus");//,"LIME");
      drawArrow(0, 1.05*max);
      DrawLabel("3.54 fb^{-1} collision data", 0.36, 0.85, 0.46, kRed+1);

      DrawCMSPrel();

      plotName = "plots/massPull_"; plotName += 10*genMasses[i]; plotName += "_"; plotName += genJESes[j]==0.96 ? "096" : genJESes[j]==0.98 ? "098" : genJESes[j]==1. ? "100" : genJESes[j]==1.02 ? "102" : genJESes[j]==1.04 ? "104" : "999"; plotName += ".eps";
      cPull->Print(plotName);

      TCanvas* cPullConstJES = new TCanvas("cPullConstJES", "cPullConstJES", 600, 600);

      max = hPullsConstJES[i][j]->GetMaximum();
      hPullsConstJES[i][j]->GetXaxis()->SetRangeUser(-10., 10.);
      hPullsConstJES[i][j]->GetYaxis()->SetRangeUser(0, 1.2*max);
      hPullsConstJES[i][j]->Draw();
      hPullsConstJES[i][j]->Fit("gaus");//,"LIME");
      drawArrow(0, 1.05*max);
      DrawLabel("3.54 fb^{-1} collision data", 0.36, 0.85, 0.46, kRed+1);

      DrawCMSPrel();

      plotName = "plots/massConstJESPull_"; plotName += 10*genMasses[i]; plotName += "_"; plotName += genJESes[j]==0.96 ? "096" : genJESes[j]==0.98 ? "098" : genJESes[j]==1. ? "100" : genJESes[j]==1.02 ? "102" : genJESes[j]==1.04 ? "104" : "999"; plotName += ".eps";
      cPullConstJES->Print(plotName);
 
      TCanvas* cPullFSig = new TCanvas("cPullFSig", "cPullFSig", 600, 600);

      max = hPullsFSig[i][j]->GetMaximum();
      hPullsFSig[i][j]->GetXaxis()->SetRangeUser(-6., 6.);
      hPullsFSig[i][j]->GetYaxis()->SetRangeUser(0, 1.2*max);
      hPullsFSig[i][j]->Draw();
      hPullsFSig[i][j]->Fit("gaus");//,"LIME");
      drawArrow(0, 1.05*max);
      DrawLabel("3.54 fb^{-1} collision data", 0.36, 0.85, 0.46, kRed+1);

      DrawCMSPrel();

      plotName = "plots/massPullFSig_"; plotName += 10*genMasses[i]; plotName += "_"; plotName += genJESes[j]==0.96 ? "096" : genJESes[j]==0.98 ? "098" : genJESes[j]==1. ? "100" : genJESes[j]==1.02 ? "102" : genJESes[j]==1.04 ? "104" : "999"; plotName += ".eps";
      cPullFSig->Print(plotName);

      TCanvas* cMass = new TCanvas("cMass", "cMass", 600, 600);

      max = hMasses[i][j]->GetMaximum();
      hMasses[i][j]->GetXaxis()->SetRangeUser(155, 190);
      hMasses[i][j]->GetYaxis()->SetRangeUser(0, 1.2*max);
      hMasses[i][j]->Draw();
      hMasses[i][j]->Fit("gaus");//,"LIME");
      drawArrow(genMasses[i], 1.05*max);
      DrawLabel("input mass from MC", 0.36, 0.85, 0.46, kRed+1);

      DrawCMSPrel();

      plotName = "plots/massMass_"; plotName += 10*genMasses[i]; plotName += "_"; plotName += genJESes[j]==0.96 ? "096" : genJESes[j]==0.98 ? "098" : genJESes[j]==1. ? "100" : genJESes[j]==1.02 ? "102" : genJESes[j]==1.04 ? "104" : "999"; plotName += ".eps";
      cMass->Print(plotName);

      TCanvas* cMassConstJES = new TCanvas("cMassConstJES", "cMassConstJES", 600, 600); plotName += genJESes[j]==0.96 ? "096" : genJESes[j]==0.98 ? "098" : genJESes[j]==1. ? "100" : genJESes[j]==1.02 ? "102" : genJESes[j]==1.04 ? "104" : "999"; plotName += ".eps";

      max = hMassesConstJES[i][j]->GetMaximum();
      hMassesConstJES[i][j]->GetXaxis()->SetRangeUser(155, 195);
      hMassesConstJES[i][j]->GetYaxis()->SetRangeUser(0, 1.2*max);
      hMassesConstJES[i][j]->Draw();
      hMassesConstJES[i][j]->Fit("gaus");//,"LIME");
      drawArrow(genMasses[i], 1.05*max);
      DrawLabel("input mass from MC", 0.36, 0.85, 0.46, kRed+1);

      DrawCMSPrel();

      plotName = "plots/massConstJESMass_"; plotName += 10*genMasses[i]; plotName += "_"; plotName += genJESes[j]==0.96 ? "096" : genJESes[j]==0.98 ? "098" : genJESes[j]==1. ? "100" : genJESes[j]==1.02 ? "102" : genJESes[j]==1.04 ? "104" : "999"; plotName += ".eps";
      cMassConstJES->Print(plotName);
 
      TCanvas* cMassFSig = new TCanvas("cMassFSig", "cMassFSig", 600, 600);

      max = hMassesFSig[i][j]->GetMaximum();
      hMassesFSig[i][j]->GetXaxis()->SetRangeUser(155, 190);
      hMassesFSig[i][j]->GetYaxis()->SetRangeUser(0, 1.2*max);
      hMassesFSig[i][j]->Draw();
      hMassesFSig[i][j]->Fit("gaus");//,"LIME");
      drawArrow(genMasses[i], 1.05*max);
      DrawLabel("input mass from MC", 0.36, 0.85, 0.46, kRed+1);

      DrawCMSPrel();

      plotName = "plots/massMassFSig_"; plotName += 10*genMasses[i]; plotName += "_"; plotName += genJESes[j]==0.96 ? "096" : genJESes[j]==0.98 ? "098" : genJESes[j]==1. ? "100" : genJESes[j]==1.02 ? "102" : genJESes[j]==1.04 ? "104" : "999"; plotName += ".eps";
      cMassFSig->Print(plotName);

    }
  }
}
