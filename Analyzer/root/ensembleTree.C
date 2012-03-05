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
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"

#include "tdrstyle.C"

enum styles          { kDown, kNominal, kUp};
int color_   [ 3 ] = { kRed+1, kBlue+1, kGreen+1};
int marker_  [ 3 ] = { 23, 20, 22};

double genMass[]      = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5, 184.5};
double genMassError[] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
double genMassN[]     = {2, 2, 2, 2, 4, 2, 2, 2, 2};
  
double genJES[]       = {0.96, 1.00, 1.04};
double genJESError[]  = {1e-6, 1e-6, 1e-6};
  
double mass[3][9];
double massBias[3][9];
double massBiasError[3][9];

double massPull[3][9];
double massPullError[3][9];
  
double JES[3][9];
double JESBias[3][9];
double JESBiasError[3][9];

double JESPull[3][9];
double JESPullError[3][9];

std::vector<TGraphErrors*> gMass;
std::vector<TGraphErrors*> gJES;

std::vector<TGraphErrors*> gMassPull;
std::vector<TGraphErrors*> gJESPull;

TTree* tree;

TF1* constFit;


void DrawLegend() {
  TLegend *legMass = new TLegend(0.74, 0.7, 0.97, 0.94);
  legMass->SetFillColor(kWhite);
  legMass->SetBorderSize(1);
  legMass->AddEntry( gMass[0], "JES=0.96", "EP");
  legMass->AddEntry( gMass[1], "JES=1.00", "EP");
  legMass->AddEntry( gMass[2], "JES=1.04", "EP");
  legMass->AddEntry( constFit, "Const. fit", "L");
  legMass->Draw();
}

void ensembleTree()
{
  //*
  TStyle *tdrStyle = setTDRStyle();
  tdrStyle->SetPadGridX(true);
  tdrStyle->SetPadGridY(true);
  tdrStyle->SetOptFit(0);
  //*/
  
  TCanvas* canvasFit = new TCanvas("canvasFit", "hadronic top h2Mass", 500, 500);
  canvasFit->cd();
  
  //// Get histos
  
  TFile* fEnsemble = new TFile("/scratch/hh/current/cms/user/mseidel/topmass_120207_1618/ensemble.root");
  
  tree = (TTree*) fEnsemble->Get("tree");
  
  TMultiGraph *mgMass = new TMultiGraph();
  mgMass->SetTitle(";m_{t,gen} [GeV];m_{t,meas}-m_{t,gen} [GeV]");
  
  TMultiGraph *mgJES = new TMultiGraph();
  mgJES->SetTitle(";m_{t,gen} [GeV];JES_{meas}-JES_{gen}");
  
  TMultiGraph *mgMassPull = new TMultiGraph();
  mgMassPull->SetTitle(";m_{t,gen} [GeV];mass pull width");
  
  TMultiGraph *mgJESPull = new TMultiGraph();
  mgJESPull->SetTitle(";m_{t,gen} [GeV];JES pull width");
  
  for (int iJES = 0; iJES < 3; iJES++) {
    for (int iMass = 0; iMass < 9; iMass++) {
      TString sel("mass>0 & JES>0 & genMass=="); sel+=genMass[iMass]; sel+=" & genJES=="; sel+=genJES[iJES];
      
      TF1* gausMassBias = new TF1("gausMassBias", "gaus");
      tree->Fit("gausMassBias", "mass", sel, "Q0");
      
      mass[iJES][iMass]          = genMass[iMass];
      massBias[iJES][iMass]      = gausMassBias->GetParameter(1) - genMass[iMass];
      massBiasError[iJES][iMass] = gausMassBias->GetParameter(2) / sqrt(genMassN[iMass]);
        
      TF1* gausJESBias = new TF1("gausJESBias", "gaus");
      tree->Fit("gausJESBias", "JES", sel, "Q0");
      
      JES[iJES][iMass]          = genJES[iJES];
      JESBias[iJES][iMass]      = gausJESBias->GetParameter(1) - genJES[iJES];
      JESBiasError[iJES][iMass] = gausJESBias->GetParameter(2) / sqrt(genMassN[iMass]);
      
      
      TF1* gausMassPull = new TF1("gausMassPull", "gaus");
      tree->Fit("gausMassPull", "massPull", sel, "Q0");
      
      massPull[iJES][iMass]      = gausMassPull->GetParameter(2);
      massPullError[iJES][iMass] = sqrt(1./2. * (1./(5144.*genMassN[iMass])+1./2000));
      
      //std::cout << sqrt(1./2. * (1./(5144.*genMassN[iMass])+1./2000)) << std::endl;
        
      TF1* gausJESPull = new TF1("gausJESPull", "gaus");
      tree->Fit("gausJESPull", "JESPull", sel, "Q0");
      
      JESPull[iJES][iMass]      = gausJESPull->GetParameter(2);
      JESPullError[iJES][iMass] = sqrt(1./2. * (1./(5144.*genMassN[iMass])+1./2000));
    }
    
    double genMassMod[9];
		for (int i = 0; i < 9; i++) {
			genMassMod[i] = genMass[i] + 0.2*(iJES-1);
		}
		
  	gMass.push_back(new TGraphErrors(9, genMassMod, massBias[iJES], genMassError, massBiasError[iJES]));
		gMass[iJES]->SetMarkerStyle(marker_[iJES]);
		gMass[iJES]->SetMarkerColor(color_ [iJES]);
		gMass[iJES]->SetLineColor  (color_ [iJES]);
		mgMass->Add(gMass[iJES]);
		
		gJES.push_back(new TGraphErrors(9, genMassMod, JESBias[iJES], genMassError, JESBiasError[iJES]));
		gJES[iJES]->SetMarkerStyle(marker_[iJES]);
		gJES[iJES]->SetMarkerColor(color_ [iJES]);
		gJES[iJES]->SetLineColor  (color_ [iJES]);
		mgJES->Add(gJES[iJES]);
		
		gMassPull.push_back(new TGraphErrors(9, genMassMod, massPull[iJES], genMassError, massPullError[iJES]));
		gMassPull[iJES]->SetMarkerStyle(marker_[iJES]);
		gMassPull[iJES]->SetMarkerColor(color_ [iJES]);
		gMassPull[iJES]->SetLineColor  (color_ [iJES]);
		mgMassPull->Add(gMassPull[iJES]);
		
		gJESPull.push_back(new TGraphErrors(9, genMassMod, JESPull[iJES], genMassError, JESPullError[iJES]));
		gJESPull[iJES]->SetMarkerStyle(marker_[iJES]);
		gJESPull[iJES]->SetMarkerColor(color_ [iJES]);
		gJESPull[iJES]->SetLineColor  (color_ [iJES]);
		mgJESPull->Add(gJESPull[iJES]);
  }
  
  constFit = new TF1("constFit", "[0]");
  constFit->SetParNames("offset");
  constFit->SetLineColor(kBlack);
  constFit->SetLineWidth(3);
  constFit->SetLineStyle(7);
  
  canvasFit->Clear();
  mgMass->Fit("constFit", "EM");
  mgMass->SetMinimum(-5);
  mgMass->SetMaximum( 5);
  mgMass->Draw("AP");
  canvasFit->Update();
  /*
  TPaveStats* statsGlobal = (TPaveStats*) mgMass->GetListOfFunctions()->FindObject("stats");
  std::cout << "got stats object" << std::endl;
  statsGlobal->SetX1NDC(0.555);
  statsGlobal->SetY1NDC(0.825);
  statsGlobal->SetX2NDC(0.95);
  statsGlobal->SetY2NDC(0.95);
  canvasFit->Modified();
  //*/
  DrawLegend(); DrawCMSSim();
  canvasFit->Print("fit_Mass.eps");
  
  canvasFit->Clear();
  mgJES->Fit("constFit", "EM");
  mgJES->SetMinimum(-0.05);
  mgJES->SetMaximum( 0.05);
  mgJES->Draw("AP");
  DrawLegend(); DrawCMSSim();
  canvasFit->Print("fit_JES.eps");
  
  canvasFit->Clear();
  mgMassPull->SetMinimum(0.5);
  mgMassPull->Fit("constFit", "EM");
  mgMassPull->SetMaximum(1.5);
  mgMassPull->Draw("AP");
  DrawLegend(); DrawCMSSim();
  canvasFit->Print("fit_MassPull.eps");
  
  canvasFit->Clear();
  mgJESPull->Fit("constFit", "EM");
  mgJESPull->SetMinimum(0.5);
  mgJESPull->SetMaximum(1.5);
  mgJESPull->Draw("AP");
  DrawLegend(); DrawCMSSim();
  canvasFit->Print("fit_JESPull.eps");
}


