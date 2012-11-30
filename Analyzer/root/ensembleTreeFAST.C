#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1F.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TBox.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TFitResult.h"

#include "tdrstyle.C"

const int nJES = 3;
enum styles          { kDown, kNominal, kUp };
int color_ [] = { kRed+1, kBlue+1, kGreen+1 };
int marker_[] = { 23, 20, 22 };
double genJES[]       = {0.96, 1.00, 1.04};
double genJESError[]  = {1e-6, 1e-6, 1e-6};
double genMass      = 172.5;
double genMassError = 1e-6;
double genMassN     = 59613991.;
double maxMCWeight  = 2.2;

double crossSection   = 164.4;
double peLumi         = 3544.844;
  
double measJES[nJES];

double offset[nJES];
double offsetError[nJES];
  
double mass[nJES];
double massBias[nJES];
double massBiasError[nJES];

double massPull[nJES];
double massPullError[nJES];
  
double JES[nJES];
double JESBias[nJES];
double JESBiasError[nJES];

double JESPull[nJES];
double JESPullError[nJES];

void ensembleTreeFAST()
{
  /*
  TStyle *tdrStyle = setTDRStyle();
  tdrStyle->SetPadGridX(true);
  tdrStyle->SetPadGridY(true);
  tdrStyle->SetOptFit(0);
  tdrStyle->SetTitleSize(0.09, "XYZ");
  tdrStyle->SetLabelSize(0.08, "XYZ");
  tdrStyle->SetPadTopMargin(0.2);
  tdrStyle->SetPadBottomMargin(0.2);
  tdrStyle->SetTitleYOffset(1.);
  //*/
  
  //// Get histos

  TString sFile("/scratch/hh/current/cms/user/eschliec/TopMass/21/ensemble/");
  //TString sFile("/scratch/hh/dust/naf/cms/user/eschliec/TopMass/");
  sFile += "ensemble_F11_Ideogram_Calibrated_1_FAST.root";
  TFile* fEnsemble = new TFile(sFile);
  
  TTree* tree = (TTree*) fEnsemble->Get("tree");

  TH1F * hMass = new TH1F("hMass", "hMass", 1000, 0.95, 1.05);
  TH1F * hMassPull = new TH1F("hMassPull", "hMassPull", 1000, 0.95, 1.05);
  
  TH1F * hJES = new TH1F("hJES", "hJES", 1000, 0.95, 1.05);
  TH1F * hJESPull = new TH1F("hJESPull", "hJESPull", 1000, 0.95, 1.05);

  for (int iJES = 0; iJES < nJES; iJES++) {
    TString sel("mass_mTop_JES>0 & JES_mTop_JES>0 & abs(mass_mTop_JES_Pull)<100 & genMass=="); sel+=genMass; sel+=" & genJES=="; sel+=genJES[iJES];
    double entries = tree->GetEntries(sel);
      
    TF1* gausMassBias = new TF1("gausMassBias", "gaus");
    tree->Fit("gausMassBias", "mass_mTop_JES", sel, "Q0");
      
    mass[iJES]          = genMass;
    massBias[iJES]      = gausMassBias->GetParameter(1) - genMass;
    massBiasError[iJES] = gausMassBias->GetParameter(2) / sqrt(genMassN/(crossSection*peLumi*maxMCWeight));
        
    TF1* gausJESBias = new TF1("gausJESBias", "gaus");
    tree->Fit("gausJESBias", "JES_mTop_JES", sel, "Q0");
      
    JES[iJES]          = genJES[iJES];
    JESBias[iJES]      = gausJESBias->GetParameter(1) - genJES[iJES];
    JESBiasError[iJES] = gausJESBias->GetParameter(2) / sqrt(genMassN/(crossSection*peLumi*maxMCWeight));
      
    // Fill calibration histos
      
    TF1* gausMassPull = new TF1("gausMassPull", "gaus");
    tree->Fit("gausMassPull", "mass_mTop_JES_Pull", sel, "Q0");
      
    double eff = 0.00135403611570646;
    massPull[iJES]      = gausMassPull->GetParameter(2);
    massPullError[iJES] = sqrt(1./2. * (maxMCWeight/(genMassN*eff) + 1./(entries-1.)));
      
    //std::cout << sqrt(1./2. * (1./(5144.*genMassN[iMass])+1./2000)) << std::endl;
        
    TF1* gausJESPull = new TF1("gausJESPull", "gaus");
    tree->Fit("gausJESPull", "JES_mTop_JES_Pull", sel, "Q0");
      
    JESPull[iJES]      = gausJESPull->GetParameter(2);
    JESPullError[iJES] = sqrt(1./2. * maxMCWeight/(genMassN*eff) + 1./entries);
  }
  
  for(int i = 0; i < nJES; ++i){
    hMass    ->SetBinContent(hMass->FindBin(genJES[i]),massBias[i]);
    hMass    ->SetBinError  (hMass->FindBin(genJES[i]),massBiasError[i]);
    hMassPull->SetBinContent(hMass->FindBin(genJES[i]),massPull[i]);
    hMassPull->SetBinError  (hMass->FindBin(genJES[i]),massPullError[i]);
    
    hJES    ->SetBinContent(hJES->FindBin(genJES[i]),JESBias[i]);
    hJES    ->SetBinError  (hJES->FindBin(genJES[i]),JESBiasError[i]);
    hJESPull->SetBinContent(hJES->FindBin(genJES[i]),JESPull[i]);
    hJESPull->SetBinError  (hJES->FindBin(genJES[i]),JESPullError[i]);
  }

  TF1* constFit = new TF1("constFit", "[0]");
  constFit->SetParNames("offset");
  constFit->SetLineColor(kBlack);
  constFit->SetLineWidth(3);
  constFit->SetLineStyle(7);

  TF1* linearFitJES = new TF1("linearFitJES", "[0]+[1]*(x-1.0)");
  linearFitJES->SetParNames("offset", "slope");
  linearFitJES->SetLineColor(kBlack);
  linearFitJES->SetLineWidth(3);
  linearFitJES->SetLineStyle(7);

  
  new TCanvas("canv1","Mass Calibration",1,1,600,600);
  std::cout << "=== calibration - Mass ===" << std::endl;
  hMass->Fit("linearFitJES", "EM");
  hMassPull->Fit("constFit", "EMN0");

  new TCanvas("canv2","JES Calibration",641,1,600,600);
  std::cout << "=== calibration - JES ===" << std::endl;
  hJES ->Fit("linearFitJES", "EM");
  hJESPull ->Fit("constFit", "EMN0");
}

