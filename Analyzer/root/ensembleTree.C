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

bool muon = true;

enum styles          { kDown, kNominal, kUp};
int color_   [ 3 ] = { kRed+1, kBlue+1, kGreen+1};
int marker_  [ 3 ] = { 23, 20, 22};

double genMass[]      = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5, 184.5};
double genMassError[] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
double genMassN[]     = {1620072, 1633197, 1669034, 1606570, 59613991, 1538301, 1648519, 1665350, 1671859};
//double genMassN[]     = {1620072, 1.5, 1.5, 1.5, 59613991, 1.5, 1.5, 1.5, 1.5};
double maxMCWeight[]  = {1.7, 2.2, 2.2, 2.2, 1.7, 2.2, 2.2, 2.2, 1.7};
double crossSection   = 164.4;
double peLumi         = 5000.;
  
double genJES[]       = {0.96, 1.00, 1.04};
double genJESError[]  = {1e-6, 1e-6, 1e-6};
double measJES[3];

double offset[3];
double offsetError[3];
  
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
TF1* linearFit096;
TF1* linearFit100;
TF1* linearFit104;
TF1* linearFitJES;
TFitResultPtr fitResult;


void DrawLegend() {
  TLegend *leg = new TLegend(0.74, 0.7, 0.97, 0.94);
  leg->SetTextSizePixels(12);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->AddEntry( gMass[0], "JES=0.96", "EP");
  leg->AddEntry( gMass[1], "JES=1.00", "EP");
  leg->AddEntry( gMass[2], "JES=1.04", "EP");
  leg->AddEntry( constFit, "Const. fit", "L");
  char chi2[6]; sprintf(chi2, "%3.1f", fitResult->Chi2());
  TString sConstFit("#chi^{2}/ndf="); sConstFit+=chi2; sConstFit+="/"; sConstFit+=fitResult->Ndf();
  leg->AddEntry((TObject*)0, sConstFit, "");
  leg->Draw();
}

void ensembleTree()
{
  //*
  TStyle *tdrStyle = setTDRStyle();
  tdrStyle->SetPadGridX(true);
  tdrStyle->SetPadGridY(true);
  tdrStyle->SetOptFit(0);
  //*/
  
  TCanvas* canvasFit = new TCanvas("canvasFit", "mt-JES measurement calibration", 500, 500);
  canvasFit->cd();
  
  //// Get histos
  TString sFile("/scratch/hh/current/cms/user/mseidel/topmass_120503_2140c/");
  if (muon) sFile += "muon"; else sFile += "electron"; sFile += "/ensemble.root";
  TFile* fEnsemble = new TFile(sFile);
  
  tree = (TTree*) fEnsemble->Get("tree");
  
  TMultiGraph *mgMass = new TMultiGraph();
  mgMass->SetTitle(";m_{t,gen} [GeV];m_{t,meas}-m_{t,gen} [GeV]");
  
  TMultiGraph *mgJES = new TMultiGraph();
  mgJES->SetTitle(";m_{t,gen} [GeV];JES_{meas}-JES_{gen}");
  
  TMultiGraph *mgMassPull = new TMultiGraph();
  mgMassPull->SetTitle(";m_{t,gen} [GeV];mass pull width");
  
  TMultiGraph *mgJESPull = new TMultiGraph();
  mgJESPull->SetTitle(";m_{t,gen} [GeV];JES pull width");
  
  TH2D* h2Mass = new TH2D("h2Mass", "h2Mass", 1000, 150, 200, 1000, 0.9, 1.1);
  TH2D* h2JES = new TH2D("h2JES", "h2JES", 1000, 150, 200, 1000, 0.9, 1.1);
  
  for (int iJES = 0; iJES < 3; iJES++) {
    for (int iMass = 0; iMass < 9; iMass++) {
      TString sel("mass>0 & JES>0 & genMass=="); sel+=genMass[iMass]; sel+=" & genJES=="; sel+=genJES[iJES];
      double entries = tree->GetEntries(sel);
      
      TF1* gausMassBias = new TF1("gausMassBias", "gaus");
      tree->Fit("gausMassBias", "mass", sel, "Q0");
      //tree->Fit("gausMassBias", "mass+2.25341e-01+0.2/0.04*(1-JES)", sel, "Q0");
      
      mass[iJES][iMass]          = genMass[iMass];
      massBias[iJES][iMass]      = gausMassBias->GetParameter(1) - genMass[iMass];
      massBiasError[iJES][iMass] = gausMassBias->GetParameter(2) / sqrt(genMassN[iMass]/(crossSection*peLumi*maxMCWeight[iMass]));
        
      TF1* gausJESBias = new TF1("gausJESBias", "gaus");
      tree->Fit("gausJESBias", "JES", sel, "Q0");
      //tree->Fit("gausJESBias", "JES + 2.66486e-03 + 0.002/0.035*(JES-1)", sel, "Q0");
      
      JES[iJES][iMass]          = genJES[iJES];
      JESBias[iJES][iMass]      = gausJESBias->GetParameter(1) - genJES[iJES];
      JESBiasError[iJES][iMass] = gausJESBias->GetParameter(2) / sqrt(genMassN[iMass]/(crossSection*peLumi*maxMCWeight[iMass]));
      
      // Fill calibration histos
      
      h2Mass->Fill(gausMassBias->GetParameter(1), gausJESBias->GetParameter(1), massBias[iJES][iMass]);
      h2Mass->SetBinError(h2Mass->FindBin(gausMassBias->GetParameter(1), gausJESBias->GetParameter(1)), massBiasError[iJES][iMass]);
      
      h2JES->Fill(gausMassBias->GetParameter(1), gausJESBias->GetParameter(1), JESBias[iJES][iMass]);
      h2JES->SetBinError(h2JES->FindBin(gausMassBias->GetParameter(1), gausJESBias->GetParameter(1)), JESBiasError[iJES][iMass]);
      
      TF1* gausMassPull = new TF1("gausMassPull", "gaus");
      tree->Fit("gausMassPull", "massPull", sel, "Q0");
      
      massPull[iJES][iMass]      = gausMassPull->GetParameter(2);
      massPullError[iJES][iMass] = sqrt(1./2. * (1./(5144.*genMassN[iMass]/(crossSection*peLumi*maxMCWeight[iMass]))+1./entries));
      
      //std::cout << sqrt(1./2. * (1./(5144.*genMassN[iMass])+1./2000)) << std::endl;
        
      TF1* gausJESPull = new TF1("gausJESPull", "gaus");
      tree->Fit("gausJESPull", "JESPull", sel, "Q0");
      
      JESPull[iJES][iMass]      = gausJESPull->GetParameter(2);
      JESPullError[iJES][iMass] = sqrt(1./2. * (1./(5144.*genMassN[iMass]/(crossSection*peLumi*maxMCWeight[iMass]))+1./entries));
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
  
  linearFit096 = new TF1("linearFit096", "[0]+[1]*(x-172.5)");
  linearFit096->SetParNames("offset", "slope");
  linearFit096->SetLineColor(color_[0]);
  linearFit096->SetLineWidth(1);
  linearFit096->SetLineStyle(7);
  
  linearFit100 = new TF1("linearFit100", "[0]+[1]*(x-172.5)");
  linearFit100->SetParNames("offset", "slope");
  linearFit100->SetLineColor(color_[1]);
  linearFit100->SetLineWidth(1);
  linearFit100->SetLineStyle(7);
  
  linearFit104 = new TF1("linearFit104", "[0]+[1]*(x-172.5)");
  linearFit104->SetParNames("offset", "slope");
  linearFit104->SetLineColor(color_[2]);
  linearFit104->SetLineWidth(1);
  linearFit104->SetLineStyle(7);
	
	linearFitJES = new TF1("linearFitJES", "[0]+[1]*(x-1.0)");
  linearFitJES->SetParNames("offset", "slope");
  
	
	canvasFit->Clear();
  
  fitResult = mgJES->Fit("constFit", "EMS");
  gJES[0]->Fit("linearFit096", "EM");
  gJES[1]->Fit("linearFit100", "EM");
  gJES[2]->Fit("linearFit104", "EM");
  mgJES->SetMinimum(-0.05);
  mgJES->SetMaximum( 0.05);
  mgJES->Draw("AP");
  DrawLegend(); if (muon) DrawCMSSimMuon(); else DrawCMSSimElectron();
  if (muon) canvasFit->Print("fit_JES.eps");
  else canvasFit->Print("fit_JES_electron.eps");
	
	canvasFit->Clear();
	
	offset[0] = linearFit096->GetParameter(0);
	offset[1] = linearFit100->GetParameter(0);
	offset[2] = linearFit104->GetParameter(0);
	offsetError[0] = linearFit096->GetParError(0);
	offsetError[1] = linearFit100->GetParError(0);
	offsetError[2] = linearFit104->GetParError(0);
	measJES[0] = linearFit096->GetParameter(0)+genJES[0];
	measJES[1] = linearFit100->GetParameter(0)+genJES[1];
	measJES[2] = linearFit104->GetParameter(0)+genJES[2];
	
	TGraphErrors* gOffset = new TGraphErrors(3, measJES, offset, genJESError, offsetError);
	gOffset->Fit("linearFitJES", "EM");
  gOffset->Draw("AP");
  
  
  canvasFit->Clear();
  fitResult = mgMass->Fit("constFit", "EMS");
  gMass[0]->Fit("linearFit096", "EM");
  gMass[1]->Fit("linearFit100", "EM");
  gMass[2]->Fit("linearFit104", "EM");
  mgMass->SetMinimum(-5);
  mgMass->SetMaximum( 5);
  mgMass->Draw("AP");
  canvasFit->Update();
  
  DrawLegend(); if (muon) DrawCMSSimMuon(); else DrawCMSSimElectron();
  if (muon) canvasFit->Print("fit_Mass.eps");
  else canvasFit->Print("fit_Mass_electron.eps");
  
  canvasFit->Clear();
	
	offset[0] = linearFit096->GetParameter(0);
	offset[1] = linearFit100->GetParameter(0);
	offset[2] = linearFit104->GetParameter(0);
	offsetError[0] = linearFit096->GetParError(0);
	offsetError[1] = linearFit100->GetParError(0);
	offsetError[2] = linearFit104->GetParError(0);
	
	gOffset = new TGraphErrors(3, measJES, offset, genJESError, offsetError);
	gOffset->Fit("linearFitJES", "EM");
	gOffset->Draw("AP");
	
	TF2* fit2D = new TF2("fit2D", "[0] + [1]*(x-172.5) + [2]*(y-1) + [3]*(x-172.5)*(y-1)");
  fit2D->SetParNames("offset", "slopeMass", "slopeJES");
  
  //*
  std::cout << "=== 2D calibration - mass ===" << std::endl;
  h2Mass->Fit("fit2D");
  h2Mass->Draw("SURF1");
  
  std::cout << "=== 2D calibration - JES ===" << std::endl;
  h2JES->Fit("fit2D");
  h2JES->Draw("SURF1");
  //*/
  
  //*
  canvasFit->Clear();
  mgMassPull->Fit("constFit", "EM");
  gMassPull[0]->Fit("linearFit096", "EM0");
  gMassPull[1]->Fit("linearFit100", "EM0");
  gMassPull[2]->Fit("linearFit104", "EM0");
  mgMassPull->SetMinimum(0.5);
  mgMassPull->SetMaximum(1.5);
  mgMassPull->Draw("AP");
  DrawLegend(); if (muon) DrawCMSSimMuon(); else DrawCMSSimElectron();
  if (muon) canvasFit->Print("fit_MassPull.eps");
  else canvasFit->Print("fit_MassPull_electron.eps");
  
  canvasFit->Clear();
  mgJESPull->Fit("constFit", "EM");
  gJESPull[0]->Fit("linearFit096", "EM0");
  gJESPull[1]->Fit("linearFit100", "EM0");
  gJESPull[2]->Fit("linearFit104", "EM0");
  mgJESPull->SetMinimum(0.5);
  mgJESPull->SetMaximum(1.5);
  mgJESPull->Draw("AP");
  DrawLegend(); if (muon) DrawCMSSimMuon(); else DrawCMSSimElectron();
  if (muon) canvasFit->Print("fit_JESPull.eps");
  else canvasFit->Print("fit_JESPull_electron.eps");
  //*/
}


