#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
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

#include "tdrstyle_new.C"
#include "CMS_lumi.C"

enum lepton           { kElectron, kMuon, kAll, kMuon_BReg};
std::string lepton_ [4] = { "electron", "muon", "lepton", "muon_BReg"};

int channel = 2;
//bool calibrated = false;
//std::string suffix = "_uncalibrated";
bool calibrated = true;
std::string suffix = "_calibrated_pull";

Long64_t nentries = 1000000000; //1000*27;
Long64_t firstentry = nentries*0 + 0;

enum styles          { kDown, kNominal, kUp};
int color_   [ 3 ] = { kRed+1, kBlue+1, kGreen+1};
int marker_  [ 3 ] = { 23, 20, 22};

double genMass[]      = {166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5};
double genMassError[] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
double genMassN[]     = {27000000, 41000000, 25000000, 62000000, 27000000, 40000000, 24000000};
//double genMassN[]     = {1620072, 1.5, 1.5, 1.5, 59613991, 1.5, 1.5, 1.5, 1.5};
double maxMCWeight[]  = {1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};
double crossSection   = 245.8;
double peLumi         = 19700.;
double eff;

std::string samples[21] = {"Summer12_TTJetsMS1665_0.98", "Summer12_TTJetsMS1715_0.98", "Summer12_TTJetsMS1735_0.98", "Summer12_TTJetsMS1785_0.98", "Summer12_TTJetsMS1665_1.00", "Summer12_TTJetsMS1715_1.00", "Summer12_TTJetsMS1735_1.00", "Summer12_TTJetsMS1785_1.00", "Summer12_TTJetsMS1665_1.02", "Summer12_TTJetsMS1715_1.02", "Summer12_TTJetsMS1735_1.02", "Summer12_TTJetsMS1785_1.02", "Summer12_TTJetsMS1695_0.98", "Summer12_TTJetsMS1725_0.98", "Summer12_TTJetsMS1755_0.98", "Summer12_TTJetsMS1695_1.00", "Summer12_TTJetsMS1725_1.00", "Summer12_TTJetsMS1755_1.00", "Summer12_TTJetsMS1695_1.02", "Summer12_TTJetsMS1725_1.02", "Summer12_TTJetsMS1755_1.02"};
  
double genJES[]       = {0.98, 1.00, 1.02};
double genJESError[]  = {1e-6, 1e-6, 1e-6};
double measJES[3];

double offset[3];
double offsetError[3];
  
double mass[3][7];
double massBias[3][7];
double massBiasError[3][7];

double massPull[3][7];
double massPullError[3][7];
  
double JES[3][7];
double JESBias[3][7];
double JESBiasError[3][7];

double JESPull[3][7];
double JESPullError[3][7];

std::vector<TGraphErrors*> gMass;
std::vector<TGraphErrors*> gJES;

std::vector<TGraphErrors*> gMassPull;
std::vector<TGraphErrors*> gJESPull;

TChain* tree;

TF1* constFit;
TF1* linearFit096;
TF1* linearFit100;
TF1* linearFit104;
TF1* linearFitJES;
TFitResultPtr fitResult;


void DrawLegend(bool bwlines = false) {
  double y1 = 0.84; double y2 = 0.9;
  int textsize = 10;
  TLegend *leg1 = new TLegend(0.25, y1, 0.92, y2);
  leg1->SetNColumns(3);
  leg1->SetTextSizePixels(textsize);
  leg1->SetFillColor(kWhite);
  leg1->SetBorderSize(0);
  if (bwlines) {
    leg1->AddEntry( linearFit096, "JSF=0.98", "LP");
    leg1->AddEntry( linearFit100, "JSF=1.00", "LP");
    leg1->AddEntry( linearFit104, "JSF=1.02", "LP");
  }
  else {
    leg1->AddEntry( linearFit096, "JSF=0.98", "P");
    leg1->AddEntry( linearFit100, "JSF=1.00", "P");
    leg1->AddEntry( linearFit104, "JSF=1.02", "P");
  }
  /*
  leg->AddEntry( constFit, "Const. fit", "L");
  char chi2[6]; sprintf(chi2, "%3.1f", fitResult->Chi2());
  TString sConstFit("#chi^{2}/ndf="); sConstFit+=chi2; sConstFit+="/"; sConstFit+=fitResult->Ndf();
  leg->AddEntry((TObject*)0, sConstFit, "");
  */
  leg1->Draw();
  
  /*
  TLegend *leg2 = new TLegend(0.44, y1, 0.68, y2);
  leg2->SetTextSizePixels(textsize);
  leg2->SetFillColor(kWhite);
  leg2->SetBorderSize(0);
  if (bwlines) leg2->AddEntry( linearFit100, "JES=1.00", "LP");
  else         leg2->AddEntry( linearFit100, "JES=1.00", "P");
  leg2->Draw();
  
  TLegend *leg3 = new TLegend(0.68, y1, 0.92, y2);
  leg3->SetTextSizePixels(textsize);
  leg3->SetFillColor(kWhite);
  leg3->SetBorderSize(0);
  if (bwlines) leg3->AddEntry( linearFit104, "JES=1.02", "LP");
  else         leg3->AddEntry( linearFit104, "JES=1.02", "P");
  leg3->Draw();
  */
}

void ensembleTreeLeptonJets(std::string pathToPE = "/nfs/dust/cms/user/mseidel/pseudoexperiments/topmass_paper")
{
  //*
  TStyle *tdrStyle = setTDRStyle();
  //tdrStyle->SetPadGridX(true);
  //tdrStyle->SetPadGridY(true);
  tdrStyle->SetOptFit(0);
  tdrStyle->SetTitleSize(0.09, "XYZ");
  tdrStyle->SetLabelSize(0.08, "XYZ");
  tdrStyle->SetPadTopMargin(0.16);
  tdrStyle->SetPadBottomMargin(0.24);
  tdrStyle->SetTitleYOffset(1.);
  //*/
  
  TCanvas* canvasFit = new TCanvas("canvasFit", "mt-JES measurement calibration", 500, 500);
  
  canvasFit->cd();
  
  //// Get histos
  // topmass_120503_2140 - e mu
  // topmass_120530_0840 - all
  // topmass_120412_2120cp
  tree = new TChain("tree");
  
  for (int s = 0; s < 21; ++s) {
    std::string pathToRootFiles = pathToPE+suffix+std::string("/")+lepton_[channel]+std::string("/")+samples[s]+std::string("/job*.root");
    std::cout << "opening root files from " << pathToRootFiles << std::endl;
    tree->Add(pathToRootFiles.c_str());
  }
  
  switch(channel) {
    case kElectron:
      eff = 0.003022831;
      break;
    case kMuon:
      eff = 0.003728336;
      break;
    case kAll:
      eff = 0.006751167;
      break;
    case kMuon_BReg:
      eff = 0.003728336;
      break;
  }
  
  TMultiGraph *mgMass = new TMultiGraph();
  mgMass->SetTitle(";m_{t,gen} [GeV];<m_{t,extr} #font[122]{\55} m_{t,gen}> [GeV]"); 

  TMultiGraph *mgJES = new TMultiGraph();
  mgJES->SetTitle(";m_{t,gen} [GeV];<JSF_{extr} #font[122]{\55} JSF>");
  
  if (calibrated) {
    mgMass->SetTitle(";m_{t,gen} [GeV];<m_{t,cal} #font[122]{\55} m_{t,gen}> [GeV]");
    mgJES->SetTitle(";m_{t,gen} [GeV];<JSF_{cal} #font[122]{\55} JSF>");
  }
  
  TMultiGraph *mgMassPull = new TMultiGraph();
  mgMassPull->SetTitle(";m_{t,gen} [GeV];Width of mass pull ");
  
  TMultiGraph *mgJESPull = new TMultiGraph();
  mgJESPull->SetTitle(";m_{t,gen} [GeV];Width of JSF pull ");
  
  TH2D* h2Mass = new TH2D("h2Mass", "h2Mass", 1000, 150, 200, 1000, 0.9, 1.1);
  TH2D* h2JES = new TH2D("h2JES", "h2JES", 1000, 150, 200, 1000, 0.9, 1.1);
  
  for (int iJES = 0; iJES < 3; iJES++) {
    for (int iMass = 0; iMass < 7; iMass++) {
      
      TString sel("mass_mTop_JES>0 & JES_mTop_JES>0 & genMass=="); sel+=genMass[iMass]; sel+=" & genJES=="; sel+=genJES[iJES];
      double entries = tree->GetEntries(sel);
      /*
      entries = nentries/27;
      if (iMass == 4) entries = entries/3;
      //*/ 
      
      TF1* gausMassBias = new TF1("gausMassBias", "gaus");
      tree->Fit("gausMassBias", "mass_mTop_JES", sel, "Q0", "", nentries, firstentry);
      //tree->Fit("gausMassBias", "mass+2.25341e-01+0.2/0.04*(1-JES)", sel, "Q0");
      
      mass[iJES][iMass]          = genMass[iMass];
      massBias[iJES][iMass]      = gausMassBias->GetParameter(1) - genMass[iMass];
      massBiasError[iJES][iMass] = gausMassBias->GetParameter(2) / sqrt(genMassN[iMass]/(crossSection*peLumi*maxMCWeight[iMass]));
        
      TF1* gausJESBias = new TF1("gausJESBias", "gaus");
      tree->Fit("gausJESBias", "JES_mTop_JES", sel, "Q0", "", nentries, firstentry);
      //tree->Fit("gausJESBias", "JES + 2.66486e-03 + 0.002/0.035*(JES-1)", sel, "Q0");
      
      JES[iJES][iMass]          = genJES[iJES];
      JESBias[iJES][iMass]      = gausJESBias->GetParameter(1) - genJES[iJES];
      JESBiasError[iJES][iMass] = gausJESBias->GetParameter(2) / sqrt(genMassN[iMass]/(crossSection*peLumi*maxMCWeight[iMass]));
      
      /* TEST
      if (iMass == 4) {
        massBiasError[iJES][iMass] = 1000.;
        JESBiasError [iJES][iMass] = 1000.;
      }
      //*/
      
      // Fill calibration histos
      //*
      h2Mass->Fill(gausMassBias->GetParameter(1), gausJESBias->GetParameter(1), massBias[iJES][iMass]);
      h2Mass->SetBinError(h2Mass->FindBin(gausMassBias->GetParameter(1), gausJESBias->GetParameter(1)), massBiasError[iJES][iMass]);
      
      h2JES->Fill(gausMassBias->GetParameter(1), gausJESBias->GetParameter(1), JESBias[iJES][iMass]);
      h2JES->SetBinError(h2JES->FindBin(gausMassBias->GetParameter(1), gausJESBias->GetParameter(1)), JESBiasError[iJES][iMass]);
      //*/
      
      TF1* gausMassPull = new TF1("gausMassPull", "gaus");
      tree->Fit("gausMassPull", "mass_mTop_JES_Pull", sel, "Q0", "", nentries, firstentry);
      
      massPull[iJES][iMass]      = gausMassPull->GetParameter(2);
      massPullError[iJES][iMass] = sqrt(1./2. * (maxMCWeight[iMass]/(genMassN[iMass]*eff) + 1./(entries-1.)));
      
      //std::cout << sqrt(1./2. * (1./(5144.*genMassN[iMass])+1./2000)) << std::endl;
        
      TF1* gausJESPull = new TF1("gausJESPull", "gaus");
      tree->Fit("gausJESPull", "JES_mTop_JES_Pull", sel, "Q0", "", nentries, firstentry);
      
      JESPull[iJES][iMass]      = gausJESPull->GetParameter(2);
      JESPullError[iJES][iMass] = sqrt(1./2. * (maxMCWeight[iMass]/(genMassN[iMass]*eff) + 1./(entries-1.)));
    }
    
    double genMassMod[7];
    for (int i = 0; i < 7; i++) {
      genMassMod[i] = genMass[i] + 0.2*(iJES-1);
    }

    gMass.push_back(new TGraphErrors(7, genMassMod, massBias[iJES], genMassError, massBiasError[iJES]));
    gMass[iJES]->SetMarkerStyle(marker_[iJES]);
    gMass[iJES]->SetMarkerColor(color_ [iJES]);
    gMass[iJES]->SetLineColor  (color_ [iJES]);
    mgMass->Add(gMass[iJES]);

    gJES.push_back(new TGraphErrors(7, genMassMod, JESBias[iJES], genMassError, JESBiasError[iJES]));
    gJES[iJES]->SetMarkerStyle(marker_[iJES]);
    gJES[iJES]->SetMarkerColor(color_ [iJES]);
    gJES[iJES]->SetLineColor  (color_ [iJES]);
    mgJES->Add(gJES[iJES]);

    gMassPull.push_back(new TGraphErrors(7, genMassMod, massPull[iJES], genMassError, massPullError[iJES]));
    gMassPull[iJES]->SetMarkerStyle(marker_[iJES]);
    gMassPull[iJES]->SetMarkerColor(color_ [iJES]);
    gMassPull[iJES]->SetLineColor  (color_ [iJES]);
    mgMassPull->Add(gMassPull[iJES]);

    gJESPull.push_back(new TGraphErrors(7, genMassMod, JESPull[iJES], genMassError, JESPullError[iJES]));
    gJESPull[iJES]->SetMarkerStyle(marker_[iJES]);
    gJESPull[iJES]->SetMarkerColor(color_ [iJES]);
    gJESPull[iJES]->SetLineColor  (color_ [iJES]);
    mgJESPull->Add(gJESPull[iJES]);
  }
  
  constFit = new TF1("constFit", "[0]");
  constFit->SetParNames("offset");
  constFit->SetLineColor(kBlack);
  constFit->SetLineWidth(2);
  constFit->SetLineStyle(1);
  
  linearFit096 = new TF1("linearFit096", "[0]+[1]*(x-172.5)");
  linearFit096->SetParNames("offset", "slope");
  linearFit096->SetLineColor(color_[0]);
  linearFit096->SetLineWidth(2);
  linearFit096->SetLineStyle(9);
  linearFit096->SetMarkerStyle(marker_[0]);
  linearFit096->SetMarkerColor(color_ [0]);
  
  linearFit100 = new TF1("linearFit100", "[0]+[1]*(x-172.5)");
  linearFit100->SetParNames("offset", "slope");
  linearFit100->SetLineColor(color_[1]);
  linearFit100->SetLineWidth(3);
  linearFit100->SetLineStyle(2);
  linearFit100->SetMarkerStyle(marker_[1]);
  linearFit100->SetMarkerColor(color_ [1]);
  
  linearFit104 = new TF1("linearFit104", "[0]+[1]*(x-172.5)");
  linearFit104->SetParNames("offset", "slope");
  linearFit104->SetLineColor(color_[2]);
  linearFit104->SetLineWidth(2);
  linearFit104->SetLineStyle(5);
  linearFit104->SetMarkerStyle(marker_[2]);
  linearFit104->SetMarkerColor(color_ [2]);
  
  linearFitJES = new TF1("linearFitJES", "[0]+[1]*(x-1.0)");
  linearFitJES->SetParNames("offset", "slope");


  canvasFit->Clear();
  canvasFit->Divide(1, 2, -1., -1.);
  canvasFit->cd(1);
  
  std::cout << "MASS, fit constant" << std::endl;
  fitResult = mgMass->Fit("constFit", "EMS");
  std::cout << "Fit linear, JES 0.98" << std::endl;
  gMass[0]->Fit("linearFit096", "EM");
  std::cout << "Fit linear" << std::endl;
  gMass[1]->Fit("linearFit100", "EM");
  std::cout << "Fit linear, JES 1.02" << std::endl;
  gMass[2]->Fit("linearFit104", "EM");
  mgMass->SetMinimum(-1.9);
  mgMass->SetMaximum( 1.9);
  if (channel==2) {
    mgMass->SetMinimum(-0.75);
    mgMass->SetMaximum( 0.75);
  }
  mgMass->Draw("AP");
  gPad->SetBottomMargin(0.005);
  mgMass->GetYaxis()->SetTitleSize(0.11);
  mgMass->GetYaxis()->SetLabelSize(0.10);
  mgMass->GetYaxis()->SetTitleOffset(0.82);
  
  /*
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
  */
  
  canvasFit->cd(2);
  
  std::cout << "JES, fit constant" << std::endl;
  fitResult = mgJES->Fit("constFit", "EMS");
  std::cout << "Fit linear, JES 0.98" << std::endl;
  gJES[0]->Fit("linearFit096", "EM");
  std::cout << "Fit linear" << std::endl;
  gJES[1]->Fit("linearFit100", "EM");
  std::cout << "Fit linear, JES 1.02" << std::endl;
  gJES[2]->Fit("linearFit104", "EM");
  mgJES->SetMinimum(-0.019);
  mgJES->SetMaximum( 0.019);
  if (channel==2) {
    mgJES->SetMinimum(-0.0115);
    mgJES->SetMaximum( 0.0115);
  }
  mgJES->Draw("AP");
  gPad->SetTopMargin(0.005);
  //canvasFit->Update();
  
  canvasFit->cd(0);
  CMS_lumi(canvasFit, 210, 0., true, "Simulation"); // lepton+jets
  DrawLegend(true);
  
  TString sFitMass("fit_Bias_"); sFitMass += lepton_[channel]; sFitMass += suffix; sFitMass += ".eps";
  canvasFit->Print(sFitMass);
  
  /*
  offset[0] = linearFit096->GetParameter(0);
  offset[1] = linearFit100->GetParameter(0);
  offset[2] = linearFit104->GetParameter(0);
  offsetError[0] = linearFit096->GetParError(0);
  offsetError[1] = linearFit100->GetParError(0);
  offsetError[2] = linearFit104->GetParError(0);

  gOffset = new TGraphErrors(3, measJES, offset, genJESError, offsetError);
  gOffset->Fit("linearFitJES", "EM");
  gOffset->Draw("AP");
  */
  
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
  
  TCanvas* canvasFitPull = new TCanvas("canvasFitPull", "mt-JES measurement pull width", 500, 500);
  canvasFitPull->Divide(1, 2, 0, 0);
  canvasFitPull->cd(1);
  
  std::cout << "MASS PULL" << std::endl;
  fitResult = mgMassPull->Fit("constFit", "EMS");
  gMassPull[0]->Fit("linearFit096", "EM0");
  gMassPull[1]->Fit("linearFit100", "EM0");
  gMassPull[2]->Fit("linearFit104", "EM0");
  mgMassPull->SetMinimum(0.81);
  mgMassPull->SetMaximum(1.19);
  mgMassPull->Draw("AP");
  gPad->SetBottomMargin(0.005);
  mgMassPull->GetYaxis()->SetTitleSize(0.11);
  mgMassPull->GetYaxis()->SetLabelSize(0.10);
  mgMassPull->GetYaxis()->SetTitleOffset(0.82);
  
  canvasFitPull->cd(2);
  std::cout << "JES PULL" << std::endl;
  fitResult = mgJESPull->Fit("constFit", "EMS");
  gJESPull[0]->Fit("linearFit096", "EM0");
  gJESPull[1]->Fit("linearFit100", "EM0");
  gJESPull[2]->Fit("linearFit104", "EM0");
  mgJESPull->SetMinimum(0.81);
  mgJESPull->SetMaximum(1.19);
  mgJESPull->Draw("AP");
  gPad->SetTopMargin(0.005);
  
  canvasFitPull->cd(0);
  CMS_lumi(canvasFitPull, 210, 0., true, "Simulation"); // lepton+jets
  DrawLegend();
  
  TString sFitJESPull("fit_Pull_"); sFitJESPull += lepton_[channel]; sFitJESPull += suffix; sFitJESPull += ".eps";
  canvasFitPull->Print(sFitJESPull);
  //*/
}


