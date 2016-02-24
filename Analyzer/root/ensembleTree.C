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

//const int nJES = 5;
//const int nMasses = 9;
//enum styles          { kDown1, kDown2, kNominal, kUp1, kUp2 };
//int color_ [] = { kRed+1, kMagenta, kBlue+1, kCyan, kGreen+1 };
//int marker_[] = { 23, 23, 20, 22, 22 };
//double genJES[]       = {0.96, 0.98, 1.00, 1.02, 1.04};
////double genJESError[]  = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
//double genMass[]      = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5, 184.5};
//double genMassError[] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
//double genMassN[]     = {5000000, 5000000, 5000000, 5000000, 90000000, 5000000, 5000000, 5000000, 5000000}; //{1620072, 1633197, 1669034, 1606570, 59613991, 1538301, 1648519, 1665350, 1671859};
//double maxMCWeight[]  = {1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8};
//
////double genMassN[]     = {1620072, 1.5, 1.5, 1.5, 59613991, 1.5, 1.5, 1.5, 1.5};
//double crossSection   = 247.0;
//double peLumi         = 18352.0;

// const int nJES = 5;
// const int nMasses = 7;
// enum styles          { kDown1, kDown2, kNominal, kUp1, kUp2 };
// int color_ [] = { kRed+1, kMagenta, kBlue+1, kCyan, kGreen+1 };
// int marker_[] = { 23, 23, 20, 22, 22 };
// double genJES[]       = {0.96, 0.98, 1.00, 1.02, 1.04};
// //double genJESError[]  = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
// double genMass[]      = {166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5};
// double genMassError[] = { 1e-6,  1e-6,  1e-6,  1e-6,  1e-6,  1e-6,  1e-6};
// double genMassN[]     = {27078777, 39518234, 24439341, 62131965, 26489020, 40244328, 24359161};
// double maxMCWeight[]  = {1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8};

const int nJES = 3;
const int nMasses = 3;
int color_ [] = { kRed+1, kBlue+1, kGreen+1 ,kOrange+1, kViolet+1};
int marker_[] = { 23, 20, 22 ,24, 25  };
double genJES[]       = {/*0.96,*/ 0.98, 1.00, 1.02/*, 1.04*/};
double genJESError[]  = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
double genMass[]      = { 169.5, 172.5,  175.5};
double genMassError[] = { 1e-6,  1e-6,  1e-6,  1e-6,  1e-6,  1e-6,  1e-6};
double genMassN[]     = {9969800, 19757190, 9904200};
double maxMCWeight[]  = {1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8};

double crossSection   = 831.7;
double peLumi         = 2192.;
  
double measJES[nJES];

double offset[nJES];
double offsetError[nJES];
  
double mass[nJES][nMasses];
double massBias[nJES][nMasses];
double massBiasError[nJES][nMasses];

double massPull[nJES][nMasses];
double massPullError[nJES][nMasses];
  
double JES[nJES][nMasses];
double JESBias[nJES][nMasses];
double JESBiasError[nJES][nMasses];

double JESPull[nJES][nMasses];
double JESPullError[nJES][nMasses];

std::vector<TGraphErrors*> gMass;
std::vector<TGraphErrors*> gJES;

std::vector<TGraphErrors*> gMassPull;
std::vector<TGraphErrors*> gJESPull;

TTree* tree;

TF1* topMass1DUncertaintyFit;
TF1* topMass2DUncertaintyFit;
TF1* jesUncertaintyFit;
TF1* constFit;
TF1* linearFit096;
TF1* linearFit098;
TF1* linearFit100;
TF1* linearFit102;
TF1* linearFit104;
TF1* linearFitJES;
TFitResultPtr fitResult;

/*
void DrawLegend() {
  TLegend *leg = new TLegend(0.74, 0.7, 0.97, 0.94);
  leg->SetTextSizePixels(12);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  if(nJES  == 5){
    leg->AddEntry( gMass[0], "JES=0.96", "EP");
    leg->AddEntry( gMass[1], "JES=0.98", "EP");
    leg->AddEntry( gMass[2], "JES=1.00", "EP");
    leg->AddEntry( gMass[3], "JES=1.02", "EP");
    leg->AddEntry( gMass[4], "JES=1.04", "EP");
    //leg->AddEntry( constFit, "Const. fit", "L");
  }
  else{
    leg->AddEntry( gMass[0], "JES=0.96", "EP");
    leg->AddEntry( gMass[1], "JES=1.00", "EP");
    leg->AddEntry( gMass[2], "JES=1.04", "EP");
    //leg->AddEntry( constFit, "Const. fit", "L");
  }
  //char chi2[6]; sprintf(chi2, "%3.1f", fitResult->Chi2());
  //TString sConstFit("#chi^{2}/ndf="); sConstFit+=chi2; sConstFit+="/"; sConstFit+=fitResult->Ndf();
  //leg->AddEntry((TObject*)0, sConstFit, "");
  leg->Draw();
}
*/
void DrawLegend() {
  double y1 = 0.915; double y2 = 0.955;
  int textsize = 7;
  
  TLegend *leg1 = new TLegend(0.170, y1, 0.334, y2);
  leg1->SetTextSizePixels(textsize);
  leg1->SetFillColor(kWhite);
  leg1->SetBorderSize(0);
  if(nJES==3) leg1->AddEntry(gMass[0], "JSF=0.98", "EP");
  if(nJES==5) leg1->AddEntry(gMass[0], "JSF=0.98", "EP");
  /*
  leg->AddEntry( constFit, "Const. fit", "L");
  char chi2[6]; sprintf(chi2, "%3.1f", fitResult->Chi2());
  TString sConstFit("#chi^{2}/ndf="); sConstFit+=chi2; sConstFit+="/"; sConstFit+=fitResult->Ndf();
  leg->AddEntry((TObject*)0, sConstFit, "");
  */
  leg1->Draw();
  
  if(nJES==5){
    TLegend *leg2 = new TLegend(0.334, y1, 0.498, y2);
    leg2->SetTextSizePixels(textsize);
    leg2->SetFillColor(kWhite);
    leg2->SetBorderSize(0);
    leg2->AddEntry(gMass[1], "JSF=0.98", "EP");
    leg2->Draw();
  }

  TLegend *leg3 = new TLegend(0.498, y1, 0.662, y2);
  leg3->SetTextSizePixels(textsize);
  leg3->SetFillColor(kWhite);
  leg3->SetBorderSize(0);
  leg3->AddEntry(gMass[nJES==5?2:1], "JSF=1.00", "EP");
  leg3->Draw();
  
  if(nJES==5){
    TLegend *leg4 = new TLegend(0.662, y1, 0.826, y2);
    leg4->SetTextSizePixels(textsize);
    leg4->SetFillColor(kWhite);
    leg4->SetBorderSize(0);
    leg4->AddEntry(gMass[3], "JSF=1.02", "EP");
    leg4->Draw();
  }
  
  TLegend *leg5 = new TLegend(0.826, y1, 0.990, y2);
  leg5->SetTextSizePixels(textsize);
  leg5->SetFillColor(kWhite);
  leg5->SetBorderSize(0);
  if(nJES==3) leg5->AddEntry(gMass[2], "JSF=1.02", "EP");
  if(nJES==5) leg5->AddEntry(gMass[4], "JSF=1.02", "EP");
  leg5->Draw();
}

void ensembleTree()
{
  //*
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
  
  TCanvas* canvasFit = new TCanvas("canvasFit", "mt-JSF measurement calibration", 500, 500);
  canvasFit->cd();

  TString sFile("/nfs/dust/cms/user/garbersc/TopMass/2015_TemplateCalibrationMerged/"); 
  sFile += "ensemble_2015D3JES_electron_Calibrated.root";
  //sFile += "ensemble_S12_Calibrated.root";

  std::cout << "Doing calibration on: " << sFile << std::endl;
  TFile* fEnsemble = new TFile(sFile);
  
  tree = (TTree*) fEnsemble->Get("tree");
  
  std::cout << "Entries: " << tree->GetEntriesFast() << std::endl;

  TMultiGraph *mgMass = new TMultiGraph();
  if(sFile.Contains("Uncalibrated")) mgMass->SetTitle(";m_{t,gen} [GeV];<m_{t,extr}-m_{t,gen}> [GeV]");
  else                               mgMass->SetTitle(";m_{t,gen} [GeV];<m_{t,cal}-m_{t,gen}> [GeV]");
  
  TMultiGraph *mgJES = new TMultiGraph();
  if(sFile.Contains("Uncalibrated")) mgJES->SetTitle(";m_{t,gen} [GeV];<JSF_{extr}-JSF>");
  else                               mgJES->SetTitle(";m_{t,gen} [GeV];<JSF_{cal}-JSF>");
  
  TMultiGraph *mgMassPull = new TMultiGraph();
  mgMassPull->SetTitle(";m_{t,gen} [GeV];Mass pull width");
  
  TMultiGraph *mgJESPull = new TMultiGraph();
  mgJESPull->SetTitle(";m_{t,gen} [GeV];JSF pull width");
  
  TH2D* h2Mass = new TH2D("h2Mass", "h2Mass", 1000, 150, 200, 1000, 0.9, 1.1);
  TH2D* h2JES = new TH2D("h2JES", "h2JES", 1000, 150, 200, 1000, 0.9, 1.1);
  


  for (int iJES = 0; iJES < nJES; iJES++) {
    for (int iMass = 0; iMass < nMasses; ++iMass) {
      //TString sel("mass_mTop_JES_fSig_fCP>0 & JES_mTop_JES_fSig_fCP>0 & JES_mTop_JES_fSig_fCP_Error>0.00 & JES_mTop_JES_fSig_fCP_Error<10.0 & mass_mTop_JES_fSig_fCP_Error>0.5 & mass_mTop_JES_fSig_fCP_Error<2.5 & abs(mass_mTop_JES_fSig_fCP_Pull)<100 & genMass=="); sel+=genMass[iMass]; sel+=" & genJES=="; sel+=genJES[iJES];
      TString sel("mass_mTop_JES>0 & JES_mTop_JES>0 & genMass=="); sel+=genMass[iMass]; sel+=" & genJES=="; sel+=genJES[iJES];
      //TString sel("mass_mTop_JES_fSig_fCP>0 & JES_mTop_JES_fSig_fCP>0 & fSig_mTop_JES_fSig_fCP>0 & fCP_mTop_JES_fSig_fCP>0 & genMass=="); sel+=genMass[iMass]; sel+=" & genJES=="; sel+=genJES[iJES];
      //TString sel("mass_mTop_JES_fSig>0 & JES_mTop_JES_fSig>0 & fSig_mTop_JES_fSig>0 & genMass=="); sel+=genMass[iMass]; sel+=" & genJES=="; sel+=genJES[iJES];
      double entries = tree->GetEntries(sel);
      
      TF1* gausMassBias = new TF1("gausMassBias", "gaus");
     // tree->Fit("gausMassBias", "mass_mTop_JES_fSig_fCP", sel, "Q0");
      //tree->Fit("gausMassBias", "mass_mTop_JES_fSig", sel, "Q0");
	tree->Fit("gausMassBias", "mass_mTop_JES", sel, "Q0");
      
      mass[iJES][iMass]          = genMass[iMass];
      massBias[iJES][iMass]      = gausMassBias->GetParameter(1) - genMass[iMass];
      massBiasError[iJES][iMass] = gausMassBias->GetParameter(2) / sqrt(genMassN[iMass]/(crossSection*peLumi*maxMCWeight[iMass]));
        
      TF1* gausJESBias = new TF1("gausJESBias", "gaus");
   //   tree->Fit("gausJESBias", "JES_mTop_JES_fSig_fCP", sel, "Q0");
      tree->Fit("gausJESBias", "JES_mTop_JES", sel, "Q0");
      //tree->Fit("gausJESBias", "JES_mTop_JES_fSig", sel, "Q0");
      
      JES[iJES][iMass]          = genJES[iJES];
      JESBias[iJES][iMass]      = gausJESBias->GetParameter(1) - genJES[iJES];
      JESBiasError[iJES][iMass] = gausJESBias->GetParameter(2) / sqrt(genMassN[iMass]/(crossSection*peLumi*maxMCWeight[iMass]));
      
      // Fill calibration histos
      
      h2Mass->Fill(gausMassBias->GetParameter(1), gausJESBias->GetParameter(1), massBias[iJES][iMass]);
      h2Mass->SetBinError(h2Mass->FindBin(gausMassBias->GetParameter(1), gausJESBias->GetParameter(1)), massBiasError[iJES][iMass]);
      
      h2JES->Fill(gausMassBias->GetParameter(1), gausJESBias->GetParameter(1), JESBias[iJES][iMass]);
      h2JES->SetBinError(h2JES->FindBin(gausMassBias->GetParameter(1), gausJESBias->GetParameter(1)), JESBiasError[iJES][iMass]);
      
      TF1* gausMassPull = new TF1("gausMassPull", "gaus");
     // tree->Fit("gausMassPull", "mass_mTop_JES_fSig_fCP_Pull", sel, "Q0");
      tree->Fit("gausMassPull", "mass_mTop_JES_Pull", sel, "Q0");
      //tree->Fit("gausMassPull", "mass_mTop_JES_fSig_Pull", sel, "Q0");
 
      double eff = 0.000838398;
      massPull[iJES][iMass]      = gausMassPull->GetParameter(2);
      massPullError[iJES][iMass] = sqrt(1./2. * (maxMCWeight[iMass]/(genMassN[iMass]*eff) + 1./(entries-1.)));
      
      TF1* gausJESPull = new TF1("gausJESPull", "gaus");
     // tree->Fit("gausJESPull", "JES_mTop_JES_fSig_fCP_Pull", sel, "Q0");
      tree->Fit("gausJESPull", "JES_mTop_JES_Pull", sel, "Q0");
      //tree->Fit("gausJESPull", "JES_mTop_JES_fSig_Pull", sel, "Q0");
      
      JESPull[iJES][iMass]      = gausJESPull->GetParameter(2);
      JESPullError[iJES][iMass] = sqrt(1./2. * maxMCWeight[iMass]/(genMassN[iMass]*eff) + 1./(entries-1.));
    }
    
    double genMassMod[nMasses];
    for (int i = 0; i < nMasses; i++) {
      genMassMod[i] = genMass[i] + 0.2*(iJES-1);
    }

    gMass.push_back(new TGraphErrors(nMasses, genMassMod, massBias[iJES], genMassError, massBiasError[iJES]));
    gMass[iJES]->SetMarkerStyle(marker_[iJES]);
    gMass[iJES]->SetMarkerColor(color_ [iJES]);
    gMass[iJES]->SetLineColor  (color_ [iJES]);
    mgMass->Add(gMass[iJES]);

    gJES.push_back(new TGraphErrors(nMasses, genMassMod, JESBias[iJES], genMassError, JESBiasError[iJES]));
    gJES[iJES]->SetMarkerStyle(marker_[iJES]);
    gJES[iJES]->SetMarkerColor(color_ [iJES]);
    gJES[iJES]->SetLineColor  (color_ [iJES]);
    mgJES->Add(gJES[iJES]);

    gMassPull.push_back(new TGraphErrors(nMasses, genMassMod, massPull[iJES], genMassError, massPullError[iJES]));
    gMassPull[iJES]->SetMarkerStyle(marker_[iJES]);
    gMassPull[iJES]->SetMarkerColor(color_ [iJES]);
    gMassPull[iJES]->SetLineColor  (color_ [iJES]);
    mgMassPull->Add(gMassPull[iJES]);

    gJESPull.push_back(new TGraphErrors(nMasses, genMassMod, JESPull[iJES], genMassError, JESPullError[iJES]));
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
  
  topMass1DUncertaintyFit= new TF1("topMass1DUncertaintyFit", "[0]+[1]*(x-172.5)");
  topMass1DUncertaintyFit->SetParNames("offset1D", "slope1D");
  topMass1DUncertaintyFit->SetLineColor(kBlack);
  topMass1DUncertaintyFit->SetLineWidth(3);
  topMass1DUncertaintyFit->SetLineStyle(7);

  topMass2DUncertaintyFit= new TF1("topMass2DUncertaintyFit", "[0]+[1]*(x-172.5)");
  topMass2DUncertaintyFit->SetParNames("offset2D", "slope2D");
  topMass2DUncertaintyFit->SetLineColor(kBlack);
  topMass2DUncertaintyFit->SetLineWidth(3);
  topMass2DUncertaintyFit->SetLineStyle(7);

  jesUncertaintyFit= new TF1("jesUncertaintyFit", "[0]+[1]*(x-172.5)");
  jesUncertaintyFit->SetParNames("offsetJES", "slopeJES");
  jesUncertaintyFit->SetLineColor(kBlack);
  jesUncertaintyFit->SetLineWidth(3);
  jesUncertaintyFit->SetLineStyle(7);

  linearFit096 = new TF1("linearFit096", "[0]+[1]*(x-172.5)");
  linearFit096->SetParNames("offset", "slope");
  if(nJES==5) linearFit096->SetLineColor(color_[0]);
  else        linearFit096->SetLineColor(color_[0]);
  linearFit096->SetLineWidth(1);
  linearFit096->SetLineStyle(7);
  
  linearFit098 = new TF1("linearFit098", "[0]+[1]*(x-172.5)");
  linearFit098->SetParNames("offset", "slope");
  if(nJES==5) linearFit098->SetLineColor(color_[1]);
  else        linearFit098->SetLineColor(color_[0]);
  linearFit098->SetLineWidth(1);
  linearFit098->SetLineStyle(7);
  
  linearFit100 = new TF1("linearFit100", "[0]+[1]*(x-172.5)");
  linearFit100->SetParNames("offset", "slope");
  if(nJES==5) linearFit100->SetLineColor(color_[2]);
  else        linearFit100->SetLineColor(color_[1]);
  linearFit100->SetLineWidth(1);
  linearFit100->SetLineStyle(7);
  
  linearFit102 = new TF1("linearFit102", "[0]+[1]*(x-172.5)");
  linearFit102->SetParNames("offset", "slope");
  if(nJES==5) linearFit102->SetLineColor(color_[3]);
  else        linearFit102->SetLineColor(color_[2]);
  linearFit102->SetLineWidth(1);
  linearFit102->SetLineStyle(7);

  linearFit104 = new TF1("linearFit104", "[0]+[1]*(x-172.5)");
  linearFit104->SetParNames("offset", "slope");
  if(nJES==5) linearFit104->SetLineColor(color_[4]);
  else        linearFit104->SetLineColor(color_[2]);
  linearFit104->SetLineWidth(1);
  linearFit104->SetLineStyle(7);

  linearFitJES = new TF1("linearFitJES", "[0]+[1]*(x-1.0)");
  linearFitJES->SetParNames("offset", "slope");

  canvasFit->Clear();
  canvasFit->Divide(1, 2, -1., -1.);
  //canvasFit->Divide(1, 2, 0, 0);
  canvasFit->cd(1);
  
  std::cout << "=== calibration - mass ===" << std::endl;
  //canvasFit->Clear();
  //fitResult = mgMass->Fit("constFit", "EMS");
  mgMass->Fit("topMass2DUncertaintyFit", "EM");
  if(nJES  == 5){
    gMass[0]->Fit("linearFit096", "EM");
    gMass[1]->Fit("linearFit098", "EM");
    gMass[2]->Fit("linearFit100", "EM");
    gMass[3]->Fit("linearFit102", "EM");
    gMass[4]->Fit("linearFit104", "EM");

    gMass[2]->Fit("topMass1DUncertaintyFit", "EM0");
  }
  else{
    gMass[0]->Fit("linearFit098", "EM");
    gMass[1]->Fit("linearFit100", "EM");
    gMass[2]->Fit("linearFit102", "EM");

    gMass[1]->Fit("topMass1DUncertaintyFit", "EM0");
  }

   mgMass->SetMinimum(-1.58);
    mgMass->SetMaximum( 1.58);

 /* if(sFile.Contains("Calibrated")){
    mgMass->SetMinimum(-1.06);
    mgMass->SetMaximum( 1.06);
  } else {
    mgMass->SetMinimum(-2.15);
    mgMass->SetMaximum( 2.15);
  }*/
  mgMass->Draw("AP");
  gPad->SetBottomMargin(0.005);
  mgMass->GetYaxis()->SetTitleSize(0.11);
  mgMass->GetYaxis()->SetLabelSize(0.10);
  mgMass->GetYaxis()->SetTitleOffset(0.82);
  //canvasFit->Update();
  
  /*
  DrawLegend();
  DrawCMSSim();
  canvasFit->Print("fit_Mass.eps");
  
  canvasFit->Clear();
  */
  if(nJES  == 5){
    offset[0] = linearFit096->GetParameter(0);
    offset[1] = linearFit098->GetParameter(0);
    offset[2] = linearFit100->GetParameter(0);
    offset[3] = linearFit102->GetParameter(0);
    offset[4] = linearFit104->GetParameter(0);
    offsetError[0] = linearFit096->GetParError(0);
    offsetError[1] = linearFit098->GetParError(0);
    offsetError[2] = linearFit100->GetParError(0);
    offsetError[3] = linearFit102->GetParError(0);
    offsetError[4] = linearFit104->GetParError(0);
  }
  else{
    offset[0] = linearFit098->GetParameter(0);
    offset[1] = linearFit100->GetParameter(0);
    offset[2] = linearFit102->GetParameter(0);
    offsetError[0] = linearFit098->GetParError(0);
    offsetError[1] = linearFit100->GetParError(0);
    offsetError[2] = linearFit102->GetParError(0);
  }

  //gOffset = new TGraphErrors(nJES, measJES, offset, genJESError, offsetError);
  //gOffset->Fit("linearFitJES", "EM");
  //gOffset->Draw("AP");
  

  canvasFit->cd(2);

  std::cout << "=== calibration - JES ===" << std::endl;
  fitResult = mgJES->Fit("constFit", "EMS");
  mgJES->Fit("jesUncertaintyFit", "EM");
  if(nJES  == 5){
    gJES[0]->Fit("linearFit096", "EM");
    gJES[1]->Fit("linearFit098", "EM");
    gJES[2]->Fit("linearFit100", "EM");
    gJES[3]->Fit("linearFit102", "EM");
    gJES[4]->Fit("linearFit104", "EM");
  }
  else{
    gJES[0]->Fit("linearFit098", "EM");
    gJES[1]->Fit("linearFit100", "EM");
    gJES[2]->Fit("linearFit102", "EM");
  }

    mgJES->SetMinimum(-0.0158);
    mgJES->SetMaximum( 0.0158);
 /* if(sFile.Contains("Calibrated")){
    mgJES->SetMinimum(-0.0106);
    mgJES->SetMaximum( 0.0106);
  } else {
    mgJES->SetMinimum(-0.0215);
    mgJES->SetMaximum( 0.0215);
  }*/
  mgJES->Draw("AP");
  gPad->SetTopMargin(0.005);
  //DrawLegend();
  //DrawCMSSim();
  //canvasFit->Print("fit_JES.eps");

  //canvasFit->Clear();

  canvasFit->cd(0);
  DrawLabel("Private Work,  #sqrt{s}=13 TeV", 0.2, 0.93, 0.9);
  DrawLegend();

  TString sFitMass("fit_Bias.eps");
  canvasFit->Print(sFitMass);

  
  if(nJES  == 5){
    offset[0] = linearFit096->GetParameter(0);
    offset[1] = linearFit098->GetParameter(0);
    offset[2] = linearFit100->GetParameter(0);
    offset[3] = linearFit102->GetParameter(0);
    offset[4] = linearFit104->GetParameter(0);
    offsetError[0] = linearFit096->GetParError(0);
    offsetError[1] = linearFit098->GetParError(0);
    offsetError[2] = linearFit100->GetParError(0);
    offsetError[3] = linearFit102->GetParError(0);
    offsetError[4] = linearFit104->GetParError(0);
    measJES[0] = linearFit096->GetParameter(0)+genJES[0];
    measJES[1] = linearFit098->GetParameter(0)+genJES[1];
    measJES[2] = linearFit100->GetParameter(0)+genJES[2];
    measJES[3] = linearFit102->GetParameter(0)+genJES[3];
    measJES[4] = linearFit104->GetParameter(0)+genJES[4];
  }
  else{
    offset[0] = linearFit098->GetParameter(0);
    offset[1] = linearFit100->GetParameter(0);
    offset[2] = linearFit102->GetParameter(0);
    offsetError[0] = linearFit098->GetParError(0);
    offsetError[1] = linearFit100->GetParError(0);
    offsetError[2] = linearFit102->GetParError(0);
    measJES[0] = linearFit098->GetParameter(0)+genJES[0];
    measJES[1] = linearFit100->GetParameter(0)+genJES[1];
    measJES[2] = linearFit102->GetParameter(0)+genJES[2];
  }
  
  //TGraphErrors* gOffset = new TGraphErrors(nJES, measJES, offset, genJESError, offsetError);
  //gOffset->Fit("linearFitJES", "EM");
  //gOffset->Draw("AP");
  
  
  TF2* fit2D = new TF2("fit2D", "[0] + [1]*(x-172.5) + [2]*(y-1) + [3]*(x-172.5)*(y-1)");
  fit2D->SetParNames("offset", "slopeMass", "slopeJES","slopeMassJES");
  
  //*
  std::cout << "=== 2D calibration - mass ===" << std::endl;
  h2Mass->Fit("fit2D");
  //h2Mass->Draw("SURF1");
  
  std::cout << "=== 2D calibration - JES ===" << std::endl;
  h2JES->Fit("fit2D");
  //h2JES->Draw("SURF1");
  //*/
  
  //*

  TCanvas* canvasFitPull = new TCanvas("canvasFitPull", "mt-JES measurement pull width", 500, 500);
  canvasFitPull->Divide(1, 2, 0, 0);
  canvasFitPull->cd(1);

  std::cout << "=== calibration - massPull ===" << std::endl;
  //canvasFit->Clear();
  fitResult = mgMassPull->Fit("constFit", "EMS");
  if(nJES  == 5){
    gMassPull[0]->Fit("linearFit096", "EM0");
    gMassPull[1]->Fit("linearFit098", "EM0");
    gMassPull[2]->Fit("linearFit100", "EM0");
    gMassPull[3]->Fit("linearFit102", "EM0");
    gMassPull[4]->Fit("linearFit104", "EM0");
  }
  else{
    gMassPull[0]->Fit("linearFit098", "EM0");
    gMassPull[1]->Fit("linearFit100", "EM0");
    gMassPull[2]->Fit("linearFit102", "EM0");
  }
  //mgMassPull->SetMinimum(0.78);
  //mgMassPull->SetMaximum(1.22);
  //mgMassPull->SetMinimum(0.96);
  //mgMassPull->SetMaximum(1.39);
  mgMassPull->SetMinimum(0.89);
  mgMassPull->SetMaximum(1.11);
  mgMassPull->Draw("AP");
  gPad->SetBottomMargin(0.005);
  mgMassPull->GetYaxis()->SetTitleSize(0.11);
  mgMassPull->GetYaxis()->SetLabelSize(0.10);
  mgMassPull->GetYaxis()->SetTitleOffset(0.82);
  //DrawLegend();
  //DrawCMSSim();
  //canvasFit->Print("fit_MassPull.eps");
  
  std::cout << "=== calibration - JESPull ===" << std::endl;
  canvasFitPull->cd(2);
  //canvasFit->Clear();
  fitResult = mgJESPull->Fit("constFit", "EMS");
  if(nJES  == 5){
    gJESPull[0]->Fit("linearFit096", "EM0");
    gJESPull[1]->Fit("linearFit098", "EM0");
    gJESPull[2]->Fit("linearFit100", "EM0");
    gJESPull[3]->Fit("linearFit102", "EM0");
    gJESPull[4]->Fit("linearFit104", "EM0");
  }
  else{
    gJESPull[0]->Fit("linearFit098", "EM0");
    gJESPull[1]->Fit("linearFit100", "EM0");
    gJESPull[2]->Fit("linearFit102", "EM0");
  }
  //mgJESPull->SetMinimum(0.68);
  //mgJESPull->SetMaximum(1.12);
  //mgJESPull->SetMinimum(0.88);
  //mgJESPull->SetMaximum(1.38);
  mgJESPull->SetMinimum(0.89);
  mgJESPull->SetMaximum(1.11);
  mgJESPull->Draw("AP");
  gPad->SetTopMargin(0.005);
  //DrawLegend();
  //DrawCMSSim();
  //canvasFit->Print("fit_JESPull.eps");
  //*/

  canvasFitPull->cd(0);
  //DrawCMSSim();
  DrawLabel("Private Work,  #sqrt{s}=13 TeV", 0.2, 0.93, 0.9);
  DrawLegend();

  TString sFitJESPull("fit_Pull.eps");
  canvasFitPull->Print(sFitJESPull);
  
  std::cout << "=== 1D calibration uncertainty (mass / pull) ===" << std::endl;
  if(nJES  == 5){
    gMass[2]    ->Fit("constFit", "EM");
    gMassPull[2]->Fit("constFit", "EM");
  }
  else{
    gMass[1]    ->Fit("constFit", "EM");
    gMassPull[1]->Fit("constFit", "EM");
  }
  std::cout << "topmass 1D calibration uncertainty: " << std::sqrt(     pow(topMass1DUncertaintyFit->GetParError(0),2) +      pow((172.30-172.50)*topMass1DUncertaintyFit->GetParError(1),2)) << std::endl;
  std::cout << "topmass 2D calibration uncertainty: " << std::sqrt(nJES*pow(topMass2DUncertaintyFit->GetParError(0),2) + nJES*pow((171.88-172.50)*topMass2DUncertaintyFit->GetParError(1),2)) << std::endl;
  std::cout << "JES        calibration uncertainty: " << std::sqrt(nJES*pow(      jesUncertaintyFit->GetParError(0),2) + nJES*pow((171.88-172.50)*      jesUncertaintyFit->GetParError(1),2)) << std::endl;
}

