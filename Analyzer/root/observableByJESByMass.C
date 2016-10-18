#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TFitResult.h"
#include "TLine.h"
#include "TPaveStats.h"

//#include "tdrstyle_new.C"
#include "tdrstyle.C"
#include "CMS_lumi.C"

const int nMass = 5;
const int nJES  = 5;

int target = 1; // 1: correct, 0: wrong, -10: unmatched
int obs    = 0; // 0: hadTopMass, 1: hadWRawMass
int lepton = 1;// 0: electron, 1: muon, 2: lepton

int iMassMin = 0;
int iMassMax = nMass;

// bool plotByMass = false; // ??
// bool pas = false;

bool plotByMass = true;
bool pas = true;

int iTarget[]     = {1, 0, -10};
TString sTarget[] = {"wp", "cp", "un"};
TString sFObs[]   = {"mt", "mW"};
TString sObs[]    = {"m_{t}", "m_{W}^{reco}"};
TString sLepton[] = {"electron", "muon", "lepton"};

TGraphErrors* gr1[nMass];
TGraphErrors* gr2[nMass];
TGraphErrors* gr3[nMass];
TGraphErrors* gr4[nMass];

TString sX[nMass] = {"m_{t,gen} = 166.5 GeV",
                     "m_{t,gen} = 169.5 GeV",
                     //"m_{t,gen} = 171.5 GeV",
                     "m_{t,gen} = 172.5 GeV",
                     //"m_{t,gen} = 173.5 GeV",
                     "m_{t,gen} = 175.5 GeV",
                     "m_{t,gen} = 178.5 GeV"};

double X  [nMass] = {166.5, 169.5, /*171.5,*/ 172.5, /*173.5,*/ 175.5, 178.5};
TString X_ [nMass] = {"166_5", "169_5", /*"171_5",*/ "172_5", /*"173_5",*/ "175_5", "178_5"};
double Y10[nMass];
double Y11[nMass];
double Y20[nMass];
double Y21[nMass];
double Y30[nMass];
double Y31[nMass];

double eX  [nMass] = { 1e-12, 1e-12, 1e-12, 1e-12};//, 1e-12, 1e-12, 1e-12, 1e-12};
double eY10[nMass];
double eY11[nMass];
double eY20[nMass];
double eY21[nMass];
double eY30[nMass];
double eY31[nMass];

double xJES[nJES] = {0.96, 0.98, 1.00, 1.02, 1.04};
TString xJES_[nJES] = {"0_96", "0_98", "1_00", "1_02", "1_04"};
double y00 [nJES];
double y01 [nJES];
double y2  [nJES];
double y3  [nJES];
double y4  [nJES];
double y5  [nJES];

double ex  [nJES] = { 1e-12, 1e-12, 1e-12, 1e-12, 1e-12};
double ey0 [nJES];
double ey1 [nJES];
double ey2 [nJES];
double ey3 [nJES];
double ey4 [nJES];
double ey5 [nJES];
double ey6 [nJES];

double params[12];

enum styles             { kDown, kNominal, kUp};
int color_   [ nJES ] = { kMagenta, kRed+1, kBlack, kGreen+1, kBlue};
int marker_  [ nJES ] = {32, 23, 20, 22, 26};
int line_    [ nJES ] = {2, 7, 1, 9, 10};
int width_   [ nJES ] = {2, 2, 3, 2, 2};

void FindParametersMass(int iMass);
TH1F* FindParameters(TString filename, int i);

namespace cb {
  double A(const double alpha, const double power) {
    return TMath::Power(power / TMath::Abs(alpha), power) * TMath::Exp(-alpha*alpha/2);
  };
  double B(const double alpha, const double power) {
    return power / TMath::Abs(alpha) - TMath::Abs(alpha);
  };
}

// 7 parameters: [0] -> [6]
double crystalBall(const double* xx, const double* p) //m_top: the probability densities for the wp and um
{
  double N     = p[0]; //N_wp = 15 ???; N_un = 3 ??
  double mu    = p[1];
  double sigma = p[2];
  double alpha = p[3];
  double power = p[4];
  double t = (xx[0] - mu) / sigma;
  
  if(t < alpha) 
    return N * TMath::Exp(-t*t/2);
  else
    return N * cb::A(alpha,power) * TMath::Power(cb::B(alpha,power) + t, -power);
}

double asymGaus(const double* xx, const double* p)
{
  double N      = p[0];
  double mu     = p[1];
  double sigma1 = p[2];
  double sigma2 = p[3];
  double t1     = (xx[0] - mu) / sigma1;
  double t2     = (xx[0] - mu) / sigma2;
  
  double N1     = 1./sqrt(2*TMath::Pi()*sigma1*sigma1);
  double N2     = 1./sqrt(2*TMath::Pi()*sigma2*sigma2);
  
  N = (N1+N2)/2. * p[0];
  
  if(xx[0] < mu)
//     if(t1 < 0)
//     return N * TMath::Exp(t1*t1/2);
    return N * TMath::Exp(-t1*t1/2);
  else
//     return N * TMath::Exp(t2*t2/2);
    return N * TMath::Exp(-t2*t2/2);
}

void DrawCutLine(double cutval, double maximum)
{
  TLine *cut = new TLine();
  cut->SetLineWidth(2);
  cut->SetLineStyle(7);
  cut->SetLineColor(1);
  cut->DrawLine(cutval, 1e-6,cutval, maximum);
}

void observableByJESByMass(int pTarget = 1, int pObs = 0, int pLepton = 1) {
  target = pTarget;
  obs    = pObs;
  lepton = pLepton;
  
  if (plotByMass || pas) {
    iMassMin = 0;
    iMassMax = 3;
  }
  
  TStyle* tdrStyle = setTDRStyle();
  TH1::SetDefaultSumw2();
  tdrStyle->SetStatX(0.);
  tdrStyle->SetStatY(0.);

  for (int i = iMassMin; i<iMassMax; i++) {
	  std::cout << "### FindParametersMass(" << i << "); ###" << std::endl;
    FindParametersMass(i);
		/*
		mypause();
		myflush ( stdin );
		//*/
	}
  
  TF1* linearFit = new TF1("linearFit", "[0]+(x-172.5)*[1]");
  linearFit->SetLineWidth(2);
  linearFit->SetLineColor(kRed+1);
  
  TCanvas* cByMass = new TCanvas("cByMass", "cByMass", 1200, 1800);
  cByMass->Divide(2,3);
  
  cByMass->cd(1);
  
  TGraphErrors* gr = new TGraphErrors(nMass,X,Y10,eX,eY10);
  gr->SetTitle("10");
  gr->Draw("A*");
  //linearFit->SetParNames("p_{#mu}^{const}", "p_{#mu}^{m_{t}}");
  gr->Fit("linearFit", "EM");
  std::cout << linearFit->GetParameter(0) << " " << linearFit->GetParameter(1) << std::endl;
  params[0] = linearFit->GetParameter(0);
  params[1] = linearFit->GetParameter(1);
  
  gr->GetXaxis()->SetTitle("m_{t,gen} [GeV]");
  gr->GetYaxis()->SetTitle("#mu const [GeV]");
  
  cByMass->Update();
  TPaveStats *st1 = (TPaveStats*)gr->FindObject("stats");
  st1->SetFillColor(kWhite);
  st1->SetX1NDC(0.5); //new x start position
  st1->SetY1NDC(0.15); //new y start position
  st1->SetX2NDC(0.95); //new x end position
  st1->SetY2NDC(0.35); //new x end position
  gr->Draw("A*");
  
  cByMass->cd(2);
  
  gr = new TGraphErrors(nMass,X,Y11,eX,eY11);
  gr->SetTitle("11");
  gr->Draw("A*");
  //linearFit->SetParNames("p_{#mu}^{JES}", "p_{#mu}^{m_{t},JES}");
  gr->Fit("linearFit", "EM");
  params[2] = linearFit->GetParameter(0);
  params[3] = linearFit->GetParameter(1);
  
  gr->GetXaxis()->SetTitle("m_{t,gen} [GeV]");
  gr->GetYaxis()->SetTitle("#mu slope [GeV]");
  
  cByMass->Update();
  TPaveStats *st2 = (TPaveStats*)gr->FindObject("stats");
  st2->SetFillColor(kWhite);
  st2->SetX1NDC(0.5); //new x start position
  st2->SetY1NDC(0.15); //new y start position
  st2->SetX2NDC(0.95); //new x end position
  st2->SetY2NDC(0.35); //new x end position
  gr->Draw("A*");
  
  cByMass->cd(3);
  
  gr = new TGraphErrors(nMass,X,Y20,eX,eY20);
  gr->SetTitle("20");
  gr->Draw("A*");
  //linearFit->SetParNames("p_{#sigma}^{const}", "p_{#sigma}^{m_{t}}");
  gr->Fit("linearFit", "EM");
  params[4] = linearFit->GetParameter(0);
  params[5] = linearFit->GetParameter(1);
  
  gr->GetXaxis()->SetTitle("m_{t,gen} [GeV]");
  switch(obs) {
    case 0: {
      gr->GetYaxis()->SetTitle("#sigma const [GeV]");
      break;
    }
    case 1: {
      gr->GetYaxis()->SetTitle("#sigma_{1} const [GeV]");
      break;
    }
  }
  
  cByMass->Update();
  TPaveStats *st3 = (TPaveStats*)gr->FindObject("stats");
  st3->SetFillColor(kWhite);
  st3->SetX1NDC(0.5); //new x start position
  st3->SetY1NDC(0.15); //new y start position
  st3->SetX2NDC(0.95); //new x end position
  st3->SetY2NDC(0.35); //new x end position
  gr->Draw("A*");
  
  cByMass->cd(4);
  
  gr = new TGraphErrors(nMass,X,Y21,eX,eY21);
  gr->SetTitle("21");
  gr->Draw("A*");
  //linearFit->SetParNames("p_{#sigma}^{JES}", "p_{#sigma}^{m_{t},JES}");
  gr->Fit("linearFit", "EM");
  params[6] = linearFit->GetParameter(0);
  params[7] = linearFit->GetParameter(1);
  
  gr->GetXaxis()->SetTitle("m_{t,gen} [GeV]");
  switch(obs) {
    case 0: {
      gr->GetYaxis()->SetTitle("#sigma slope [GeV]");
      break;
    }
    case 1: {
      gr->GetYaxis()->SetTitle("#sigma_{1} slope [GeV]");
      break;
    }
  }
  
  cByMass->Update();
  TPaveStats *st4 = (TPaveStats*)gr->FindObject("stats");
  st4->SetFillColor(kWhite);
  st4->SetX1NDC(0.5); //new x start position
  st4->SetY1NDC(0.15); //new y start position
  st4->SetX2NDC(0.95); //new x end position
  st4->SetY2NDC(0.35); //new x end position
  gr->Draw("A*");
  
  cByMass->cd(5);
  
  gr = new TGraphErrors(nMass,X,Y30,eX,eY30);
  gr->SetTitle("30");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  params[8] = linearFit->GetParameter(0);
  params[9] = linearFit->GetParameter(1);
  
  gr->GetXaxis()->SetTitle("m_{t,gen} [GeV]");
  switch(obs) {
    case 0: {
      gr->GetYaxis()->SetTitle("#alpha const");
      break;
    }
    case 1: {
      gr->GetYaxis()->SetTitle("#sigma_{2} const");
      break;
    }
  }
  
  cByMass->Update();
  TPaveStats *st5 = (TPaveStats*)gr->FindObject("stats");
  st5->SetFillColor(kWhite);
  st5->SetX1NDC(0.5); //new x start position
  st5->SetY1NDC(0.15); //new y start position
  st5->SetX2NDC(0.95); //new x end position
  st5->SetY2NDC(0.35); //new x end position
  gr->Draw("A*");
  
  cByMass->cd(6);
  
  gr = new TGraphErrors(nMass,X,Y31,eX,eY31);
  gr->SetTitle("31");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  params[10] = linearFit->GetParameter(0);
  params[11] = linearFit->GetParameter(1);
  
  gr->GetXaxis()->SetTitle("m_{t,gen} [GeV]");
  switch(obs) {
    case 0: {
      gr->GetYaxis()->SetTitle("#alpha slope");
      break;
    }
    case 1: {
      gr->GetYaxis()->SetTitle("#sigma_{2} slope");
      break;
    }
  }
  
  cByMass->Update();
  TPaveStats *st6 = (TPaveStats*)gr->FindObject("stats");
  st6->SetFillColor(kWhite);
  st6->SetX1NDC(0.5); //new x start position
  st6->SetY1NDC(0.15); //new y start position
  st6->SetX2NDC(0.95); //new x end position
  st6->SetY2NDC(0.35); //new x end position
  gr->Draw("A*");
  
  if (!pas && !plotByMass) {
    TString pathByMass("../plot/template/"); pathByMass+= sFObs[obs]; pathByMass += "_bymass";
    pathByMass += "_"; pathByMass += sTarget[abs(target%8)];
    pathByMass += "_"; pathByMass += sLepton[lepton]; pathByMass += ".eps";
    cByMass->Print(pathByMass);
  }
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  TCanvas* cSetOfCurves = new TCanvas("cSetOfCurves", "cSetOfCurves", 1800, 600);
  cSetOfCurves->Divide(3,1);
  
  cSetOfCurves->cd(1);
  
  TMultiGraph *mg1 = new TMultiGraph();
  for (int i = iMassMin; i<iMassMax; i++) {
    mg1->Add(gr1[i]);
  }
  
  mg1->SetTitle("#mu [GeV];JSF;#mu [GeV]");
  mg1->Draw("AP");
  mg1->GetXaxis()->SetLimits(0.945, 1.055);
  mg1->GetXaxis()->SetLabelSize(1.25*mg1->GetXaxis()->GetLabelSize());
  mg1->GetYaxis()->SetLabelSize(1.25*mg1->GetYaxis()->GetLabelSize());
  mg1->GetXaxis()->SetTitleSize(1.25*mg1->GetXaxis()->GetTitleSize());
  mg1->GetYaxis()->SetTitleSize(1.25*mg1->GetYaxis()->GetTitleSize());
  //mg1->GetYaxis()->SetTitleOffset(1./1.25*mg1->GetYaxis()->GetTitleOffset());
  mg1->Draw("AP");
  cSetOfCurves->Update();
  
  cSetOfCurves->cd(2);
  
  TMultiGraph *mg2 = new TMultiGraph();
  for (int i = iMassMin; i<iMassMax; i++) {
    mg2->Add(gr2[i]);
  }
  
  mg2->SetTitle("#sigma [GeV];JSF;#sigma [GeV]");
  mg2->Draw("AP");
  mg2->GetXaxis()->SetLimits(0.945, 1.055);
  mg2->GetXaxis()->SetLabelSize(1.25*mg2->GetXaxis()->GetLabelSize());
  mg2->GetYaxis()->SetLabelSize(1.25*mg2->GetYaxis()->GetLabelSize());
  mg2->GetXaxis()->SetTitleSize(1.25*mg2->GetXaxis()->GetTitleSize());
  mg2->GetYaxis()->SetTitleSize(1.25*mg2->GetYaxis()->GetTitleSize());
  mg2->GetYaxis()->SetTitleOffset(1./1.25*mg2->GetYaxis()->GetTitleOffset());
  mg2->Draw("AP");
  cSetOfCurves->Update();
  
  cSetOfCurves->cd(3);
  
  TLegend *leg1 = new TLegend(0.1, 0.1, 0.9, 0.9);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  for (int i = iMassMin; i<iMassMax; i++) {
    leg1->AddEntry( gr1[i], sX[i], "PL" );
  }
  
  leg1->Draw();
  
  if (!pas && !plotByMass) {
    TString pathByJES("../plot/template/"); pathByJES+= sFObs[obs]; pathByJES += "_byjes";
    pathByJES += "_"; pathByJES += sTarget[abs(target%8)];
    pathByJES += "_"; pathByJES += sLepton[lepton]; pathByJES += ".eps";
    cSetOfCurves->Print(pathByJES);
  }
  
  std::cout << "Parameters for IdeogramCombLikelihood.cxx:" << std::endl;
  ofstream parFile;
  parFile.open ("../plot/template/parameters.txt", ios::app);
  parFile   << sFObs[obs] << " " << sTarget[abs(target%8)] << " " << sLepton[lepton] << "\n";
  for (int i=0; i<12; i++) {
    std::cout << params[i] << "|";
    parFile   << params[i] << "|";
  }
  std::cout << std::endl;
  parFile   << "\n";
  parFile.close();
}

void FindParametersMass(int iMass)
{
  setTDRStyle();
  gStyle->SetOptFit(0);
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0); 
    
  TCanvas* cObservable = new TCanvas("cObservable", "cObservable", 600, 600);
  
  cObservable->cd();
  
  TF1* linearFit = new TF1("linearFit", "[0]+(x-1)*[1]");
  linearFit->SetLineWidth(1);
  linearFit->SetLineColor(kRed+1);
	linearFit->SetParLimits(1, -50, 200);
	linearFit->SetParameters(100, 50);
	
  TH1F* h096 = new TH1F();
  TH1F* h098 = new TH1F();
  TH1F* h100 = new TH1F();
  TH1F* h102 = new TH1F();
  TH1F* h104 = new TH1F();
  TH1F* h166 = new TH1F();
  TH1F* h178 = new TH1F();
  
  
  if (!plotByMass) {
    switch(iMass) {
      case 0: {
	h096 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1665_0.96", 0);
        h098 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1665_0.98", 1);
        h100 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1665_1.00", 2);
	h102 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1665_1.02", 3);
        h104 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1665_1.04", 4);
        break;
      }
      case 1: {
	h096 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1695_0.96", 0);
        h098 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1695_0.98", 1);
        h100 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1695_1.00", 2);
	h102 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1695_1.02", 3);
        h104 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1695_1.04", 4);
        break;
      }
      case 2: {
	h096 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1725_0.96", 0);
        h098 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1725_0.98", 1);
        h100 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1725_1.00", 2);
	h102 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1725_1.02", 3);
        h104 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1725_1.04", 4);
        break;
      }
      case 3: {
	h096 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1755_0.96", 0);
        h098 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1755_0.98", 1);
        h100 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1755_1.00", 2);
	h102 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1755_1.02", 3);
        h104 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1755_1.04", 4);
        break;
      }
      case 4: {
	h096 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1785_0.96", 0);
        h098 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1785_0.98", 1);
        h100 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1785_1.00", 2);
	h102 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1785_1.02", 3);
        h104 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1785_1.04", 4);
        break;
      }
		  case 5: {
		    		    std::cout<<"ERROR: atm only 5 masses valid!"<<std::endl;
        h096 = FindParameters("/nfs/dust/cms/user/mseidel/trees_paper/Summer12_TTJetsMS1755_0.96", 0);
        h100 = FindParameters("/nfs/dust/cms/user/mseidel/trees_paper/Summer12_TTJetsMS1755_1.00", 1);
        h104 = FindParameters("/nfs/dust/cms/user/mseidel/trees_paper/Summer12_TTJetsMS1755_1.04", 2);
        break;
      }
		  case 6: {
        h096 = FindParameters("/nfs/dust/cms/user/mseidel/trees_paper/Summer12_TTJetsMS1785_0.96", 0);
        h100 = FindParameters("/nfs/dust/cms/user/mseidel/trees_paper/Summer12_TTJetsMS1785_1.00", 1);
        h104 = FindParameters("/nfs/dust/cms/user/mseidel/trees_paper/Summer12_TTJetsMS1785_1.04", 2);
        break;
      }
    }
  }
  else if (plotByMass) {
    std::cout << "PLOT BY MASS" << std::endl;
    h096 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1665_1.00", 0);
    h098 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1695_1.00", 1);
    h100 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1725_1.00", 2);
    h104 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1755_1.00", 3);
    h166 = FindParameters("/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/Trees_80X/2016B_JESVariations_ICHEPdata/TT_TuneCUETP8M1_13TeV_powheg_pythia8_Spring16_MS1785_1.00", 3);
  }
  else  {std::cout<<"wtf, how did he got here???"<<std::endl;
    if (obs==0) h166 = FindParameters("/nfs/dust/cms/user/mseidel/trees/Summer12_TTJets1665_1.00", 3);
    h096 = FindParameters("/nfs/dust/cms/user/mseidel/trees/Summer12_TTJets1725_0.96", 0);
    h100 = FindParameters("/nfs/dust/cms/user/mseidel/trees/Summer12_TTJets1725_1.00", 1);
    h104 = FindParameters("/nfs/dust/cms/user/mseidel/trees/Summer12_TTJets1725_1.04", 2);
    if (obs==0) h178 = FindParameters("/nfs/dust/cms/user/mseidel/trees/Summer12_TTJets1785_1.00", 4); 
  }
  
  h096->Draw();
  if (obs == 0) {
    h096->GetXaxis()->SetRangeUser(100, 240);
    h098->GetXaxis()->SetRangeUser(100, 240);
    h100->GetXaxis()->SetRangeUser(100, 240);
    h102->GetXaxis()->SetRangeUser(100, 240);
    h104->GetXaxis()->SetRangeUser(100, 240);
  }
  if (obs == 1) {
    h096->GetXaxis()->SetRangeUser(60, 110);
    h098->GetXaxis()->SetRangeUser(60, 110);
    h100->GetXaxis()->SetRangeUser(60, 110);
    h102->GetXaxis()->SetRangeUser(60, 110);
    h104->GetXaxis()->SetRangeUser(60, 110);

  }
  h098->Draw("SAME");
  h100->Draw("SAME");
  h102->Draw("SAME");
  h104->Draw("SAME"); 
  // ---
  //    create legend
  // ---
  TLegend *leg0 = new TLegend(0.45, 0.7, 0.9, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  if (!plotByMass) {
    leg0->AddEntry( h096, "JSF = 0.96", "PL");
    leg0->AddEntry( h098, "JSF = 0.98", "PL");
    leg0->AddEntry( h100, "JSF = 1.00", "PL");
    leg0->AddEntry( h102, "JSF = 1.02", "PL");
    leg0->AddEntry( h104, "JSF = 1.04", "PL");
    leg0->AddEntry((TObject*)0, sX[iMass], "");
  }
  else if (plotByMass) {
    leg0->AddEntry( h096, "m_{t,gen} = 166.5 GeV", "PL");
    leg0->AddEntry( h098, "m_{t,gen} = 169.5 GeV", "PL");
    leg0->AddEntry( h100, "m_{t,gen} = 172.5 GeV", "PL");
    leg0->AddEntry( h102, "m_{t,gen} = 175.5 GeV", "PL");
    leg0->AddEntry( h104, "m_{t,gen} = 178.5 GeV", "PL");
    leg0->AddEntry((TObject*)0, "JSF = 1.00", "");
  }
  else {
    if (obs==0) leg0->AddEntry( h166, "m_{t,gen} = 169.5 GeV", "PL");
    leg0->AddEntry( h096, "m_{t,gen} = 172.5 GeV, JSF -4%", "PL");
    leg0->AddEntry( h100, "m_{t,gen} = 172.5 GeV", "PL");
    leg0->AddEntry( h104, "m_{t,gen} = 172.5 GeV, JSF +4%", "PL");
    if (obs==0) leg0->AddEntry( h178, "m_{t,gen} = 175.5 GeV", "PL");
  }
  
  leg0->Draw();

  CMS_lumi(cObservable, 4, 0., true, "Privat work");
  
  if (obs==0) DrawCutLine(172.5, h096->GetMaximum()*1.);
  else DrawCutLine(80.4, h096->GetMaximum()*1.);
  h096->GetYaxis()->SetRangeUser(0, h096->GetMaximum()*1.5);
  
  gStyle->SetOptFit(1);
  
  TString path("../plot/template/"); path+= sFObs[obs]; path += "_";
  if (plotByMass) {
    path += "mass";
  }
  else {
    path += "jes";
    path += X_[iMass];
  }
  path += "_"; path += sTarget[abs(target%8)];
  path += "_"; path += sLepton[lepton]; path += ".eps";
  cObservable->Print(path);
  
  TCanvas* cObservablePar = new TCanvas("cObservablePar", "cObservablePar", 1600, 400);
  cObservablePar->Divide(4,1);
  
  int masscolor[]  = {kRed+1, kMagenta+1, kViolet+1, kBlue+1, kAzure+1, kCyan+1, kGreen+1};
  int massmarker[] = { 23, 32, 32, 20, 26, 26, 22};
  double masssize[]   = { 1, 1, 0.5, 1, 0.5, 1, 1};
  
  cObservablePar->cd(1);
  gr1[iMass] = new TGraphErrors(nJES,xJES,y01,ex,ey1);
  gr1[iMass]->SetTitle("p1");
  gr1[iMass]->SetLineColor(masscolor[iMass]);
  gr1[iMass]->SetMarkerColor(masscolor[iMass]);
  gr1[iMass]->SetMarkerStyle(massmarker[iMass]);
  gr1[iMass]->SetMarkerSize(1.5*masssize[iMass]);
  linearFit->SetLineColor(masscolor[iMass]);
  gr1[iMass]->SetLineStyle(iMass+1);
  linearFit->SetLineStyle(iMass+1);
  gr1[iMass]->Draw("AP");
  gr1[iMass]->Fit("linearFit", "EM");
  
  gr1[iMass]->GetXaxis()->SetTitle("JSF");
  gr1[iMass]->GetYaxis()->SetTitle("#mu");
  
  Y10 [iMass] = linearFit->GetParameter(0);
  eY10[iMass] = linearFit->GetParError(0);
  
  Y11 [iMass] = linearFit->GetParameter(1);
  eY11[iMass] = linearFit->GetParError(1);
  
  
  cObservablePar->cd(2);
  gr2[iMass] = new TGraphErrors(nJES,xJES,y2,ex,ey2);
  gr2[iMass]->SetTitle("p2");
  gr2[iMass]->SetLineColor(masscolor[iMass]);
  gr2[iMass]->SetMarkerColor(masscolor[iMass]);
  gr2[iMass]->SetMarkerStyle(massmarker[iMass]);
  gr2[iMass]->SetMarkerSize(1.5*masssize[iMass]);
  linearFit->SetLineColor(masscolor[iMass]);
  gr2[iMass]->SetLineStyle(iMass+1);
  linearFit->SetLineStyle(iMass+1);
  gr2[iMass]->Draw("AP");
  linearFit->SetParameters(100, 50);
  gr2[iMass]->Fit("linearFit", "EM");
  
  gr2[iMass]->GetXaxis()->SetTitle("JSF");
  gr2[iMass]->GetYaxis()->SetTitle("#sigma");
  
  Y20 [iMass] = linearFit->GetParameter(0);
  eY20[iMass] = linearFit->GetParError(0);
  
  Y21 [iMass] = linearFit->GetParameter(1);
  eY21[iMass] = linearFit->GetParError(1);
  
  
  cObservablePar->cd(3);
  gr3[iMass] = new TGraphErrors(nJES,xJES,y3,ex,ey3);
  gr3[iMass]->SetTitle("p3");
  gr3[iMass]->Draw("A*");
  linearFit->SetParameters(100, 50);
  gr3[iMass]->Fit("linearFit", "EM");
  
  gr3[iMass]->GetXaxis()->SetTitle("JSF");
  if (obs == 0) gr3[iMass]->GetYaxis()->SetTitle("#alpha");
  if (obs == 1) gr3[iMass]->GetYaxis()->SetTitle("#sigma2");
  
  Y30 [iMass] = linearFit->GetParameter(0);
  eY30[iMass] = linearFit->GetParError(0);
  
  Y31 [iMass] = linearFit->GetParameter(1);
  eY31[iMass] = linearFit->GetParError(1);
  
  if (target!=1) {
    cObservablePar->cd(4);
    gr4[iMass] = new TGraphErrors(nJES,xJES,y4,ex,ey4);
    gr4[iMass]->SetTitle("p4");
    gr4[iMass]->Draw("A*");
    linearFit->SetParameters(100, 50);
    gr4[iMass]->Fit("linearFit", "EM");
    
    gr4[iMass]->GetXaxis()->SetTitle("JSF");
    gr4[iMass]->GetYaxis()->SetTitle("power");
  }
}



TH1F* FindParameters(TString filename, int i)
{
  TChain* eventTree = new TChain("analyzeHitFit/eventTree");
  TString filenameEle = filename; filenameEle += "_electron/*.root";
  TString filenameMu  = filename; filenameMu  += "_muon/*.root";
  if (lepton != 1) {
    std::cout << filenameEle << std::endl;
    eventTree->Add(filenameEle);
  }
  if (lepton != 0) {
    std::cout << filenameMu << std::endl;
    eventTree->Add(filenameMu);
  }
  
  TF1* fit = new TF1();
  
  if (obs == 0) {
    switch(target) {
      case   1: {
        fit = new TF1("fit", "[0]*TMath::Voigt(x-[1], [2], [3])", 100, 400);
        fit->SetLineColor(kBlack);
        fit->SetLineWidth(width_[i]);
        fit->SetParameters(100000, 170, 10, 2);
        
        fit->SetParLimits(0, 1, 1000000);
        fit->SetParLimits(1, 150, 200);
        fit->SetParLimits(2, 5, 15);
        fit->SetParLimits(3, 2, 2);
        
        break;
      }
      case   0: {
        fit = new TF1("fit", crystalBall, 100, 400, 5);
        fit->SetLineColor(kBlack);
        fit->SetLineWidth(width_[i]);
        
        double power = 15.;
        
        fit->SetParNames("N", "#mu", "#sigma", "#alpha", "power");
        fit->SetParameters(1, 170, 25, 0.45, power);
        
        fit->SetParLimits(0, 0, 1000000);
        fit->SetParLimits(1, 150, 200);
        fit->SetParLimits(2, 15, 50);
        fit->SetParLimits(3, 0.1, 0.95);
        fit->SetParLimits(4, power, power);
        
        break;
      }
      case -10: {
        fit = new TF1("fit", crystalBall, 100, 400, 5);
        fit->SetLineColor(kBlack);
        fit->SetLineWidth(width_[i]);
        
        double power = 3.; //TODO oder 5?
        
        fit->SetParNames("N", "#mu", "#sigma", "#alpha", "power");
        fit->SetParameters(1, 170, 15, 0.45, power);
        
        fit->SetParLimits(0, 0, 1000000);
        fit->SetParLimits(1, 150, 200);
        fit->SetParLimits(2, 10, 30);
        fit->SetParLimits(3, 0.5, 1.95);
        fit->SetParLimits(4, power, power);
        
        break;
      }
    }
  }
  
  else if (obs == 1) {
    fit = new TF1("fit", asymGaus, 0, 1000, 4);
    fit->SetLineColor(kBlack);
    fit->SetLineWidth(width_[i]);
    
    fit->SetParNames("N", "#mu", "#sigma1", "#sigma2");
    fit->SetParameters(1000, 80, 5, 5);
    
    fit->SetParLimits(0, 100, 1000000);
    fit->SetParLimits(1, 50, 150);
    fit->SetParLimits(2, 0, 10);
    fit->SetParLimits(3, 0, 15);
  }
  
  else if (obs == 4) {
    fit = new TF1("fit", crystalBall, 0, 1000, 5);
    fit->SetLineColor(kBlack);
    fit->SetLineWidth(width_[i]);
    
    fit->SetParNames("N", "#mu", "#sigma", "#alpha", "power");
    fit->SetParameters(1, 0, 1, 2, 3);
    
    fit->SetParLimits(0, 0, 1000000);
    fit->SetParLimits(1, 0.5, 2);
    fit->SetParLimits(2, 0.2, 2);
    fit->SetParLimits(3, 0, 2);
    fit->SetParLimits(4, 0, 20);
  }
  
  else if (obs == 5) {
    fit = new TF1("fit", "[0]*TMath::Voigt(x-[1], [2], [3])");

    fit->SetLineColor(kBlack);
    fit->SetLineWidth(width_[i]);
    fit->SetParameters(100000, 0.0, 1e-12, 1.35);
    
    fit->SetParLimits(0, 1, 1000000);
    fit->SetParLimits(1, 0, 2);
    fit->SetParLimits(2, 1e-12, 1e-12);
    fit->SetParLimits(3, 0, 10);
  }
  
  fit->SetNpx(300);
  
  TString sObservable;
  switch(obs) {
    case 0: {
      switch(target) {
        case   1: {
          sObservable = "top.fitTop1.M() >> htemp(30, 100, 250)";
          break;
        }
        case   0: {
          sObservable = "top.fitTop1.M() >> htemp(30, 100, 400)";
          break;
        }
        case -10: {
          sObservable = "top.fitTop1.M() >> htemp(30, 100, 400)";
          break;
        }
      }
      break;
    }
    case 1: {
      sObservable = "top.recoW1.M() >> htemp(30, 60, 120)";
      break;
    }
    case 4: {
      sObservable = "deltaRHadWHadB/deltaRHadQHadQBar >> htemp(40, 0, 4)";
      break;
    }
    case 5: {
      sObservable = "-(hadWPt-lepWPt)/(hadBPt-lepBPt) >> htemp(50, -5, 5)";
      break;
    }
  }
  TString sCutAndWeight("(top.fitProb*weight.combinedWeight)*(");
  switch(target) {
    case   1: {
      sCutAndWeight += "top.combinationType==1";
      break;
    }
    case   0: {
      sCutAndWeight += "top.combinationType>=2 & top.combinationType<=4";
      break;
    }
    case -10: {
      sCutAndWeight += "(top.combinationType==6 | top.combinationType<0)";
      break;
    }
  }
  sCutAndWeight += " & top.fitProb > 0.2 & Min$(sqrt(pow(top.recoW2Prod1[0].Eta()-jet.jet.Eta(),2) + pow(TVector2::Phi_mpi_pi(top.recoW2Prod1[0].Phi()-jet.jet.Phi()),2)))>0.3)";

  std::cout << sCutAndWeight << std::endl;
  
  // Get observable
  eventTree->Draw(sObservable, sCutAndWeight, "");
  
  TH1F *h1 = (TH1F*)gDirectory->Get("htemp")->Clone();
  
  switch(obs) {
    case 0: {
      switch(target) {
        case   1: {
          h1->GetYaxis()->SetTitle("Fraction of entries / 5 GeV");
          break;
        }
        case   0: {
          h1->GetYaxis()->SetTitle("Fraction of entries / 10 GeV");
          break;
        }
        case -10: {
          h1->GetYaxis()->SetTitle("Fraction of entries / 10 GeV");
          break;
        }
      }
      h1->GetXaxis()->SetTitle("m_{t," + sTarget[abs(target%8)] + "}^{fit} [GeV]");
      break;
    }
    case 1: {
      h1->GetYaxis()->SetTitle("Fraction of entries / 2 GeV");
      h1->GetXaxis()->SetTitle("m_{W," + sTarget[abs(target%8)] + "}^{reco} [GeV]");
      break;
    }
  }
  
  h1->SetFillColor(color_[i]);
  h1->SetLineColor(color_[i]);
  h1->SetMarkerColor(color_[i]);
  h1->SetMarkerStyle(marker_[i]);
  h1->SetMarkerSize(1.5);
  fit->SetLineColor(color_[i]);
  fit->SetLineStyle(line_[i]);
  
  TFitResultPtr r = h1->Fit("fit","WEMSR");
  
  std::cout << "chi2/ndf = " << r->Chi2() << "/" << r->Ndf() << std::endl;

  y00[i] = fit->GetParameter(0);
  ey0[i] = 2*fit->GetParError(0);
 
  y01[i] = fit->GetParameter(1);
  ey1[i] = 2*fit->GetParError(1);
  
  y2 [i] = fit->GetParameter(2);
  ey2[i] = 2*fit->GetParError(2);
  if (ey2[i] == 0) ey2[i] = 1e-12;
  
  y3 [i] = fit->GetParameter(3);
  ey3[i] = 2*fit->GetParError(3);
  
  y4 [i] = fit->GetParameter(4);
  ey4[i] = 2*fit->GetParError(4);
  
  TF1 *fitted = h1->GetFunction("fit");

  fitted->SetParameter(0, fitted->GetParameter(0)/h1->Integral());
  h1->Scale(1/h1->Integral());
  
  return h1;
}

void batchObservableByJESByMass() { //that the "real" start
  pas = false;
  
  for (int t = 0; t < 3; ++t) { // target
    for (int o = 0; o < 2; ++o) { // obs
     // for (int l = 1; l < 2; ++l) { // lepton  atm only muon
        observableByJESByMass(iTarget[t], o, 1);
     // }
    }
  }
  
  plotByMass = true;
  
  for (int t = 0; t < 3; ++t) { // target
    for (int o = 0; o < 2; ++o) { // obs
      //for (int l = 1; l < 2; ++l) { // lepton  atm only muon
        observableByJESByMass(iTarget[t], o, 1);
     // }
    }
  }
}
