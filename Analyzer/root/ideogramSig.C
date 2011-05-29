#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"

double x[3] = {166.5, 172.5, 178.5};
double y0[3];
double y1[3];
double y2[3];
double y3[3];
double y4[3];
double y5[3];
double ex[3] = {0.001, 0.001, 0.001};
double ey0[3];
double ey1[3];
double ey2[3];
double ey3[3];
double ey4[3];
double ey5[3];
double ey6[3];

void ideogramSig()
{
  // open input file/directory
//  file1725 = new TFile("TopMass/Analyzer/root/analyzeTop_1725.root");
//  file1785 = new TFile("TopMass/Analyzer/root/analyzeTop_1785.root");
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1);
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 900, 900);
  canvas->Divide(3,2);
  
  canvas->cd(1);
  FindParameters("analyzeTop_1665.root", 0);
  canvas->cd(2);
  FindParameters("analyzeTop_1725.root", 1);
  canvas->cd(3);
  FindParameters("analyzeTop_1785.root", 2);
  
  /*
  canvas->cd(4);
  gr = new TGraphErrors(3,x,y0,ex,ey0);
  gr->SetTitle("p0");
  gr->Draw("A*");
  gr->Fit("pol1", "M");
  */
  
  canvas->cd(4);
  gr = new TGraphErrors(3,x,y1,ex,ey1);
  gr->SetTitle("p1");
  gr->Draw("A*");
  gr->Fit("pol1", "EM");
  
  canvas->cd(5);
  gr = new TGraphErrors(3,x,y2,ex,ey2);
  gr->SetTitle("p2");
  gr->Draw("A*");
  gr->Fit("pol1", "EM");
  
  canvas->cd(6);
  gr = new TGraphErrors(3,x,y3,ex,ey3);
  gr->SetTitle("p3");
  gr->Draw("A*");
  gr->Fit("pol1", "EM");
  
  /*
  canvas->cd(6);
  gr = new TGraphErrors(3,x,y4,ex,ey4);
  gr->SetTitle("p4");
  gr->Draw("A*");
  gr->Fit("pol1", "M");
  
  /*
  canvas->cd(9);
  gr = new TGraphErrors(3,x,y5,ex,ey5);
  gr->SetTitle("p5");
  gr->Draw("A*");
  gr->Fit("pol1");
  */
}

void FindParameters(TString filename, int i)
{

  TFile* file = new TFile(filename);
  analyzeHitFit->cd();

  TF1* voigt = new TF1("voigt", "[0]*TMath::Voigt(x-[1], [2], [3])");
  voigt->SetLineColor(kRed+1);
  voigt->SetParameters(100000, 170, 10, 2);
  
  voigt->SetParLimits(0, 20, 1000000);
  voigt->SetParLimits(1, 150, 200);
  voigt->SetParLimits(2, 5, 20);
  voigt->SetParLimits(3, 2, 2);
  
  TH1F* hSig;
  
  eventTree->Draw("hadTopMass >> hSig(50, 100, 250)", "(bProbSSV*hitFitProb)*(target==1 & (bProbSSV * hitFitProb) > 0.05)");
  
  TH1F *hSig = (TH1F*)gDirectory->Get("hSig");
  
  //hSig->Scale(1/hSig->Integral());
  
  hSig->Fit("voigt","LEM");

  y0 [i] = voigt->GetParameter(0);
  ey0[i] = voigt->GetParError(0);
 
  y1 [i] = voigt->GetParameter(1);
  ey1[i] = voigt->GetParError(1);
  
  y2 [i] = voigt->GetParameter(2);
  ey2[i] = voigt->GetParError(2);
  
  y3 [i] = voigt->GetParameter(3);
  ey3[i] = voigt->GetParError(3);

}
