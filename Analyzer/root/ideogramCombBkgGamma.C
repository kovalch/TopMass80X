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

void ideogramCombBkgGamma()
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
  
  canvas->cd(4);
  gr = new TGraphErrors(3,x,y0,ex,ey0);
  gr->SetTitle("p0");
  gr->Draw("A*");
  gr->Fit("pol1");
  
  canvas->cd(5);
  gr = new TGraphErrors(3,x,y1,ex,ey1);
  gr->SetTitle("p1");
  gr->Draw("A*");
  gr->Fit("pol1");
  
  canvas->cd(6);
  gr = new TGraphErrors(3,x,y2,ex,ey2);
  gr->SetTitle("p2");
  gr->Draw("A*");
  gr->Fit("pol1");
}

void FindParameters(TString filename, int i)
{

  TFile* file = new TFile(filename);
  analyzeKinFit->cd();

  TF1 *gamma = new TF1("gamma", "[0]*TMath::GammaDist(x, [3], [1], [2])");
  gamma->SetParameters(2000, 130, 20, 3);
  
  /*
  gamma->SetParLimits(0, 20, 100000);
  gamma->SetParLimits(1, 150, 200);
  gamma->SetParLimits(2, 10, 100);
  gamma->SetParLimits(3, 200, 300);
  gamma->SetParLimits(4, 20, 100);
  */
  
  TH1F* hCombBkg;
  
  eventTree->Draw("hadTopMass >> hCombBkg","fitProb*(target==0)");
  
  TH1F *hCombBkg = (TH1F*)gDirectory->Get("hCombBkg");
  
  //hCombBkg->Scale(1/hCombBkg->Integral());
  
  hCombBkg->Fit("gamma","L");

  y0 [i] = gamma->GetParameter(0);
  ey0[i] = gamma->GetParError(0);
  //ey0[i] = sqrt(pow(1/gamma->GetParameter(3)*gamma->GetParError(0), 2)
  //          + pow(gamma->GetParameter(0)/pow(gamma->GetParameter(3), 2)*gamma->GetParError(3), 2));
  
  y1 [i] = gamma->GetParameter(1);
  ey1[i] = gamma->GetParError(1);
  
  y2 [i] = gamma->GetParameter(2);
  ey2[i] = gamma->GetParError(2);


}
