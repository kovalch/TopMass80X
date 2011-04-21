#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"

double x[5] = {166.5, 172.5, 172.5, 172.5, 178.5};
double y0[5];
double y1[5];
double y2[5];
double y3[5];
double y4[5];
double y5[5];
double ex[5] = {0.001, 0.001, 0.001, 0.001, 0.001};
double ey0[5];
double ey1[5];
double ey2[5];
double ey3[5];
double ey4[5];
double ey5[5];
double ey6[5];

void ideogramCombBkgLogNormal()
{
  // open input file/directory
//  file1725 = new TFile("TopMass/Analyzer/root/analyzeTop_1725.root");
//  file1785 = new TFile("TopMass/Analyzer/root/analyzeTop_1785.root");
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1);
  
  /*TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 900, 900);
  TF1 *lognormal = new TF1("lognormal", "TMath::LogNormal(x, [0], [1], [2])", 0, 500);
  
  double xx = 266.5;
  
  double p0 =  1.27437 - 0.00409434 * xx;
  double p1 =  206.205 - 0.67529 * xx;
  double p2 = -187.656 + 1.65207 * xx;
  
  lognormal->SetParameters(p0, p1, p2);
  
  lognormal->Draw();*/
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 900, 900);
  canvas->Divide(4,3);
  
  canvas->cd(5);
  FindParameters("analyzeTop_1665.root", 0);
  canvas->cd(2);
  FindParameters("analyzeTop_1725_jes_up.root", 1);
  canvas->cd(6);
  FindParameters("analyzeTop_1725.root", 2);
  canvas->cd(10);
  FindParameters("analyzeTop_1725_jes_down.root", 3);
  canvas->cd(7);
  FindParameters("analyzeTop_1785.root", 4);
  
  canvas->cd(4);
  gr = new TGraphErrors(5,x,y0,ex,ey0);
  gr->SetTitle("p0");
  gr->Draw("A*");
  gr->Fit("pol1");
  
  canvas->cd(8);
  gr = new TGraphErrors(5,x,y1,ex,ey1);
  gr->SetTitle("p1");
  gr->Draw("A*");
  gr->Fit("pol1");
  
  canvas->cd(12);
  gr = new TGraphErrors(5,x,y2,ex,ey2);
  gr->SetTitle("p2");
  gr->Draw("A*");
  gr->Fit("pol1");
  
  
  
  
}

void FindParameters(TString filename, int i)
{

  TFile* file = new TFile(filename);
  analyzeKinFit->cd();

  TF1 *lognormal = new TF1("lognormal", "[3]*TMath::LogNormal(x, 0.5, [1], [2])", 100, 300);
  lognormal->SetParameters(0.5, 100, 100, 10000);
  
  lognormal->SetParLimits(0, 0.3, 1);
  lognormal->SetParLimits(1, 50, 200);
  lognormal->SetParLimits(2, 50, 200);
  
  TH1F* hCombBkg;
  
  eventTree->Draw("hadTopMass >> hCombBkg(100, 100, 300)","fitProb*(target==0)");
  
  TH1F *hCombBkg = (TH1F*)gDirectory->Get("hCombBkg");
  
  //hCombBkg->Scale(1/hCombBkg->Integral());
  
  hCombBkg->Fit("lognormal","LEMR");

  y0 [i] = lognormal->GetParameter(0);
  ey0[i] = lognormal->GetParError(0);
  //ey0[i] = sqrt(pow(1/lognormal->GetParameter(3)*lognormal->GetParError(0), 2)
  //          + pow(lognormal->GetParameter(0)/pow(lognormal->GetParameter(3), 2)*lognormal->GetParError(3), 2));
  
  y1 [i] = lognormal->GetParameter(1);
  ey1[i] = lognormal->GetParError(1);
  
  y2 [i] = lognormal->GetParameter(2);
  ey2[i] = lognormal->GetParError(2);


}
