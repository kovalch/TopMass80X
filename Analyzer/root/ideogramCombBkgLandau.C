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

void ideogramCombBkgLandau()
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
  gr->Fit("pol1", "M");
  
  canvas->cd(5);
  gr = new TGraphErrors(3,x,y1,ex,ey1);
  gr->SetTitle("p1");
  gr->Draw("A*");
  gr->Fit("pol1", "M");
  
  /*
  canvas->cd(6);
  gr = new TGraphErrors(3,x,y2,ex,ey2);
  gr->SetTitle("p2");
  gr->Draw("A*");
  gr->Fit("pol1");
  
  canvas->cd(7);
  gr = new TGraphErrors(3,x,y3,ex,ey3);
  gr->SetTitle("p3");
  gr->Draw("A*");
  gr->Fit("pol1");
  */
  
  canvas->cd(6);
  gr = new TGraphErrors(3,x,y2,ex,ey2);
  gr->SetTitle("p2");
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
  analyzeKinFit->cd();

  TF1 *tLandau = new TF1("tLandau", "[0] * TMath::Landau(x, [1], [2], 1)");
  tLandau->SetParameters(200, 160, 20);
  
  /*
  if (i == 1) tLandau->SetParLimits(1, y1[0] + 6, 200);
  if (i == 2) tLandau->SetParLimits(1, y1[0] + 12, 200);
  */
  /*
  tLandau->SetParLimits(0, 0, 10);
  tLandau->SetParLimits(1, 150, 200);
  tLandau->SetParLimits(2, 25, 25);
  */
  
  TH1F* hCombBkg;
  
  eventTree->Draw("hadTopMass >> hCombBkg(20, 100, 300)","(hadBProb*exp(-fitChi2))*(target==0 & hadBProb > 1e-3 & exp(-fitChi2) > 1e-3)");
  
  TH1F *hCombBkg = (TH1F*)gDirectory->Get("hCombBkg");
  
  //hCombBkg->Scale(1/hCombBkg->Integral());
  
  hCombBkg->Fit("tLandau","LEM");

  y0 [i] = tLandau->GetParameter(0);
  ey0[i] = tLandau->GetParError(0);
  /*
  ey0[i] = sqrt(pow(1/tLandau->GetParameter(3)*tLandau->GetParError(0), 2)
            + pow(tLandau->GetParameter(0)/pow(tLandau->GetParameter(3), 2)*tLandau->GetParError(3), 2));
  */
  
  y1 [i] = tLandau->GetParameter(1);
  ey1[i] = tLandau->GetParError(1);
  
  y2 [i] = tLandau->GetParameter(2);
  ey2[i] = tLandau->GetParError(2);

}
