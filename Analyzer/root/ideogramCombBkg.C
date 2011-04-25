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

void ideogramCombBkg()
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
  gr->SetTitle("p0/p3");
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
  analyzeKinFit->cd();

  TF1 *gaus2 = new TF1("gaus2", "[0]*TMath::Gaus(x, [1], [2], 1)+[3]*TMath::Gaus(x, [4], [5], 1)");
  gaus2->SetParameters(100000,170,20,1000,250,50);
  
  gaus2->SetParLimits(0, 20, 1000000);
  gaus2->SetParLimits(1, 150, 200);
  gaus2->SetParLimits(2, 30, 30);
  gaus2->SetParLimits(3, 20, 1000000);
  gaus2->SetParLimits(4, 170, 300);
  gaus2->SetParLimits(5, 50, 50);
  
  TH1F* hCombBkg;
  
  eventTree->Draw("hadTopMass >> hCombBkg(20, 100, 300)","(bProb*fitProb)*(target==0)");
  
  TH1F *hCombBkg = (TH1F*)gDirectory->Get("hCombBkg");
  
  //hCombBkg->Scale(1/hCombBkg->Integral());
  
  hCombBkg->Fit("gaus2","LEM");

  y0 [i] = gaus2->GetParameter(0)/gaus2->GetParameter(3);;
  //ey0[i] = gaus2->GetParError(0);
  ey0[i] = sqrt(pow(1/gaus2->GetParameter(3)*gaus2->GetParError(0), 2)
            + pow(gaus2->GetParameter(0)/pow(gaus2->GetParameter(3), 2)*gaus2->GetParError(3), 2));
  
  y1 [i] = gaus2->GetParameter(1);
  ey1[i] = gaus2->GetParError(1);
  
  y2 [i] = gaus2->GetParameter(2);
  ey2[i] = gaus2->GetParError(2);
  
  y3 [i] = gaus2->GetParameter(3);
  ey3[i] = gaus2->GetParError(3);
  
  y4 [i] = gaus2->GetParameter(4);
  ey4[i] = gaus2->GetParError(4);
  
  y5 [i] = gaus2->GetParameter(5);
  ey5[i] = gaus2->GetParError(5);

}
