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
double y6[3];
double ex[3] = {0.001, 0.001, 0.001};
double ey0[3];
double ey1[3];
double ey2[3];
double ey3[3];
double ey4[3];
double ey5[3];
double ey6[3];

namespace cb {
  double A(const double alpha, const double power) {
    return TMath::Power(power / TMath::Abs(alpha), power) * TMath::Exp(-alpha*alpha/2);
  };
  double B(const double alpha, const double power) {
    return power / TMath::Abs(alpha) - TMath::Abs(alpha);
  };
}

// 7 parameters: [0] -> [6]
double crystalBall(const double* x, const double* p)
{
  double N     = p[0];
  double mu    = p[1];
  double sigma = p[2];
  double alpha = p[3];
  double power = p[4];
  double t = (x[0] - mu) / sigma;
  
  /*
  double N1 = -sqrt(TMath::PiOver2()) * sigma * (-1+TMath::Erf(-sigma*alpha)/(sqrt(2)*sigma));
  std::cout << "N1: " << N1 << std::endl;
  double N2 = cb::A(alpha,power) * sigma/(-1.+power) * ( 0. + //lim(x->inf)
                pow(((cb::B(alpha,power) * sigma - sigma*alpha)/sigma), (1. - power))
              );
  std::cout << "N2: " << N2 << std::endl;
  N = N1 + N2;
  //*/
  
  if(t < alpha)
    return N * TMath::Exp(-t*t/2);
  else
    return N * cb::A(alpha,power) * TMath::Power(cb::B(alpha,power) + t, -power);
}


void ideogramCombBkgCB()
{
  // open input file/directory
//  file1725 = new TFile("TopMass/Analyzer/root/analyzeTop_1725.root");
//  file1785 = new TFile("TopMass/Analyzer/root/analyzeTop_1785.root");
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1);
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 900, 900);
  canvas->Divide(3,3);
  
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
  
  canvas->cd(6);
  gr = new TGraphErrors(3,x,y2,ex,ey2);
  gr->SetTitle("p2");
  gr->Draw("A*");
  gr->Fit("pol1", "M");
  
  canvas->cd(7);
  gr = new TGraphErrors(3,x,y3,ex,ey3);
  gr->SetTitle("p3");
  gr->Draw("A*");
  gr->Fit("pol1", "M");
  
  canvas->cd(8);
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

  TF1* cb = new TF1("cb", crystalBall, 0, 1000, 5);
  cb->SetLineColor(kRed+1);
  
  cb->SetParameters(1, 170, 25, 0.45, 15);
  
  cb->SetParLimits(0, 0, 1000000);
  cb->SetParLimits(1, 150, 200);
  cb->SetParLimits(2, 20, 30);
  cb->SetParLimits(3, 0.05, 0.95);
  cb->SetParLimits(4, 5, 5);
  
  TH1F* hCombBkg;
  
  eventTree->Draw("hadTopMass >> hCombBkg(80, 100, 500)", "(bProbSSV*hitFitProb)*((target==0 | target==-2) & (bProbSSV * hitFitProb) > 0.01)");
  
  TH1F *hCombBkg = (TH1F*)gDirectory->Get("hCombBkg");
  
  double integral = hCombBkg->Integral();
  integral = 1;
  hCombBkg->Scale(1/integral);
  
  hCombBkg->Fit("cb","WEM");

  y0 [i] = cb->GetParameter(0);
  ey0[i] = cb->GetParError(0)/sqrt(integral);
  /*
  ey0[i] = sqrt(pow(1/cb->GetParameter(3)*cb->GetParError(0), 2)
            + pow(cb->GetParameter(0)/pow(cb->GetParameter(3), 2)*cb->GetParError(3), 2));
  */
  
  y1 [i] = cb->GetParameter(1);
  ey1[i] = cb->GetParError(1)/sqrt(integral);
  
  y2 [i] = cb->GetParameter(2);
  ey2[i] = cb->GetParError(2)/sqrt(integral);
  
  y3 [i] = cb->GetParameter(3);
  ey3[i] = cb->GetParError(3)/sqrt(integral);
  
  y4 [i] = cb->GetParameter(4);
  ey4[i] = cb->GetParError(4)/sqrt(integral);

}
