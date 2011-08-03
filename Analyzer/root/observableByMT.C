#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"

#include "tdrstyle.C"

int target = 1;

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
  double M     = p[0];
  double mu    = p[1];
  double sigma = p[2];
  double alpha = p[3];
  double power = p[4];
  double t = (x[0] - mu) / sigma;
  
  double N1 = -sqrt(TMath::PiOver2()) * sigma * (TMath::Erf(-alpha/sqrt(2)) - TMath::Erf(mu/(sqrt(2)*sigma)));
  double N2 = cb::A(alpha,power) / (-1.+power) * (
                (-cb::B(alpha,power)*sigma+mu-10000) * TMath::Power(cb::B(alpha,power)+(-mu+10000)/sigma, -power)
               -(-cb::B(alpha,power)*sigma-sigma*alpha) * TMath::Power(cb::B(alpha,power)+(sigma*alpha)/sigma, -power)
              );
  double N = N1 + N2;
  
  if(t < alpha)
    return M * 1./N * TMath::Exp(-t*t/2);
  else
    return M * 1./N * cb::A(alpha,power) * TMath::Power(cb::B(alpha,power) + t, -power);
}

void observableByMT()
{
  setTDRStyle();
  gStyle->SetOptFit(0); 
    
  TCanvas* cObservable = new TCanvas("cObservable", "cObservable", 600, 600);
  
  cObservable->cd();
  
  TF1* linearFit = new TF1("linearFit", "[0]+(x-172.5)*[1]");
  linearFit->SetLineWidth(2);
  linearFit->SetLineColor(kRed+1);
  
  TH1F* h1665 = FindParameters("/scratch/hh/current/cms/user/mseidel/Spring11_TTJets1665_1.00_1.00/analyzeTop.root", 0);
  TH1F* h1725 = FindParameters("/scratch/hh/current/cms/user/mseidel/Spring11_TTJets1725_1.00_1.00/analyzeTop.root", 1);
  TH1F* h1785 = FindParameters("/scratch/hh/current/cms/user/mseidel/Spring11_TTJets1785_1.00_1.00/analyzeTop.root", 2);
  
  h1665->Draw();
  h1665->GetXaxis()->SetRangeUser(100, 250);
  h1725->Draw("SAME");
  h1785->Draw("SAME");
  
  // ---
  //    create legend
  // ---
  TLegend *leg0 = new TLegend(0.65, 0.7, 0.95, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  leg0->AddEntry( h1665 , "m_{t,gen} = 166.5 GeV", "F");
  leg0->AddEntry( h1725 , "m_{t,gen} = 172.5 GeV", "F" );
  leg0->AddEntry( h1785 , "m_{t,gen} = 178.5 GeV", "F" );
  
  leg0->Draw();
  
  gStyle->SetOptFit(1);
  
  TCanvas* cObservablePar = new TCanvas("cObservablePar", "cObservablePar", 1200, 400);
  cObservablePar->Divide(3,1);
  
  cObservablePar->cd(1);
  gr = new TGraphErrors(3,x,y1,ex,ey1);
  gr->SetTitle("p1");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle("m_{t,gen}");
  gr->GetYaxis()->SetTitle("#mu");
  
  cObservablePar->cd(2);
  gr = new TGraphErrors(3,x,y2,ex,ey2);
  gr->SetTitle("p2");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle("m_{t,gen}");
  gr->GetYaxis()->SetTitle("#sigma");
  
  cObservablePar->cd(3);
  gr = new TGraphErrors(3,x,y3,ex,ey3);
  gr->SetTitle("p3");
  gr->Draw("A*");
  gr->Fit("linearFit", "M");
  
  gr->GetXaxis()->SetTitle("m_{t,gen}");
  gr->GetYaxis()->SetTitle("#alpha");
}

TH1F* FindParameters(TString filename, int i)
{

  TFile* file = new TFile(filename);
  analyzeHitFit->cd();
  
  TF1* fit;
  
  switch(target) {
    case   1: {
      fit = new TF1("fit", "[0]*TMath::Voigt(x-[1], [2], [3])");
      fit->SetLineColor(kBlack);
      fit->SetLineWidth(2);
      fit->SetParameters(100000, 170, 10, 2);
      
      fit->SetParLimits(0, 1, 1000000);
      fit->SetParLimits(1, 150, 200);
      fit->SetParLimits(2, 5, 20);
      fit->SetParLimits(3, 2, 2);
      
      break;
    }
    case   0: {
      fit = new TF1("fit", crystalBall, 0, 1000, 5);
      fit->SetLineColor(kBlack);
      fit->SetLineWidth(2);
      
      fit->SetParNames("N", "#mu", "#sigma", "#alpha", "power");
      fit->SetParameters(1, 170, 25, 0.45, 15);
      
      fit->SetParLimits(0, 0, 1000000);
      fit->SetParLimits(1, 150, 200);
      fit->SetParLimits(2, 15, 30);
      fit->SetParLimits(3, 0.05, 0.95);
      fit->SetParLimits(4, 15, 15);
      
      break;
    }
    case -10: {
      fit = new TF1("fit", crystalBall, 0, 1000, 5);
      fit->SetLineColor(kBlack);
      fit->SetLineWidth(2);
      
      fit->SetParNames("N", "#mu", "#sigma", "#alpha", "power");
      fit->SetParameters(1, 170, 25, 0.45, 5);
      
      fit->SetParLimits(0, 0, 1000000);
      fit->SetParLimits(1, 150, 200);
      fit->SetParLimits(2, 15, 30);
      fit->SetParLimits(3, 0.05, 0.95);
      fit->SetParLimits(4, 5, 5);
      
      break;
    }
  }

  
  TH1F* hSig;
  
  TString sCutAndWeight("(bProbSSV*hitFitProb)*(target=="); sCutAndWeight += target; sCutAndWeight += " & (bProbSSV * hitFitProb) > 0.05)";
  //sCutAndWeight = "target==1";
  
  // Get observable
  eventTree->Draw("hadTopMass >> hSig(80, 100, 500)", sCutAndWeight);
  
  TH1F *hSig = (TH1F*)gDirectory->Get("hSig");
  
  hSig->GetXaxis()->SetTitle("m_{i}");
  hSig->GetYaxis()->SetTitle("Fraction of entries / 5 GeV");
  hSig->SetFillColor(kRed-8-i);
  
  hSig->Fit("fit","LEM");

  y0 [i] = fit->GetParameter(0);
  ey0[i] = fit->GetParError(0);
 
  y1 [i] = fit->GetParameter(1);
  ey1[i] = fit->GetParError(1);
  
  y2 [i] = fit->GetParameter(2);
  ey2[i] = fit->GetParError(2);
  
  y3 [i] = fit->GetParameter(3);
  ey3[i] = fit->GetParError(3);
  
  TF1 *fitted = hSig->GetFunction("fit");

  fitted->SetParameter(0, fitted->GetParameter(0)/hSig->Integral());
  hSig->Scale(1/hSig->Integral());
  
  return hSig;
}
