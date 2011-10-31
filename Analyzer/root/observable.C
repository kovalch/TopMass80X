#include <iostream>

#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"

#include "tdrstyle.C"

int target = 1;
int obs    = 4; // 0: hadTopMass, 1: hadWRawMass, 2: hadTopPt-lepTopPt
int test   = 0; // 0: mt,gen, 1: JES, 2: bJES

TString sX0[3];
double x0[3];
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

TString sObservable;
TString sObservableShort;

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
  
  if(t < alpha)
    return N * TMath::Exp(-t*t/2);
  else
    return N * cb::A(alpha,power) * TMath::Power(cb::B(alpha,power) + t, -power);
}

double asymGaus(const double* x, const double* p)
{
  double N      = p[0];
  double mu     = p[1];
  double sigma1 = p[2];
  double sigma2 = p[3];
  double t1     = (x[0] - mu) / sigma1;
  double t2     = (x[0] - mu) / sigma2;
  
  double N1     = sqrt(2*sigma1*sigma1)/2.;
  double N2     = sqrt(2*sigma2*sigma2)/2.;
  
  N = N1+N2;
  
  if(t1 < 0)
    return 1./N * p[0] * TMath::Exp(-t1*t1/2);
  else
    return 1./N * p[0] * TMath::Exp(-t2*t2/2);
}

void observable()
{
  setTDRStyle();
  gStyle->SetOptFit(0); 
  
  TH1::SetDefaultSumw2();
    
  TCanvas* cObservable = new TCanvas("cObservable", "cObservable", 600, 600);
  
  cObservable->cd();
  
  TF1* linearFit;
  switch(test) {
    case 0: {
      linearFit = new TF1("linearFit", "[0]+(x-172.5)*[1]");
      break;
    }
    default: {
      linearFit = new TF1("linearFit", "[0]+(x-1)*[1]");
    }
  }
  linearFit->SetLineWidth(2);
  linearFit->SetLineColor(kRed+1);
  
  switch(obs) {
    case 0: {
      sObservableShort = "m_{t}";
      sObservable = "hadTopMass >> h1(400, 100, 500)";
      break;
    }
    case 1: {
      sObservableShort = "m_{W}";
      sObservable = "hadWRawMass >> h1(120, 60, 120)";
      break;
    }
    case 2: {
      sObservableShort = "A_{p_{T}}";
      sObservable = "hadTopPt-lepTopPt >> h1(25, -100, 100)";
      break;
    }
    case 3: {
      sObservableShort = "A_{p_{T}}";
      sObservable = "(hadTopPt-lepTopPt)/(hadTopPt+lepTopPt) >> h1(1000, -0.2, 0.2)";
      break;
    }
    case 4: {
      sObservableShort = "#DeltaR";
      sObservable = "deltaRHadWHadB/deltaRHadQHadQBar >> h1(40, 0, 4)";
      break;
    }
    case 5: {
      sObservableShort = "S'";
      sObservable = "-(hadWPt-lepWPt)/(hadBPt-lepBPt) >> h1(800, -3, 5)";
      break;
    }
    case 6: {
      sObservableShort = "sin(#theta*)";
      sObservable = "( ( sqrt( pow(sin(deltaThetaHadWHadB),2)*pow(172.5,2) / pow(1/2*(pow(172.5,2)-pow(80.4,2))/(hadWE-sqrt(pow(hadWE,2) - pow(80.4,2))) - (hadWE-sqrt(pow(hadWE,2) - pow(80.4,2))*cos(deltaThetaHadWHadB))*cos(deltaThetaHadWHadB),2) + 1/pow(((hadWE-(sqrt(pow(hadWE,2) - pow(80.4,2)))*cos(deltaThetaHadWHadB)))/172.5,2)) ) ) >> h1(25, 0, 5)";
      break;
    }
  }
  
  TH1F* h0;
  TH1F* h1;
  TH1F* h2;
  
  switch(test) {
    case 0: {
      h0 = FindParameters("/scratch/hh/current/cms/user/mseidel/Spring11_TTJets1665_1.00_1.00_2b/analyzeTop.root", 0);
      h1 = FindParameters("/scratch/hh/current/cms/user/mseidel/Spring11_TTJets1725_1.00_1.00_2b/analyzeTop.root", 1);
      h2 = FindParameters("/scratch/hh/current/cms/user/mseidel/Spring11_TTJets1785_1.00_1.00_2b/analyzeTop.root", 2);
      break;
    }
    case 1: {
      h0 = FindParameters("/scratch/hh/current/cms/user/mseidel/Spring11_TTJets1725_0.96_1.00_2b/analyzeTop.root", 0);
      h1 = FindParameters("/scratch/hh/current/cms/user/mseidel/Spring11_TTJets1725_1.00_1.00_2b/analyzeTop.root", 1);
      h2 = FindParameters("/scratch/hh/current/cms/user/mseidel/Spring11_TTJets1725_1.04_1.00_2b/analyzeTop.root", 2);
      break;
    }
    case 2: {
      h0 = FindParameters("/scratch/hh/current/cms/user/mseidel/Spring11_TTJets1725_1.00_0.96/analyzeTop.root", 0);
      h1 = FindParameters("/scratch/hh/current/cms/user/mseidel/Spring11_TTJets1725_1.00_1.00/analyzeTop.root", 1);
      h2 = FindParameters("/scratch/hh/current/cms/user/mseidel/Spring11_TTJets1725_1.00_1.04/analyzeTop.root", 2);
    }
  }
  
  switch(test) {
    case 0: {
      sX0[0] = " = 166.5 GeV";
      sX0[1] = " = 172.5 GeV";
      sX0[2] = " = 178.5 GeV";
      
      x0[0] = 166.5;
      x0[1] = 172.5;
      x0[2] = 178.5;
      break;
    }
    default: {
      sX0[0] = " -4%";
      sX0[1] = "";
      sX0[2] = " +4%";
      
      x0[0] = 0.96;
      x0[1] = 1.00;
      x0[2] = 1.04;
    }
  }
  
  TString sTest;
  switch(test) {
    case 0: {
      sTest = "m_{t,gen}";
      break;
    }
    case 1: {
      sTest = "JES";
      break;
    }
    case 2: {
      sTest = "bJES";
      break;
    }
  }
  
  /* ISR/FSR test
  TH1F* h096 = FindParameters("/scratch/hh/current/cms/user/mseidel/Fall10_TTJets1725RS_1.00_1.00/analyzeTop.root", 0);
  TH1F* h100 = FindParameters("/scratch/hh/current/cms/user/mseidel/Fall10_TTJets1725_1.00_1.00/analyzeTop.root", 1);
  TH1F* h104 = FindParameters("/scratch/hh/current/cms/user/mseidel/Fall10_TTJets1725RL_1.00_1.00/analyzeTop.root", 2);
  //*/
  
  h0->Draw();
  if (obs == 0) h0->GetXaxis()->SetRangeUser(100, 250);
  if (obs == 1) h0->GetXaxis()->SetRangeUser(20, 160);
  h1->Draw("SAME");
  h2->Draw("SAME");
  
  // ---
  //    create legend
  // ---
  TLegend *leg0 = new TLegend(0.65, 0.7, 0.95, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  
  TString sLeg0(sTest); sLeg0 += sX0[0];
  TString sLeg1(sTest); sLeg1 += sX0[1];
  TString sLeg2(sTest); sLeg2 += sX0[2];
  
  leg0->AddEntry( h0, sLeg0, "F");
  leg0->AddEntry( h1, sLeg1, "F");
  leg0->AddEntry( h2, sLeg2, "F");
  
  leg0->Draw();
  
  gStyle->SetOptFit(1);
  
  TCanvas* cObservablePar = new TCanvas("cObservablePar", "cObservablePar", 1600, 400);
  cObservablePar->Divide(4,1);
  
  cObservablePar->cd(1);
  gr = new TGraphErrors(3,x0,y1,ex,ey1);
  gr->SetTitle("p1");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle(sTest);
  gr->GetYaxis()->SetTitle("#mu");
  
  cObservablePar->cd(2);
  gr = new TGraphErrors(3,x0,y2,ex,ey2);
  gr->SetTitle("p2");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle(sTest);
  gr->GetYaxis()->SetTitle("#sigma");
  
  cObservablePar->cd(3);
  gr = new TGraphErrors(3,x0,y3,ex,ey3);
  gr->SetTitle("p3");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle(sTest);
  if (obs == 0) gr->GetYaxis()->SetTitle("#alpha");
  if (obs == 1) gr->GetYaxis()->SetTitle("#sigma2");
	gr->GetYaxis()->SetTitle("#alpha");
  
  //if (target!=1 || obs==2) {
    cObservablePar->cd(4);
    gr = new TGraphErrors(3,x0,y4,ex,ey4);
    gr->SetTitle("p4");
    gr->Draw("A*");
    gr->Fit("linearFit", "EM");
    
  gr->GetXaxis()->SetTitle(sTest);
    gr->GetYaxis()->SetTitle("power");
  //}
}

TH1F* FindParameters(TString filename, int i)
{

  TFile* file = new TFile(filename);
  analyzeHitFit->cd();
  
  TF1* fit;
  
  if (obs == 0) {
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
        fit->SetParLimits(2, 15, 50);
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
  }
  
  else if (obs == 1) {
    fit = new TF1("fit", asymGaus, -4, 4, 4);
    fit->SetLineColor(kBlack);
    fit->SetLineWidth(2);
    
    fit->SetParNames("N", "#mu", "#sigma1", "#sigma2");
    fit->SetParameters(1, 80, 5, 5);
    
    fit->SetParLimits(0, 0, 1000000);
    fit->SetParLimits(1, 60, 100);
    fit->SetParLimits(2, 0, 10);
    fit->SetParLimits(3, 0, 10);
  }
  
  else if (obs == 2) {
    /*
    fit = new TF1("fit", asymGaus, -4, 4, 4);
    fit->SetLineColor(kBlack);
    fit->SetLineWidth(2);
    
    fit->SetParNames("N", "#mu", "#sigma1", "#sigma2");
    fit->SetParameters(1, 0, 5, 5);
    
    fit->SetParLimits(0, 0, 1000000);
    fit->SetParLimits(1, -10, 20);
    fit->SetParLimits(2, 0, 30);
    fit->SetParLimits(3, 0, 30);
    */
    
    fit = new TF1("fit", crystalBall, -20, 20, 5);
    fit->SetLineColor(kBlack);
    fit->SetLineWidth(2);
    
    fit->SetParNames("N", "#mu", "#sigma", "#alpha", "power");
    fit->SetParameters(1, 0, 25, 0.45, 3);
    
    fit->SetParLimits(0, 0, 1000000);
    fit->SetParLimits(1, -10, 20);
    fit->SetParLimits(2, 10, 30);
    fit->SetParLimits(3, 0.05, 2);
    fit->SetParLimits(4, 3, 3);
  }
  
  else if (obs == 3) {
    fit = new TF1("fit", "[0]*TMath::Voigt(x-[1], [2], [3])");
    //fit = new TF1("fit", "gaus");
    fit->SetLineColor(kBlack);
    fit->SetLineWidth(2);
    fit->SetParameters(100000, 0.5, 0.1, 0.09);
    
    fit->SetParLimits(0, 1, 1000000);
    fit->SetParLimits(1, -0.1, 0.1);
    fit->SetParLimits(2, 0.02, 0.2);
    fit->SetParLimits(3, 0.01, 0.5);
  }
  
  else if (obs == 4) {
    fit = new TF1("fit", crystalBall, 0, 1000, 5);
    fit->SetLineColor(kBlack);
    fit->SetLineWidth(2);
    
    fit->SetParNames("N", "#mu", "#sigma", "#alpha", "power");
    fit->SetParameters(1, 0, 1, 2, 3);
    
    fit->SetParLimits(0, 0, 1000000);
    fit->SetParLimits(1, 0.5, 2);
    fit->SetParLimits(2, 0.2, 2);
    fit->SetParLimits(3, 2, 2);
    fit->SetParLimits(4, 0, 20);
  }
  
  else if (obs == 5) {
    fit = new TF1("fit", "[0]*TMath::Voigt(x-[1], [2], [3])");

    fit->SetLineColor(kBlack);
    fit->SetLineWidth(2);
    fit->SetParameters(100000, 0.0, 0.001, 1.55);
    
    fit->SetParLimits(0, 1, 1000000);
    fit->SetParLimits(1, 0, 2);
    fit->SetParLimits(2, 0.001, 0.001);
    fit->SetParLimits(3, 1., 1.5);
  }
  
  else if (obs == 6) {
    fit = new TF1("fit", crystalBall, 0, 1000, 5);
    fit->SetLineColor(kBlack);
    fit->SetLineWidth(2);
    
    fit->SetParNames("N", "#mu", "#sigma", "#alpha", "power");
    fit->SetParameters(1, 0, 1, 0.45, 10);
    
    fit->SetParLimits(0, 0, 1000000);
    fit->SetParLimits(1, 0.5, 2);
    fit->SetParLimits(2, 0.2, 2);
    fit->SetParLimits(3, 0.05, 10);
    fit->SetParLimits(4, 10, 10);
  }

  
  TH1F* hSig;
  

  TString sCutAndWeight("(hitFitProb*MCWeight)*(target=="); sCutAndWeight += target; sCutAndWeight += " & (hitFitProb) > 0.2)";
  //TString sCutAndWeight("(bProbSSV*hitFitProb*PUWeight)*(target=="); sCutAndWeight += target; sCutAndWeight += " & (bProbSSV*hitFitProb) > 0.05 & sumBPt > 50)";
  
  std::cout << sObservableShort << std::endl;
  std::cout << sObservable << std::endl;
  std::cout << sCutAndWeight << std::endl;
  
  // Get observable
  eventTree->Draw(sObservable, sCutAndWeight);
  
  TH1F *h1 = (TH1F*)gDirectory->Get("h1");
  
  h1->GetXaxis()->SetTitle(sObservableShort);
  h1->GetYaxis()->SetTitle("Fraction of entries");
  h1->SetFillColor(kRed-4-i);
  h1->SetLineColor(kRed-4-i);
  h1->SetMarkerColor(kRed-4-i);
  h1->SetMarkerStyle(2);
  //fit->SetLineColor(kRed-4-i);
  
  h1->Fit("fit","WLEM");

  y0 [i] = fit->GetParameter(0);
  ey0[i] = fit->GetParError(0);
 
  y1 [i] = fit->GetParameter(1);
  ey1[i] = fit->GetParError(1);
  
  y2 [i] = fit->GetParameter(2);
  ey2[i] = 0.01;
  
  y3 [i] = fit->GetParameter(3);
  ey3[i] = fit->GetParError(3);
  
  y4 [i] = fit->GetParameter(4);
  ey4[i] = fit->GetParError(4);
  
  TF1 *fitted = h1->GetFunction("fit");

  fitted->SetParameter(0, fitted->GetParameter(0)/h1->Integral());
  h1->Scale(1/h1->Integral());
  
  return h1;
}
