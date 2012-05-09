#include <stdio.h>

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
#include "TString.h"
#include "TFitResult.h"
#include "TLine.h"

#include "tdrstyle.C"

TString sObs[] = {"m_{t}", "m_{W}^{raw}"};

enum styles          { kDown, kNominal, kUp};
int color_   [ 3 ] = { kRed+1, kBlue+1, kGreen+1};
int marker_  [ 3 ] = { 23, 20, 22};
int line_    [ 3 ] = { 7, 1, 9};

//TH1F* FindParameters(TString filename, int obs, int target);

namespace cb {
  double A(const double alpha, const double power) {
    return TMath::Power(power / TMath::Abs(alpha), power) * TMath::Exp(-alpha*alpha/2);
  };
  double B(const double alpha, const double power) {
    return power / TMath::Abs(alpha) - TMath::Abs(alpha);
  };
}

// 7 parameters: [0] -> [6]
double crystalBall(const double* xx, const double* p)
{
  double N     = p[0];
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
  
  double N1     = 1./sqrt(2*sigma1*sigma1);
  double N2     = 1./sqrt(2*sigma2*sigma2);
  
  N = (N1+N2)/2. * p[0];
  
  if(xx[0] < mu)
    return N * TMath::Exp(-t1*t1/2);
  else
    return N * TMath::Exp(-t2*t2/2);
}

TH1F* FindParameters(int iFile = 1, int obs = 0, int target = 1)
{
  setTDRStyle();
  TH1::SetDefaultSumw2();
  
  TFile* file;
  
  switch (iFile) {
    case 0: {
      file = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_1.00_muon/analyzeTop.root");
      break;
    }
    case 1: {
      file = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_P11_muon/analyzeTop.root");
      break;
    }
  }
  
  TTree* eventTree = (TTree*) file->Get("analyzeHitFit/eventTree");
  
  TF1* fit;
  
  if (obs == 0) {
    switch(target) {
      case   1: {
        fit = new TF1("fit", "[0]*TMath::Voigt(x-[1], [2], [3])", 100, 400);
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
        fit = new TF1("fit", crystalBall, 100, 400, 5);
        fit->SetLineColor(kBlack);
        fit->SetLineWidth(2);
        
        double power = 15.;
        
        fit->SetParNames("N", "#mu", "#sigma", "#alpha", "power");
        fit->SetParameters(1, 170, 25, 0.45, power);
        
        fit->SetParLimits(0, 0, 1000000);
        fit->SetParLimits(1, 150, 200);
        fit->SetParLimits(2, 15, 40);
        fit->SetParLimits(3, 0.3, 0.95);
        fit->SetParLimits(4, power, power);
        
        break;
      }
      case -10: {
        fit = new TF1("fit", crystalBall, 100, 400, 5);
        fit->SetLineColor(kBlack);
        fit->SetLineWidth(2);
        
        double power = 3.;
        
        fit->SetParNames("N", "#mu", "#sigma", "#alpha", "power");
        fit->SetParameters(1, 170, 15, 0.45, power);
        
        fit->SetParLimits(0, 0, 1000000);
        fit->SetParLimits(1, 150, 200);
        fit->SetParLimits(2, 10, 30);
        fit->SetParLimits(3, 0.05, 1.95);
        fit->SetParLimits(4, power, power);
        
        break;
      }
    }
  }
  
  else if (obs == 1) {
    fit = new TF1("fit", asymGaus, 0, 1000, 4);
    fit->SetLineColor(kBlack);
    fit->SetLineWidth(2);
    
    fit->SetParNames("N", "#mu", "#sigma1", "#sigma2");
    fit->SetParameters(1000, 80, 5, 5);
    
    fit->SetParLimits(0, 800, 1000000);
    fit->SetParLimits(1, 50, 150);
    fit->SetParLimits(2, 0, 10);
    fit->SetParLimits(3, 0, 15);
  }
  
  fit->SetNpx(300);
  
  TString sObservable;
  switch(obs) {
    case 0: {
      switch(target) {
        case   1: {
          sObservable = "hadTopMass >> h1(30, 100, 250)";
          break;
        }
        case   0: {
          sObservable = "hadTopMass >> h1(30, 100, 400)";
          break;
        }
        case -10: {
          sObservable = "hadTopMass >> h1(30, 100, 400)";
          break;
        }
      }
      break;
    }
    case 1: {
      sObservable = "hadWRawMass >> h1(30, 60, 120)";
      break;
    }
  }
  TString sCutAndWeight("(hitFitProb*MCWeight)*(target=="); sCutAndWeight += target; sCutAndWeight += " & hitFitProb > 0.2 & leptonPt > 30 & hadQBCSV<0.679 & hadQBarBCSV<0.679 & hadBBCSV>0.679 & lepBBCSV>0.679)";

  std::cout << sCutAndWeight << std::endl;
  
  // Get observable
  eventTree->Draw(sObservable, sCutAndWeight);
  
  TH1F *h1 = (TH1F*)gDirectory->Get("h1");
  
  double p[12];
  
  switch(obs) {
    case 0: {
      switch(target) {
        case   1: {
          h1->GetYaxis()->SetTitle("Fraction of entries / 5 GeV");
          double q[12] = {171.073, 0.996416, 86.5831, 0.778731, 9.69888, 0.0752241, 9.06047, -0.170575, 1, 1, 1, 1};
          for (int i = 0; i < 12; ++i) p[i] = q[i];
          break;
        }
        case   0: {
          h1->GetYaxis()->SetTitle("Fraction of entries / 10 GeV");
          double q[12] = {171.55, 0.857602, 107.067, -1.00712, 29.0538, 0.397468, 40.5463, -1.16334, 0.410231, 0.00464444, 0.360723, -0.0268234};
          for (int i = 0; i < 12; ++i) p[i] = q[i];
          break;
        }
        case -10: {
          h1->GetYaxis()->SetTitle("Fraction of entries / 10 GeV");
          double q[12] = {168.905, 0.902044, 88.8572, 0.454322, 18.5913, 0.185929, 11.1974, -0.143464, 0.860393, 0.00876846, 0.401755, -0.00365658};
          for (int i = 0; i < 12; ++i) p[i] = q[i];
          break;
        }
      }
      h1->GetXaxis()->SetTitle("m_{t,i}^{fit} [GeV]");
      break;
    }
    case 1: {
      h1->GetYaxis()->SetTitle("Fraction of entries / 2 GeV");
      h1->GetXaxis()->SetTitle("m_{W,i}^{reco} [GeV]");
      switch(target) {
        case   1: {
          double q[12] = {83.1044, 0.0105928, 51.3068, 0.179922, 6.17959, -0.0156169, 14.5686, 0.0232517, 7.07706, 0.000679047, -2.79713, -0.20092};
          for (int i = 0; i < 12; ++i) p[i] = q[i];
          break;
        }
        case   0: {
          double q[12] = {82.8984, -0.027674, 48.3411, -0.0550518, 6.21915, -0.0496059, 13.6744, -0.55865, 7.33281, 0.0217102, -4.74054, -0.161831};
          for (int i = 0; i < 12; ++i) p[i] = q[i];
          break;
        }
        case -10: {
          double q[12] = {82.6481, 0.0123294, 33.5792, 0.176464, 6.80263, 0.00166665, 13.651, 0.137787, 8.86199, 0.00685281, -10.6243, -0.309378};
          for (int i = 0; i < 12; ++i) p[i] = q[i];
          break;
        }
      }
      break;
    }
  }
  
  //*
  h1->SetFillColor(color_[target%10]);
  h1->SetLineColor(color_[target%10]);
  h1->SetMarkerColor(color_[target%10]);
  h1->SetMarkerStyle(marker_[target%10]);
  fit->SetLineColor(color_[target%10]);
  fit->SetLineStyle(line_[target%10]);
  //*/
  
  TFitResultPtr r = h1->Fit("fit","WLEMSR");
  
  std::cout << "chi2/ndf = " << r->Chi2() << "/" << r->Ndf() << std::endl;

  h1->Draw("E");
  
  double f[] = {fit->GetParameter(1), fit->GetParameter(2), fit->GetParameter(3)};
  double e[] = {fit->GetParError(1), fit->GetParError(2), fit->GetParError(3)};
  
  for (int i = 0; i < 3; ++i) {
    //std::cout << f[i] << " " << e[i] << " " << p[4*i+0] << " " << p[4*i+1] << " " << std::endl;
    if (obs == 0) {
      std::cout << i << ": " << (f[i] - p[4*i+0]) / p[4*i+1] + 172.5 << " +/- " << e[i] / p[4*i+1] << std::endl;
    }
    else {
      std::cout << i << ": " << (f[i] - p[4*i+0]) / p[4*i+2] + 1 << " +/- " << e[i] / p[4*i+2] << std::endl;
    }
  }
  
  return h1;
}
