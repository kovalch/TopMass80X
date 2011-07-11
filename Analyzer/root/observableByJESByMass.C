#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"

#include "tdrstyle.C"

int target = 1;
int obs    = 1; // 0: hadTopMass, 1: hadWRawMass

double X[3] = {166.5, 172.5, 178.5};
double Y10[3];
double Y11[3];
double Y20[3];
double Y21[3];
double Y30[3];
double Y31[3];

double eX[3] = {0.001, 0.001, 0.001};
double eY10[3];
double eY11[3];
double eY20[3];
double eY21[3];
double eY30[3];
double eY31[3];

double x[3] = {0.95, 1.0, 1.05};
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
  
  double N1     = 1./sqrt(2*sigma1*sigma1);
  double N2     = 1./sqrt(2*sigma2*sigma2);
  
  N = (N1+N2)/2. * p[0];
  
  if(x[0] < mu)
    return N * TMath::Exp(-t1*t1/2);
  else
    return N * TMath::Exp(-t2*t2/2);
}

void observableByJESByMass() {

  std::cout << "### FindParametersMass(0); ###" << std::endl;
  FindParametersMass(0);
  std::cout << "### FindParametersMass(1); ###" << std::endl;
  FindParametersMass(1);
  std::cout << "### FindParametersMass(2); ###" << std::endl;
  FindParametersMass(2);
  
  TF1* linearFit = new TF1("linearFit", "[0]+(x-172.5)*[1]");
  linearFit->SetLineWidth(2);
  linearFit->SetLineColor(kRed+1);
  
  TCanvas* cByMass = new TCanvas("cByMass", "cByMass", 800, 800);
  cByMass->Divide(2,3);
  
  cByMass->cd(1);
  
  gr = new TGraphErrors(3,X,Y10,eX,eY10);
  gr->SetTitle("10");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle("m_{t}");
  gr->GetYaxis()->SetTitle("#mu const");
  
  cByMass->cd(2);
  
  gr = new TGraphErrors(3,X,Y11,eX,eY11);
  gr->SetTitle("11");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle("m_{t}");
  gr->GetYaxis()->SetTitle("#mu slope");
  
  cByMass->cd(3);
  
  gr = new TGraphErrors(3,X,Y20,eX,eY20);
  gr->SetTitle("20");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle("m_{t}");
  gr->GetYaxis()->SetTitle("#sigma const");
  
  cByMass->cd(4);
  
  gr = new TGraphErrors(3,X,Y21,eX,eY21);
  gr->SetTitle("21");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle("m_{t}");
  gr->GetYaxis()->SetTitle("#sigma slope");
  
  cByMass->cd(5);
  
  gr = new TGraphErrors(3,X,Y30,eX,eY30);
  gr->SetTitle("30");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle("m_{t}");
  gr->GetYaxis()->SetTitle("#sigma2 const");
  
  cByMass->cd(6);
  
  gr = new TGraphErrors(3,X,Y31,eX,eY31);
  gr->SetTitle("31");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle("m_{t}");
  gr->GetYaxis()->SetTitle("#sigma2 slope");
}

void FindParametersMass(int iMass)
{
  setTDRStyle();
  gStyle->SetOptFit(0); 
    
  TCanvas* cObservable = new TCanvas("cObservable", "cObservable", 600, 600);
  
  cObservable->cd();
  
  TF1* linearFit = new TF1("linearFit", "[0]+(x-1)*[1]");
  linearFit->SetLineWidth(2);
  linearFit->SetLineColor(kRed+1);
  
  switch(iMass) {
    case 0: {
      TH1F* h095 = FindParameters("/scratch/hh/current/cms/user/mseidel/TTJets1665_0.95/analyzeTop.root", 0);
      TH1F* h100 = FindParameters("/scratch/hh/current/cms/user/mseidel/TTJets1665_1.0/analyzeTop.root", 1);
      TH1F* h105 = FindParameters("/scratch/hh/current/cms/user/mseidel/TTJets1665_1.05/analyzeTop.root", 2);
      break;
    }
    case 1: {
      TH1F* h095 = FindParameters("/scratch/hh/current/cms/user/mseidel/TTJets1725_0.95/analyzeTop.root", 0);
      TH1F* h100 = FindParameters("/scratch/hh/current/cms/user/mseidel/TTJets1725_1.0/analyzeTop.root", 1);
      TH1F* h105 = FindParameters("/scratch/hh/current/cms/user/mseidel/TTJets1725_1.05/analyzeTop.root", 2);
      break;
    }
    case 2: {
      TH1F* h095 = FindParameters("/scratch/hh/current/cms/user/mseidel/TTJets1785_0.95/analyzeTop.root", 0);
      TH1F* h100 = FindParameters("/scratch/hh/current/cms/user/mseidel/TTJets1785_1.0/analyzeTop.root", 1);
      TH1F* h105 = FindParameters("/scratch/hh/current/cms/user/mseidel/TTJets1785_1.05/analyzeTop.root", 2);
      break;
    }
  }
  
  h095->Draw();
  if (obs == 0) h095->GetXaxis()->SetRangeUser(100, 250);
  if (obs == 1) h095->GetXaxis()->SetRangeUser(50, 150);
  h100->Draw("SAME");
  h105->Draw("SAME");
  
  // ---
  //    create legend
  // ---
  TLegend *leg0 = new TLegend(0.65, 0.7, 0.95, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  leg0->AddEntry( h095 , "JES = 0.95", "F");
  leg0->AddEntry( h100 , "JES = 1.00", "F");
  leg0->AddEntry( h105 , "JES = 1.05", "F");
  
  leg0->Draw();
  
  gStyle->SetOptFit(1);
  
  TCanvas* cObservablePar = new TCanvas("cObservablePar", "cObservablePar", 1600, 400);
  cObservablePar->Divide(4,1);
  
  cObservablePar->cd(1);
  gr = new TGraphErrors(3,x,y1,ex,ey1);
  gr->SetTitle("p1");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle("JES");
  gr->GetYaxis()->SetTitle("#mu");
  
  Y10 [iMass] = linearFit->GetParameter(0);
  eY10[iMass] = linearFit->GetParError(0);
  
  Y11 [iMass] = linearFit->GetParameter(1);
  eY11[iMass] = linearFit->GetParError(1);
  
  
  cObservablePar->cd(2);
  gr = new TGraphErrors(3,x,y2,ex,ey2);
  gr->SetTitle("p2");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle("JES");
  gr->GetYaxis()->SetTitle("#sigma");
  
  Y20 [iMass] = linearFit->GetParameter(0);
  eY20[iMass] = linearFit->GetParError(0);
  
  Y21 [iMass] = linearFit->GetParameter(1);
  eY21[iMass] = linearFit->GetParError(1);
  
  
  cObservablePar->cd(3);
  gr = new TGraphErrors(3,x,y3,ex,ey3);
  gr->SetTitle("p3");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle("JES");
  if (obs == 0) gr->GetYaxis()->SetTitle("#alpha");
  if (obs == 1) gr->GetYaxis()->SetTitle("#sigma2");
  
  Y30 [iMass] = linearFit->GetParameter(0);
  eY30[iMass] = linearFit->GetParError(0);
  
  Y31 [iMass] = linearFit->GetParameter(1);
  eY31[iMass] = linearFit->GetParError(1);
  
  if (target!=1) {
    cObservablePar->cd(4);
    gr = new TGraphErrors(3,x,y4,ex,ey4);
    gr->SetTitle("p4");
    gr->Draw("A*");
    gr->Fit("linearFit", "EM");
    
    gr->GetXaxis()->SetTitle("JES");
    gr->GetYaxis()->SetTitle("power");
  }
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
  }
  
  else if (obs == 1) {
    fit = new TF1("fit", asymGaus, 0, 1000, 4);
    fit->SetLineColor(kBlack);
    fit->SetLineWidth(2);
    
    fit->SetParNames("N", "#mu", "#sigma1", "#sigma2");
    fit->SetParameters(1, 80, 5, 5);
    
    fit->SetParLimits(0, 0, 1000000);
    fit->SetParLimits(1, 50, 150);
    fit->SetParLimits(2, 0, 10);
    fit->SetParLimits(3, 0, 10);
  }

  
  TH1F* hSig;
  
  TString sObservable;
  switch(obs) {
    case 0: {
      sObservable = "hadTopMass >> h1(80, 100, 500)";
      break;
    }
    case 1: {
      sObservable = "hadWRawMass >> h1(80, 50, 150)";
      break;
    }
  }
  TString sCutAndWeight("(bProbSSV*hitFitProb)*(target=="); sCutAndWeight += target; sCutAndWeight += " & (bProbSSV * hitFitProb) > 0.05)";

  std::cout << sCutAndWeight << std::endl;
  
  // Get observable
  eventTree->Draw(sObservable, sCutAndWeight);
  
  TH1F *h1 = (TH1F*)gDirectory->Get("h1");
  
  h1->GetXaxis()->SetTitle("m_{i}");
  h1->GetYaxis()->SetTitle("Fraction of entries / 5 GeV");
  h1->SetFillColor(kRed-6-i);
  
  h1->Fit("fit","LEM");

  y0 [i] = fit->GetParameter(0);
  ey0[i] = fit->GetParError(0);
 
  y1 [i] = fit->GetParameter(1);
  ey1[i] = fit->GetParError(1);
  
  y2 [i] = fit->GetParameter(2);
  ey2[i] = fit->GetParError(2);
  
  y3 [i] = fit->GetParameter(3);
  ey3[i] = fit->GetParError(3);
  
  y4 [i] = fit->GetParameter(4);
  ey4[i] = fit->GetParError(4);
  
  TF1 *fitted = h1->GetFunction("fit");

  fitted->SetParameter(0, fitted->GetParameter(0)/h1->Integral());
  h1->Scale(1/h1->Integral());
  
  return h1;
}
