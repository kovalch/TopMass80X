#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"

#include "tdrstyle.C"

int target = 1;
int obs    = 4; // 0: hadTopMass, 1: hadWRawMass

bool plotByMass = false;

TString sObs[] = {"m_{t}", "m_{W}^{raw}"};

double X[9] = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5, 184.5};
double Y10[9];
double Y11[9];
double Y20[9];
double Y21[9];
double Y30[9];
double Y31[9];

double eX[9] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001};
double eY10[9];
double eY11[9];
double eY20[9];
double eY21[9];
double eY30[9];
double eY31[9];

double x[3] = {0.96, 1.0, 1.04};
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

#include <stdio.h>

enum styles          { kDown, kNominal, kUp};
int color_   [ 3 ] = { kRed+1, kBlue+1, kGreen+1};
int marker_  [ 3 ] = { 23, 20, 22};

void myflush ( FILE *in )
{
  int ch;

  do
    ch = fgetc ( in ); 
  while ( ch != EOF && ch != '\n' ); 

  clearerr ( in );
}

void mypause ( void ) 
{ 
  printf ( "Press [Enter] to continue . . ." );
  fflush ( stdout );
  getchar();
} 

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
  setTDRStyle();
  TH1::SetDefaultSumw2();

  for (int i=4; i<5; i++) {
	  std::cout << "### FindParametersMass(i); ###" << std::endl;
    FindParametersMass(i);
		/*
		mypause();
		myflush ( stdin );
		//*/
	}
  
  TF1* linearFit = new TF1("linearFit", "[0]+(x-172.5)*[1]");
  linearFit->SetLineWidth(2);
  linearFit->SetLineColor(kRed+1);
  
  TCanvas* cByMass = new TCanvas("cByMass", "cByMass", 800, 800);
  cByMass->Divide(2,3);
  
  cByMass->cd(1);
  
  gr = new TGraphErrors(9,X,Y10,eX,eY10);
  gr->SetTitle("10");
  gr->Draw("A*");
  //linearFit->SetParNames("p_{#mu}^{const}", "p_{#mu}^{m_{t}}");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle("m_{t}");
  gr->GetYaxis()->SetTitle("#mu const");
  
  cByMass->cd(2);
  
  gr = new TGraphErrors(9,X,Y11,eX,eY11);
  gr->SetTitle("11");
  gr->Draw("A*");
  //linearFit->SetParNames("p_{#mu}^{JES}", "p_{#mu}^{m_{t},JES}");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle("m_{t}");
  gr->GetYaxis()->SetTitle("#mu slope");
  
  cByMass->cd(3);
  
  gr = new TGraphErrors(9,X,Y20,eX,eY20);
  gr->SetTitle("20");
  gr->Draw("A*");
  //linearFit->SetParNames("p_{#sigma}^{const}", "p_{#sigma}^{m_{t}}");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle("m_{t}");
  switch(obs) {
    case 0: {
      gr->GetYaxis()->SetTitle("#sigma const");
      break;
    }
    case 1: {
      gr->GetYaxis()->SetTitle("#sigma_{1} const");
      break;
    }
  }
  
  cByMass->cd(4);
  
  gr = new TGraphErrors(9,X,Y21,eX,eY21);
  gr->SetTitle("21");
  gr->Draw("A*");
  //linearFit->SetParNames("p_{#sigma}^{JES}", "p_{#sigma}^{m_{t},JES}");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle("m_{t}");
  switch(obs) {
    case 0: {
      gr->GetYaxis()->SetTitle("#sigma slope");
      break;
    }
    case 1: {
      gr->GetYaxis()->SetTitle("#sigma_{1} slope");
      break;
    }
  }
  
  cByMass->cd(5);
  
  gr = new TGraphErrors(9,X,Y30,eX,eY30);
  gr->SetTitle("30");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle("m_{t}");
  switch(obs) {
    case 0: {
      gr->GetYaxis()->SetTitle("#alpha const");
      break;
    }
    case 1: {
      gr->GetYaxis()->SetTitle("#sigma_{2} const");
      break;
    }
  }
  
  cByMass->cd(6);
  
  gr = new TGraphErrors(9,X,Y31,eX,eY31);
  gr->SetTitle("31");
  gr->Draw("A*");
  gr->Fit("linearFit", "EM");
  
  gr->GetXaxis()->SetTitle("m_{t}");
  switch(obs) {
    case 0: {
      gr->GetYaxis()->SetTitle("#alpha slope");
      break;
    }
    case 1: {
      gr->GetYaxis()->SetTitle("#sigma_{2} slope");
      break;
    }
  }
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
	linearFit->SetParLimits(1, -50, 200);
	linearFit->SetParameters(100, 50);
  
  if (!plotByMass) {
    switch(iMass) {
      case 0: {
        TH1F* h096 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1615_0.96_2b/analyzeTop.root", 0);
        TH1F* h100 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1615_1.00_2b/analyzeTop.root", 1);
        TH1F* h104 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1615_1.04_2b/analyzeTop.root", 2);
        break;
      }
      case 1: {
        TH1F* h096 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1635_0.96_2b/analyzeTop.root", 0);
        TH1F* h100 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1635_1.00_2b/analyzeTop.root", 1);
        TH1F* h104 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1635_1.04_2b/analyzeTop.root", 2);
        break;
      }
		  case 2: {
        TH1F* h096 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1665_0.96_2b/analyzeTop.root", 0);
        TH1F* h100 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1665_1.00_2b/analyzeTop.root", 1);
        TH1F* h104 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1665_1.04_2b/analyzeTop.root", 2);
        break;
      }
		  case 3: {
        TH1F* h096 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1695_0.96_2b/analyzeTop.root", 0);
        TH1F* h100 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1695_1.00_2b/analyzeTop.root", 1);
        TH1F* h104 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1695_1.04_2b/analyzeTop.root", 2);
        break;
      }
		  case 4: {
        TH1F* h096 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_0.96_2b/analyzeTop.root", 0);
        TH1F* h100 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_1.00_2b/analyzeTop.root", 1);
        TH1F* h104 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_1.04_2b/analyzeTop.root", 2);
        break;
      }
		  case 5: {
        TH1F* h096 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1755_0.96_2b/analyzeTop.root", 0);
        TH1F* h100 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1755_1.00_2b/analyzeTop.root", 1);
        TH1F* h104 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1755_1.04_2b/analyzeTop.root", 2);
        break;
      }
		  case 6: {
        TH1F* h096 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1785_0.96_2b/analyzeTop.root", 0);
        TH1F* h100 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1785_1.00_2b/analyzeTop.root", 1);
        TH1F* h104 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1785_1.04_2b/analyzeTop.root", 2);
        break;
      }
		  case 7: {
        TH1F* h096 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1815_0.96_2b/analyzeTop.root", 0);
        TH1F* h100 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1815_1.00_2b/analyzeTop.root", 1);
        TH1F* h104 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1815_1.04_2b/analyzeTop.root", 2);
        break;
      }
		  case 8: {
        TH1F* h096 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1845_0.96_2b/analyzeTop.root", 0);
        TH1F* h100 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1845_1.00_2b/analyzeTop.root", 1);
        TH1F* h104 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1845_1.04_2b/analyzeTop.root", 2);
        break;
      }
    }
  }
  else {
    TH1F* h096 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1665_1.00_2b/analyzeTop.root", 0);
    TH1F* h100 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_1.00_2b/analyzeTop.root", 1);
    TH1F* h104 = FindParameters("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1785_1.00_2b/analyzeTop.root", 2); 
  }
  
  h096->Draw();
  if (obs == 0) h096->GetXaxis()->SetRangeUser(100, 250);
  if (obs == 1) h096->GetXaxis()->SetRangeUser(50, 150);
  h100->Draw("SAME");
  h104->Draw("SAME");
  
  // ---
  //    create legend
  // ---
  TLegend *leg0 = new TLegend(0.65, 0.7, 0.95, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  if (!plotByMass) {
    leg0->AddEntry( h096 , "JES = 0.96", "F");
    leg0->AddEntry( h100 , "JES = 1.00", "F");
    leg0->AddEntry( h104 , "JES = 1.04", "F");
  }
  else {
    leg0->AddEntry( h096 , "m_{t,gen} = 166.5 GeV", "F");
    leg0->AddEntry( h100 , "m_{t,gen} = 172.5 GeV", "F");
    leg0->AddEntry( h104 , "m_{t,gen} = 178.5 GeV", "F");
  }
  
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
        fit->SetParLimits(2, 15, 40);
        fit->SetParLimits(3, 0.3, 0.95);
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
        fit->SetParLimits(3, 0.05, 1.95);
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
    fit->SetParameters(1000, 80, 5, 5);
    
    fit->SetParLimits(0, 800, 1000000);
    fit->SetParLimits(1, 50, 150);
    fit->SetParLimits(2, 0, 10);
    fit->SetParLimits(3, 0, 10);
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
    fit->SetParLimits(3, 0, 2);
    fit->SetParLimits(4, 0, 20);
  }
  
  else if (obs == 5) {
    fit = new TF1("fit", "[0]*TMath::Voigt(x-[1], [2], [3])");

    fit->SetLineColor(kBlack);
    fit->SetLineWidth(2);
    fit->SetParameters(100000, 0.0, 0.001, 1.35);
    
    fit->SetParLimits(0, 1, 1000000);
    fit->SetParLimits(1, 0, 2);
    fit->SetParLimits(2, 0.001, 0.001);
    fit->SetParLimits(3, 0, 10);
  }

  
  TH1F* hSig;
  
  TString sObservable;
  switch(obs) {
    case 0: {
      sObservable = "hadTopMass >> h1(250, 100, 500)";
      break;
    }
    case 1: {
      sObservable = "hadWRawMass >> h1(100, 50, 150)";
      break;
    }
    case 4: {
      sObservableShort = "#DeltaR";
      sObservable = "deltaRHadWHadB/deltaRHadQHadQBar >> h1(40, 0, 4)";
      break;
    }
    case 5: {
      sObservable = "-(hadWPt-lepWPt)/(hadBPt-lepBPt) >> h1(50, -5, 5)";
      break;
    }
  }
  TString sCutAndWeight("(hitFitProb*MCWeight)*(target=="); sCutAndWeight += target; sCutAndWeight += " & (hitFitProb) > 0.2)";

  std::cout << sCutAndWeight << std::endl;
  
  // Get observable
  eventTree->Draw(sObservable, sCutAndWeight);
  
  TH1F *h1 = (TH1F*)gDirectory->Get("h1");
  
  h1->GetXaxis()->SetTitle("m_{i}");
  //h1->GetXaxis()->SetTitle("#DeltaR^{decay}_{t}/#DeltaR^{decay}_{W}");
  h1->GetYaxis()->SetTitle("Fraction of entries");
  h1->SetFillColor(color_[i]);
  h1->SetLineColor(color_[i]);
  h1->SetMarkerColor(color_[i]);
  h1->SetMarkerStyle(2);
  fit->SetLineColor(color_[i]);
  
  h1->Fit("fit","WLEM");

  y0 [i] = fit->GetParameter(0);
  ey0[i] = fit->GetParError(0);
 
  y1 [i] = fit->GetParameter(1);
  ey1[i] = fit->GetParError(1);
  
  y2 [i] = fit->GetParameter(2);
  ey2[i] = fit->GetParError(2);
  if (ey2[i] == 0) ey2[i] = 0.001;
  
  y3 [i] = fit->GetParameter(3);
  ey3[i] = fit->GetParError(3);
  
  y4 [i] = fit->GetParameter(4);
  ey4[i] = fit->GetParError(4);
  
  TF1 *fitted = h1->GetFunction("fit");

  fitted->SetParameter(0, fitted->GetParameter(0)/h1->Integral());
  h1->Scale(1/h1->Integral());
  
  return h1;
}
