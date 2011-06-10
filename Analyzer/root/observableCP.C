#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"

#include "tdrstyle.C"

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

void observableCP()
{
  setTDRStyle();
  gStyle->SetOptFit(0); 
    
  TCanvas* cObservableCP = new TCanvas("cObservableCP", "cObservableCP", 600, 600);
  
  cObservableCP->cd();
  
  TH1F* h1665 = FindParameters("analyzeTop_1665.root", 0);
  TH1F* h1725 = FindParameters("analyzeTop_1725.root", 1);
  TH1F* h1785 = FindParameters("analyzeTop_1785.root", 2);
  
  h1665->Draw();
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
  
  TCanvas* cObservableCPPar = new TCanvas("cObservableCPPar", "cObservableCPPar", 1200, 600);
  cObservableCPPar->Divide(2,1);
  
  cObservableCPPar->cd(1);
  gr = new TGraphErrors(3,x,y1,ex,ey1);
  gr->SetTitle("p1");
  gr->Draw("A*");
  gr->Fit("pol1", "EM");
  
  func = gr->GetFunction("pol1");
  func->SetLineColor(kRed+1);
  func->SetLineWidth(2);
  gr->GetXaxis()->SetTitle("m_{t,gen}");
  gr->GetYaxis()->SetTitle("#mu");
  
  cObservableCPPar->cd(2);
  gr = new TGraphErrors(3,x,y2,ex,ey2);
  gr->SetTitle("p2");
  gr->Draw("A*");
  gr->Fit("pol1", "EM");
  
  func = gr->GetFunction("pol1");
  func->SetLineColor(kRed+1);
  func->SetLineWidth(2);
  gr->GetXaxis()->SetTitle("m_{t,gen}");
  gr->GetYaxis()->SetTitle("#sigma");
}

TH1F* FindParameters(TString filename, int i)
{

  TFile* file = new TFile(filename);
  analyzeHitFit->cd();

  TF1* voigt = new TF1("voigt", "[0]*TMath::Voigt(x-[1], [2], [3])");
  voigt->SetLineColor(kBlack);
  voigt->SetLineWidth(2);
  voigt->SetParameters(100000, 170, 10, 2);
  
  voigt->SetParLimits(0, 1, 1000000);
  voigt->SetParLimits(1, 150, 200);
  voigt->SetParLimits(2, 5, 20);
  voigt->SetParLimits(3, 2, 2);
  
  TH1F* hSig;
  
  eventTree->Draw("hadTopMass >> hSig(50, 100, 250)", "(bProbSSV*hitFitProb)*(target==1 & (bProbSSV * hitFitProb) > 0.05)");
  
  TH1F *hSig = (TH1F*)gDirectory->Get("hSig");
  
  hSig->GetXaxis()->SetTitle("m_{i}");
  hSig->GetYaxis()->SetTitle("Fraction of entries / 3 GeV");
  hSig->SetFillColor(kRed-8-i);
  
  hSig->Fit("voigt","LEM");

  y0 [i] = voigt->GetParameter(0);
  ey0[i] = voigt->GetParError(0);
 
  y1 [i] = voigt->GetParameter(1);
  ey1[i] = voigt->GetParError(1);
  
  y2 [i] = voigt->GetParameter(2);
  ey2[i] = voigt->GetParError(2);
  
  y3 [i] = voigt->GetParameter(3);
  ey3[i] = voigt->GetParError(3);
  
  hSig->Scale(1/hSig->Integral());
  hSig->Fit("voigt","L");
  
  return hSig;
}
