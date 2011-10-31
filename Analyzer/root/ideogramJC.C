#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"

void ideogramJC()
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TFile* file = new TFile("/scratch/hh/current/cms/user/mseidel/TTJets1725-S4_1.00_1.00/analyzeTop.root");
  analyzeHitFit->cd();
  
  /*
  TF1 *gaus2 = new TF1("gaus2", "[0]*TMath::Gaus(x, [1], [2], 1)+[3]*TMath::Gaus(x, [4], [5], 1)");
  gaus2->SetParameters(100000,170,20,1000,250,50);
  
  gaus2->SetParLimits(0, 20, 1000000);
  gaus2->SetParLimits(1, 150, 200);
  gaus2->SetParLimits(2, 10, 30);
  gaus2->SetParLimits(3, 20, 1000000);
  gaus2->SetParLimits(4, 170, 300);
  gaus2->SetParLimits(5, 30, 100);
  */
  
  //TH1F* hhadQJC;
  
  int bins = 60;
  
  eventTree->Draw("abs(leptonC+lepBJC) >> h1(60, 0, 3)", "target==1");
  eventTree->Draw("abs(leptonC+lepBJC) >> hS(60, 0, 3)");
  
  TH1F *h1 = (TH1F*)gDirectory->Get("h1");
  TH1F *hS = (TH1F*)gDirectory->Get("hS");
   
  TH1F *hQProb = new TH1F("hQProb", "hQProb", bins, 0, 3);
  hQProb->Divide(h1, hS);
  
  hQProb->Draw();

  TF1 *exponential = new TF1("exponential", "exp([0]+[1]*x)", 1.25, 5);
  exponential->SetLineColor(kRed+1);
  
  hQProb->Fit("exponential", "WREM");
  
  hQProb->GetXaxis()->SetTitle("B_{SSVHE}");
  hQProb->GetYaxis()->SetTitle("P_{lq}(B)");
  
  //hCombBkg->Scale(1/hCombBkg->Integral());
  /*
  hCombBkg->Fit("gaus2","LEM");

  y0 [i] = gaus2->GetParameter(0)/gaus2->GetParameter(3);;
  //ey0[i] = gaus2->GetParError(0);
  ey0[i] = sqrt(pow(1/gaus2->GetParameter(3)*gaus2->GetParError(0), 2)
            + pow(gaus2->GetParameter(0)/pow(gaus2->GetParameter(3), 2)*gaus2->GetParError(3), 2));
  */


}
