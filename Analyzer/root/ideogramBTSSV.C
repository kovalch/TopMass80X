#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"

void ideogramBTSSV()
{

  TFile* file = new TFile("analyzeTop_1725.root");
  analyzeGenMatch->cd();
  
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
  
  //TH1F* hHadQBSSV;
  
  int bins = 60;
  
  eventTree->Draw("hadQBSSV >> hHadQBSSV(60, -1, 5)");
  eventTree->Draw("hadQbarBSSV >> hHadQbarBSSV(60, -1, 5)");
  eventTree->Draw("hadBBSSV >> hHadBBSSV(60, -1, 5)");
  eventTree->Draw("lepBBSSV >> hLepBBSSV(60, -1, 5)");
  
  TH1F *hHadQBSSV = (TH1F*)gDirectory->Get("hHadQBSSV");
  TH1F *hHadQbarBSSV = (TH1F*)gDirectory->Get("hHadQbarBSSV");
  TH1F *hHadBBSSV = (TH1F*)gDirectory->Get("hHadBBSSV");
  TH1F *hLepBBSSV = (TH1F*)gDirectory->Get("hLepBBSSV");
  
  TH1F *hQSum = new TH1F("hQSum", "hQSum", bins, -1, 5);
  hQSum->Add(hHadQBSSV, hHadQbarBSSV);
  
  TH1F *hBSum = new TH1F("hBSum", "hBSum", bins, -1, 5);
  hBSum->Add(hHadBBSSV, hLepBBSSV);
  
  TH1F *hQBSum = new TH1F("hQBSum", "hQBSum", bins, -1, 5);
  hQBSum->Add(hQSum, hBSum);
  
  TH1F *hQProb = new TH1F("hQProb", "hQProb", bins, -1, 5);
  hQProb->Divide(hQSum, hQBSum);
  
  hQProb->Draw();
  
  TF1 *voigt = new TF1("voigt", "[0]*TMath::Voigt(x, [1], [2])", 0, 50);
  voigt->SetParameters(7, 2, 3);

  TF1 *exponential = new TF1("exponential", "exp([0]+[1]*x)", 1.25, 5);
  
  hQProb->Fit("exponential", "WREM");
  
  //hCombBkg->Scale(1/hCombBkg->Integral());
  /*
  hCombBkg->Fit("gaus2","LEM");

  y0 [i] = gaus2->GetParameter(0)/gaus2->GetParameter(3);;
  //ey0[i] = gaus2->GetParError(0);
  ey0[i] = sqrt(pow(1/gaus2->GetParameter(3)*gaus2->GetParError(0), 2)
            + pow(gaus2->GetParameter(0)/pow(gaus2->GetParameter(3), 2)*gaus2->GetParError(3), 2));
  */


}
