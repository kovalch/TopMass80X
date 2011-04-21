#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"

void ideogramBT()
{

  TFile* file = new TFile("analyzeTop_1725.root");
  analyzeKinFit->cd();
  
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
  
  //TH1F* hHadQBTCHE;
  
  int bins = 150;
  
  eventTree->Draw("hadQBTCHE >> hHadQBTCHE(150, -100, 50)","target==1");
  eventTree->Draw("hadQbarBTCHE >> hHadQbarBTCHE(150, -100, 50)","target==1");
  eventTree->Draw("hadBBTCHE >> hHadBBTCHE(150, -100, 50)","target==1");
  eventTree->Draw("lepBBTCHE >> hLepBBTCHE(150, -100, 50)","target==1");
  
  TH1F *hHadQBTCHE = (TH1F*)gDirectory->Get("hHadQBTCHE");
  TH1F *hHadQbarBTCHE = (TH1F*)gDirectory->Get("hHadQbarBTCHE");
  TH1F *hHadBBTCHE = (TH1F*)gDirectory->Get("hHadBBTCHE");
  TH1F *hLepBBTCHE = (TH1F*)gDirectory->Get("hLepBBTCHE");
  
  TH1F *hQSum = new TH1F("hQSum", "hQSum", bins, -100, 50);
  hQSum->Add(hHadQBTCHE, hHadQbarBTCHE);
  
  TH1F *hBSum = new TH1F("hBSum", "hBSum", bins, -100, 50);
  hBSum->Add(hHadBBTCHE, hLepBBTCHE);
  
  TH1F *hQBSum = new TH1F("hQBSum", "hQBSum", bins, -100, 50);
  hQBSum->Add(hQSum, hBSum);
  
  TH1F *hQProb = new TH1F("hQProb", "hQProb", bins, -100, 50);
  hQProb->Divide(hQSum, hQBSum);
  
  hQProb->Draw();
  
  TF1 *voigt = new TF1("voigt", "[0]*TMath::Voigt(x, [1], [2])", 0, 50);
  voigt->SetParameters(7, 2, 3);
  
  hQProb->Fit("voigt", "REM");
  
  //hCombBkg->Scale(1/hCombBkg->Integral());
  /*
  hCombBkg->Fit("gaus2","LEM");

  y0 [i] = gaus2->GetParameter(0)/gaus2->GetParameter(3);;
  //ey0[i] = gaus2->GetParError(0);
  ey0[i] = sqrt(pow(1/gaus2->GetParameter(3)*gaus2->GetParError(0), 2)
            + pow(gaus2->GetParameter(0)/pow(gaus2->GetParameter(3), 2)*gaus2->GetParError(3), 2));
  */


}
