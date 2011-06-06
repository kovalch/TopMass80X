#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "THStack.h"

void drawcutline(double cutval, double maximum)
{
  TLine *cut = new TLine();
  cut->SetLineWidth(2);
  cut->SetLineStyle(7);
  cut->SetLineColor(1);
  cut->DrawLine(cutval, 0.,cutval, maximum);
}

void hitFit()
{  
  gStyle->SetOptStat(0);
  
  enum styles             {kSig,   kBkg,   kWjets,   kData,  };
  int color_      [ 4 ] = {kRed+1, kRed-7, kGreen-3, kBlack, };
  int markerStyle_[ 4 ] = {20,     22,     23,       20,     };

  
  // ---
  //    open input file
  // ---
  TFile* fTTJets = new TFile("./analyzeTop_1725.root");
  
  // ---
  //    Get trees
  // ---
  TTree* tH = (TTree*) fTTJets->Get("analyzeHitFit/eventTree");
  TTree* tK = (TTree*) fTTJets->Get("analyzeKinFit/eventTree");
  
  // target
  
  TCanvas* cTarget = new TCanvas("cTarget", "cTarget", 1000, 400);
  cTarget->Divide(3,1);
  
  cTarget->cd(1);
  
  tH->Draw("target >> hHTarget(30, -12, 3)");
  tK->Draw("target >> hKTarget(30, -12, 3)");
  
  hHTarget->Scale(1./2);
  
  hHTarget->SetLineColor(color_[kSig]);
  hKTarget->SetLineColor(color_[kWjets]);
  
  hHTarget->Draw();
  hKTarget->Draw("SAME");
  
  cTarget->cd(2);
  
  tH->Draw("target >> hHTargetB1(30, -12, 3)", "bProbSSV > 0.1");
  tK->Draw("target >> hKTargetB1(30, -12, 3)", "bProbSSV > 0.1");
  
  hHTargetB1->Scale(1./2);
  
  hHTargetB1->SetLineColor(color_[kSig]);
  hKTargetB1->SetLineColor(color_[kWjets]);
  
  hHTargetB1->Draw();
  hKTargetB1->Draw("SAME");
  
  cTarget->cd(3);
  
  tH->Draw("target >> hHTargetB2(30, -12, 3)", "bProbSSV > 0.3");
  tK->Draw("target >> hKTargetB2(30, -12, 3)", "bProbSSV > 0.3");
  
  hHTargetB2->Scale(1./2);
  
  hHTargetB2->SetLineColor(color_[kSig]);
  hKTargetB2->SetLineColor(color_[kWjets]);
  
  hHTargetB2->Draw();
  hKTargetB2->Draw("SAME");
  
  // target, best permutation
  
  TCanvas* cTargetCombi0 = new TCanvas("cTargetCombi0", "cTargetCombi0", 1000, 400);
  cTargetCombi0->Divide(3,1);
  
  cTargetCombi0->cd(1);
  
  tH->Draw("target >> hHTargetCombi0(30, -12, 3)", "combi==0");
  tK->Draw("target >> hKTargetCombi0(30, -12, 3)", "combi==0");
  
  hHTargetCombi0->SetLineColor(color_[kSig]);
  hKTargetCombi0->SetLineColor(color_[kWjets]);
  
  hHTargetCombi0->Draw();
  hKTargetCombi0->Draw("SAME");
  
  cTargetCombi0->cd(2);
  
  tH->Draw("target >> hHTargetCombi0B1(30, -12, 3)", "combi==0 & bProbSSV > 0.1");
  tK->Draw("target >> hKTargetCombi0B1(30, -12, 3)", "combi==0 & bProbSSV > 0.1");
  
  hHTargetCombi0B1->SetLineColor(color_[kSig]);
  hKTargetCombi0B1->SetLineColor(color_[kWjets]);
  
  hHTargetCombi0B1->Draw();
  hKTargetCombi0B1->Draw("SAME");
  
  cTargetCombi0->cd(3);
  
  tH->Draw("target >> hHTargetCombi0B2(30, -12, 3)", "combi==0 & bProbSSV > 0.3");
  tK->Draw("target >> hKTargetCombi0B2(30, -12, 3)", "combi==0 & bProbSSV > 0.3");
  
  hHTargetCombi0B2->SetLineColor(color_[kSig]);
  hKTargetCombi0B2->SetLineColor(color_[kWjets]);
  
  hHTargetCombi0B2->Draw();
  hKTargetCombi0B2->Draw("SAME");
  
  // target, best permutation, cut on fit prob
  
  TCanvas* cTargetCombi0FitProb = new TCanvas("cTargetCombi0FitProb", "cTargetCombi0FitProb", 1000, 400);
  cTargetCombi0FitProb->Divide(3,1);
  
  cTargetCombi0FitProb->cd(1);
  
  tH->Draw("target >> hHTargetCombi0FitProb(30, -12, 3)", "combi==0 & hitFitProb > 0.05");
  tK->Draw("target >> hKTargetCombi0FitProb(30, -12, 3)", "combi==0 & fitProb > 0.05");
  
  hHTargetCombi0FitProb->SetLineColor(color_[kSig]);
  hKTargetCombi0FitProb->SetLineColor(color_[kWjets]);
  
  hHTargetCombi0FitProb->GetYaxis()->SetRangeUser(0, 20000);
  
  hHTargetCombi0FitProb->Draw();
  hKTargetCombi0FitProb->Draw("SAME");
  
  cTargetCombi0FitProb->cd(2);
  
  tH->Draw("target >> hHTargetCombi0FitProbB1(30, -12, 3)", "combi==0 & bProbSSV > 0.1 & hitFitProb > 0.05");
  tK->Draw("target >> hKTargetCombi0FitProbB1(30, -12, 3)", "combi==0 & bProbSSV > 0.1 & fitProb > 0.05");
  
  hHTargetCombi0FitProbB1->SetLineColor(color_[kSig]);
  hKTargetCombi0FitProbB1->SetLineColor(color_[kWjets]);
  
  hHTargetCombi0FitProbB1->GetYaxis()->SetRangeUser(0, 8000);
  
  hHTargetCombi0FitProbB1->Draw();
  hKTargetCombi0FitProbB1->Draw("SAME");
  
  cTargetCombi0FitProb->cd(3);
  
  tH->Draw("target >> hHTargetCombi0FitProbB2(30, -12, 3)", "combi==0 & bProbSSV > 0.3 & hitFitProb > 0.05");
  tK->Draw("target >> hKTargetCombi0FitProbB2(30, -12, 3)", "combi==0 & bProbSSV > 0.3 & fitProb > 0.05");
  
  hHTargetCombi0FitProbB2->SetLineColor(color_[kSig]);
  hKTargetCombi0FitProbB2->SetLineColor(color_[kWjets]);
  
  hHTargetCombi0FitProbB2->Draw();
  hKTargetCombi0FitProbB2->Draw("SAME");
  
  
  /*
  // fitProb
  
  TCanvas* cFitProb = new TCanvas("cFitProb", "cFitProb", 1000, 400);
  cFitProb->Divide(4,1);
  
  cFitProb->cd(1);
  
  tH->Draw("hitFitProb >> hHFitProbCombi0(30, 0, 1)", "combi==0");
  tK->Draw("fitProb >> hKFitProbCombi0(30, 0, 1)", "combi==0");
  
  hHFitProbCombi0->SetLineColor(color_[kSig]);
  hKFitProbCombi0->SetLineColor(color_[kWjets]);
  
  hHFitProbCombi0->Draw();
  hKFitProbCombi0->Draw("SAME");
  
  cFitProb->cd(2);
  
  tH->Draw("hitFitProb >> hHFitProbCombi0B1(30, 0, 1)", "combi==0 & bProbSSV > 0.1");
  tK->Draw("fitProb >> hKFitProbCombi0B1(30, 0, 1)", "combi==0 & bProbSSV > 0.1");
  
  hHFitProbCombi0B1->SetLineColor(color_[kSig]);
  hKFitProbCombi0B1->SetLineColor(color_[kWjets]);
  
  hHFitProbCombi0B1->Draw();
  hKFitProbCombi0B1->Draw("SAME");
  
  cFitProb->cd(3);
  
  tH->Draw("hitFitProb >> hHFitProbCombi0B2(30, 0, 1)", "combi==0 & bProbSSV > 0.3");
  tK->Draw("fitProb >> hKFitProbCombi0B2(30, 0, 1)", "combi==0 & bProbSSV > 0.3");
  
  hHFitProbCombi0B2->SetLineColor(color_[kSig]);
  hKFitProbCombi0B2->SetLineColor(color_[kWjets]);
  
  hHFitProbCombi0B2->Draw();
  hKFitProbCombi0B2->Draw("SAME");
  
  cFitProb->cd(4);
  
  tH->Draw("hitFitProb >> hHFitProbTarget1(30, 0, 1)", "target==1");
  tK->Draw("fitProb >> hKFitProbTarget1(30, 0, 1)", "target==1");
  
  hHFitProbTarget1->Scale(1./2);
  
  hHFitProbTarget1->SetLineColor(color_[kSig]);
  hKFitProbTarget1->SetLineColor(color_[kWjets]);
  
  hHFitProbTarget1->Draw();
  hKFitProbTarget1->Draw("SAME");
  */
  
  /*
  // top pt
  
  TCanvas* cHadTopPt = new TCanvas("cHadTopPt", "cHadTopPt", 1000, 400);
  cHadTopPt->Divide(4,1);
  
  cHadTopPt->cd(1);
  
  tH->Draw("hadTopPt >> hHHadTopPtCombi0(30, 0, 300)", "combi==0 & hitFitProb > 0.05");
  tK->Draw("hadTopPt >> hKHadTopPtCombi0(30, 0, 300)", "combi==0 & fitProb > 0.05");
  
  hHHadTopPtCombi0->Scale(1./hHHadTopPtCombi0->Integral());
  hKHadTopPtCombi0->Scale(1./hKHadTopPtCombi0->Integral());
  
  hHHadTopPtCombi0->SetLineColor(color_[kSig]);
  hKHadTopPtCombi0->SetLineColor(color_[kWjets]);
  
  hHHadTopPtCombi0->Draw();
  hKHadTopPtCombi0->Draw("SAME");
  
  cHadTopPt->cd(2);
  
  tH->Draw("hadTopPt >> hHHadTopPtCombi0B1(30, 0, 300)", "combi==0 & bProbSSV > 0.1 & hitFitProb > 0.05");
  tK->Draw("hadTopPt >> hKHadTopPtCombi0B1(30, 0, 300)", "combi==0 & bProbSSV > 0.1 & fitProb > 0.05");
  
  hHHadTopPtCombi0B1->Scale(1./hHHadTopPtCombi0B1->Integral());
  hKHadTopPtCombi0B1->Scale(1./hKHadTopPtCombi0B1->Integral());
  
  hHHadTopPtCombi0B1->SetLineColor(color_[kSig]);
  hKHadTopPtCombi0B1->SetLineColor(color_[kWjets]);
  
  hHHadTopPtCombi0B1->Draw();
  hKHadTopPtCombi0B1->Draw("SAME");
  
  cHadTopPt->cd(3);
  
  tH->Draw("hadTopPt >> hHHadTopPtCombi0B2(30, 0, 300)", "combi==0 & bProbSSV > 0.3 & hitFitProb > 0.05");
  tK->Draw("hadTopPt >> hKHadTopPtCombi0B2(30, 0, 300)", "combi==0 & bProbSSV > 0.3 & fitProb > 0.05");
  
  hHHadTopPtCombi0B2->Scale(1./hHHadTopPtCombi0B2->Integral());
  hKHadTopPtCombi0B2->Scale(1./hKHadTopPtCombi0B2->Integral());
  
  hHHadTopPtCombi0B2->SetLineColor(color_[kSig]);
  hKHadTopPtCombi0B2->SetLineColor(color_[kWjets]);
  
  hHHadTopPtCombi0B2->Draw();
  hKHadTopPtCombi0B2->Draw("SAME");
  
  cHadTopPt->cd(4);
  
  tH->Draw("hadTopPt >> hHHadTopPtTarget1(30, 0, 300)", "target==1");
  tK->Draw("hadTopPt >> hKHadTopPtTarget1(30, 0, 300)", "target==1");
  
  hHHadTopPtTarget1->Scale(1./hHHadTopPtTarget1->Integral());
  hKHadTopPtTarget1->Scale(1./hKHadTopPtTarget1->Integral());
  
  hHHadTopPtTarget1->SetLineColor(color_[kSig]);
  hKHadTopPtTarget1->SetLineColor(color_[kWjets]);
  
  hHHadTopPtTarget1->Draw();
  hKHadTopPtTarget1->Draw("SAME");
  
  // top eta
  
  TCanvas* cHadTopPt = new TCanvas("cHadTopEta", "cHadTopEta", 1000, 400);
  cHadTopEta->Divide(4,1);
  
  cHadTopEta->cd(1);
  
  tH->Draw("hadTopEta >> hHHadTopEtaCombi0(30, -3, 3)", "combi==0 & hitFitProb > 0.05");
  tK->Draw("hadTopEta >> hKHadTopEtaCombi0(30, -3, 3)", "combi==0 & fitProb > 0.05");
  
  hHHadTopEtaCombi0->SetLineColor(color_[kSig]);
  hKHadTopEtaCombi0->SetLineColor(color_[kWjets]);
  
  hHHadTopEtaCombi0->Draw();
  hKHadTopEtaCombi0->Draw("SAME");
  
  cHadTopEta->cd(2);
  
  tH->Draw("hadTopEta >> hHHadTopEtaCombi0B1(30, -3, 3)", "combi==0 & bProbSSV > 0.1 & hitFitProb > 0.05");
  tK->Draw("hadTopEta >> hKHadTopEtaCombi0B1(30, -3, 3)", "combi==0 & bProbSSV > 0.1 & fitProb > 0.05");
  
  hHHadTopEtaCombi0B1->SetLineColor(color_[kSig]);
  hKHadTopEtaCombi0B1->SetLineColor(color_[kWjets]);
  
  hHHadTopEtaCombi0B1->Draw();
  hKHadTopEtaCombi0B1->Draw("SAME");
  
  cHadTopEta->cd(3);
  
  tH->Draw("hadTopEta >> hHHadTopEtaCombi0B2(30, -3, 3)", "combi==0 & bProbSSV > 0.3 & hitFitProb > 0.05");
  tK->Draw("hadTopEta >> hKHadTopEtaCombi0B2(30, -3, 3)", "combi==0 & bProbSSV > 0.3 & fitProb > 0.05");
  
  hHHadTopEtaCombi0B2->SetLineColor(color_[kSig]);
  hKHadTopEtaCombi0B2->SetLineColor(color_[kWjets]);
  
  hHHadTopEtaCombi0B2->Draw();
  hKHadTopEtaCombi0B2->Draw("SAME");
  
  cHadTopEta->cd(4);
  
  tH->Draw("hadTopEta >> hHHadTopEtaTarget1(30, -3, 3)", "target==1");
  tK->Draw("hadTopEta >> hKHadTopEtaTarget1(30, -3, 3)", "target==1");
  
  hHHadTopEtaTarget1->Scale(1./2);
  
  hHHadTopEtaTarget1->SetLineColor(color_[kSig]);
  hKHadTopEtaTarget1->SetLineColor(color_[kWjets]);
  
  hHHadTopEtaTarget1->Draw();
  hKHadTopEtaTarget1->Draw("SAME");
  */
  
  // top mass
  
  TCanvas* cHadTopMass = new TCanvas("cHadTopMass", "cHadTopMass", 1000, 400);
  cHadTopMass->Divide(4,1);
  
  cHadTopMass->cd(1);
  
  tH->Draw("hadTopMass >> hHHadTopMassCombi0(50, 100, 250)", "combi==0 & hitFitProb > 0.05");
  tK->Draw("hadTopMass >> hKHadTopMassCombi0(50, 100, 250)", "combi==0 & fitProb > 0.05");
  
  hHHadTopMassCombi0->Scale(1./hHHadTopMassCombi0->Integral());
  hKHadTopMassCombi0->Scale(1./hKHadTopMassCombi0->Integral());
  
  hHHadTopMassCombi0->SetLineColor(color_[kSig]);
  hKHadTopMassCombi0->SetLineColor(color_[kWjets]);
  
  hHHadTopMassCombi0->Draw();
  hKHadTopMassCombi0->Draw("SAME");
  
  cHadTopMass->cd(2);
  
  tH->Draw("hadTopMass >> hHHadTopMassCombi0B1(50, 100, 250)", "combi==0 & bProbSSV > 0.1 & hitFitProb > 0.05");
  tK->Draw("hadTopMass >> hKHadTopMassCombi0B1(50, 100, 250)", "combi==0 & bProbSSV > 0.1 & fitProb > 0.05");
  
  hHHadTopMassCombi0B1->Scale(1./hHHadTopMassCombi0B1->Integral());
  hKHadTopMassCombi0B1->Scale(1./hKHadTopMassCombi0B1->Integral());
  
  hHHadTopMassCombi0B1->SetLineColor(color_[kSig]);
  hKHadTopMassCombi0B1->SetLineColor(color_[kWjets]);
  
  hHHadTopMassCombi0B1->Draw();
  hKHadTopMassCombi0B1->Draw("SAME");
  
  cHadTopMass->cd(3);
  
  tH->Draw("hadTopMass >> hHHadTopMassCombi0B2(50, 100, 250)", "combi==0 & bProbSSV > 0.3 & hitFitProb > 0.05");
  tK->Draw("hadTopMass >> hKHadTopMassCombi0B2(50, 100, 250)", "combi==0 & bProbSSV > 0.3 & fitProb > 0.05");
  
  hHHadTopMassCombi0B2->Scale(1./hHHadTopMassCombi0B2->Integral());
  hKHadTopMassCombi0B2->Scale(1./hKHadTopMassCombi0B2->Integral());
  
  hHHadTopMassCombi0B2->SetLineColor(color_[kSig]);
  hKHadTopMassCombi0B2->SetLineColor(color_[kWjets]);
  
  hHHadTopMassCombi0B2->Draw();
  hKHadTopMassCombi0B2->Draw("SAME");
  
  cHadTopMass->cd(4);
  
  tH->Draw("hadTopMass >> hHHadTopMassTarget1(50, 100, 250)", "target==1");
  tK->Draw("hadTopMass >> hKHadTopMassTarget1(50, 100, 250)", "target==1");
  
  hHHadTopMassTarget1->Scale(1./2);
  
  hHHadTopMassTarget1->SetLineColor(color_[kSig]);
  hKHadTopMassTarget1->SetLineColor(color_[kWjets]);
  
  hHHadTopMassTarget1->Draw();
  hKHadTopMassTarget1->Draw("SAME");

  
  /*
  tT->Draw("log10(bProbSSV*hitFitProb) >> hTSig(40, -10, 0)", "target==1");
  tT->Draw("log10(bProbSSV*hitFitProb) >> hTBkg(40, -10, 0)", "target!=1");
  tW->Draw("log10(bProbSSV*hitFitProb) >> hW   (40, -10, 0)");
  tD->Draw("log10(bProbSSV*hitFitProb) >> hD   (40, -10, 0)");
  
  TH1F* hNull = new TH1F("null", "", 40, -10, 0);
  hNull->GetYaxis()->SetRangeUser(0.1, 1000);
  hNull->GetXaxis()->SetTitle("log(w_{i})");
  
  TH1F* hTSig = (TH1F*) gDirectory->Get("hTSig");
  TH1F* hTBkg = (TH1F*) gDirectory->Get("hTBkg");
  TH1F* hW    = (TH1F*) gDirectory->Get("hW");
  TH1F* hD    = (TH1F*) gDirectory->Get("hD");
  
  hTBkg->Scale(1/normalization * 188/35.9*luminosity);
  hTSig->Scale(1/normalization * 188/35.9*luminosity);
  hW   ->Scale(1/hW->Integral() *  21/35.9*luminosity);
  
  hTSig->SetFillColor(color_[kSig]);
  hTBkg->SetFillColor(color_[kBkg]);
  hW   ->SetFillColor(color_[kWjets]);
  hD   ->SetMarkerStyle(markerStyle_[kData]);
  
  THStack* stack = new THStack("stack", "");
  stack->Add(hW);
  stack->Add(hTSig);
  stack->Add(hTBkg);
  
  // ---
  //    create legend
  // ---
  // samples: separate canvas
  TLegend *leg0 = new TLegend(0.6, 0.7, 0.85, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  leg0->AddEntry( hD , "Data (23 pb ^{-1})"       , "PL");
  leg0->AddEntry( hW , "W#rightarrowl#nu"              , "F" );
  leg0->AddEntry( hTBkg , "t#bar{t} wp"                , "F" );
  leg0->AddEntry( hTSig , "t#bar{t} cp"  , "F" );

  hNull->Draw();
  stack->Draw("SAME");
  hD   ->Draw("E,SAME");
  leg0->Draw("");
  
  drawcutline(-1.3, 100);
  
  */
}
