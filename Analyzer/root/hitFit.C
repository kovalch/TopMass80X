#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "THStack.h"

#include "tdrstyle.C"

void drawcutline(double cutval, double maximum)
{
  TLine *cut = new TLine();
  cut->SetLineWidth(2);
  cut->SetLineStyle(7);
  cut->SetLineColor(1);
  cut->DrawLine(cutval, 0.,cutval, maximum);
}

void setHistStyle(TH1* hist, int kStyle) {
  enum styles             { kHitFit, kKinFit  };
  int color_      [ 2 ] = { kRed+1,  kGreen-3 };
  int fillStyle_  [ 2 ] = { 3654,    3645     };
  
  if (kStyle == kHitFit) hist->Scale(1./2);
  hist->SetLineColor(color_[kStyle]);
  hist->SetLineWidth(2);
  /*
  hist->SetFillColor(color_[kStyle]);
  hist->SetFillStyle(fillStyle_[kStyle]);
  */
  hist->GetYaxis()->SetTitle("Number of permutations");
  hist->GetXaxis()->SetTitle("Permutation type");
}

void hitFit()
{
  enum styles             { kHitFit, kKinFit  };
  int color_      [ 2 ] = { kRed+1,  kGreen-3 };
  int fillStyle_  [ 2 ] = { 3654,    3645     };
  
  setTDRStyle();
  //gStyle->SetPadLeftMargin(0.2);
  //gStyle->SetTitleYOffset(1.75);

  // ---
  //    open input file
  // ---
  TFile* fTTJets = new TFile("/scratch/hh/current/cms/user/mseidel/TTJets1725_abs_test/analyzeTop_1725.root");
  
  // ---
  //    Get trees
  // ---
  TTree* tH = (TTree*) fTTJets->Get("analyzeHitFit/eventTree");
  TTree* tK = (TTree*) fTTJets->Get("analyzeKinFit/eventTree");
  
  //*
  TCanvas* cEfficiency = new TCanvas("cEfficiency", "cEfficiency", 450, 450);
  cEfficiency->cd();
  cEfficiency->SetLogy(1);
  
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("Kinematic Fit Efficiency;Correct permutation efficiency;Bad permutation efficiency");
  
  double aHSigEfficiency[100];
  double aHBkgEfficiency[100];
  double aKSigEfficiency[100];
  double aKBkgEfficiency[100];
  
  double bins = 100;
  
  for (int i = 1; i < bins; i++) {
    TString sHSig("target == 1 & hitFitProb > "); sHSig+= 1./bins * (double)i;
    TString sHBkg("target != 1 & hitFitProb > "); sHBkg+= 1./bins * (double)i;
    TString sKSig("target == 1 & fitProb > ");    sKSig+= 1./bins * (double)i;
    TString sKBkg("target != 1 & fitProb > ");    sKBkg+= 1./bins * (double)i;
    
    std::cout << sHSig << std::endl;
    
    aHSigEfficiency[i-1] = (double)tH->GetEntries(sHSig)/(double)tH->GetEntries("target == 1");
    aHBkgEfficiency[i-1] = (double)tH->GetEntries(sHBkg)/(double)tH->GetEntries("target != 1");
    
    aKSigEfficiency[i-1] = (double)tK->GetEntries(sKSig)/(double)tK->GetEntries("target == 1");
    aKBkgEfficiency[i-1] = (double)tK->GetEntries(sKBkg)/(double)tK->GetEntries("target != 1");
  }
  
  gH = new TGraph(bins-1, aHSigEfficiency, aHBkgEfficiency);
  gH->SetLineColor(color_[kHitFit]);
  gH->SetLineWidth(2);
  mg->Add(gH);
  gK = new TGraph(bins-1, aKSigEfficiency, aKBkgEfficiency);
  gK->SetLineColor(color_[kKinFit]);
  gK->SetLineWidth(2);
  mg->Add(gK);
  
  mg->Draw("AC");
  
  TLegend *legEfficiency = new TLegend(0.25, 0.75, 0.55, 0.9);
  legEfficiency->SetFillStyle(0);
  legEfficiency->SetBorderSize(0);
  legEfficiency->AddEntry( gH, "HitFit"   , "L");
  legEfficiency->AddEntry( gK, "KinFitter", "L");
  legEfficiency->Draw();
  //*/
    
  /*
  // target
  
  TCanvas* cTarget = new TCanvas("cTarget", "cTarget", 450, 450);
  
  cTarget->cd();
  
  tH->Draw("target >> hHTargetB1(30, -12, 3)", "bProbSSV > 0.1");
  tK->Draw("target >> hKTargetB1(30, -12, 3)", "bProbSSV > 0.1");
  
  setHistStyle(hHTargetB1, kHitFit);
  setHistStyle(hKTargetB1, kKinFit);
  
  hKTargetB1->Draw();
  hHTargetB1->Draw("SAME");
  
  TLegend *leg0 = new TLegend(0.7, 0.8, 0.95, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  leg0->AddEntry( hHTargetB1 , "HitFit"   , "F");
  leg0->AddEntry( hKTargetB1 , "KinFitter", "F");
  
  leg0->Draw();
  
  // target, best permutation
  
  TCanvas* cTargetCombi0 = new TCanvas("cTargetCombi0", "cTargetCombi0", 450, 450);
  
  cTargetCombi0->cd();
  
  tH->Draw("target >> hHTargetCombi0B1(30, -12, 3)", "(combi==0 || combi==1) & bProbSSV > 0.1");
  tK->Draw("target >> hKTargetCombi0B1(30, -12, 3)", "combi==0 & bProbSSV > 0.1");
  
  setHistStyle(hHTargetCombi0B1, kHitFit);
  setHistStyle(hKTargetCombi0B1, kKinFit);
  
  hKTargetCombi0B1->Draw();
  hHTargetCombi0B1->Draw("SAME");
  leg0->Draw();
  
  // target, best permutation, cut on fit prob
  
  TCanvas* cTargetCombi0FitProb = new TCanvas("cTargetCombi0FitProb", "cTargetCombi0FitProb", 450, 450);
  
  cTargetCombi0FitProb->cd();
  
  tH->Draw("target >> hHTargetCombi0FitProbB1(30, -12, 3)", "(combi==0 || combi==1) & bProbSSV > 0.1 & hitFitProb > 0.05");
  tK->Draw("target >> hKTargetCombi0FitProbB1(30, -12, 3)", "combi==0 & bProbSSV > 0.1 & fitProb > 0.40");
  
  setHistStyle(hHTargetCombi0FitProbB1, kHitFit);
  setHistStyle(hKTargetCombi0FitProbB1, kKinFit);
  
  hKTargetCombi0FitProbB1->Draw();
  hHTargetCombi0FitProbB1->Draw("SAME");
  leg0->Draw();
  
  // target, weight and cut on fit prob
  
  TCanvas* cTargetFitProbWeighted = new TCanvas("cTargetFitProbWeighted", "cTargetFitProbWeighted", 450, 450);
  
  cTargetFitProbWeighted->cd();
  
  tH->Draw("target >> hHTargetFitProbWeightedB1(30, -12, 3)", "(hitFitProb)*(bProbSSV > 0.1 & hitFitProb > 0.05)");
  tK->Draw("target >> hKTargetFitProbWeightedB1(30, -12, 3)", "(fitProb)*(bProbSSV > 0.1 & fitProb > 0.4)");
  
  setHistStyle(hHTargetFitProbWeightedB1, kHitFit);
  setHistStyle(hKTargetFitProbWeightedB1, kKinFit);
  
  //hKTargetFitProbWeightedB1->GetYaxis()->SetRangeUser(0, 0.6);
  
  hKTargetFitProbWeightedB1->Draw();
  hHTargetFitProbWeightedB1->Draw("SAME");
  leg0->Draw();
  */
  
  // top mass
  
  TCanvas* cHadTopMass = new TCanvas("cHadTopMass", "cHadTopMass", 450, 450);
  
  cHadTopMass->cd();
  
  tH->Draw("hadTopMass >> hHHadTopMass(100, 100, 250)", "target==1");
  tK->Draw("hadTopMass >> hKHadTopMass(100, 100, 250)", "target==1");
  
  setHistStyle(hHHadTopMass, kHitFit);
  setHistStyle(hKHadTopMass, kKinFit);
  
  hHHadTopMass->GetXaxis()->SetTitle("m_{i}");
  hKHadTopMass->GetXaxis()->SetTitle("m_{i}");
  
  //hKHadTopMass->GetYaxis()->SetRangeUser(0, 0.085);
  
  hKHadTopMass->Draw();
  hHHadTopMass->Draw("SAME");
  
  TLegend *leg0 = new TLegend(0.7, 0.8, 0.95, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  leg0->AddEntry( hHHadTopMass , "HitFit"   , "F");
  leg0->AddEntry( hKHadTopMass , "KinFitter", "F");
  leg0->Draw();
  
  
  /*
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
