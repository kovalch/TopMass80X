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

enum styles             { kHitFit, kKinFit, kGenMatch, kMVADisc};
int color_      [ 4 ] = { kRed,    kGreen,  kBlack,    kGray};

void drawcutline(double cutval, double maximum)
{
  TLine *cut = new TLine();
  cut->SetLineWidth(2);
  cut->SetLineStyle(7);
  cut->SetLineColor(1);
  cut->DrawLine(cutval, 0.,cutval, maximum);
}

void setHistStyle(TString titleX, TH1* hist, int kStyle, int col) {
  hist->Scale(1./hist->Integral());
  
  hist->SetLineColor(color_[kStyle]+col);
  hist->SetLineWidth(2);
  
  //hist->GetYaxis()->SetTitle("Fraction of permutations");
  hist->GetXaxis()->SetTitle(titleX);
}

void hitFitKinFitter()
{
  setTDRStyle();
  //gStyle->SetPadLeftMargin(0.2);
  //gStyle->SetTitleYOffset(1.75);

  // ---
  //    open input file
  // ---
  TFile* fTTJets = new TFile("/scratch/hh/current/cms/user/mseidel/Spring11_TTJets1725_1.00_1.00_FitComparison/analyzeTop.root");
  
  // ---
  //    Get trees
  // ---
  TTree* tH = (TTree*) fTTJets->Get("analyzeHitFit/eventTree");
  TTree* tK = (TTree*) fTTJets->Get("analyzeKinFit/eventTree");
  TTree* tG = (TTree*) fTTJets->Get("analyzeGenMatch/eventTree");
  TTree* tM = (TTree*) fTTJets->Get("analyzeMVADisc/eventTree");
  
  //* efficiency
  
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
    TString sHSig("target == 1 & bProbSSV > 0.1 & hitFitProb > "); sHSig+= 1./bins * (double)i;
    TString sHBkg("target != 1 & bProbSSV > 0.1 & hitFitProb > "); sHBkg+= 1./bins * (double)i;
    TString sKSig("target == 1 & bProbSSV > 0.1 & fitProb > ");    sKSig+= 1./bins * (double)i;
    TString sKBkg("target != 1 & bProbSSV > 0.1 & fitProb > ");    sKBkg+= 1./bins * (double)i;
    
    std::cout << sHSig << std::endl;
    
    aHSigEfficiency[i-1] = (double)tH->GetEntries(sHSig)/(double)tH->GetEntries("target == 1 & bProbSSV > 0.1");
    aHBkgEfficiency[i-1] = (double)tH->GetEntries(sHBkg)/(double)tH->GetEntries("target != 1 & bProbSSV > 0.1");
    
    aKSigEfficiency[i-1] = (double)tK->GetEntries(sKSig)/(double)tK->GetEntries("target == 1 & bProbSSV > 0.1");
    aKBkgEfficiency[i-1] = (double)tK->GetEntries(sKBkg)/(double)tK->GetEntries("target != 1 & bProbSSV > 0.1");
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
  
  //* top mass resolution
  
  TCanvas* cHadTopMassRes = new TCanvas("cHadTopMassRes", "cHadTopMassRes", 450, 450);
  
  cHadTopMassRes->cd();
  
  tH->Draw("hadTopMass-genHadTopMass >> hHHadTopMassRes(25, -100, 100)", "target==1");
  tK->Draw("hadTopMass-genHadTopMass >> hKHadTopMassRes(25, -100, 100)", "target==1");
  tG->Draw("hadTopMass-genHadTopMass >> hGHadTopMassRes(25, -100, 100)", "");
  
  setHistStyle("m_{t,rec}-m_{t,gen}", hHHadTopMassRes, kHitFit, 0);
  setHistStyle("m_{t,rec}-m_{t,gen}", hKHadTopMassRes, kKinFit, 0);
  setHistStyle("m_{t,rec}-m_{t,gen}", hGHadTopMassRes, kGenMatch, 0);
  
  //hKHadTopMass->Draw();
  hHHadTopMassRes->Draw("C");
  hKHadTopMassRes->Draw("C,SAME");
  hGHadTopMassRes->Draw("C,SAME");
  
  TLegend *leg0 = new TLegend(0.65, 0.55, 0.95, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  leg0->AddEntry( hHHadTopMassRes , "HitFit", "F");
  leg0->AddEntry( hKHadTopMassRes , "KinFitter", "F");
  leg0->AddEntry( hGHadTopMassRes , "GenMatch", "F");
  leg0->Draw();
  //*/
  
  //* top mass bias
  
  TCanvas* cHadTopMass = new TCanvas("cHadTopMass", "cHadTopMass", 450, 450);
  
  cHadTopMass->cd();
  
  tH->Draw("hadTopMass-genHadTopMass >> hHTHadTopMass(25, -100, 100)", "target==1 & bProbSSV>0.1");
  tH->Draw("hadTopMass-genHadTopMass >> hHPHadTopMass(25, -100, 100)", "hitFitProb>0.3 & bProbSSV>0.1");
  tH->Draw("hadTopMass-genHadTopMass >> hHBHadTopMass(25, -100, 100)", "combi==0 & bProbSSV>0.1");
  tH->Draw("hadTopMass-genHadTopMass >> hHAHadTopMass(25, -100, 100)", "bProbSSV>0.1");
  
  tK->Draw("hadTopMass-genHadTopMass >> hKTHadTopMass(25, -100, 100)", "target==1 & bProbSSV>0.1");
  tK->Draw("hadTopMass-genHadTopMass >> hKPHadTopMass(25, -100, 100)", "fitProb>0.4 & bProbSSV>0.1");
  tK->Draw("hadTopMass-genHadTopMass >> hKBHadTopMass(25, -100, 100)", "combi==0 & bProbSSV>0.1");
  tK->Draw("hadTopMass-genHadTopMass >> hKAHadTopMass(25, -100, 100)", "bProbSSV>0.1");
  
  setHistStyle("m_{t,rec}-m_{t,gen}", hHTHadTopMass, kHitFit, 0);
  setHistStyle("m_{t,rec}-m_{t,gen}", hHPHadTopMass, kHitFit, 1);
  setHistStyle("m_{t,rec}-m_{t,gen}", hHBHadTopMass, kHitFit, 2);
  setHistStyle("m_{t,rec}-m_{t,gen}", hHAHadTopMass, kHitFit, 3);
  
  setHistStyle("m_{t,rec}-m_{t,gen}", hKTHadTopMass, kKinFit, 0);
  setHistStyle("m_{t,rec}-m_{t,gen}", hKPHadTopMass, kKinFit, 1);
  setHistStyle("m_{t,rec}-m_{t,gen}", hKBHadTopMass, kKinFit, 2);
  setHistStyle("m_{t,rec}-m_{t,gen}", hKAHadTopMass, kKinFit, 3);
  
  hHTHadTopMass->Draw("C");
  hHPHadTopMass->Draw("C,SAME");
  hHBHadTopMass->Draw("C,SAME");
  hHAHadTopMass->Draw("C,SAME");
  
  hKTHadTopMass->Draw("C,SAME");
  hKPHadTopMass->Draw("C,SAME");
  hKBHadTopMass->Draw("C,SAME");
  hKAHadTopMass->Draw("C,SAME");
  
  TLegend *leg0 = new TLegend(0.65, 0.55, 0.95, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  
  leg0->AddEntry( hHTHadTopMass , "HitFit correct", "F");
  leg0->AddEntry( hHPHadTopMass , "HitFit cut50", "F");
  leg0->AddEntry( hHBHadTopMass , "HitFit 'best'", "F");
  leg0->AddEntry( hHAHadTopMass , "HitFit all", "F");
  
  leg0->AddEntry( hKTHadTopMass , "KinFit correct", "F");
  leg0->AddEntry( hKPHadTopMass , "KinFit cut50", "F");
  leg0->AddEntry( hKBHadTopMass , "KinFit 'best'", "F");
  leg0->AddEntry( hKAHadTopMass , "KinFit all", "F");
  
  leg0->Draw();
  //*/
  
  //* top pt bias
  
  TCanvas* cHadTopPt = new TCanvas("cHadTopPt", "cHadTopPt", 450, 450);
  
  cHadTopPt->cd();
  
  tH->Draw("hadTopPt-genHadTopPt >> hHTHadTopPt(25, -100, 100)", "target==1 & bProbSSV>0.1");
  tH->Draw("hadTopPt-genHadTopPt >> hHPHadTopPt(25, -100, 100)", "hitFitProb>0.3 & bProbSSV>0.1");
  tH->Draw("hadTopPt-genHadTopPt >> hHBHadTopPt(25, -100, 100)", "combi==0 & bProbSSV>0.1");
  tH->Draw("hadTopPt-genHadTopPt >> hHAHadTopPt(25, -100, 100)", "bProbSSV>0.1");
  
  tK->Draw("hadTopPt-genHadTopPt >> hKTHadTopPt(25, -100, 100)", "target==1 & bProbSSV>0.1");
  tK->Draw("hadTopPt-genHadTopPt >> hKPHadTopPt(25, -100, 100)", "fitProb>0.4 & bProbSSV>0.1");
  tK->Draw("hadTopPt-genHadTopPt >> hKBHadTopPt(25, -100, 100)", "combi==0 & bProbSSV>0.1");
  tK->Draw("hadTopPt-genHadTopPt >> hKAHadTopPt(25, -100, 100)", "bProbSSV>0.1");
  
  setHistStyle("p_{T,t,rec}-p_{T,t,gen}", hHTHadTopPt, kHitFit, 0);
  setHistStyle("p_{T,t,rec}-p_{T,t,gen}", hHPHadTopPt, kHitFit, 1);
  setHistStyle("p_{T,t,rec}-p_{T,t,gen}", hHBHadTopPt, kHitFit, 2);
  setHistStyle("p_{T,t,rec}-p_{T,t,gen}", hHAHadTopPt, kHitFit, 3);
  
  setHistStyle("p_{T,t,rec}-p_{T,t,gen}", hKTHadTopPt, kKinFit, 0);
  setHistStyle("p_{T,t,rec}-p_{T,t,gen}", hKPHadTopPt, kKinFit, 1);
  setHistStyle("p_{T,t,rec}-p_{T,t,gen}", hKBHadTopPt, kKinFit, 2);
  setHistStyle("p_{T,t,rec}-p_{T,t,gen}", hKAHadTopPt, kKinFit, 3);
  
  hHTHadTopPt->Draw("C");
  hHPHadTopPt->Draw("C,SAME");
  hHBHadTopPt->Draw("C,SAME");
  hHAHadTopPt->Draw("C,SAME");
  
  hKTHadTopPt->Draw("C,SAME");
  hKPHadTopPt->Draw("C,SAME");
  hKBHadTopPt->Draw("C,SAME");
  hKAHadTopPt->Draw("C,SAME");
  
  TLegend *leg0 = new TLegend(0.65, 0.55, 0.95, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  
  leg0->AddEntry( hHTHadTopPt , "HitFit correct", "F");
  leg0->AddEntry( hHPHadTopPt , "HitFit cut50", "F");
  leg0->AddEntry( hHBHadTopPt , "HitFit 'best'", "F");
  leg0->AddEntry( hHAHadTopPt , "HitFit all", "F");
  
  leg0->AddEntry( hKTHadTopPt , "KinFit correct", "F");
  leg0->AddEntry( hKPHadTopPt , "KinFit cut50", "F");
  leg0->AddEntry( hKBHadTopPt , "KinFit 'best'", "F");
  leg0->AddEntry( hKAHadTopPt , "KinFit all", "F");
  
  leg0->Draw();
  //*/
  
  //* top eta bias
  
  TCanvas* cHadTopEta = new TCanvas("cHadTopEta", "cHadTopEta", 450, 450);
  
  cHadTopEta->cd();
  
  tH->Draw("hadTopEta-genHadTopEta >> hHTHadTopEta(25, -0.5, 0.5)", "target==1 & bProbSSV>0.1");
  tH->Draw("hadTopEta-genHadTopEta >> hHPHadTopEta(25, -0.5, 0.5)", "hitFitProb>0.3 & bProbSSV>0.1");
  tH->Draw("hadTopEta-genHadTopEta >> hHBHadTopEta(25, -0.5, 0.5)", "combi==0 & bProbSSV>0.1");
  tH->Draw("hadTopEta-genHadTopEta >> hHAHadTopEta(25, -0.5, 0.5)", "bProbSSV>0.1");
  
  tK->Draw("hadTopEta-genHadTopEta >> hKTHadTopEta(25, -0.5, 0.5)", "target==1 & bProbSSV>0.1");
  tK->Draw("hadTopEta-genHadTopEta >> hKPHadTopEta(25, -0.5, 0.5)", "fitProb>0.4 & bProbSSV>0.1");
  tK->Draw("hadTopEta-genHadTopEta >> hKBHadTopEta(25, -0.5, 0.5)", "combi==0 & bProbSSV>0.1");
  tK->Draw("hadTopEta-genHadTopEta >> hKAHadTopEta(25, -0.5, 0.5)", "bProbSSV>0.1");
  
  setHistStyle("#eta_{t,rec}-#eta_{t,gen}", hHTHadTopEta, kHitFit, 0);
  setHistStyle("#eta_{t,rec}-#eta_{t,gen}", hHPHadTopEta, kHitFit, 1);
  setHistStyle("#eta_{t,rec}-#eta_{t,gen}", hHBHadTopEta, kHitFit, 2);
  setHistStyle("#eta_{t,rec}-#eta_{t,gen}", hHAHadTopEta, kHitFit, 3);
  
  setHistStyle("#eta_{t,rec}-#eta_{t,gen}", hKTHadTopEta, kKinFit, 0);
  setHistStyle("#eta_{t,rec}-#eta_{t,gen}", hKPHadTopEta, kKinFit, 1);
  setHistStyle("#eta_{t,rec}-#eta_{t,gen}", hKBHadTopEta, kKinFit, 2);
  setHistStyle("#eta_{t,rec}-#eta_{t,gen}", hKAHadTopEta, kKinFit, 3);
  
  hHTHadTopEta->Draw("C");
  hHPHadTopEta->Draw("C,SAME");
  hHBHadTopEta->Draw("C,SAME");
  hHAHadTopEta->Draw("C,SAME");
  
  hKTHadTopEta->Draw("C,SAME");
  hKPHadTopEta->Draw("C,SAME");
  hKBHadTopEta->Draw("C,SAME");
  hKAHadTopEta->Draw("C,SAME");
  
  TLegend *leg0 = new TLegend(0.65, 0.55, 0.95, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  
  leg0->AddEntry( hHTHadTopEta , "HitFit correct", "F");
  leg0->AddEntry( hHPHadTopEta , "HitFit cut50", "F");
  leg0->AddEntry( hHBHadTopEta , "HitFit 'best'", "F");
  leg0->AddEntry( hHAHadTopEta , "HitFit all", "F");
  
  leg0->AddEntry( hKTHadTopEta , "KinFit correct", "F");
  leg0->AddEntry( hKPHadTopEta , "KinFit cut50", "F");
  leg0->AddEntry( hKBHadTopEta , "KinFit 'best'", "F");
  leg0->AddEntry( hKAHadTopEta , "KinFit all", "F");
  
  leg0->Draw();
  //*/


}
