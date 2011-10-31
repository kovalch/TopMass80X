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
  
  hist->SetLineColor(col);
  hist->SetLineWidth(2);
  
  if (kStyle == kGenMatch) {
    hist->SetLineStyle(7);
    hist->SetLineWidth(4);
  }
  
  //hist->GetYaxis()->SetTitle("Fraction of permutations");
  hist->GetXaxis()->SetTitle(titleX);
  hist->GetYaxis()->SetTitle("Fraction of permutations");
}

void hitFit()
{
  setTDRStyle();
  //gStyle->SetPadLeftMargin(0.2);
  //gStyle->SetTitleYOffset(1.75);

  // ---
  //    open input file
  // ---
  TFile* fTTJets = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_1.00_2b/analyzeTop.root");
  
  // ---
  //    Get trees
  // ---
  TTree* tH = (TTree*) fTTJets->Get("analyzeHitFit/eventTree");
  TTree* tG = (TTree*) fTTJets->Get("analyzeGenMatch/eventTree");
  
  //* efficiency
  
  TCanvas* cEfficiency = new TCanvas("cEfficiency", "cEfficiency", 450, 450);
  cEfficiency->cd();
  cEfficiency->SetLogy(1);
  
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("Kinematic Fit Efficiency;Correct permutation efficiency;Bad permutation efficiency");
  
  double aHSigEfficiency[10];
  double aHBkgEfficiency[10];
  
  double bins = 10;
  
  for (int i = 1; i < bins; i++) {
    TString sHSig("target == 1 & hitFitProb > "); sHSig+= 1./bins * (double)i;
    TString sHBkg("target != 1 & hitFitProb > "); sHBkg+= 1./bins * (double)i;
    
    std::cout << sHSig << std::endl;
    
    aHSigEfficiency[i-1] = (double)tH->GetEntries(sHSig)/(double)tH->GetEntries("target == 1");
    aHBkgEfficiency[i-1] = (double)tH->GetEntries(sHBkg)/(double)tH->GetEntries("target != 1");
  }
  
  gH = new TGraph(bins-1, aHSigEfficiency, aHBkgEfficiency);
  gH->SetLineColor(color_[kHitFit]);
  gH->SetLineWidth(2);
  mg->Add(gH);
  
  mg->Draw("AC");
  
  TLegend *legEfficiency = new TLegend(0.25, 0.75, 0.55, 0.9);
  legEfficiency->SetFillStyle(0);
  legEfficiency->SetBorderSize(0);
  legEfficiency->AddEntry( gH, "HitFit"   , "L");
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
  
  /* top mass resolution
  
  TCanvas* cHadTopMassRes = new TCanvas("cHadTopMassRes", "cHadTopMassRes", 450, 450);
  
  cHadTopMassRes->cd();
  
  tH->Draw("hadTopMass-genHadTopMass >> hHHadTopMassRes(25, -100, 100)", "target==1");
  tG->Draw("hadTopMass-genHadTopMass >> hGHadTopMassRes(25, -100, 100)", "");
  
  setHistStyle("m_{t,rec}-m_{t,gen}", hHHadTopMassRes, kHitFit, 0);
  setHistStyle("m_{t,rec}-m_{t,gen}", hGHadTopMassRes, kGenMatch, 0);
  
  hHHadTopMassRes->Draw("C");
  hGHadTopMassRes->Draw("C,SAME");
  
  TLegend *leg0 = new TLegend(0.65, 0.55, 0.95, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  leg0->AddEntry( hHHadTopMassRes , "HitFit", "F");
  leg0->AddEntry( hGHadTopMassRes , "GenMatch", "F");
  leg0->Draw();
  //*/
  
  //* top mass bias
  
  TCanvas* cHadTopMass = new TCanvas("cHadTopMass", "cHadTopMass", 450, 450);
  
  cHadTopMass->cd();
  
  tH->Draw("hadTopMass-genHadTopMass >> hHTHadTopMass(25, -100, 100)", "target==1");
  tH->Draw("hadTopMass-genHadTopMass >> hHPHadTopMass(25, -100, 100)", "hitFitProb>0.2");
  tH->Draw("hadTopMass-genHadTopMass >> hHBHadTopMass(25, -100, 100)", "combi==0");
  tH->Draw("hadTopMass-genHadTopMass >> hHAHadTopMass(25, -100, 100)", "");
  tH->Draw("hadTopRawMass-genHadTopMass >> hHRHadTopMass(25, -100, 100)", "");
  tG->Draw("hadTopMass-genHadTopMass >> hGHadTopMass(25, -100, 100)", "");
  
  setHistStyle("m_{t,rec}-m_{t,gen}", hHTHadTopMass, kHitFit, kRed+1);
  setHistStyle("m_{t,rec}-m_{t,gen}", hHPHadTopMass, kHitFit, kMagenta+1);
  setHistStyle("m_{t,rec}-m_{t,gen}", hHBHadTopMass, kHitFit, kBlue+1);
  setHistStyle("m_{t,rec}-m_{t,gen}", hHAHadTopMass, kHitFit, kGreen+1);
  setHistStyle("m_{t,rec}-m_{t,gen}", hHRHadTopMass, kHitFit, kYellow+1);
  setHistStyle("m_{t,rec}-m_{t,gen}", hGHadTopMass, kGenMatch, kBlack);
  
  hHTHadTopMass->Draw("C");
  hHPHadTopMass->Draw("C,SAME");
  hHBHadTopMass->Draw("C,SAME");
  hHAHadTopMass->Draw("C,SAME");
  hHRHadTopMass->Draw("C,SAME");
  hGHadTopMass->Draw("C,SAME");
  
  TLegend *leg0 = new TLegend(0.65, 0.55, 0.95, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  
  leg0->AddEntry( hHTHadTopMass , "HitFit correct", "F");
  leg0->AddEntry( hHPHadTopMass , "HitFit P_{fit}>0.2", "F");
  leg0->AddEntry( hHBHadTopMass , "HitFit lowest #chi^{2}", "F");
  leg0->AddEntry( hHAHadTopMass , "HitFit all", "F");
  leg0->AddEntry( hHRHadTopMass , "No kinematic fit", "F");
  leg0->AddEntry( hGHadTopMass  , "GenMatch", "F");
  
  leg0->Draw();
  //*/
  
  //* top pt bias
  
  TCanvas* cHadTopPt = new TCanvas("cHadTopPt", "cHadTopPt", 450, 450);
  
  cHadTopPt->cd();
  
  tH->Draw("hadTopPt-genHadTopPt >> hHTHadTopPt(25, -100, 100)", "target==1");
  tH->Draw("hadTopPt-genHadTopPt >> hHPHadTopPt(25, -100, 100)", "hitFitProb>0.2");
  tH->Draw("hadTopPt-genHadTopPt >> hHBHadTopPt(25, -100, 100)", "combi==0");
  tH->Draw("hadTopPt-genHadTopPt >> hHAHadTopPt(25, -100, 100)", "");
  tH->Draw("hadTopRawPt-genHadTopPt >> hHRHadTopPt(25, -100, 100)", "");
  tG->Draw("hadTopPt-genHadTopPt >> hGHadTopPt(25, -100, 100)", "");
  
  setHistStyle("p_{T,t,rec}-p_{T,t,gen}", hHTHadTopPt, kHitFit, kRed+1);
  setHistStyle("p_{T,t,rec}-p_{T,t,gen}", hHPHadTopPt, kHitFit, kMagenta+1);
  setHistStyle("p_{T,t,rec}-p_{T,t,gen}", hHBHadTopPt, kHitFit, kBlue+1);
  setHistStyle("p_{T,t,rec}-p_{T,t,gen}", hHAHadTopPt, kHitFit, kGreen+1);
  setHistStyle("p_{T,t,rec}-p_{T,t,gen}", hHRHadTopPt, kHitFit, kYellow+1);
  setHistStyle("p_{T,t,rec}-p_{T,t,gen}", hGHadTopPt , kGenMatch, kBlack);
  
  hHTHadTopPt->Draw("C");
  hHPHadTopPt->Draw("C,SAME");
  hHBHadTopPt->Draw("C,SAME");
  hHAHadTopPt->Draw("C,SAME");
  hHRHadTopPt->Draw("C,SAME");
  hGHadTopPt->Draw("C,SAME");
  
  TLegend *leg0 = new TLegend(0.65, 0.55, 0.95, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  
  leg0->AddEntry( hHTHadTopPt , "HitFit correct", "F");
  leg0->AddEntry( hHPHadTopPt , "HitFit P_{fit}>0.2", "F");
  leg0->AddEntry( hHBHadTopPt , "HitFit lowest #chi^{2}", "F");
  leg0->AddEntry( hHAHadTopPt , "HitFit all", "F");
  leg0->AddEntry( hHRHadTopPt , "No kinematic fit", "F");
  leg0->AddEntry( hGHadTopPt  , "GenMatch", "F");
  
  leg0->Draw();
  //*/
  
  //* top eta bias
  
  TCanvas* cHadTopEta = new TCanvas("cHadTopEta", "cHadTopEta", 450, 450);
  
  cHadTopEta->cd();
  
  tH->Draw("hadTopEta-genHadTopEta >> hHTHadTopEta(25, -0.5, 0.5)", "target==1");
  tH->Draw("hadTopEta-genHadTopEta >> hHPHadTopEta(25, -0.5, 0.5)", "hitFitProb>0.2");
  tH->Draw("hadTopEta-genHadTopEta >> hHBHadTopEta(25, -0.5, 0.5)", "combi==0");
  tH->Draw("hadTopEta-genHadTopEta >> hHAHadTopEta(25, -0.5, 0.5)", "");
  tH->Draw("hadTopRawEta-genHadTopEta >> hHRHadTopEta(25, -0.5, 0.5)", "");
  tG->Draw("hadTopEta-genHadTopEta >> hGHadTopEta(25, -0.5, 0.5)", "");
  
  setHistStyle("#eta_{t,rec}-#eta_{t,gen}", hHTHadTopEta, kHitFit, kRed+1);
  setHistStyle("#eta_{t,rec}-#eta_{t,gen}", hHPHadTopEta, kHitFit, kMagenta+1);
  setHistStyle("#eta_{t,rec}-#eta_{t,gen}", hHBHadTopEta, kHitFit, kBlue+1);
  setHistStyle("#eta_{t,rec}-#eta_{t,gen}", hHAHadTopEta, kHitFit, kGreen+1);
  setHistStyle("#eta_{t,rec}-#eta_{t,gen}", hHRHadTopEta, kHitFit, kYellow+1);
  setHistStyle("#eta_{t,rec}-#eta_{t,gen}", hGHadTopEta, kGenMatch, kBlack);
  
  hHTHadTopEta->Draw("C");
  hHPHadTopEta->Draw("C,SAME");
  hHBHadTopEta->Draw("C,SAME");
  hHAHadTopEta->Draw("C,SAME");
  hHRHadTopEta->Draw("C,SAME");
  hGHadTopEta->Draw("C,SAME");
  
  TLegend *leg0 = new TLegend(0.65, 0.55, 0.95, 0.9);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  
  leg0->AddEntry( hHTHadTopEta , "HitFit correct", "F");
  leg0->AddEntry( hHPHadTopEta , "HitFit P_{fit}>0.2", "F");
  leg0->AddEntry( hHBHadTopEta , "HitFit lowest #chi^{2}", "F");
  leg0->AddEntry( hHAHadTopEta , "HitFit all", "F");
  leg0->AddEntry( hHRHadTopEta , "No kinematic fit", "F");
  leg0->AddEntry( hGHadTopEta  , "GenMatch", "F");
  
  leg0->Draw();
  //*/


}
