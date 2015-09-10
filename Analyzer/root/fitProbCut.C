#include <vector>
#include <iostream>

#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TGraph.h"

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

void fitProbCut()
{
  setTDRStyle();
  //gStyle->SetPadLeftMargin(0.2);
  //gStyle->SetTitleYOffset(1.75);

  // ---
  //    open input file
  // ---
  TChain* tH = new TChain("analyzeHitFit/eventTree");
  tH->Add("/nfs/dust/cms/user/mseidel/trees/Summer12_TTJetsMS1725_1.00_muon/job*.root");
  
  //* efficiency
  
  TCanvas* cEfficiency = new TCanvas("cEfficiency", "cEfficiency", 450, 450);
  cEfficiency->cd();
  //cEfficiency->SetLogy(1);
  
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("#chi^{2} Cut;#chi^{2} cut;");
  
  double aSB[25];
  double aCut[25];
  double aEffSig[25];
  double aPur[25];
  
  double bins = 25;
  
  double nSigTotal = tH->GetEntries("(weight.combinedWeight)*(top.combinationType==1)");
  double nBkgTotal = tH->GetEntries("(weight.combinedWeight)*(top.combinationType!=1)");
  
  for (int i = 1; i < bins; i++) {
    aCut[i-1] = 1./bins*i;
    
    for (int n = 0; n < chain->GetEntries(); ++n) {
    
    }
    
    TString sHSig("(weight.combinedWeight)*(top.combinationType==1 & top.fitProb > "); sHSig+= aCut[i-1]; sHSig+= ")";
    TString sHBkg("(weight.combinedWeight)*(top.combinationType!=1 & top.fitProb > "); sHBkg+= aCut[i-1]; sHBkg+= ")";
    
    std::cout << sHSig << std::endl;
    
    double nSig = tH->GetEntries(sHSig);
    double nBkg = tH->GetEntries(sHBkg);
    
    aSB[i-1] = nSig/sqrt(nBkg);
    aEffSig[i-1] = nSig/nSigTotal * 100;
    aPur[i-1] = nSig/(nSig+nBkg) * 100;
    std::cout << "nSig: " << nSig << std::endl;
    std::cout << "nBkg: " << nBkg << std::endl;
    std::cout << "S/sqrt(B): " << aSB[i-1] << std::endl;
  }
  
  TGraph* gH = new TGraph(bins-1, aCut, aSB);
  gH->SetLineColor(color_[kHitFit]);
  gH->SetLineWidth(2);
  mg->Add(gH);
  
  TGraph* gH2 = new TGraph(bins-1, aCut, aEffSig);
  gH2->SetLineColor(color_[kKinFit]);
  gH2->SetLineWidth(2);
  mg->Add(gH2);
  
  TGraph* gH3 = new TGraph(bins-1, aCut, aPur);
  gH3->SetLineColor(kBlue);
  gH3->SetLineWidth(2);
  mg->Add(gH3);
  
  mg->Draw("AC");
  mg->SetMinimum(0);
  
  drawcutline(3.218875825, 100.);
  
  //*
  TLegend *legEfficiency = new TLegend(0.6, 0.25, 0.9, 0.4);
  legEfficiency->SetFillStyle(0);
  legEfficiency->SetBorderSize(0);
  legEfficiency->AddEntry( gH, "S/#sqrt{B}", "L");
  legEfficiency->AddEntry( gH3, "S/(S+B) * 100", "L");
  legEfficiency->AddEntry( gH2, "#varepsilon_{S} * 100", "L");
  legEfficiency->Draw();
  //*/
}
