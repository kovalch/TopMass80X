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

void fitProbCut()
{
  setTDRStyle();
  //gStyle->SetPadLeftMargin(0.2);
  //gStyle->SetTitleYOffset(1.75);

  // ---
  //    open input file
  // ---
  TFile* fTTJets = new TFile("/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_1.00/analyzeTop.root");
  
  // ---
  //    Get trees
  // ---
  TTree* tH = (TTree*) fTTJets->Get("analyzeHitFit/eventTree");
  
  //* efficiency
  
  TCanvas* cEfficiency = new TCanvas("cEfficiency", "cEfficiency", 450, 450);
  cEfficiency->cd();
  //cEfficiency->SetLogy(1);
  
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("#chi^{2} Cut;#chi^{2} cut;S / #sqrt{S+B}");
  
  double aSB[10];
  double aCut[10];
  
  double bins = 10;
  
  for (int i = 1; i < bins; i++) {
    aCut[i-1] = 2*i;
    
    TString sHSig("(MCWeight)*(target==1 & leptonPt > 30 & hadBBSSV>1.74 & lepBBSSV>1.74 & hadQBSSV<1.74 & hadQBarBSSV<1.74 & hitFitChi2 < "); sHSig+= aCut[i-1]; sHSig+= ")";
    TString sHBkg("(MCWeight)*(target!=1 & leptonPt > 30 & hadBBSSV>1.74 & lepBBSSV>1.74 & hadQBSSV<1.74 & hadQBarBSSV<1.74 & hitFitChi2 < "); sHBkg+= aCut[i-1]; sHBkg+= ")";
    
    std::cout << sHSig << std::endl;
    
    double nSig = tH->GetEntries(sHSig);
    double nBkg = tH->GetEntries(sHBkg);
    
    aSB[i-1] = nSig/sqrt(nSig+nBkg);
    std::cout << "nSig: " << nSig << std::endl;
    std::cout << "nBkg: " << nBkg << std::endl;
    std::cout << "S/sqrt(S+B): " << aSB[i-1] << std::endl;
  }
  
  gH = new TGraph(bins-1, aCut, aSB);
  gH->SetLineColor(color_[kHitFit]);
  gH->SetLineWidth(2);
  mg->Add(gH);
  
  mg->Draw("AC");
  
  /*
  TLegend *legEfficiency = new TLegend(0.25, 0.75, 0.55, 0.9);
  legEfficiency->SetFillStyle(0);
  legEfficiency->SetBorderSize(0);
  legEfficiency->AddEntry( gH, "HitFit"   , "L");
  legEfficiency->Draw();
  //*/
}
