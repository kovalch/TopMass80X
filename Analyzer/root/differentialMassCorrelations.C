#include <string>
#include <vector>
#include <iostream>

#include "TFile.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TF1.h"
#include "TLegend.h"
#include "TChain.h"
#include "TGraph.h"

#include "tdrstyle.C"

struct obs {
  TString name;
  TString formula;
  obs(TString n, TString f) : name(n), formula(f) {}
};

std::vector<obs> vobs;

void differentialMassCorrelations() {
  TH1::SetDefaultSumw2();
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  vobs.push_back(obs("top pt", "top.fitTop1.Pt()")); // 1 x
  vobs.push_back(obs("top eta", "top.fitTop1.Eta()")); // 2
  vobs.push_back(obs("W pt", "top.fitW1.Pt()")); // 3
  vobs.push_back(obs("W eta", "top.fitW1.Eta()")); // 4
  vobs.push_back(obs("b pt", "top.fitB1.Pt()")); // 5 x
  vobs.push_back(obs("b eta", "top.fitB1.Eta()")); // 6 x
  vobs.push_back(obs("q pt", "top.fitW1Prod1.Pt()")); // 7
  vobs.push_back(obs("q eta", "top.fitW1Prod1.Eta()")); // 8
  vobs.push_back(obs("ttbar mass", "top.fitTTBar.M()")); // 9 x
  vobs.push_back(obs("ttbar pt", "top.fitTTBar.Pt()")); // )10 x
  vobs.push_back(obs("dR qq", "sqrt(pow(top.fitW1Prod1.Eta()-top.fitW1Prod2.Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitW1Prod1.Phi()-top.fitW1Prod2.Phi()),2))")); // 11 dRqq x
  vobs.push_back(obs("dR bb", "sqrt(pow(top.fitB1.Eta()-top.fitB2.Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1.Phi()-top.fitB2.Phi()),2))")); // 12 dRbb x
  vobs.push_back(obs("HT", "top.fitB1.Pt()+top.fitB2.Pt()+top.fitW1Prod1.Pt()+top.fitW1Prod2.Pt()")); // 13 HT
  vobs.push_back(obs("N(jet)", "Max$(Iteration$*(jet.jet.Pt()>30)) + 1")); // 15 x
  vobs.push_back(obs("alpha", "jet.jet[4].Pt()/(top.recoB1[0].Pt()+top.recoB2[0].Pt())*2"));
  
  TLegend *leg = new TLegend(0.25, 0.7, 0.55, 0.9);
  leg->SetTextSize(0.03);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  
  TH2D* corr = new TH2D("corr", "corr", vobs.size(), 0, vobs.size(), vobs.size(), 0, vobs.size());
  
  TChain* chain = new TChain("analyzeHitFit/eventTree");
  chain->Add("/nfs/dust/cms/user/mseidel/trees_paper/Summer12_TTJetsMS1725_1.00_muon/job*.root");
  
  for (unsigned int i = 0; i < vobs.size(); ++i) {
    for (unsigned int j = 0; j < vobs.size(); ++j) {
      std::cout << vobs[i].name << ":" << vobs[j].name << std::endl;
      chain->Draw(vobs[i].formula+":"+vobs[j].formula, "(top.fitProb>0.2)*weight.combinedWeight", "", 10000);
      TGraph *graph = (TGraph*)gPad->GetPrimitive("Graph");
      
      corr->SetBinContent(i+1, j+1, graph->GetCorrelationFactor());
    }
    corr->GetXaxis()->SetBinLabel(i+1, vobs[i].name);
    corr->GetYaxis()->SetBinLabel(i+1, vobs[i].name);
  }
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
  corr->GetZaxis()->SetRangeUser(-1., 1.);
  corr->Draw("COLZ");
  
  canvas->Print("differentialMassCorrelations.eps");
}
