#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TChain.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"

#include "tdrstyle.C"

TChain* tData;
TChain* tTTJets096;
TChain* tTTJets100;
TChain* tTTJets104;
TChain* tTTJetsFlavorUp;
TChain* tTTJetsFlavorDown;
TChain* tTTJetsP11;
TChain* tTTJetsP11noCR;

enum lepton           { kElectron, kMuon, kAll};
TString lepton_ [3] = { "electron", "muon", "all"};

int channel = 1;

int nbins = 20;
double xfirst = 0;
double xlast  = 3.1416;

int color_ [] =  {kRed+1, kRed-7, kRed-10, kGray , kGreen+1, kRed+1 , kAzure-2, kGreen-3, 
                 kYellow , kMagenta, 10      , kBlack  , 
                 kYellow , kYellow , kYellow , kYellow , kYellow , kYellow ,
                 10      , 10      , 10      , 
                 kMagenta, kMagenta, kMagenta, kMagenta, kMagenta, kMagenta };

int marker_[] = {20, 22, 22, 20, 22, 23, 29, 23, 
                 21, 27, 28, 20, 
                 21, 21, 21, 21, 21, 21,
                 28, 28, 28,
                 27, 27, 27, 27, 27, 27};

struct sample {
  TChain* chain;
  TH1F* profile;
  TGraphErrors* graph;
  const char* label;
  int color;
  double scale;
  sample(TChain* ch, const char* l, int c = kBlack, double s = 1.)
  : chain(ch), label(l), color(c), scale(s) {}
};


std::vector<sample> samples;
std::vector<sample>::iterator it;

TH1F* responseProfile(TString sDraw, TChain* chain, TString label) {
  sDraw += " >> h2"; sDraw += label; sDraw += sDraw(0, 3); sDraw += "("; sDraw += nbins; sDraw += ","; sDraw += xfirst; sDraw += ","; sDraw += xlast; sDraw += ")";
  chain->Draw(sDraw, "(top.fitProb>0.2)*(weight.combinedWeight)");
  TString sH2("h2"); sH2 += label; sH2 += sDraw(0, 3);
  TH1F* h2 = (TH1F*) gDirectory->Get(sH2);
  /*
  for (int i = 0; i < nbins; ++i) {
    h2->SetBinError(h2->GetBinError()*2);
  }
  */
  //h2->Fit("gaus");
  //h2->Scale(1./h2->Integral());
  return h2;
}

TGraphErrors* ratioGraph(TH1F* p1, TH1F* p2, int color) {
  double position[nbins];
  double posError[nbins];
  double ratio[nbins];
  double ratioError[nbins];
  
  for (int i = 0; i < nbins; ++i) {
    position[i] = 1./(2.*nbins) + 1./nbins*i;
    posError[i] = 1./(4.*nbins);
    
    if (p2->GetBinContent(i+1)>0) {
      ratio[i]      = p1->GetBinContent(i+1)/p2->GetBinContent(i+1) * p2->Integral()/p1->Integral();
      ratioError[i] = sqrt(pow(p1->GetBinError(i+1)/p2->GetBinContent(i+1), 2) + pow(p1->GetBinContent(i+1)/pow(p2->GetBinContent(i+1), 2) * p2->GetBinError(i+1), 2)) * p2->Integral()/p1->Integral();
    }
    /*
    else {
      ratio[i] = 1;
      ratioError[i] = 0.00001;
    }
    */
  }
  
  TGraphErrors* graph = new TGraphErrors(nbins, position, ratio, posError, ratioError);
  
  TF1* constFit = new TF1("constFit", "[0]+x*[1]");
  constFit->SetParNames("offset");
  constFit->SetLineColor(color);
  constFit->SetLineWidth(1);
  constFit->SetLineStyle(7);
  graph->Fit("constFit");
  
  return graph;
}

void CRStudy()
{
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  
  // Define chains
  TString sTree("analyzeHitFit/eventTree");
  tData             = new TChain(sTree);
  tTTJets100        = new TChain(sTree);
  tTTJetsP11        = new TChain(sTree);
  tTTJetsP11noCR    = new TChain(sTree);
  
  
  // ---
  //    open input files
  // ---
  if (channel == kMuon || channel == kAll) {
    tData             ->Add("/scratch/hh/dust/naf/cms/user/mseidel/trees/Run2012_muon/job_*.root");
    tTTJets100        ->Add("/scratch/hh/dust/naf/cms/user/mseidel/trees/Summer12_TTJets1725_1.00_muon/job_*.root");
    tTTJetsP11        ->Add("/scratch/hh/dust/naf/cms/user/mseidel/trees/Summer12_TTJets1725_SemiLept_P11_muon/job_*.root");
    tTTJetsP11noCR    ->Add("/scratch/hh/dust/naf/cms/user/mseidel/trees/Summer12_TTJets1725_SemiLept_P11noCR_muon/job_*.root");
  }
  if (channel == kElectron || channel == kAll) {
    tData             ->Add("/scratch/hh/dust/naf/cms/user/mseidel/trees/Run2011_CRStudy_electron/analyzeTop.root");
    tTTJetsP11        ->Add("/scratch/hh/dust/naf/cms/user/mseidel/trees/Fall11_TTJets1725_P11_CRStudy_electron/analyzeTop.root");
    tTTJetsP11noCR    ->Add("/scratch/hh/dust/naf/cms/user/mseidel/trees/Fall11_TTJets1725_P11noCR_CRStudy_electron/analyzeTop.root");
  }
  
  samples.push_back(sample(tData, "Data", kBlack));
  samples.push_back(sample(tTTJets100, "Z2", kRed+1, 19.7));
  samples.push_back(sample(tTTJetsP11, "P11", kMagenta+1, 19.7*9./4.));
  samples.push_back(sample(tTTJetsP11noCR, "P11noCR", kCyan+1, 19.7*9./4.));
    
  //"[0] + exp([1]+[2]*x)"
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    
    //it->profile  = responseProfile("hadTopMass", it->chain, it->label);
    it->profile  = responseProfile("abs(TVector2::Phi_mpi_pi(jet.pullCharged[top.recoJetIdxW1Prod1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod2.Phi()-top.fitW1Prod1.Phi()),-(top.fitW1Prod2.Eta()-top.fitW1Prod1.Eta())))))", it->chain, it->label);
    //it->profile  = responseProfile("abs(TVector2::Phi_mpi_pi(jet.pull[top.recoJetIdxW1Prod2].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitW1Prod1.Phi()-top.fitW1Prod2.Phi()),-(top.fitW1Prod1.Eta()-top.fitW1Prod2.Eta())))))", it->chain, it->label);
    //it->profile  = responseProfile("abs(TVector2::Phi_mpi_pi(jet.pull[top.recoJetIdxB1].Phi()-(TMath::Pi()+TMath::ATan2(-(top.fitB2.Phi()-top.fitB1.Phi()),-(top.fitB2.Eta()-top.fitB1.Eta())))))", it->chain, it->label);
    std::cout << it->profile->GetEntries() << std::endl;
    //it->profile  = responseProfile("abs(TVector2::Phi_mpi_pi(jet[recoJetIdxB1].pull.Phi()-(TMath::Pi()+TMath::ATan2(-(fitB2.Phi()-fitB1.Phi()),-(fitB2.Eta()-fitB1.Eta())))))", it->chain, it->label);
    
    //it->profile  = responseProfile("(hadQRaw.Pt()+hadQBarRaw.Pt())/(hadQGen.Pt()+hadQBarGen.Pt()):(hadQGen.Pt()+hadQBarGen.Pt())/2", it->chain, it->label);
    //it->profileB = responseProfile("(lepBRaw.Pt()/lepBGen.Pt() + hadBRaw.Pt()/hadBGen.Pt())/2:(hadBGen.Pt()+lepBGen.Pt())/2", it->chain, it->label);
    
    //it->profile  = responseProfile("hadQRaw.Pt()/hadQGen.Pt():hadQGen.Pt()", it->chain, it->label);
    //it->profile  = responseProfile("hadWRaw.E()/hadWGen.E():hadWGen.Pt()", it->chain, it->label);
    //it->profileB = responseProfile("(hadBRaw.Pt()/hadBGen.Pt()):hadBGen.Pt()", it->chain, it->label);
    //it->profileB = responseProfile("(lepBRaw.Pt()/lepBGen.Pt()):lepBGen.Pt()", it->chain, it->label);
  }
  
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle(";[bin];ratio. to P11");
  
  TLegend *leg1 = new TLegend(0.6, 0.75, 0.9, 0.925);
  leg1->SetTextSize(0.03);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    it->profile->SetLineColor(it->color);
    if (it == samples.begin()) { // DATA
      it->profile->SetMarkerStyle(20);
      it->profile->GetYaxis()->SetRangeUser(1, it->profile->GetMaximum()*1.5);
      it->profile->Draw("E");
    }
    else {
      it->profile->SetMarkerStyle(1);
      it->profile->Scale(it->scale);
      it->profile->Draw("SAME");
    }
    leg1->AddEntry(it->profile, it->label, "PL");
  }
  
  /*
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle(";[bin];ratio. to P11");
  
  TLegend *leg1 = new TLegend(0.6, 0.75, 0.9, 0.925);
  leg1->SetTextSize(0.03);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    it->graph = ratioGraph(it->profile, samples.begin()->profile, it->color);
    it->graph->SetLineColor(it->color);
    it->graph->SetMarkerStyle(1);
    mg->Add(it->graph);
    leg1->AddEntry(it->graph, it->label, "PL");
  }
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
  
  mg->Draw("AP");
  //*
  mg->SetMinimum(0.6);
  mg->SetMaximum(1.4);
  //*/
  leg1->Draw();
}
