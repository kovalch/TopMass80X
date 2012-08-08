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
TChain* tTTJets100;
TChain* tTTJetsQ2Up;
TChain* tTTJetsQ2Down;
TChain* tTTJetsMatchingUp;
TChain* tTTJetsMatchingDown;
TChain* tTTFast100;
TChain* tTTFastMatchingDown;

enum lepton           { kElectron, kMuon, kAll};
TString lepton_ [3] = { "electron", "muon", "all"};

int channel = 2;

int nbins = 6;
double xfirst = 4;
double xlast  = 10;
TString sObs("jetMultiplicity");

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
  bool fast;
  sample(TChain* ch, const char* l, int c = kBlack, bool f = false)
  : chain(ch), label(l), color(c), fast(f) {}
};


std::vector<sample> samples;
std::vector<sample>::iterator it;

TH1F* responseProfile(TString sDraw, sample its) {
  sDraw += " >> h2"; sDraw += its.label; sDraw += sDraw(1, 7); sDraw += "("; sDraw += nbins; sDraw += ","; sDraw += xfirst; sDraw += ","; sDraw += xlast; sDraw += ")";
  TString sCut("(leptonPt > 30 & hadQBCSV<0.679 & hadQBarBCSV<0.679 & hadBBCSV>0.679 & lepBBCSV>0.679 & combi==0)*(MCWeight)");
  //if (its.line != 2) sCut += "*(MCWeight)";
  its.chain->Draw(sDraw, sCut);
  TString sH2("h2"); sH2 += its.label; sH2 += sDraw(1, 7);
  TH1F* h2 = (TH1F*) gDirectory->Get(sH2);
  /*
  for (int i = 0; i < nbins; ++i) {
    h2->SetBinError(h2->GetBinError()*2);
  }
  */
  //h2->Fit("gaus");
  //h2->Scale(1./h2->Integral());
  h2->SetLineColor(its.color);
  if (its.fast) h2->SetLineStyle(2);
  
  return h2;
}

TGraphErrors* ratioGraph(TH1F* p1, TH1F* p2, int color) {
  double position[nbins];
  double posError[nbins];
  double ratio[nbins];
  double ratioError[nbins];
  
  for (int i = 0; i < nbins; ++i) {
    position[i] = p2->GetBinCenter(i+1);
    posError[i] = (xlast-xfirst)/(2*nbins);
    
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
  //graph->Fit("constFit");
  
  return graph;
}

void theoryUncertainties()
{
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  
  TH1::SetDefaultSumw2(true);
  
  // Define chains
  TString sTree("analyzeHitFit/eventTree");
  tData               = new TChain(sTree);
  tTTJets100          = new TChain(sTree);
  tTTJetsQ2Up         = new TChain(sTree);
  tTTJetsQ2Down       = new TChain(sTree);
  tTTJetsMatchingUp   = new TChain(sTree);
  tTTJetsMatchingDown = new TChain(sTree);
  tTTFast100          = new TChain(sTree);
  tTTFastMatchingDown = new TChain(sTree);
  
  // ---
  //    open input files
  // ---
  if (channel == kMuon || channel == kAll) {
    tData               ->Add("/scratch/hh/current/cms/user/mseidel/Run2011_CRStudy_muon/analyzeTop.root");
    tTTJets100          ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_1.00_muon/analyzeTop.root");
    tTTJetsQ2Up         ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_scaleup_muon/analyzeTop.root");
    tTTJetsQ2Down       ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_scaledown_muon/analyzeTop.root");
    tTTJetsMatchingUp   ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_matchingup_muon/analyzeTop.root");
    tTTJetsMatchingDown ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_matchingdown_muon/analyzeTop.root");
    tTTFast100          ->Add("/scratch/hh/current/cms/user/mseidel/TT_FASTSIM_TuneZ2_madgraph_muon/analyzeTop.root");
    tTTFastMatchingDown ->Add("/scratch/hh/current/cms/user/mseidel/TT_FASTSIM_Z2_matchingdown_madgraph_muon/analyzeTop.root");
  }
  if (channel == kElectron || channel == kAll) {
    tData               ->Add("/scratch/hh/current/cms/user/mseidel/Run2011_CRStudy_electron/analyzeTop.root");
    tTTJets100          ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_1.00_electron/analyzeTop.root");
    tTTJetsQ2Up         ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_scaleup_electron/analyzeTop.root");
    tTTJetsQ2Down       ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_scaledown_electron/analyzeTop.root");
    tTTJetsMatchingUp   ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_matchingup_electron/analyzeTop.root");
    tTTJetsMatchingDown ->Add("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1725_matchingdown_electron/analyzeTop.root");
    tTTFast100          ->Add("/scratch/hh/current/cms/user/mseidel/TT_FASTSIM_TuneZ2_madgraph_electron/analyzeTop.root");
    tTTFastMatchingDown ->Add("/scratch/hh/current/cms/user/mseidel/TT_FASTSIM_Z2_matchingdown_madgraph_electron/analyzeTop.root");
  }
  
  samples.push_back(sample(tTTJets100, "Z2", kRed+1));
  samples.push_back(sample(tTTJetsQ2Up, "Q2 up", kGreen+1));
  samples.push_back(sample(tTTJetsQ2Down, "Q2 down", kBlue+1));
  samples.push_back(sample(tTTJetsMatchingUp, "MEPS up", kMagenta+1));
  samples.push_back(sample(tTTJetsMatchingDown, "MEPS down", kCyan+1));
  samples.push_back(sample(tTTFast100, "Fast Z2", kRed-7, 2));
  samples.push_back(sample(tTTFastMatchingDown, "Fast MEPS down", kCyan-7, 2));
  samples.push_back(sample(tData, "Data", kBlack));
    
  //"[0] + exp([1]+[2]*x)"
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    
    //it->profile  = responseProfile("hadTopMass", it->chain, it->label);
    it->profile  = responseProfile(sObs, *it);
    //it->profile  = responseProfile("(hadQRaw.Pt()+hadQBarRaw.Pt())/(hadQGen.Pt()+hadQBarGen.Pt()):(hadQGen.Pt()+hadQBarGen.Pt())/2", it->chain, it->label);
    //it->profileB = responseProfile("(lepBRaw.Pt()/lepBGen.Pt() + hadBRaw.Pt()/hadBGen.Pt())/2:(hadBGen.Pt()+lepBGen.Pt())/2", it->chain, it->label);
    
    //it->profile  = responseProfile("hadQRaw.Pt()/hadQGen.Pt():hadQGen.Pt()", it->chain, it->label);
    //it->profile  = responseProfile("hadWRaw.E()/hadWGen.E():hadWGen.Pt()", it->chain, it->label);
    //it->profileB = responseProfile("(hadBRaw.Pt()/hadBGen.Pt()):hadBGen.Pt()", it->chain, it->label);
    //it->profileB = responseProfile("(lepBRaw.Pt()/lepBGen.Pt()):lepBGen.Pt()", it->chain, it->label);
  }
  
  TMultiGraph *mg = new TMultiGraph();
  TString mgtitle(";"); mgtitle += sObs; mgtitle += ";ratio to Z2";
  mg->SetTitle(mgtitle);
  
  TLegend *leg1 = new TLegend(0.6, 0.725, 0.9, 0.925);
  leg1->SetTextSize(0.03);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    if (!it->fast) it->graph = ratioGraph(it->profile, samples.begin()->profile, it->color);
    else it->graph = ratioGraph(it->profile, samples[5].profile, it->color);
    it->graph->SetLineColor(it->color);
    if (it->fast) it->graph->SetLineStyle(2);
    it->graph->SetMarkerStyle(1);
    mg->Add(it->graph);
    leg1->AddEntry(it->graph, it->label, "PL");
  }
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    it->profile->Scale(1./it->profile->Integral());
  }
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
  
  mg->Draw("AP");
  //*
  mg->SetMinimum(0.6);
  mg->SetMaximum(1.4);
  //*/
  leg1->Draw();
  
  TString path("controlPlots_upDown/gen_"); path+= sObs; path += "_"; path += lepton_[channel]; path += "_ratio.eps";
  canvas->Print(path);
  
  TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 600, 600);
  canvas2->cd();
  canvas2->SetLogy(true);
  
  TH1F* hNull = new TH1F("null", "", nbins, xfirst, xlast);
  hNull->GetYaxis()->SetRangeUser(1e-3, samples[0].profile->GetMaximum()*10);
  hNull->GetXaxis()->SetTitle(sObs);
  hNull->Draw();
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    if (it->color == kBlack) it->profile->Draw("SAME,E");
    else it->profile->Draw("SAME,HIST");
  }
  
  leg1->Draw();
  
  TString path2("controlPlots_upDown/gen_"); path2 += sObs; path2 += "_"; path2 += lepton_[channel]; path2 += ".eps";
  canvas2->Print(path2);
}
