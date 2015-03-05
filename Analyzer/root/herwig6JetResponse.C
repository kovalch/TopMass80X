#include <string>
#include <vector>
#include <iostream>

#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"

#include "tdrstyle.C"

struct sample {
  TFile* file;
  TH2F* ud;
  TH2F* s;
  TH2F* c;
  TH2F* b;
  TH2F* glu;
  TH2F* nob;
  TH2F* all;
  TH1D* mean_ud;
  TH1D* mean_s;
  TH1D* mean_c;
  TH1D* mean_b;
  TH1D* mean_glu;
  TH1D* mean_nob;
  TH1D* mean_all;
  std::string filename;
  const char* label;
  int color;
  int linestyle;
  sample(const char* f, const char* l, int co = kBlack, int st = 1)
  : filename(f), label(l), color(co), linestyle(st) {}
};

struct flavor {
  const char* label;
  TH1D* mean_pythia;
  TH1D* mean_herwig;
  TH1D* ratio;
  TF1* ratioFit;
  TH1D* ratioAll;
  TF1* ratioAllFit;
  int color;
  int markerstyle;
  flavor(const char* l, TH1D* mc, TH1D* ml, int c = kBlack, int m = 20)
  : label(l), mean_pythia(mc), mean_herwig(ml), color(c), markerstyle(m) {}
};

std::vector<sample> samples;
std::vector<sample>::iterator it;
std::vector<flavor> flavors;
std::vector<flavor>::iterator itf;

TH1D* getProfile(TH2F* h, int c, int s) {
  TH1D* profile = h->ProfileX()->ProjectionX();
  profile->SetLineColor(c);
  profile->SetLineStyle(s);
  
  return profile;
}

void treeToResponseHist(std::string sample, Long64_t nentries = 1000000000) {
  TFile* histFile = new TFile((sample+std::string("_responseHist.root")).c_str(), "RECREATE");
  
  TChain* chain = new TChain("analyzeHitFit/eventTree");
  chain->Add((sample+std::string("_muon/job*.root")).c_str());
  chain->Add((sample+std::string("_electron/job*.root")).c_str());
  
  chain->Draw("jet.jet.Pt()/jet.genJet.Pt() : jet.genJet.Pt() >>  ud(50,0,250,99,0.02,2)", "top.fitProb[0]>0.2 & jet.jet.Pt()>30 & jet.bTagCSV<0.679 & Iteration$<4", "weight.combinedWeight", nentries);
  chain->Draw("jet.jet.Pt()/jet.genJet.Pt() : jet.genJet.Pt() >>   b(50,0,250,99,0.02,2)", "top.fitProb[0]>0.2 & jet.jet.Pt()>30 & jet.bTagCSV>0.679 & Iteration$<4", "weight.combinedWeight", nentries);
  
  /*
  chain->Draw("jet.jet.Pt()/jet.genJet.Pt() : jet.genJet.Pt() >>  ud(50,0,250,99,0.02,2)", "jet.jet.Pt()>30 & abs(jet.flavour)==1 || abs(jet.flavour)==2", "", nentries);
  chain->Draw("jet.jet.Pt()/jet.genJet.Pt() : jet.genJet.Pt() >>   s(50,0,250,99,0.02,2)", "jet.jet.Pt()>30 & abs(jet.flavour)==3", "", nentries);
  chain->Draw("jet.jet.Pt()/jet.genJet.Pt() : jet.genJet.Pt() >>   c(50,0,250,99,0.02,2)", "jet.jet.Pt()>30 & abs(jet.flavour)==4", "", nentries);
  chain->Draw("jet.jet.Pt()/jet.genJet.Pt() : jet.genJet.Pt() >>   b(50,0,250,99,0.02,2)", "jet.jet.Pt()>30 & abs(jet.flavour)==5", "", nentries);
  chain->Draw("jet.jet.Pt()/jet.genJet.Pt() : jet.genJet.Pt() >> glu(50,0,250,99,0.02,2)", "jet.jet.Pt()>30 & abs(jet.flavour)==21", "", nentries);
  */
  
  histFile->Write();
}

void herwig6JetResponse() {
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  // Samples
  //samples.push_back(sample("Summer12_TTJets1725_powheg_responseHist.root", "Powheg+Pythia6 Z2*", kRed+1));
  //samples.push_back(sample("Summer12_TTJets1725_powheg_herwig_responseHist.root", "Powheg+Herwig6 AUET2", kBlue+1));
  //samples.push_back(sample("Summer12_TTJets1725_FSIM_1.00_responseHist.root", "MadGraph+Pythia6 Z2*", kRed+1));
  samples.push_back(sample("Summer12_TT1725_amcatnlo_pythia8_responseHist.root", "aMC@NLO+Pythia8 4C", kGreen+1));
  //samples.push_back(sample("Summer12_TT1725_amcatnlo_herwigpp_responseHist.root", "aMC@NLO+Herwig++", kBlue+1));
  samples.push_back(sample("Summer12_TT1725_amcatnlo_herwigpp:pythia8_responseHist.root", "aMC@NLO+Herwig++", kBlue+1));
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    it->file = new TFile((std::string("/nfs/dust/cms/user/mseidel/trees_paper/")+it->filename).c_str());
    
    it->ud   = (TH2F*) it->file->Get("ud");
    //it->s    = (TH2F*) it->file->Get("s");
    //it->c    = (TH2F*) it->file->Get("c");
    it->b    = (TH2F*) it->file->Get("b");
    //it->glu  = (TH2F*) it->file->Get("glu");
    
    it->ud  ->Rebin2D(2,1);
    //it->s   ->Rebin2D(2,1);
    //it->c   ->Rebin2D(2,1);
    it->b   ->Rebin2D(2,1);
    //it->glu ->Rebin2D(2,1);
    
    //it->nob = (TH2F*) it->ud->Clone(); it->nob->SetName("nob");
    //it->nob->Add(it->s);
    //it->nob->Add(it->c);
    //it->nob->Add(it->glu);
    
    // TODO: Mixture buggy?
    it->all = (TH2F*) it->ud->Clone(); it->all->SetName("all");
    it->all->Add(it->b);
    
    it->mean_ud   = getProfile(it->ud  , it->color, it->linestyle);
    //it->mean_s    = getProfile(it->s   , it->color, it->linestyle);
    //it->mean_c    = getProfile(it->c   , it->color, it->linestyle);
    it->mean_b    = getProfile(it->b   , it->color, it->linestyle);
    //it->mean_glu  = getProfile(it->glu , it->color, it->linestyle);
    //it->mean_nob  = getProfile(it->nob , it->color, it->linestyle);
    it->mean_all  = getProfile(it->all , it->color, it->linestyle);
  }
  
  // Flavors
  flavors.push_back(flavor("t#bar{t} mixture", samples[0].mean_all, samples[1].mean_all, kRed-4, 20));
  flavors.push_back(flavor("light",      samples[0].mean_ud,   samples[1].mean_ud,   kBlue-4, 22));
  //flavors.push_back(flavor("strange", samples[0].mean_s,    samples[1].mean_s,    kOrange-4, 22));
  //flavors.push_back(flavor("charm",   samples[0].mean_c,    samples[1].mean_c,    kViolet-4, 22));
  flavors.push_back(flavor("bottom",  samples[0].mean_b,    samples[1].mean_b,    kBlack, 23));
  //flavors.push_back(flavor("gluon",   samples[0].mean_glu,  samples[1].mean_glu,  kGreen-3, 24));
  //flavors.push_back(flavor("W mix",   samples[0].mean_nob,  samples[1].mean_nob,  kGray, 26));
  
  TString fitFunction = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4";
  //TString fitFunction = "[0]";
  
  // Ratio Cluster/String
  TCanvas* canvas = new TCanvas("canvas", "canvas", 500, 500);
  TLegend *leg = new TLegend(0.6, 0.65, 0.9, 0.9);
  leg->SetTextSize(0.05);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  
  for (itf = flavors.begin(); itf != flavors.end(); ++itf) {
    itf->ratio = (TH1D*) itf->mean_herwig->Clone();
    itf->ratio->Divide(itf->mean_pythia);
    itf->ratio->SetLineColor(itf->color);
    itf->ratio->SetFillColor(itf->color);
    itf->ratio->SetMarkerColor(itf->color);
    itf->ratio->SetMarkerStyle(itf->markerstyle);
    itf->ratio->GetXaxis()->SetRangeUser(20, 245);
    itf->ratio->GetYaxis()->SetRangeUser(0.95, 1.1);
    itf->ratio->SetTitle("Ratio;GenJet p_{T};Jet response, Herwig/Pythia");
    
    leg->AddEntry(itf->ratio, itf->label, "LP");
    
    // FIT
    itf->ratioFit = new TF1((std::string(itf->label)+std::string("Fit")).c_str(), fitFunction);
    itf->ratioFit->SetLineColor(itf->color);
    itf->ratioFit->SetLineWidth(2);
    itf->ratioFit->SetRange(20, 250);
    
    itf->ratio->Fit((std::string(itf->label)+std::string("Fit")).c_str(), "FEMR");
  }
  
  for (itf = flavors.begin(); itf != flavors.end(); ++itf) {
    if (itf == flavors.begin()) {
      itf->ratio->Draw();
    }
    else itf->ratio->Draw("SAME");
  }
  
  leg->Draw();
  canvas->Print("amcatnloGenJetResponse.eps");
  
  // Ratio Flavor/mixture
  TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 500, 500);
  
  for (itf = flavors.begin(); itf != flavors.end(); ++itf) {
    itf->ratioAll = (TH1D*) itf->ratio->Clone();
    itf->ratioAll->Divide(flavors[0].ratio);
    itf->ratioAll->GetXaxis()->SetRangeUser(20, 245);
    itf->ratioAll->GetYaxis()->SetRangeUser(0.95, 1.1);
    itf->ratioAll->SetTitle("Ratio;GenJet p_{T};Rel. Jet response, Herwig/Pythia");
    
    // FIT
    itf->ratioAllFit = new TF1((std::string(itf->label)+std::string("AllFit")).c_str(), fitFunction);
    itf->ratioAllFit->SetLineColor(itf->color);
    itf->ratioAllFit->SetLineWidth(2);
    itf->ratioAllFit->SetRange(20, 250);
    
    itf->ratioAll->Fit((std::string(itf->label)+std::string("AllFit")).c_str(), "FEMR");
  }
  
  for (itf = flavors.begin(); itf != flavors.end(); ++itf) {
    if (itf == flavors.begin()) {
      itf->ratioAll->Draw();
    }
    else itf->ratioAll->Draw("SAME");
  }
  
  leg->Draw();
  canvas2->Print("amcatnloJetResponse_ratio.eps");
}
