#include <string>
#include <vector>
#include <iostream>

#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"

#include "tdrstyle.C"

struct sample {
  TFile* file;
  TH2F* udsc;
  TH2F* b;
  TH2F* glu;
  TH2F* all;
  TH1D* mean_udsc;
  TH1D* mean_b;
  TH1D* mean_glu;
  TH1D* mean_all;
  std::string filename;
  const char* label;
  int color;
  int linestyle;
  sample(const char* f, const char* l, int c = kBlack, int s = 1)
  : filename(f), label(l), color(c), linestyle(s) {}
};

struct flavor {
  const char* label;
  TH1D* mean_cluster;
  TH1D* mean_string;
  TH1D* ratio;
  TF1* ratioFit;
  TH1D* ratioAll;
  TF1* ratioAllFit;
  int color;
  int markerstyle;
  flavor(const char* l, TH1D* mc, TH1D* ml, int c = kBlack, int m = 20)
  : label(l), mean_cluster(mc), mean_string(ml), color(c), markerstyle(m) {}
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

void sherpaHadResponse() {
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  // Samples
  samples.push_back(sample("TT_Cluster05_8TeV-sherpa.root", "Sherpa 2.1, Cluster", kRed+1));
  samples.push_back(sample("TT_Lund05_8TeV-sherpa.root", "Sherpa 2.1, Lund", kBlue+1));
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    it->file = new TFile((std::string("/nfs/dust/cms/user/mseidel/batch_output/analyzeSherpaGenEvent/")+it->filename).c_str());
    
    it->udsc = (TH2F*) it->file->Get("udsc");
    it->b    = (TH2F*) it->file->Get("b");
    it->glu  = (TH2F*) it->file->Get("glu");
    
    it->udsc->Rebin2D(2,1);
    it->b   ->Rebin2D(2,1);
    it->glu ->Rebin2D(2,1);
    
    it->all  = (TH2F*) it->udsc->Clone(); it->all->SetName("all");
    it->all->Add(it->b); it->all->Add(it->glu);
    
    it->mean_udsc = getProfile(it->udsc, it->color, it->linestyle);
    it->mean_b    = getProfile(it->b   , it->color, it->linestyle);
    it->mean_glu  = getProfile(it->glu , it->color, it->linestyle);
    it->mean_all  = getProfile(it->all , it->color, it->linestyle);
  }
  
  // Flavors
  flavors.push_back(flavor("t#bar{t} mixture", samples[0].mean_all, samples[1].mean_all, kRed-4, 20));
  flavors.push_back(flavor("udsc",   samples[0].mean_udsc, samples[1].mean_udsc, kBlue-4, 22));
  flavors.push_back(flavor("bottom", samples[0].mean_b, samples[1].mean_b, kBlack, 23));
  flavors.push_back(flavor("gluon",  samples[0].mean_glu, samples[1].mean_glu, kGreen-3, 24));
  
  // Ratio Cluster/String
  TCanvas* canvas = new TCanvas("canvas", "canvas", 500, 500);
  TLegend *leg = new TLegend(0.6, 0.65, 0.9, 0.9);
  leg->SetTextSize(0.05);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  
  for (itf = flavors.begin(); itf != flavors.end(); ++itf) {
    itf->ratio = (TH1D*) itf->mean_cluster->Clone();
    itf->ratio->Divide(itf->mean_string);
    itf->ratio->SetLineColor(itf->color);
    itf->ratio->SetFillColor(itf->color);
    itf->ratio->SetMarkerColor(itf->color);
    itf->ratio->SetMarkerStyle(itf->markerstyle);
    itf->ratio->GetXaxis()->SetRangeUser(20, 135);
    itf->ratio->GetYaxis()->SetRangeUser(0.9985, 1.005);
    itf->ratio->SetTitle("Ratio;GenJet p_{T};GenJet response, Cluster/String");
    
    leg->AddEntry(itf->ratio, itf->label, "LP");
    
    // FIT
    itf->ratioFit = new TF1((std::string(itf->label)+std::string("Fit")).c_str(), "[0]+[1]*x+[2]*x**2+[3]*x**3");
    itf->ratioFit->SetLineColor(itf->color);
    itf->ratioFit->SetLineWidth(2);
    itf->ratioFit->SetRange(20, 140);
    
    itf->ratio->Fit((std::string(itf->label)+std::string("Fit")).c_str(), "FEMR");
  }
  
  for (itf = flavors.begin(); itf != flavors.end(); ++itf) {
    if (itf == flavors.begin()) {
      itf->ratio->Draw();
    }
    else itf->ratio->Draw("SAME");
  }
  
  leg->Draw();
  canvas->Print("sherpaGenJetResponse_cluster_string.eps");
  
  // Ratio Flavor/mixture
  TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 500, 500);
  
  for (itf = flavors.begin(); itf != flavors.end(); ++itf) {
    itf->ratioAll = (TH1D*) itf->ratio->Clone();
    itf->ratioAll->Divide(flavors[0].ratio);
    itf->ratioAll->GetXaxis()->SetRangeUser(20, 135);
    itf->ratioAll->GetYaxis()->SetRangeUser(0.998, 1.003);
    itf->ratioAll->SetTitle("Ratio;GenJet p_{T};Rel. GenJet response, Cluster/String");
    
    // FIT
    itf->ratioAllFit = new TF1((std::string(itf->label)+std::string("AllFit")).c_str(), "[0]+[1]*x+[2]*x**2+[3]*x**3");
    itf->ratioAllFit->SetLineColor(itf->color);
    itf->ratioAllFit->SetLineWidth(2);
    itf->ratioAllFit->SetRange(20, 140);
    
    itf->ratioAll->Fit((std::string(itf->label)+std::string("AllFit")).c_str(), "FEMR");
  }
  
  for (itf = flavors.begin(); itf != flavors.end(); ++itf) {
    if (itf == flavors.begin()) {
      itf->ratioAll->Draw();
    }
    else itf->ratioAll->Draw("SAME");
  }
  
  leg->Draw();
  canvas2->Print("sherpaGenJetResponse_flavor_ratio.eps");
}
