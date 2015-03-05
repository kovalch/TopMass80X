#include <string>
#include <vector>
#include <iostream>

#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TF1.h"
#include "TLegend.h"

#include "tdrstyle.C"

struct sample {
  TFile* file;
  TH1F* w;
  TH1F* t;
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
  TH1F* mean_cluster;
  TH1F* mean_string;
  TH1F* ratio;
  TF1* ratioFit;
  TH1F* ratioAll;
  TF1* ratioAllFit;
  int color;
  int markerstyle;
  flavor(const char* l, TH1F* mc, TH1F* ml, int c = kBlack, int m = 20)
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

void treeToHist(std::string sample, Long64_t nentries = 1000000000) {
  TFile* histFile = new TFile((sample+std::string("_massHist.root")).c_str(), "RECREATE");
  
  TChain* chain = new TChain("analyzeSherpaGenEvent/eventTree");
  chain->Add((sample+std::string("/job*.root")).c_str());
  
  chain->Draw("sherpa.genJetW[0].M() >> w(500,0,250)", "sherpa.genJet[3].Pt()>30 & sherpa.genJet[4].Pt()>30", "", nentries);
  chain->Draw("sherpa.genJetTop[0].M() >> t(500,0,500)", "sherpa.genJet[1].Pt()>30 & sherpa.genJet[3].Pt()>30 & sherpa.genJet[4].Pt()>30 & sherpa.genJetW[0].M()>70 & sherpa.genJetW[0].M()<90", "", nentries);
  
  histFile->Write();
}

void sherpaMass() {
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TH1::SetDefaultSumw2();
  
  // Samples
  samples.push_back(sample("TT_Cluster_8TeV-sherpa_massHist.root", "Cluster", kRed+1));
  samples.push_back(sample("TT_Lund_8TeV-sherpa_massHist.root", "Lund", kBlue+1));
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    it->file = new TFile((std::string("/nfs/dust/cms/user/mseidel/batch_output/analyzeSherpaGenEvent/")+it->filename).c_str());
    
    it->w = (TH1F*) it->file->Get("w");
    it->t = (TH1F*) it->file->Get("t");
    
    it->w->Sumw2();
    it->t->Sumw2();
    
    for (int i = 0; i<100; ++i) {
      it->w->SetBinContent(i, 0);
    }
    for (int i = 500; i<502; ++i) {
      it->w->SetBinContent(i, 0);
    }
    
    for (int i = 0; i<100; ++i) {
      it->t->SetBinContent(i, 0);
    }
    
    it->w->ComputeIntegral();
    it->t->ComputeIntegral();
    
    //it->w->Scale(1./it->w->Integral());
    //it->t->Scale(1./it->t->Integral());
    
    it->w->Rebin(1);
    it->w->GetXaxis()->SetRangeUser(50,300);
    it->w->GetXaxis()->SetTitle("m_{W/t}^{particle} [GeV]");
    
    
    it->t->Rebin(1);
    it->t->GetXaxis()->SetRangeUser(50,300);
    it->t->GetXaxis()->SetTitle("m_{W/t}^{particle} [GeV]");
    
    //it->w->Fit("gaus", "", "", 70, 90);
    //it->t->Fit("gaus", "", "", 155, 190);

    std::cout << "Mean mW: " << it->w->GetMean() << std::endl;    
    std::cout << "Mean mt: " << it->t->GetMean() << std::endl;
    
    it->w->SetLineColor(it->color);
    it->w->SetMarkerColor(it->color);
    
    it->t->SetLineColor(it->color);
    it->t->SetMarkerColor(it->color);
  }
  
  // Flavors
  flavors.push_back(flavor("top", samples[0].t, samples[1].t, kBlue-4, 22));
  flavors.push_back(flavor("W",   samples[0].w, samples[1].w, kRed-4, 20));
  
  // Ratio Cluster/String
  TCanvas* canvas = new TCanvas("canvas", "canvas", 500, 500);
  TLegend *leg = new TLegend(0.6, 0.75, 0.9, 0.9);
  leg->SetTextSize(0.05);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  
  for (itf = flavors.begin(); itf != flavors.end(); ++itf) {
    itf->ratio = (TH1F*) itf->mean_cluster->Clone();
    itf->ratio->Divide(itf->mean_string);
    itf->ratio->SetLineColor(itf->color);
    itf->ratio->SetFillColor(itf->color);
    itf->ratio->SetMarkerColor(itf->color);
    itf->ratio->SetMarkerStyle(itf->markerstyle);
    itf->ratio->GetXaxis()->SetRangeUser(0, 500);
    itf->ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
    itf->ratio->SetTitle("Ratio;Mass [GeV];Cluster/String");
    
    //leg->AddEntry(itf->ratio, itf->label, "LP");
    
    // FIT
    itf->ratioFit = new TF1((std::string(itf->label)+std::string("Fit")).c_str(), "[0]+[1]*(x-172.5)");
    itf->ratioFit->SetLineColor(itf->color);
    itf->ratioFit->SetLineWidth(2);
    itf->ratioFit->SetRange(50, 200);
    
    itf->ratio->Fit((std::string(itf->label)+std::string("Fit")).c_str(), "FEMR");
  }
  
  /*
  for (itf = flavors.begin(); itf != flavors.end(); ++itf) {
    if (itf == flavors.begin()) {
      itf->w->Draw("HIST");
    }
    else itf->w->Draw("SAME,HIST");
  }
  */
  
  samples[1].w->SetLineWidth(2);
  samples[1].w->GetYaxis()->SetRangeUser(0, 100000);
  samples[1].w->Draw("HIST");
  samples[0].w->SetLineWidth(2);
  samples[0].w->Draw("HIST,SAME");
  
  samples[1].t->SetLineWidth(2);
  samples[1].t->Draw("HIST,SAME");
  samples[0].t->SetLineWidth(2);
  samples[0].t->Draw("HIST,SAME");

  
  
  leg->AddEntry(samples[0].t, samples[0].label, "LP");
  leg->AddEntry(samples[1].t, samples[1].label, "LP");
  
  //flavors[0].ratio->Draw();
  //flavors[1].ratio->Draw("SAME");
  
  leg->Draw();
  canvas->Print("sherpaMass_cluster_string.eps");
  
  /*
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
  canvas2->Print("sherpaMass_ratio.eps");
  */
}
