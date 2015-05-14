#include <string>
#include <vector>
#include <iostream>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"

#include "tdrstyle.C"

struct sample {
  TFile* file;
  TH1F* hist;
  std::string filename;
  const char* label;
  int color;
  int linestyle;
  sample(const char* f, const char* l, int c = kBlack, int s = 1)
  : filename(f), label(l), color(c), linestyle(s) {}
};

std::vector<sample> samples;
std::vector<sample>::iterator it;

void bFragmentationPlot() {
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  
  samples.push_back(sample("MC_BJES_TuneZ2star.root", "Pythia6, Z2*", kRed+1));
  samples.push_back(sample("MC_BJES_TuneZ2star_rbLEP.root", "Pythia6, Z2*rbLEP", kBlue+1));
  samples.push_back(sample("MC_BJES_TuneZ2star_rbLEPhard.root", "Pythia6, Z2*rbLEP hard", kGreen+1, 7));
  samples.push_back(sample("MC_BJES_TuneZ2star_rbLEPsoft.root", "Pythia6, Z2*rbLEP soft", kMagenta+1, 7));
  samples.push_back(sample("MC_BJES_TuneP12.root", "Pythia6, P12", kBlack));
  
  TLegend *leg = new TLegend(0.25, 0.7, 0.55, 0.9);
  leg->SetTextSize(0.03);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  
  for (it = samples.begin(); it != samples.end(); ++it) {
    it->file = new TFile((std::string("/afs/desy.de/user/m/mseidel/xxl/CMSSW_5_3_11/src/TopAnalysis/TopUtils/data/")+it->filename).c_str());
    it->hist = (TH1F*) it->file->Get("EventWeightBJES/genBHadronPtFraction");
    
    it->hist->Scale(1./it->hist->Integral());
    it->hist->SetLineColor(it->color);
    it->hist->SetLineStyle(it->linestyle);
    
    it->hist->GetXaxis()->SetTitle("p_{T}^{B hadron}/p_{T}^{genJet}");
    it->hist->GetXaxis()->SetRangeUser(0, 1.25);
    it->hist->GetYaxis()->SetRangeUser(0, it->hist->GetMaximum()*1.5);
    
    leg->AddEntry(it->hist, it->label, "L");
  }
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 500, 500);
  for (it = samples.begin(); it != samples.end(); ++it) {
    if (it == samples.begin()) it->hist->Draw();
    else it->hist->Draw("SAME");
  }
  leg->Draw();
  
  canvas->Print("bFragmentationTTJets.eps");
}
