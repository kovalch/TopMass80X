#include <vector>
#include <iostream>
#include <iomanip>

#include "TCanvas.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TLine.h"

#include "tdrstyle.C"

struct uncertainty {
  const char* name;
  double syst;
  double stat;
  uncertainty(const char* n, double _syst, double _stat)
  : name(n), syst(_syst), stat(_stat) {}
};

void drawcutline(double cutval, double maximum)
{
  TLine *cut = new TLine();
  cut->SetLineWidth(2);
  cut->SetLineStyle(7);
  cut->SetLineColor(kBlue+1);
  cut->DrawLine(cutval, 0.,cutval, maximum);
}

void toySystematics() {
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(1);
  gStyle->SetErrorX(0.5);

  std::vector<uncertainty> uncertainties;
  std::vector<uncertainty>::iterator it;

  double correlated = 0.;

  uncertainties.push_back(uncertainty("Calibration", 0.06, correlated));
  uncertainties.push_back(uncertainty("PDF", 0.07, correlated));
  uncertainties.push_back(uncertainty("UE", 0.153, 0.187));
  uncertainties.push_back(uncertainty("Background", 0.126, correlated));
  uncertainties.push_back(uncertainty("bJES", 0.61, correlated));
  uncertainties.push_back(uncertainty("JES", 0.28, correlated));
  uncertainties.push_back(uncertainty("JER", 0.22, correlated));
  uncertainties.push_back(uncertainty("Matching", 0.17, 0.26));
  uncertainties.push_back(uncertainty("Scale", 0.24, 0.27));
  uncertainties.push_back(uncertainty("CR", 0.54, 0.25));
  uncertainties.push_back(uncertainty("PU", 0.07, correlated));
  uncertainties.push_back(uncertainty("b-tag", 0.12, correlated));
  uncertainties.push_back(uncertainty("LES", 0.02, correlated));
  uncertainties.push_back(uncertainty("MET", 0.06, correlated));
  
  TH1F* hTotal = new TH1F("hTotal", ";Total systematic uncertainty;Number of pseudo-experiments", 100, 0, 2);
  TRandom3* random = new TRandom3(0);
  double trueUncertainty2;
  
  for (int i = 0; i < 1000000; i++) {
    double totalUncertainty2 = 0;
    trueUncertainty2 = 0;
    for (it = uncertainties.begin(); it != uncertainties.end(); ++it) {
      totalUncertainty2 += pow(random->Gaus(it->syst, it->stat), 2);
      trueUncertainty2  += pow(it->syst, 2);
    }
    hTotal->Fill(sqrt(totalUncertainty2));
  }

  TCanvas* canvas = new TCanvas("canvas", "canvas", 500, 500);
  canvas->cd();
  
  hTotal->Draw();
  drawcutline(sqrt(trueUncertainty2), hTotal->GetMaximum());
};
