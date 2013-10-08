#ifndef TOPMASSCONTROLPLOTS_H
#define TOPMASSCONTROLPLOTS_H

#include "TChain.h"
#include "TTreeFormula.h"
#include "TH1F.h"

#include <vector>

class TopMassControlPlots{
public:
  TopMassControlPlots();
  ~TopMassControlPlots() {}

private:

  class MyHistogram{
  public:
    MyHistogram(std::string name, std::string formula, TChain* chain, std::string title, int nBins, double min, double max) :
      var(new TTreeFormula((std::string("f")+name).c_str(), formula.c_str(), chain)),
      data(new TH1F((std::string("hD")+name).c_str(), title.c_str(), nBins, min, max))
    {
      data->SetLineWidth(2);
      data->SetLineColor(kBlack);
      data->SetMarkerStyle(20);
      data->SetMarkerColor(kBlack);
    }

    void AddSignal()
    {
      sig.push_back((TH1F*)data->Clone());
      sig.back()->SetLineWidth(2);
      sig.back()->SetLineColor(kRed);
      sig.back()->SetLineStyle(9);
    }
    void AddBackground()
    {
      bkg.push_back((TH1F*)data->Clone());
      bkg.back()->SetLineWidth(2);
      bkg.back()->SetLineColor(kBlue);
      bkg.back()->SetLineStyle(2);
    }
    TTreeFormula* var;
    TH1F *data;
    std::vector<TH1F*> sig, bkg;
  };

  void doPlots();

  std::vector<MyHistogram> hists;
};

#endif /* TOPMASSCONTROLPLOTS_H */
