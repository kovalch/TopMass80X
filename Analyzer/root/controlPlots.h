#ifndef TOPMASSCONTROLPLOTS_H
#define TOPMASSCONTROLPLOTS_H

#include "TChain.h"
#include "TTreeFormula.h"
#include "TH1F.h"

#include <vector>

enum sampleType {kData, kSig, kBkg, kSigVar};

class TopMassControlPlots{
public:
  TopMassControlPlots();
  ~TopMassControlPlots() {}

private:

  class MySample{
  public:
    MySample(std::string name_, std::string file_, int type_, int color_, int line_ = 1, double scale_ = 1.) :
      name(name_), file(file_), type(type_), color(color_), line(line_), scale(scale_) {}
    
    std::string name, file;
    int type, color, line;
    double scale;
    TChain* chain;
  };

  class MyHistogram{
  public:
    MyHistogram(std::string name_, std::string formula_, std::string selection_, std::string title, int nBins, double min, double max) :
      name(name_), formula(formula_), selection(selection_),
      data(new TH1F((std::string("hD")+name).c_str(), (std::string("Data")+title).c_str(), nBins, min, max))
    {
      data->SetLineWidth(1);
      data->SetLineColor(kBlack);
      data->SetMarkerStyle(20);
      data->SetMarkerColor(kBlack);
    }

    void Init(TChain* chain)
    {
      var = new TTreeFormula((std::string("f")+name).c_str(), formula.c_str(), chain);
      if (selection.size() > 0) sel = new TTreeFormula((std::string("s")+name).c_str(), selection.c_str(), chain);
    }
    void AddSignal(MySample* sample)
    {
      int colorShift[] =  {-11, -8, 0};
      std::string combinationType[] = {" unmatched", " wrong", " correct"};
      for(int i = 0; i < 3; ++i) {
        sig.push_back((TH1F*)data->Clone());
        sig.back()->Reset();
        sig.back()->SetLineWidth(1);
        sig.back()->SetFillColor(sample->color+colorShift[i]);
        sig.back()->SetTitle((sample->name + combinationType[i]).c_str());
      }
    }
    void AddBackground(MySample* sample)
    {
      bkg.push_back((TH1F*)data->Clone());
      bkg.back()->Reset();
      bkg.back()->SetLineWidth(1);
      bkg.back()->SetFillColor(sample->color);
      bkg.back()->SetTitle(sample->name.c_str());
    }
    void AddSignalVariation(MySample* sample)
    {
      sigvar.push_back((TH1F*)data->Clone());
      sigvar.back()->Reset();
      sigvar.back()->SetLineWidth(2);
      sigvar.back()->SetLineColor(sample->color);
      sigvar.back()->SetLineStyle(sample->line);
      sigvar.back()->SetTitle(sample->name.c_str());
    }
    
    TTreeFormula* var;
    TTreeFormula* sel;
    std::vector<TH1F*> sig, bkg, sigvar;
    std::string name, formula, selection;
    TH1F *data;
  };

  void doPlots();

  std::vector<MyHistogram> hists;
  std::vector<MySample>    samples;
};

#endif /* TOPMASSCONTROLPLOTS_H */
