#include <iostream>
#include <list>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TList.h"

#include "RooAbsData.h"
#include "RooAbsPdf.h"
//#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

int color_ [] = { kBlue+1, kGreen+1, kRed+1, kBlack };

std::vector<TH1F*> hists(4,NULL);

void modifyAndSaveCanvas(TCanvas* canvas)
{
  TList* list = canvas->GetListOfPrimitives();
  for(int i = 0, l = list->GetSize(); i < l; ++i){
    TObject* obj = list->At(i);
    //std::cout << obj->ClassName() << ": "<< obj->GetName() << std::endl;
    if(obj->InheritsFrom("TH1")){
      TH1* hist = (TH1*)obj;
      hist->SetTitle("");
      hist->GetYaxis()->SetLabelSize(0);
      hist->GetYaxis()->SetTitle("Probability Density (a.u.)");
      hist->GetYaxis()->SetRangeUser(0,hist->GetMaximum()*1.1);
    }
  }
  canvas->Print((std::string(canvas->GetName())+std::string(".eps")).c_str());
}

void plotter(std::string prefix, std::vector<std::string> samples, std::vector<std::string> pdfsTMP)
{
  bool isBKG = false;
  if(samples.at(0).find("BKG") != std::string::npos) isBKG = true;
  TFile* file = 0;
  //if(isBKG)
  //  file = TFile::Open("Calibration_AllJets_bkg.root", "READ");
  //else
  //  file = TFile::Open("Calibration_AllJets_sig.root", "READ");
  file = TFile::Open("Calibration_AllJets_new.root", "READ");

  RooWorkspace* ws = (RooWorkspace*)file->Get("workspaceMtop");

  //RooArgSet set1 = ws->allPdfs();
  //set1.dump();
  //
  //RooArgSet set2 = ws->allVars();
  //set2.dump();

  //std::list<RooAbsData*> set3 = ws->allData();
  //for(const auto& dat : set3)
  //  std::cout << dat->GetName() << std::endl;

  TCanvas* canvas = new TCanvas("canvas", "canvas", 10, 10, 600, 600);
  canvas->cd();
  for(unsigned int j = 0; j < samples.size(); ++j){
    RooPlot* frame = 0;
    if(prefix.substr(0,4) == "mTop"){
      RooRealVar topMass = RooRealVar("topMass","m_{t}^{fit}",100.,550.,"GeV");
      if(isBKG) frame = topMass.frame(RooFit::Range(100., 550.));
      else{
        if     (j == 0) frame = topMass.frame(RooFit::Range(120., 270.));
        else if(j == 1) frame = topMass.frame(RooFit::Range(100., 350.));
      }
    }
    else if(prefix.substr(0,2) == "mW"){
      RooRealVar topMass = RooRealVar("meanWMass","m_{W}^{reco}",50.,300.,"GeV");
      if(isBKG) frame = topMass.frame(RooFit::Range(60., 150.));
      else      frame = topMass.frame(RooFit::Range(60., 130.));
    }
    TLegend* leg = new TLegend(0.5,0.7,0.89,0.89);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    std::vector<std::string> pdfs;
    for(const auto& pdf : pdfsTMP){
      pdfs.push_back(samples[j]+pdf);
      if(isBKG)
        leg->AddEntry(hists[3], "Background", "L");
      else{
        if     (pdf == "_jes100mass1665"){
          leg->AddEntry(hists[0], "m_{t} = 166.5 GeV", "L");
          leg->AddEntry(hists[1], "m_{t} = 172.5 GeV", "L");
          leg->AddEntry(hists[2], "m_{t} = 178.5 GeV", "L");
        }
        else if(pdf == "_jes096mass1725"){
          leg->SetX1(0.6);
          leg->AddEntry(hists[0], "JES=0.96", "L");
          leg->AddEntry(hists[1], "JES=1.00", "L");
          leg->AddEntry(hists[2], "JES=1.04", "L");
        }
      }
    }
    for(unsigned int i = 0; i < pdfs.size(); ++i){
      if(isBKG){
        ws->data("BKG")->plotOn(frame, RooFit::MarkerColor(color_[3]));
        ws->pdf(pdfs[i].c_str())->plotOn(frame, RooFit::LineColor(color_[3]));
      }
      else{
        std::string datasetName = "dataset_"; datasetName += samples[j].substr(4,1);
        std::cout << datasetName << std::endl;
        if     (pdfs[i].find("_jes100mass1725") != std::string::npos) ws->data((datasetName+"_17").c_str())->plotOn(frame, RooFit::MarkerColor(color_[i]),RooFit::Rescale(5.65));
        else if(pdfs[i].find("_jes100mass1665") != std::string::npos) ws->data((datasetName+"_2" ).c_str())->plotOn(frame, RooFit::MarkerColor(color_[i]));
        else if(pdfs[i].find("_jes100mass1785") != std::string::npos) ws->data((datasetName+"_32").c_str())->plotOn(frame, RooFit::MarkerColor(color_[i]),RooFit::Rescale(1.05));
        else if(pdfs[i].find("_jes096mass1725") != std::string::npos) ws->data((datasetName+"_15").c_str())->plotOn(frame, RooFit::MarkerColor(color_[i]),RooFit::Rescale(5.65));
        else if(pdfs[i].find("_jes104mass1725") != std::string::npos) ws->data((datasetName+"_19").c_str())->plotOn(frame, RooFit::MarkerColor(color_[i]),RooFit::Rescale(5.85));
        ws->pdf(pdfs[i].c_str())->plotOn(frame, RooFit::LineColor(color_[i]));
      }
    }
    frame->Draw();
    leg->Draw("same");
    canvas->SetName((prefix+"_"+samples[j]).c_str());
    modifyAndSaveCanvas(canvas);
  }
}

void function2()
{
  TLegend* leg = new TLegend(0.5,0.7,0.89,0.89);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hists[3], "Data", "L");
  leg->AddEntry(hists[2], "Variation", "L");

  TFile* file1 = TFile::Open("Calibration_AllJets_new.root", "READ");
  TFile* file2 = TFile::Open("Calibration_AllJets_bkg_from_sig.root", "READ");
  RooWorkspace* ws1 = (RooWorkspace*)file1->Get("workspaceMtop");
  RooWorkspace* ws2 = (RooWorkspace*)file2->Get("workspaceMtop");

  TCanvas* canvas1 = new TCanvas("canvas1", "canvas1", 10, 10, 600, 600);
  canvas1->cd();
  RooRealVar topMass = RooRealVar("topMass","m_{t}^{fit}",100.,550.,"GeV");
  RooPlot* frame = topMass.frame(RooFit::Range(100., 550.));
  ws1->pdf("topBKG")->plotOn(frame, RooFit::LineColor(kBlack));
  ws2->pdf("topBKG")->plotOn(frame, RooFit::LineColor(kRed));
  frame->Draw();
  leg->Draw();
  canvas1->SetName("mTop_shape_uncertainty");
  modifyAndSaveCanvas(canvas1);

  TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 660, 10, 600, 600);
  canvas2->cd();
  RooRealVar wMass = RooRealVar("meanWMass","m_{W}^{reco}",50.,300.,"GeV");
  RooPlot* frame2 = wMass.frame(RooFit::Range(60., 130.));
  ws1->pdf("wBKG")->plotOn(frame2, RooFit::LineColor(kBlack));
  ws2->pdf("wBKG")->plotOn(frame2, RooFit::LineColor(kRed));
  frame2->Draw();
  leg->Draw();
  canvas2->SetName("mW_shape_uncertainty");
  modifyAndSaveCanvas(canvas2);

}

void PlotsFromRooWorkspace()
{
  hists[0] = new TH1F();
  hists[0]->SetLineWidth(2);
  for(unsigned int i = 0, l = hists.size(); i < l; ++i){
    if(i > 0) hists[i] = (TH1F*)hists[0]->Clone();
    hists[i]->SetLineColor(color_[i]);
  }

  //plotter("mTop1", {"sig_0", "sig_1"}, {"_jes100mass1665", "_jes100mass1725", "_jes100mass1785"});
  //plotter("mTop2", {"sig_0", "sig_1"}, {"_jes096mass1725", "_jes100mass1725", "_jes104mass1725"});
  //plotter("mW1"  , {"sig_2", "sig_3"}, {"_jes100mass1665", "_jes100mass1725", "_jes100mass1785"});
  //plotter("mW2"  , {"sig_2", "sig_3"}, {"_jes096mass1725", "_jes100mass1725", "_jes104mass1725"});
  //plotter("mTop", {"topBKG"}, {""});
  //plotter("mW"  , {"wBKG"  }, {""});
  function2();
}

