#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLatex.h"
#include "TBox.h"
#include "TPad.h"
#include "TPaveLabel.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"

#include "tdrstyle.C"

enum styles             { kStat,  kOverall, kTev};
int color_      [ 4 ] = { kRed+1, kBlack  , kGray};

struct measurement {
  double value;
  double stat;
  double syst;
  const char* label1;
  const char* label2;
  int color;
  measurement(double v, double st, double sy, const char* l1, const char* l2, int c = kBlack)
  : value(v), stat(st), syst(sy), label1(l1), label2(l2), color(c) {}
};


std::vector<measurement> measurements;
std::vector<measurement>::iterator it;


void DrawLabel(TString text, double size, const double x1, const double y1, const double x2, Color_t color = kBlack)
{
  // function to directly draw a label into the active canvas
  double y2 = y1 + 0.05;
  TPaveLabel *label = new TPaveLabel(x1, y1, x2, y2, text, "NDC");
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  label->SetTextSize(size);
  label->SetTextAlign(12);
  label->SetTextColor(color);
  label->Draw("same");
}


void CRComparison()
{
  TStyle* tdrStyle = setTDRStyle();
  tdrStyle->SetPadLeftMargin(0.05);
  tdrStyle->SetPadRightMargin(0.05);
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.12);
  tdrStyle->SetNdivisions(505, "XYZ");
  
  TCanvas* measurementComparison = new TCanvas("measurementComparison", "measurementComparison", 550, 650);
  measurementComparison->cd();
  
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle(";m_{t} [GeV];");
  
  measurements.push_back(measurement(172.856, 0.277397, 0., "Tune P0, no CR, 4M", "PARP(77)=0, PARP(78)=0"));
  measurements.push_back(measurement(172.633, 0.278354, 0., "Tune P0, CR, 4M", ""));
  measurements.push_back(measurement(172.693, 0.270812, 0., "Tune Z2, no CR, 4M (1)", "MSTP(95)=0"));
  measurements.push_back(measurement(172.575, 0.273624, 0., "Tune Z2, no CR, 4M (2)", "PARP(77)=0, PARP(78)=0"));
  measurements.push_back(measurement(172.633, 0.195448, 0., "Tune Z2, no CR, 8M", "Combined (1) and (2))"));
  measurements.push_back(measurement(172.47, 0.193422, 0., "Tune Z2, CR, 8M", ""));
  
  int nBox = -1; // set 'reference' value, -1 for last entry
  int n = measurements.size();
  if (nBox == -1) nBox = n-1;
  
  double xBegin = 171.5;
  double xEnd   = 173.5;
  
  double centerX[] = {(xEnd+xBegin)/2.};
  double centerY[] = {-1};
  double rangeX[]  = {(xEnd-xBegin)/2.};
  double rangeY[]  = {0.0001};
  
  double position[n];
  double mass[n];
  double stat[n];
  double syst[n];
  double overall[n];
  double posError[n];
  
  
  
  for(it = measurements.begin(); it != measurements.end(); ++it) {
    int i = it - measurements.begin();
    
    position[i] = n - i;
    mass[i] = it->value;
    stat[i] = it->stat;
    syst[i] = it->syst;
    overall[i]  = sqrt(pow(it->stat,2)+pow(it->syst,2));
    posError[i] = 0.0001;
  }
  
  TGraphErrors* gStat = new TGraphErrors(n, mass, position, stat, posError);
  gStat->SetMarkerColor(color_[kStat]);
  gStat->SetLineColor(color_[kStat]);
  gStat->SetMarkerStyle(20);
  gStat->SetLineWidth(4);
  gStat->SetMarkerSize(1.5);
  
  TGraphErrors* gOverall = new TGraphErrors(n, mass, position, overall, posError);
  gOverall->SetMarkerColor(color_[kStat]);
  gOverall->SetLineColor(color_[kStat]);
  gOverall->SetMarkerStyle(20);
  gOverall->SetLineWidth(2);
  gOverall->SetMarkerSize(1.5);
  
  TGraphErrors* gRange = new TGraphErrors(1, centerX, centerY, rangeX, rangeY);
  gRange->SetLineColor(kWhite);
  gRange->SetLineWidth(0);
  gRange->SetMarkerStyle(20);
  gRange->SetMarkerColor(kWhite);
  
  mg->Add(gRange);
  mg->Add(gOverall);
  mg->Add(gStat);
  
  mg->Draw("AP");
  
  mg->SetMinimum(0);
  mg->SetMaximum(n * 1.1);
  mg->GetYaxis()->SetLabelColor(kWhite);
  mg->GetYaxis()->SetTickLength(0.);
  
  for(it = measurements.begin(); it != measurements.end(); ++it) {
    int i = it - measurements.begin();
    
    DrawLabel(it->label1, 0.45, 0.07,
              tdrStyle->GetPadBottomMargin() + (1. - tdrStyle->GetPadBottomMargin() - tdrStyle->GetPadTopMargin()) * 1./(n*1.1) - 0.04
            + (n - i - 1)
            * (1. - tdrStyle->GetPadBottomMargin() - tdrStyle->GetPadTopMargin()) / (n*1.1),
              0.5, it->color);
              
    DrawLabel(it->label2, 0.4, 0.07,
              tdrStyle->GetPadBottomMargin() + (1. - tdrStyle->GetPadBottomMargin() - tdrStyle->GetPadTopMargin()) * 1./(n*1.1) - 0.065
            + (n - i - 1)
            * (1. - tdrStyle->GetPadBottomMargin() - tdrStyle->GetPadTopMargin()) / (n*1.1),
              0.5, it->color);
    
    char value[6]; sprintf(value, "%3.2f", it->value);
    char stat[6]; sprintf(stat, "%3.2f", it->stat);
    char syst[6]; sprintf(syst, "%3.2f", it->syst);
    
    TString numbers = "";
    numbers += value; numbers += " #pm ";
    numbers += stat;// numbers += " #pm ";
    //numbers += syst;
    numbers += " GeV";
    
    DrawLabel(numbers, 0.45, 0.68,
              tdrStyle->GetPadBottomMargin() + (1. - tdrStyle->GetPadBottomMargin() - tdrStyle->GetPadTopMargin()) * 1./(n*1.1) - 0.04
            + (n - i - 1)
            * (1. - tdrStyle->GetPadBottomMargin() - tdrStyle->GetPadTopMargin()) / (n*1.1),
              0.93, it->color);
    
    DrawLabel("(value #pm stat)", 0.4, 0.68,
              tdrStyle->GetPadBottomMargin() + (1. - tdrStyle->GetPadBottomMargin() - tdrStyle->GetPadTopMargin()) * 1./(n*1.1) - 0.065
            + (n - i - 1)
            * (1. - tdrStyle->GetPadBottomMargin() - tdrStyle->GetPadTopMargin()) / (n*1.1),
              0.93, it->color);
  }
  
  DrawLabel("FASTSIM, PYTHIA6 only", 0.4, 0.07, 0.0, 0.5);
  
  TH1F* h1 = new TH1F("h1", "h1", 1, measurements[nBox].value - sqrt(pow(measurements[nBox].stat,2)+pow(measurements[nBox].syst,2)), measurements[nBox].value + sqrt(pow(measurements[nBox].stat,2)+pow(measurements[nBox].syst,2)));
  h1->SetBinContent(1, n+1);
  h1->SetLineColor(kGray);
  h1->SetFillColor(kGray);
  h1->Draw("same");
  
  mg->Draw("AP,same");
}
