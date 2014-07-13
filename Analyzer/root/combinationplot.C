#include <vector>
#include <iostream>

#include "TCanvas.h"
#include "TH1F.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"

#include "tdrstyle.C"

enum styles        { kStat,  kOverall, kTev };
int color_ [ 4 ] = { kRed+1, kBlack  , kGray};

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

void combinationplot()
{
  TStyle* tdrStyle = setTDRStyle();
  tdrStyle->SetPadLeftMargin(0.05);
  tdrStyle->SetPadRightMargin(0.05);
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.12);
  tdrStyle->SetNdivisions(505, "XYZ");
  
  TCanvas* measurementComparison = new TCanvas("topmass_cms", "topmass_cms", 550, 650);
  measurementComparison->cd();
  
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle(";m_{t} [GeV];");
  
  measurements.push_back(measurement(175.5, 4.6, 4.6, "CMS 2010, dilepton", "JHEP 07 (2011) 049, 36 pb^{-1}"));
  measurements.push_back(measurement(173.10, 2.10, 2.63, "CMS 2010, lepton+jets", "PAS TOP-10-009, 36 pb^{-1}"));
  measurements.push_back(measurement(172.50, 0.43, 1.43, "CMS 2011, dilepton", "EPJC 72 (2012) 2202, 5.0 fb^{-1}"));
  measurements.push_back(measurement(173.49, 0.43, 0.98, "CMS 2011, lepton+jets", "JHEP 12 (2012) 105, 5.0 fb^{-1}"));
  measurements.push_back(measurement(173.49, 0.69, 1.21, "CMS 2011, all-hadronic", "arXiv:1307.4617, 3.5 fb^{-1}"));
  
  measurements.push_back(measurement(172.04, 0.188, 0.75, "CMS 2012, lepton+jets", "PAS TOP-14-001, 19.7 fb^{-1}"));
  measurements.push_back(measurement(172.03, 0.36 , 0.82, "CMS 2012, all+jets", "PAS TOP-14-002, 18.2 fb^{-1}"));
  
  measurements.push_back(measurement(172.21, 0.13, 0.73, "CMS combination", "May 2014", kRed+1));
  
  measurements.push_back(measurement(173.18, 0.56, 0.75, "Tevatron combination", "Phys. Rev. D86 (2012) 092003", kGray+3));
  
  measurements.push_back(measurement(173.34, 0.27, 0.71, "World combination 2014", "ATLAS, CDF, CMS, D0", kGray+3)); //CMS PAS TOP-13-014
  
  int nBox = 8; // set 'reference' value, -1 for last entry
  int n = measurements.size();
  if (nBox == -1) nBox = n-1;
  
  double xBegin = 164.;
  double xEnd   = 182.;
  
  double centerX[] = {(xEnd+xBegin)/2.};
  double centerY[] = {-1};
  double rangeX[]  = {(xEnd-xBegin)/2.};
  double rangeY[]  = {0.0001};
  
  double position[n];
  double mass[n];
  double stat[n];
  //double syst[n];
  double overall[n];
  double posError[n];
  
  for(it = measurements.begin(); it != measurements.end(); ++it) {
    int i = it - measurements.begin();
    
    position[i] = n - i;
    if (position[i] < 3) position[i] = n - i - 0.;
    mass[i] = it->value;
    stat[i] = it->stat;
    //syst[i] = it->syst;
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
    
    char value[6];   sprintf(value  , "%3.1f", it->value);
    char statVal[6]; sprintf(statVal, "%3.1f", it->stat);
    char systVal[6]; sprintf(systVal, "%3.1f", it->syst);
    
    TString numbers = "";
    numbers += value; numbers += " #pm ";
    numbers += statVal; numbers += " #pm ";
    numbers += systVal; numbers += " GeV";
    
    DrawLabel(numbers, 0.45, 0.68,
              tdrStyle->GetPadBottomMargin() + (1. - tdrStyle->GetPadBottomMargin() - tdrStyle->GetPadTopMargin()) * 1./(n*1.1) - 0.04
            + (n - i - 1)
            * (1. - tdrStyle->GetPadBottomMargin() - tdrStyle->GetPadTopMargin()) / (n*1.1),
              0.93, it->color);
    
    DrawLabel("(value #pm stat #pm syst)", 0.4, 0.68,
              tdrStyle->GetPadBottomMargin() + (1. - tdrStyle->GetPadBottomMargin() - tdrStyle->GetPadTopMargin()) * 1./(n*1.1) - 0.065
            + (n - i - 1)
            * (1. - tdrStyle->GetPadBottomMargin() - tdrStyle->GetPadTopMargin()) / (n*1.1),
              0.93, it->color);
  }
  
  /*
  DrawLabel("Methods", 0.4, 0.07, 0.025, 0.5);
  DrawLabel("T: Template, I: Ideogram, M: Matrix Element, J: with JES", 0.4, 0.07, 0.00, 0.5);
  */
  
  TH1F* h0 = new TH1F("h0", "h0", 1, measurements[nBox-1].value - sqrt(pow(measurements[nBox-1].stat,2)+pow(measurements[nBox-1].syst,2)), measurements[nBox-1].value + sqrt(pow(measurements[nBox-1].stat,2)+pow(measurements[nBox-1].syst,2)));
  h0->SetBinContent(1, n+1);
  h0->SetLineColor(kRed-9);
  h0->SetFillColor(kRed-9);
  h0->SetFillStyle(3345);
  h0->Draw("same");
  
  TH1F* h1 = new TH1F("h1", "h1", 1, measurements[nBox].value - sqrt(pow(measurements[nBox].stat,2)+pow(measurements[nBox].syst,2)), measurements[nBox].value + sqrt(pow(measurements[nBox].stat,2)+pow(measurements[nBox].syst,2)));
  h1->SetBinContent(1, n+1);
  h1->SetLineColor(kBlue-9);
  h1->SetFillColor(kBlue-9);
  h1->SetFillStyle(3354);
  //h1->Draw("same");
  
  mg->Draw("AP,same");
  
  double lineY = 2.2;
  TLine *line = new TLine();
  line->SetLineStyle(2);
  line->DrawLine(164., lineY, 182., lineY);
  
  //DrawCMSPrel();
  
  measurementComparison->Print("topmass_cms.eps");
  measurementComparison->Print("topmass_cms.pdf");
}
