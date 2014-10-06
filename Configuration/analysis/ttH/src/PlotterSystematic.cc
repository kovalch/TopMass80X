#include <fstream>
#include <iostream>
#include <cstdio>
#include <sstream>
#include <cmath>
#include <iomanip>

#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TExec.h>
#include <TStyle.h>
#include <TMath.h>
#include <TROOT.h>
#include <THStack.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGaxis.h>
#include <TPaveText.h>
#include <TClass.h>
#include <TError.h>

#include "PlotterSystematic.h"
#include "higgsUtils.h"
#include "Samples.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"
#include "../../common/include/plotterUtils.h"





PlotterSystematic::PlotterSystematic(const char* outputDir,
                                     const std::map<Channel::Channel, std::map<Systematic::Systematic, std::map<TString, std::pair<TString, TString> > > >& inputFileLists):
outputDir_(outputDir),
inputFileLists_(inputFileLists),
fileReader_(RootFileReader::getInstance()),
name_("defaultName"),
rangemin_(0),
rangemax_(3),
ymin_(0),
ymax_(0),
YAxis_(""),
XAxis_(""),
logX_(false),
logY_(false)
{
    // Suppress default info that canvas is printed
    gErrorIgnoreLevel = 1001;
}



void PlotterSystematic::setOptions(const TString& name, const TString&,
                         const TString& YAxis, const TString& XAxis,
                         const int , const bool ,
                         const bool logX, const bool logY,
                         const double& ymin, const double& ymax,
                         const double& rangemin, const double& rangemax)
{
    name_ = name; //Histogram name to be plotted
    YAxis_ = YAxis; //Y-axis title
    XAxis_ = XAxis; //X-axis title
    logX_ = logX; //Draw X-axis in Log scale
    logY_ = logY; //Draw Y-axis in Log scale
    ymin_ = ymin; //Min. value in Y-axis
    ymax_ = ymax; //Max. value in Y-axis
    rangemin_ = rangemin; //Min. value in X-axis
    rangemax_ = rangemax; //Max. value in X-axis
    
    printf("Running plot for: %s\n", name.Data());
    
    processNames_.push_back("ttbb");
    processNames_.push_back("tt2b");
    processNames_.push_back("ttOther");
}



void PlotterSystematic::producePlots()
{
    //std::cout<<"--- Beginning of plot production\n\n";
    
    prepareStyle();
    
    for(TString processName : processNames_) {
        for(auto channelCollection : inputFileLists_) {
            Channel::Channel channel = channelCollection.first;
            SystematicHistoMap systematicHistos;
            // Getting up/down variation histograms for each systematic
            for(auto systematicCollection : channelCollection.second) {
                Systematic::Systematic systematic = systematicCollection.first;
                if(systematicCollection.second.count(name_) < 1) {
                    printf("### Warning! No file (%s) for %s : %s\n\n", name_.Data(), Systematic::convertType(systematic.type()).Data(), Channel::convert(channel).Data());
                    continue;
                }
                std::pair<TString, TString> upDownFile = systematicCollection.second.at(name_);
                TH1* histoUp = fileReader_->GetClone<TH1>(upDownFile.first, name_+"_"+processName, true, false);
                if(!histoUp) {
                    printf("### Warning! No histogram for process (%s) in file: %s\n\n", processName.Data(), upDownFile.first.Data());
                }
                TH1* histoDown = fileReader_->GetClone<TH1>(upDownFile.second, name_+"_"+processName, true, false);
                if(!histoDown) {
                    printf("### Warning! No histogram for process (%s) in file: %s\n\n", processName.Data(), upDownFile.second.Data());
                }
                HistoPair pair(histoUp, histoDown);
                systematicHistos[systematic.type()] = pair;
            }
            writeVariations(systematicHistos, channel, processName);
            if(logY_) writeVariations(systematicHistos, channel, processName, true);
        }
    }
    
    //std::cout<<"\n=== Finishing of plot production\n\n";
}


void PlotterSystematic::writeVariations(const SystematicHistoMap& histoCollection, const Channel::Channel channel, const TString processName, const bool logY)
{
    // Prepare canvas and legend
    TCanvas* canvas = new TCanvas("","");
    TLegend* legend = new TLegend(0.67,0.55,0.92,0.85);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.25);
    legend->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 - legend->GetNRows()*0.05);
    legend->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength());
    legend->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength());
    legend->Clear();
    canvas->Clear();
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    canvas->SetName("");
    canvas->SetTitle("");
    if(logY) canvas->SetLogy();
    
    // Getting nominal histogram
    TH1* h_nominal = histoCollection.at(Systematic::nominal).first;
    if(h_nominal) {
        updateHistoAxis(h_nominal, logY);
        setHistoStyle(h_nominal, 1, 1, 3, 0, 0, 20, 1, 0.75);
        h_nominal->Draw("HISTE");
    } else {
        printf("No Nominal histogram for sample: ");
    }
    
    drawRatioPad(canvas, 0.75, 1.25, h_nominal);
    TVirtualPad* ratioPad = gPad;
    
    // Looping over all available systematics
    size_t orderId = 0;
    for(auto systematicHistos : histoCollection) {
        canvas->cd(0);
        const Systematic::Type systematicType(systematicHistos.first);
        if(systematicHistos.first == Systematic::nominal) continue;
        TH1* h_up = systematicHistos.second.first;
//         std::cout << "UP: " << h_up << std::endl;
        if(h_up) {
            updateHistoAxis(h_up, logY);
            setHistoStyle(h_up, lineStyles_.first, lineColors_.at(orderId), 2);
            h_up->Draw("sameHIST");
        }
        TH1* h_down = systematicHistos.second.second;
//         std::cout << "DOWN: " << h_down << std::endl;
        if(h_down) {
            updateHistoAxis(h_down, logY);
            setHistoStyle(h_down, lineStyles_.second, lineColors_.at(orderId), 2);
            h_down->Draw("sameHIST");
        }
        
        legend->AddEntry(h_up, Systematic::convertType(systematicType), "l");
        ratioPad->cd();
        TH1* h_ratio_up = ratioHistogram(h_up, h_nominal, 0);
        h_ratio_up->Draw("sameHIST");
        TH1* h_ratio_down = ratioHistogram(h_down, h_nominal, 0);
        h_ratio_down->Draw("sameHIST");

        if(orderId<lineColors_.size()-1) orderId++;
    }
    // Drawing statistical error bars of the nominal histogram
    ratioPad->cd();
    TH1* h_nominal_stat = ratioHistogram(h_nominal, h_nominal, 1);
    setHistoStyle(h_nominal_stat, 1, 1, 3, 3004, 1, 20, 1, 0.75);
    h_nominal_stat->Draw("sameE1");
    canvas->cd(0);
    h_nominal->Draw("sameE1");
    legend->Draw();
    
    this->drawCmsLabels(2, 8);
    this->drawDecayChannelLabel(channel);
    
    TString eventFileString = common::assignFolder(outputDir_, channel, Systematic::Systematic("Nominal"))+name_+"_"+processName+"_systematics";
    if(logY) eventFileString.Append("_logY");
    canvas->Print(eventFileString+".eps");
    
    delete canvas;
}


void PlotterSystematic::prepareStyle()
{
    // Set style
    common::setHHStyle(*gStyle);

    // Different line styles for up/down variations
    lineStyles_.first = 7;
    lineStyles_.second = 4;
    
    // Different colors for up to 11 systematics
    lineColors_.push_back(2);
    lineColors_.push_back(kAzure+2);
    lineColors_.push_back(kTeal+4);
    lineColors_.push_back(kOrange-3);
    lineColors_.push_back(kViolet-4);
    lineColors_.push_back(kOrange+2);
    lineColors_.push_back(kSpring-9);
    lineColors_.push_back(kPink+1);
    lineColors_.push_back(kAzure+10);
    lineColors_.push_back(17);
    lineColors_.push_back(29);
    
}


void PlotterSystematic::setHistoStyle(TH1* hist, Style_t line, Color_t lineColor, Size_t lineWidth, 
                                  Style_t fill, Color_t fillColor, 
                                  Style_t marker, Color_t markerColor, Size_t markerSize)const
{
    hist->SetLineStyle(line);
    hist->SetLineColor(lineColor);
    hist->SetLineWidth(lineWidth);
    
    hist->SetFillStyle(fill);
    hist->SetFillColor(fillColor);
    
    hist->SetMarkerStyle(marker);
    hist->SetMarkerColor(markerColor);
    hist->SetMarkerSize(markerSize);
    
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleOffset(1.7);
    hist->GetXaxis()->SetTitleOffset(1.25);
}


void PlotterSystematic::updateHistoAxis(TH1* histo, const bool logY)const
{
    // Set x and y axis ranges
    if(logY){
      // Setting minimum to >0 value
      // FIXME: Should we automatically calculate minimum value instead of the fixed value?
      histo->SetMinimum(1e-1);
      if(ymin_>0) histo->SetMinimum(ymin_);
    }
    else histo->SetMinimum(ymin_);

    if(rangemin_!=0. || rangemax_!=0.) {histo->SetAxisRange(rangemin_, rangemax_, "X");}
    
    if(ymax_==0.){
        // Determining the highest Y value that is plotted
        float yMax = histo ? histo->GetBinContent(histo->GetMaximumBin()) + histo->GetBinError(histo->GetMaximumBin()) : ymax_;
        
        // Scaling the Y axis
        if(logY){histo->SetMaximum(18.*yMax);}
        else{histo->SetMaximum(1.35*yMax);}
    }
    else{histo->SetMaximum(ymax_);}

    histo->GetXaxis()->SetNoExponent(kTRUE);
    // Set axis titles
    if(XAxis_ != "-") histo->GetYaxis()->SetTitle(YAxis_);
    if(YAxis_ != "-") histo->GetXaxis()->SetTitle(XAxis_);
}


void PlotterSystematic::drawDecayChannelLabel(const Channel::Channel& channel, const double& textSize)const
{
    TPaveText* decayChannel = new TPaveText();

    decayChannel->AddText(Channel::label(channel));

    decayChannel->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
    decayChannel->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );
    decayChannel->SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
    decayChannel->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength()        );

    decayChannel->SetFillStyle(0);
    decayChannel->SetBorderSize(0);
    if (textSize!=0) decayChannel->SetTextSize(textSize);
    decayChannel->SetTextAlign(12);
    decayChannel->Draw("same");
}


void PlotterSystematic::drawCmsLabels(const int cmsprelim, const double& energy, const double& textSize)const
{
    const char* text;
    if(cmsprelim == 2) text = "Private Work, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV"; // Private work for PhDs students
    else if (cmsprelim == 1) text = "CMS Preliminary, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV"; // CMS preliminary label
    else text = "CMS, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV"; // CMS label

    TPaveText* label = new TPaveText();
    label->SetX1NDC(gStyle->GetPadLeftMargin());
    label->SetY1NDC(1.0 - gStyle->GetPadTopMargin());
    label->SetX2NDC(1.0 - gStyle->GetPadRightMargin());
    label->SetY2NDC(1.0);
    label->SetTextFont(42);
//     label->AddText(Form(text, samples_.luminosityInInversePb()/1000., energy));
    label->AddText(Form(text, 19.7, energy));
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    if(textSize != 0) label->SetTextSize(textSize);
    label->SetTextAlign(32);
    label->Draw("same");
}


void PlotterSystematic::setGraphStyle( TGraph* graph, Style_t marker, Color_t markerColor, Size_t markerSize, 
                             Style_t line, Color_t lineColor, Size_t lineWidth)const 
{
    graph->SetMarkerStyle(marker);
    graph->SetMarkerColor(markerColor);
    graph->SetMarkerSize(markerSize);
    graph->SetLineStyle(line);
    graph->SetLineColor(lineColor);
    graph->SetLineWidth(lineWidth);
}


void PlotterSystematic::normalize(TH1* histo)const
{
    histo->Scale(1./histo->Integral(1,-1));
}


TH1* PlotterSystematic::drawRatioPad(TPad* pad, const double yMin, const double yMax, TH1* axisHisto, 
                                 const double fraction, const TString title)const
{
    // y:x size ratio for canvas
    double canvAsym = (pad->GetY2() - pad->GetY1())/(pad->GetX2()-pad->GetX1());
    Double_t left  = pad->GetLeftMargin();
    Double_t right = pad->GetRightMargin();
    // change old pad
    pad->SetBottomMargin(fraction);
    pad->SetRightMargin(right);
    pad->SetLeftMargin(left);
    pad->SetBorderMode(0);
    pad->SetBorderSize(0);
    pad->SetFillColor(10);
    // create new pad for ratio plot
    TPad *rPad;
    rPad = new TPad("rPad","",0,0,1,fraction+0.001);
#ifdef DILEPTON_MACRO
    rPad->SetFillColor(10);
#else
    rPad->SetFillStyle(0);
    rPad->SetFillColor(0);
#endif
    rPad->SetBorderSize(0);
    rPad->SetBorderMode(0);
    rPad->Draw();
    rPad->cd();
    rPad->SetLogy(0);
    rPad->SetTicky(1);
    // configure ratio plot
    double scaleFactor = 1./(canvAsym*fraction);
    TH1* ratio = (TH1*)axisHisto->Clone("axis_ratio");
    ratio->SetStats(kFALSE);
    ratio->SetTitle("");
    ratio->SetName("axis_ratio");
    ratio->SetMaximum(yMax);
    ratio->SetMinimum(yMin);
    // configure axis of ratio plot
    ratio->GetXaxis()->SetTitleSize(axisHisto->GetXaxis()->GetTitleSize()*scaleFactor*1.3);
    ratio->GetXaxis()->SetTitleOffset(axisHisto->GetXaxis()->GetTitleOffset()*0.9);
    ratio->GetXaxis()->SetLabelSize(axisHisto->GetXaxis()->GetLabelSize()*scaleFactor*1.4);
    ratio->GetXaxis()->SetTitle(axisHisto->GetXaxis()->GetTitle());
    ratio->GetXaxis()->SetNdivisions(axisHisto->GetNdivisions());
    ratio->GetYaxis()->CenterTitle();
    ratio->GetYaxis()->SetTitle(title);
    ratio->GetYaxis()->SetTitleSize(axisHisto->GetYaxis()->GetTitleSize()*scaleFactor);
    ratio->GetYaxis()->SetTitleOffset(axisHisto->GetYaxis()->GetTitleOffset()/scaleFactor);
    ratio->GetYaxis()->SetLabelSize(axisHisto->GetYaxis()->GetLabelSize()*scaleFactor);
    ratio->GetYaxis()->SetLabelOffset(axisHisto->GetYaxis()->GetLabelOffset()*3.3);
    ratio->GetYaxis()->SetTickLength(0.03);
    ratio->GetYaxis()->SetNdivisions(405);
    ratio->GetXaxis()->SetRange(axisHisto->GetXaxis()->GetFirst(), axisHisto->GetXaxis()->GetLast());
    // delete axis of initial plot
    axisHisto->GetXaxis()->SetLabelSize(0);
    axisHisto->GetXaxis()->SetTitleSize(0);
    //This is frustrating and stupid but apparently necessary...
    TExec *setex1 { new TExec("setex1","gStyle->SetErrorX(0.5)") };
    setex1->Draw();
    TExec *setex2 { new TExec("setex2","gStyle->SetErrorX(0.)") };
    setex2->Draw();
//     ratio->SetMarkerSize(0.8);
//     ratio->SetLineWidth(2);
    ratio->Draw("axis");
    rPad->SetTopMargin(0.0);
    rPad->SetBottomMargin(0.15*scaleFactor);
    rPad->SetRightMargin(right);
    pad->SetLeftMargin(left);
    pad->RedrawAxis();
    // draw grid
    rPad->SetGrid(0,1);
    rPad->cd();
    
    // draw a horizontal line on the pad
    Double_t xmin = ratio->GetXaxis()->GetXmin();
    Double_t xmax = ratio->GetXaxis()->GetXmax();
    TString height = ""; height += 1;
    TF1 *f = new TF1("f", height, xmin, xmax);
    f->SetLineStyle(1);
    f->SetLineWidth(1);
    f->SetLineColor(kBlack);
    f->Draw("L same");
    
    axisHisto->GetXaxis()->SetLabelSize(0);
    axisHisto->GetXaxis()->SetTitleSize(0);
    
    
    return ratio;
    
}


TH1* PlotterSystematic::ratioHistogram(const TH1* h_nominator, const TH1* h_denominator, const int errorId)const
{
    TH1* h_ratio = (TH1*)h_nominator->Clone();
    for(int iBin = 1; iBin<=h_ratio->GetNbinsX(); ++iBin) {
        const ValueError nominator(h_nominator->GetBinContent(iBin), h_nominator->GetBinError(iBin));
        const ValueError denominator(h_denominator->GetBinContent(iBin), h_denominator->GetBinError(iBin));
        ValueError ratio;
        ratio.v = nominator.v / denominator.v;
        if(errorId == 0) ratio.e = 0.;
        else if(errorId == 1) ratio.e = nominator.e/nominator.v;
        else if(errorId == 2) ratio.e = denominator.e/denominator.v;
        else if(errorId == 3) ratio.e = ratio.v * sqrt(nominator.eOv2() + denominator.eOv2());
            
        
        
        if(ratio.v != ratio.v) ratio.v = 1.;
        if(ratio.e != ratio.e) ratio.e = 0.;
        
        h_ratio->SetBinContent(iBin, ratio.v);
        h_ratio->SetBinError(iBin, ratio.e);
    }
    
    return h_ratio;
}