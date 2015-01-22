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
normalizeToNominal_(false),
logX_(false),
logY_(false)
{
    // Suppress default info that canvas is printed
    gErrorIgnoreLevel = 1001;
}



void PlotterSystematic::setOptions(const TString& name, const TString& processNames,
                         const TString& YAxis, const TString& XAxis,
                         const int, const bool normalizeToNominal,
                         const bool logX, const bool logY,
                         const double& ymin, const double& ymax,
                         const double& rangemin, const double& rangemax)
{
    name_ = name; //Histogram name to be plotted
    YAxis_ = YAxis; //Y-axis title
    XAxis_ = XAxis; //X-axis title
    normalizeToNominal_ = normalizeToNominal; //
    logX_ = logX; //Draw X-axis in Log scale
    logY_ = logY; //Draw Y-axis in Log scale
    ymin_ = ymin; //Min. value in Y-axis
    ymax_ = ymax; //Max. value in Y-axis
    rangemin_ = rangemin; //Min. value in X-axis
    rangemax_ = rangemax; //Max. value in X-axis
    
    // Extracting the list of process names
    TObjArray* strings = processNames.Tokenize(" ");
    processNames_.clear();
    for(int i=0; i<strings->GetEntries(); ++i) {
        processNames_.push_back(((TObjString*)strings->At(i))->GetString());
    }
}



void PlotterSystematic::producePlots()
{
    prepareStyle();
    
    
    for(auto channelCollection : inputFileLists_) {
        // Map of nominal histograms for different processes
        std::map<TString, TH1*> m_processHistograms;
            
        Channel::Channel channel = channelCollection.first;
        SystematicHistoMap systematicHistos;
        
        for(TString processName : processNames_) {
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
                    printf("### Warning! No UP histogram for process (%s) in file: %s\n\n", processName.Data(), upDownFile.first.Data());
                }
                TH1* histoDown = fileReader_->GetClone<TH1>(upDownFile.second, name_+"_"+processName, true, false);
                if(!histoDown) {
                    printf("### Warning! No DOWN histogram for process (%s) in file: %s\n\n", processName.Data(), upDownFile.second.Data());
                }
                HistoPair pair(histoUp, histoDown);
                systematicHistos[systematic.type()] = pair;
                
                // Storing the nominal shape of the process to compare processes
                if(systematic.type() != Systematic::nominal) continue;
                m_processHistograms[processName] = (TH1*)histoUp->Clone();
            }
            // Plotting different systematic shapes of each process
            writeVariations(systematicHistos, channel, processName, logY_);
            if(logY_) writeVariations(systematicHistos, channel, processName, true);
        }
        if(processNames_.size()<2) continue;
        // Plotting nominal shapes of different processes
        writeNominalShapes(m_processHistograms, channel);
        if(logY_) writeNominalShapes(m_processHistograms, channel, true);
    }
    
    
    //std::cout<<"\n=== Finishing of plot production\n\n";
}


void PlotterSystematic::writeVariations(const SystematicHistoMap& histoCollection, const Channel::Channel channel, const TString processName, 
                                        const bool logY)
{
    // Getting nominal histogram
    if(histoCollection.count(Systematic::nominal) < 1) {
        std::cout << "No nominal histogram for process: " << processName << " in " << name_ << std::endl;
        return;
    }
    
    // Prepare canvas and legend
    TCanvas* canvas(0);
    canvas = new TCanvas("","");
    canvas->SetName("");
    canvas->SetTitle("");
    TLegend* legend = new TLegend(0.67,0.55,0.92,0.85);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.25);
    legend->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 - legend->GetNRows()*0.05);
    legend->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength());
    legend->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength());
    legend->Clear();
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    canvas->cd(0);
    canvas->Clear();
    if(logY) canvas->SetLogy();
    
    TH1* h_nominal = (TH1*)histoCollection.at(Systematic::nominal).first->Clone();
    if(h_nominal) {
        common::setHistoStyle(h_nominal, 1,1,3, 20,1,0.75, 0,0);
        h_nominal->Draw("HIST E X0");
        legend->AddEntry(h_nominal, processName, "l");
    }
    
    TPad* ratioPad = common::drawRatioPad(canvas, 0.5, 2.0);
    
    // Looping over all available systematics
    size_t orderId = 0;
    for(auto systematicHistos : histoCollection) {
        canvas->cd(0);
        const Systematic::Type systematicType(systematicHistos.first);
        if(systematicHistos.first == Systematic::nominal) continue;
        TH1* h_up = (TH1*)systematicHistos.second.first->Clone();
        if(h_up) {
            if(normalizeToNominal_) common::normalize(h_up, h_nominal->Integral());
            common::setHistoStyle(h_up, lineStyles_.first, lineColors_.at(orderId), 2, 0,-1,-1, 0,0);
            h_up->Draw("same HIST");
        }
        TH1* h_down = (TH1*)systematicHistos.second.second->Clone();
        if(h_down) {
            if(normalizeToNominal_) common::normalize(h_down, h_nominal->Integral());
            common::setHistoStyle(h_down, lineStyles_.second, lineColors_.at(orderId), 2, 0,-1,-1, 0,0);
            h_down->Draw("same HIST");
        }
        
        legend->AddEntry(h_up, Systematic::convertType(systematicType), "l");
        ratioPad->cd();
        TH1* h_ratio_up = common::ratioHistogram(h_up, h_nominal, 0);
        h_ratio_up->Draw("same HIST");
        TH1* h_ratio_down = common::ratioHistogram(h_down, h_nominal, 0);
        h_ratio_down->Draw("same HIST");
        
        // Calculating deviations due to the systematic
        printf(" Variation of the %s:\n", Systematic::convertType(systematicType).Data());
        systematicDeviation(h_nominal, h_up, h_down);

        if(orderId<lineColors_.size()-1) orderId++;
    }
    // Drawing statistical error bars of the nominal histogram
    ratioPad->cd();
    TH1* h_nominal_stat = common::ratioHistogram(h_nominal, h_nominal, 1);
    common::setHistoStyle(h_nominal_stat, 1,1,3, 20,1,0.75);
    h_nominal_stat->Draw("same E1 X0");
    canvas->cd(0);
    h_nominal->Draw("same E1 X0");
    legend->Draw();

    // Updating axis
    updateHistoAxis(canvas);
    
    this->drawCmsLabels(2, 8);
    this->drawDecayChannelLabel(channel);
    
    TString eventFileString = common::assignFolder(outputDir_, channel, Systematic::Systematic("Nominal"))+name_+"_"+processName+"_systematics";
    if(logY) eventFileString.Append("_logY");
    canvas->Print(eventFileString+".eps");
    
    delete canvas;
    delete legend;
}



void PlotterSystematic::systematicDeviation(const TH1* h_nominal, const TH1* h_up, const TH1* h_down)const
{
    if(!h_nominal) return;
    const double inclusive_up = h_up ? h_up->Integral("width") - h_nominal->Integral("width") : 0.;
    const double inclusive_down = h_down ? h_down->Integral("width") - h_nominal->Integral("width") : 0.;
    const double inclusive_error = 0.5*std::max(std::fabs(inclusive_up), std::fabs(inclusive_down));
    std::cout<< "  nominal: " << h_nominal->Integral("width") << "  up: " << h_up->Integral("width") << "  down: " << h_down->Integral("width") << std::endl;
    printf("  Inclusive: %.5f  (%.5f%%)\n", inclusive_error, 100.*inclusive_error/h_nominal->Integral("width"));
}


void PlotterSystematic::writeNominalShapes(const std::map<TString, TH1*>& processHistograms, const Channel::Channel channel, const bool logY)
{
    if(processHistograms.size() < processNames_.size()) {
        std::cout << "Nominal histograms (" << processHistograms.size() << ") not available for all processes (" << processNames_.size() << ")\n";
        return;
    }
    
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
    
    char tempChar[100];
    int histoId = 0;
    // Calculating the total integral of all histograms combined
    double integralTotal = 0.;
    for(TString processName : processNames_) integralTotal += processHistograms.at(processName)->Integral();
    // Plotting histogram for each process
    for(TString processName : processNames_) {
        TH1* histo = (TH1*)processHistograms.at(processName)->Clone();
        double integral = histo->Integral();
        common::normalize(histo);
        common::setHistoStyle(histo, -1, lineColors_.at(histoId), 2, 20,lineColors_.at(histoId),0.5, 0,0);
        
        sprintf(tempChar, "%s [%.1f%%]", processName.Data(), integral/integralTotal*100.);
        legend->AddEntry(histo, tempChar, "l");
        
        if(histoId == 0) histo->Draw("E");
        else histo->Draw("same E");
        
        histoId++;
    }
    legend->Draw();
    
    updateHistoAxis(canvas);
    
    this->drawCmsLabels(2, 8);
    this->drawDecayChannelLabel(channel);
    
    // Saving the plot
    TString eventFileString = common::assignFolder(outputDir_, channel, Systematic::Systematic("Nominal"))+name_+"_nominalShapes";
    if(logY) eventFileString.Append("_logY");
    canvas->Print(eventFileString+".eps");
    
    delete canvas;
    delete legend;
}


void PlotterSystematic::prepareStyle()
{
    // Set style
    common::setHHStyle(*gStyle, false);

    // Different line styles for up/down variations
    lineStyles_.first = 1;
    lineStyles_.second = 7;
    
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


void PlotterSystematic::updateHistoAxis(TPad* pad)const
{

    TH1* histo = common::updatePadYAxisRange(pad, 0.35);

    // Applying the configured X axis range
    if(rangemin_!=0. || rangemax_!=0.) {histo->SetAxisRange(rangemin_, rangemax_, "X");}
    histo->GetXaxis()->SetNoExponent(kTRUE);
    
    // Applying the configured Y axis range
    if(ymax_ != 0.) histo->SetMaximum(ymax_);
    if(ymin_ != 0.) histo->SetMinimum(ymin_);

    // Applying the configured axis titles
    if(XAxis_ != "-") histo->GetYaxis()->SetTitle(YAxis_);
    if(YAxis_ != "-") histo->GetXaxis()->SetTitle(XAxis_);
    
    // Redrawing the axis
    pad->RedrawAxis();
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
    label->AddText(Form(text, 19.7, energy));
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    if(textSize != 0) label->SetTextSize(textSize);
    label->SetTextAlign(32);
    label->Draw("same");
}

