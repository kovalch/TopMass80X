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
#include <TGraphAsymmErrors.h>
#include <TGaxis.h>
#include <TPaveText.h>
#include <TClass.h>
#include <TError.h>

#include "PlotterDiffXSSystematic.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"
#include "../../common/include/plotterUtils.h"





PlotterDiffXSSystematic::PlotterDiffXSSystematic(const char* outputDir,
                                     const std::map<Channel::Channel, std::map<Systematic::Systematic, std::map<TString, std::pair<TString, TString> > > >& inputFileLists):
outputDir_(outputDir),
inputFileLists_(inputFileLists),
fileReader_(RootFileReader::getInstance()),
name_("defaultName"),
rangemin_(0),
rangemax_(3),
ymin_(0),
ymax_(0),
nRatio_max_(3.),
YAxis_(""),
XAxis_(""),
logX_(false),
logY_(false)
{
    // Suppress default info that canvas is printed
    gErrorIgnoreLevel = 1001;
    
    // Setting the list of systematics that are included in the tt+HF template fit (only shape variations should be taken into account)
    normalisedSystematicTypes_.push_back(Systematic::jes);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrPurity);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrBstat1);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrBstat2);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrLstat1);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrLstat2);
    normalisedSystematicTypes_.push_back(Systematic::xsec_ttcc);
}



void PlotterDiffXSSystematic::setOptions(const TString& name, const TString&,
                         const TString& YAxis, const TString& XAxis,
                         const int nRatio_max, const bool,
                         const bool logX, const bool logY,
                         const double& ymin, const double& ymax,
                         const double& rangemin, const double& rangemax)
{
    name_ = name; //Histogram name to be plotted
    YAxis_ = YAxis; //Y-axis title
    XAxis_ = XAxis; //X-axis title
    nRatio_max_ = nRatio_max == 0. ? nRatio_max_ : nRatio_max; //Upper bound of the Y axis in the ratio plot
    logX_ = logX; //Draw X-axis in Log scale
    logY_ = logY; //Draw Y-axis in Log scale
    ymin_ = ymin; //Min. value in Y-axis
    ymax_ = ymax; //Max. value in Y-axis
    rangemin_ = rangemin; //Min. value in X-axis
    rangemax_ = rangemax; //Max. value in X-axis
}



void PlotterDiffXSSystematic::producePlots()
{
    common::setHHStyle(*gStyle, false);
    
    for(auto channelCollection : inputFileLists_) {
        Channel::Channel channel = channelCollection.first;
        plotXSection(channel);
    }
}


PlotterDiffXSSystematic::SystematicHistoMap PlotterDiffXSSystematic::readSystematicHistograms(TString histoName, const Channel::Channel& channel)
{
    SystematicHistoMap m_systematicHistos;
    if(inputFileLists_.count(channel) < 1) {
        std::cerr << "No histograms (" << histoName << ") available for the channel: " << Channel::convert(channel) << std::endl;
        return m_systematicHistos;
    }
    
    for(auto systematicCollection : inputFileLists_.at(channel)) {
        Systematic::Systematic systematic = systematicCollection.first;
        if(systematicCollection.second.count(name_) < 1) {
            printf("### Warning! No file (%s) for %s : %s\n\n", name_.Data(), Systematic::convertType(systematic.type()).Data(), Channel::convert(channel).Data());
            continue;
        }
        std::pair<TString, TString> upDownFile = systematicCollection.second.at(name_);
        TH1* histoUp = fileReader_->GetClone<TH1>(upDownFile.first, histoName, true, false);
        if(!histoUp) {
            printf("### Warning! No UP histogram (%s) in file: %s\n\n", histoName.Data(), upDownFile.first.Data());
        }
        TH1* histoDown = fileReader_->GetClone<TH1>(upDownFile.second, histoName, true, false);
        if(!histoDown) {
            printf("### Warning! No DOWN histogram (%s) in file: %s\n\n", histoName.Data(), upDownFile.second.Data());
        }
        HistoPair pair(histoUp, histoDown);
        m_systematicHistos[systematic.type()] = pair;
    }
    
    return m_systematicHistos;
}


std::vector<PlotterDiffXSSystematic::ErrorMap> PlotterDiffXSSystematic::extractVariations(const SystematicHistoMap& m_systematicHistos)
{
    std::vector<ErrorMap> errors;
    
    // Getting nominal histogram
    if(m_systematicHistos.count(Systematic::nominal) < 1) {
        std::cerr << "No nominal histogram for " << name_ << std::endl;
        return errors;
    }
    
    TH1* h_nominal = (TH1*)m_systematicHistos.at(Systematic::nominal).first->Clone();
    if(!h_nominal) {
        std::cerr << "Nominal histogram is NULL for " << name_ << std::endl;
        return errors;
    }
    const size_t nBins = h_nominal->GetNbinsX();
    // Initialising the empty errors for each bin
    errors = std::vector<ErrorMap>(nBins);
    
    // Extracting up/down systematic relative variations in each bin
    for(auto systematicHistoPair : m_systematicHistos) {
        Systematic::Type systematicType = systematicHistoPair.first;
        
        TH1* h_up = systematicHistoPair.second.first;
        TH1* h_down= systematicHistoPair.second.second;
        // Normalising systematic histograms to Nominal if only shape variation should be taken into account
        if(std::find(normalisedSystematicTypes_.begin(), normalisedSystematicTypes_.end(), systematicType) != normalisedSystematicTypes_.end()) {
            common::normalize(h_up, h_nominal->Integral());
            common::normalize(h_down, h_nominal->Integral());
        }
        // Calculating difference between nominal and systematic variation in each bin
        for(size_t iBin = 0; iBin < nBins; ++iBin) {
            UpDown errorPair;
            const double nominal = h_nominal->GetBinContent(iBin+1);
            if(h_up) errorPair.u = (h_up->GetBinContent(iBin+1) - nominal)/nominal;
            if(h_down) errorPair.d = (h_down->GetBinContent(iBin+1) - h_nominal->GetBinContent(iBin+1))/nominal;
            errors.at(iBin)[systematicType] = errorPair;
        }
    }
    // Storing the statistical uncertainty of the nominal distribution as a Nominal systematic
    for(size_t iBin = 0; iBin < nBins; ++iBin) {
        UpDown errorPair;
        const double nominal = h_nominal->GetBinContent(iBin+1);
        const double error = h_nominal->GetBinError(iBin+1) / nominal;
        errorPair.u = error;
        errorPair.d = error;
        errors.at(iBin)[Systematic::nominal] = errorPair;
    }
    
    return errors;
}


void PlotterDiffXSSystematic::plotXSection(const Channel::Channel& channel)
{
    TString fileName_nominal = inputFileLists_.at(channel).at(Systematic::nominalSystematic()).at(name_).first;
    // Getting histograms of the measured cross section for all systematic variations
    SystematicHistoMap m_systematicHistos_measured = readSystematicHistograms(name_+"_xs_data", channel);
    
    // Extracting up/down systematic/statistical errors for each bin of the distribution
    std::vector<ErrorMap> e_measured = extractVariations(m_systematicHistos_measured);
    
    TH1* h_nominal= m_systematicHistos_measured.at(Systematic::nominal).first;
    const size_t nBins = h_nominal->GetNbinsX();
    
    // Calculating final uncertainties in each bin
    std::vector<UpDown> e_measured_stat;
    std::vector<UpDown> e_measured_total;
    for(size_t iBin = 0; iBin < nBins; ++iBin) {
        const ErrorMap& m_errors = e_measured.at(iBin);
        const UpDown e_stat = binUncertaintyOfType(m_errors, ErrorType::stat);
        const UpDown e_total = binUncertaintyOfType(m_errors, ErrorType::total);
        e_measured_stat.push_back(e_stat);
        e_measured_total.push_back(e_total);
    }
    
    TGraphAsymmErrors* g_measured_stat = errorGraphFromHisto(h_nominal, e_measured_stat);
    TGraphAsymmErrors* g_measured_total = errorGraphFromHisto(h_nominal, e_measured_total);
    common::setGraphStyle(g_measured_stat, 1,1,1, 20,1,1);
    common::setGraphStyle(g_measured_total, 1,1,1, 20,1,1);
    
    TH1* h_madgraph = fileReader_->GetClone<TH1>(fileName_nominal, name_+"_xs_madgraph", true, false);
    common::setHistoStyle(h_madgraph, 1,2,1);

    
    // Prepare canvas and legend
    TCanvas* canvas = new TCanvas("","");
    TLegend* legend = new TLegend(0.67,0.7,0.92,0.85);
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
    
    // Drawing axis and xsections
    h_madgraph->Draw("hist");
    g_measured_total->Draw("same PZ");
    g_measured_stat->Draw("same ||");
    
    // Drawing legend
    legend->AddEntry(g_measured_total, "Data", "LPE");
    legend->AddEntry(h_madgraph, "Madgraph+Pythia", "L");
    legend->Draw("same");
    
    updateHistoAxis(canvas);
    
    this->drawCmsLabels(2, 8);
    this->drawDecayChannelLabel(channel);
    
    // Drawing ratios
    common::drawRatioPad(canvas, 0., double(nRatio_max_), "#frac{MC}{Data}");
    
    // Plotting the statistics band of the Data
    TGraph* g_ratioData_stat = common::ratioGraph(g_measured_stat, g_measured_stat, 1);
    common::setGraphStyle(g_ratioData_stat, 1,1,1, -1,-1,-1, 1001,18);
    TGraph* g_ratioData_total = common::ratioGraph(g_measured_total, g_measured_total, 1);
    common::setGraphStyle(g_ratioData_total, 1,1,1, -1,-1,-1, 1001,16);
    TH1* h_ratioMadgraphData = common::ratioHistogram(h_madgraph, h_nominal);
    common::setHistoStyle(h_ratioMadgraphData, 1,2,1);
    
    // Performing hardcoded blinding of Mjj
    if(name_.Contains("Mjj")) {
        h_ratioMadgraphData->SetBinContent(3, 1.);
        double x, y;
        g_measured_total->GetPoint(2, x, y);
        g_measured_total->SetPoint(2, x, 1e-10);
        g_measured_total->SetPointEYhigh(2, 0.);
        g_measured_total->SetPointEYlow(2, 0.);
        g_measured_stat->SetPoint(2, x, 1e-10);
        g_measured_stat->SetPointEYhigh(2, 0.);
        g_measured_stat->SetPointEYlow(2, 0.);
        g_ratioData_total->SetPoint(2, x, -1e10);
    }

    h_ratioMadgraphData->Draw("same hist");
    g_ratioData_total->Draw("same PZ");
    g_ratioData_stat->Draw("same ||");
    gPad->RedrawAxis();
    
    // Saving the plot
    TString eventFileString = common::assignFolder(outputDir_, channel, Systematic::Systematic("Nominal"))+name_+"_systematic";
    canvas->Print(eventFileString+".eps");
    
    delete canvas;
    delete legend;
}


PlotterDiffXSSystematic::UpDown PlotterDiffXSSystematic::binUncertaintyOfType(const ErrorMap& m_errors, const ErrorType type)const
{
    UpDown error(0.,0.);
    
    // Statistical uncertainty is stored under Nominal systematic
    if(type == ErrorType::stat) {
        if(m_errors.count(Systematic::nominal) < 1) return error;
        else return m_errors.at(Systematic::nominal);
    } else 
    // All or only systematic uncertainties added in quadrature
    if(type == ErrorType::syst || type == ErrorType::total) {
        for(auto systematicValuePair : m_errors) {
            Systematic::Type systematicType = systematicValuePair.first;
            if(systematicType == Systematic::nominal && type == ErrorType::syst) continue;
            error.addInQuadrature(systematicValuePair.second);
        }
        error.sqrt();
        return error;
    }
    
    return error;
}


TGraphAsymmErrors* PlotterDiffXSSystematic::errorGraphFromHisto(const TH1* histo, const std::vector<UpDown>& errors)const
{
    const int nBins = errors.size();
    if(nBins != histo->GetNbinsX()) {
        std::cerr << "setErrorsToHisto: Number of bins errors and histogram do not match. Stopping..." << std::endl;
        exit(102);
    }
    
    TGraphAsymmErrors* graph = new TGraphAsymmErrors(histo);
    
    // Updating the up/down Y errors of the graph
    for(int iBin = 0; iBin < nBins; ++iBin) {
        const double value = histo->GetBinContent(iBin+1);
        // Multiplying error by bin content (errors are relative)
        graph->SetPointEYhigh(iBin, errors.at(iBin).u*value);
        graph->SetPointEYlow(iBin, errors.at(iBin).d*value);
    }
    
    return graph;
}


void PlotterDiffXSSystematic::updateHistoAxis(TPad* pad)const
{

    TH1* histo = common::updatePadYAxisRange(pad, logY_, 0.35);

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


void PlotterDiffXSSystematic::drawDecayChannelLabel(const Channel::Channel& channel, const double& textSize)const
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


void PlotterDiffXSSystematic::drawCmsLabels(const int cmsprelim, const double& energy, const double& textSize)const
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

