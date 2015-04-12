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
    
    // Setting the list of systematics that are included in the tt+HF template fit (only their shape variations are taken into account)
    normalisedSystematicTypes_.push_back(Systematic::jes);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrPurity);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrBstat1);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrBstat2);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrLstat1);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrLstat2);
    normalisedSystematicTypes_.push_back(Systematic::xsec_ttcc);
    
    // Setting the list of systematics that should be combined 
    // Combination types: 0 - in quadrature; 1 - maximum variation;
    combinedSystematicTypes_[Systematic::btag] = SystematicCombination(0);
    combinedSystematicTypes_.at(Systematic::btag).addSystematic(Systematic::btagDiscrPurity);
    combinedSystematicTypes_.at(Systematic::btag).addSystematic(Systematic::btagDiscrBstat1);
    combinedSystematicTypes_.at(Systematic::btag).addSystematic(Systematic::btagDiscrBstat2);
    combinedSystematicTypes_.at(Systematic::btag).addSystematic(Systematic::btagDiscrLstat1);
    combinedSystematicTypes_.at(Systematic::btag).addSystematic(Systematic::btagDiscrLstat2);
    
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


PlotterDiffXSSystematic::SystematicHistoMap PlotterDiffXSSystematic::readSystematicHistograms(TString histoName, 
                                                                                              const Channel::Channel& channel)const
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


std::vector<PlotterDiffXSSystematic::ErrorMap> PlotterDiffXSSystematic::extractVariations(const SystematicHistoMap& m_systematicHistos)const
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
    const std::vector<ErrorMap> e_measured = extractVariations(m_systematicHistos_measured);
    
    TH1* h_nominal = m_systematicHistos_measured.at(Systematic::nominal).first;
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
    
    // Getting systematic and statistical variations of the inclusive cross section
    SystematicHistoMap m_systematicHistos_measured_inclusive = readSystematicHistograms(name_+"_xs_inclusive_data", channel);
    const std::vector<ErrorMap> e_measured_inclusive = extractVariations(m_systematicHistos_measured_inclusive);
    const double xsec_inclusive = m_systematicHistos_measured_inclusive.at(Systematic::nominal).first->GetBinContent(1);
    const UpDown e_measured_stat_inclusive = binUncertaintyOfType(e_measured_inclusive.at(0), ErrorType::stat);
    const UpDown e_measured_syst_inclusive = binUncertaintyOfType(e_measured_inclusive.at(0), ErrorType::syst);
    const UpDown e_measured_total_inclusive = binUncertaintyOfType(e_measured_inclusive.at(0), ErrorType::total);
    printf("-----------------------------------------------------------------\n");
    printf("Inclusive xsection: %.3f  +- %.3f (stat)  + %.3f /- %.3f (syst)   [ + %.3f /- %.3f (total)]\n", xsec_inclusive, e_measured_stat_inclusive.u,  e_measured_syst_inclusive.u, e_measured_syst_inclusive.d, e_measured_total_inclusive.u, e_measured_total_inclusive.d);
    printf("-----------------------------------------------------------------\n");
    
    // Printing the list of uncertainties in each bin and the median
    printAllUncertainties(e_measured, e_measured_inclusive.at(0));
    
    TGraphAsymmErrors* g_measured_stat = errorGraphFromHisto(h_nominal, e_measured_stat);
    TGraphAsymmErrors* g_measured_total = errorGraphFromHisto(h_nominal, e_measured_total);
    common::setGraphStyle(g_measured_stat, 1,1,1, 20,1,1);
    common::setGraphStyle(g_measured_total, 1,1,1, 20,1,1);
    
    TH1* h_madgraph = fileReader_->GetClone<TH1>(fileName_nominal, name_+"_xs_madgraph", true, false);
    common::setHistoStyle(h_madgraph, 1,2,1);

    
    // Prepare canvas and legend
    TCanvas* canvas = new TCanvas("","");
    canvas->Clear();
    TLegend* legend = common::createLegend(0.6, 0.6, 1, 2, 0.05);
    
    // Drawing axis and xsections
    h_madgraph->Draw("hist");
    g_measured_total->Draw("same PZ");
    g_measured_stat->Draw("same ||");
    
    // Drawing legend
    legend->AddEntry(g_measured_total, "Data", "LPE");
    legend->AddEntry(h_madgraph, "Madgraph+Pythia", "L");
    legend->Draw("same");
    
    updateHistoAxis(canvas);
    
    common::drawCmsLabels(0, 8, 19.7);
    
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
            UpDown errorPair = systematicValuePair.second;
            if(systematicType == Systematic::nominal && type == ErrorType::syst) continue;
            error.addInQuadrature(errorPair);
        }
        error.sqrt();
        return error;
    }
    
    
    return error;
}


PlotterDiffXSSystematic::ErrorMap PlotterDiffXSSystematic::binUncertaintiesOfType(const ErrorMap& m_errors, const ErrorType type)const
{
    ErrorMap errors;
    // Looping over all available systematics which should not be merged
    for(auto systematicValuePair : m_errors) {
        Systematic::Type systematicType = systematicValuePair.first;
        // Skipping irrelevant variations
        if(type == ErrorType::stat && systematicType != Systematic::nominal) continue;
        if(type == ErrorType::syst && systematicType == Systematic::nominal) continue;
        // Skipping this systematic if it should be combined
        bool toBeCombined = false;
        for(auto combinedSystematics : combinedSystematicTypes_) {
            if(!combinedSystematics.second.hasSystematic(systematicType)) continue;
            toBeCombined = true;
            break;
        }
        if(toBeCombined) continue;
        
        errors[systematicType] = systematicValuePair.second;
    }
    // Combining systematics if necessary
    if(type == ErrorType::syst || type == ErrorType::total) {
        for(auto systematicCombinationPair : combinedSystematicTypes_) {
            Systematic::Type combinedSystematic = systematicCombinationPair.first;
            SystematicCombination systematicsCombination= systematicCombinationPair.second;
            UpDown combinedError(0., 0.);
            if(systematicsCombination.type == 0) {
                // Adding in quadrature
                for(Systematic::Type systematicType : systematicsCombination.systematics) {
                    if(m_errors.count(systematicType) < 1) continue;
                    combinedError.addInQuadrature(m_errors.at(systematicType));
                }
                combinedError.sqrt();
            } else if(systematicsCombination.type == 1) {
                // Taking the envelope of absolute variations
                std::vector<double> ups, downs;
                for(Systematic::Type systematicType : systematicsCombination.systematics) {
                    if(m_errors.count(systematicType) < 1) continue;
                    ups.push_back(m_errors.at(systematicType).u);
                    downs.push_back(m_errors.at(systematicType).d);
                }
                // Sorting to have the largest absolute values first
                std::sort(ups.begin(), ups.end(), uncertaintySortingFunction);
                std::sort(downs.begin(), downs.end(), uncertaintySortingFunction);
                if(ups.size() > 0 && downs.size() > 0) combinedError = UpDown(ups.at(0), downs.at(0));
            }
            errors[combinedSystematic] = combinedError;
        }
    }
    
    return errors;
}


void PlotterDiffXSSystematic::printAllUncertainties(const std::vector<ErrorMap>& errorMaps, const ErrorMap& errorMap_inclusive)const
{
    if(errorMaps.size()<1) return;
    const ErrorMap& systematicsMap = errorMaps.at(0);
    const int numbersIndentation = 20;
    // Looping over available systematics
    for(auto systematics : systematicsMap) {
        const Systematic::Type systematicType = systematics.first;
        TString systematicName = systematicType == Systematic::nominal ? "STATISTICAL" : Systematic::convertType(systematicType);
        const int systematicLength = systematicName.Length();
        std::vector<double> ups, downs;
        printf("%s", systematicName.Data());
        for(int iSpace = 0; iSpace < numbersIndentation - systematicLength; ++iSpace) printf(" ");
        printf("| ");
        // Looping over all bins
        for(auto errorMap : errorMaps) {
            printf("%.3f  ", errorMap.at(systematicType).maxAbsVariation());
            ups.push_back(errorMap.at(systematicType).u);
            downs.push_back(errorMap.at(systematicType).d);
        }
        printf(" | ");
        // Sorting the up/down variations to have the largest absolute values first
        std::sort(ups.begin(), ups.end(), uncertaintySortingFunction);
        std::sort(downs.begin(), downs.end(), uncertaintySortingFunction);
        const int nBins = errorMaps.size();
        const int binMedian = nBins/2;
        printf("Median: %.3f ", std::max(std::fabs(ups.at(binMedian)), std::fabs(downs.at(binMedian))) );
        printf(" | ");
        printf("Incl: %.3f\n", errorMap_inclusive.at(systematicType).maxAbsVariation());
        
//         printf("-----------------------------------------------------------------\n");
//         printf("Total: %.3f\t");
    }
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

bool PlotterDiffXSSystematic::uncertaintySortingFunction(const double a, const double b) {
    return (std::fabs(a) < std::fabs(b));
}

