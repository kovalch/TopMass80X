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
    TH1::AddDirectory(false);
    
    // Setting the list of systematics that are included in the tt+HF template fit (only their shape variations are taken into account)
    normalisedSystematicTypes_.push_back(Systematic::jes);
    normalisedSystematicTypes_.push_back(Systematic::pu);
    normalisedSystematicTypes_.push_back(Systematic::lept);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrBpurity);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrLpurity);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrBstat1);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrBstat2);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrLstat1);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrLstat2);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrCerr1);
    normalisedSystematicTypes_.push_back(Systematic::btagDiscrCerr2);
    normalisedSystematicTypes_.push_back(Systematic::xsec_ttcc);
    normalisedSystematicTypes_.push_back(Systematic::xsec_ttother);
    normalisedSystematicTypes_.push_back(Systematic::xsec_ttH);
    normalisedSystematicTypes_.push_back(Systematic::xsec_ttZ);
    normalisedSystematicTypes_.push_back(Systematic::powheg);
    normalisedSystematicTypes_.push_back(Systematic::mass);
    normalisedSystematicTypes_.push_back(Systematic::match);
    
    // Setting the list of systematics that should be compared at generator level without reco selection and sample weights
    // Required due a lower statistics in the varied samples compared to nominal after reco selection
    generatorUnweightedSystematicTypes_.push_back(Systematic::powheg);
    generatorUnweightedSystematicTypes_.push_back(Systematic::mass);
    generatorUnweightedSystematicTypes_.push_back(Systematic::match);
    
    // Setting the list of systematics that should not be included in the total systematic uncertainty
    ignoredSystematicTypes_.push_back(Systematic::powhegHerwig);
    ignoredSystematicTypes_.push_back(Systematic::mcatnlo);
    ignoredSystematicTypes_.push_back(Systematic::perugia11);
    ignoredSystematicTypes_.push_back(Systematic::perugia11NoCR);
    
    // Setting the list of systematics that should be combined 
    // Combination types: 0 - in quadrature; 1 - maximum absolute variation;
    combinedSystematicTypes_[Systematic::btagBeff] = SystematicCombination(0);
    combinedSystematicTypes_.at(Systematic::btagBeff).addSystematic(Systematic::btagDiscrBpurity);
    combinedSystematicTypes_.at(Systematic::btagBeff).addSystematic(Systematic::btagDiscrBstat1);
    combinedSystematicTypes_.at(Systematic::btagBeff).addSystematic(Systematic::btagDiscrBstat2);
    combinedSystematicTypes_.at(Systematic::btagBeff).addSystematic(Systematic::btagDiscrLpurity);
    combinedSystematicTypes_.at(Systematic::btagBeff).addSystematic(Systematic::btagDiscrLstat1);
    combinedSystematicTypes_.at(Systematic::btagBeff).addSystematic(Systematic::btagDiscrLstat2);
    combinedSystematicTypes_.at(Systematic::btagBeff).addSystematic(Systematic::btagDiscrCerr1);
    combinedSystematicTypes_.at(Systematic::btagBeff).addSystematic(Systematic::btagDiscrCerr2);
    combinedSystematicTypes_[Systematic::xsec_ttother] = SystematicCombination(0);
    combinedSystematicTypes_.at(Systematic::xsec_ttother).addSystematic(Systematic::xsec_ttother);
    combinedSystematicTypes_.at(Systematic::xsec_ttother).addSystematic(Systematic::xsec_tt2b);
    combinedSystematicTypes_.at(Systematic::xsec_ttother).addSystematic(Systematic::xsec_ttcc);
    combinedSystematicTypes_.at(Systematic::xsec_ttother).addSystematic(Systematic::xsec_ttH);
    combinedSystematicTypes_.at(Systematic::xsec_ttother).addSystematic(Systematic::xsec_ttZ);
    combinedSystematicTypes_.at(Systematic::xsec_ttother).addSystematic(Systematic::frac_tthf);
    combinedSystematicTypes_.at(Systematic::xsec_ttother).addSystematic(Systematic::frac_ttother);
    
    // Setting the uncertainties that should be forced to particular values
    overridenSystematics_[Systematic::scale] = UpDown(0.08, 0.08);
    
    // Setting the list of systematics that should be used to obtain MC predictions that should be plotted
    predictionSystematicLegends_[Systematic::nominal] = PredictionEntry("Madgraph+Pythia", 2, 1);
    predictionSystematicLegends_[Systematic::powheg] = PredictionEntry("Powheg+Pythia", kAzure+2, 1);
    predictionSystematicLegends_[Systematic::powhegHerwig] = PredictionEntry("Powheg+Herwig", kTeal+4, 7);
    predictionSystematicLegends_[Systematic::mcatnlo] = PredictionEntry("MC@NLO+Herwig", kOrange+7, 7);
    
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
    
    // Getting up/down histograms for each systematic
    for(auto systematicCollection : inputFileLists_.at(channel)) {
        Systematic::Systematic systematic = systematicCollection.first;
        // Skipping PDF variations (need special treatment: later in the code)
        if(systematic.type() == Systematic::pdf) continue;
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
    
    // Getting two up/down histograms out of 53 PDF variations
    TH1* histoUp_pdf = getPdfHisto(histoName, channel, 1, 0);
    TH1* histoDown_pdf = getPdfHisto(histoName, channel, -1, 0);
    if(histoUp_pdf && histoDown_pdf) m_systematicHistos[Systematic::pdf] = HistoPair(histoUp_pdf, histoDown_pdf);
    
    
    return m_systematicHistos;
}


TH1* PlotterDiffXSSystematic::getPdfHisto(TString histoName, const Channel::Channel& channel, const int variation, const int combinationType)const
{
    TH1* h_nominal(0);
    TH1* h_central(0);
    // Getting the Nominal and Central PDF histograms
    for(auto systematicCollection : inputFileLists_.at(channel)) {
        Systematic::Systematic systematic = systematicCollection.first;
        std::pair<TString, TString> upDownFile = systematicCollection.second.at(name_);
        // Getting true nominal
        if(systematic.type() == Systematic::nominal) h_nominal = fileReader_->GetClone<TH1>(upDownFile.first, histoName, true, false);
        // Skipping non-PDF variations
        if(systematic.type() != Systematic::pdf) continue;
        // Getting central PDF variation
        if(systematic.variationNumber() != 0) continue;
        h_central = fileReader_->GetClone<TH1>(upDownFile.first, histoName, true, false);
        break;
    }
    if(!h_nominal || !h_central) return 0;
    
    const int nBins = h_nominal->GetNbinsX();
    TH1* h_varied = (TH1*)h_central->Clone("h_pdf_varied");
    
    // Reading every PDF variation and comparing to the Central
    for(auto systematicCollection : inputFileLists_.at(channel)) {
        Systematic::Systematic systematic = systematicCollection.first;
       // Skipping non-PDF variations
        if(systematic.type() != Systematic::pdf) continue;
        // Skipping central PDF variation
        if(systematic.variationNumber() == 0) continue;
        std::pair<TString, TString> upDownFile = systematicCollection.second.at(name_);
        TH1* h_up = fileReader_->GetClone<TH1>(upDownFile.first, histoName, true, false);
        TH1* h_down = fileReader_->GetClone<TH1>(upDownFile.second, histoName, true, false);
        
        if(combinationType == 1) {
            // Updating every bin of the corresponding histogram if it goes in the proper direction 
            for(int iBin = 1; iBin <= nBins; ++iBin) {
                const double up = h_up->GetBinContent(iBin);
                const double down = h_down->GetBinContent(iBin);
                const double mostVaried = variation > 0 ? std::max(up, down) : std::min(up, down);
                const double lastVaried = h_varied->GetBinContent(iBin);
                const double newContent = variation > 0 ? std::max(lastVaried, mostVaried) : std::min(lastVaried, mostVaried);
                
                h_varied->SetBinContent(iBin, newContent);
            }
        } else if(combinationType == 0) {
            const double upIntegral = h_up->Integral("width");
            const double downIntegral = h_down->Integral("width");
            TH1* h_mostVaried(0);
            if(variation > 0) {
                h_mostVaried = upIntegral > downIntegral ? h_up : h_down;
            } else {
                h_mostVaried = upIntegral < downIntegral ? h_up : h_down;
            }
//             printf("Variation: %d \tup: %.4f  down: %.4f  last: %.4f  most: %.4f  num: %d\n", variation, upIntegral, downIntegral, h_varied->Integral("width"), h_mostVaried->Integral("width"), systematic.variationNumber());
            if(variation > 0 && h_mostVaried->Integral("width") > h_varied->Integral("width")) h_varied = h_mostVaried;
            else if(variation < 0 && h_mostVaried->Integral("width") < h_varied->Integral("width")) h_varied = h_mostVaried;
        }
    }
    // Applying the difference from Central to the Nominal
    for(int iBin = 1; iBin <= nBins; ++iBin) {
	const double factor = h_varied->GetBinContent(iBin)/h_central->GetBinContent(iBin);
	h_nominal->SetBinContent(iBin, factor*h_nominal->GetBinContent(iBin));
    }
    
    delete h_varied;
    delete h_central;
    return h_nominal;
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
            common::normalize(h_up, h_nominal->Integral("width"), false, "width");
            common::normalize(h_down, h_nominal->Integral("width"), false, "width");
        }
        // Calculating difference between nominal and systematic variation in each bin
        for(size_t iBin = 0; iBin < nBins; ++iBin) {
            UpDown errorPair;
            const double nominal = h_nominal->GetBinContent(iBin+1);
            errorPair.c = nominal;
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
        errorPair.c = nominal;
        errorPair.u = error;
        errorPair.d = error;
        errors.at(iBin)[Systematic::nominal] = errorPair;
    }
    
    return errors;
}


void PlotterDiffXSSystematic::plotXSection(const Channel::Channel& channel)
{
    if(inputFileLists_.at(channel).at(Systematic::nominalSystematic()).count(name_) < 1) {
        std::cerr << "  ERROR!!! No nominal histogram for " << name_ << std::endl;
        return;
    }
    TString fileName_nominal = inputFileLists_.at(channel).at(Systematic::nominalSystematic()).at(name_).first;
    // Getting histograms of the measured cross section for all systematic variations
    SystematicHistoMap m_systematicHistos_measured = readSystematicHistograms(name_+"_xs_data", channel);
    SystematicHistoMap m_systematicHistos_generated = readSystematicHistograms(name_+"_xs_madgraph", channel);
    
    // Extracting up/down systematic/statistical errors for each bin of the distribution
    const std::vector<ErrorMap> em_measured = extractVariations(m_systematicHistos_measured);
    const std::vector<ErrorMap> em_generated = extractVariations(m_systematicHistos_generated);
    // The list of uncertainties for each including proper treatment of uncertainties estimated on reco- and gen-level
    std::vector<ErrorMap> em_combined;
    
    TH1* h_nominal = m_systematicHistos_measured.at(Systematic::nominal).first;
    const size_t nBins = h_nominal->GetNbinsX();
    
    // Calculating final uncertainties in each bin
    std::vector<UpDown> e_measured_stat;
    std::vector<UpDown> e_measured_total;
    for(size_t iBin = 0; iBin < nBins; ++iBin) {
        const ErrorMap& m_errors = em_measured.at(iBin);
        const ErrorMap& m_errors_gen = em_generated.at(iBin);
        const ErrorMap m_combined = binUncertaintiesRecoAndGen(m_errors, m_errors_gen);
        const UpDown e_stat = binUncertaintyOfType(m_combined, ErrorType::stat);
        const UpDown e_total = binUncertaintyOfType(m_combined, ErrorType::total);
        e_measured_stat.push_back(e_stat);
        e_measured_total.push_back(e_total);
        em_combined.push_back(m_combined);
    }
    
    // Getting systematic and statistical variations of the inclusive cross section
    SystematicHistoMap m_systematicHistos_measured_inclusive = readSystematicHistograms(name_+"_xs_inclusive_data", channel);
    SystematicHistoMap m_systematicHistos_generated_inclusive = readSystematicHistograms(name_+"_xs_inclusive_madgraph", channel);
    const std::vector<ErrorMap> em_measured_inclusive = extractVariations(m_systematicHistos_measured_inclusive);
    const std::vector<ErrorMap> em_generated_inclusive = extractVariations(m_systematicHistos_generated_inclusive);
    const ErrorMap m_combined_inclusive = binUncertaintiesRecoAndGen(em_measured_inclusive.at(0), em_generated_inclusive.at(0));
    
    // Printing the list of uncertainties in each bin and the median
    // FIXME: Make the BOOL parameter steerable in the HistoDiffXSSystematic (-l to list all systematics)
    printAllUncertainties(em_combined, m_combined_inclusive, true);
    
    TGraphAsymmErrors* g_measured_stat = errorGraphFromHisto(h_nominal, e_measured_stat);
    TGraphAsymmErrors* g_measured_total = errorGraphFromHisto(h_nominal, e_measured_total);
    common::setGraphStyle(g_measured_stat, 1,1,2, 20,1,1);
    common::setGraphStyle(g_measured_total, 1,1,2, 20,1,1);
    
    std::vector<TH1*> prediction_histograms;
    std::vector<TH1*> prediction_ratioHistograms;
    std::vector<TString> prediction_legends;
    for(auto systematicLegend : predictionSystematicLegends_) {
        const Systematic::Type& systematicType = systematicLegend.first;
        const PredictionEntry legendColorStyle = systematicLegend.second;
        const TString fileName = inputFileLists_.at(channel).at(Systematic::Systematic(systematicType, Systematic::undefinedVariation, -1)).at(name_).first;
        // Getting the histogram of prediction
        TH1* histo = fileReader_->GetClone<TH1>(fileName, name_+"_xs_madgraph", true, false);
        // Normalising histogram
        common::normalize(histo, h_nominal->Integral("width"), false, "width");
        common::setHistoStyle(histo, legendColorStyle.style, legendColorStyle.color, 2);
        prediction_histograms.push_back(histo);
        prediction_legends.push_back(legendColorStyle.legend);
        // Adding ratio histogram
        TH1* ratioHisto = common::ratioHistogram(histo, h_nominal);
        if(name_.Contains("Mjj")) ratioHisto->SetBinContent(3, 1.);
        common::setHistoStyle(ratioHisto, legendColorStyle.style, legendColorStyle.color, 2);
        prediction_ratioHistograms.push_back(ratioHisto);
    }

    
    // Prepare canvas and legend
    TCanvas* canvas = new TCanvas("","");
    canvas->Clear();
    TLegend* legend = common::createLegend(0.6, 0.7, 1, 1+prediction_histograms.size(), 0.05);
    
    // Drawing axis and xsections
    for(size_t iHisto = 0; iHisto < prediction_histograms.size(); ++iHisto) {
        if(iHisto == 0) prediction_histograms.at(iHisto)->Draw("hist");
        else prediction_histograms.at(iHisto)->Draw("hist same");
    }
    g_measured_total->Draw("same PZ");
    g_measured_stat->Draw("same ||");
    
    // Drawing legend
    legend->AddEntry(g_measured_total, "Data", "LPE");
    for(size_t iLegend = 0; iLegend < prediction_legends.size(); ++iLegend) {
        legend->AddEntry(prediction_histograms.at(iLegend), prediction_legends.at(iLegend), "L");
    }
    legend->Draw("same");
    
    updateHistoAxis(canvas);
    
    common::drawCmsLabels(-1, 8, 19.7);
    
    // Drawing ratios
    common::drawRatioPad(canvas, 0., double(nRatio_max_), "#frac{Theory}{Data}");
    
    // Plotting the statistics band of the Data
    TGraph* g_ratioData_stat = common::ratioGraph(g_measured_stat, g_measured_stat, 1);
    common::setGraphStyle(g_ratioData_stat, 1,1,1, -1,-1,-1, 1001,18);
    TGraph* g_ratioData_total = common::ratioGraph(g_measured_total, g_measured_total, 1);
    common::setGraphStyle(g_ratioData_total, 1,1,1, -1,-1,-1, 1001,16);
    
    // Performing hardcoded blinding of Mjj
    if(name_.Contains("Mjj")) {
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

    for(TH1* ratioHisto : prediction_ratioHistograms) ratioHisto->Draw("same hist");
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
    
    // Adding all reco-level uncertainties in quadrature
    for(auto systematicValuePair : m_errors) {
        const Systematic::Type& systematic = systematicValuePair.first;
        if(systematic == Systematic::nominal && type == ErrorType::syst) continue;
        if(systematic != Systematic::nominal && type == ErrorType::stat) continue;
        const UpDown& errorPair = systematicValuePair.second;
        error.addInQuadrature(errorPair);
        // Setting the central value as well
        if(type == ErrorType::stat) error.c = errorPair.c;
    }
    error.sqrt();

    return error;
}


PlotterDiffXSSystematic::ErrorMap PlotterDiffXSSystematic::binUncertaintiesRecoAndGen(const ErrorMap& m_errors, const ErrorMap& m_errors_gen)const
{
    ErrorMap m_combined;
    const ErrorMap m_grouped_reco = binUncertainties(m_errors);
    const ErrorMap m_grouped_gen = binUncertainties(m_errors_gen);
    // Adding all reco-level uncertainties in quadrature
    for(auto systematicValuePair : m_grouped_reco) {
        const Systematic::Type& systematic = systematicValuePair.first;
        if(std::find(generatorUnweightedSystematicTypes_.begin(), generatorUnweightedSystematicTypes_.end(), systematic) != 
            generatorUnweightedSystematicTypes_.end()) continue;
        const UpDown& errorPair = systematicValuePair.second;
        m_combined[systematic] = errorPair;
    }
    // Adding gen-level systematics
    for(auto systematicValuePair : m_grouped_gen) {
        const Systematic::Type& systematic = systematicValuePair.first;
        if(std::find(generatorUnweightedSystematicTypes_.begin(), generatorUnweightedSystematicTypes_.end(), systematic) == 
            generatorUnweightedSystematicTypes_.end()) continue;
        const UpDown& errorPair = systematicValuePair.second;
        m_combined[systematic] = errorPair;
    }
    
    return m_combined;
}


PlotterDiffXSSystematic::ErrorMap PlotterDiffXSSystematic::binUncertainties(const ErrorMap& m_errors)const
{
    ErrorMap errors;
    // Looping over all available systematics which should not be merged
    for(auto systematicValuePair : m_errors) {
        Systematic::Type systematicType = systematicValuePair.first;
        // Skipping variation if it should be ignored in systematic uncertainties
        if(std::find(ignoredSystematicTypes_.begin(), ignoredSystematicTypes_.end(), systematicType) != ignoredSystematicTypes_.end()) continue;
        // Skipping this systematic if it should be combined
        bool toBeCombined = false;
        for(auto combinedSystematics : combinedSystematicTypes_) {
            if(!combinedSystematics.second.hasSystematic(systematicType)) continue;
            toBeCombined = true;
            break;
        }
        if(toBeCombined) continue;
        
        errors[systematicType] = systematicValuePair.second;
        
        // Overriding value for the systematic if needed
        if(overridenSystematics_.count(systematicType) < 1) continue;
        errors[systematicType] = overridenSystematics_.at(systematicType);
    }
    // Combining systematics
    for(auto systematicCombinationPair : combinedSystematicTypes_) {
        Systematic::Type combinedSystematic = systematicCombinationPair.first;
        const SystematicCombination& systematicsCombination = systematicCombinationPair.second;
        UpDown combinedError(0., 0.);
        if(systematicsCombination.type == 0) {
            // Adding in quadrature
            for(Systematic::Type systematicType : systematicsCombination.systematics) {
                if(m_errors.count(systematicType) < 1) continue;
                // Skipping variation if it should be ignored in systematic uncertainties
                if(std::find(ignoredSystematicTypes_.begin(), ignoredSystematicTypes_.end(), systematicType) != ignoredSystematicTypes_.end()) continue;
                combinedError.addInQuadrature(m_errors.at(systematicType));
            }
            combinedError.sqrt();
        } else if(systematicsCombination.type == 1) {
            // Taking the envelope of absolute variations
            std::vector<double> ups, downs;
            for(Systematic::Type systematicType : systematicsCombination.systematics) {
                if(m_errors.count(systematicType) < 1) continue;
                // Skipping variation if it should be ignored in systematic uncertainties
                if(std::find(ignoredSystematicTypes_.begin(), ignoredSystematicTypes_.end(), systematicType) != ignoredSystematicTypes_.end()) continue;
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
    
    return errors;
}


void PlotterDiffXSSystematic::printAllUncertainties(const std::vector<ErrorMap>& errorMaps, const ErrorMap& errorMap_inclusive, const bool listSystematics)const
{
    printf("Bin Id \t& cent \t& stat \t& syst \t& total  // highest asymetric uncertainty [ %%x100 ]\n");
    printf("-----------------------------------------------------------------------------------\n");
    // Printing inclusive cross-section results
    const UpDown e_stat_inclusive = binUncertaintyOfType(errorMap_inclusive, ErrorType::stat);
    const UpDown e_syst_inclusive = binUncertaintyOfType(errorMap_inclusive, ErrorType::syst);
    const UpDown e_total_inclusive = binUncertaintyOfType(errorMap_inclusive, ErrorType::total);
    const double xsec_inclusive = e_stat_inclusive.c;
    printf("Incl. \t& %.2f \t& %.2f \t& %.2f \t& %.2f\n", xsec_inclusive, e_stat_inclusive.maxAbsVariation(),  e_syst_inclusive.maxAbsVariation(), e_total_inclusive.maxAbsVariation());
    printf("-----------------------------------------------------------------------------------\n");

    if(errorMaps.size()<1) return;
    // Printing differential cross-section results for each bin
    int binId = 0;
    for(auto errorMap : errorMaps) {
	binId++;
	const UpDown e_stat = binUncertaintyOfType(errorMap, ErrorType::stat);
    	const UpDown e_syst = binUncertaintyOfType(errorMap, ErrorType::syst);
    	const UpDown e_total = binUncertaintyOfType(errorMap, ErrorType::total);
	printf("Bin %d \t& %.2f \t& %.2f \t& %.2f \t& %.2f\n", binId, e_stat.c, e_stat.maxAbsVariation(),  e_syst.maxAbsVariation(), e_total.maxAbsVariation());
    }
    printf("-----------------------------------------------------------------------------------\n");
    
    // Printing the list of systematic uncertainties if required
    if(!listSystematics) return;
    const ErrorMap& systematicsMap = errorMaps.at(0);
    const int numbersIndentation = 20;
    const int nBins = errorMaps.size();
    // Looping over available systematics
    for(auto systematics : systematicsMap) {
        const Systematic::Type systematicType = systematics.first;
        TString systematicName = systematicType == Systematic::nominal ? "STATISTICAL" : Systematic::convertType(systematicType);
        const int systematicLength = systematicName.Length();
        std::vector<double> maximums(0);
        printf("%s", systematicName.Data());
        for(int iSpace = 0; iSpace < numbersIndentation - systematicLength; ++iSpace) printf(" ");
        printf("| ");
        // Looping over all bins
        for(auto errorMap : errorMaps) {
            printf("%.2f  ", errorMap.at(systematicType).maxAbsVariation());
            maximums.push_back(errorMap.at(systematicType).maxAbsVariation());
        }
        printf(" | ");
        // Sorting the up/down variations to have the largest absolute values first
        std::sort(maximums.begin(), maximums.end(), uncertaintySortingFunction);
        const int binMedian = nBins%2 == 0 ? nBins/2 - 1 : nBins/2;
        printf("Median: %.2f ", maximums.at(binMedian) );
        printf(" | ");
        printf("Max: %.2f ", maximums.at(0) );
        printf(" | ");
        printf("Incl: %.2f\n", errorMap_inclusive.at(systematicType).maxAbsVariation());
        if(systematicType == Systematic::nominal) printf("-----------------------------------------------------------------------------------\n");
    }
    // Calculating total uncertainties
    printf("-----------------------------------------------------------------------------------\n");
    printf("Total:");
    for(int iSpace = 0; iSpace < numbersIndentation - 6; ++iSpace) printf(" ");
    printf("| ");
    std::vector<double> maximums(0);
    for(auto errorMap : errorMaps) {
        const UpDown& e_total = binUncertaintyOfType(errorMap, ErrorType::total);
        printf("%.2f  ", e_total.maxAbsVariation());
        maximums.push_back(e_total.maxAbsVariation());
    }
    printf(" | ");
    // Sorting the up/down variations to have the largest absolute values first
    std::sort(maximums.begin(), maximums.end(), uncertaintySortingFunction);
    const int binMedian = nBins%2 == 0 ? nBins/2 - 1 : nBins/2;
    printf("Median: %.2f ", maximums.at(binMedian) );
    printf(" | ");
    printf("Max: %.2f ", maximums.at(0) );
    printf(" | ");
    printf("Incl: %.2f\n", e_total_inclusive.maxAbsVariation());
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
    histo->GetXaxis()->SetLabelSize(histo->GetXaxis()->GetLabelSize()*0.7);;
    histo->GetXaxis()->SetTitleSize(histo->GetXaxis()->GetTitleSize()*0.8);;

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
    return (std::fabs(a) > std::fabs(b));
}

