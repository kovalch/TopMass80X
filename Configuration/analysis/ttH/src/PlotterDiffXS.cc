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

#include "PlotterDiffXS.h"
#include "higgsUtils.h"
#include "Samples.h"
#include "DilepSVDFunctions.h"
#include "TheoryTopHistoReader.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"
#include "../../common/include/plotterUtils.h"





PlotterDiffXS::PlotterDiffXS(const char* outputDir,
                             const Samples& samples,
                             const double luminosity):
outputDir_(outputDir),
samples_(samples),
fileReader_(RootFileReader::getInstance()),
inputDirTheoryTop_(tth::DATA_PATH_TTH() + "/" + "theoryPredictions"),
hasPseudodata_(false),
name_("defaultName"),
nameGen_("defaultName"),
nameGenEventBased_("events_weighted_step0b"),
signalType_(0),
plotResponse_(false),
normalizeXS_(false),
rangemin_(0),
rangemax_(3),
ymin_(0),
ymax_(0),
sampleTypesData_(std::vector<Sample::SampleType>(0)),
sampleTypesSignal_(std::vector<Sample::SampleType>(0)),
sampleTypesBackground_(std::vector<Sample::SampleType>(0)),
sampleTypesTtbar_(std::vector<Sample::SampleType>(0)),
YAxis_(""),
XAxis_(""),
logX_(false),
logY_(false),
luminosity_(luminosity)
{
    // Suppress default info that canvas is printed
    gErrorIgnoreLevel = 1001;
}



void PlotterDiffXS::setOptions(const TString& name, const TString& namesGen,
                         const TString& YAxis, const TString& XAxis,
                         const int signalType, const bool plotResponse,
                         const bool normalizeXS, const bool logY,
                         const double& ymin, const double& ymax,
                         const double& rangemin, const double& rangemax)
{
    name_ = name; //Variable name for the reconstructed level
    if(name_.Contains("#")) name_.Remove(name_.Last('#')); //Remove possible #TEXT in the end (to allow multiple entries of the same histogram)
    YAxis_ = YAxis; //Y-axis title
    XAxis_ = XAxis; //X-axis title
    signalType_ = signalType; //Type of the signal (0-ttbb, 1-ttbb+ttb+tt2b)
    plotResponse_ = plotResponse; //Whether the response matrix should be plotted instead of the cross section
    normalizeXS_ = normalizeXS; // Whether normalized cross sections should be plotted (for shape comparison)
//     logX_ = logX; //Draw X-axis in Log scale
    logY_ = logY; //Draw Y-axis in Log scale
    ymin_ = ymin; //Min. value in Y-axis
    ymax_ = ymax; //Max. value in Y-axis
    rangemin_ = rangemin; //Min. value in X-axis
    rangemax_ = rangemax; //Max. value in X-axis
    
    // Extracting all generator level histogram names to be used
    namesTheory_.clear();
    TObjArray* genNames = namesGen.Tokenize("|");
    for(int nameId = 0; nameId < genNames->GetEntries(); ++nameId) {
        // The first name is for the generator level histogram used for XS calculation
        if(nameId == 0) nameGen_ = ((TObjString*)genNames->At(nameId))->GetString();
        // Other names are for theory predictions in *.top format
        else namesTheory_.push_back(((TObjString*)genNames->At(nameId))->GetString());
    }
    
    // Clearing the sample collections from previous run of the plotter
    sampleTypesData_.clear();
    sampleTypesSignal_.clear();
    sampleTypesBackground_.clear();
    sampleTypesTtbar_.clear();
    
    // Setting the lists of sample types which should be treated as data, signal or background
    for(int type = Sample::data; type <= Sample::dummy; ++type) {
        if(type==Sample::data || type==Sample::pseudodata) sampleTypesData_.push_back(Sample::SampleType(type));
        else if(type==Sample::ttbb && (signalType_==0 || signalType_==1)) sampleTypesSignal_.push_back(Sample::SampleType(type));
        else if(type==Sample::ttb && signalType_==1) sampleTypesSignal_.push_back(Sample::SampleType(type));
        else if(type==Sample::tt2b && signalType_==1) sampleTypesSignal_.push_back(Sample::SampleType(type));
//         else if(type==Sample::ttbb || type==Sample::ttb || type==Sample::tt2b || type==Sample::ttcc || type==Sample::ttother) sampleTypesSignal_.push_back(Sample::SampleType(type));
        else sampleTypesBackground_.push_back(Sample::SampleType(type));
        
        if(type==Sample::ttbb || type==Sample::ttb || type==Sample::tt2b || type==Sample::ttcc || type==Sample::ttother)
            sampleTypesTtbar_.push_back(Sample::SampleType(type));
    }
    
}



void PlotterDiffXS::producePlots()
{
    //std::cout<<"--- Beginning of plot production\n\n";
    
    // Set style
    common::setHHStyle(*gStyle, false);
    gROOT->ForceStyle();
    
    // Access correction factors
    const SystematicChannelFactors globalWeights = this->scaleFactors();
    const SystematicChannelFactors globalWeights_noScale = this->scaleFactors(false, false);
    
    
    // Loop over all channels and systematics and produce plots
    const SystematicChannelSamples& m_systematicChannelSample(samples_.getSystematicChannelSamples());
    for(const auto& systematicChannelSamples : m_systematicChannelSample){
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        for(const auto& channelSample : systematicChannelSamples.second){
            const Channel::Channel& channel(channelSample.first);
            if(channel != Channel::combined) continue;
            const std::vector<Sample>& v_sample(channelSample.second);
            const std::vector<double>& v_weight(globalWeights.at(systematic).at(channel));
            if(!this->prepareDataset(v_sample, v_weight, v_sampleHistPair_, name_)){
                std::cout<<"WARNING! Cannot find histograms for all datasets, for (channel/systematic): "
                         << Channel::convert(channel) << "/" << systematic.name()
                         <<"\n... skip this plot\n";
                return;
            }
            if(plotResponse_) {
                this->writeResponseMatrix(channel, systematic);
                continue;
            }
            const std::vector<double>& v_weight_noScale(globalWeights_noScale.at(systematic).at(channel));
            this->prepareDataset(v_sample, v_weight_noScale, v_sampleHistPair_noWeight_, name_);
            if(nameGen_ == "-") {
                std::cout<<"WARNING! No name of the generator level histogram provided, for (channel/systematic): "
                         << Channel::convert(channel) << "/" << systematic.name()
                         <<"\n... skip this plot\n";
                return;
            }
            // Getting generator level histograms for the specified quantity (single channel to avoid double-counting)
            const std::vector<double>& v_weight_noScale_ee(globalWeights_noScale.at(systematic).at(Channel::ee));
            const std::vector<double>& v_weight_ee(globalWeights.at(systematic).at(Channel::ee));
            const std::vector<Sample>& v_sample_ee = m_systematicChannelSample.at(systematic).at(Channel::ee);
            this->prepareDataset(v_sample_ee, v_weight_noScale_ee, v_sampleHistPairGen_noWeight_, nameGen_);
            this->prepareDataset(v_sample_ee, v_weight_ee, v_sampleHistPairGen_, nameGen_);
            
            // Getting generator level histograms based on events (to get ttbb/tt ratio from MC prediction)
            if(!this->prepareDataset(v_sample_ee, v_weight_noScale_ee, v_sampleHistPairGenEventBased_, nameGenEventBased_)) {
                std::cout<<"WARNING! Cannot find generator level event based histograms for all datasets, for (channel/systematic): "
                         << Channel::convert(channel) << "/" << systematic.name()
                         <<"\n... skip this plot\n";
                return;
            }
            // Getting response 2D histograms for the cross section
            TString nameResponse(name_);
            nameResponse.ReplaceAll("_step", "VsGen_step");
            if(!this->prepareDataset(v_sample, v_weight, v_sampleHistPairResponse_, nameResponse)) {
                std::cout<<"WARNING! Cannot find response histogram for all datasets, for (channel/systematic): "
                         << Channel::convert(channel) << "/" << systematic.name()
                         <<"\n... skip this plot\n";
                return;
            }
            
            std::cout << "Systematic: " << systematic.name() << "  Channel: " << Channel::convert(channel) << std::endl;
            this->writeDiffXS(channel, systematic);
        }
    }
    
    //std::cout<<"\n=== Finishing of plot production\n\n";
}



const SystematicChannelFactors PlotterDiffXS::scaleFactors(const bool dyScale, const bool hfScale)
{
    const TString fullStepname = tth::extractSelectionStepAndJetCategory(name_);
    
    return samples_.globalWeights(fullStepname, dyScale, hfScale).first;
}



bool PlotterDiffXS::prepareDataset(const std::vector<Sample>& v_sample,
                                   const std::vector<double>& v_weight,
                                   std::vector<SampleHistPair>& v_sampleHistPair,
                                   const TString name)
{
    bool allHistosAvailable(true);
    
    // Associate histogram to dataset if histogram can be found
    v_sampleHistPair.clear();
    TH1::AddDirectory(kFALSE);
    for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
        const auto& sample(v_sample.at(iSample));
        SampleHistPair p_sampleHist;
        TH1 *hist = fileReader_->GetClone<TH1>(sample.inputFile(), name, true, false);
        if(!hist){
            // Print message only for one histogram
            if(allHistosAvailable)
                std::cout<<"Histogram ("<<name<<") not found e.g. in file ("<<sample.inputFile()<<")\n";
            p_sampleHist = SampleHistPair(sample, 0);
            allHistosAvailable = false;
        } else {
            if(hist->GetSumw2N() < 1) hist->Sumw2();
            // Apply weights
            if(sample.sampleType() != Sample::data){
                const double& weight = v_weight.at(iSample);
                hist->Scale(weight);
            }

            // Clone histogram directly here
            TH1* histClone = (TH1*) hist->Clone();
            p_sampleHist = SampleHistPair(sample, histClone);
        }
        if(sample.sampleType() == Sample::pseudodata) hasPseudodata_ = true;
        v_sampleHistPair.push_back(p_sampleHist);
    }
    //std::cout<<"Number of samples used for histogram ("<<name<<"): "<<v_sampleHistPair.size()<<"\n";

    return allHistosAvailable;
}




void PlotterDiffXS::writeDiffXS(const Channel::Channel& channel, const Systematic::Systematic& systematic)
{
    // Prepare canvas and legend
    TCanvas* canvas = new TCanvas("canvas_diffXS","");
    canvas->Clear();
    
    
    // Generic information about the histogram being plotted
    TString histo_name = v_sampleHistPair_.size() > 0 ? TString(v_sampleHistPair_.at(0).second->GetName()) : "";
    
    
    // Loop over all samples and add those with identical legendEntry
    // And sort them into the categories data, signal, background
    
    // Setting up a map of histograms used to calculate the xsection
    std::map<TString, TH1*> m_inputHistos;
    
    // Create histogram corresponding to the sum of all stacked data histograms
    std::vector<LegendHistPair> dataHists = legendHistPairsForSamples(sampleTypesData_, v_sampleHistPair_);
    m_inputHistos["data"] = sumOfHistos(dataHists, "h_data");
    
    // Create histogram corresponding to the sum of all stacked signal histograms
    std::vector<LegendHistPair> signalHists = legendHistPairsForSamples(sampleTypesSignal_, v_sampleHistPair_);
    m_inputHistos["signal"] = sumOfHistos(signalHists, "h_signal_reco");
    
    // Create histogram corresponding to the sum of all stacked background histograms
    std::vector<LegendHistPair> backgroundHists = legendHistPairsForSamples(sampleTypesBackground_, v_sampleHistPair_);
    m_inputHistos["background"] = sumOfHistos(backgroundHists, "h_background_reco");
    
    // Create histogram corresponding to the sum of all stacked signal histograms without sample level weights
    std::vector<LegendHistPair> signalHists_noWeight = legendHistPairsForSamples(sampleTypesSignal_, v_sampleHistPair_noWeight_);
    m_inputHistos["signal_noWeight"] = sumOfHistos(signalHists_noWeight, "h_signal_reco_noWeight");
    
        
    // Get all signal samples at generator level
    std::vector<LegendHistPair> signalHistsGen = legendHistPairsForSamples(sampleTypesSignal_, v_sampleHistPairGen_);
    m_inputHistos["signal_gen"] = sumOfHistos(signalHistsGen, "h_signal_gen");
    
    // Get all signal samples at generator level without sample level weights
    std::vector<LegendHistPair> signalHistsGen_noWeight = legendHistPairsForSamples(sampleTypesSignal_, v_sampleHistPairGen_noWeight_);
    m_inputHistos["signal_gen_noWeight"] = sumOfHistos(signalHistsGen_noWeight, "h_signal_gen_noWeight");
    
    // Get all ttbar dileptonic samples at generator level (to get fraction of ttbb in ttbar from the generator)
    std::vector<LegendHistPair> ttbarHistsGen = legendHistPairsForSamples(sampleTypesTtbar_, v_sampleHistPairGenEventBased_);
    m_inputHistos["ttbar_gen"] = sumOfHistos(ttbarHistsGen, "h_ttbar_gen");
    
    // Get response matrix
    std::vector<LegendHistPair> signalHists_response = legendHistPairsForSamples(sampleTypesSignal_, v_sampleHistPairResponse_);
    m_inputHistos["signal_response"] = (TH2*)sumOfHistos(signalHists_response, "h_signal_response");
    
    // Calculating actual differential cross section and getting corresponding histograms
    std::map<TString, TH1*> m_xs = calculateDiffXS(m_inputHistos, true);
    std::map<TString, TH1*> m_xs_reweighted;
    // Calculating differential cross section again using the reweighted generator MC
    if(hasPseudodata_) {
        // Get all signal samples at generator level that are reweighted for pseudodata (only for closure tests)
        std::vector<LegendHistPair> signalHistsGenReweighted_noWeight = legendHistPairsForSamples(sampleTypesData_, v_sampleHistPairGen_noWeight_);
        m_inputHistos["signal_gen_noWeight"] = sumOfHistos(signalHistsGenReweighted_noWeight, "h_signal_gen_reweighted_noWeight");
        printf("Integral gen reweighted: %.4f\n", m_inputHistos.at("signal_gen_noWeight")->Integral());
        m_xs_reweighted = calculateDiffXS(m_inputHistos, true);
    }
    
    
    // Obtaining the theory curve
    std::vector<TString> theoryLegends({"PowHel"});
    std::map<TString, TH1*> m_theoryHisto;
    TheoryTopHistoReader theoryReader;
    for(size_t nameId = 0; nameId < namesTheory_.size(); ++nameId) {
        TString fileHistoName = namesTheory_.at(nameId);
        TString fileName = fileHistoName(0, fileHistoName.First(':'));
        TString histoName = fileHistoName(fileHistoName.First(':')+1, fileHistoName.Length());
        TString histoLegend = theoryLegends.at(nameId);
        
        TH1* histo = theoryReader.getHisto1D(inputDirTheoryTop_+"/"+fileName, histoName);
        // Matching the histogram's binning and normalization
        histo = common::rebinHistoToHisto(histo, m_xs.at("data"));
        common::normalizeToBinWidth(histo, true);
        double norm_scale = common::normalize(histo, m_xs.at("inclusive_data")->Integral(), false);
        common::normalizeToBinWidth(histo, false);
        histo->SetBinContent(0, norm_scale);
        m_theoryHisto[histoLegend] = histo;
    }
    
    

    if(!m_xs.size()) {
        std::cerr<<"ERROR in Plotter! Histograms for the cross section not produced\n...break\n"<<std::endl;
        exit(237);
    }

    TLegend* legend = common::createLegend(0.7, 0.6, 1, 2+(int)hasPseudodata_+m_theoryHisto.size(), 0.035);
    char legendEntry[100];
    // Add entries to legend
    legend->AddEntry(m_xs.at("data"), "Data","pe");
    
    const double norm_scale_mc = m_xs.at("madgraph")->GetBinContent(0);
    sprintf(legendEntry, "Madgraph x%.2f", norm_scale_mc);
    legend->AddEntry(m_xs.at("madgraph"), legendEntry,"l");
    
    // Adding entries for the theory curves
    for(auto nameHisto : m_theoryHisto) {
        const double norm_scale = nameHisto.second->GetBinContent(0);
        sprintf(legendEntry, nameHisto.first+" x%.2f", norm_scale);
        legend->AddEntry(nameHisto.second, legendEntry,"l");
    }

    if(logY_) canvas->SetLogy();
    common::setHistoStyle(m_xs.at("data"), 1,1,1, 20,1,1.5, 0,0);
    common::setHistoStyle(m_xs.at("madgraph"), 1,2,1, 0,2,1.5, 0,0);
    for(auto nameHisto : m_theoryHisto) {
        common::setHistoStyle(nameHisto.second, 1,kAzure+2,1, 0,kAzure+2,1.5, 0,0);
    }
    
    // Draw cross section histograms
    canvas->cd();
    // Canceling 0 X error, that is globally set in the common::setHistoStyle
    gStyle->SetErrorX();
    m_xs.at("data")->Draw("E1 X0");
    m_xs.at("madgraph")->Draw("hist same");
    for(auto nameHisto : m_theoryHisto) nameHisto.second->Draw("hist same");
    if(hasPseudodata_) {
        sprintf(legendEntry, "Madgraph (reweighted) x%.2f", norm_scale_mc);
        legend->AddEntry(m_xs_reweighted.at("madgraph"), legendEntry,"l");
        // Applying the same scale for reweighted MC shape
        m_xs_reweighted.at("madgraph")->Scale(norm_scale_mc);
        common::setHistoStyle(m_xs_reweighted.at("madgraph"), 7,kAzure+2,1, 0,2,1.5, 0,0);
        m_xs_reweighted.at("madgraph")->Draw("hist same");
    }
    m_xs.at("data")->Draw("same E1 X0");
    
    TH1* axisHisto = common::updatePadYAxisRange(canvas, logY_);
    updateHistoAxis(axisHisto);
    
    // Put additional stuff to histogram
    common::drawCmsLabels(0, 8, samples_.luminosityInInversePb()/1000.);
    legend->Draw("same");
    common::drawRatioPad(canvas, 0., 2.4, "#frac{MC}{Data}");
    
    // Plotting the statistics band of the Data
    TH1* h_ratioData_statBand = common::ratioHistogram(m_xs.at("data"), m_xs.at("data"), 1);
    common::setHistoStyle(h_ratioData_statBand, 1,1,1, 0,0,0, 3354,13);
    h_ratioData_statBand->Draw("same E2");
//     h_ratioData_statBand->Draw("L same");
    
    TH1* h_ratioMadgraphData = common::ratioHistogram(m_xs.at("madgraph"), m_xs.at("data"));
    if(hasPseudodata_) {
        common::setHistoStyle(h_ratioMadgraphData, 1,2,1, 0,2,1.5, 0,0);
        TH1* h_ratioMadgraphData_reweighted = common::ratioHistogram(m_xs_reweighted.at("madgraph"), m_xs.at("data"));
        common::setHistoStyle(h_ratioMadgraphData_reweighted, 7,kAzure+2,1, 0,2,1.5, 0,0);
        h_ratioMadgraphData->Draw("hist same");
        h_ratioMadgraphData_reweighted->Draw("hist same");
    } else {
        h_ratioMadgraphData->Draw("hist same");
    }
    for(auto nameHisto : m_theoryHisto) {
        TH1* h_ratioDataTheory = common::ratioHistogram(nameHisto.second, m_xs.at("data"));
        h_ratioDataTheory->Draw("hist E0 same");
    }

    // Create Directory for Output Plots and write them
    TString xsAddName = "_diffXS";
    if(normalizeXS_) xsAddName.Append("_norm");
    TString eventFileString = common::assignFolder(outputDir_, channel, systematic)+name_+xsAddName;
    canvas->Print(eventFileString+".eps");
    
    //Save input and cross section histograms to the root file
    TFile out_root(eventFileString+"_source.root", "RECREATE");
    // Storing input histograms
    for(auto nameHisto : m_inputHistos) {
        if(!nameHisto.second) continue;
        nameHisto.second->Write(name_+xsAddName+"_input_"+nameHisto.first);
    }
    // Storing measured XS histograms
    for(auto nameHisto : m_xs) {
        if(!nameHisto.second) continue;
        nameHisto.second->Write(name_+xsAddName+"_xs_"+nameHisto.first);
    }
    // Storing XS histograms from theory predictions
    for(auto nameHisto : m_theoryHisto) {
        if(!nameHisto.second) continue;
        nameHisto.second->Write(name_+xsAddName+"_xsTheory_"+nameHisto.first);
    }
    
    out_root.Close();
    
    // Removing created objects
    m_xs.at("data")->Delete();
    m_xs.at("madgraph")->Delete();
    canvas->Clear();
    legend->Clear();
    
    delete canvas;
    delete legend;
}




std::map<TString, TH1*> PlotterDiffXS::calculateDiffXS(const std::map<TString, TH1*> m_inputHistos, const bool normalizeMcToData)const
{
    std::map<TString, TH1*> result;
    
    // Extracting histograms from the input map
    TH1* h_data = m_inputHistos.at("data");
    TH1* h_signal = m_inputHistos.at("signal");
    TH1* h_background = m_inputHistos.at("background");
    TH1* h_signal_noWeight = m_inputHistos.at("signal_noWeight");
    TH1* h_signal_gen = m_inputHistos.at("signal_gen");
    TH1* h_signal_gen_noWeight = m_inputHistos.at("signal_gen_noWeight");
    TH1* h_ttbar_gen = m_inputHistos.at("ttbar_gen");
    
    // Initial values for the cross section calculation
    const ValueError luminosity(luminosity_, 0.);
    const ValueError br_dilepton(0.06391, 0.00063);
    const ValueError xsection_tt_dilepton(241.5*br_dilepton.v, 8.5*br_dilepton.v);
//     const ValueError xsection_tt_dilepton(234000.0*br_dilepton.v, 14290.0*br_dilepton.v);
    
    TH1* h_dataMinusBackground = (TH1*)h_data->Clone("h_dataMinusBackground");
    h_dataMinusBackground->Add(h_background, -1);
    
    TH1* h_mc = (TH1*)h_signal->Clone("h_signalPlusBackground_reco");
    h_mc->Add(h_background, 1.0);
    
    
    const int nBins = h_data->GetNbinsX();
    ////////////////////////////////////////////////////////////////////// CALCULATING THE INCLUSIVE CROSS SECTION: ttbb
    ValueError N_data_reco(1.,1.);
    N_data_reco.v = h_data->IntegralAndError(1,-1, N_data_reco.e);
    ValueError N_dataMinusBackground(1.,1.);
    N_dataMinusBackground.v = h_dataMinusBackground->IntegralAndError(1,-1, N_dataMinusBackground.e);
    ValueError N_signal_reco(1.,1.);
    N_signal_reco.v = h_signal->IntegralAndError(1,-1, N_signal_reco.e);
    ValueError N_background_reco(1.,1.);
    N_background_reco.v = h_background->IntegralAndError(1,-1, N_background_reco.e);
    ValueError N_mc_reco(1.,1.);
    N_mc_reco.v = h_mc->IntegralAndError(1,-1, N_mc_reco.e);
    ValueError N_signal_reco_noWeight(1.,1.);
    N_signal_reco_noWeight.v = h_signal_noWeight->IntegralAndError(1,-1, N_signal_reco_noWeight.e);
    ValueError N_signal_gen(1.,1.);
    N_signal_gen.v = h_signal_gen->IntegralAndError(1,-1, N_signal_gen.e);
    ValueError N_signal_gen_noWeight(1.,1.);
    N_signal_gen_noWeight.v = h_signal_gen_noWeight->IntegralAndError(1,-1, N_signal_gen_noWeight.e);
    ValueError N_ttbar_gen_events(1.,1.);
    N_ttbar_gen_events.v = h_ttbar_gen->IntegralAndError(1,-1, N_ttbar_gen_events.e);
    
    printf("\n### [Full reco selection, generator weights * PU weights * reco weights * sample SFs (DY, tt+HF)]\n");
    printf("Data:   \t%.3f \t+- %.2f\n", N_data_reco.v, N_data_reco.e);
    printf("Signal: \t%.3f \t+- %.2f\t(RECO signal)\n", N_signal_reco.v, N_signal_reco.e);
    printf("Background: \t%.3f \t+- %.2f\t(RECO rest MC)\n", N_background_reco.v, N_background_reco.e);
    printf("-------------------------------\n");
    printf("MC total: \t%.3f \t+- %.2f\n", N_mc_reco.v, N_mc_reco.e);
    printf("\nData-Bkg: \t%.3f \t+- %.2f\n\n", N_dataMinusBackground.v, N_dataMinusBackground.e);

    printf("### [Full reco selection, generator weights * PU weights * reco weights]\n");
    printf("Signal[reco]: \t%.3f \t+- %.2f\n", N_signal_reco_noWeight.v, N_signal_reco_noWeight.e);
    printf("### [No reco selection, generator weights * PU weights * reco weights]\n");
    printf("Signal[gen]: \t%.3f \t+- %.2f\n", N_signal_gen.v, N_signal_gen.e);
    printf("### [No reco selection, generator weights * PU weights]\n");
    printf("Signal[gen]: \t%.3f \t+- %.2f\n", N_signal_gen_noWeight.v, N_signal_gen_noWeight.e);
    printf("-----------------------\n");
    printf("Acceptance: \t%.3f\n\n", N_signal_reco_noWeight.v/N_signal_gen_noWeight.v);
    printf("tt+jets[gen]: \t%.3f \t+- %.2f\n", N_ttbar_gen_events.v, N_ttbar_gen_events.e);
    printf("-----------------------\n");
    
    ValueError ttbbFraction_Madgraph;
    ttbbFraction_Madgraph.v = N_signal_gen_noWeight.v/N_ttbar_gen_events.v;
    ttbbFraction_Madgraph.e = uncertaintyBinomial(N_signal_gen.v, N_ttbar_gen_events.v);
    
    printf("ttbb/tt: \t%.3f \t+- %.2f\n\n", ttbbFraction_Madgraph.v, ttbbFraction_Madgraph.e);
    
    printf("Luminosity: \t%.3f \t+- %.2f\n", luminosity.v, luminosity.e);
//     printf("BR(tt->ll): \t%.3f \t+- %.2f\n\n", br_dilepton.v, br_dilepton.e);
    
    ValueError xsection_inclusive;
    xsection_inclusive.v = N_dataMinusBackground.v / (luminosity.v * N_signal_reco.v / N_signal_gen.v);
    xsection_inclusive.e = N_dataMinusBackground.e / (luminosity.v * N_signal_reco.v / N_signal_gen.v);
    ValueError xsection_inclusive_mc;
    xsection_inclusive_mc.v = N_signal_reco.v / (luminosity.v * N_signal_reco.v / N_signal_gen.v);
    xsection_inclusive_mc.e = N_signal_reco.e / (luminosity.v * N_signal_reco.v / N_signal_gen.v);
    ValueError xsection_inclusive_madgraph;
    xsection_inclusive_madgraph.v = xsection_tt_dilepton.v * N_signal_gen_noWeight.v / N_ttbar_gen_events.v;
    xsection_inclusive_madgraph.e = xsection_tt_dilepton.e * N_signal_gen_noWeight.v / N_ttbar_gen_events.v;
    
    printf("###############################\n");
    printf("Inclusive x-section (Data):   \t%.3f \t+- %.3f\n", xsection_inclusive.v, xsection_inclusive.e);
    printf("Inclusive x-section (MC reco):\t%.3f \t+- %.3f\n", xsection_inclusive_mc.v, xsection_inclusive_mc.e);
    printf("Inclusive x-section (MC gen): \t%.3f \t+- %.3f\n\n", xsection_inclusive_madgraph.v, xsection_inclusive_madgraph.e);
    
    // Creating histograms with a single bin for inclusive cross sections
    TH1* h_incXS_data = new TH1D("h_incXS_data", ";Data;#sigma_{inclusive}", 1, 0, 1);
    h_incXS_data->SetBinContent(1, xsection_inclusive.v);
    h_incXS_data->SetBinError(1, xsection_inclusive.e);
    TH1* h_incXS_madgraph = new TH1D("h_incXS_madgraph", ";Madgraph;#sigma_{inclusive}", 1, 0, 1);
    h_incXS_madgraph->SetBinContent(1, xsection_inclusive_madgraph.v);
    h_incXS_madgraph->SetBinError(1, xsection_inclusive_madgraph.e);
    
    result["inclusive_data"] = h_incXS_data;
    result["inclusive_madgraph"] = h_incXS_madgraph;
    
    
    ////////////////////////////////////////////////////////////////////// CALCULATING THE DIFFERENTIAL CROSS SECTION
    // Calculating the cross section from data
    TH1* h_diffXS_data = unfoldedHistogram(m_inputHistos, 1);
    
    for(int iBin = 1; iBin <= nBins; ++iBin) {
        const ValueError N_signal_reco( h_diffXS_data->GetBinContent(iBin), 
                                        h_diffXS_data->GetBinError(iBin));
        ValueError xsection(1.,1.);
        xsection.v = N_signal_reco.v / (luminosity.v);
        xsection.e = N_signal_reco.e / (luminosity.v);
        if(xsection.v != xsection.v) xsection.v = 0.;
        if(xsection.e != xsection.e) xsection.e = 0.;
        
        h_diffXS_data->SetBinContent(iBin, xsection.v);
        h_diffXS_data->SetBinError(iBin, xsection.e);
    }
    
    // Calculating the cross section from MC
    TH1* h_diffXS_madgraph = (TH1*)h_data->Clone("h_diffXS_madgraph");
    for(int iBin = 0; iBin <= nBins+1; ++iBin) {
        ValueError signalFraction;
        signalFraction.v = h_signal_gen_noWeight->GetBinContent(iBin)/N_ttbar_gen_events.v;
        signalFraction.e = uncertaintyBinomial(h_signal_gen_noWeight->GetBinContent(iBin), N_ttbar_gen_events.v);
        
        ValueError xsection(0., 0.);
        xsection.v = xsection_tt_dilepton.v * signalFraction.v;
        xsection.e = sqrt(xsection_tt_dilepton.eOv2() * signalFraction.v2() + signalFraction.eOv2());
        if(xsection.v != xsection.v) xsection.v = 0.;
        if(xsection.e != xsection.e) xsection.e = 0.;
        
        h_diffXS_madgraph->SetBinContent(iBin, xsection.v);
        h_diffXS_madgraph->SetBinError(iBin, xsection.e);
    }
    // Scaling to the measured inclusive cross section (for shape comparison)
    double normScale = 1.0;
    if(normalizeMcToData) {
        normScale = xsection_inclusive.v/xsection_inclusive_madgraph.v;
        h_diffXS_madgraph->Scale(normScale);
    }
    
    // Normalising to unity
    if(normalizeXS_) {
        common::normalize(h_diffXS_data);
        common::normalize(h_diffXS_madgraph);
    }
    
    // Normalising to bin width
    common::normalizeToBinWidth(h_diffXS_data);
    common::normalizeToBinWidth(h_diffXS_madgraph);
    
    
    // Setting the normalisation scale to the underflow bin
    h_diffXS_madgraph->SetBinContent(0, normScale);
    
    // Blinding the Mjj data histogram
    if(name_.Contains("Mjj") && !hasPseudodata_) {
        h_diffXS_data->SetBinContent(3, 1e-10);
        h_diffXS_data->SetBinError(3, 0);
    }
    
    result["data"] = h_diffXS_data;
    result["madgraph"] = h_diffXS_madgraph;
    
    
    // Removing created objects
    h_dataMinusBackground->Delete();
    h_mc->Delete();
    
    return result;
}



TH1* PlotterDiffXS::unfoldedHistogram(const std::map<TString, TH1*> m_inputHistos, const int unfoldingType)const
{
    TH1* histo_unfolded(0);
    
    // Custom bin-by-bin unfolding
    if(unfoldingType == 0) {
        TH1* histo = (TH1*)m_inputHistos.at("data")->Clone("h_data_unfolded");
        histo->Reset("ICESM");
        TH1* h_dataMinusBackground = (TH1*)m_inputHistos.at("data")->Clone();
        h_dataMinusBackground->Add(m_inputHistos.at("background"), -1);
        
        for(int iBin = 1; iBin <= histo->GetNbinsX(); ++iBin) {
            const ValueError N_signal_reco( h_dataMinusBackground->GetBinContent(iBin), 
                                            h_dataMinusBackground->GetBinError(iBin));
            const ValueError acceptance( m_inputHistos.at("signal_noWeight")->GetBinContent(iBin)/m_inputHistos.at("signal_gen_noWeight")->GetBinContent(iBin), 
                                         uncertaintyBinomial( m_inputHistos.at("signal_noWeight")->GetBinContent(iBin), 
                                                              m_inputHistos.at("signal_gen_noWeight")->GetBinContent(iBin) ) );
            ValueError content(0., 0.);
            content.v = N_signal_reco.v / acceptance.v;
            content.e = N_signal_reco.e / acceptance.v;
            if(content.v != content.v) content.v = 0.;
            if(content.e != content.e) content.e = 0.;
            
            histo->SetBinContent(iBin, content.v);
            histo->SetBinError(iBin, content.e);
        }
        
        return histo;
    }
    // Unfolding as defined in the DilepSVDFunctions.cc
    else if(unfoldingType == 1) {
        // Setting the unfolding function
        TH1D* histogram(0);
        TH1D* histogram_norm(0);
        // Getting the binning
        const int nBins = m_inputHistos.at("data")->GetNbinsX();
        double bins[nBins+1];
        for(int iBin=0; iBin<nBins; ++iBin) {
            bins[iBin] = m_inputHistos.at("data")->GetBinLowEdge(iBin+1);
            // Ensuring that data always has more entries than background
            if(m_inputHistos.at("data")->GetBinContent(iBin+1) < m_inputHistos.at("background")->GetBinContent(iBin+1)) {
                printf("####################### WARNING! More background than data in bin starting from %.2f\n", bins[iBin]);
                printf("####################### Differential cross-section histogram is reset\n\n");
                TH1* histo = (TH1*)m_inputHistos.at("data")->Clone("h_data_unfolded");
                histo->Reset("ICESM");
                return histo;
            }
        }
        bins[nBins] = bins[nBins-1] + m_inputHistos.at("data")->GetBinWidth(nBins);
        double* bins_ptr = bins;
        
        // Setting the particle and quantity name
        TString particleName(name_), quantityName(name_);
        if(particleName.Contains("1st")) particleName = "1st"; 
        else if(particleName.Contains("2nd")) particleName = "2nd";
        else particleName = "dijet";
        if(quantityName.Contains("Pt")) quantityName = "Pt";
        else if(quantityName.Contains("Eta")) quantityName = "Eta";
        else if(quantityName.Contains("M")) quantityName = "M";
        else if(quantityName.Contains("dR")) quantityName = "dR";
        
        
        TH1D* h_backgroundSignal = new TH1D("h_backgroundSignal", "", nBins, bins);
        h_backgroundSignal->Reset();
        
        DilepSVDFunctions mySVDFunctions;
        mySVDFunctions.SetOutputPath(outputDir_);
        
//         printf("Bins of the response matrix:\n");
//         TH2* h_response = (TH2*)m_inputHistos.at("signal_response");
//         for(int bin1 = 0; bin1 <= 1+h_response->GetNbinsX(); ++bin1) {
//             for(int bin2 = 0; bin2 <= 1+h_response->GetNbinsY(); ++bin2) {
//                 printf("  (%d,%d):  %.3f \t+/-  %.4f\n", bin1, bin2, h_response->GetBinContent(bin1, bin2), h_response->GetBinError(bin1, bin2));
//             }
//         }
        
        
        mySVDFunctions.SVD_DoUnfold(
            (TH1D*)m_inputHistos.at("data"),
            (TH1D*)m_inputHistos.at("background"),
            (TH1D*)h_backgroundSignal,
            (TH1D*)m_inputHistos.at("signal_gen"),
            (TH1D*)m_inputHistos.at("signal"),
            (TH2D*)m_inputHistos.at("signal_response"),
            NULL, NULL, NULL, NULL, NULL,
            bins_ptr, nBins,
            histogram, histogram_norm,
            0, "combined", particleName, quantityName, name_, ""
        );
        
        // Updating bin contents since unfolded histogram has outside bins shifted into visible area
        TH1* histo = (TH1*)m_inputHistos.at("data")->Clone("h_data_unfolded");
        histo->Reset("ICESM");
        for(int iBin = 1; iBin <= histo->GetNbinsX(); ++iBin) {
            histo->SetBinContent(iBin, histogram->GetBinContent(iBin+1));
            histo->SetBinError(iBin, histogram->GetBinError(iBin+1));
        }
        
        delete h_backgroundSignal;
        
        return histo;
    }
    
    return histo_unfolded;
}


void PlotterDiffXS::writeResponseMatrix(const Channel::Channel& channel, const Systematic::Systematic& systematic)
{
    // Prepare canvas and legend
    TCanvas* canvas = new TCanvas("canvas_response","");
    canvas->Clear();
    canvas->SetName("");
    canvas->SetTitle("");
    
    // Loop over all specified samples and add those with identical legendEntry
    std::vector<LegendHistPair> signalHists = legendHistPairsForSamples(sampleTypesSignal_, v_sampleHistPair_);
    
    // Create histogram corresponding to the sum of all specified samples
    TH1* h_response = sumOfHistos(signalHists, "h_ttbar_reco");

    if(!h_response){
        std::cerr<<"ERROR in PlotterDiffXS! No single sample for drawing exists\n...break\n"<<std::endl;
        exit(237);
    }

    if(logY_) canvas->SetLogy();    
    updateHistoAxis(h_response);
    

    // Scaling the content to number of entries for the response matrix
    h_response->Scale(h_response->GetEntries()/((TH2*)h_response)->Integral(0,-1,0,-1));

    // Draw data histogram and stack and error bars
    h_response->SetLineColor(0);
    h_response->Draw("colz");
    // Moving the pad to the left to fit the color palette
    if(ymax_ != 0) h_response->SetMaximum(ymax_);
    gPad->SetRightMargin(0.15);
    
    // Put additional stuff to histogram
    common::drawCmsLabels(0, 8, samples_.luminosityInInversePb()/1000.);

    // Create Directory for Output Plots and write them
    const TString eventFileString = common::assignFolder(outputDir_, channel, systematic);
    canvas->Print(eventFileString+name_+".eps");
    
    // Plotting Purity/Stability if this is a response matrix
    drawPurityStability((TH2*)h_response, eventFileString+name_+"_PurStab"+".eps");
    
    h_response->Delete();
    
    canvas->Clear();
    
    delete canvas;
    
}



double PlotterDiffXS::uncertaintyBinomial(const double pass, const double all)const 
{
    return (1./all)*sqrt(pass - pass*pass/all);
}


TGraphErrors* PlotterDiffXS::purityStabilityGraph(TH2* h2d, const int type)const
{

    const int nBins = h2d->GetNbinsX();

    TGraphErrors* graph = new TGraphErrors(nBins);

    // Calculating each point of graph for each diagonal bin
    for(int iBin = 1; iBin<=nBins; ++iBin) {
        const double diag = h2d->GetBinContent(iBin, iBin);
        const double reco = h2d->Integral(iBin, iBin, 1, -1)+1e-30;
        const double gen = h2d->Integral(1, -1, iBin, iBin)+1e-30;
        
        const double value = type == 0 ? diag/reco : diag/gen;
        const double error = type == 0 ? uncertaintyBinomial(diag, reco) : uncertaintyBinomial(diag, gen);

        const double bin = h2d->GetXaxis()->GetBinCenter(iBin);
        const double binW = h2d->GetXaxis()->GetBinWidth(iBin);
        
        graph->SetPoint(iBin-1, bin, value);
        graph->SetPointError(iBin-1, binW/2., error);
    }

    return graph;

}


std::vector<PlotterDiffXS::LegendHistPair> PlotterDiffXS::legendHistPairsForSamples(const std::vector<Sample::SampleType> allowedSampleTypes, 
                                                                                    const std::vector<SampleHistPair> samples,
                                                                                    const bool addMCtoPseudodata)const
{
    std::vector<LegendHistPair> result;
    TH1* tmpHist(0);
    
    Sample dataSample;
    bool isPseudodata = false;
    
    for(std::vector<SampleHistPair>::const_iterator i_sampleHistPair = samples.begin();
    i_sampleHistPair != samples.end(); ++i_sampleHistPair)
    {
        const Sample& sample = i_sampleHistPair->first;
        const Sample::SampleType& sampleType(sample.sampleType());
        if((!addMCtoPseudodata || !isPseudodata) && std::find(allowedSampleTypes.begin(), allowedSampleTypes.end(), sampleType) == allowedSampleTypes.end()) continue;
        
        const TString& legendEntry(sample.legendEntry());
        const TH1* hist = i_sampleHistPair->second;
        
        if(!hist) {
            std::cerr << "ERROR!!! Histogram not found for sample: " << legendEntry << std::endl;
            exit(1);
        }
        
        // Detecting whether pseudodata is used instead of data
        if(sampleType == Sample::pseudodata) {
            dataSample = sample;
            isPseudodata = true;
        }
        // Skipping other MC samples if it should be added
        if(isPseudodata && !addMCtoPseudodata && sampleType!=Sample::pseudodata) continue;
        // Skipping files from MC samples that are already included in pseudodata
        if(isPseudodata && sampleType!=Sample::pseudodata && dataSample.containsFilenamesOfSample(sample, true)) continue;
        
        std::vector<SampleHistPair>::const_iterator incrementIterator(i_sampleHistPair);
        ++incrementIterator;
        const bool lastHist(incrementIterator == samples.end());
        const bool newHist(!tmpHist);

        if(newHist)tmpHist = (TH1*)hist->Clone();

        if(lastHist || (legendEntry!=incrementIterator->first.legendEntry())){
            if(!newHist)tmpHist->Add(hist);
            // Adding histogram to the list
            result.push_back(LegendHistPair(legendEntry, tmpHist));
            tmpHist = 0;
        }
        else{
            if(!newHist)tmpHist->Add(hist);
        }
    }
    
    return result;
}


TH1* PlotterDiffXS::sumOfHistos(const std::vector<LegendHistPair> histos, const TString name)const
{
    TH1* sum(0);
    for(size_t i = 0; i < histos.size(); ++i){
        if(i==0) sum = (TH1*)histos.at(0).second->Clone(name);
        else sum->Add((TH1*)histos.at(i).second);
    }
    
    return sum;
}


void PlotterDiffXS::drawPurityStability(TH2* histo2d, TString name)const 
{
    TGraphErrors* g_purity = purityStabilityGraph(histo2d, 0);
    TGraphErrors* g_stability = purityStabilityGraph(histo2d, 1);
    
    
    // Styling axis
    TH1* histo = histo2d->ProjectionX(TString(histo2d->GetName())+"_px");
    histo->SetAxisRange(0., 1.1, "Y");
    histo->GetYaxis()->SetTitle("Values");
    gPad->SetRightMargin(0.05);
    
    histo->Draw("AXIS");
    
    // Styling graphs
    common::setGraphStyle(g_purity, 1,1,2, 21,1,1.5);
    common::setGraphStyle(g_stability, 1,2,2, 20,2,1.5);
    
    // Drawing graphs
    g_purity->Draw("P0same");
    g_stability->Draw("P0same");
    gPad->Update();
    
    
    // Adding a legend
    TLegend* leg = common::createLegend(0.6, 0.6, 1, 2, 0.05);
    leg->AddEntry(g_purity, "Purity", "p");
    leg->AddEntry(g_stability, "Stability", "p");
    leg->Draw();
    
    common::drawCmsLabels(0, 8, samples_.luminosityInInversePb()/1000.);
    
    // Storing the same canvas with a different name
    gPad->Print(name);
    
    delete g_purity;
    delete g_stability;
    delete leg;
    delete histo;
}


void PlotterDiffXS::updateHistoAxis(TH1* histo)const
{
    // Set x and y axis ranges
    if(rangemin_!=0. || rangemax_!=0.) {histo->SetAxisRange(rangemin_, rangemax_, "X");}
    
    if(ymin_!=0) histo->SetMinimum(ymin_);
    if(ymax_!=0.) histo->SetMaximum(ymax_);

    histo->GetXaxis()->SetNoExponent(kTRUE);
    // Set axis titles
    if(XAxis_ != "-") histo->GetYaxis()->SetTitle(YAxis_);
    if(YAxis_ != "-") histo->GetXaxis()->SetTitle(XAxis_);
}

