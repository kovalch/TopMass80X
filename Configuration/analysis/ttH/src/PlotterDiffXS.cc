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
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"
#include "../../common/include/plotterUtils.h"





PlotterDiffXS::PlotterDiffXS(const char* outputDir,
                             const Samples& samples,
                             const double luminosity):
outputDir_(outputDir),
samples_(samples),
fileReader_(RootFileReader::getInstance()),
name_("defaultName"),
nameGen_("defaultName"),
nameGenEventBased_("events_weighted_step0b"),
mode_(0),
rebin_(1),
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



void PlotterDiffXS::setOptions(const TString& name, const TString& nameGen,
                         const TString& YAxis, const TString& XAxis,
                         const int rebin, const bool mode,
                         const bool logX, const bool logY,
                         const double& ymin, const double& ymax,
                         const double& rangemin, const double& rangemax)
{
    name_ = name; //Variable name for the reconstructed level
    nameGen_ = nameGen; //Variable name for the generator level
    YAxis_ = YAxis; //Y-axis title
    XAxis_ = XAxis; //X-axis title
    rebin_ = rebin; //Nr. of bins to be merged together
    logX_ = logX; //Draw X-axis in Log scale
    logY_ = logY; //Draw Y-axis in Log scale
    ymin_ = ymin; //Min. value in Y-axis
    ymax_ = ymax; //Max. value in Y-axis
    rangemin_ = rangemin; //Min. value in X-axis
    rangemax_ = rangemax; //Max. value in X-axis
    mode_ = mode; //Mode for the histogram (0 - 1D cross section, 1 - 2D response matrix)
    
    // Setting the lists of sample types which should be treated as data, signal or background
    for(int type = Sample::data; type != Sample::dummy; ++type) {
        if(type==Sample::data) sampleTypesData_.push_back(Sample::SampleType(type));
        else if(type==Sample::ttbb ) sampleTypesSignal_.push_back(Sample::SampleType(type));
        else sampleTypesBackground_.push_back(Sample::SampleType(type));
        
        if(type==Sample::ttbb || type==Sample::ttb || type==Sample::tt2b || type==Sample::ttcc || type==Sample::ttother)
            sampleTypesTtbar_.push_back(Sample::SampleType(type));
    }
    
}



void PlotterDiffXS::producePlots()
{
    //std::cout<<"--- Beginning of plot production\n\n";
    
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
            if(mode_==1) {
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
            // Getting generator level histograms for the specified quantity
            const std::vector<double>& v_weight_noScale_ee(globalWeights_noScale.at(systematic).at(Channel::ee));
            const std::vector<Sample>& v_sample_ee = m_systematicChannelSample.at(systematic).at(Channel::ee);
            if(!this->prepareDataset(v_sample_ee, v_weight_noScale_ee, v_sampleHistPairGen_, nameGen_)) {
                std::cout<<"WARNING! Cannot find generator level histograms for all datasets, for (channel/systematic): "
                         << Channel::convert(channel) << "/" << systematic.name()
                         <<"\n... skip this plot\n";
                return;
            }
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
        }
        else{
            // Apply weights
            if(sample.sampleType() != Sample::data){
                const double& weight = v_weight.at(iSample);
                hist->Scale(weight);
            }

            // Set style
            common::setHHStyle(*gStyle);

            // Clone histogram directly here
            TH1* histClone = (TH1*) hist->Clone();
            p_sampleHist = SampleHistPair(sample, histClone);
        }
        v_sampleHistPair.push_back(p_sampleHist);
    }
    //std::cout<<"Number of samples used for histogram ("<<name<<"): "<<v_sampleHistPair.size()<<"\n";


    return allHistosAvailable;
}




void PlotterDiffXS::writeDiffXS(const Channel::Channel& channel, const Systematic::Systematic& systematic)
{
    // Prepare canvas and legend
    TCanvas* canvas = new TCanvas("","");
    TLegend* legend = new TLegend(0.70,0.7,0.92,0.85);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.25);
    legend->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 - legend->GetNRows()*0.04);
    legend->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength());
    legend->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength());
    legend->Clear();
    canvas->Clear();
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    canvas->SetName("");
    canvas->SetTitle("");
    
    
    // Here fill colors and line width are adjusted, and potentially rebinning applied
    if(rebin_>1) {
        for(auto sampleHistPair : v_sampleHistPair_) sampleHistPair.second->Rebin(rebin_);
    }
    
    
    // Generic information about the histogram being plotted
    TString histo_name = v_sampleHistPair_.size() > 0 ? TString(v_sampleHistPair_.at(0).second->GetName()) : "";
    
    
    // Loop over all samples and add those with identical legendEntry
    // And sort them into the categories data, signal, background
    
    // Create histogram corresponding to the sum of all stacked signal histograms
    std::vector<LegendHistPair> dataHists = legendHistPairsForSamples(sampleTypesData_, v_sampleHistPair_);
    TH1* h_data = sumOfHistos(dataHists, "h_data");
    
    // Create histogram corresponding to the sum of all stacked signal histograms
    std::vector<LegendHistPair> signalHists = legendHistPairsForSamples(sampleTypesSignal_, v_sampleHistPair_);
    TH1* h_signal = sumOfHistos(signalHists, "h_signal_reco");
    
    // Create histogram corresponding to the sum of all stacked background histograms
    std::vector<LegendHistPair> backgroundHists = legendHistPairsForSamples(sampleTypesBackground_, v_sampleHistPair_);
    TH1* h_background = sumOfHistos(backgroundHists, "h_background_reco");
    
    // Create histogram corresponding to the sum of all stacked signal histograms without sample level weights
    std::vector<LegendHistPair> signalHists_noWeight = legendHistPairsForSamples(sampleTypesSignal_, v_sampleHistPair_noWeight_);
    TH1* h_signal_noWeight = sumOfHistos(signalHists_noWeight, "h_signal_reco_noWeight");
    
        
    // Get all signal samples at generator level
    std::vector<LegendHistPair> signalHistsGen = legendHistPairsForSamples(sampleTypesSignal_, v_sampleHistPairGen_);
    TH1* h_signal_gen = sumOfHistos(signalHistsGen, "h_signal_gen");
    
    // Get all ttbar dileptonic samples at generator level (to get fraction of ttbb in ttbar from the generator)
    std::vector<LegendHistPair> ttbarHistsGen = legendHistPairsForSamples(sampleTypesTtbar_, v_sampleHistPairGenEventBased_);
    TH1* h_ttbar_gen = sumOfHistos(ttbarHistsGen, "h_ttbar_gen");
    
    // Get response matrix
    std::vector<LegendHistPair> ttbarHists_response = legendHistPairsForSamples(sampleTypesTtbar_, v_sampleHistPairResponse_);
    TH1* h_ttbar_response = sumOfHistos(ttbarHists_response, "h_ttbar_response");
    
    // Calculating actual differential cross section and getting corresponding histograms
    std::map<TString, TH1*> xsectionSet= calculateDiffXS(h_data, h_signal, h_background, h_signal_noWeight, h_signal_gen, h_ttbar_response, h_ttbar_gen);
    

    if(!xsectionSet.size()) {
        std::cerr<<"ERROR in Plotter! Histograms for the cross section not produced\n...break\n"<<std::endl;
        exit(237);
    }

    
    // Obtaining the proper histogram to be plotted
    TH1* h_xs_data = xsectionSet.at("diffXS_data");
    TH1* h_xs_mc = xsectionSet.at("diffXS_madgraph");
    const double norm_scale_mc = h_xs_mc->GetBinContent(0);
    char legendEntry[100];
    sprintf(legendEntry, "Madgraph x%.1f", norm_scale_mc);
    
    
    // Add entries to legend
    legend->AddEntry(h_xs_data, "Data","pe");
    legend->AddEntry(h_xs_mc, legendEntry,"l");

    if(logY_) canvas->SetLogy();
    setHistoStyle(h_xs_data, 1, 1, 1, 0, 0, 20, 1, 1.5);
    setHistoStyle(h_xs_mc, 1, 2, 1, 0, 0, 0, 2, 1.5);
    
    updateHistoAxis(h_xs_data);
    updateHistoAxis(h_xs_mc);
    

    // Draw cross section histograms
    h_xs_data->Draw("E1");
    h_xs_mc->Draw("hist same");
    h_xs_data->Draw("sameE1");
    
    
    // Put additional stuff to histogram
    this->drawCmsLabels(2, 8);
    this->drawDecayChannelLabel(channel);
    legend->Draw("same");
    drawRatioPad(canvas, 0.0, 3.0, h_xs_data);
    
    TH1* h_ratioDataMadgraph = ratioHistogram(h_xs_data, h_xs_mc);
    h_ratioDataMadgraph->Draw("sameE1");

    // Create Directory for Output Plots and write them
    const TString eventFileString = common::assignFolder(outputDir_, channel, systematic);
    canvas->Print(eventFileString+name_+".eps");
    
    return;
    
//     // Prepare additional histograms for root-file
//     TH1* sumMC(0);
//     TH1* sumSignal(0);
//     for(const auto& sampleHistPair : v_sampleHistPair_){
//         if(sampleHistPair.first.sampleType() == Sample::SampleType::ttHbb || sampleHistPair.first.sampleType() == Sample::SampleType::ttHother){
//             if (sumSignal) sumSignal->Add(sampleHistPair.second);
//             else sumSignal = static_cast<TH1*>(sampleHistPair.second->Clone());
//         }
//         if(sampleHistPair.first.sampleType() != Sample::SampleType::data){
//             if(sumMC) sumMC->Add(sampleHistPair.second);
//             else sumMC = static_cast<TH1*>(sampleHistPair.second->Clone());
//         }
//     }
//     if(sumMC) sumMC->SetName(name_);
//     
//     
//     //save Canvas AND sources in a root file
//     TFile out_root(eventFileString+name_+"_source.root", "RECREATE");
//     if(dataHist.second) dataHist.second->Write(name_+"_data");
//     if(sumSignal) sumSignal->Write(name_+"_signalmc");
//     if(sumMC) sumMC->Write(name_+"_allmc");
//     canvas->Write(name_ + "_canvas");
//     out_root.Close();
    
    // Removing created objects
    h_data->Delete();
    h_background->Delete();
    h_signal->Delete();
    h_signal_noWeight->Delete();
    h_signal_gen->Delete();
    h_ttbar_gen->Delete();
    
    h_xs_data->Delete();
    canvas->Clear();
    legend->Clear();
    canvas->Delete();
    legend->Delete();
}




std::map<TString, TH1*> PlotterDiffXS::calculateDiffXS(const TH1* h_data, const TH1* h_signal, const TH1* h_bkg, const TH1* h_signal_noWeight, 
                                                             const TH1* h_signal_gen, const TH1* , const TH1* h_ttbar_gen)const
{
    std::map<TString, TH1*> result;
    
    // Initial values for the cross section calculation
    const ValueError luminosity(luminosity_/1000., 0.);     // Converted to fb-1
    const ValueError br_dilepton(0.06391, 0.00063);
    const ValueError xsection_tt_dilepton(245100.0, 14290.0);
    
    TH1* h_dataMinusBackground = (TH1*)h_data->Clone("h_dataMinusBackground");
    h_dataMinusBackground->Add(h_bkg, -1);
    
    TH1* h_mc = (TH1*)h_signal->Clone("h_signalPlusBackground_reco");
    h_mc->Add(h_bkg, 1.0);
    
    
    ////////////////////////////////////////////////////////////////////// CALCULATING THE INCLUSIVE CROSS SECTION: ttbb
    ValueError N_data_reco(1.,1.);
    N_data_reco.v = h_data->IntegralAndError(0,-1, N_data_reco.e);
    ValueError N_dataMinusBackground(1.,1.);
    N_dataMinusBackground.v = h_dataMinusBackground->IntegralAndError(0,-1, N_dataMinusBackground.e);
    ValueError N_signal_reco(1.,1.);
    N_signal_reco.v = h_signal->IntegralAndError(0,-1, N_signal_reco.e);
    ValueError N_background_reco(1.,1.);
    N_background_reco.v = h_bkg->IntegralAndError(0,-1, N_background_reco.e);
    ValueError N_mc_reco(1.,1.);
    N_mc_reco.v = h_mc->IntegralAndError(0,-1, N_mc_reco.e);
    ValueError N_signal_reco_noWeight(1.,1.);
    N_signal_reco_noWeight.v = h_signal_noWeight->IntegralAndError(0,-1, N_signal_reco_noWeight.e);
    ValueError N_signal_gen(1.,1.);
    N_signal_gen.v = h_signal_gen->IntegralAndError(0,-1, N_signal_gen.e);
    ValueError N_ttbar_gen_events(1.,1.);
    N_ttbar_gen_events.v = h_ttbar_gen->IntegralAndError(0,-1, N_ttbar_gen_events.e);
    
    printf("\n### [Full reco selection, generator weights * PU weights * reco weights * sample SFs (DY, tt+HF)]\n");
    printf("Data:   \t%.3f \t+- %.2f\n", N_data_reco.v, N_data_reco.e);
    printf("Signal: \t%.3f \t+- %.2f\t(RECO tt+bb)\n", N_signal_reco.v, N_signal_reco.e);
    printf("Background: \t%.3f \t+- %.2f\t(RECO rest MC)\n", N_background_reco.v, N_background_reco.e);
    printf("-------------------------------\n");
    printf("MC total: \t%.3f \t+- %.2f\n", N_mc_reco.v, N_mc_reco.e);
    printf("\nData-Bkg: \t%.3f \t+- %.2f\n\n", N_dataMinusBackground.v, N_dataMinusBackground.e);

    printf("### [Full reco selection, generator weights * PU weights * reco weights]\n");
    printf("Signal[reco]: \t%.3f \t+- %.2f\n", N_signal_reco_noWeight.v, N_signal_reco_noWeight.e);
    printf("### [No reco selection, generator weights * PU weights]\n");
    printf("Signal[gen]: \t%.3f \t+- %.2f\n", N_signal_gen.v, N_signal_gen.e);
    printf("-----------------------\n");
    printf("Acceptance: \t%.3f\n\n", N_signal_reco_noWeight.v/N_signal_gen.v);
    printf("tt+jets[gen]: \t%.3f \t+- %.2f\n", N_ttbar_gen_events.v, N_ttbar_gen_events.e);
    printf("-----------------------\n");
    
    ValueError ttbbFraction_Madgraph;
    ttbbFraction_Madgraph.v = N_signal_gen.v/N_ttbar_gen_events.v;
    ttbbFraction_Madgraph.e = uncertaintyBinomial(N_signal_gen.v, N_ttbar_gen_events.v);
    
    printf("ttbb/tt: \t%.3f \t+- %.2f\n\n", ttbbFraction_Madgraph.v, ttbbFraction_Madgraph.e);
    
    printf("Luminosity: \t%.3f \t+- %.2f\n", luminosity.v, luminosity.e);
    printf("BR(tt->ll): \t%.3f \t+- %.2f\n\n", br_dilepton.v, br_dilepton.e);
    
    ValueError xsection_inclusive;
    xsection_inclusive.v = N_dataMinusBackground.v / (br_dilepton.v * luminosity.v * N_signal_reco_noWeight.v / N_signal_gen.v);
    xsection_inclusive.e = N_dataMinusBackground.e / (br_dilepton.v * luminosity.v * N_signal_reco_noWeight.v / N_signal_gen.v);
    ValueError xsection_inclusive_mc;
    xsection_inclusive_mc.v = N_signal_reco.v / (br_dilepton.v * luminosity.v * N_signal_reco_noWeight.v / N_signal_gen.v);
    xsection_inclusive_mc.e = N_signal_reco.e / (br_dilepton.v * luminosity.v * N_signal_reco_noWeight.v / N_signal_gen.v);
    ValueError xsection_inclusive_madgraph;
    xsection_inclusive_madgraph.v = xsection_tt_dilepton.v * N_signal_gen.v / N_ttbar_gen_events.v;
    xsection_inclusive_madgraph.e = xsection_tt_dilepton.e * N_signal_gen.v / N_ttbar_gen_events.v;
    
    printf("###############################\n");
    printf("Inclusive x-seciton(data):     \t%.3f \t+- %.2f\n", xsection_inclusive.v, xsection_inclusive.e);
    printf("Inclusive x-seciton(mc):     \t%.3f \t+- %.2f\n", xsection_inclusive_mc.v, xsection_inclusive_mc.e);
    printf("Inclusive x-seciton(Madgraph): \t%.3f \t+- %.2f\n\n", xsection_inclusive_madgraph.v, xsection_inclusive_madgraph.e);
    
    
    ////////////////////////////////////////////////////////////////////// CALCULATING THE DIFFERENTIAL CROSS SECTION: ttbb
    // Calculating the cross section from data
    const int nBins = h_data->GetNbinsX();
    TH1* h_diffXS_data = (TH1*)h_data->Clone("h_diffXS_data");
    for(int iBin = 1; iBin <= nBins; ++iBin) {
        const ValueError N_signal_reco( h_dataMinusBackground->GetBinContent(iBin), 
                                        h_dataMinusBackground->GetBinError(iBin));
        const ValueError acceptance( h_signal_noWeight->GetBinContent(iBin)/h_signal_gen->GetBinContent(iBin), 
                                     uncertaintyBinomial(h_signal_noWeight->GetBinContent(iBin), h_signal_gen->GetBinContent(iBin)));
        ValueError xsection(1.,1.);
        xsection.v = N_signal_reco.v / (acceptance.v*br_dilepton.v*luminosity.v);
        xsection.e = N_signal_reco.e / (acceptance.v*br_dilepton.v*luminosity.v);
        if(xsection.v != xsection.v) xsection.v = 0.;
        if(xsection.e != xsection.e) xsection.e = 0.;
        // Dividing by the bin width
        xsection.v /= h_diffXS_data->GetBinWidth(iBin);
        xsection.e /= h_diffXS_data->GetBinWidth(iBin);
        
        h_diffXS_data->SetBinContent(iBin, xsection.v);
        h_diffXS_data->SetBinError(iBin, xsection.e);
    }
    result["diffXS_data"] = h_diffXS_data;
    
    // Calculating the cross section from MC
    TH1* h_diffXS_madgraph = (TH1*)h_data->Clone("h_diffXS_madgraph");
    for(int iBin = 1; iBin <= nBins; ++iBin) {
        ValueError ttbbFraction;
        ttbbFraction.v = h_signal_gen->GetBinContent(iBin)/N_ttbar_gen_events.v;
        ttbbFraction.e = uncertaintyBinomial(h_signal_gen->GetBinContent(iBin), N_ttbar_gen_events.v);
        
        ValueError xsection(1.,1.);
        xsection.v = xsection_tt_dilepton.v * ttbbFraction.v;
        xsection.e = sqrt(xsection_tt_dilepton.eOv2() * ttbbFraction.v2() + ttbbFraction.eOv2());
        if(xsection.v != xsection.v) xsection.v = 0.;
        if(xsection.e != xsection.e) xsection.e = 0.;
        // Dividing by the bin width
        xsection.v /= h_diffXS_madgraph->GetBinWidth(iBin);
        xsection.e /= h_diffXS_madgraph->GetBinWidth(iBin);
        
        h_diffXS_madgraph->SetBinContent(iBin, xsection.v);
        h_diffXS_madgraph->SetBinError(iBin, xsection.e);
    }
    // Scaling to the measured inclusive cross section (for shape comparison)
    const double normScale = h_diffXS_data->Integral() / h_diffXS_madgraph->Integral();
    h_diffXS_madgraph->Scale(normScale);
    // Setting the normalisation scale to the underflow bin
    h_diffXS_madgraph->SetBinContent(0, normScale);
    
    result["diffXS_madgraph"] = h_diffXS_madgraph;
    
    
    // Removing created objects
    h_dataMinusBackground->Delete();
    h_mc->Delete();
    
    return result;
}




void PlotterDiffXS::writeResponseMatrix(const Channel::Channel& channel, const Systematic::Systematic& systematic)
{
    // Prepare canvas and legend
    TCanvas* canvas = new TCanvas("","");
    canvas->SetName("");
    canvas->SetTitle("");
    
    
    // Loop over all specified samples and add those with identical legendEntry
    std::vector<LegendHistPair> ttbarHists = legendHistPairsForSamples(sampleTypesTtbar_, v_sampleHistPair_);
    
    // Create histogram corresponding to the sum of all specified samples
    TH1* h_ttbar = sumOfHistos(ttbarHists, "h_ttbar_reco");

    
    // Obtaining the proper histogram to be plotted
    TH1* firstHistToDraw(0);
    if(h_ttbar){
        firstHistToDraw = (TH1*)h_ttbar->Clone();
    }
    else{
        std::cerr<<"ERROR in Plotter! No single sample for drawing exists\n...break\n"<<std::endl;
        exit(237);
    }

    if(logY_) canvas->SetLogy();    
    updateHistoAxis(firstHistToDraw);
    

    // Scaling the content to number of entries for the response matrix
    firstHistToDraw->Scale(firstHistToDraw->GetEntries()/((TH2*)firstHistToDraw)->Integral(0,-1,0,-1));

    // Draw data histogram and stack and error bars
    firstHistToDraw->SetLineColor(0);
    firstHistToDraw->Draw("colz");
    // Moving the pad to the left to fit the color palette
    if(ymax_ != 0) firstHistToDraw->SetMaximum(ymax_);
    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.1);      
    
    
    // Put additional stuff to histogram
    this->drawCmsLabels(2, 8);
    this->drawDecayChannelLabel(channel);

    // Create Directory for Output Plots and write them
    const TString eventFileString = common::assignFolder(outputDir_, channel, systematic);
    canvas->Print(eventFileString+name_+".eps");
    
    // Plotting Purity/Stability if this is a response matrix
    drawPurityStability((TH2*)firstHistToDraw, eventFileString+name_+"_PurStab"+".eps");
    
    return;
    
    // Prepare additional histograms for root-file
    TH1* sumMC(0);
    TH1* sumSignal(0);
    for(const auto& sampleHistPair : v_sampleHistPair_){
        if(sampleHistPair.first.sampleType() == Sample::SampleType::ttHbb || sampleHistPair.first.sampleType() == Sample::SampleType::ttHother){
            if (sumSignal) sumSignal->Add(sampleHistPair.second);
            else sumSignal = static_cast<TH1*>(sampleHistPair.second->Clone());
//             // Do not add Higgs samples to sumMC (all MC samples) in case of overlaid drawing
//             if(drawHiggsOverlaid)continue;
        }
        if(sampleHistPair.first.sampleType() != Sample::SampleType::data){
            if(sumMC) sumMC->Add(sampleHistPair.second);
            else sumMC = static_cast<TH1*>(sampleHistPair.second->Clone());
        }
    }
    if(sumMC) sumMC->SetName(name_);
    
    
    //save Canvas AND sources in a root file
    TFile out_root(eventFileString+name_+"_source.root", "RECREATE");
    if(sumSignal) sumSignal->Write(name_+"_signalmc");
    if(sumMC) sumMC->Write(name_+"_allmc");
    canvas->Write(name_ + "_canvas");
    out_root.Close();
    
    firstHistToDraw->Delete();
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
                                                                                    const std::vector<SampleHistPair> samples)const
{
    std::vector<LegendHistPair> result;
    TH1* tmpHist(0);
    
    for(std::vector<SampleHistPair>::const_iterator i_sampleHistPair = samples.begin();
    i_sampleHistPair != samples.end(); ++i_sampleHistPair)
    {
        const Sample::SampleType& sampleType(i_sampleHistPair->first.sampleType());
        const TString& legendEntry(i_sampleHistPair->first.legendEntry());
        const TH1* hist = i_sampleHistPair->second;

        std::vector<SampleHistPair>::const_iterator incrementIterator(i_sampleHistPair);
        ++incrementIterator;
        const bool lastHist(incrementIterator == samples.end());
        const bool newHist(!tmpHist);

        if(newHist)tmpHist = (TH1*)hist->Clone();

        if(lastHist || (legendEntry!=incrementIterator->first.legendEntry())){
            if(!newHist)tmpHist->Add(hist);
            if(std::find(allowedSampleTypes.begin(), allowedSampleTypes.end(), sampleType) != allowedSampleTypes.end()) 
            {           // Adding histogram to the signal stack
                result.push_back(LegendHistPair(legendEntry, tmpHist));
            }
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


void PlotterDiffXS::drawPurityStability(TH2* histo2d, TString name)const {
    TGraphErrors* g_purity = purityStabilityGraph(histo2d, 0);
    TGraphErrors* g_stability = purityStabilityGraph(histo2d, 1);
    
    
    // Styling axis
    TH1* histo = histo2d->ProjectionX();
    histo->Draw("AXIS");
    histo->SetAxisRange(0., 1.1, "Y");
    histo->GetYaxis()->SetTitle("Values");
    // Styling graphs
    setGraphStyle(g_purity, 21, 1, 1.5, 1, 1);
    setGraphStyle(g_stability, 20, 2, 1.5, 1, 2);
    gPad->SetRightMargin(0.05);
    
    // Drawing graphs
    g_purity->Draw("P0same");
    g_stability->Draw("P0same");
    gPad->Update();
    
    
    // Adding a legend
    TLegend *leg;
    leg = new TLegend ( 0.7,0.84,0.98,0.95,NULL,"brNDC" );
    leg->SetFillColor ( 0 );
    leg->AddEntry(g_purity, "Purity", "p");
    leg->AddEntry(g_stability, "Stability", "p");
    leg->Draw();
    
    // Storing the same canvas with a different name
    gPad->Print(name);
    
    delete g_purity;
    delete g_stability;
}


void PlotterDiffXS::setHistoStyle(TH1* hist, Style_t line, Color_t lineColor, Size_t lineWidth, 
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


void PlotterDiffXS::updateHistoAxis(TH1* histo)const
{
    // Set x and y axis ranges
    if(logY_){
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
        if(logY_){histo->SetMaximum(18.*yMax);}
        else{histo->SetMaximum(1.35*yMax);}
    }
    else{histo->SetMaximum(ymax_);}

    histo->GetXaxis()->SetNoExponent(kTRUE);
    // Set axis titles
    if(XAxis_ != "-") histo->GetYaxis()->SetTitle(YAxis_);
    if(YAxis_ != "-") histo->GetXaxis()->SetTitle(XAxis_);
}


void PlotterDiffXS::drawDecayChannelLabel(const Channel::Channel& channel, const double& textSize)const
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


void PlotterDiffXS::drawCmsLabels(const int cmsprelim, const double& energy, const double& textSize)const
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
    label->AddText(Form(text, samples_.luminosityInInversePb()/1000., energy));
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    if(textSize != 0) label->SetTextSize(textSize);
    label->SetTextAlign(32);
    label->Draw("same");
}


TPaveText* PlotterDiffXS::drawSignificance(TH1* signal, TH1* sigBkg, float min, float max, float yOffset, std::string sLabel, const int type)const
{
    TPaveText *label = new TPaveText();
    if(max<=min) {
        std::cout<<"Wrong range for signal significance in  histogram ("<<name_<<")\n";
        return label;
    }
    
    if(!signal || !sigBkg) {
        std::cout<<"Some input histogram for the signal significance calculation not available\n";
        return label;
    }

    // Finding the bin range corresponding to [min;max]
    int bin1 = signal->FindBin(min);
    int bin2 = signal->FindBin(max);

    // Calculating integral of the signal and background
    float sigInt = signal->Integral(bin1,bin2);
    float sigBkgInt = sigBkg->Integral(bin1,bin2);

    float sigSign = 0.f;
    if(type == 0) sigSign = sigInt/sqrt(sigBkgInt);
    else if(type == 1) sigSign = sigInt/(sigBkgInt - sigInt);
    

    char text[40];
    sprintf(text,"%s = %.2f", sLabel.c_str(), sigSign);
    
    label->SetX1NDC(gStyle->GetPadLeftMargin()+0.4);
    label->SetX2NDC(label->GetX1NDC()+0.1);
    label->SetY2NDC(1.0-gStyle->GetPadTopMargin()-0.13 - yOffset);
    label->SetY1NDC(label->GetY2NDC()-0.05);
    label->SetTextFont(42);
    label->AddText(text);
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->SetTextSize(0.025);
    label->SetTextAlign(32);

    return label;
}


void PlotterDiffXS::setGraphStyle( TGraph* graph, Style_t marker, Color_t markerColor, Size_t markerSize, 
                             Style_t line, Color_t lineColor, Size_t lineWidth)const 
{
    graph->SetMarkerStyle(marker);
    graph->SetMarkerColor(markerColor);
    graph->SetMarkerSize(markerSize);
    graph->SetLineStyle(line);
    graph->SetLineColor(lineColor);
    graph->SetLineWidth(lineWidth);
}


TH1* PlotterDiffXS::drawRatioPad(TPad* pad, const double yMin, const double yMax, TH1* axisHisto, 
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


TH1* PlotterDiffXS::ratioHistogram(const TH1* h_nominator, const TH1* h_denominator)const
{
    TH1* h_ratio = (TH1*)h_nominator->Clone();
    for(int iBin = 1; iBin<=h_ratio->GetNbinsX(); ++iBin) {
        const ValueError nominator(h_nominator->GetBinContent(iBin), h_nominator->GetBinError(iBin));
        const ValueError denominator(h_denominator->GetBinContent(iBin), h_denominator->GetBinError(iBin));
        ValueError ratio;
        ratio.v = nominator.v / denominator.v;
        ratio.e = ratio.v * sqrt(nominator.eOv2() + denominator.eOv2());
        
        h_ratio->SetBinContent(iBin, ratio.v);
        h_ratio->SetBinError(iBin, ratio.e);
    }
    
    return h_ratio;
}