#include <fstream>
#include <iostream>
#include <cstdio>
#include <sstream>
#include <cmath>
#include <iomanip>

#include <TCanvas.h>
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
nameGen_(""),
mode_(0),
rebin_(1),
rangemin_(0),
rangemax_(3),
ymin_(0),
ymax_(0),
sampleTypesData_(std::vector<Sample::SampleType>(0)),
sampleTypesSignal_(std::vector<Sample::SampleType>(0)),
sampleTypesBackground_(std::vector<Sample::SampleType>(0)),
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
        else if(mode_ == 0 && type==Sample::ttbb ) sampleTypesSignal_.push_back(Sample::SampleType(type));
        else if(mode_ == 1 && (type==Sample::ttbb || type==Sample::ttb || type==Sample::tt2b || type==Sample::ttcc || type==Sample::ttother )) 
            sampleTypesSignal_.push_back(Sample::SampleType(type));
        else sampleTypesBackground_.push_back(Sample::SampleType(type));
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
            // Getting generator level histograms
            const std::vector<double>& v_weight_noScale_ee(globalWeights_noScale.at(systematic).at(Channel::ee));
            const std::vector<Sample>& v_sample_ee = m_systematicChannelSample.at(systematic).at(Channel::ee);
            if(!this->prepareDataset(v_sample_ee, v_weight_noScale_ee, v_sampleHistPairGen_, nameGen_)) {
                std::cout<<"WARNING! Cannot find generator level histograms for all datasets, for (channel/systematic): "
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
    TLegend* legend = new TLegend(0.70,0.55,0.92,0.85);
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
    for(auto sampleHistPair : v_sampleHistPair_){
        TH1* tmp_hist = sampleHistPair.second;
        if(rebin_>1) tmp_hist->Rebin(rebin_);
        setStyle(sampleHistPair);
    }
    
    
    // Generic information about the histogram being plotted
    TString histo_name = v_sampleHistPair_.size() > 0 ? TString(v_sampleHistPair_.at(0).second->GetName()) : "";
    
    
    // Loop over all samples and add those with identical legendEntry
    // And sort them into the categories data, signal, background
    
    // Create histogram corresponding to the sum of all stacked signal histograms
    std::vector<LegendHistPair> dataHists = legendHistPairsForSamples(sampleTypesData_, v_sampleHistPair_);
    TH1* h_data = sumOfHistos(dataHists);
    
    // Create histogram corresponding to the sum of all stacked signal histograms
    std::vector<LegendHistPair> signalHists = legendHistPairsForSamples(sampleTypesSignal_, v_sampleHistPair_);
    TH1* h_signal = sumOfHistos(signalHists);
    
    // Create histogram corresponding to the sum of all stacked background histograms
    std::vector<LegendHistPair> backgroundHists = legendHistPairsForSamples(sampleTypesBackground_, v_sampleHistPair_);
    TH1* h_background = sumOfHistos(backgroundHists);
    
    // Create histogram corresponding to the sum of all stacked signal histograms without sample level weights
    std::vector<LegendHistPair> signalHists_noWeight = legendHistPairsForSamples(sampleTypesSignal_, v_sampleHistPair_noWeight_);
    TH1* h_signal_noWeight = sumOfHistos(signalHists_noWeight);
    
        
    // Get all signal samples at generator level
    std::vector<LegendHistPair> signalHistsGen = legendHistPairsForSamples(sampleTypesSignal_, v_sampleHistPairGen_);
    TH1* h_signal_gen = sumOfHistos(signalHistsGen);
    
    printf("N histos Gen: %d   N events: %.1f   Integral: %.3f\n", (int)signalHistsGen.size(), h_signal_gen->GetEntries(), h_signal_gen->Integral());
    
    // Calculating actual differential cross section and getting corresponding histograms
    LegendHistPair xsectionSet_data = calculateDiffXS(h_data, h_signal, h_background, h_signal_noWeight, h_signal_gen);
    
    return;

    
//     // Add entries to legend
//     if(dataHist.second) legend->AddEntry(dataHist.second, dataHist.first,"pe");

    
    // Obtaining the proper histogram to be plotted
    TH1* firstHistToDraw(0);
    if(!firstHistToDraw)
    {
        std::cerr<<"ERROR in Plotter! Histogram for drawing doesn't exist\n...break\n"<<std::endl;
        exit(237);
    }
    
//     // Set x and y axis ranges
//     if(logY_){
//       // Setting minimum to >0 value
//       // FIXME: Should we automatically calculate minimum value instead of the fixed value?
//       firstHistToDraw->SetMinimum(1e-1);
//       if(ymin_>0) firstHistToDraw->SetMinimum(ymin_);
//       canvas->SetLogy();
//     }
//     else firstHistToDraw->SetMinimum(ymin_);
// 
//     if(rangemin_!=0. || rangemax_!=0.) {firstHistToDraw->SetAxisRange(rangemin_, rangemax_, "X");}
//     
//     if(ymax_==0.){
//         // Determining the highest Y value that is plotted
//         float yMax = stacksumSignal ? stacksumSignal->GetBinContent(stacksumSignal->GetMaximumBin()) : 0.f;
//         
//         // Scaling the Y axis
//         if(logY_){firstHistToDraw->SetMaximum(18.*yMax);}
//         else{firstHistToDraw->SetMaximum(1.35*yMax);}
//     }
//     else{firstHistToDraw->SetMaximum(ymax_);}
// 
//     firstHistToDraw->GetXaxis()->SetNoExponent(kTRUE);
// 
// 
//     // Draw data histogram and stack and error bars
//     firstHistToDraw->SetLineColor(0);
//     firstHistToDraw->Draw();
//     
//     // Put additional stuff to histogram
//     this->drawCmsLabels(2, 8);
//     this->drawDecayChannelLabel(channel);
//     if(mode_!=1) legend->Draw("SAME");
//     if(mode_==0){
//         common::drawRatio(dataHist.second, stacksumSignal, 0, 0.5, 1.7);
//         firstHistToDraw->GetXaxis()->SetLabelSize(0);
//         firstHistToDraw->GetXaxis()->SetTitleSize(0);
//     }

    // Create Directory for Output Plots and write them
    const TString eventFileString = common::assignFolder(outputDir_, channel, systematic);
    canvas->Print(eventFileString+name_+".eps");
    
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
    
    firstHistToDraw->Delete();
    canvas->Clear();
    legend->Clear();
    delete canvas;
    delete legend;
}




void PlotterDiffXS::writeResponseMatrix(const Channel::Channel& channel, const Systematic::Systematic& systematic)
{
    // Prepare canvas and legend
    TCanvas* canvas = new TCanvas("","");
    canvas->SetName("");
    canvas->SetTitle("");
    
    
    // Loop over all samples and add those with identical legendEntry
    // And sort them into the categories data, Higgs, other
    std::vector<LegendHistPair> signalHists;
    TH1* tmpHist(0);
    for(std::vector<SampleHistPair>::iterator i_sampleHistPair = v_sampleHistPair_.begin();
        i_sampleHistPair != v_sampleHistPair_.end(); ++i_sampleHistPair)
    {
        const Sample::SampleType& sampleType(i_sampleHistPair->first.sampleType());
        const TString& legendEntry(i_sampleHistPair->first.legendEntry());
        const TH1* hist = i_sampleHistPair->second;

        if(sampleType == Sample::data) continue;
        
        std::vector<SampleHistPair>::iterator incrementIterator(i_sampleHistPair);
        ++incrementIterator;
        const bool lastHist(incrementIterator == v_sampleHistPair_.end());
        const bool newHist(!tmpHist);

        if(newHist)tmpHist = (TH1*)hist->Clone();

        if(lastHist || (legendEntry!=incrementIterator->first.legendEntry())){
            if(!newHist)tmpHist->Add(hist);
            if(std::find(sampleTypesSignal_.begin(), sampleTypesSignal_.end(), sampleType) != sampleTypesSignal_.end()) 
            {           // Adding histogram to the signal stack
                signalHists.push_back(LegendHistPair(legendEntry, tmpHist));
            }
            tmpHist = 0;
            
        }
        else{
            if(!newHist)tmpHist->Add(hist);
        }
    }
    
    // Create histogram corresponding to the sum of all stacked signal histograms
    TH1* stacksumSignal = sumOfHistos(signalHists);

    
    // Obtaining the proper histogram to be plotted
    TH1* firstHistToDraw(0);
    if(stacksumSignal){
        firstHistToDraw = (TH1*)stacksumSignal->Clone();
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



void PlotterDiffXS::setStyle(SampleHistPair& sampleHistPair)
{
    TH1* hist(sampleHistPair.second);

    hist->SetFillColor(sampleHistPair.first.color());
    hist->SetLineColor(sampleHistPair.first.color());
    hist->SetLineWidth(1);
    
    if(XAxis_ == "-") XAxis_ = hist->GetXaxis()->GetTitle();
    if(YAxis_ == "-") YAxis_ = hist->GetYaxis()->GetTitle();

    if(sampleHistPair.first.sampleType() == Sample::SampleType::data){
        hist->SetFillColor(0);
        hist->SetMarkerStyle(20);
        hist->SetMarkerSize(1.);
        hist->SetLineWidth(1);
        hist->GetXaxis()->SetLabelFont(42);
        hist->GetYaxis()->SetLabelFont(42);
        hist->GetXaxis()->SetTitleFont(42);
        hist->GetYaxis()->SetTitleFont(42);
        hist->GetYaxis()->SetTitleOffset(1.7);
        hist->GetXaxis()->SetTitleOffset(1.25);
    }
    hist->GetYaxis()->SetTitle(YAxis_);
    hist->GetXaxis()->SetTitle(XAxis_);
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
        float yMax = histo ? histo->GetBinContent(histo->GetMaximumBin()) : 0.f;
        
        // Scaling the Y axis
        if(logY_){histo->SetMaximum(18.*yMax);}
        else{histo->SetMaximum(1.35*yMax);}
    }
    else{histo->SetMaximum(ymax_);}

    histo->GetXaxis()->SetNoExponent(kTRUE);
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


TH1* PlotterDiffXS::sumOfHistos(const std::vector<LegendHistPair> histos)const
{
    TH1* sum(0);
    for(size_t i = 0; i < histos.size(); ++i){
        if(i==0) sum = (TH1*)histos.at(0).second->Clone();
        else sum->Add((TH1*)histos.at(i).second);
    }
    
    return sum;
}


PlotterDiffXS::LegendHistPair PlotterDiffXS::calculateDiffXS(TH1* h_data, TH1* h_signal, TH1* h_bkg, TH1* h_signal_noWeight, TH1* h_signal_gen)const
{
    LegendHistPair result;
    
    // Initial values for the cross section calculation
    ValueError luminosity(luminosity_, 0.);
    ValueError br_dilepton(0.6391, 0.00063);
    
    TH1* h_dataMinusBackground = (TH1*)h_data->Clone("h_dataMinusBackground");
    h_dataMinusBackground->Add(h_bkg, -1);
    
    TH1* h_mc = (TH1*)h_signal->Clone("mc");
    h_mc->Add(h_bkg, 1.0);
    
    printf("\n[Full reco selection, generator weights * PU weights * reco weights * sample SFs (DY, tt+HF)]\n");
    printf("Data:   \t%.3f\n", h_data->Integral());
    printf("Signal: \t%.3f\t(RECO tt+bb)\n", h_signal->Integral());
    printf("Background: \t%.3f\t(RECO rest MC)\n", h_bkg->Integral());
    printf("-------------------------------\n");
    printf("MC total: \t%.3f\n", h_mc->Integral());
    printf("\nData-Bkg: \t%.3f\n", h_dataMinusBackground->Integral());

    printf("\n[Full reco selection, generator weights * PU weights * reco weights]\n");
    printf("Signal[reco]: \t%.3f\n", h_signal_noWeight->Integral());
    printf("[No reco selection, generator weights * PU weights]\n");
    printf("Signal[gen]: \t%.3f\n", h_signal_gen->Integral());
    printf("-----------------------\n");
    printf("Acceptance: \t%.3f\n", h_signal_noWeight->Integral()/h_signal_gen->Integral());
    
    return result;
}


TGraphErrors* PlotterDiffXS::purityStabilityGraph(TH2* h2d, const int type)const
{

    const int nBins = h2d->GetNbinsX();

    TGraphErrors* graph = new TGraphErrors(nBins);

    // Calculating each point of graph for each diagonal bin
    for(int iBin = 1; iBin<=nBins; ++iBin) {
        const double diag = h2d->GetBinContent(iBin, iBin);
        const double reco = h2d->Integral(iBin, iBin, 0, -1)+1e-30;
        const double gen = h2d->Integral(0, -1, iBin, iBin)+1e-30;
        
        const double value = type == 0 ? diag/reco : diag/gen;
        const double error = type == 0 ? uncertaintyBinomial(diag, reco) : uncertaintyBinomial(diag, gen);

        const double bin = h2d->GetXaxis()->GetBinCenter(iBin);
        const double binW = h2d->GetXaxis()->GetBinWidth(iBin);
        
        graph->SetPoint(iBin-1, bin, value);
        graph->SetPointError(iBin-1, binW/2., error);
    }

    return graph;

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


double PlotterDiffXS::uncertaintyBinomial(const double pass, const double all)const 
{
    return (1./all)*sqrt(pass - pass*pass/all);
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
