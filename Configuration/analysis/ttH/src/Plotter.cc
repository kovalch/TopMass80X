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
#include <TH1D.h>
#include <TPaveText.h>
#include <TError.h>

#include "Plotter.h"
#include "higgsUtils.h"
#include "Samples.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"
#include "../../common/include/plotterUtils.h"





Plotter::Plotter(const char* outputDir,
                 const Samples& samples,
                 const DrawMode::DrawMode& drawMode):
outputDir_(outputDir),
topLeftLabelId_(1),
samples_(samples),
drawMode_(drawMode),
fileReader_(RootFileReader::getInstance()),
hasPseudodata_(false),
name_("defaultName"),
bins_(0),
rebin_(1),
rangemin_(0),
rangemax_(3),
ymin_(0),
ymax_(0),
YAxis_(""),
XAxis_(""),
logX_(false),
logY_(false),
normalizeToData_(false)
{
    // Suppress default info that canvas is printed
    gErrorIgnoreLevel = 1001;
}



void Plotter::setOptions(const TString& name, const TString&,
                         const TString& YAxis, const TString& XAxis,
                         const int rebin, const bool normalizeToData,
                         const bool logX, const bool logY,
                         const double& ymin, const double& ymax,
                         const double& rangemin, const double& rangemax,
                         const int bins,
                         const std::vector<double>& XAxisbins, const std::vector<double>& XAxisbinCenters)
{
    name_ = name; //Variable name you want to plot
    YAxis_ = YAxis; //Y-axis title
    XAxis_ = XAxis; //X-axis title
    rebin_ = rebin; //Nr. of bins to be merged together
    logX_ = logX; //Draw X-axis in Log scale
    logY_ = logY; //Draw Y-axis in Log scale
    ymin_ = ymin; //Min. value in Y-axis
    ymax_ = ymax; //Max. value in Y-axis
    rangemin_ = rangemin; //Min. value in X-axis
    rangemax_ = rangemax; //Max. value in X-axis
    bins_ = bins; //Number of bins to plot
    XAxisbins_.clear();
    XAxisbins_ = XAxisbins; // Bins edges=bins+1
    XAxisbinCenters_.clear();
    XAxisbinCenters_ = XAxisbinCenters; //Central point for BinCenterCorrection=bins
    normalizeToData_ = normalizeToData;	// Whether stacksum should be normalized to data
}



void Plotter::producePlots()
{
    //std::cout<<"--- Beginning of plot production\n\n";
    
    // Set style
    common::setHHStyle(*gStyle, false);
    
    // Access correction factors
    const SystematicChannelFactors globalWeights = this->scaleFactors();
    
    // Loop over all channels and systematics and produce plots
    const SystematicChannelSamples& m_systematicChannelSample(samples_.getSystematicChannelSamples());
    for(const auto& systematicChannelSamples : m_systematicChannelSample){
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        for(const auto& channelSample : systematicChannelSamples.second){
            const Channel::Channel& channel(channelSample.first);
            const std::vector<Sample>& v_sample(channelSample.second);
            if(v_sample.size() < 1) continue;
            const std::vector<double>& v_weight(globalWeights.at(systematic).at(channel));
            if(!this->prepareDataset(v_sample, v_weight)){
                std::cout<<"WARNING! Cannot find histograms for all datasets, for (channel/systematic): "
                         << Channel::convert(channel) << "/" << systematic.name()
                         <<"\n... skip this plot\n";
                return;
            }
            this->write(channel, systematic);
        }
    }
    
    //std::cout<<"\n=== Finishing of plot production\n\n";
}



const SystematicChannelFactors& Plotter::scaleFactors()
{
    const TString fullStepname = tth::extractSelectionStepAndJetCategory(name_);
    
    // Check if map contains already scale factors for this step, else access and fill them
    if(m_stepFactors_.find(fullStepname) == m_stepFactors_.end()){
        m_stepFactors_[fullStepname] = samples_.globalWeights(fullStepname).first;
    }
    
    return m_stepFactors_.at(fullStepname);
}



bool Plotter::prepareDataset(const std::vector<Sample>& v_sample,
                             const std::vector<double>& v_weight)
{
    bool allHistosAvailable(true);
    
    // Associate histogram to dataset if histogram can be found
    v_sampleHistPair_.clear();
    TH1::AddDirectory(kFALSE);
    for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
        const auto& sample(v_sample.at(iSample));
        SampleHistPair p_sampleHist;
        TH1D* hist = fileReader_->GetClone<TH1D>(sample.inputFile(), name_, true, false);
        if(!hist){
            // Print message only for one histogram
            if(allHistosAvailable)
                std::cout<<"Histogram ("<<name_<<") not found e.g. in file ("<<sample.inputFile()<<")\n";
            p_sampleHist = SampleHistPair(sample, 0);
            allHistosAvailable = false;
        } else {
            if(hist->GetSumw2N()<1) hist->Sumw2();
            // Apply weights
            if(sample.sampleType() != Sample::data){
                const double& weight = v_weight.at(iSample);
                hist->Scale(weight);
            }

            // Clone histogram directly here
            TH1D* histClone = (TH1D*) hist->Clone();
            p_sampleHist = SampleHistPair(sample, histClone);
        }
        if(sample.sampleType() == Sample::pseudodata) hasPseudodata_ = true;
        v_sampleHistPair_.push_back(p_sampleHist);
    }
    //std::cout<<"Number of samples used for histogram ("<<name_<<"): "<<v_sampleHistPair_.size()<<"\n";


    return allHistosAvailable;
}



void Plotter::write(const Channel::Channel& channel, const Systematic::Systematic& systematic)
{
    // Prepare canvas and legend
    TCanvas* canvas = new TCanvas("","");
    canvas->Clear();
    
    
    // Here fill colors and line width are adjusted, and potentially rebinning applied
    for(auto& sampleHistPair : v_sampleHistPair_){
        sampleHistPair.second = common::rebinnedHistoInRange(sampleHistPair.second, rebin_, rangemin_, rangemax_);
        this->setStyle(sampleHistPair, true);
    }
    
    
    // Check whetherHiggs sample should be drawn overlaid and/or scaled
    bool drawHiggsOverlaid(false);
    bool drawHiggsScaled(false);
    if(drawMode_ == DrawMode::overlaid){drawHiggsOverlaid = true;}
    else if(drawMode_ == DrawMode::scaledoverlaid){drawHiggsOverlaid = true; drawHiggsScaled = true;}
    
    
    // Specific histograms for ttbb and ttHbb to calculate signal significances
    TH1* ttbbHist(0);
    TH1* ttHbbHist(0);
    
    
    // Loop over all samples and add those with identical legendEntry
    // And sort them into the categories data, Higgs, other
    LegendHistPair dataHist;
    Sample dataSample;
    
    std::vector<LegendHistPair> higgsHists;
    std::vector<LegendHistPair> stackHists;
    TH1D* tmpHist(0);
    for(std::vector<SampleHistPair>::iterator i_sampleHistPair = v_sampleHistPair_.begin();
        i_sampleHistPair != v_sampleHistPair_.end(); ++i_sampleHistPair)
    {
        const Sample& sample = i_sampleHistPair->first;
        const Sample::SampleType& sampleType(sample.sampleType());
        const TString& legendEntry(sample.legendEntry());
        const TH1* hist = i_sampleHistPair->second;
        
        std::vector<SampleHistPair>::iterator incrementIterator(i_sampleHistPair);
        ++incrementIterator;
        const bool lastHist(incrementIterator == v_sampleHistPair_.end());
        const bool newHist(!tmpHist);
        
        if(newHist) tmpHist = (TH1D*)hist->Clone();
        
        if(lastHist || (legendEntry!=incrementIterator->first.legendEntry())){
            if(!newHist) tmpHist->Add(hist);
            if(sampleType == Sample::data || sampleType == Sample::pseudodata) dataHist = LegendHistPair(legendEntry, tmpHist);
            else if(sampleType == Sample::ttHbb) higgsHists.push_back(LegendHistPair(legendEntry, tmpHist));
            else if(sampleType == Sample::ttHother) higgsHists.push_back(LegendHistPair(legendEntry, tmpHist));
            else stackHists.push_back(LegendHistPair(legendEntry, tmpHist));
            if(sampleType == Sample::ttHbb) ttHbbHist = tmpHist;
            else if(sampleType == Sample::ttbb) ttbbHist = tmpHist;
            tmpHist = 0;
        }
        else{
            if(!newHist) tmpHist->Add(hist);
        }
    }
    
    
    // Create histogram corresponding to the sum of all stacked histograms
    TH1D* stacksum(0);
    for(size_t iHist = 0; iHist < stackHists.size(); ++iHist){
        if(iHist == 0) stacksum = (TH1D*)stackHists.at(0).second->Clone();
        else stacksum->Add((TH1D*)stackHists.at(iHist).second);
    }
    if(!drawHiggsOverlaid){
        bool stackDefined(stacksum);
        for(auto higgsHist : higgsHists){
            if(!stackDefined){
                stacksum = (TH1D*)higgsHist.second->Clone();
                stackDefined = true;
            }
            else{
                stacksum->Add((TH1D*)higgsHist.second);
            }
        }
    }

    // Normalizing MC to data if required
    if(normalizeToData_) {
        const double normScale = common::normalize(stacksum, dataHist.second->Integral());
        for(auto stackHisto : stackHists) stackHisto.second->Scale(normScale);
        if(drawHiggsOverlaid) for(auto higgsHisto : higgsHists) higgsHisto.second->Scale(normScale);
        ttbbHist->Scale(normScale);
        ttHbbHist->Scale(normScale);
    }
    
    
    // Draw signal significance for dijet_mass H->bb
    std::vector<TPaveText*> significanceLabels;
    const TString histoName(stackHists.at(0).second->GetName());
    if(!drawHiggsOverlaid && stackHists.size() && histoName.Contains("dijet_dijet_mass")){
//         significanceLabels.push_back(this->drawSignificance(ttbbHist, stacksum, 100.1, 139.9, 0.f, "#frac{S_{ttbb}}{#sqrt{S_{ttbb}+B}}", 0));
        significanceLabels.push_back(this->drawSignificance(ttHbbHist, stacksum, 100.1, 139.9, 0.1f, "#frac{N_{t#bar{t}H}}{#sqrt{N_{t#bar{t}H}+N_{bkg}}}", 0));
        significanceLabels.push_back(this->drawSignificance(ttHbbHist, stacksum, 100.1, 139.9, 0.2f, "#frac{N_{t#bar{t}H}}{N_{bkg}}", 1));
    }
    
    // Blinding dijet mass plots around the Higgs mass region
    if(( (TString(stackHists.at(0).second->GetName())).Contains("dijet_mass") || (TString(stackHists.at(0).second->GetName())).Contains("Mjj") ) 
        && dataHist.second && !hasPseudodata_) {
        int bin1 = dataHist.second->FindBin(100.1);
        int bin2 = dataHist.second->FindBin(139.9);
        while(bin1 <= bin2) {
            dataHist.second->SetBinContent(bin1, 1e-10);
            dataHist.second->SetBinError(bin1, 0);
            bin1++;
        }
    }
    
    
    // If Higgs signal scaled: scale sample and modify legend entry
    if(drawHiggsScaled){
        for(auto& higgsHist : higgsHists){
            std::stringstream ss_scaleFactor;
            const double signalScaleFactor = stacksum ? stacksum->Integral()/higgsHist.second->Integral() : 1.;
            if(isnan(signalScaleFactor) || isinf(signalScaleFactor)) continue;
            const int precision(signalScaleFactor>=100 ? 0 : 1);
            ss_scaleFactor<<" x"<<std::fixed<<std::setprecision(precision)<<signalScaleFactor;
            higgsHist.second->Scale(signalScaleFactor);
            higgsHist.first.Append(ss_scaleFactor.str());
        }
    }
    
    
    // If Higgs samples should be stacked, add them to stack histograms and clear Higgs vector
    if(!drawHiggsOverlaid){
        stackHists.insert(stackHists.end(), higgsHists.begin(), higgsHists.end());
        higgsHists.clear();
    }
    
    // Creating a legend
    const int nLegendEntries = stackHists.size() + higgsHists.size();
    TLegend* legend = common::createLegend(0.62, 0.55, 2, nLegendEntries, 0.03);
    
    // Create the stack and add entries to legend
    THStack* stack(0);
    if(dataHist.second) legend->AddEntry(dataHist.second, dataHist.first,"pe");
    if(stackHists.size()) stack = new THStack("def", "def");
    for(std::vector<LegendHistPair>::reverse_iterator higgsHist = higgsHists.rbegin(); higgsHist!=higgsHists.rend(); ++higgsHist){
        higgsHist->second->SetFillStyle(0);
        higgsHist->second->SetLineWidth(2);
        legend->AddEntry(higgsHist->second, higgsHist->first,"l");
    }
    for(std::vector<LegendHistPair>::reverse_iterator stackHist = stackHists.rbegin(); stackHist!=stackHists.rend(); ++stackHist){
        stackHist->second->SetLineColor(1);
        legend->AddEntry(stackHist->second, stackHist->first,"f");
    }
    for(std::vector<LegendHistPair>::const_iterator stackHist = stackHists.begin(); stackHist!=stackHists.end(); ++stackHist){
        stack->Add(stackHist->second);
    }

    
    // Preparing a statistical uncertainty band for MC stack
    TH1* stacksum_statBand(0);
    if(stacksum){
        stacksum_statBand = (TH1*)stacksum->Clone();
        common::setHistoStyle(stacksum_statBand, 1,1,1, -1,-1,-1, 3354, 1);
    }
    
    
    // Set x and y axis
    TH1* firstHistToDraw(0);
    if(dataHist.second) firstHistToDraw = (TH1*)dataHist.second->Clone();
    else if(stacksum){
        firstHistToDraw = (TH1*)stacksum->Clone();
    }
    else if(higgsHists.size()){
        firstHistToDraw = (TH1*)higgsHists.at(0).second->Clone();
    }
    else{
        std::cerr<<"ERROR in Plotter::write()! No single sample for drawing exists\n...break\n"<<std::endl;
        exit(237);
    }

    firstHistToDraw->GetXaxis()->SetNoExponent(kTRUE);

    
    // Add binwidth including unit to y-axis
    common::addBinDivisionToYaxis(firstHistToDraw, {"_multiplicity_", "jetCategories"});
    
    
    // Draw data histogram and stack
    firstHistToDraw->SetLineColor(0);
    firstHistToDraw->Draw();
    if(dataHist.second) dataHist.second->Draw("same E1 X0");
    if(stack) stack->Draw("same HIST");
    gPad->RedrawAxis();
    stacksum_statBand->Draw("same E2");  // error bars for stack (which, stat or combined with syst ?)
    if(dataHist.second) dataHist.second->Draw("same E1 X0");
    for(const auto& higgsHist : higgsHists){
        higgsHist.second->Draw("same HIST");
    }
    // Adding the uncertainty band to the legend
    const int nColumns = legend->GetNColumns();
    const int nEntries = legend->GetListOfPrimitives()->GetEntries();
    const int nEntriesToInsert = nColumns - nEntries%nColumns - 1;
    // Adding empty entries to ensure that the uncertainty is in the rightmost column
    for(int iEn = 0; iEn < nEntriesToInsert; ++iEn) legend->AddEntry((TObject*)0, "", "");
    legend->AddEntry(stacksum_statBand, "Uncertainty", "F");
    
    // Updating the histo axis
    TH1* axisHisto = common::updatePadYAxisRange(canvas, logY_, 0.35);
    
    // Applying the configured Y axis range
    if(ymax_ != 0.) axisHisto->SetMaximum(ymax_);
    if(ymin_ != 0.) axisHisto->SetMinimum(ymin_);
    axisHisto->GetXaxis()->SetNoExponent(kTRUE);

    // Redrawing the axis
    canvas->RedrawAxis();
    
    // Put additional stuff to histogram
    common::drawCmsLabels(0, 8, samples_.luminosityInInversePb()/1000.);
//     if(topLeftLabelId_==1) this->drawTopLeftLabel(Channel::label(channel));
//     else if(topLeftLabelId_==2) this->drawTopLeftLabel(Systematic::convertType(systematic.type())+=Systematic::convertVariation(systematic.variation()));
    for(TPaveText* label : significanceLabels) if(label) label->Draw("same");
    legend->Draw("SAME");
    if(dataHist.second && stacksum){
        common::drawRatioPad(canvas, 0.5, 1.7, "#frac{Data}{MC}");
        
        TH1* ratio_histo = common::ratioHistogram(dataHist.second, stacksum, 1);
        common::setHistoStyle(ratio_histo, -1,-1,2);
        ratio_histo->Draw("same E");
        
        TH1* ratio_statBand = common::ratioHistogram(stacksum, stacksum, 1);
        common::setHistoStyle(ratio_statBand, -1,-1,-1, -1,-1,-1, 3354,1);
        ratio_statBand->Draw("same E2");
    }

    // Create Directory for Output Plots and write them
    const TString eventFileString = common::assignFolder(outputDir_, channel, systematic);
    canvas->Print(eventFileString+name_+".eps");
    
    
    // Prepare additional histograms for root-file
    TH1* sumMC(0);
    TH1* sumSignal(0);
    TH1* sum_ttbb(0);
    TH1* sum_ttb(0);
    TH1* sum_tt2b(0);
    TH1* sum_ttOther(0);
    TH1* sum_bkg(0);
    for(const auto& sampleHistPair : v_sampleHistPair_){
        if(sampleHistPair.first.sampleType() == Sample::SampleType::ttHbb || sampleHistPair.first.sampleType() == Sample::SampleType::ttHother){
            sumSignal = addOrCreateHisto(sumSignal, sampleHistPair.second);
            // Do not add Higgs samples to sumMC (all MC samples) in case of overlaid drawing
            if(drawHiggsOverlaid) continue;
        }
        if(sampleHistPair.first.sampleType() != Sample::SampleType::data){
            sumMC = addOrCreateHisto(sumMC, sampleHistPair.second);
        }
        if(sampleHistPair.first.sampleType() == Sample::SampleType::ttbb) {
            sum_ttbb = addOrCreateHisto(sum_ttbb, sampleHistPair.second);
        } else if(sampleHistPair.first.sampleType() == Sample::SampleType::ttb) {
            sum_ttb = addOrCreateHisto(sum_ttb, sampleHistPair.second);
        } else if(sampleHistPair.first.sampleType() == Sample::SampleType::tt2b) {
            sum_tt2b = addOrCreateHisto(sum_tt2b, sampleHistPair.second);
        } else if(sampleHistPair.first.sampleType() == Sample::SampleType::ttcc || sampleHistPair.first.sampleType() == Sample::SampleType::ttother) {
            sum_ttOther = addOrCreateHisto(sum_ttOther, sampleHistPair.second);
        } else {
            sum_bkg = addOrCreateHisto(sum_bkg, sampleHistPair.second);
        }
            
    }
    if(sumMC) sumMC->SetName(name_);
    
    
    // Save Canvas AND sources in a root file
    TFile out_root(eventFileString+name_+"_source.root", "RECREATE");
    if(dataHist.second) dataHist.second->Write(name_+"_data");
    if(sumSignal) sumSignal->Write(name_+"_signalmc");
    if(sumMC) sumMC->Write(name_+"_allmc");
    if(sum_ttbb) sum_ttbb->Write(name_+"_ttbb");
    if(sum_ttb) sum_ttb->Write(name_+"_ttb");
    if(sum_tt2b) sum_tt2b->Write(name_+"_tt2b");
    if(sum_ttOther) sum_ttOther->Write(name_+"_ttOther");
    if(sum_bkg) sum_bkg->Write(name_+"_bkg");
    canvas->Write(name_ + "_canvas");
    out_root.Close();
    
    firstHistToDraw->Delete();
    for(TPaveText* label : significanceLabels) if(label) label->Delete();
    canvas->Clear();
    legend->Clear();
    delete canvas;
    delete legend;
}



void Plotter::setStyle(SampleHistPair& sampleHistPair, const bool isControlPlot)
{
    TH1* hist(sampleHistPair.second);
    const Sample& sample = sampleHistPair.first;

    hist->SetFillColor(sample.color());
    hist->SetLineColor(sample.color());
    hist->SetLineWidth(1);

    if(XAxis_ == "-") XAxis_ = hist->GetXaxis()->GetTitle();
    if(YAxis_ == "-") YAxis_ = hist->GetYaxis()->GetTitle();

    if(sample.sampleType() == Sample::data || sample.sampleType() == Sample::pseudodata){
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
        if ((name_.Contains("pT") || name_.Contains("Mass")) && !name_.Contains("Rapidity")) {
            hist->GetXaxis()->SetTitle(XAxis_+" #left[GeV#right]");
            hist->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d"+XAxis_+"}"+" #left[GeV^{-1}#right]");
        } else {
            hist->GetXaxis()->SetTitle(XAxis_);
            hist->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d"+XAxis_+"}");
        }
        if (isControlPlot) hist->GetYaxis()->SetTitle(YAxis_);
    }
}



void Plotter::drawTopLeftLabel(const TString& text, const double& textSize)const
{
    TPaveText* decayChannel = new TPaveText();

    decayChannel->AddText(text);

    decayChannel->SetX1NDC(      gStyle->GetPadLeftMargin());
    decayChannel->SetY1NDC(1.0 - gStyle->GetPadTopMargin() );
    decayChannel->SetX2NDC(      gStyle->GetPadLeftMargin() + 0.15 );
    decayChannel->SetY2NDC(1.0 );

    decayChannel->SetFillStyle(0);
    decayChannel->SetBorderSize(0);
    if (textSize!=0) decayChannel->SetTextSize(textSize);
    decayChannel->SetTextAlign(12);
    decayChannel->Draw("same");
}



TPaveText* Plotter::drawSignificance(const TH1* const signal, const TH1* const signalPlusBackground,
                                     const float xMin, const float xMax,
                                     const float yOffset, const std::string& sLabel, const int type)const
{
    if(xMax <= xMin){
        std::cout<<"Wrong range for signal significance in  histogram ("<<name_<<")\n";
        return 0;
    }
    
    if(!signal || !signalPlusBackground){
        std::cout<<"Some input histogram for the signal significance calculation not available\n";
        return 0;
    }
    
    // Find the bin range corresponding to [min;max]
    const int bin1 = signal->FindFixBin(xMin);
    const int bin2 = signal->FindFixBin(xMax);
    
    // Calculating integral of the signal and background
    const float sigInt = signal->Integral(bin1, bin2);
    Double_t sigBkgIntErr(0.0);
    float sigBkgInt = signalPlusBackground->IntegralAndError(bin1, bin2, sigBkgIntErr);
    // Varying the background maximally up to have the lower limit on the ratio
    sigBkgInt += sigBkgIntErr;
    
    float sigSign = 0.f;
    if(type == 0) sigSign = sigInt/sqrt(sigBkgInt);
    else if(type == 1) sigSign = sigInt/(sigBkgInt - sigInt);
    
    char text[40];
    sprintf(text,"%s = %.2f", sLabel.c_str(), sigSign);
    
    TPaveText* label = new TPaveText();
    label->SetX1NDC(gStyle->GetPadLeftMargin()+0.55);
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



TH1* Plotter::addOrCreateHisto(TH1* base, TH1* add_histo) 
{
    if(!base) base = static_cast<TH1*>(add_histo->Clone());
    else base->Add(add_histo);
    
    return base;
}

