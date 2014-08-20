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

#include "Plotter.h"
#include "higgsUtils.h"
#include "Samples.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"
#include "../../common/include/plotterUtils.h"





/// Folder for Control Plots output
constexpr const char* ControlPlotDIR = "Plots";





Plotter::Plotter(const Samples& samples,
                 const DrawMode::DrawMode& drawMode):
samples_(samples),
drawMode_(drawMode),
fileReader_(RootFileReader::getInstance()),
name_("defaultName"),
drawOpt_(""),
bins_(0),
rebin_(1),
stackToNEntries_(1),
rangemin_(0),
rangemax_(3),
ymin_(0),
ymax_(0),
sampleTypesToStack2D_(std::vector<Sample::SampleType>(0)),
XAxis_(""),
logX_(false),
logY_(false)
{
    // Suppress default info that canvas is printed
    gErrorIgnoreLevel = 1001;
}



void Plotter::setOptions(const TString& name, const TString& drawOpt,
                         const TString& samplesToStack2D, const TString& XAxis,
                         const int rebin, const bool stackToNEntries,
                         const bool logX, const bool logY,
                         const double& ymin, const double& ymax,
                         const double& rangemin, const double& rangemax,
                         const int bins,
                         const std::vector<double>& XAxisbins, const std::vector<double>& XAxisbinCenters)
{
    name_ = name; //Variable name you want to plot
    drawOpt_ = drawOpt; // Option used for a Draw command
    XAxis_ = XAxis; //X-axis title
    rebin_ = rebin; //Nr. of bins to be merged together
    stackToNEntries_ = stackToNEntries; // For 2D histograms fill by N entries
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
    
    // Setting the list of sample types which should be used in stack for 2D histograms
    for(int type = Sample::data; type != Sample::dummy; ++type) {
        if(samplesToStack2D == "ttJets") {
            if(type!=Sample::ttbb && type!=Sample::ttb && type!=Sample::tt2b && type!=Sample::ttcc && type!=Sample::ttother ) continue;
        }
        if(samplesToStack2D == "ttbb") {
            if(type!=Sample::ttbb ) continue;
        }
        if(samplesToStack2D == "ttLightJets") {
            if(type!=Sample::ttother ) continue;
        }
        sampleTypesToStack2D_.push_back(Sample::SampleType(type));
    }
    
}



void Plotter::producePlots()
{
    //std::cout<<"--- Beginning of plot production\n\n";
    
    // Access correction factors
    const SystematicChannelFactors globalWeights = this->scaleFactors();
    
    // Loop over all channels and systematics and produce plots
    const SystematicChannelSamples& m_systematicChannelSample(samples_.getSystematicChannelSamples());
    for(const auto& systematicChannelSamples : m_systematicChannelSample){
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        for(const auto& channelSample : systematicChannelSamples.second){
            const Channel::Channel& channel(channelSample.first);
            const std::vector<Sample>& v_sample(channelSample.second);
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
        TH1 *hist = fileReader_->GetClone<TH1>(sample.inputFile(), name_, true, false);
        if(!hist){
            // Print message only for one histogram
            if(allHistosAvailable)
                std::cout<<"Histogram ("<<name_<<") not found e.g. in file ("<<sample.inputFile()<<")\n";
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
        v_sampleHistPair_.push_back(p_sampleHist);
    }
    //std::cout<<"Number of samples used for histogram ("<<name_<<"): "<<v_sampleHistPair_.size()<<"\n";


    return allHistosAvailable;
}



void Plotter::write(const Channel::Channel& channel, const Systematic::Systematic& systematic)
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
    
    
    // Check whether Higgs sample should be drawn overlaid and/or scaled
    bool drawHiggsOverlaid(false);
    bool drawHiggsScaled(false);
    if(drawMode_ == DrawMode::overlaid){drawHiggsOverlaid = true;}
    else if(drawMode_ == DrawMode::scaledoverlaid){drawHiggsOverlaid = true; drawHiggsScaled = true;}
    
    // Specific histograms for ttbb and ttHbb to calculate signal significances
    TH1* ttbbHist(0);
    TH1* ttHbbHist(0);
    
        // Generic information about the histogram being plotted
    TString histo_name = v_sampleHistPair_.size() > 0 ? TString(v_sampleHistPair_.at(0).second->GetName()) : "";
    const bool isTH2 = TString(v_sampleHistPair_.at(0).second->ClassName()).Contains("TH2");
    
    
    // Loop over all samples and add those with identical legendEntry
    // And sort them into the categories data, Higgs, other
    LegendHistPair dataHist;
    std::vector<LegendHistPair> higgsHists;
    std::vector<LegendHistPair> stackHists;
    TH1* tmpHist(0);
    for(std::vector<SampleHistPair>::iterator i_sampleHistPair = v_sampleHistPair_.begin();
        i_sampleHistPair != v_sampleHistPair_.end(); ++i_sampleHistPair)
    {
        const Sample::SampleType& sampleType(i_sampleHistPair->first.sampleType());
        const TString& legendEntry(i_sampleHistPair->first.legendEntry());
        const TH1* hist = i_sampleHistPair->second;

        std::vector<SampleHistPair>::iterator incrementIterator(i_sampleHistPair);
        ++incrementIterator;
        const bool lastHist(incrementIterator == v_sampleHistPair_.end());
        const bool newHist(!tmpHist);

        if(newHist)tmpHist = (TH1*)hist->Clone();

        if(lastHist || (legendEntry!=incrementIterator->first.legendEntry())){
            if(!newHist)tmpHist->Add(hist);
            if(sampleType == Sample::data) dataHist = LegendHistPair(legendEntry, tmpHist);
            else if(sampleType == Sample::ttHbb) higgsHists.push_back(LegendHistPair(legendEntry, tmpHist));
            else if(sampleType == Sample::ttHother) higgsHists.push_back(LegendHistPair(legendEntry, tmpHist));
            else if(std::find(sampleTypesToStack2D_.begin(), sampleTypesToStack2D_.end(), sampleType) != sampleTypesToStack2D_.end() || !isTH2) 
            {   // Adding histogram to the stack (for 2D only if the sample was selected via samplesToStack2D option)
                stackHists.push_back(LegendHistPair(legendEntry, tmpHist));
            }
            // Assigning histograms for specific samples
            if(sampleType == Sample::ttHbb) ttHbbHist = tmpHist;
            else if(sampleType == Sample::ttbb) ttbbHist = tmpHist;
            tmpHist = 0;
            
        }
        else{
            if(!newHist)tmpHist->Add(hist);
        }
    }
    
    // Create histogram corresponding to the sum of all stacked histograms
    TH1* stacksum(0);
    for(size_t i = 0; i < stackHists.size(); ++i){
        if(i==0) stacksum = (TH1*)stackHists.at(0).second->Clone();
        else stacksum->Add((TH1*)stackHists.at(i).second);
    }
    if(!drawHiggsOverlaid){
        bool stackDefined(stacksum);
        for(auto higgsHist : higgsHists){
            if(!stackDefined){
                stacksum = (TH1*)higgsHist.second->Clone();
                stackDefined = true;
            }
            else{
                stacksum->Add((TH1*)higgsHist.second);
            }
        }
    }
    

    // Drawing signal significance for dijet_mass H->bb
    std::vector<TPaveText*> significanceLabels;
    if(!drawHiggsOverlaid && histo_name.Contains("dijet_dijet_mass")){
        significanceLabels.push_back( this->drawSignificance(ttbbHist,stacksum,85,140,0.f,"#frac{S_{ttbb}}{#sqrt{S_{ttbb}+B}}", 0) );
        significanceLabels.push_back( this->drawSignificance(ttHbbHist,stacksum,85,140,0.08,"#frac{S_{ttH}}{#sqrt{S_{ttH}+B}}", 0) );
        significanceLabels.push_back( this->drawSignificance(ttHbbHist,stacksum,85,140,0.16,"#frac{S_{ttH}}{B}", 1) );
    }
    
    // If Higgs signal scaled: scale sample and modify legend entry
    if(drawHiggsScaled){
        for(auto& higgsHist : higgsHists){
            std::stringstream ss_scaleFactor;
            const double signalScaleFactor = stacksum ? stacksum->Integral()/higgsHist.second->Integral() : 1.;
            if(isnan(signalScaleFactor) || isinf(signalScaleFactor))continue;
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
    
    // Create the stack and add entries to legend
    THStack* stack(0);
    if(dataHist.second) legend->AddEntry(dataHist.second, dataHist.first,"pe");
    if(stackHists.size()) stack = new THStack("def", "def");
    for(std::vector<LegendHistPair>::reverse_iterator higgsHist = higgsHists.rbegin(); higgsHist!=higgsHists.rend(); ++higgsHist)
    {
        higgsHist->second->SetFillStyle(0);
        higgsHist->second->SetLineWidth(2);
        legend->AddEntry(higgsHist->second, higgsHist->first,"l");
    }
    for(std::vector<LegendHistPair>::reverse_iterator stackHist = stackHists.rbegin(); stackHist!=stackHists.rend(); ++stackHist)
    {
        stackHist->second->SetLineColor(1);
        legend->AddEntry(stackHist->second, stackHist->first,"f");
    }
    for(std::vector<LegendHistPair>::const_iterator stackHist = stackHists.begin(); stackHist!=stackHists.end(); ++stackHist)
    {
        stack->Add(stackHist->second);
    }

    
    // FIXME: is this histo for error band on stack? but it is commented out ?!
    TH1* syshist(0);
    if(stacksum){
        syshist = (TH1*)stacksum->Clone();
        for(Int_t i=0; i<=syshist->GetNbinsX(); ++i){
            Double_t binc = 0;
            binc += stacksum->GetBinContent(i);
            syshist->SetBinContent(i, binc);
        }
        syshist->SetFillStyle(3004);
        syshist->SetFillColor(kBlack);
    }
    
    
    // Set x and y axis
    TH1* firstHistToDraw(0);
    if(dataHist.second && !isTH2) firstHistToDraw = (TH1*)dataHist.second->Clone();
    else if(stacksum){
        firstHistToDraw = (TH1*)stacksum->Clone();
    }
    else if(higgsHists.size()){
        firstHistToDraw = (TH1*)higgsHists.at(0).second->Clone();
    }
    else{
        std::cerr<<"ERROR in Plotter! No single sample for drawing exists\n...break\n"<<std::endl;
        exit(237);
    }
    
    if(logY_){
      // Setting minimum to >0 value
      // FIXME: Should we automatically calculate minimum value instead of the fixed value?
      firstHistToDraw->SetMinimum(1e-1);
      if(ymin_>0) firstHistToDraw->SetMinimum(ymin_);
      canvas->SetLogy();
    }
    else firstHistToDraw->SetMinimum(ymin_);

//     if(rangemin_!=0. || rangemax_!=0.) {firstHistToDraw->SetAxisRange(rangemin_, rangemax_, "X");}
    
    if(ymax_==0.){
        // Determining the highest Y value that is plotted
        float yMax = dataHist.second ? dataHist.second->GetBinContent(dataHist.second->GetMaximumBin()) : 0.f;
        if(stacksum){
            const float maxTmp = stacksum->GetBinContent(stacksum->GetMaximumBin());
            if(maxTmp > yMax) yMax = maxTmp;
        }
        for(const auto& legendHistPair : higgsHists){
            const TH1* hist = legendHistPair.second;
            const float maxTmp = hist->GetBinContent(hist->GetMaximumBin());
            if(maxTmp > yMax) yMax = maxTmp;
        }
        
        // Scaling the Y axis
        if(logY_){firstHistToDraw->SetMaximum(18.*yMax);}
        else{firstHistToDraw->SetMaximum(1.35*yMax);}
    }
    else{firstHistToDraw->SetMaximum(ymax_);}

    firstHistToDraw->GetXaxis()->SetNoExponent(kTRUE);
    
    // Blinding the data histogram region
    if((rangemin_!=0. || rangemax_!=0.) && dataHist.second) {
        if(rangemax_ < rangemin_) {
            std::cerr<<"ERROR in Plotter! Blinding region has right border with lower value than left border\n...break\n"<<std::endl;
            exit(238);
        }
        int bin1 = dataHist.second->FindBin(rangemin_);
        int bin2 = dataHist.second->FindBin(rangemax_);
        while(bin1 <= bin2) {
            dataHist.second->SetBinContent(bin1, -1e10);
            dataHist.second->SetBinError(bin1, 0);
            bin1++;
        }
    }

    
//     //Add the binwidth to the yaxis in yield plots (FIXME: works only correctly for equidistant bins)
//     TString ytitle = TString(firstHistToDraw->GetYaxis()->GetTitle()).Copy();
//     const double binWidth = firstHistToDraw->GetXaxis()->GetBinWidth(1);
//     std::ostringstream width;
//     width<<binWidth;
//     if(name_.Contains("Rapidity") || name_.Contains("Eta")){ytitle.Append(" / ").Append(width.str());}
//     else if(name_.Contains("pT") || name_.Contains("Mass") || name_.Contains("mass") || name_.Contains("MET") || name_.Contains("HT")){ytitle.Append(" / ").Append(width.str()).Append(" GeV");};
//     firstHistToDraw->GetYaxis()->SetTitle(ytitle);

    // Draw data histogram and stack and error bars
    firstHistToDraw->SetLineColor(0);
    firstHistToDraw->Draw(drawOpt_);
    if(!isTH2) {
        if(dataHist.second) dataHist.second->Draw(drawOpt_+"same e1");
        if(stack) stack->Draw(drawOpt_+"same HIST");
        gPad->RedrawAxis();
        TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0.5)");//this is frustrating and stupid but apparently necessary...
        setex1->Draw();  // error bars for data
        if(syshist) syshist->SetMarkerStyle(0);
        //syshist->Draw("same,E2");  // error bars for stack (which, stat or combined with syst ?)
        TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.)");
        setex2->Draw();  // remove error bars for data in x-direction
        if(dataHist.second) dataHist.second->Draw(drawOpt_+"same,e1");
        for(const auto& higgsHist : higgsHists){
            higgsHist.second->Draw(drawOpt_+"same");
        }
    } else {
        // Moving the pad to the left to fit the color palette
        if(stackToNEntries_) firstHistToDraw->Scale(firstHistToDraw->GetEntries()/((TH2*)firstHistToDraw)->Integral(0,-1,0,-1));
        if(ymax_ != 0) firstHistToDraw->SetMaximum(ymax_);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.1);      
    }
    
    
    // Put additional stuff to histogram
    this->drawCmsLabels(2, 8);
    this->drawDecayChannelLabel(channel);
    if(!isTH2) {
        legend->Draw("SAME");
        
        // Drawing significance labels
        for(TPaveText* label : significanceLabels) {
            label->Draw("same");
        }
    }
    if(dataHist.second && stacksum && !isTH2){
        common::drawRatio(dataHist.second, stacksum, 0, 0.5, 1.7);
        firstHistToDraw->GetXaxis()->SetLabelSize(0);
        firstHistToDraw->GetXaxis()->SetTitleSize(0);
    }

    // Create Directory for Output Plots and write them
    const TString eventFileString = common::assignFolder(ControlPlotDIR, channel, systematic);
    canvas->Print(eventFileString+name_+".eps");
    
    // Plotting Purity/Stability if this is a 2D histogram
    if(isTH2) {
        TGraphErrors* g_purity = purityStabilityGraph((TH2*)firstHistToDraw, 0);
        TGraphErrors* g_stability = purityStabilityGraph((TH2*)firstHistToDraw, 1);
        
        
        // Styling axis
        firstHistToDraw = ((TH2*)firstHistToDraw)->ProjectionX();
        firstHistToDraw->Draw("AXIS");
        firstHistToDraw->SetAxisRange(0., 1.1, "Y");
        firstHistToDraw->GetYaxis()->SetTitle("Values");
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
        canvas->Print(eventFileString+name_+"_PurStab"+".eps");
        
        delete g_purity;
        delete g_stability;
    }
    
    // Prepare additional histograms for root-file
    TH1* sumMC(0);
    TH1* sumSignal(0);
    for(const auto& sampleHistPair : v_sampleHistPair_){
        if(sampleHistPair.first.sampleType() == Sample::SampleType::ttHbb || sampleHistPair.first.sampleType() == Sample::SampleType::ttHother){
            if (sumSignal) sumSignal->Add(sampleHistPair.second);
            else sumSignal = static_cast<TH1*>(sampleHistPair.second->Clone());
            // Do not add Higgs samples to sumMC (all MC samples) in case of overlaid drawing
            if(drawHiggsOverlaid)continue;
        }
        if(sampleHistPair.first.sampleType() != Sample::SampleType::data){
            if(sumMC) sumMC->Add(sampleHistPair.second);
            else sumMC = static_cast<TH1*>(sampleHistPair.second->Clone());
        }
    }
    if(sumMC) sumMC->SetName(name_);
    
    
    //save Canvas AND sources in a root file
    TFile out_root(eventFileString+name_+"_source.root", "RECREATE");
    if(dataHist.second) dataHist.second->Write(name_+"_data");
    if(sumSignal) sumSignal->Write(name_+"_signalmc");
    if(sumMC) sumMC->Write(name_+"_allmc");
    canvas->Write(name_ + "_canvas");
    out_root.Close();
    
    firstHistToDraw->Delete();
    canvas->Clear();
    legend->Clear();
    delete canvas;
    delete legend;
}



void Plotter::setStyle(SampleHistPair& sampleHistPair)
{
    TH1* hist(sampleHistPair.second);

    hist->SetFillColor(sampleHistPair.first.color());
    hist->SetLineColor(sampleHistPair.first.color());
    hist->SetLineWidth(1);

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
}



void Plotter::drawDecayChannelLabel(const Channel::Channel& channel, const double& textSize)const
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



void Plotter::drawCmsLabels(const int cmsprelim, const double& energy, const double& textSize)const
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



TPaveText* Plotter::drawSignificance(TH1* signal, TH1* sigBkg, float min, float max, float yOffset, std::string sLabel, const int type)const
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


TGraphErrors* Plotter::purityStabilityGraph(TH2* h2d, const int type)const {

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


double Plotter::uncertaintyBinomial(const double pass, const double all)const 
{
    return (1./all)*sqrt(pass - pass*pass/all);
}



void Plotter::setGraphStyle( TGraph* graph, Style_t marker, Color_t markerColor, Size_t markerSize, 
                             Style_t line, Color_t lineColor, Size_t lineWidth)const 
{
    graph->SetMarkerStyle(marker);
    graph->SetMarkerColor(markerColor);
    graph->SetMarkerSize(markerSize);
    graph->SetLineStyle(line);
    graph->SetLineColor(lineColor);
    graph->SetLineWidth(lineWidth);
}
