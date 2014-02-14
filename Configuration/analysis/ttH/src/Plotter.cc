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
#include <TGaxis.h>
#include <TPaveText.h>
#include <TError.h>

#include "Plotter.h"
#include "higgsUtils.h"
#include "Samples.h"
#include "DyScaleFactors.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"
#include "../../common/include/plotterUtils.h"





/// Folder for Control Plots output
constexpr const char* ControlPlotDIR = "Plots";





Plotter::Plotter(const Samples& samples,
                 const double& luminosity,
                 const DyScaleFactors& dyScaleFactors,
                 const DrawMode::DrawMode& drawMode):
samples_(samples),
luminosity_(luminosity),
dyScaleFactors_(dyScaleFactors),
drawMode_(drawMode),
fileReader_(RootFileReader::getInstance()),
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
doDYScale_(false),
scaleMCtoData_(false),
ttbbScale_(1.0)
{
    // Suppress default info that canvas is printed
    gErrorIgnoreLevel = 1001;
}



void Plotter::setOptions(const TString& name, const TString&,
                         const TString& YAxis, const TString& XAxis,
                         const int rebin, const bool doDYScale,
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
    doDYScale_ = doDYScale; //Apply DY scale factor?
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
}



void Plotter::producePlots()
{
    //std::cout<<"--- Beginning of plot production\n\n";

    const SystematicChannelSamples& m_systematicChannelSample(samples_.getSystematicChannelSamples());
    for(const auto& systematicChannelSamples : m_systematicChannelSample){
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        for(const auto& channelSample : systematicChannelSamples.second){
            const Channel::Channel& channel(channelSample.first);
            const std::vector<Sample>& v_sample(channelSample.second);
            if(!this->prepareDataset(v_sample, systematic)){
                std::cout<<"WARNING! Cannot find histograms for all datasets, for (channel/systematic): "
                         << Channel::convertChannel(channel) << "/" << Systematic::convertSystematic(systematic)
                         <<"\n... skip this plot\n";
                return;
            }
            this->write(channel, systematic);
        }
    }

    //std::cout<<"\n=== Finishing of plot production\n\n";
}



bool Plotter::prepareDataset(const std::vector<Sample>& v_sample, const Systematic::Systematic& systematic)
{
    bool allHistosAvailable(true);


    // Associate histogram to dataset if histogram can be found
    v_sampleHistPair_.clear();
    TH1::AddDirectory(kFALSE);
    for(const auto& sample : v_sample){
        SampleHistPair p_sampleHist;
        TH1D *hist = fileReader_->GetClone<TH1D>(sample.inputFile(), name_, true, false);
        if (!hist){
            // Print message only for one histogram
            if(allHistosAvailable)
                std::cout<<"Histogram ("<<name_<<") not found e.g. in file ("<<sample.inputFile()<<")\n";
            p_sampleHist = SampleHistPair(sample, 0);
            allHistosAvailable = false;
        }
        else{
            //Rescaling to the data luminosity
            if(sample.sampleType() != Sample::data){
                const double lumiWeight = Tools::luminosityWeight(sample, luminosity_, fileReader_);
                hist->Scale(lumiWeight);
            }

            // Drell-Yan reweighting
            if(doDYScale_){
                dyScaleFactors_.applyDyScaleFactor(hist, sample, systematic);
            }

            // Set style
            common::setHHStyle(*gStyle);

            // Clone histogram directly here
            TH1D* histClone = (TH1D*) hist->Clone();
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
        setStyle(sampleHistPair, true);
    }
    
    
    // Check whether Higgs sample should be drawn overlaid and/or scaled
    bool drawHiggsOverlaid(false);
    bool drawHiggsScaled(false);
    if(drawMode_ == DrawMode::overlaid){drawHiggsOverlaid = true;}
    else if(drawMode_ == DrawMode::scaledoverlaid){drawHiggsOverlaid = true; drawHiggsScaled = true;}
    
    
    // Loop over all samples and add those with identical legendEntry
    // And sort them into the categories data, Higgs, other
    LegendHistPair dataHist;
    std::vector<LegendHistPair> higgsHists;
    std::vector<LegendHistPair> stackHists;
    TH1D* tmpHist(0);
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

        if(newHist)tmpHist = (TH1D*)hist->Clone();

        if(lastHist || (legendEntry!=incrementIterator->first.legendEntry())){
            if(!newHist)tmpHist->Add(hist);
            if(sampleType == Sample::data) dataHist = LegendHistPair(legendEntry, tmpHist);
            else if(sampleType == Sample::higgssignal) higgsHists.push_back(LegendHistPair(legendEntry, tmpHist));
            else stackHists.push_back(LegendHistPair(legendEntry, tmpHist));
            tmpHist = 0;
        }
        else{
            if(!newHist)tmpHist->Add(hist);
        }
    }


    // Scaling the tt+bb component of the tt+jets sample
//     double ttbbInt = 0.;
//     double ttOtherInt = 0.;
//     double ttScale = 1.;
//     TH1D* ttbbHist(0);
//     TH1D* ttOtherHist(0);
//     for(size_t i = 0; i < stackHists.size(); ++i){
//         if(stackHists.at(i).first == "t#bar{t}b#bar{b}") ttbbHist = (TH1D*)stackHists.at(i).second;
//         if(stackHists.at(i).first == "t#bar{t}Other") ttOtherHist = (TH1D*)stackHists.at(i).second;
//     }
//     if(ttbbHist) ttbbInt = ttbbHist->Integral();
//     if(ttOtherHist) ttOtherInt = ttOtherHist->Integral();
//     double ttInt = ttbbInt + ttOtherInt;
//     double ttIntScaled = ttbbInt*ttbbScale_ + ttOtherInt;
//     if(ttInt>0. && ttIntScaled>0.) ttScale = ttInt/ttIntScaled;
//     ttbbHist->Scale(ttbbScale_*ttScale);
//     ttOtherHist->Scale(ttScale);
    
    
    // Create histogram corresponding to the sum of all stacked histograms
    TH1D* stacksum(0);
    for(size_t i = 0; i < stackHists.size(); ++i){
        if(i==0) stacksum = (TH1D*)stackHists.at(0).second->Clone();
        else stacksum->Add((TH1D*)stackHists.at(i).second);
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


    // Scale the summed MC sample to the data
//     if(scaleMCtoData_) {
//         // Estimate the scaling factor for MC histos to have the same integral as Data
//         float MCtoData = (float)dataHist.second->Integral()/(float)stacksum->Integral();
//         // Scaling MC to match data
//         for (const auto& stackHist : stackHists) { stackHist.second->Scale(MCtoData); }
//         for (const auto& higgsHist : higgsHists) { higgsHist.second->Scale(MCtoData); }
//         stacksum->Scale(MCtoData);
//     }


    // Drawing signal significance for dijet_mass H->bb
    TPaveText* significanceLabelTTH(0);
    TPaveText* significanceLabelTTbb(0);
//     if(higgsHists.size()){
//         const TString histo_name = TString(higgsHists.at(0).second->GetName());
//         if(histo_name.Contains("dijet_dijet_mass")){
//             TH1D* stacksum_ttH = (TH1D*)stacksum->Clone();
//             stacksum_ttH->Add((TH1D*)higgsHists.at(0).second);
//             TH1D* ttHbbHist = 0;
//             for (const auto& higgsHist : higgsHists) { if(higgsHist.first == "t#bar{t}H (b#bar{b})") ttHbbHist = (TH1D*)higgsHist.second; }
//             if(ttbbHist) significanceLabelTTbb = this->drawSignificance(ttbbHist,stacksum,85,140,0.f,"ttbb");
//             if(ttHbbHist) significanceLabelTTH = this->drawSignificance(ttHbbHist,stacksum_ttH,85,140,0.1,"ttH");
//         }
//     }
    
    
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
    for(const auto& stackHist : stackHists){
        stackHist.second->SetLineColor(1);
        legend->AddEntry(stackHist.second, stackHist.first,"f");
        stack->Add(stackHist.second);
    }
    for(const auto& higgsHist : higgsHists){
        higgsHist.second->SetFillStyle(0);
        higgsHist.second->SetLineWidth(2);
        legend->AddEntry(higgsHist.second, higgsHist.first,"l");
    }

    // In fact the legend is filled oppositely than the stack, so this is used for turning the order (but not completely)
//    legend = ControlLegend(dataHist, stackHists, higgsHists, legend);

    
    // FIXME: is this histo for error band on stack? but it is commented out ?!
    TH1D* syshist(0);
    if(stacksum){
        syshist = (TH1D*)stacksum->Clone();
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
    if(dataHist.second) firstHistToDraw = (TH1*)dataHist.second->Clone();
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

    if(rangemin_!=0. || rangemax_!=0.) {firstHistToDraw->SetAxisRange(rangemin_, rangemax_, "X");}

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

    
    //Add the binwidth to the yaxis in yield plots (FIXME: works only correctly for equidistant bins)
    TString ytitle = TString(firstHistToDraw->GetYaxis()->GetTitle()).Copy();
    const double binWidth = firstHistToDraw->GetXaxis()->GetBinWidth(1);
    std::ostringstream width;
    width<<binWidth;
    if(name_.Contains("Rapidity") || name_.Contains("Eta")){ytitle.Append(" / ").Append(width.str());}
    else if(name_.Contains("pT") || name_.Contains("Mass") || name_.Contains("mass") || name_.Contains("MET") || name_.Contains("HT")){ytitle.Append(" / ").Append(width.str()).Append(" GeV");};
    firstHistToDraw->GetYaxis()->SetTitle(ytitle);

    // Draw data histogram and stack and error bars
    firstHistToDraw->SetLineColor(0);
    firstHistToDraw->Draw();
    if(dataHist.second) dataHist.second->Draw("same e1");
    if(stack) stack->Draw("same HIST");
    gPad->RedrawAxis();
    TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0.5)");//this is frustrating and stupid but apparently necessary...
    setex1->Draw();  // error bars for data
    if(syshist) syshist->SetMarkerStyle(0);
    //syshist->Draw("same,E2");  // error bars for stack (which, stat or combined with syst ?)
    TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.)");
    setex2->Draw();  // remove error bars for data in x-direction
    if(dataHist.second) dataHist.second->Draw("same,e1");
    for(const auto& higgsHist : higgsHists){
        higgsHist.second->Draw("same");
    }
    
    
    // Put additional stuff to histogram
    this->drawCmsLabels(2, 8);
    this->drawDecayChannelLabel(channel);
    if(significanceLabelTTH) significanceLabelTTH->Draw("same");
    if(significanceLabelTTbb) significanceLabelTTbb->Draw("same");
    legend->Draw("SAME");
    if(dataHist.second && stacksum){
        common::drawRatio(dataHist.second, stacksum, 0, 0.5, 1.7);
        firstHistToDraw->GetXaxis()->SetLabelSize(0);
        firstHistToDraw->GetXaxis()->SetTitleSize(0);
    }

    // Create Directory for Output Plots and write them
    const TString eventFileString = common::assignFolder(ControlPlotDIR, channel, systematic);
    canvas->Print(eventFileString+name_+".eps");
    
    
    // Prepare additional histograms for root-file
    TH1* sumMC(0);
    TH1* sumSignal(0);
    for(const auto& sampleHistPair : v_sampleHistPair_){
        if(sampleHistPair.first.sampleType() == Sample::SampleType::higgssignal){
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



void Plotter::setStyle(SampleHistPair& sampleHistPair, const bool isControlPlot)
{
    TH1* hist(sampleHistPair.second);

    hist->SetFillColor(sampleHistPair.first.color());
    hist->SetLineColor(sampleHistPair.first.color());
    hist->SetLineWidth(1);

    if(XAxis_ == "-")XAxis_ = hist->GetXaxis()->GetTitle();
    if(YAxis_ == "-")YAxis_ = hist->GetYaxis()->GetTitle();

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



TLegend* Plotter::controlLegend(const LegendHistPair& dataHist, const std::vector<LegendHistPair>& stackHists,
                                const std::vector<LegendHistPair>& higgsHists, TLegend* leg)
{
    //hardcoded ControlPlot legend
    std::vector<TString> v_orderedLegend;
    v_orderedLegend.push_back("Data");
    v_orderedLegend.push_back("t#bar{t} Signal");
    v_orderedLegend.push_back("t#bar{t} Other");
    v_orderedLegend.push_back("Single Top");
    v_orderedLegend.push_back("W+Jets");
    v_orderedLegend.push_back("Z / #gamma* #rightarrow ee/#mu#mu");
    v_orderedLegend.push_back("Z / #gamma* #rightarrow #tau#tau");
    v_orderedLegend.push_back("Diboson");
    v_orderedLegend.push_back("QCD Multijet");
    v_orderedLegend.push_back("t#bar{t}H (incl.)");
    v_orderedLegend.push_back("t#bar{t}H (b#bar{b})");

    leg->Clear();
    leg->SetX1NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength()-0.25);
    leg->SetY1NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength()-0.30);
    leg->SetX2NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength());
    leg->SetY2NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength());
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextAlign(12);
    for(const auto& orderedLegend : v_orderedLegend){
        if(dataHist.first.Contains(orderedLegend)){
            leg->AddEntry(dataHist.second, dataHist.first, "pe");
            continue;
        }
        for(const auto& stackHist : stackHists){
            if(stackHist.first.Contains(orderedLegend)){
                leg->AddEntry(stackHist.second, stackHist.first, "f");
                break;
            }
        }
        for(const auto& higgsHist : higgsHists){
            if(higgsHist.first.Contains(orderedLegend)){
                leg->AddEntry(higgsHist.second, higgsHist.first, "l");
                break;
            }
        }
    }
    return leg;
}



void Plotter::drawDecayChannelLabel(const Channel::Channel& channel, const double& textSize)
{
    TPaveText* decayChannel = new TPaveText();

    decayChannel->AddText(Channel::channelLabel(channel).c_str());

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



void Plotter::drawCmsLabels(const int cmsprelim, const double& energy, const double& textSize)
{
    const char* text;
    if(cmsprelim ==2 ) {//Private work for PhDs students
        text = "Private Work, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV";
    } else if (cmsprelim==1) {//CMS preliminary label
        text = "CMS Preliminary, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV";
    } else {//CMS label
        text = "CMS, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV";
    }

    TPaveText *label = new TPaveText();
    label->SetX1NDC(gStyle->GetPadLeftMargin());
    label->SetY1NDC(1.0-gStyle->GetPadTopMargin());
    label->SetX2NDC(1.0-gStyle->GetPadRightMargin());
    label->SetY2NDC(1.0);
    label->SetTextFont(42);
    label->AddText(Form(text, luminosity_/1000, energy));
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    if (textSize!=0) label->SetTextSize(textSize);
    label->SetTextAlign(32);
    label->Draw("same");
}



TPaveText* Plotter::drawSignificance(TH1* signal, TH1* sigBkg, float min, float max, float yOffset, std::string sLabel)
{
    if(max<=min) {
        std::cout<<"Wrong range for signal significance in  histogram ("<<name_<<")\n";
        return 0;
    }

    // Finding the bin range corresponding to [min;max]
    int bin1 = signal->FindBin(min);
    int bin2 = signal->FindBin(max);

    // Calculating integral of the signal and background
    float sigInt = signal->Integral(bin1,bin2);
    float sigBkgInt = sigBkg->Integral(bin1,bin2);

    float sigSign = sigInt/sqrt(sigBkgInt);

    char text[40];
    sprintf(text,"#frac{S_{%s}}{#sqrt{S_{%s}+B}} = %.2f", sLabel.c_str(), sLabel.c_str(), sigSign);

    TPaveText *label = new TPaveText();
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







