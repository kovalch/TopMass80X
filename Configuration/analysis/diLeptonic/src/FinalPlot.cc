
#include <fstream>
#include <iostream>
// #include <cstdio>
#include <sstream>
// #include <cmath>
// #include <iomanip>

//#include <TCanvas.h>
//#include <TLegend.h>
//#include <TExec.h>
//#include <TStyle.h>
//#include <TSystem.h>
//#include <TMath.h>
//#include <TROOT.h>
//#include <THStack.h>
//#include <TDirectory.h>
//#include <TFile.h>
//#include <TString.h>
//#include <TTree.h>
//#include <TH1.h>
//#include <TH2.h>
//#include <TGraph.h>

// #include "utils.h"
#include "FinalPlot.h"
//#include "ttbarUtils.h"
#include "Samples.h"
//#include "Sample.h"
#include "../../common/include/utils.h"
//#include "../../common/include/plotterUtils.h"
//#include "UsefulTools.h"
//#include "Output.h"

//#include "TUnfold.h"
//#include <TLine.h>
//#include <TMath.h>
//#include "TUnfoldBinning.h"
//#include "TUnfoldDensity.h"





FinalPlot::FinalPlot(const Samples& samples, const double lumi, const double topxsec):

samples_(samples),
lumi_(lumi),
topxsec_(topxsec),
nD_(-999),
v_plotName_(std::vector<TString>(0))
{
    
}



void FinalPlot::setOptions(const std::vector<TString> v_plotName)
{
    
    //std::cout<<"[FinalPlot]:--- Beginning of plot options settings\n\n";
    
    //Clearing previous declarations  
    this->clearMemory();
    
    v_plotName_ = v_plotName;
    nD_ = (int)v_plotName_.size();
    TString nDflag;
    if(nD_==1) nDflag="1d";
    else if(nD_==2) nDflag="2d";
    else{
         std::cout<<"ERROR in Plotter! nD_ is not equal to 1 or 2.\n...break\n"<<std::endl;
         exit(237);
     }
    
    //std::cout<<"[FinalPlot]:--- Finishing of plot options settings\n\n";
}


void FinalPlot::clearMemory()
{
    v_plotName_.clear();
}


void FinalPlot::producePlots()
{
    
        //std::cout<<"[FinalPlot]:--- Beginning of plot production\n\n";
    
    // Loop over all channels and systematics and produce plots
    const SystematicChannelSamples& m_systematicChannelSample(samples_.getSystematicChannelSamples());
    for(const auto& systematicChannelSamples : m_systematicChannelSample){                     //systematic loop
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        for(const auto& channelSample : systematicChannelSamples.second){                     //channel loop
            ///this->clearMemoryPerSystematicChannel();
            const Channel::Channel& channel(channelSample.first);
            
            //Set plotsFolder_
            plotsFolder_ = common::assignFolder("Plots",channel,systematic,v_plotName_.at(0)+"-"+v_plotName_.at(1));
            
            //std::cout << "channel: " << channel << "systematic: " << systematic << std::endl;
            
        }//channel loop
        
    }//systematics loop
    
    //std::cout<<"\n[FinalPlot]:=== Finishing of plot production\n\n";
    
}





