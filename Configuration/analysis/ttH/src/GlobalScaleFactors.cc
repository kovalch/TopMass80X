#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>

#include <TString.h>
#include <TColorWheel.h>
#include <TH1.h>

#include "GlobalScaleFactors.h"
#include "DyScaleFactors.h"
#include "higgsUtils.h"
#include "Samples.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"











GlobalScaleFactors::GlobalScaleFactors(const std::vector<Channel::Channel>& v_channel,
                                       const std::vector<Systematic::Systematic>& v_systematic,
                                       const double& luminosityInInverseFb,
                                       const bool dyCorrection,
                                       const bool ttbbCorrection):
dyScaleFactors_(0),
ttbbScaleFactors_(0)
{
    std::cout<<"--- Beginning to set up global scale factors\n";
    
    // Need to check which samples are needed for calculating the requested scale factors
    // Drell-Yan scale factors require Samples for channels "ee" "emu" "mumu"
    std::vector<Channel::Channel> v_channelForCorrections = v_channel;
    if(dyCorrection){
        if(std::find(v_channel.begin(), v_channel.end(), Channel::ee) == v_channel.end()) v_channelForCorrections.push_back(Channel::ee);
        if(std::find(v_channel.begin(), v_channel.end(), Channel::emu) == v_channel.end()) v_channelForCorrections.push_back(Channel::emu);
        if(std::find(v_channel.begin(), v_channel.end(), Channel::mumu) == v_channel.end()) v_channelForCorrections.push_back(Channel::mumu);
    }
    
    Samples scalingSamples(v_channelForCorrections, v_systematic, 0);
    
    
    // Produce map for luminosity weights
    for(const auto& systematicChannelSamples : scalingSamples.getSystematicChannelSamples()){
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        for(const auto& channelSample : systematicChannelSamples.second){
            const Channel::Channel& channel(channelSample.first);
            const std::vector<Sample>& samples(channelSample.second);
            for(const Sample& sample : samples){
                double weight(-999.);
                if(sample.sampleType() != Sample::data) weight = sample.luminosityWeight(luminosityInInverseFb);
                //std::cout<<"\n\nLuminosity weight: "<<weight<<"\n\n";
                m_luminosityWeight_[systematic][channel].push_back(weight);
            }
        }
    }
    scalingSamples.setGlobalWeights(this);
    
    // Produce Drell-Yan scale factors
    if(dyCorrection){
        dyScaleFactors_ = new DyScaleFactors(scalingSamples);
        scalingSamples.setGlobalWeights(this);
    }
    
    // Produce ttbb scale factors
    if(ttbbCorrection){
        // FIXME: implement scaling here
        ttbbScaleFactors_ = 0;
        scalingSamples.setGlobalWeights(this);
    }
    
    std::cout<<"=== Finishing to set up global scale factors\n\n";
}



std::pair<SystematicChannelFactors, bool> GlobalScaleFactors::scaleFactors(const Samples& samples, const TString& step,
                                                                           const bool dyCorrection, const bool ttbbCorrection)const
{
    bool anyCorrectionApplied(false);
    SystematicChannelFactors systematicChannelFactors;
    
    for(const auto& systematicChannelSamples : samples.getSystematicChannelSamples()){
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        for(const auto& channelSample : systematicChannelSamples.second){
            const Channel::Channel& channel(channelSample.first);
            const std::vector<Sample>& samples(channelSample.second);
            for(size_t iSample = 0; iSample < samples.size(); ++iSample){
                const Sample& sample(samples.at(iSample));
                double weight(-999.);
                if(sample.sampleType() != Sample::data){
                    weight = m_luminosityWeight_.at(systematic).at(channel).at(iSample);
                    
                    const int dyStatus = (dyCorrection && dyScaleFactors_) ? dyScaleFactors_->applyScaleFactor(weight, step, sample, systematic) : 0;
                    if(dyStatus == 1) anyCorrectionApplied = true;
                    
                    // FIXME: implement tt+HF correction
                    const int ttbbStatus = (ttbbCorrection && ttbbScaleFactors_) ? 0 : 0;
                    if(ttbbStatus == 1) anyCorrectionApplied = true;
                }
                
                systematicChannelFactors[systematic][channel].push_back(weight);
            }
        }
    }
    
    return std::make_pair(systematicChannelFactors, anyCorrectionApplied);
}



std::pair<SystematicChannelFactors, bool> GlobalScaleFactors::scaleFactors(const Samples& samples, const TString& step)const
{
    return this->scaleFactors(samples, step, dyScaleFactors_, ttbbScaleFactors_);
}



GlobalScaleFactors::ScaleFactorStruct::ScaleFactorStruct(const TString& step, const bool dyCorrection, const bool ttbbCorrection):
step_(step),
dyCorrection_(dyCorrection),
ttbbCorrection_(ttbbCorrection),
anyCorrectionApplied_(false)
{}



bool GlobalScaleFactors::ScaleFactorStruct::scaleFactorsExist(const TString& step, const bool dyCorrection, const bool ttbbCorrection)const
{
    return step_==step && dyCorrection_==dyCorrection && ttbbCorrection_==ttbbCorrection;
}












