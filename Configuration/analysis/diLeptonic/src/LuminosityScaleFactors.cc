#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>

#include <TH1.h>

#include "LuminosityScaleFactors.h"
#include "Samples.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"






LuminosityScaleFactors::LuminosityScaleFactors(const Samples& samples,
                                               const double& luminosityInInversePb,
                                               RootFileReader* const rootFileReader):
rootFileReader_(rootFileReader)
{
    std::cout<<"--- Beginning production of luminosity scale factors\n\n";
    
    this->produceScaleFactors(samples, luminosityInInversePb);
    
    std::cout<<"\n=== Finishing production of luminosity scale factors\n\n";
}




void LuminosityScaleFactors::produceScaleFactors(const Samples& samples, const double& luminosityInInversePb)
{
    // Produce map for luminosity weights
    for(const auto& systematicChannelSamples : samples.getSystematicChannelSamples()){
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        for(const auto& channelSamples : systematicChannelSamples.second){
            const Channel::Channel& channel(channelSamples.first);
            const std::vector<Sample>& samples(channelSamples.second);
            for(const Sample& sample : samples){
                double weight(-999.);
                if(sample.sampleType() != Sample::data) weight = this->luminosityWeight(sample, luminosityInInversePb);
                //std::cout<<"\n\nLuminosity weight: "<<weight<<"\n\n";
                
                //systematic, bg
                    if(systematic.type() == Systematic::bg && sample.sampleType() != Sample::data && sample.sampleType() != Sample::ttother && 
                        sample.sampleType() != Sample::dyee && sample.sampleType() != Sample::dymumu && sample.sampleType() != Sample::dytautau )
                    {
                            if(systematic.variation() == Systematic::up)weight = weight*1.3;
                            if(systematic.variation() == Systematic::down)weight = weight*0.7;
                    }
                // systematic, dy
                    if(systematic.type() == Systematic::dy && (sample.sampleType() == Sample::dyee || sample.sampleType() == Sample::dymumu || sample.sampleType() == Sample::dytautau) )
                    {
                            if(systematic.variation() == Systematic::up)weight = weight*1.3;
                            if(systematic.variation() == Systematic::down)weight = weight*0.7;
                    }
                // ...
                    
                m_weight_[systematic][channel].push_back(weight);
            }
        }
    }
}



double LuminosityScaleFactors::luminosityWeight(const Sample& sample, const double& luminosityInInversePb)const
{
    const double luminosityWeightPerInversePb = this->luminosityWeightPerInversePb(sample);
    
    if(luminosityWeightPerInversePb < 0.){
        std::cerr<<"ERROR in LuminosityScaleFactors::luminosityWeight()! Value is negative: "
                 <<luminosityWeightPerInversePb<<"\n...break\n"<<std::endl;
        exit(886);
    }
    return luminosityInInversePb*luminosityWeightPerInversePb;
}



double LuminosityScaleFactors::luminosityWeightPerInversePb(const Sample& sample)const
{
    const double& crossSection = sample.crossSection();
    if(crossSection <= 0.){
        std::cerr<<"ERROR in LuminosityScaleFactors::calculateLuminosityWeight()! Sample XSection is <0 (sample, xsec): "
                 <<sample.inputFile()<<" , "<<crossSection<<"\n...break\n"<<std::endl;
        exit(887);
    }
    
    const TH1* const h_weightedEvents = rootFileReader_->Get<TH1>(sample.inputFile(), "weightedEvents");
    const double weightedEvents(h_weightedEvents->GetBinContent(1));
    
    //std::cout<<"Input file: "<<inputFileName_<<std::endl;
    //std::cout<<"Xsection, weighted events, lumi weight per inverse fb: "
    //         <<crossSection<<" , "<<weightedEvents<<" , "<<crossSection/weightedEvents*1000.<<std::endl;
    return crossSection/weightedEvents;
}



const SystematicChannelFactors& LuminosityScaleFactors::scaleFactorMap()const
{
    return m_weight_;
}








