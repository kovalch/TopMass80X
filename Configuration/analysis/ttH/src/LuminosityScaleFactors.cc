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
                                               const double& luminosityUncertaintyRelative,
                                               RootFileReader* const rootFileReader):
rootFileReader_(rootFileReader)
{
    std::cout<<"--- Beginning production of luminosity scale factors\n\n";
    
    this->produceScaleFactors(samples, luminosityInInversePb, luminosityUncertaintyRelative);
    
    std::cout<<"\n=== Finishing production of luminosity scale factors\n\n";
}




void LuminosityScaleFactors::produceScaleFactors(const Samples& samples, const double& luminosityInInversePb,
                                                 const double& luminosityUncertaintyRelative)
{
    // Produce map for luminosity weights
    for(const auto& systematicChannelSamples : samples.getSystematicChannelSamples()){
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        double factor(1.);
        if(systematic.type() == Systematic::lumi){
            if(systematic.variation() == Systematic::up) factor += luminosityUncertaintyRelative;
            else if (systematic.variation() == Systematic::down) factor -= luminosityUncertaintyRelative;
            else{
                std::cerr<<"ERROR in LuminosityScaleFactors::produceScaleFactors()! "
                         <<"Variation of luminosity systematic is neither up nor down\n...break\n"<<std::endl;
                exit(731);
            }
        }
        for(const auto& channelSamples : systematicChannelSamples.second){
            const Channel::Channel& channel(channelSamples.first);
            const std::vector<Sample>& samples(channelSamples.second);
            for(const Sample& sample : samples){
                double weight(-999.);
                if(sample.sampleType() != Sample::data) weight = this->luminosityWeight(sample, systematic, luminosityInInversePb*factor);
                //std::cout<<"\n\nLuminosity weight: "<<weight<<"\n\n";
                m_weight_[systematic][channel].push_back(weight);
            }
        }
    }
}



double LuminosityScaleFactors::luminosityWeight(const Sample& sample, const Systematic::Systematic& systematic, const double& luminosityInInversePb)const
{
    const double luminosityWeightPerInversePb = this->luminosityWeightPerInversePb(sample, systematic);
    
    if(luminosityWeightPerInversePb < 0.){
        std::cerr<<"ERROR in LuminosityScaleFactors::luminosityWeight()! Value is negative: "
                 <<luminosityWeightPerInversePb<<"\n...break\n"<<std::endl;
        exit(886);
    }
    return luminosityInInversePb*luminosityWeightPerInversePb;
}



double LuminosityScaleFactors::luminosityWeightPerInversePb(const Sample& sample, const Systematic::Systematic& systematic)const
{
    const double& crossSection = sample.crossSection(systematic);
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








