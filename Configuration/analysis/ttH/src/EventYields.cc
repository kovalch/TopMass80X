#include <fstream>
#include <sstream>

#include <TH1D.h>
#include <TString.h>

#include "EventYields.h"
#include "higgsUtils.h"
#include "Samples.h"
#include "DyScaleFactors.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"





/// Folder for event yields output
constexpr const char* YieldDIR = "EventYields";





EventYields::EventYields(const Samples& samples, const double& luminosity, const DyScaleFactors& dyScaleFactors):
luminosity_(luminosity),
dyScaleFactors_(dyScaleFactors),
fileReader_(RootFileReader::getInstance())
{
    this->produceYields(samples);
}



void EventYields::produceYields(const Samples& samples)const
{
    std::cout<<"--- Beginning event yield table processing\n\n";
    
    // Find all histograms containing information for cutflow table (in systematic Nominal and first channel, first file)
    const std::vector<std::pair<TString, TString> > v_nameStepPair =
        tth::nameStepPairs(samples.getSamples(samples.getSystematicChannelSamples().at(Systematic::nominal).begin()->first, Systematic::nominal).at(0).inputFile(), "events_weighted_step");
    
    // Loop over systematics (exclude all but Nominal - so outer loop could be removed) and channels
    for(const auto& systematicChannelSamples : samples.getSystematicChannelSamples()){
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        if(systematic != Systematic::nominal)continue;
        for(const auto& channelSample : systematicChannelSamples.second){
            const Channel::Channel& channel(channelSample.first);
            this->writeYields(channel, channelSample.second, v_nameStepPair);
            this->writeYields(channel, channelSample.second, v_nameStepPair, true);
        }
    }
    std::cout<<"\n=== Finishing event yield table processing\n\n";
}



void EventYields::writeYields(const Channel::Channel& channel, const std::vector<Sample>& v_sample,
                              const std::vector<std::pair<TString, TString> >& v_nameStepPair,
                              const bool useCorrections)const
{
    // Loop over all selection steps writing out event yields
    for(std::vector<std::pair<TString, TString> >::const_iterator i_nameStepPair = v_nameStepPair.begin(); i_nameStepPair != v_nameStepPair.end(); ++i_nameStepPair){
        
        // Assign histogram to sample and weight it to luminosity
        std::vector<SampleHistPair> v_numhist;
        for(const auto& sample : v_sample){
            TH1D* temp_hist = fileReader_->GetClone<TH1D>(sample.inputFile(), i_nameStepPair->first);
            if(sample.sampleType() != Sample::data){
                const double& lumiWeight = sample.luminosityWeight(luminosity_);
                temp_hist->Scale(lumiWeight);
            }
            v_numhist.push_back(SampleHistPair(sample, temp_hist));
        }
        
        // Apply possible corrections to yields
        bool dyScalingApplied(false);
        if(useCorrections){
            // Apply Drell-Yan scaling
            for(auto sampleHistPair : v_numhist){
                const int dyScaleStatus = dyScaleFactors_.applyDyScaleFactor(sampleHistPair.second, sampleHistPair.first,
                                                                             Systematic::nominal, true);
                if(dyScaleStatus==-1){
                    dyScalingApplied = false;
                    break;
                }
                else if(dyScaleStatus==1){
                    dyScalingApplied = true;
                }
            }
        }
        
        // If useCorrections is set, but no corrections are applied, no need to store output
        if(useCorrections && !dyScalingApplied) continue;
        
        // Prepare output folder and text file
        // At present, event yields are only possible for nominal systematic
        std::ofstream eventFile;
        TString eventFileString = common::assignFolder(YieldDIR, channel, Systematic::nominal);
        if(useCorrections) eventFileString.Append("corrected_");
        eventFileString.Append("events" + i_nameStepPair->second + ".txt");
        eventFile.open(eventFileString.Data());
        
        // Make output for tables
        double tmp_num = 0;
        double bg_num = 0;
        for(auto i_numhist = v_numhist.begin(); i_numhist != v_numhist.end(); ++i_numhist){
            auto iterator = i_numhist;
            ++iterator;
            tmp_num += i_numhist->second->Integral();
            if(i_numhist == --(v_numhist.end())){
                eventFile<<i_numhist->first.legendEntry()<<": "<<tmp_num<<std::endl;
                bg_num += tmp_num;
                tmp_num = 0;
            }
            else if(i_numhist->first.legendEntry() != iterator->first.legendEntry()){
                eventFile<<i_numhist->first.legendEntry()<<": "<<tmp_num<<std::endl;
                if(i_numhist->first.sampleType() != Sample::SampleType::data){
                    bg_num+=tmp_num;
                }
                tmp_num = 0;
            }
        }
        eventFile<<"Total background: "<<bg_num<<std::endl;
        
        // Close text file
        eventFile.close();
        //std::cout<<"\nEvent yields saved in "<<eventFilestring<<"\n"<<std::endl;
     }
}















