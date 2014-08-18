#include <fstream>
#include <sstream>

#include <TH1D.h>
#include <TString.h>

#include "EventYields.h"
#include "higgsUtils.h"
#include "Samples.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"





EventYields::EventYields(const char* outputDirectory, const Samples& samples):
fileReader_(RootFileReader::getInstance())
{
    this->produceYields(outputDirectory, samples);
}



void EventYields::produceYields(const char* outputDirectory, const Samples& samples)const
{
    std::cout<<"--- Beginning event yield table processing\n\n";
    
    // Find all histograms containing information for cutflow table (in first systematic and first channel, first file)
    const std::vector<std::pair<TString, TString> > v_nameStepPair =
        tth::nameStepPairs(samples.getSamples(samples.getSystematicChannelSamples().begin()->second.begin()->first, samples.getSystematicChannelSamples().begin()->first).at(0).inputFile(), "events_weighted_step");
    
    // Loop over steps and write yields
    for(const auto& nameStepPair : v_nameStepPair){
        this->writeYields(outputDirectory, samples, nameStepPair);
        this->writeYields(outputDirectory, samples, nameStepPair, true);
    }
    
    std::cout<<"\n=== Finishing event yield table processing\n\n";
}



void EventYields::writeYields(const char* outputDirectory,
                              const Samples& samples,
                              const std::pair<TString, TString>& nameStepPair,
                              const bool useCorrections)const
{
    // Get the scale factors from the samples
    const std::pair<SystematicChannelFactors, bool> globalWeightsPair = samples.globalWeights(nameStepPair.second, useCorrections, useCorrections);
    const SystematicChannelFactors& globalWeights = globalWeightsPair.first;
    const bool anyCorrectionApplied = globalWeightsPair.second;
    
    // Loop over systematics and channels
    for(auto systematicChannelSamples : samples.getSystematicChannelSamples()){
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        const auto& channelSamples(systematicChannelSamples.second);
        for(const auto& channelSample : systematicChannelSamples.second){
            const Channel::Channel& channel(channelSample.first);
            const std::vector<Sample>& v_sample(channelSamples.at(channel));
            
            // Do not produce corrected yields in case there are no corrections applied for this step
            if(useCorrections && !anyCorrectionApplied) continue;
            
            // Assign histogram to sample and weight it
            std::vector<SampleHistPair> v_numhist;
            for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
                const Sample& sample = v_sample.at(iSample);
                TH1D* temp_hist = fileReader_->GetClone<TH1D>(sample.inputFile(), nameStepPair.first);
                if(sample.sampleType() != Sample::data){
                    const double& weight = globalWeights.at(systematic).at(channel).at(iSample);
                    temp_hist->Scale(weight);
                }
                v_numhist.push_back(SampleHistPair(sample, temp_hist));
            }
            
            // Prepare output folder and text file
            std::ofstream eventFile;
            TString eventFileString = common::assignFolder(outputDirectory, channel, systematic);
            if(useCorrections) eventFileString.Append("corrected_");
            eventFileString.Append("events" + nameStepPair.second + ".txt");
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
                    if(i_numhist->first.sampleType() != Sample::data){
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
}















