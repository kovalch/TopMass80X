#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>

#include <TString.h>

#include "Sample.h"
#include "SampleDefinitions.h"
#include "Samples.h"
#include "GlobalScaleFactors.h"
#include "higgsUtils.h"
#include "../../common/include/sampleHelpers.h"





Samples::Samples():
globalScaleFactors_(0)
{}



Samples::Samples(const TString& filelistDirectory,
                 const std::vector<Channel::Channel>& v_channel,
                 const std::vector<Systematic::Systematic>& v_systematic,
                 const GlobalScaleFactors* globalScaleFactors):
globalScaleFactors_(globalScaleFactors)
{
    std::cout<<"--- Beginning to set up the samples\n\n";
    
    for(const auto& systematic : v_systematic){
        for(const auto& channel : v_channel){
            this->addSamples(filelistDirectory, channel, systematic);
        }
    }
    
    std::cout<<"\n=== Finishing to set up the samples\n\n";
}



void Samples::setGlobalWeights(const GlobalScaleFactors* globalScaleFactors)
{
    globalScaleFactors_ = globalScaleFactors;
}



std::vector<std::pair<TString, Sample> > Samples::setSamples(const std::vector<TString>& v_filename,
                                                             const std::map<TString, Sample>& m_samples,
                                                             const std::vector<TString>& v_sampleIdentifier)const
{
    std::vector<std::pair<TString, Sample> > result;
    
    std::vector<size_t> v_selectedFileIndex;
    for(const TString& sampleIdentifier : v_sampleIdentifier){
        // Check if specified sampleIdentifier fits with any defined sample, and access it
        if(m_samples.find(sampleIdentifier) == m_samples.end()){
            std::cerr<<"ERROR in Samples::setSamples()! Specified sample identifier not defined in samples: "
                     <<sampleIdentifier<<"\n...break\n"<<std::endl;
            exit(58);
        }
        const Sample& sample = m_samples.at(sampleIdentifier);
        
        // Loop over all filenames and check if it fits with the sample
        for(size_t iFile = 0; iFile < v_filename.size(); ++iFile){
            const TString& filename = v_filename.at(iFile);
            
            if(sample.checkFilename(filename)){
                // Check that no file is associated to two different samples
                if(std::find(v_selectedFileIndex.begin(), v_selectedFileIndex.end(), iFile) != v_selectedFileIndex.end()){
                    std::cerr<<"ERROR in Samples::setSamples()! Input file selected twice for different samples: "
                             <<filename<<"\n...break\n"<<std::endl;
                    exit(513);
                }
                v_selectedFileIndex.push_back(iFile);
                
                result.push_back(std::make_pair(filename, sample));
            }
        }
    }
    
    if(result.size() != v_filename.size()){
        std::cout<<"WARNING: Not all samples of input file list will be used (# input, # usage): "
                 <<v_filename.size()<<" , "<<result.size()<<"\n";
    }
    
    return result;
}



void Samples::addSamples(const TString& filelistDirectory,
                         const Channel::Channel& channel,
                         const Systematic::Systematic& systematic)
{
    // Read the full input filenames from the FileList
    const auto& v_filename = this->readFilelist(filelistDirectory, channel, systematic);
    
    // Add all samples as they are defined
    std::vector<std::pair<TString, Sample> > v_filenameSamplePair =
        this->setSamples(v_filename, SampleDefinitions::samples8TeV(), SampleDefinitions::selectAndOrderSamples8TeV());
    
    // Set sample options via filename
    std::vector<Sample> v_sample(this->setSampleOptions(systematic, v_filenameSamplePair));
    
    // Order files by legendEntry
    this->orderByLegend(v_sample);
    
    // Create map of maps, containing Sample per channel per systematic
    m_systematicChannelSample_[systematic][channel] = v_sample;
    
    //for(auto systematicChannelSample : m_systematicChannelSample_){
    //    for(auto channelSample : systematicChannelSample.second){
    //        for(auto sample : channelSample.second){
    //            std::cout<<"We have samples: "<<Systematic::convertSystematic(systematic)
    //                <<" , "<<Channel::convertChannel(channel)
    //                <<" , "<<sample.inputFile()
    //                <<" , "<<Channel::convertChannel(sample.finalState())<<"\n";
    //        }
    //    }
    //}
}



std::vector<TString> Samples::readFilelist(const TString& filelistDirectory,
                                           const Channel::Channel& channel,
                                           const Systematic::Systematic& systematic)const
{
    // Access fileList containing list of input root files
    const TString histoListName(filelistDirectory + "/HistoFileList_" + Systematic::convertSystematic(systematic) + "_" + Channel::convertChannel(channel) + ".txt");
    std::cout << "Reading file: " << histoListName << std::endl;
    ifstream fileList(histoListName);
    if (fileList.fail()) {
        std::cerr << "Error reading file: " << histoListName << std::endl;
        exit(1);
    }
    
    // Read in fileList to a vector
    std::vector<TString> v_filename;
    while(!fileList.eof()){
        TString filename;
        fileList>>filename;
        // Skip empty lines
        if(filename == "") continue;
        // Comment lines in FileList with '#'
        if(filename.BeginsWith("#")) continue;
        v_filename.push_back(filename);
    }
    return v_filename;
}



std::vector<Sample> Samples::setSampleOptions(const Systematic::Systematic& systematic,
                                              const std::vector<std::pair<TString, Sample> >& v_filenameSamplePair)
{
    std::vector<Sample> v_sample;

    for(auto filenameSamplePair : v_filenameSamplePair){
        TString filename(filenameSamplePair.first);
        Sample sample(filenameSamplePair.second);
        
        // Assign dilepton final state to each sample
        sample.setFinalState(this->assignFinalState(filename));
        
        // Assign specific systematic to each sample and adjust filename accordingly
        sample.setSystematic(this->assignSystematic(filename, systematic));
        
        // Check if input file really exists and set it
        sample.setInputFile(filename);
        
        v_sample.push_back(sample);
    }

    return v_sample;
}



void Samples::orderByLegend(std::vector<Sample>& v_sample)
{
    // Associate vector constituents to legend entry
    // and store all unequal legend entries
    std::vector<std::pair<TString, Sample> > v_legendSamplePair;
    std::vector<TString> v_legendEntry;
    for(auto sample : v_sample){
        const TString& legendEntry(sample.legendEntry());
        //std::cout<<"Legends before: "<<legendEntry<<std::endl;
        v_legendSamplePair.push_back(std::pair<TString, Sample>(legendEntry, sample));
        if(std::find(v_legendEntry.begin(), v_legendEntry.end(), legendEntry) != v_legendEntry.end()) continue;
        else v_legendEntry.push_back(legendEntry);
    }

    // Clear vector and fill it again in correct order
    v_sample.clear();
    for(auto legendEntry : v_legendEntry){
        int testColor(-999);
        for(auto legendSamplePair : v_legendSamplePair){
             if(legendSamplePair.first != legendEntry) continue;

             //std::cout<<"Legends after: "<<legendEntry<<std::endl;
             v_sample.push_back(legendSamplePair.second);

             // Check if all with same legend do have same colour
             int color = legendSamplePair.second.color();
             if(testColor == -999) testColor = color;
             else if(testColor != color){
                 std::cerr<<"ERROR! Samples have different colors but same legends: "<<legendEntry
                          <<"\n...break\n";
                 exit(999);
             }
         }
    }
}



Channel::Channel Samples::assignFinalState(const TString& filename)const
{
    return common::finalState(filename);
}



Systematic::Systematic Samples::assignSystematic(TString&, const Systematic::Systematic&)
//Systematic::Systematic Samples::assignSystematic(TString& filename, const Systematic::Systematic& systematic)
{
    // FIXME: adjust filename corresponding to specific systematic

    return Systematic::undefined;
}



const SystematicChannelSamples& Samples::getSystematicChannelSamples()const
{
    return m_systematicChannelSample_;
}



const std::vector<Sample>& Samples::getSamples(const Channel::Channel& channel, const Systematic::Systematic& systematic)const
{
    if(m_systematicChannelSample_.find(systematic) == m_systematicChannelSample_.end()){
        std::cerr<<"ERROR in getSamples! No samples available for requested systematic: "<<Systematic::convertSystematic(systematic)
                 <<"\n...break\n"<<std::endl;
        exit(321);
    }
    if(m_systematicChannelSample_.at(systematic).find(channel) == m_systematicChannelSample_.at(systematic).end()){
        std::cerr<<"ERROR in getSamples! No samples available for requested channel: "<<Channel::convertChannel(channel)
                 <<"\n...break\n"<<std::endl;
        exit(322);
    }
    return m_systematicChannelSample_.at(systematic).at(channel);
}



std::pair<SystematicChannelFactors, bool> Samples::globalWeights(const TString& objectname,
                                                                 const bool dyCorrection,
                                                                 const bool ttbbCorrection)const
{
    const TString fullStepname = tth::extractSelectionStepAndJetCategory(objectname);
    
    return globalScaleFactors_->scaleFactors(*this, fullStepname, dyCorrection, ttbbCorrection);
}




std::pair<SystematicChannelFactors, bool> Samples::globalWeights(const TString& objectname)const
{
    const TString fullStepname = tth::extractSelectionStepAndJetCategory(objectname);
    
    return globalScaleFactors_->scaleFactors(*this, fullStepname);
}



double Samples::luminosityInInversePb()const
{
    if(!globalScaleFactors_){
        std::cerr<<"ERROR in Samples::luminosityInInversePb()! GlobalScaleFactors do not exist, cannot return value\n...break\n"<<std::endl;
    }
    
    return globalScaleFactors_->luminosityInInversePb();
}








