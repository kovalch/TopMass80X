#include <vector>
#include <map>
#include <set>
#include <iostream>
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
                                                             const std::vector<TString>& v_sampleIdentifier,
                                                             const bool hasPseudodata)const
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
            
            // Adding the filename if it suits the current sample
            if(sample.checkFilename(filename)){
                // Check that no file is associated to two different samples
                if(std::find(v_selectedFileIndex.begin(), v_selectedFileIndex.end(), iFile) != v_selectedFileIndex.end()){
                    std::cerr<<"ERROR in Samples::setSamples()! Input file selected twice for different samples: "
                             <<filename<<"\n...break\n"<<std::endl;
                    exit(513);
                }
                // Adding the file index to check for double-counting
                // Ignoring files for pseudodata: to allow the same file with different normalisation in pseudodata
                if(sample.sampleType() != Sample::pseudodata) v_selectedFileIndex.push_back(iFile);
                
                result.push_back(std::make_pair(filename, sample));
            } else continue;
            
            // Adding the file to pseudodata sample if it exists
            if(hasPseudodata) {
                bool hasSampleInPseudodata = false;
                Sample samplePseudodata_reference;
                // Checking whether this file is already included in any of modified pseudodata files
                for(auto nameSamplePair : m_samples) {
                    Sample& theSample = nameSamplePair.second;
                    if(theSample.sampleType() != Sample::pseudodata) continue;
                    if(std::find(v_sampleIdentifier.begin(), v_sampleIdentifier.end(), nameSamplePair.first) == v_sampleIdentifier.end()) continue;
                    samplePseudodata_reference = theSample;
                    if(theSample.containsFilenamesOfSample(sample, true)) {
                        hasSampleInPseudodata = true;
                        break;
                    }
                }
                if(!hasSampleInPseudodata) {
                    // Making a copy of the sample with visual properties and type of 0.the pseudodata sample
                    Sample samplePseudodata(sample);
                    samplePseudodata.setSampleType(samplePseudodata_reference.sampleType());
                    samplePseudodata.setLegendEntry(samplePseudodata_reference.legendEntry());
                    samplePseudodata.setColor(samplePseudodata_reference.color());
                    result.push_back(std::make_pair(filename, samplePseudodata));
                }
            }
        }
    }
    
    if(result.size() < v_filename.size()){
        std::cout<<"WARNING: Not all samples of input file list will be used (# input, # usage): "
                 <<v_filename.size()<<" , "<<result.size()<<"\n";
    }
    
    return result;
}



void Samples::addSamples(const TString& filelistDirectory,
                         const Channel::Channel& channel,
                         const Systematic::Systematic& systematic)
{
    // Full input filenames from the systematic and the nominal FileList
    std::vector<std::pair<TString, Sample> > v_filenameSamplePair;
    std::vector<std::pair<TString, Sample> > v_filenameSamplePairNominal;
    const bool hasPseudodata = SampleDefinitions::usingPseudodata(SampleDefinitions::samples8TeV(), SampleDefinitions::selectAndOrderSamples8TeV());
    
    // Add all samples as they are defined for given systematic (except of systematics which only scale nominal samples)
    if(systematic.type() != Systematic::lumi && 
        std::find(Systematic::crossSectionTypes.begin(), Systematic::crossSectionTypes.end(), systematic.type()) == Systematic::crossSectionTypes.end() &&
        std::find(Systematic::tthfFractionTypes.begin(), Systematic::tthfFractionTypes.end(), systematic.type()) == Systematic::tthfFractionTypes.end()
    ) {
        const auto& v_filename = common::readFilelist(filelistDirectory, channel, systematic);
        v_filenameSamplePair =
            this->setSamples(v_filename, SampleDefinitions::samples8TeV(), SampleDefinitions::selectAndOrderSamples8TeV(), hasPseudodata);
    }
    
    // Add nominal samples for those not varied (i.e. not found in systematic FileList)
    if(systematic.type() != Systematic::nominal){
        const auto& v_filenameNominal = common::readFilelist(filelistDirectory, channel, Systematic::nominalSystematic());
        v_filenameSamplePairNominal =
            this->setSamples(v_filenameNominal, SampleDefinitions::samples8TeV(), SampleDefinitions::selectAndOrderSamples8TeV(), hasPseudodata);
    }
    
    // Set sample options via filename
    std::vector<Sample> v_sample(this->setSampleOptions(systematic, v_filenameSamplePair, v_filenameSamplePairNominal));
    
    // If nominal samples will be merged: reordering samples based on legends
    if(v_filenameSamplePairNominal.size()>0 || hasPseudodata) {
        // Getting a list of legends in a proper order
        std::vector<TString> v_legend = SampleDefinitions::legendList(SampleDefinitions::samples8TeV(), SampleDefinitions::selectAndOrderSamples8TeV());
        // Order files by legendEntry
        this->orderByLegend(v_sample, v_legend);
    }
    
    // Create map of maps, containing Sample per channel per systematic
    m_systematicChannelSample_[systematic][channel] = v_sample;
}



std::vector<Sample> Samples::setSampleOptions(const Systematic::Systematic& systematic,
                                              const std::vector<std::pair<TString, Sample> >& v_filenameSamplePair,
                                              const std::vector<std::pair<TString, Sample> >& v_filenameSamplePairNominal)
{
    std::vector<Sample> v_sample;
    
    // Set all samples processed for given systematic
    for(auto filenameSamplePair : v_filenameSamplePair){
        const TString& filename(filenameSamplePair.first);
        Sample sample(filenameSamplePair.second);
        
        // Assign real final state to each sample, ie. only "ee", "emu", "mumu", but not "combined"
        sample.setFinalState(common::finalState(filename));
        
        // Assign specific systematic to each sample
        sample.setSystematic(systematic);
        
        // Check if input file really exists and set it
        sample.setInputFile(filename);
        
        v_sample.push_back(sample);
    }
    
    // Set all samples not processed for given systematic to nominal one
    for(auto filenameSamplePair : v_filenameSamplePairNominal){
        const TString& filename(filenameSamplePair.first);
        Sample sample(filenameSamplePair.second);
        
        // Check if sample is already contained in systematic-specific samples
        bool inSystematicSamples(false);
        for(const auto& filenameSamplePairSystematic : v_filenameSamplePair){
            if(sample.sampleType() != filenameSamplePairSystematic.second.sampleType()) continue;
            inSystematicSamples = true;
            break;
        }
        if(inSystematicSamples) continue;
        
        // Set real final state, nominal systematic, and corresponding filename
        sample.setFinalState(common::finalState(filename));
        sample.setSystematic(Systematic::nominalSystematic());
        sample.setInputFile(filename);
        
        v_sample.push_back(sample);
    }
    
    return v_sample;
}



void Samples::orderByLegend(std::vector<Sample>& v_sample, const std::vector<TString>& v_legend)const
{
    // Associate vector constituents to legend entry
    // and store all unequal legend entries
    std::vector<Sample> v_sample_unordered;
    for(auto sample : v_sample) v_sample_unordered.push_back(sample);

    // Clear vector and fill it again in correct order
    v_sample.clear();
    for(auto legendEntry : v_legend){
        int testColor(-999);
        for(auto sample : v_sample_unordered){
             if(sample.legendEntry() != legendEntry) continue;

             v_sample.push_back(sample);

             // Check if all with same legend do have same colour
             int color = sample.color();
             if(testColor == -999) testColor = color;
             else if(testColor != color){
                 std::cerr<<"ERROR! Samples have different colors but same legends: "<<legendEntry
                          <<"\n...break\n";
                 exit(999);
             }
         }
    }
}



const SystematicChannelSamples& Samples::getSystematicChannelSamples()const
{
    return m_systematicChannelSample_;
}



const std::vector<Sample>& Samples::getSamples(const Channel::Channel& channel, const Systematic::Systematic& systematic)const
{
    if(m_systematicChannelSample_.find(systematic) == m_systematicChannelSample_.end()){
        std::cerr<<"ERROR in getSamples! No samples available for requested systematic: "<<systematic.name()
                 <<"\n...break\n"<<std::endl;
        exit(321);
    }
    if(m_systematicChannelSample_.at(systematic).find(channel) == m_systematicChannelSample_.at(systematic).end()){
        std::cerr<<"ERROR in getSamples! No samples available for requested channel: "<<Channel::convert(channel)
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








