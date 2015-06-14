#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

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
        this->writeYields(outputDirectory, samples, nameStepPair, plain);
        this->writeYields(outputDirectory, samples, nameStepPair, lumiCorrected);
        this->writeYields(outputDirectory, samples, nameStepPair, fullyCorrected);
    }
    
    std::cout<<"\n=== Finishing event yield table processing\n\n";
}



void EventYields::writeYields(const char* outputDirectory,
                              const Samples& samples,
                              const std::pair<TString, TString>& nameStepPair,
                              const WeightType weightType)const
{
    // Access global weights for samples if requested
    SystematicChannelFactors globalWeights;
    if(weightType != plain){
        const bool useCorrections = weightType==fullyCorrected;
        
        // Get the scale factors from the samples
        const std::pair<SystematicChannelFactors, bool> globalWeightsPair = samples.globalWeights(nameStepPair.second, useCorrections, useCorrections);
        globalWeights = globalWeightsPair.first;
        const bool anyCorrectionApplied = globalWeightsPair.second;
        
        // Do not produce corrected yields in case there are no corrections applied for this step
        if(useCorrections && !anyCorrectionApplied) return;
    }
    
    // Loop over systematics and channels
    for(auto systematicChannelSamples : samples.getSystematicChannelSamples()){
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        const auto& channelSamples(systematicChannelSamples.second);
        for(const auto& channelSample : systematicChannelSamples.second){
            const Channel::Channel& channel(channelSample.first);
            const std::vector<Sample>& v_sample(channelSamples.at(channel));
            
            // Assign histogram to sample and weight it
            std::vector<SampleHistPair> v_numhist;
            for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
                const Sample& sample = v_sample.at(iSample);
                TH1D* temp_hist = fileReader_->GetClone<TH1D>(sample.inputFile(), nameStepPair.first);
                if(weightType!=plain && sample.sampleType()!=Sample::data){
                    const double& weight = globalWeights.at(systematic).at(channel).at(iSample);
                    temp_hist->Scale(weight);
                }
                v_numhist.push_back(SampleHistPair(sample, temp_hist));
            }
            
            // Prepare output folder and text file
            std::ofstream eventFile;
            TString eventFileString = common::assignFolder(outputDirectory, channel, systematic);
            if(weightType == plain) eventFileString.Append("plain_");
            else if(weightType == fullyCorrected) eventFileString.Append("fullyCorrected_");
            else if(weightType == lumiCorrected) eventFileString.Append("lumiCorrected_");
            eventFileString.Append("events" + nameStepPair.second + ".txt");
            eventFile.open(eventFileString.Data());
            
            // Make output for tables
            double tmp_num = 0.;
            double tmp_err2 = 0.;
            double bg_num = 0.;
            double bg_err2 = 0.;
            double mc_num = 0.;
            double mc_err2 = 0.;
            std::streamsize sSize = std::cout.precision();
            for(std::vector<SampleHistPair>::const_iterator i_numhist = v_numhist.begin(); i_numhist != v_numhist.end(); ++i_numhist){
                std::vector<SampleHistPair>::const_iterator iterator = i_numhist;
                ++iterator;
                Sample::SampleType sampleType = i_numhist->first.sampleType();
                double statError(0.);
                tmp_num += i_numhist->second->IntegralAndError(0, -1, statError);
                tmp_err2 += statError*statError;
                if(i_numhist == --(v_numhist.end())){
                    const double tmp_err = std::sqrt(tmp_err2);
                    const double tmp_errRel = tmp_num>0. ? tmp_err/tmp_num*100. : 0.;
                    eventFile<<i_numhist->first.legendEntry()<<": "<<tmp_num<<"     +- "<<tmp_err<<" ("<<std::setprecision(2)<<tmp_errRel<<std::setprecision(sSize)<<"%)"<<std::endl;
                    if(sampleType != Sample::ttHbb){
                        bg_num += tmp_num;
                        bg_err2 += tmp_err2;
                    }
                    mc_num += tmp_num;
                    mc_err2 += tmp_err2;
                    tmp_num = 0.;
                    tmp_err2 = 0.;
                }
                else if(i_numhist->first.legendEntry() != iterator->first.legendEntry()){
                    const double tmp_err = std::sqrt(tmp_err2);
                    const double tmp_errRel = tmp_num>0. ? tmp_err/tmp_num*100. : 0.;
                    eventFile<<i_numhist->first.legendEntry()<<": "<<tmp_num<<"     +- "<<tmp_err<<" ("<<std::setprecision(2)<<tmp_errRel<<std::setprecision(sSize)<<"%)"<<std::endl;
                    if(sampleType!=Sample::data && sampleType!=Sample::pseudodata){
                        if(sampleType != Sample::ttHbb){
                            bg_num += tmp_num;
                            bg_err2 += tmp_err2;
                        }
                        mc_num += tmp_num;
                        mc_err2 += tmp_err2;
                    }
                    tmp_num = 0.;
                    tmp_err2 = 0.;
                }
            }
            const double bg_err = std::sqrt(bg_err2);
            const double bg_errRel = bg_num>0. ? bg_err/bg_num*100. : 0.;
            const double mc_err = std::sqrt(mc_err2);
            const double mc_errRel = mc_num>0. ? mc_err/mc_num*100. : 0.;
            eventFile<<"Total background: "<<bg_num<<"     +- "<<bg_err<<" ("<<std::setprecision(2)<<bg_errRel<<std::setprecision(sSize)<<"%)"<<std::endl;
            eventFile<<"Total MC: "<<mc_num<<"     +- "<<mc_err<<" ("<<std::setprecision(2)<<mc_errRel<<std::setprecision(sSize)<<"%)"<<std::endl;
            
            // Cleanup
            eventFile.close();
        }
    }
}















