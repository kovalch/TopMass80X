#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <algorithm>

#include <TMVA/Types.h>
#include <TString.h>
#include <TList.h>
#include <TFile.h>
#include <TTree.h>

#include "../../common/include/sampleHelpers.h"
#include "mvaSetup.h"
#include "higgsUtils.h"







mvaSetup::MvaConfig::MvaConfig(const TString& options,
                               const TString& methodAppendix,
                               const TMVA::Types::EMVA& mvaType):
mvaType_(mvaType),
methodAppendix_(methodAppendix),
options_(options)
{}



mvaSetup::MvaSet::MvaSet(const std::vector<Channel::Channel>& v_channel,
                         const TString& step,
                         const std::vector<int>& v_category,
                         const std::vector<mvaSetup::MvaConfig>& v_mvaConfigCorrect,
                         const std::vector<mvaSetup::MvaConfig>& v_mvaConfigSwapped):
v_channel_(v_channel),
step_(step),
v_category_(v_category),
v_mvaConfigCorrect_(v_mvaConfigCorrect),
v_mvaConfigSwapped_(v_mvaConfigSwapped)
{
    for(const int category : v_category_){
        if(category<0){
            std::cerr<<"ERROR in constructor of MvaSet! jetCategories contain value <0: "<<category
                     <<"\n...break\n"<<std::endl;
            exit(136);
        }
    }
}



std::vector<mvaSetup::MvaSet> mvaSetup::cleanForChannel(const std::vector<mvaSetup::MvaSet>& v_mvaSet, const Channel::Channel& channel)
{
    std::vector<MvaSet> result;
    for(const MvaSet& mvaSet : v_mvaSet){
        // Check if MVA set is valid for this channel
        if(std::find(mvaSet.v_channel_.begin(), mvaSet.v_channel_.end(), channel) == mvaSet.v_channel_.end()) continue;
        
        result.push_back(mvaSet);
    }
    return result;
}



mvaSetup::SystematicChannelFileNames mvaSetup::mergeTrees(const char* mvaInputDir,
                                                          const mvaSetup::SystematicChannelFileNames& m_systematicChannelFileNamesTraining,
                                                          const mvaSetup::SystematicChannelFileNames& m_systematicChannelFileNamesTesting,
                                                          const std::vector<mvaSetup::MvaSet>& v_mvaSet)
{
    SystematicChannelFileNames result;
    
//     // Find all trees of all steps/categories containing MVA input variables from first input file
//     std::cout<<"The following steps exist in input files:\n";
//     const std::vector<std::pair<TString, TString> > v_nameStepPair =
//         tth::nameStepPairs(m_systematicChannelFileNamesTraining.begin()->second.begin()->second.at(0),
//                            "mvaInputTopJets_step");
    
    // Loop over all channels and systematics
    for(const auto& systematicChannelFileNamesTraining : m_systematicChannelFileNamesTraining){
        const Systematic::Systematic& systematic = systematicChannelFileNamesTraining.first;
        for(const auto& channelFileNamesTraining : systematicChannelFileNamesTraining.second){
            const Channel::Channel& channel = channelFileNamesTraining.first;
            const std::vector<TString>& v_fileNameTraining = channelFileNamesTraining.second;
            const std::vector<TString>& v_fileNameTesting = m_systematicChannelFileNamesTesting.at(systematic).at(channel);
            std::cout<<"\nProcessing (Channel, Systematic): "<<Channel::convert(channel)<<" , "<<systematic.name()<<"\n\n";
            
            // Access sets which are valid for the channel
            const std::vector<MvaSet> v_cleanSet = cleanForChannel(v_mvaSet, channel);
            std::cout<<"Merge trees for the following MVA sets (step, merged categories):\n";
            for(const MvaSet& mvaSet : v_cleanSet){
                std::cout<<"\t"<<mvaSet.step_<<"\t\t";
                for(const auto& category : mvaSet.v_category_) std::cout<<category<<" ";
                std::cout<<std::endl;
            }
            
            // Open the input files and access the MVA input training trees
            std::map<TString, TList*> m_stepListTraining;
            for(const auto& fileName : v_fileNameTraining){
                std::cout<<"File for training: "<<fileName<<std::endl;
                // FIXME: need to check whether input file really exists
                TFile* inputFile(0);
                inputFile = TFile::Open(fileName);
                for(const MvaSet& mvaSet : v_cleanSet){
                    const TString mergedStepName = tth::stepName(mvaSet.step_, mvaSet.v_category_);
                    // Clone MVA set to handle steps without categories
                    MvaSet mvaSetClone = mvaSet;
                    if(!mvaSetClone.v_category_.size()) mvaSetClone.v_category_ = {-1};
                    for(const int& category : mvaSetClone.v_category_){
                        const TString singleStepName = tth::stepName(mvaSet.step_, category);
                        //std::cout<<"Step name (merged, single): "<<mergedStepName<<" , "<<singleStepName<<"\n";
                        // FIXME: need to check whether input tree really exists
                        TTree* inputTree = (TTree*)inputFile->Get("topJets_mvaVariables"+singleStepName);
                        if(m_stepListTraining.find(mergedStepName) == m_stepListTraining.end()){
                            m_stepListTraining[mergedStepName] = new TList;
                        }
                        m_stepListTraining.at(mergedStepName)->Add(inputTree);
                    }
                }
                
            }
            std::cout<<std::endl;
            
            // Open the input files and access the MVA input testing trees
            std::map<TString, TList*> m_stepListTesting;
            for(const auto& fileName : v_fileNameTesting){
                std::cout<<"File for testing: "<<fileName<<std::endl;
                // FIXME: need to check whether input file really exists
                TFile* inputFile(0);
                inputFile = TFile::Open(fileName);
                for(const MvaSet& mvaSet : v_cleanSet){
                    const TString mergedStepName = tth::stepName(mvaSet.step_, mvaSet.v_category_);
                    // Clone MVA set to handle steps without categories
                    MvaSet mvaSetClone = mvaSet;
                    if(!mvaSetClone.v_category_.size()) mvaSetClone.v_category_ = {-1};
                    for(const int& category : mvaSetClone.v_category_){
                        const TString singleStepName = tth::stepName(mvaSet.step_, category);
                        //std::cout<<"Step name (merged, single): "<<mergedStepName<<" , "<<singleStepName<<"\n";
                        // FIXME: need to check whether input tree really exists
                        TTree* inputTree = (TTree*)inputFile->Get("topJets_mvaVariables"+singleStepName);
                        
                        if(m_stepListTesting.find(mergedStepName) == m_stepListTesting.end()){
                            m_stepListTesting[mergedStepName] = new TList;
                        }
                        m_stepListTesting.at(mergedStepName)->Add(inputTree);
                    }
                }
            }
            std::cout<<std::endl;
            
            
            
            
            // Write merged trees to file for given channel and systematic
            TString mergedTreesFileName = common::assignFolder(mvaInputDir, channel, systematic);
            mergedTreesFileName.Append("/");
            mergedTreesFileName.Append("mergedTrees.root");
            TFile* mergedTreesFile = new TFile(mergedTreesFileName, "RECREATE");
            for(const auto& stepList : m_stepListTraining){
                TTree* tree = TTree::MergeTrees(stepList.second);
                tree->SetName("training"+stepList.first);
                tree->Write();
            }
            for(const auto& stepList : m_stepListTesting){
                TTree* tree = TTree::MergeTrees(stepList.second);
                tree->SetName("testing"+stepList.first);
                tree->Write();
            }
            mergedTreesFile->Close();
            result[systematic][channel].push_back(mergedTreesFileName);
        }
    }
    
    return result;
}











// FIXME: for now systematic is not used to study systematic variations which modify Higgs samples,
// FIXME: but running on all Higgs masses
namespace Systematic{
    const std::vector<Type> allowedSystematics = {
        nominal, mH110, mH115, mH120, mH1225, mH1275, mH130, mH135, mH140
    };
}



mvaSetup::SystematicChannelFileNames mvaSetup::systematicChannelFileNames(const char* filelistDirectory,
                                                                          const std::vector<Channel::Channel>& v_channel,
                                                                          const std::vector<Systematic::Systematic>& v_systematic,
                                                                          const bool forTraining)
{
    SystematicChannelFileNames m_systematicChannelFileNames;
    
    for(const auto& systematic : v_systematic){
        std::map<Channel::Channel, std::vector<TString> >& m_channelFileNames = m_systematicChannelFileNames[systematic];
        for(const auto& channel : v_channel){
            std::vector<TString>& v_filename = m_channelFileNames[channel];
            // FIXME: for now systematic is not used to study systematic variations which modify Higgs samples,
            // FIXME: but running on all Higgs masses
            for(const auto& systematicMass : Systematic::allowedSystematicsAnalysis(Systematic::allowedSystematics)){
                // Read the full input filenames from the FileList
                if(forTraining) v_filename = common::readFilelist(filelistDirectory, channel, systematicMass, {"ttbarH", "inclusiveBbbar"});
                else v_filename = common::readFilelist(filelistDirectory, channel, systematicMass, {"ttbarH", "tobbbar"});
            }
        }
    }
    
    return m_systematicChannelFileNames;
}







