#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>

#include "TString.h"

#include "MvaFactoryEventClassification.h"
#include "MvaTreeHandlerEventClassification.h"
#include "MvaTreePlotterBase.h"
#include "MvaWeights2d.h"
#include "mvaSetup.h"
#include "higgsUtils.h"
#include "Samples.h"
#include "GlobalScaleFactors.h"
#include "plotterHelpers.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/CommandLineParameters.h"
#include "../../common/include/utils.h"





/// Set data luminosity in pb-1
constexpr double Luminosity = 19712.;

/// Set relative luminosity uncertainty
constexpr double LuminosityUNCERTAINTY = 0.026;



/// The MVA input base folder
constexpr const char* MvaInputDIR = "mvaInput_mvaEvent";

/// The MVA output base folder
constexpr const char* MvaOutputDIR = "mvaOutput_mvaEvent";

/// The MVA output sub-folder for weights (1D and 2D weights)
constexpr const char* MvaWeightFileDIR = "weights";



/// Input base for the file lists containing the samples to be processed
constexpr const char* FileListBASE = "FileLists_mvaInput_mvaEvent/HistoFileList_";



/// The output file name for the control and separation power plots
constexpr const char* PlotOutputFILE = "plots.root";

/// The output file name for the MVA weights histograms (2D)
constexpr const char* WeightsOutputFILE = "weights2d.root";





void trainMvaEventClassification(const std::vector<Channel::Channel>& v_channel,
                                 const std::vector<Systematic::Systematic>& v_systematic,
                                 const std::vector<GlobalCorrection::GlobalCorrection> v_globalCorrection,
                                 const std::vector<std::string>& v_mode)
{
    // Set up scale factors
    const bool dyCorrection = std::find(v_globalCorrection.begin(), v_globalCorrection.end(), GlobalCorrection::dy) != v_globalCorrection.end();
    const bool ttbbCorrection = std::find(v_globalCorrection.begin(), v_globalCorrection.end(), GlobalCorrection::ttbb) != v_globalCorrection.end();
    const GlobalScaleFactors* globalScaleFactors = new GlobalScaleFactors(v_channel, v_systematic, Luminosity, LuminosityUNCERTAINTY, dyCorrection, ttbbCorrection);
    
    // Access all samples
    const Samples samples("FileLists_mvaInput_mvaEvent", v_channel, v_systematic, globalScaleFactors);
    
    
    
    // Set up several parameters for automated application of zillions of test trainings
    std::vector<TString> v_NTrees = {"100", "200", "400", "600", "800", "1000"};
    std::vector<TString> v_MaxDepth = {"3", "4", "5"};
    //std::vector<TString> v_BoostType = {"AdaBoost:AdaBoostBeta=0.45", "Grad"};
    std::vector<TString> v_BoostType = {"AdaBoost:AdaBoostBeta=0.45"};
    std::vector<TString> v_UseNVars = {"4", "6", "8", "10", "12"};
    std::vector<TString> v_nCuts = {"100", "1000", "5000"};
    std::vector<mvaSetup::MvaConfig> v_config;
    int counter(1);
    for(const auto& NTrees : v_NTrees){
        for(const auto& MaxDepth : v_MaxDepth){
            for(const auto& BoostType : v_BoostType){
                for(const auto& UseNVars : v_UseNVars){
                    for(const auto& nCuts : v_nCuts){
                        std::stringstream ss_counter;
                        ss_counter<<counter;
                        const TString commands = "!H:!V:NTrees="+NTrees+":nEventsMin=400:MaxDepth="+MaxDepth+":BoostType="+BoostType+":UseRandomisedTrees=True:UseNVars="+UseNVars+":nCuts="+nCuts+":PruneMethod=CostComplexity:PruneStrength=-1";
                        TString name = "d";
                        name.Append(ss_counter.str().c_str());
                        mvaSetup::MvaConfig config(commands, name);
                        v_config.push_back(config);
                        ++counter;
                    }
                }
            }
        }
    }
    
    // Define MVA sets, i.e. for which merged categories of which step to apply MVA (also separated by channels)
    // First set is for correct combinations, second for swapped combinations
    const std::vector<Channel::Channel> allChannels = {Channel::ee, Channel::emu, Channel::mumu, Channel::combined};
    std::vector<mvaSetup::MvaSet> mvaSets;
    mvaSets.push_back(mvaSetup::MvaSet(
        allChannels,
        "7",
        {0},
        {mvaSetup::c1, mvaSetup::c2, mvaSetup::c3}));
    mvaSets.push_back(mvaSetup::MvaSet(
        allChannels,
        "7",
        {1},
        {mvaSetup::c1, mvaSetup::c2, mvaSetup::c3}));
    mvaSets.push_back(mvaSetup::MvaSet(
        allChannels,
        "7",
        {2},
        {mvaSetup::c1, mvaSetup::c2, mvaSetup::c3}));
    
    const MvaFactoryEventClassification mvaFactory(MvaOutputDIR, MvaWeightFileDIR, samples);
    mvaFactory.train2(mvaSets);
    
    
    
/*    
    // Access all MVA input file names for all systematics and channels
    mvaSetup::SystematicChannelFileNames m_systematicChannelFileNamesTraining =
        mvaSetup::systematicChannelFileNames(FileListBASE, v_channel, v_systematic);
    mvaSetup::SystematicChannelFileNames m_systematicChannelFileNamesTesting =
        mvaSetup::systematicChannelFileNames(FileListBASE, v_channel, v_systematic, false);
    
    // Loop over all channels and systematics and merge trees
    mvaSetup::SystematicChannelFileNames m_systematicChannelMergedFiles = 
        mvaSetup::mergeTrees(MvaInputDIR,
                             m_systematicChannelFileNamesTraining,
                             m_systematicChannelFileNamesTesting,
                             mvaSets);
    
    // Loop over all channels and systematics and run tools
    for(const auto& systematicChannelMergedFiles : m_systematicChannelMergedFiles){
        const Systematic::Systematic& systematic = systematicChannelMergedFiles.first;
        for(const auto& channelMergedFiles : systematicChannelMergedFiles.second){
            const Channel::Channel& channel = channelMergedFiles.first;
            const TString outputFolder = common::assignFolder(MvaOutputDIR, channel, systematic);
            const TString& fileName = channelMergedFiles.second.at(0);
            
            // Print all separation power plots
            if(std::find(v_mode.begin(), v_mode.end(), "cp") != v_mode.end()){
                TString outputFile = outputFolder;
                outputFile.Append(PlotOutputFILE);
                MvaTreeHandlerEventClassification mvaTreeHandler("", {});
                mvaTreeHandler.importTrees(fileName.Data(), "training");
                MvaTreePlotterBase* mvaTreePlotter = mvaTreeHandler.setPlotter(mvaTreeHandler.stepMvaVariablesMap(), true);
                mvaTreePlotter->plotVariables(outputFile.Data());
                mvaTreePlotter->clear();
                delete mvaTreePlotter;
                mvaTreeHandler.clear();
            }
            
            // Produce MVA weights
            if(std::find(v_mode.begin(), v_mode.end(), "mva") != v_mode.end()){
                
                // Clean the MVA sets in case they are not selected for training
                std::vector<mvaSetup::MvaSet> cleanSets = mvaSetup::cleanForChannel(mvaSets, channel);
                
                // Run the MVA training for correct and swapped combinations
                const MvaFactoryEventClassification mvaFactory(outputFolder, MvaWeightFileDIR, fileName);
                mvaFactory.train(cleanSets);
                
                // Build 2D histograms of MVA weights for correct and swapped combinations
                MvaTreeHandlerEventClassification mvaTreeHandler("", {});
                mvaTreeHandler.importTrees(fileName.Data(), "training");
                TString weightsFolder(outputFolder);
                weightsFolder.Append(MvaWeightFileDIR).Append("/");
                TString outputFile = weightsFolder;
                outputFile.Append(WeightsOutputFILE);
                MvaWeights2d mvaWeights2d(mvaTreeHandler.stepMvaVariablesMap(),
                                          weightsFolder.Data(), mvaSets, true);
                mvaWeights2d.plotVariables(outputFile.Data());
                mvaTreeHandler.clear();
            }
        }
    }
*/    
    std::cout<<"MVA program successfully finished\n";
}



/// All systematics allowed as steering parameter
namespace Systematic{
    const std::vector<Type> allowedSystematics = {
        nominal,
    };
}



int main(int argc, char** argv)
{
    CLParameter<std::string> opt_channel("c", "Specify channel(s), valid: emu, ee, mumu, combined. Default: combined", false, 1, 4,
        common::makeStringCheck(Channel::convert(Channel::allowedChannelsPlotting)));
    CLParameter<std::string> opt_systematic("s", "Systematic variation - default is Nominal", false, 1, 100,
        common::makeStringCheckBegin(Systematic::convertType(Systematic::allowedSystematics)));
    CLParameter<std::string> opt_globalCorrection("g", "Specify global correction, valid: empty argument for none, Drell-Yan (dy), tt+HF (ttbb). Default: none", false, 0, 2,
        common::makeStringCheck(GlobalCorrection::convert(GlobalCorrection::allowedGlobalCorrections)));
    CLParameter<std::string> opt_mode("m", "Mode: separation plots (cp), run MVA (mva). Default is cp", false, 1, 2,
        //[](const std::string& m){return m=="" || m=="cp" || m=="mva";});
        common::makeStringCheck({"", "cp", "mva"}));
    CLAnalyser::interpretGlobal(argc, argv);
    
    // Set up channels
    std::vector<Channel::Channel> v_channel({Channel::combined});
    if(opt_channel.isSet()) v_channel = Channel::convert(opt_channel.getArguments());
    std::cout<<"Processing channels: ";
    for(auto channel : v_channel) std::cout<<Channel::convert(channel)<<" ";
    std::cout<<"\n\n";
    
    // Set up systematics
    std::vector<Systematic::Systematic> v_systematic({Systematic::nominalSystematic()});
    if(opt_systematic.isSet()) v_systematic = Systematic::setSystematics(opt_systematic.getArguments());
    std::cout<<"Processing systematics: "; 
    for(auto systematic : v_systematic) std::cout<<systematic.name()<<" ";
    std::cout<<"\n\n";
    
    // Set up global corrections
    std::vector<GlobalCorrection::GlobalCorrection> v_globalCorrection({});
    if(opt_globalCorrection.isSet()) v_globalCorrection = GlobalCorrection::convert(opt_globalCorrection.getArguments());
    std::cout<<"\n";
    std::cout<<"Using global corrections: ";
    for(auto globalCorrection : v_globalCorrection) std::cout<<GlobalCorrection::convert(globalCorrection)<<" ";
    std::cout<<"\n\n";
    
    // Set up modes
    const std::vector<std::string> v_mode = opt_mode.isSet() ? opt_mode.getArguments() : std::vector<std::string>{"cp"};
    std::cout<<"Running following modes: ";
    for(auto mode : v_mode) std::cout << mode << " ";
    std::cout<<"\n\n";
    
    trainMvaEventClassification(v_channel, v_systematic, v_globalCorrection, v_mode);
}





