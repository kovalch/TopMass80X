#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>

#include "TString.h"

#include "AnalysisConfig.h"
#include "MvaFactoryEventClassification.h"
#include "mvaSetup.h"
#include "Samples.h"
#include "GlobalScaleFactors.h"
#include "plotterHelpers.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/CommandLineParameters.h"
#include "../../common/include/utils.h"





void trainMvaEventClassification(const std::vector<Channel::Channel>& v_channel,
                                 const std::vector<Systematic::Systematic>& v_systematic,
                                 const std::vector<GlobalCorrection::GlobalCorrection> v_globalCorrection)
{
    // Read analysis config from text file
    const AnalysisConfig analysisConfig;
    
    // Set up scale factors
    const bool dyCorrection = std::find(v_globalCorrection.begin(), v_globalCorrection.end(), GlobalCorrection::dy) != v_globalCorrection.end();
    const bool ttbbCorrection = std::find(v_globalCorrection.begin(), v_globalCorrection.end(), GlobalCorrection::ttbb) != v_globalCorrection.end();
    const GlobalScaleFactors* globalScaleFactors = new GlobalScaleFactors(analysisConfig, v_channel, v_systematic, dyCorrection, ttbbCorrection);
    
    // Access all samples
    const Samples samples("FileLists_mvaInput_mvaEvent", analysisConfig, v_channel, v_systematic, globalScaleFactors);
    
    // MVA training configurations
    const mvaSetup::MvaConfig d1("!H:!V:NTrees=300:nEventsMin=40:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.15:UseRandomisedTrees=False:UseNVars=13:nCuts=30:PruneMethod=CostComplexity:PruneStrength=-1",
                                 "d1");
    
    const mvaSetup::MvaConfig e1("!H:!V:NTrees=100:nEventsMin=40:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.10:UseRandomisedTrees=False:UseNVars=13:nCuts=30:PruneMethod=CostComplexity:PruneStrength=-1",
                                 "e1");
    const mvaSetup::MvaConfig e2("!H:!V:NTrees=100:nEventsMin=40:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.15:UseRandomisedTrees=False:UseNVars=13:nCuts=30:PruneMethod=CostComplexity:PruneStrength=-1",
                                 "e2");
    const mvaSetup::MvaConfig e3("!H:!V:NTrees=100:nEventsMin=40:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.20:UseRandomisedTrees=False:UseNVars=13:nCuts=30:PruneMethod=CostComplexity:PruneStrength=-1",
                                 "e3");
    const mvaSetup::MvaConfig e4("!H:!V:NTrees=200:nEventsMin=40:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.10:UseRandomisedTrees=False:UseNVars=13:nCuts=30:PruneMethod=CostComplexity:PruneStrength=-1",
                                 "e4");
    const mvaSetup::MvaConfig e5("!H:!V:NTrees=200:nEventsMin=40:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.15:UseRandomisedTrees=False:UseNVars=13:nCuts=30:PruneMethod=CostComplexity:PruneStrength=-1",
                                 "e5");
    const mvaSetup::MvaConfig e6("!H:!V:NTrees=200:nEventsMin=40:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.20:UseRandomisedTrees=False:UseNVars=13:nCuts=30:PruneMethod=CostComplexity:PruneStrength=-1",
                                 "e6");
    const mvaSetup::MvaConfig e7("!H:!V:NTrees=300:nEventsMin=40:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.10:UseRandomisedTrees=False:UseNVars=13:nCuts=30:PruneMethod=CostComplexity:PruneStrength=-1",
                                 "e7");
    const mvaSetup::MvaConfig e8("!H:!V:NTrees=300:nEventsMin=40:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.15:UseRandomisedTrees=False:UseNVars=13:nCuts=30:PruneMethod=CostComplexity:PruneStrength=-1",
                                 "e8");
    const mvaSetup::MvaConfig e9("!H:!V:NTrees=300:nEventsMin=40:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.20:UseRandomisedTrees=False:UseNVars=13:nCuts=30:PruneMethod=CostComplexity:PruneStrength=-1",
                                 "e9");
    
    const mvaSetup::MvaConfig f1("!H:!V:NTrees=200:nEventsMin=40:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.20:UseRandomisedTrees=False:UseNVars=13:nCuts=30:PruneMethod=CostComplexity:PruneStrength=-1",
                                 "f1");
    
    const mvaSetup::MvaConfig g1("!H:!V:NTrees=100:nEventsMin=40:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.10:UseRandomisedTrees=False:UseNVars=13:nCuts=30:PruneMethod=CostComplexity:PruneStrength=-1",
                                 "g1");
    const mvaSetup::MvaConfig g2("!H:!V:NTrees=100:nEventsMin=40:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.15:UseRandomisedTrees=False:UseNVars=13:nCuts=30:PruneMethod=CostComplexity:PruneStrength=-1",
                                 "g2");
    
    const mvaSetup::MvaConfig h1("!H:!V:NTrees=600:nEventsMin=40:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.45:UseRandomisedTrees=False:UseNVars=13:nCuts=1000:PruneMethod=CostComplexity:PruneStrength=-1",
                                 "h1");
    const mvaSetup::MvaConfig h2("!H:!V:NTrees=400:nEventsMin=40:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.45:UseRandomisedTrees=False:UseNVars=13:nCuts=1000:PruneMethod=CostComplexity:PruneStrength=-1",
                                 "h2");
    
    
    // Define MVA sets, i.e. for which (merged) categories of which step to apply MVA, also separated by channels
    const std::vector<Channel::Channel> allChannels = {Channel::ee, Channel::emu, Channel::mumu, Channel::combined};
    std::vector<mvaSetup::MvaSet> mvaSets;
    mvaSets.push_back(mvaSetup::MvaSet(
        allChannels,
        "7",
        {0},
        {d1}));
    mvaSets.push_back(mvaSetup::MvaSet(
        allChannels,
        "7",
        {1},
        {e1, e2, e3, e4, e5, e6, e7, e8, e9}));
    mvaSets.push_back(mvaSetup::MvaSet(
        allChannels,
        "7",
        {2},
        {f1}));
    mvaSets.push_back(mvaSetup::MvaSet(
        allChannels,
        "7",
        {3},
        {g1, g2}));
//     mvaSets.push_back(mvaSetup::MvaSet(
//         allChannels,
//         "7",
//         {0, 1, 2, 3},
//         {h1, h2}));
    
    // Set up and run factory
    const MvaFactoryEventClassification mvaFactory("mvaOutput_mvaEvent", "weights", samples);
    mvaFactory.train2(mvaSets);
    
    
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
    
    trainMvaEventClassification(v_channel, v_systematic, v_globalCorrection);
}





