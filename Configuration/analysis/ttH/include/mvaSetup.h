#ifndef mvaSetup_h
#define mvaSetup_h

#include <vector>
#include <string>

#include <TMVA/Types.h>

#include "../../common/include/sampleHelpers.h"






namespace mvaSetup{
    
    
    /// Struct keeping steering parameters for MVA training
    struct MvaConfig{
        MvaConfig(const TString& options,
                  const TString& methodAppendix,
                  const TMVA::Types::EMVA& mvaType =TMVA::Types::kBDT);
        ~MvaConfig(){}
        
        TMVA::Types::EMVA mvaType_;
        TString methodAppendix_;
        TString options_;
    };
    
    
    
    /// Struct for setting up same set of MvaConfigs in different event selections
    struct MvaSet{
        MvaSet(const std::vector<Channel::Channel>& channel,
               const TString& step,
               const std::vector<int>& v_category,
               const std::vector<MvaConfig>& v_mvaConfigCorrect,
               const std::vector<MvaConfig>& v_mvaConfigSwapped);
        ~MvaSet(){}
        
        std::vector<Channel::Channel> v_channel_;
        TString step_;
        std::vector<int> v_category_;
        std::vector<MvaConfig> v_mvaConfigCorrect_;
        std::vector<MvaConfig> v_mvaConfigSwapped_;
    };
    
    
    /// Vector of file names for MVA, separated by channels and by systematics
    typedef std::map<Systematic::Systematic, std::map<Channel::Channel, std::vector<TString> > > SystematicChannelFileNames;
    
    /// Find all files from file list for MVA, separated by channels and by systematics
    /// Flag forTraining is used to separate training and testing files
    SystematicChannelFileNames systematicChannelFileNames(const char* filelistDirectory,
                                                          const std::vector<Channel::Channel>& v_channel,
                                                          const std::vector<Systematic::Systematic>& v_systematic,
                                                          const bool forTraining =true);
    
    /// Merge training and testing trees as they are sprecified in MVA sets (for given steps/categories and channels)
    SystematicChannelFileNames mergeTrees(const char* mvaInputDir,
                                          const SystematicChannelFileNames& m_systematicChannelFileNamesTraining,
                                          const SystematicChannelFileNames& m_systematicChannelFileNamesTesting,
                                          const std::vector<MvaSet>& v_mvaSet);
    
    /// Remove MVA sets which are not valid for given channel
    std::vector<MvaSet> cleanForChannel(const std::vector<MvaSet>& v_mvaSet, const Channel::Channel& channel);
    
    
    
    /// MVA training configurations
    const MvaConfig c1("!H:!V:NTrees=600:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.45:UseRandomisedTrees=False:UseNVars=13:nCuts=1000:PruneMethod=CostComplexity:PruneStrength=-1",
                       "c1");
    const MvaConfig c2("!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.45:UseRandomisedTrees=False:UseNVars=13:nCuts=1000:PruneMethod=CostComplexity:PruneStrength=-1",
                       "c2");
    const MvaConfig c3("!H:!V:NTrees=200:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.45:UseRandomisedTrees=False:UseNVars=13:nCuts=1000:PruneMethod=CostComplexity:PruneStrength=-1",
                       "c3");
}





#endif







