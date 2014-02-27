#ifndef MvaVariablesTopJets_h
#define MvaVariablesTopJets_h

#include <map>
#include <vector>
#include <string>

#include "MvaVariablesBase.h"
#include "analysisStructsFwd.h"
#include "../../common/include/classesFwd.h"

class RecoObjects;
namespace tth{
    class RecoObjectIndices;
    class GenObjectIndices;
}









/// Class holding the input variables for one entry for MVA for top system jet identification,
/// i.e. one entry of each quantity per selected jet combination
class MvaVariablesTopJets : public MvaVariablesBase{
    
public:
    
    /// Empty constructor
    MvaVariablesTopJets();
    
    /// Constructor setting up input variables from physics objects
    MvaVariablesTopJets(const LV& lepton, const LV& antiLepton,
                        const LV& bJet, const LV& antiBJet,
                        const double& bJetBtagDiscriminator, const double& antiBJetBtagDiscriminator,
                        const double& jetChargeDiff,
                        const LV& jetRecoil, const LV& met,
                        const bool bQuarkRecoJetMatched,
                        const bool correctCombination, const bool swappedCombination,
                        const bool lastInEvent,
                        const double& eventWeight);
    
    /// Destructor
    ~MvaVariablesTopJets(){}
    
    /// Fill the MVA input structs for all jet combinations of one event
    static std::vector<MvaVariablesBase*> fillVariables(const tth::RecoObjectIndices& recoObjectIndices,
                                                        const tth::GenObjectIndices& genObjectIndices,
                                                        const RecoObjects& recoObjects,
                                                        const double& eventWeight);
    
    /// Calculate the jet recoil for a given jet pair, i.e. vector sum of all jets except selected jet pair
    static VLV recoilForJetPairs(const tth::IndexPairs& jetIndexPairs,
                                 const std::vector<int>& jetIndices,
                                 const VLV& jets);
    
    
    
    // The variables needed for MVA
    
    /// Is it the last dijet combination in the event
    MvaVariableInt lastInEvent_;
    /// Could b quark and anti-b quark be matched to reco jets
    MvaVariableInt bQuarkRecoJetMatched_;
    /// Is it the true correct jet combination
    MvaVariableInt correctCombination_;
    /// Is it the true but swapped jet combination
    MvaVariableInt swappedCombination_;
    
    /// Difference of the jet charges for (anti-b jet - b jet), i.e. it is within [0,2]
    MvaVariableFloat jetChargeDiff_;
    // FIXME: describe each variable in doxygen
    /// Variables for MVA
    MvaVariableFloat meanDeltaPhi_b_met_;
    MvaVariableFloat massDiff_recoil_bbbar_;
    MvaVariableFloat pt_b_antiLepton_;
    MvaVariableFloat pt_antiB_lepton_;
    MvaVariableFloat deltaR_b_antiLepton_;
    MvaVariableFloat deltaR_antiB_lepton_;
    MvaVariableFloat btagDiscriminatorSum_;
    MvaVariableFloat deltaPhi_antiBLepton_bAntiLepton_;
    MvaVariableFloat massDiff_fullBLepton_bbbar_;
    MvaVariableFloat meanMt_b_met_;
    MvaVariableFloat massSum_antiBLepton_bAntiLepton_;
    MvaVariableFloat massDiff_antiBLepton_bAntiLepton_;
    
    
    
private:
    
    // The names associated to the variables
    
    static constexpr const char* name_lastInEvent_ = "lastInEvent";
    static constexpr const char* name_bQuarkRecoJetMatched_ = "bQuarkRecoJetMatched";
    static constexpr const char* name_correctCombination_ = "correctCombination";
    static constexpr const char* name_swappedCombination_ = "swappedCombination";
    
    static constexpr const char* name_jetChargeDiff_ = "jetChargeDiff";
    static constexpr const char* name_meanDeltaPhi_b_met_ = "meanDeltaPhi_b_met";
    static constexpr const char* name_massDiff_recoil_bbbar_ = "massDiff_recoil_bbbar";
    static constexpr const char* name_pt_b_antiLepton_ = "pt_b_antiLepton";
    static constexpr const char* name_pt_antiB_lepton_ = "pt_antiB_lepton";
    static constexpr const char* name_deltaR_b_antiLepton_ = "deltaR_b_antiLepton";
    static constexpr const char* name_deltaR_antiB_lepton_ = "deltaR_antiB_lepton";
    static constexpr const char* name_btagDiscriminatorSum_ = "btagDiscriminatorSum";
    static constexpr const char* name_deltaPhi_antiBLepton_bAntiLepton_ = "deltaPhi_antiBLepton_bAntiLepton";
    static constexpr const char* name_massDiff_fullBLepton_bbbar_ = "massDiff_fullBLepton_bbbar";
    static constexpr const char* name_meanMt_b_met_ = "meanMt_b_met";
    static constexpr const char* name_massSum_antiBLepton_bAntiLepton_ = "massSum_antiBLepton_bAntiLepton";
    static constexpr const char* name_massDiff_antiBLepton_bAntiLepton_ = "massDiff_antiBLepton_bAntiLepton";
};




/// Class holding the input variables for all jet combinations of one event for MVA,
/// potentially together with estimated MVA weights for correct and swapped combinations
class MvaVariablesTopJetsPerEvent{
    
public:
    
    /// Constructor setting up mvaInputVariables
    MvaVariablesTopJetsPerEvent(const std::vector<MvaVariablesBase*>& v_mvaVariables);
    
    /// Constructor setting up mvaInputVariables together with MVA weights for correct and for swapped combinations
    MvaVariablesTopJetsPerEvent(const std::vector<MvaVariablesBase*>& v_mvaVariables,
                                const std::map<std::string, std::vector<float> >& m_weightCorrect,
                                const std::map<std::string, std::vector<float> >& m_weightSwapped,
                                const std::map<std::string, std::map<std::string, std::vector<float> > >& m_weightCombined);
    
    /// Destructor
    ~MvaVariablesTopJetsPerEvent(){}
    
    
    
    /// Index of dijet combination with maximum correct MVA weight for given config name
    size_t maxWeightCorrectIndex(const std::string& mvaConfigName)const;
    
    /// Index of dijet combination with maximum swapped MVA weight for given config name
    size_t maxWeightSwappedIndex(const std::string& mvaConfigName)const;
    
    /// Index of dijet combination with maximum combined MVA weight for given config names
    size_t maxWeightCombinedIndex(const std::string& mvaConfigNameCorrect, const std::string& mvaConfigNameSwapped)const;
    
    /// Returns the maximum correct MVA weight per event for given config name
    float maxWeightCorrect(const std::string& mvaConfigName)const;
    
    /// Returns the maximum swapped MVA weight per event for given config name
    float maxWeightSwapped(const std::string& mvaConfigName)const;
    
    /// Returns the maximum combined MVA weight per event for given config names
    float maxWeightCombined(const std::string& mvaConfigNameCorrect, const std::string& mvaConfigNameSwapped)const;
    
    /// Does the same dijet combination have the maximum correct and the maximum swapped weight per event for given config names
    bool isSameMaxCombination(const std::string& mvaConfigNameCorrect, const std::string& mvaConfigNameSwapped)const;
    
    /// Get the vector of MVA variables
    std::vector<MvaVariablesBase*> variables()const;
    
    /// Get the vector of MVA correct weights for given config name
    std::vector<float> mvaWeightsCorrect(const std::string& mvaConfigName)const;
    
    /// Get the vector of MVA swapped weights for given config name
    std::vector<float> mvaWeightsSwapped(const std::string& mvaConfigName)const;
    
    /// Get the vector of MVA combined weights for given config names
    std::vector<float> mvaWeightsCombined(const std::string& mvaConfigNameCorrect, const std::string& mvaConfigNameSwapped)const;
    
    /// Get the vector of MVA correct weights
    std::map<std::string, std::vector<float> > mvaWeightsCorrectMap()const;
    
    /// Get the vector of MVA swapped weights
    std::map<std::string, std::vector<float> > mvaWeightsSwappedMap()const;
    
    /// Get the vector of MVA combined weights
    std::map<std::string, std::map<std::string, std::vector<float> > > mvaWeightsCombinedMap()const;
    
    
    
private:
    
    /// Vector containing MVA input variables for all dijet pairs per event
    std::vector<MvaVariablesBase*> v_mvaVariables_;
    
    /// Map MVA config name to the vector containing correct MVA weights for all dijet pairs per event
    std::map<std::string, std::vector<float> > m_weightCorrect_;
    
    /// Map MVA config name to the vector containing swapped MVA weights for all dijet pairs per event
    std::map<std::string, std::vector<float> > m_weightSwapped_;
    
    /// Map MVA config name to the vector containing swapped MVA weights for all dijet pairs per event
    std::map<std::string, std::map<std::string, std::vector<float> > > m_weightCombined_;
};





#endif







