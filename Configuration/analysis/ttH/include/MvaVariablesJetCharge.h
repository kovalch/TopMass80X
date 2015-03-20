#ifndef MvaVariablesJetCharge_h
#define MvaVariablesJetCharge_h

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









/// Class holding the input variables for one entry for MVA for jet charge,
/// i.e. one entry of each quantity per selected jet
class MvaVariablesJetCharge : public MvaVariablesBase{
    
public:
    
    /// Empty constructor
    MvaVariablesJetCharge();
    
    /// Constructor setting up input variables from physics objects
    MvaVariablesJetCharge(const LV& lepton, const LV& antiLepton,
                        const LV& bJet, const LV& antiBJet,
                        const double& bJetBtagDiscriminator, const double& antiBJetBtagDiscriminator,
                        const double& jetChargeDiff,
                        const LV& jetRecoil, const LV& met,
                        const bool bQuarkRecoJetMatched,
                        const bool correctCombination, const bool swappedCombination,
                        const bool lastInEvent,
                        const double& eventWeight);
    
    /// Destructor
    ~MvaVariablesJetCharge(){}
    
    /// Fill the MVA input structs for all jet combinations of one event
    static std::vector<MvaVariablesBase*> fillVariables(const tth::RecoObjectIndices& recoObjectIndices,
                                                        const tth::GenObjectIndices& genObjectIndices,
                                                        const RecoObjects& recoObjects,
                                                        const double& eventWeight);
    
    
    
    // The variables needed for MVA
    
    // FIXME: Describe each variable in doxygen
    // FIXME: Review the names
    /// Variables for MVA
    // FIXME: Why are these floats?
    MvaVariableFloat trueBJetId_;
    MvaVariableFloat thereIsALeadingLepton_;
    MvaVariableFloat thereIsALeadingMuon_;
    MvaVariableFloat thereIsASecondaryVertex_;

    MvaVariableFloat longChargeJet_;
    MvaVariableFloat relChargeJet_;
    MvaVariableFloat leadingTrackPtWeightedCharge_;
    MvaVariableFloat subleadingTrackPtWeightedCharge_;
    MvaVariableFloat thirdleadingTrackPtWeightedCharge_;
    MvaVariableFloat leadingMuonPtWeightedCharge_;
    MvaVariableFloat leadingElectronPtWeightedCharge_;
    MvaVariableFloat trackNumberWeightedJetPt_;
    MvaVariableFloat chargeWeightedTrackId_;
    MvaVariableFloat svChargeWeightedFlightDistance_;
    MvaVariableFloat secondaryVertexCharge_;
    MvaVariableFloat ipSignificanceLeadingTrack_;
    
    
    
private:
    
    // The names associated to the variables
    static constexpr const char* name_trueBJetId_ = "trueBJetId";
    static constexpr const char* name_thereIsALeadingLepton_ = "thereIsALeadingLepton";
    static constexpr const char* name_thereIsALeadingMuon_ = "thereIsALeadingMuon";
    static constexpr const char* name_thereIsASecondaryVertex_ = "thereIsASecondaryVertex";
    
    static constexpr const char* name_longChargeJet_ = "longChargeJet";
    static constexpr const char* name_relChargeJet_ = "relChargeJet";
    static constexpr const char* name_leadingTrackPtWeightedCharge_ = "leadingTrackPtWeightedCharge";
    static constexpr const char* name_subleadingTrackPtWeightedCharge_ = "subleadingTrackPtWeightedCharge";
    static constexpr const char* name_thirdleadingTrackPtWeightedCharge_ = "thirdleadingTrackPtWeightedCharge";
    static constexpr const char* name_leadingMuonPtWeightedCharge_ = "leadingMuonPtWeightedCharge";
    static constexpr const char* name_leadingElectronPtWeightedCharge_ = "leadingElectronPtWeightedCharge";
    static constexpr const char* name_trackNumberWeightedJetPt_ = "trackNumberWeightedJetPt";
    static constexpr const char* name_chargeWeightedTrackId_ = "chargeWeightedTrackId";
    static constexpr const char* name_svChargeWeightedFlightDistance_ = "svChargeWeightedFlightDistance";
    static constexpr const char* name_secondaryVertexCharge_ = "secondaryVertexCharge";
    static constexpr const char* name_ipSignificanceLeadingTrack_ = "ipSignificanceLeadingTrack";
};




/// Class holding the input variables for all jets of one event for jet charge MVA,
/// potentially together with estimated MVA weights
class MvaVariablesJetChargePerEvent{
    
public:
    
    /// Constructor setting up mvaInputVariables
    MvaVariablesJetChargePerEvent(const std::vector<MvaVariablesBase*>& v_mvaVariables);
    
    /// Constructor setting up mvaInputVariables together with MVA weights for correct and for swapped combinations
    MvaVariablesJetChargePerEvent(const std::vector<MvaVariablesBase*>& v_mvaVariables,
                                 const std::map<std::string, std::vector<float> >& m_weight);
    
    /// Destructor
    ~MvaVariablesJetChargePerEvent(){}
    
    
    
    /// Index of dijet combination with maximum correct MVA weight for given config name
    size_t maxWeightIndex(const std::string& mvaConfigName)const;
    
    /// Returns the maximum correct MVA weight per event for given config name
    float maxWeight(const std::string& mvaConfigName)const;
    
    /// Get the vector of MVA variables
    std::vector<MvaVariablesBase*> variables()const;
    
    /// Get the vector of MVA correct weights for given config name
    std::vector<float> mvaWeights(const std::string& mvaConfigName)const;
    
    /// Get the vector of MVA correct weights
    std::map<std::string, std::vector<float> > mvaWeightsMap()const;
    
    
    
private:
    
    /// Vector containing MVA input variables for all dijet pairs per event
    std::vector<MvaVariablesBase*> v_mvaVariables_;
    
    /// Map MVA config name to the vector containing correct MVA weights for all dijet pairs per event
    std::map<std::string, std::vector<float> > m_weight_;
};





#endif







