#include <iostream>
#include <cstdlib>

#include <Math/VectorUtil.h>

#include "MvaVariablesJetCharge.h"
#include "analysisStructs.h"
#include "../../common/include/classes.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/analysisObjectStructs.h"





// ----------------------------------------- Methods for MvaVariablesJetCharge -----------------------------------------------------------



MvaVariablesJetCharge::MvaVariablesJetCharge():
MvaVariablesBase(),
trueBJetId_(MvaVariableFloat(name_trueBJetId_)),
thereIsALeadingLepton_(MvaVariableFloat(name_thereIsALeadingLepton_)),
thereIsALeadingMuon_(MvaVariableFloat(name_thereIsALeadingMuon_)),
thereIsASecondaryVertex_(MvaVariableFloat(name_thereIsASecondaryVertex_)),
longChargeJet_(MvaVariableFloat(name_longChargeJet_)),
relChargeJet_(MvaVariableFloat(name_relChargeJet_)),
leadingTrackPtWeightedCharge_(MvaVariableFloat(name_leadingTrackPtWeightedCharge_)),
subleadingTrackPtWeightedCharge_(MvaVariableFloat(name_subleadingTrackPtWeightedCharge_)),
thirdleadingTrackPtWeightedCharge_(MvaVariableFloat(name_thirdleadingTrackPtWeightedCharge_)),
leadingMuonPtWeightedCharge_(MvaVariableFloat(name_leadingMuonPtWeightedCharge_)),
leadingElectronPtWeightedCharge_(MvaVariableFloat(name_leadingElectronPtWeightedCharge_)),
trackNumberWeightedJetPt_(MvaVariableFloat(name_trackNumberWeightedJetPt_)),
chargeWeightedTrackId_(MvaVariableFloat(name_chargeWeightedTrackId_)),
svChargeWeightedFlightDistance_(MvaVariableFloat(name_svChargeWeightedFlightDistance_)),
secondaryVertexCharge_(MvaVariableFloat(name_secondaryVertexCharge_)),
ipSignificanceLeadingTrack_(MvaVariableFloat(name_ipSignificanceLeadingTrack_))
{}



MvaVariablesJetCharge::MvaVariablesJetCharge(const LV& lepton, const LV& antiLepton,
                                         const LV& bJet, const LV& antiBJet,
                                         const double& bjetBtagDiscriminator, const double& antiBjetBtagDiscriminator,
                                         const double& jetChargeDiff,
                                         const LV& jetRecoil, const LV& met,
                                         const bool bQuarkRecoJetMatched,
                                         const bool correctCombination, const bool swappedCombination,
                                         const bool lastInEvent,
                                         const double& eventWeight):
MvaVariablesBase(eventWeight),
trueBJetId_(MvaVariableFloat(name_trueBJetId_)),
thereIsALeadingLepton_(MvaVariableFloat(name_thereIsALeadingLepton_)),
thereIsALeadingMuon_(MvaVariableFloat(name_thereIsALeadingMuon_)),
thereIsASecondaryVertex_(MvaVariableFloat(name_thereIsASecondaryVertex_)),
longChargeJet_(MvaVariableFloat(name_longChargeJet_)),
relChargeJet_(MvaVariableFloat(name_relChargeJet_)),
leadingTrackPtWeightedCharge_(MvaVariableFloat(name_leadingTrackPtWeightedCharge_)),
subleadingTrackPtWeightedCharge_(MvaVariableFloat(name_subleadingTrackPtWeightedCharge_)),
thirdleadingTrackPtWeightedCharge_(MvaVariableFloat(name_thirdleadingTrackPtWeightedCharge_)),
leadingMuonPtWeightedCharge_(MvaVariableFloat(name_leadingMuonPtWeightedCharge_)),
leadingElectronPtWeightedCharge_(MvaVariableFloat(name_leadingElectronPtWeightedCharge_)),
trackNumberWeightedJetPt_(MvaVariableFloat(name_trackNumberWeightedJetPt_)),
chargeWeightedTrackId_(MvaVariableFloat(name_chargeWeightedTrackId_)),
svChargeWeightedFlightDistance_(MvaVariableFloat(name_svChargeWeightedFlightDistance_)),
secondaryVertexCharge_(MvaVariableFloat(name_secondaryVertexCharge_)),
ipSignificanceLeadingTrack_(MvaVariableFloat(name_ipSignificanceLeadingTrack_))
{
    // FIXME: calculate variables
}



std::vector<MvaVariablesBase*> MvaVariablesJetCharge::fillVariables(const tth::RecoObjectIndices& recoObjectIndices,
                                                                    const tth::GenObjectIndices& genObjectIndices,
                                                                    const RecoObjects& recoObjects,
                                                                    const double& eventWeight)
{
    std::vector<MvaVariablesBase*> result;
    
    const std::vector<int>& jetIndices = recoObjectIndices.jetIndices_;
    
    
    // Loop over all jet pairs
    for(int index : jetIndices){
        // Is it the last entry of an event
        const bool lastInEvent = index == static_cast<int>(jetIndices.size()-1);
        
        MvaVariablesBase* mvaVariables = new MvaVariablesJetCharge();
        
        result.push_back(mvaVariables);
    }
    
    return result;
}







// ----------------------------------------- Methods for MvaVariablesJetChargePerEvent -----------------------------------------------------------



MvaVariablesJetChargePerEvent::MvaVariablesJetChargePerEvent(const std::vector<MvaVariablesBase*>& v_mvaVariables):
v_mvaVariables_(v_mvaVariables)
{}



MvaVariablesJetChargePerEvent::MvaVariablesJetChargePerEvent(const std::vector<MvaVariablesBase*>& v_mvaVariables,
                                                         const std::map<std::string, std::vector<float> >& m_mvaWeight):
v_mvaVariables_(v_mvaVariables),
m_weight_(m_mvaWeight)
{
    for(const auto& weights : m_weight_){
        if(weights.second.size() != v_mvaVariables_.size()){
            std::cerr<<"ERROR in constructor of MvaVariablesJetChargePerEvent! Vector sizes do not match for weights (variables, weights): "
                     <<v_mvaVariables.size()<<" , "<<weights.second.size()<<"\n...break\n"<<std::endl;
            exit(47);
        }
    }
}



size_t MvaVariablesJetChargePerEvent::maxWeightIndex(const std::string& mvaConfigName)const
{
    const std::vector<float>& v_weight = this->mvaWeights(mvaConfigName);
    return common::extremumIndex(v_weight);
}



float MvaVariablesJetChargePerEvent::maxWeight(const std::string& mvaConfigName)const
{
    const size_t maxWeightIndex = this->maxWeightIndex(mvaConfigName);
    return m_weight_.at(mvaConfigName).at(maxWeightIndex);
}



std::vector<MvaVariablesBase*> MvaVariablesJetChargePerEvent::variables()const
{
    return v_mvaVariables_;
}



std::vector<float> MvaVariablesJetChargePerEvent::mvaWeights(const std::string& mvaConfigName)const
{
    if(m_weight_.find(mvaConfigName) == m_weight_.end()){
        std::cerr<<"ERROR in mvaWeights()! No weights found for mvaConfigName: "
                 <<mvaConfigName<<"\n...break\n"<<std::endl;
        exit(50);
    }
    
    return m_weight_.at(mvaConfigName);
}



std::map<std::string, std::vector<float> > MvaVariablesJetChargePerEvent::mvaWeightsMap()const
{
    return m_weight_;
}







