#include <algorithm>

#include <Math/VectorUtil.h>

#include "MvaVariablesEventClassification.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/classes.h"




MvaVariablesEventClassification::MvaVariablesEventClassification():
MvaVariablesBase(),
multiplicity_jets_(MvaVariableInt(name_multiplicity_jets_)),
btagDiscriminatorAverage_tagged_(MvaVariableFloat(name_btagDiscriminatorAverage_tagged_)),
btagDiscriminatorAverage_untagged_(MvaVariableFloat(name_averageBtagDiscriminatorUntagged_)),
minDeltaR_jet_jet_(MvaVariableFloat(name_minDeltaR_jet_jet_)),
ptSum_jets_leptons_(MvaVariableFloat(name_ptSum_jets_leptons_))
{}



MvaVariablesEventClassification::MvaVariablesEventClassification(
    const int multiplicity_jets,
    const double& btagDiscriminatorAverage_tagged, const double& btagDiscriminatorAverage_untagged,
    const double& minDeltaR_jet_jet, const double& ptSum_jets_leptons,
    const double& eventWeight):
MvaVariablesBase(eventWeight),
multiplicity_jets_(MvaVariableInt(name_multiplicity_jets_)),
btagDiscriminatorAverage_tagged_(MvaVariableFloat(name_btagDiscriminatorAverage_tagged_)),
btagDiscriminatorAverage_untagged_(MvaVariableFloat(name_averageBtagDiscriminatorUntagged_)),
minDeltaR_jet_jet_(MvaVariableFloat(name_minDeltaR_jet_jet_)),
ptSum_jets_leptons_(MvaVariableFloat(name_ptSum_jets_leptons_))
{
    // Fill the variables for MVA TTree
    multiplicity_jets_.value_ = multiplicity_jets;
    btagDiscriminatorAverage_tagged_.value_ = btagDiscriminatorAverage_tagged;
    btagDiscriminatorAverage_untagged_.value_ = btagDiscriminatorAverage_untagged;
    minDeltaR_jet_jet_.value_ = minDeltaR_jet_jet;
    ptSum_jets_leptons_.value_ = ptSum_jets_leptons;
}



MvaVariablesEventClassification* MvaVariablesEventClassification::fillVariables(const tth::RecoObjectIndices& recoObjectIndices,
                                                                                const RecoObjects& recoObjects,
                                                                                const double& eventWeight)
{
    using ROOT::Math::VectorUtil::DeltaR;
    
    // Access relevant objects and indices
    const std::vector<double>& jetBtag(*recoObjects.jetBTagCSV_);
    const VLV& leptons(*recoObjects.allLeptons_);
    const VLV& jets(*recoObjects.jets_);
    
    // Calculate several jet-dependent quantities
    double btagDiscriminatorSumTagged(0.);
    double btagDiscriminatorSumUntagged(0.);
    double minDeltaR_jet_jet(999.);
    double ptSumJets(0.);
    for(auto i_index = recoObjectIndices.jetIndices_.begin(); i_index != recoObjectIndices.jetIndices_.end(); ++i_index){
        // Calculate the btag-discriminator averages, setting values<0. to 0.
        const double& btagDiscriminator(jetBtag.at(*i_index));
        const double btagDiscriminatorPositive(btagDiscriminator>=0. ? btagDiscriminator : 0.);
        if(std::find(recoObjectIndices.bjetIndices_.begin(), recoObjectIndices.bjetIndices_.end(), *i_index) != recoObjectIndices.bjetIndices_.end()){
            btagDiscriminatorSumTagged += btagDiscriminatorPositive;
        }
        else{
            btagDiscriminatorSumUntagged += btagDiscriminatorPositive;
        }
        
        // Minimum deltaR between jets
        for(auto j_index = i_index+1; j_index != recoObjectIndices.jetIndices_.end(); ++j_index){
            const double deltaR = DeltaR(jets.at(*i_index), jets.at(*j_index));
            if(deltaR < minDeltaR_jet_jet) minDeltaR_jet_jet = deltaR;
        }
        
        // Scalar sum of pt of all jets
        ptSumJets += jets.at(*i_index).pt();
    }
    const int multiplicity_jets(recoObjectIndices.jetIndices_.size());
    const int numberOfTaggedJets(recoObjectIndices.bjetIndices_.size());
    const int numberOfUntaggedJets(multiplicity_jets - numberOfTaggedJets);
    const double btagDiscriminatorAverage_tagged = numberOfTaggedJets>0 ? btagDiscriminatorSumTagged/static_cast<double>(numberOfTaggedJets) : 0.;
    const double btagDiscriminatorAverage_untagged = numberOfUntaggedJets>0 ? btagDiscriminatorSumUntagged/static_cast<double>(numberOfUntaggedJets) : 0.;
    const double ptSum_jets_leptons = ptSumJets + leptons.at(recoObjectIndices.leptonIndex_).pt() + leptons.at(recoObjectIndices.antiLeptonIndex_).pt();
    
    return new MvaVariablesEventClassification(multiplicity_jets,
                                               btagDiscriminatorAverage_tagged, btagDiscriminatorAverage_untagged,
                                               minDeltaR_jet_jet, ptSum_jets_leptons,
                                               eventWeight);
}








