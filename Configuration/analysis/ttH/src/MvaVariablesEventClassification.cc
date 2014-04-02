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
ptSum_jets_leptons_(MvaVariableFloat(name_ptSum_jets_leptons_)),
multiplicity_higgsLikeDijet15_(MvaVariableInt(name_multiplicity_higgsLikeDijet15_)),
mass_higgsLikeDijet_(MvaVariableFloat(name_mass_higgsLikeDijet_)),
mass_higgsLikeDijet2_(MvaVariableFloat(name_mass_higgsLikeDijet2_))
{}



MvaVariablesEventClassification::MvaVariablesEventClassification(
    const int multiplicity_jets,
    const double& btagDiscriminatorAverage_tagged, const double& btagDiscriminatorAverage_untagged,
    const double& minDeltaR_jet_jet, const double& ptSum_jets_leptons,
    const int multiplicity_higgsLikeDijet15,
    const double& mass_higgsLikeDijet, const double& mass_higgsLikeDijet2,
    const double& eventWeight):
MvaVariablesBase(eventWeight),
multiplicity_jets_(MvaVariableInt(name_multiplicity_jets_)),
btagDiscriminatorAverage_tagged_(MvaVariableFloat(name_btagDiscriminatorAverage_tagged_)),
btagDiscriminatorAverage_untagged_(MvaVariableFloat(name_averageBtagDiscriminatorUntagged_)),
minDeltaR_jet_jet_(MvaVariableFloat(name_minDeltaR_jet_jet_)),
ptSum_jets_leptons_(MvaVariableFloat(name_ptSum_jets_leptons_)),
multiplicity_higgsLikeDijet15_(MvaVariableInt(name_multiplicity_higgsLikeDijet15_)),
mass_higgsLikeDijet_(MvaVariableFloat(name_mass_higgsLikeDijet_)),
mass_higgsLikeDijet2_(MvaVariableFloat(name_mass_higgsLikeDijet2_))
{
    // Fill the variables for MVA TTree
    multiplicity_jets_.value_ = multiplicity_jets;
    btagDiscriminatorAverage_tagged_.value_ = btagDiscriminatorAverage_tagged;
    btagDiscriminatorAverage_untagged_.value_ = btagDiscriminatorAverage_untagged;
    minDeltaR_jet_jet_.value_ = minDeltaR_jet_jet;
    ptSum_jets_leptons_.value_ = ptSum_jets_leptons;
    multiplicity_higgsLikeDijet15_.value_ = multiplicity_higgsLikeDijet15;
    mass_higgsLikeDijet_.value_ = mass_higgsLikeDijet;
    mass_higgsLikeDijet2_.value_ = mass_higgsLikeDijet2;
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
    double minDeltaRJetJet(999.);
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
            if(deltaR < minDeltaRJetJet) minDeltaRJetJet = deltaR;
        }
        
        // Scalar sum of pt of all jets
        ptSumJets += jets.at(*i_index).pt();
    }
    const int numberOfJets(recoObjectIndices.jetIndices_.size());
    const int numberOfTaggedJets(recoObjectIndices.bjetIndices_.size());
    const int numberOfUntaggedJets(numberOfJets - numberOfTaggedJets);
    const double btagDiscriminatorAverage_tagged = numberOfTaggedJets>0 ? btagDiscriminatorSumTagged/static_cast<double>(numberOfTaggedJets) : 0.;
    const double btagDiscriminatorAverage_untagged = numberOfUntaggedJets>0 ? btagDiscriminatorSumUntagged/static_cast<double>(numberOfUntaggedJets) : 0.;
    const double ptSumJetsLeptons = ptSumJets + leptons.at(recoObjectIndices.leptonIndex_).pt() + leptons.at(recoObjectIndices.antiLeptonIndex_).pt();
    
    // Calculate several dijet dependent quantities
    int numberOfHiggsLikeDijet15(0);
    double higgsLikeDijetMass(-999.);
    double higgsLikeDijetMass2(-999.);
    for(const auto& indexPair : recoObjectIndices.jetIndexPairs_){
        const bool hasBtag = (std::find(recoObjectIndices.bjetIndices_.begin(), recoObjectIndices.bjetIndices_.end(), indexPair.first) != recoObjectIndices.bjetIndices_.end()) || 
                             (std::find(recoObjectIndices.bjetIndices_.begin(), recoObjectIndices.bjetIndices_.end(), indexPair.second) != recoObjectIndices.bjetIndices_.end());
        constexpr double higgsMass(125.);
        const double dijetMass = (jets.at(indexPair.first) + jets.at(indexPair.second)).M();
        if(std::abs(higgsMass - dijetMass) < std::abs(higgsMass - higgsLikeDijetMass)) higgsLikeDijetMass = dijetMass;
        if(hasBtag){
            if(std::abs(higgsMass - dijetMass) < std::abs(higgsMass - higgsLikeDijetMass2)) higgsLikeDijetMass2 = dijetMass;
            if(std::abs(dijetMass - higgsMass) < 15.) ++numberOfHiggsLikeDijet15;
        }
    }
    
    
    return new MvaVariablesEventClassification(numberOfJets,
                                               btagDiscriminatorAverage_tagged, btagDiscriminatorAverage_untagged,
                                               minDeltaRJetJet, ptSumJetsLeptons,
                                               numberOfHiggsLikeDijet15,
                                               higgsLikeDijetMass, higgsLikeDijetMass2,
                                               eventWeight);
}








