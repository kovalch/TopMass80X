#include <algorithm>

#include "MvaVariablesEventClassification.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"





MvaVariablesEventClassification::MvaVariablesEventClassification():
MvaVariablesBase(),
btagDiscriminatorAverageTagged_(MvaVariableFloat(name_btagDiscriminatorAverageTagged_)),
btagDiscriminatorAverageUntagged_(MvaVariableFloat(name_averageBtagDiscriminatorUntagged_))
{}



MvaVariablesEventClassification::MvaVariablesEventClassification(const double& btagDiscriminatorAverageTagged,
                                                                 const double& btagDiscriminatorAverageUntagged,
                                                                 const double& eventWeight):
MvaVariablesBase(eventWeight),
btagDiscriminatorAverageTagged_(MvaVariableFloat(name_btagDiscriminatorAverageTagged_)),
btagDiscriminatorAverageUntagged_(MvaVariableFloat(name_averageBtagDiscriminatorUntagged_))
{
    // Fill the variables for MVA TTree
    btagDiscriminatorAverageTagged_.value_ = btagDiscriminatorAverageTagged;
    btagDiscriminatorAverageUntagged_.value_ = btagDiscriminatorAverageUntagged;
}



MvaVariablesEventClassification MvaVariablesEventClassification::fillVariables(const tth::RecoObjectIndices& recoObjectIndices,
                                                                               const RecoObjects& recoObjects,
                                                                               const double& eventWeight)
{
    // Access relevant objects and indices
    const std::vector<double>& jetBtag(*recoObjects.jetBTagCSV_);
    
    // Calculate the btag-discriminator averages, setting values<0. to 0.
    double btagDiscriminatorSumTagged(0.);
    double btagDiscriminatorSumUntagged(0.);
    for(const int index : recoObjectIndices.jetIndices_){
        const double& btagDiscriminator(jetBtag.at(index));
        const double btagDiscriminatorPositive(btagDiscriminator>=0. ? btagDiscriminator : 0.);
        if(std::find(recoObjectIndices.bjetIndices_.begin(), recoObjectIndices.bjetIndices_.end(), index) != recoObjectIndices.bjetIndices_.end()){
            btagDiscriminatorSumTagged += btagDiscriminatorPositive;
        }
        else{
            btagDiscriminatorSumUntagged += btagDiscriminatorPositive;
        }
    }
    const int numberOfTaggedJets(recoObjectIndices.bjetIndices_.size());
    const int numberOfUntaggedJets(recoObjectIndices.jetIndices_.size() - numberOfTaggedJets);
    const double btagDiscriminatorAverageTagged = numberOfTaggedJets>0 ? btagDiscriminatorSumTagged/static_cast<double>(numberOfTaggedJets) : 0.;
    const double btagDiscriminatorAverageUntagged = numberOfUntaggedJets>0 ? btagDiscriminatorSumUntagged/static_cast<double>(numberOfUntaggedJets) : 0.;
    
    return MvaVariablesEventClassification(btagDiscriminatorAverageTagged, btagDiscriminatorAverageUntagged, eventWeight);
}








