#include <iostream>
#include <cstdlib>

#include <TTree.h>
#include <TString.h>
#include <Rtypes.h>

#include "MvaTreeHandlerEventClassification.h"
#include "MvaVariablesBase.h"
#include "MvaVariablesEventClassification.h"
#include "MvaTreePlotterBase.h"
#include "MvaTreePlotterEventClassification.h"
#include "JetCategories.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"





MvaTreeHandlerEventClassification::MvaTreeHandlerEventClassification(const char* mvaInputDir,
                                                                     const std::vector<TString>& selectionStepsNoCategories,
                                                                     const std::vector<TString>& stepsForCategories,
                                                                     const JetCategories* jetCategories):
MvaTreeHandlerBase("event_", mvaInputDir, selectionStepsNoCategories, stepsForCategories, jetCategories)
{
    std::cout<<"--- Beginning setting up MVA tree handler for event classification\n";
    std::cout<<"=== Finishing setting up MVA tree handler for event classification\n\n";
}



void MvaTreeHandlerEventClassification::fillVariables(const RecoObjects& recoObjects, const CommonGenObjects&,
                                                      const TopGenObjects&, const HiggsGenObjects&,
                                                      const KinRecoObjects&,
                                                      const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices&,
                                                      const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                                      const double& weight, const TString&,
                                                      std::vector<MvaVariablesBase*>& mvaVariables)const
{
    // Get MVA input variables
    MvaVariablesBase* mvaVariablesEventClassification = MvaVariablesEventClassification::fillVariables(recoObjectIndices, recoObjects, weight);
    
    // Fill the MVA variables
    mvaVariables.push_back(mvaVariablesEventClassification);
}



void MvaTreeHandlerEventClassification::createAndFillBranches(TTree* tree, const std::vector<MvaVariablesBase*>& v_mvaVariables)const
{
    MvaVariablesEventClassification* const mvaVariablesEventClassification = new MvaVariablesEventClassification();
    
    this->createBranch(tree, mvaVariablesEventClassification->eventWeight_);
    this->createBranch(tree, mvaVariablesEventClassification->multiplicity_jets_);
    this->createBranch(tree, mvaVariablesEventClassification->btagDiscriminatorAverage_tagged_);
    this->createBranch(tree, mvaVariablesEventClassification->btagDiscriminatorAverage_untagged_);
    this->createBranch(tree, mvaVariablesEventClassification->minDeltaR_jet_jet_);
    this->createBranch(tree, mvaVariablesEventClassification->ptSum_jets_leptons_);
    
    for(const MvaVariablesBase* mvaVariablesTmp : v_mvaVariables){
        const MvaVariablesEventClassification* mvaVariablesEventClassificationTmp = dynamic_cast<const MvaVariablesEventClassification*>(mvaVariablesTmp);
        if(!mvaVariablesEventClassificationTmp){
            std::cerr<<"ERROR in MvaTreeHandlerEventClassification::createAndFillBranches()! MvaVariables are of wrong type, cannot typecast\n"
                     <<"...break\n"<<std::endl;
            exit(395);
        }
        
        *mvaVariablesEventClassification = *mvaVariablesEventClassificationTmp;
        tree->Fill();
    }
    
    delete mvaVariablesEventClassification;
}



void MvaTreeHandlerEventClassification::importBranches(TTree* tree, std::vector<MvaVariablesBase*>& v_mvaVariables)const
{
    // Set up variables struct
    MvaVariablesEventClassification mvaVariablesEventClassification;
    
    // Set branch addresses
    this->importBranch(tree, mvaVariablesEventClassification.eventWeight_);
    this->importBranch(tree, mvaVariablesEventClassification.multiplicity_jets_);
    this->importBranch(tree, mvaVariablesEventClassification.btagDiscriminatorAverage_tagged_);
    this->importBranch(tree, mvaVariablesEventClassification.btagDiscriminatorAverage_untagged_);
    this->importBranch(tree, mvaVariablesEventClassification.minDeltaR_jet_jet_);
    this->importBranch(tree, mvaVariablesEventClassification.ptSum_jets_leptons_);
    
    // Loop over all tree entries and fill vector of structs
    for(Long64_t iEntry = 0; iEntry < tree->GetEntries(); ++iEntry){
        tree->GetEntry(iEntry);
        MvaVariablesEventClassification* const mvaVariablesEventClassificationPtr = new MvaVariablesEventClassification();
        *mvaVariablesEventClassificationPtr = mvaVariablesEventClassification;
        v_mvaVariables.push_back(mvaVariablesEventClassificationPtr);
    }
}



MvaTreePlotterBase* MvaTreeHandlerEventClassification::setPlotter(const std::map<TString, std::vector<MvaVariablesBase*> >& m_stepMvaVariables,
                                                                  const bool separationPowerPlots)const
{
    return new MvaTreePlotterEventClassification(m_stepMvaVariables, separationPowerPlots);
}








