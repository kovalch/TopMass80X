#include <iostream>
#include <cstdlib>

#include <TTree.h>
#include <TString.h>
#include <Rtypes.h>

#include "MvaTreeHandlerTopJets.h"
#include "MvaVariablesBase.h"
#include "MvaVariablesTopJets.h"
#include "MvaTreePlotterBase.h"
#include "MvaTreePlotterTopJets.h"
#include "JetCategories.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/KinematicReconstructionSolution.h"





MvaTreeHandlerTopJets::MvaTreeHandlerTopJets(const char* mvaInputDir,
                                             const std::vector<TString>& selectionStepsNoCategories,
                                             const std::vector<TString>& stepsForCategories,
                                             const JetCategories* jetCategories):
MvaTreeHandlerBase("topJets_", mvaInputDir, selectionStepsNoCategories, stepsForCategories, jetCategories)
{
    std::cout<<"--- Beginning setting up MVA tree handler for top jets assignment\n";
    std::cout<<"=== Finishing setting up MVA tree handler for top jets assignment\n\n";
}



void MvaTreeHandlerTopJets::fillVariables(const EventMetadata&,
                                          const RecoObjects& recoObjects, const CommonGenObjects&,
                                          const TopGenObjects&, const HiggsGenObjects&,
                                          const KinematicReconstructionSolutions&,
                                          const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                                          const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                          const double& weight, const TString&,
                                          std::vector<MvaVariablesBase*>& mvaVariables)const
{
    // Loop over all jet combinations and get MVA input variables
    const std::vector<MvaVariablesBase*> v_mvaVariablesTopJets = 
            MvaVariablesTopJets::fillVariables(recoObjectIndices, genObjectIndices, recoObjects, weight);
    
    // Fill the MVA variables
    mvaVariables.insert(mvaVariables.end(), v_mvaVariablesTopJets.begin(), v_mvaVariablesTopJets.end());
}



void MvaTreeHandlerTopJets::createAndFillBranches(TTree* tree, const std::vector<MvaVariablesBase*>& v_mvaVariables)const
{
    MvaVariablesTopJets* const mvaVariablesTopJets = new MvaVariablesTopJets();
    
    this->createBranch(tree, mvaVariablesTopJets->lastInEvent_);
    this->createBranch(tree, mvaVariablesTopJets->eventWeight_);
    this->createBranch(tree, mvaVariablesTopJets->bQuarkRecoJetMatched_);
    this->createBranch(tree, mvaVariablesTopJets->correctCombination_);
    this->createBranch(tree, mvaVariablesTopJets->swappedCombination_);
    this->createBranch(tree, mvaVariablesTopJets->jetChargeDiff_);
    this->createBranch(tree, mvaVariablesTopJets->meanDeltaPhi_b_met_);
    this->createBranch(tree, mvaVariablesTopJets->massDiff_recoil_bbbar_);
    this->createBranch(tree, mvaVariablesTopJets->pt_b_antiLepton_);
    this->createBranch(tree, mvaVariablesTopJets->pt_antiB_lepton_);
    this->createBranch(tree, mvaVariablesTopJets->deltaR_b_antiLepton_);
    this->createBranch(tree, mvaVariablesTopJets->deltaR_antiB_lepton_);
    this->createBranch(tree, mvaVariablesTopJets->btagDiscriminatorSum_);
    this->createBranch(tree, mvaVariablesTopJets->deltaPhi_antiBLepton_bAntiLepton_);
    this->createBranch(tree, mvaVariablesTopJets->massDiff_fullBLepton_bbbar_);
    this->createBranch(tree, mvaVariablesTopJets->meanMt_b_met_);
    this->createBranch(tree, mvaVariablesTopJets->massSum_antiBLepton_bAntiLepton_);
    this->createBranch(tree, mvaVariablesTopJets->massDiff_antiBLepton_bAntiLepton_);
    
    for(const MvaVariablesBase* mvaVariablesTmp : v_mvaVariables){
        const MvaVariablesTopJets* mvaVariablesTopJetsTmp = dynamic_cast<const MvaVariablesTopJets*>(mvaVariablesTmp);
        if(!mvaVariablesTopJetsTmp){
            std::cerr<<"ERROR in MvaTreeHandlerTopJets::createAndFillBranches()! MvaVariables are of wrong type, cannot typecast\n"
                     <<"...break\n"<<std::endl;
            exit(395);
        }
        
        *mvaVariablesTopJets = *mvaVariablesTopJetsTmp;
        tree->Fill();
    }
    
    delete mvaVariablesTopJets;
}



void MvaTreeHandlerTopJets::importBranches(TTree* tree, std::vector<MvaVariablesBase*>& v_mvaVariables)const
{
    // Set up variables struct
    MvaVariablesTopJets mvaVariablesTopJets;
    
    // Set branch addresses
    this->importBranch(tree, mvaVariablesTopJets.lastInEvent_);
    this->importBranch(tree, mvaVariablesTopJets.eventWeight_);
    this->importBranch(tree, mvaVariablesTopJets.bQuarkRecoJetMatched_);
    this->importBranch(tree, mvaVariablesTopJets.correctCombination_);
    this->importBranch(tree, mvaVariablesTopJets.swappedCombination_);
    this->importBranch(tree, mvaVariablesTopJets.jetChargeDiff_);
    this->importBranch(tree, mvaVariablesTopJets.meanDeltaPhi_b_met_);
    this->importBranch(tree, mvaVariablesTopJets.massDiff_recoil_bbbar_);
    this->importBranch(tree, mvaVariablesTopJets.pt_b_antiLepton_);
    this->importBranch(tree, mvaVariablesTopJets.pt_antiB_lepton_);
    this->importBranch(tree, mvaVariablesTopJets.deltaR_b_antiLepton_);
    this->importBranch(tree, mvaVariablesTopJets.deltaR_antiB_lepton_);
    this->importBranch(tree, mvaVariablesTopJets.btagDiscriminatorSum_);
    this->importBranch(tree, mvaVariablesTopJets.deltaPhi_antiBLepton_bAntiLepton_);
    this->importBranch(tree, mvaVariablesTopJets.massDiff_fullBLepton_bbbar_);
    this->importBranch(tree, mvaVariablesTopJets.meanMt_b_met_);
    this->importBranch(tree, mvaVariablesTopJets.massSum_antiBLepton_bAntiLepton_);
    this->importBranch(tree, mvaVariablesTopJets.massDiff_antiBLepton_bAntiLepton_);
    
    // Loop over all tree entries and fill vector of structs
    for(Long64_t iEntry = 0; iEntry < tree->GetEntries(); ++iEntry){
        tree->GetEntry(iEntry);
        MvaVariablesTopJets* const mvaVariablesTopJetsPtr = new MvaVariablesTopJets();
        *mvaVariablesTopJetsPtr = mvaVariablesTopJets;
        v_mvaVariables.push_back(mvaVariablesTopJetsPtr);
    }
}



MvaTreePlotterBase* MvaTreeHandlerTopJets::setPlotter(const std::map<TString, std::vector<MvaVariablesBase*> >& m_stepMvaVariables,
                                                      const bool separationPowerPlots)const
{
    return new MvaTreePlotterTopJets(m_stepMvaVariables, separationPowerPlots);
}








