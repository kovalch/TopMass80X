#include <iostream>
#include <cstdlib>

#include <TTree.h>
#include <TString.h>
#include <Rtypes.h>

#include "TreeHandlerBoostedTop.h"
#include "VariablesBase.h"
#include "VariablesBoostedTop.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/KinematicReconstructionSolution.h"





TreeHandlerBoostedTop::TreeHandlerBoostedTop(const char* inputDir,
                                             const std::vector<TString>& selectionStepsNoCategories):
TreeHandlerBase("bTop_", inputDir, selectionStepsNoCategories, new VariablesBoostedTop())
{
    std::cout<<"--- Beginning setting up MVA tree handler for top jets assignment\n";
    std::cout<<"=== Finishing setting up MVA tree handler for top jets assignment\n\n";
}



void TreeHandlerBoostedTop::fillVariables(const EventMetadata& eventMetadata,
                                     const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                          const TopGenObjects& topGenObjects,
                                          const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                                          const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                                          const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                                          const double& weight, const TString&,
                                          std::vector<VariablesBase*>& variables)const
{
    // Loop over all jet combinations and get MVA input variables
    const std::vector<VariablesBase*> v_variablesBoostedTop = 
            VariablesBoostedTop::fillVariables(eventMetadata, recoObjects, commonGenObjects,topGenObjects,kinematicReconstructionSolutions, recoObjectIndices,  genObjectIndices,genLevelWeights,recoLevelWeights,weight);
    
    // Fill the MVA variables
    variables.insert(variables.end(), v_variablesBoostedTop.begin(), v_variablesBoostedTop.end());
}



void TreeHandlerBoostedTop::bookBranches(TTree* tree, VariablesBase* const variables_)const
{
    
    VariablesBoostedTop* const variablesBoostedTop = dynamic_cast<VariablesBoostedTop*>(variables_);

    
    this->createBranch(tree, variablesBoostedTop->eventWeight_);
    
    this->createBranch(tree, variablesBoostedTop->isTopGen_);
    this->createBranch(tree, variablesBoostedTop->entry_);
    this->createBranch(tree, variablesBoostedTop->isKinReco_);
    this->createBranch(tree, variablesBoostedTop->trueLevelWeight_);
    
    this->createBranch(tree, variablesBoostedTop->lep_pt_);
    this->createBranch(tree, variablesBoostedTop->anti_lep_pt_);
    
    this->createBranch(tree, variablesBoostedTop->top_pt_);
    this->createBranch(tree, variablesBoostedTop->topbar_pt_);
    this->createBranch(tree, variablesBoostedTop->ttbar_delta_phi_);
    this->createBranch(tree, variablesBoostedTop->ttbar_pt_);
    this->createBranch(tree, variablesBoostedTop->top_rapidity_);
    this->createBranch(tree, variablesBoostedTop->ttbar_delta_eta_);
    this->createBranch(tree, variablesBoostedTop->ttbar_rapidity_);
    this->createBranch(tree, variablesBoostedTop->ttbar_mass_);
    this->createBranch(tree, variablesBoostedTop->jet_multiplicity_);
    this->createBranch(tree, variablesBoostedTop->x1_);
    this->createBranch(tree, variablesBoostedTop->x2_);
    this->createBranch(tree, variablesBoostedTop->mlblbmet_);
    
    this->createBranch(tree, variablesBoostedTop->gen_top_pt_);
    this->createBranch(tree, variablesBoostedTop->gen_topbar_pt_);
    this->createBranch(tree, variablesBoostedTop->gen_ttbar_delta_phi_);
    this->createBranch(tree, variablesBoostedTop->gen_ttbar_pt_);
    this->createBranch(tree, variablesBoostedTop->gen_top_rapidity_);
    this->createBranch(tree, variablesBoostedTop->gen_ttbar_delta_eta_);
    this->createBranch(tree, variablesBoostedTop->gen_ttbar_rapidity_);
    this->createBranch(tree, variablesBoostedTop->gen_ttbar_mass_);
    this->createBranch(tree, variablesBoostedTop->gen_jet_multiplicity_);
    this->createBranch(tree, variablesBoostedTop->gen_x1_);
    this->createBranch(tree, variablesBoostedTop->gen_x2_);
}



void TreeHandlerBoostedTop::fillBranches(TTree* tree, const std::vector<VariablesBase*>& v_variables)
{

    
    for(const VariablesBase* variablesTmp : v_variables){
        const VariablesBoostedTop* variablesBoostedTopTmp = dynamic_cast<const VariablesBoostedTop*>(variablesTmp);
        if(!variablesBoostedTopTmp){
            std::cerr<<"ERROR in TreeHandlerBoostedTop::fillBranches()! variables are of wrong type, cannot typecast\n"
                     <<"...break\n"<<std::endl;
            exit(395);
        }
        
        *(dynamic_cast<VariablesBoostedTop*>(variables_)) = *variablesBoostedTopTmp;
        
        tree->Fill();
    }
    
}



void TreeHandlerBoostedTop::importBranches(TTree* tree, std::vector<VariablesBase*>& v_variables)const
{
    // Set up variables struct
    VariablesBoostedTop variablesBoostedTop;
    
    // Set branch addresses
    this->importBranch(tree, variablesBoostedTop.eventWeight_);
    
    this->importBranch(tree, variablesBoostedTop.isTopGen_);
    this->importBranch(tree, variablesBoostedTop.entry_);
    this->importBranch(tree, variablesBoostedTop.isKinReco_);
    this->importBranch(tree, variablesBoostedTop.trueLevelWeight_);
    
    this->importBranch(tree, variablesBoostedTop.lep_pt_);
    this->importBranch(tree, variablesBoostedTop.anti_lep_pt_);
    
    this->importBranch(tree, variablesBoostedTop.top_pt_);
    this->importBranch(tree, variablesBoostedTop.topbar_pt_);
    this->importBranch(tree, variablesBoostedTop.ttbar_delta_phi_);
    this->importBranch(tree, variablesBoostedTop.ttbar_pt_);
    this->importBranch(tree, variablesBoostedTop.top_rapidity_);
    this->importBranch(tree, variablesBoostedTop.ttbar_delta_eta_);
    this->importBranch(tree, variablesBoostedTop.ttbar_rapidity_);
    this->importBranch(tree, variablesBoostedTop.ttbar_mass_);
    this->importBranch(tree, variablesBoostedTop.jet_multiplicity_);
    this->importBranch(tree, variablesBoostedTop.x1_);
    this->importBranch(tree, variablesBoostedTop.x2_);
    this->createBranch(tree, variablesBoostedTop.mlblbmet_);
    
    this->importBranch(tree, variablesBoostedTop.gen_top_pt_);
    this->importBranch(tree, variablesBoostedTop.gen_topbar_pt_);
    this->importBranch(tree, variablesBoostedTop.gen_ttbar_delta_phi_);
    this->importBranch(tree, variablesBoostedTop.gen_ttbar_pt_);
    this->importBranch(tree, variablesBoostedTop.gen_top_rapidity_);
    this->importBranch(tree, variablesBoostedTop.gen_ttbar_delta_eta_);
    this->importBranch(tree, variablesBoostedTop.gen_ttbar_rapidity_);
    this->importBranch(tree, variablesBoostedTop.gen_ttbar_mass_);
    this->importBranch(tree, variablesBoostedTop.gen_jet_multiplicity_);
    this->importBranch(tree, variablesBoostedTop.gen_x1_);
    this->importBranch(tree, variablesBoostedTop.gen_x2_);
    
    // Loop over all tree entries and fill vector of structs
    for(Long64_t iEntry = 0; iEntry < tree->GetEntries(); ++iEntry){
        tree->GetEntry(iEntry);
        VariablesBoostedTop* const variablesBoostedTopPtr = new VariablesBoostedTop();
        *variablesBoostedTopPtr = variablesBoostedTop;
        v_variables.push_back(variablesBoostedTopPtr);
    }
}








