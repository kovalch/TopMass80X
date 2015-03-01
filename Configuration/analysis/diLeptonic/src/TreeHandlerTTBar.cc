#include <iostream>
#include <cstdlib>

#include <TTree.h>
#include <TString.h>
#include <Rtypes.h>

#include "TreeHandlerTTBar.h"
#include "VariablesBase.h"
#include "VariablesTTBar.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/KinematicReconstructionSolution.h"





TreeHandlerTTBar::TreeHandlerTTBar(const char* inputDir,
                                             const std::vector<TString>& selectionStepsNoCategories):
TreeHandlerBase("ttBar_", inputDir, selectionStepsNoCategories, new VariablesTTBar())
{
    std::cout<<"--- Beginning setting up MVA tree handler for top jets assignment\n";
    std::cout<<"=== Finishing setting up MVA tree handler for top jets assignment\n\n";
}



void TreeHandlerTTBar::fillVariables(const EventMetadata& eventMetadata,
                                     const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                          const TopGenObjects& topGenObjects,
                                          const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                                          const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                                          const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                                          const double& weight, const TString&,
                                          std::vector<VariablesBase*>& variables)const
{
    // Loop over all jet combinations and get MVA input variables
    const std::vector<VariablesBase*> v_variablesTTbar = 
            VariablesTTBar::fillVariables(eventMetadata, recoObjects, commonGenObjects,topGenObjects,kinematicReconstructionSolutions, recoObjectIndices,  genObjectIndices,genLevelWeights,recoLevelWeights,weight);
    
    // Fill the MVA variables
    variables.insert(variables.end(), v_variablesTTbar.begin(), v_variablesTTbar.end());
}



void TreeHandlerTTBar::bookBranches(TTree* tree, VariablesBase* const variables_)const
{
    VariablesTTBar* const variablesTTBar = dynamic_cast<VariablesTTBar*>(variables_);
    
    this->createBranch(tree, variablesTTBar->eventWeight_);
    
    this->createBranch(tree, variablesTTBar->isTopGen_);
    this->createBranch(tree, variablesTTBar->entry_);
    this->createBranch(tree, variablesTTBar->isKinReco_);
    this->createBranch(tree, variablesTTBar->trueLevelWeight_);
    
    this->createBranch(tree, variablesTTBar->top_pt_);
    this->createBranch(tree, variablesTTBar->ttbar_delta_phi_);
    this->createBranch(tree, variablesTTBar->ttbar_pt_);
    this->createBranch(tree, variablesTTBar->top_rapidity_);
    this->createBranch(tree, variablesTTBar->ttbar_delta_eta_);
    this->createBranch(tree, variablesTTBar->ttbar_rapidity_);
    this->createBranch(tree, variablesTTBar->ttbar_mass_);
    this->createBranch(tree, variablesTTBar->jet_multiplicity_);
    this->createBranch(tree, variablesTTBar->x1_);
    this->createBranch(tree, variablesTTBar->x2_);
    
    this->createBranch(tree, variablesTTBar->gen_top_pt_);
    this->createBranch(tree, variablesTTBar->gen_ttbar_delta_phi_);
    this->createBranch(tree, variablesTTBar->gen_ttbar_pt_);
    this->createBranch(tree, variablesTTBar->gen_top_rapidity_);
    this->createBranch(tree, variablesTTBar->gen_ttbar_delta_eta_);
    this->createBranch(tree, variablesTTBar->gen_ttbar_rapidity_);
    this->createBranch(tree, variablesTTBar->gen_ttbar_mass_);
    this->createBranch(tree, variablesTTBar->gen_jet_multiplicity_);
    this->createBranch(tree, variablesTTBar->gen_x1_);
    this->createBranch(tree, variablesTTBar->gen_x2_);

}



void TreeHandlerTTBar::fillBranches(TTree* tree, const std::vector<VariablesBase*>& v_variables)
{
    
    for(const VariablesBase* variablesTmp : v_variables){
        const VariablesTTBar* variablesTTBarTmp = dynamic_cast<const VariablesTTBar*>(variablesTmp);
        if(!variablesTTBarTmp){
            std::cerr<<"ERROR in TreeHandlerTTBar::fillBranches()! variables are of wrong type, cannot typecast\n"
                     <<"...break\n"<<std::endl;
            exit(395);
        }
        
        *(dynamic_cast<VariablesTTBar*>(variables_)) = *variablesTTBarTmp;
        
        tree->Fill();
    }
}



void TreeHandlerTTBar::importBranches(TTree* tree, std::vector<VariablesBase*>& v_variables)const
{
    // Set up variables struct
    VariablesTTBar variablesTTBar;
    
    // Set branch addresses
    this->importBranch(tree, variablesTTBar.eventWeight_);
    
    this->importBranch(tree, variablesTTBar.isTopGen_);
    this->importBranch(tree, variablesTTBar.entry_);
    this->importBranch(tree, variablesTTBar.isKinReco_);
    this->importBranch(tree, variablesTTBar.trueLevelWeight_);
    
    this->importBranch(tree, variablesTTBar.top_pt_);
    this->importBranch(tree, variablesTTBar.ttbar_delta_phi_);
    this->importBranch(tree, variablesTTBar.ttbar_pt_);
    this->importBranch(tree, variablesTTBar.top_rapidity_);
    this->importBranch(tree, variablesTTBar.ttbar_delta_eta_);
    this->importBranch(tree, variablesTTBar.ttbar_rapidity_);
    this->importBranch(tree, variablesTTBar.ttbar_mass_);
    this->importBranch(tree, variablesTTBar.jet_multiplicity_);
    this->importBranch(tree, variablesTTBar.x1_);
    this->importBranch(tree, variablesTTBar.x2_);
    
    this->importBranch(tree, variablesTTBar.gen_top_pt_);
    this->importBranch(tree, variablesTTBar.gen_ttbar_delta_phi_);
    this->importBranch(tree, variablesTTBar.gen_ttbar_pt_);
    this->importBranch(tree, variablesTTBar.gen_top_rapidity_);
    this->importBranch(tree, variablesTTBar.gen_ttbar_delta_eta_);
    this->importBranch(tree, variablesTTBar.gen_ttbar_rapidity_);
    this->importBranch(tree, variablesTTBar.gen_ttbar_mass_);
    this->importBranch(tree, variablesTTBar.gen_jet_multiplicity_);
    this->importBranch(tree, variablesTTBar.gen_x1_);
    this->importBranch(tree, variablesTTBar.gen_x2_);
    
    // Loop over all tree entries and fill vector of structs
    for(Long64_t iEntry = 0; iEntry < tree->GetEntries(); ++iEntry){
        tree->GetEntry(iEntry);
        VariablesTTBar* const variablesTTBarPtr = new VariablesTTBar();
        *variablesTTBarPtr = variablesTTBar;
        v_variables.push_back(variablesTTBarPtr);
    }
}








