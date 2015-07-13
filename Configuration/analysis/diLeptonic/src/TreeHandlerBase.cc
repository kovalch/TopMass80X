#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <utility>

#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TSelectorList.h>
#include <TIterator.h>
#include <TObject.h>
#include <Rtypes.h>

#include "TreeHandlerBase.h"
#include "VariablesBase.h"
#include "analysisStructs.h"
#include "ttbarUtils.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/KinematicReconstructionSolution.h"





TreeHandlerBase::TreeHandlerBase(const TString& prefix,
                                       const char* inputDir,
                                       const std::vector<TString>& selectionStepsNoCategories,
                                       VariablesBase* const variables):
variables_(variables),
prefix_(prefix),
selectorList_(0),
selectionSteps_(selectionStepsNoCategories),
inputDir_(inputDir)
{}



void TreeHandlerBase::book()
{
    TSelectorList* output = new TSelectorList();
    selectorList_ = output;
    for(const auto& stepShort : selectionSteps_){
        const TString step = ttbar::stepName(stepShort);
        this->addStep(step);
    }
}



void TreeHandlerBase::addStep(const TString& step)
{
    // Check whether step already exists
    if(this->checkExistence(step)){
        std::cout<<"Warning in addStep()! Selection step already contained: "<<step
                 <<"\n...skip this one\n";
        return;
    }

    // Book the step
    m_stepVariables_[step] = std::vector<VariablesBase*>();
    m_stepTree_[step] = new TTree(prefix_+"treeVariables"+step, prefix_+"treeVariables");
    this->store(m_stepTree_[step]);
    this->bookBranches(m_stepTree_[step],variables_);
}



bool TreeHandlerBase::checkExistence(const TString& step)const
{
    return m_stepVariables_.find(step) != m_stepVariables_.end();
}



void TreeHandlerBase::fill(const EventMetadata& eventMetadata,
                           const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                              const TopGenObjects& topGenObjects,
                              const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                              const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                              const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                              const double& weight, const TString& stepShort)
{
    // Set up step name and check if step exists
    const TString step = ttbar::stepName(stepShort);
    const bool stepExists(this->checkExistence(step));
    if(!stepExists) return;
    
    // Fill the variables of the specific treeHandler
    std::vector<VariablesBase*> variables;
    this->fillVariables(eventMetadata,
                        recoObjects, commonGenObjects,
                        topGenObjects, kinematicReconstructionSolutions,
                        recoObjectIndices, genObjectIndices,
                        genLevelWeights, recoLevelWeights,
                        weight, step,
                        variables);
    this->fillBranches(m_stepTree_.at(step),variables);
    VariablesBase::clearVariables(variables);
    variables.clear();
}



void TreeHandlerBase::writeTrees(const TString& outputFilename,
                                    const Channel::Channel& channel, const Systematic::Systematic& systematic)
{
    // Create output file for tree
    TString f_savename = common::assignFolder(inputDir_, channel, systematic);
    f_savename.Append(outputFilename);
    TFile outputFile(f_savename, "RECREATE");
    std::cout<<"\nOutput file for input trees: "<<f_savename<<"\n";
    
    // Write file and cleanup
    TIterator* it = selectorList_->MakeIterator();
    while(TObject* obj = it->Next()){
        obj->Write();
    }
    
    std::cout<<"tree variables multiplicity per step (step, no. of entries):\n";
    for(auto stepTree : m_stepTree_){
        std::cout<<"\t"<<stepTree.first<<" , "<< stepTree.second->GetEntries() <<"\n";
    }
    
    outputFile.Close();
    selectorList_->SetOwner();
    selectorList_->Clear();
}



void TreeHandlerBase::writeTrees(TSelectorList* output)
{
    std::cout<<"--- Beginning production of input trees\n";
    
    // Set pointer to output, so that TTrees are owned by it
    selectorList_ = output;
    
    std::cout<<"tree variables multiplicity per step (step, no. of entries):\n";
    for(auto stepTree : m_stepTree_){
        std::cout<<"\t"<<stepTree.first<<" , "<< (this->store(stepTree.second))->GetEntries() <<"\n";
    }
    
    std::cout<<"=== Finishing production of input trees\n\n";
}



void TreeHandlerBase::bookBranches(TTree*, VariablesBase* const )const
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method bookBranches() in TreeHandlerBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



void TreeHandlerBase::fillBranches(TTree*, const std::vector<VariablesBase*>&)
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method fillBranches() in TreeHandlerBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



void TreeHandlerBase::createBranch(TTree* tree, const VariableInt& variable)const
{
    std::string name(variable.name());
    std::string nameType(name);
    nameType.append("/").append(variable.type());
    tree->Branch(name.data(), (Long_t)&variable.value_, nameType.data());
}



void TreeHandlerBase::createBranch(TTree* tree, const VariableFloat& variable)const
{
    std::string name(variable.name());
    std::string nameType(name);
    nameType.append("/").append(variable.type());
    tree->Branch(name.c_str(), (Long_t)&variable.value_, nameType.data());
}



void TreeHandlerBase::importTrees(const TString& f_savename, const TString& prefix)
{
    std::cout<<"--- Beginning import of TTrees with variables\n";
    
    // Open input file
    TFile* inputFile = TFile::Open(f_savename);
    if(inputFile->IsZombie()){
        std::cerr<<"ERROR in importTrees()! Cannot open input file to import TTrees, filename is: "
                 <<f_savename<<"\n...break\n"<<std::endl;
        exit(77);
    }
    
    // Find all trees of all steps/categories containing MVA input variables
    const TString& treeName = prefix;
    const std::vector<std::pair<TString, TString> > v_nameStepPair =
        ttbar::nameStepPairs(f_savename, treeName);
    
    // Loop over steps and import trees
    this->clear();
    for(const auto& nameStepPair : v_nameStepPair){
        TTree* tree(0);
        tree = dynamic_cast<TTree*>(inputFile->Get(nameStepPair.first));
        if(!tree){
            std::cerr<<"ERROR in importTrees()! TTree not found in file, tree name is: "
                     <<treeName<<"\n...break\n"<<std::endl;
            exit(78);
        }
        this->importBranches(tree, m_stepVariables_[nameStepPair.second]);
        tree = 0;
        
        std::cout<<"Step, Number of entries: "<<nameStepPair.second<<" , "<<m_stepVariables_[nameStepPair.second].size()<<"\n";
    }
    inputFile->Close();
    
    std::cout<<"=== Finishing import of TTrees with variables\n\n";
}



void TreeHandlerBase::importBranches(TTree*, std::vector<VariablesBase*>&)const
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method importBranches() in MvaTreeHandlerBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



void TreeHandlerBase::importBranch(TTree* tree, VariableInt& variable)const
{
    tree->SetBranchAddress(variable.name().data(), &variable.value_, &variable.branch_);
}



void TreeHandlerBase::importBranch(TTree* tree, VariableFloat& variable)const
{
    tree->SetBranchAddress(variable.name().data(), &variable.value_, &variable.branch_);
}



const std::map<TString, std::vector<VariablesBase*> >& TreeHandlerBase::stepVariablesMap()const
{
    return m_stepVariables_;
}



void TreeHandlerBase::clear()
{
    selectorList_ = 0;
    
    for(auto& stepVariables : m_stepVariables_){
        VariablesBase::clearVariables(stepVariables.second);
    }
    m_stepVariables_.clear();
    
    for(auto& stepTree : m_stepTree_){
        if(0/*somesing exists in stepTree.second*/)stepTree.second->Delete();//FIXME
    }
    m_stepTree_.clear();
}



void TreeHandlerBase::fillVariables(const EventMetadata& ,
                                    const RecoObjects&, const CommonGenObjects&,
                                       const TopGenObjects&,
                                       const KinematicReconstructionSolutions&,
                                       const ttbar::RecoObjectIndices&, const ttbar::GenObjectIndices&,
                                       const ttbar::GenLevelWeights&, const ttbar::RecoLevelWeights&,
                                       const double&, const TString&,
                                       std::vector<VariablesBase*>&)const
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method fillVariables() in MvaTreeHandlerBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}











