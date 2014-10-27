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

#include "MvaTreeHandlerBase.h"
#include "MvaVariablesBase.h"
#include "MvaTreePlotterBase.h"
#include "JetCategories.h"
#include "analysisStructs.h"
#include "higgsUtils.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/KinematicReconstructionSolution.h"





MvaTreeHandlerBase::MvaTreeHandlerBase(const TString& prefix,
                                       const char* mvaInputDir,
                                       const std::vector<TString>& selectionStepsNoCategories,
                                       const std::vector<TString>& stepsForCategories,
                                       const JetCategories* jetCategories):
prefix_(prefix),
selectorList_(0),
selectionSteps_(selectionStepsNoCategories),
stepsForCategories_(stepsForCategories),
jetCategories_(jetCategories),
mvaInputDir_(mvaInputDir)
{
    if(!jetCategories_ && stepsForCategories_.size()>0){
        std::cerr<<"ERROR in constructor for MvaTreeHandlerBase! "
                 <<"No jet categories passed, but request for category-wise selection steps\n...break\n"<<std::endl;
        exit(236);
    }
    if(jetCategories_ && stepsForCategories_.size()==0){
        std::cerr<<"ERROR in constructor for MvaTreeHandlerBase! "
                 <<"Jet categories passed, but no category-wise selection steps defined\n...break\n"<<std::endl;
        exit(237);
    }
}



void MvaTreeHandlerBase::book()
{
    for(const auto& stepShort : selectionSteps_){
        const TString step = tth::stepName(stepShort);
        this->addStep(step);
    }
    
    for(const auto& stepShort : stepsForCategories_){
        for(int category = 0; category<jetCategories_->numberOfCategories(); ++category){
            const TString step = tth::stepName(stepShort, category);
            this->addStep(step);
        }
    }
}



void MvaTreeHandlerBase::addStep(const TString& step)
{
    // Check whether step already exists
    if(this->checkExistence(step)){
        std::cout<<"Warning in addStep()! Selection step already contained: "<<step
                 <<"\n...skip this one\n";
        return;
    }

    // Book the step
    m_stepMvaVariables_[step] = std::vector<MvaVariablesBase*>();
}



bool MvaTreeHandlerBase::checkExistence(const TString& step)const
{
    return m_stepMvaVariables_.find(step) != m_stepMvaVariables_.end();
}



void MvaTreeHandlerBase::fill(const EventMetadata& eventMetadata,
                              const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                              const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                              const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                              const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                              const tth::GenLevelWeights& genLevelWeights, const tth::RecoLevelWeights& recoLevelWeights,
                              const double& weight, const TString& stepShort)
{
    // Number of selected jets and bjets
    const int numberOfJets = recoObjectIndices.jetIndices_.size();
    const int numberOfBjets = recoObjectIndices.bjetIndices_.size();
    
    // Set up step name and check if step exists
    const bool stepInCategory = stepShort.Contains("_cate");
    const TString step = stepInCategory ? stepShort : tth::stepName(stepShort);
    const bool stepExists(this->checkExistence(step));
    if(!stepInCategory && jetCategories_){
        const int categoryId = jetCategories_->categoryId(numberOfJets, numberOfBjets);
        
        // Here check the individual jet categories
        const TString fullStepName = tth::stepName(stepShort, categoryId);
        this->fill(eventMetadata,
                   recoObjects, commonGenObjects,
                   topGenObjects, higgsGenObjects,
                   kinematicReconstructionSolutions,
                   recoObjectIndices, genObjectIndices,
                   genLevelWeights, recoLevelWeights,
                   weight, fullStepName);
    }
    if(!stepExists) return;
    
    // Fill the MVA variables of the specific treeHandler
    std::vector<MvaVariablesBase*>& mvaVariables = m_stepMvaVariables_.at(step);
    this->fillVariables(eventMetadata,
                        recoObjects, commonGenObjects,
                        topGenObjects, higgsGenObjects,
                        kinematicReconstructionSolutions,
                        recoObjectIndices, genObjectIndices,
                        genLevelWeights, recoLevelWeights,
                        weight, step,
                        mvaVariables);
}



void MvaTreeHandlerBase::writeTrees(const TString& outputFilename,
                                    const Channel::Channel& channel, const Systematic::Systematic& systematic)
{
    // Create output file for MVA tree
    TString f_savename = common::assignFolder(mvaInputDir_, channel, systematic);
    f_savename.Append(outputFilename);
    TFile outputFile(f_savename, "RECREATE");
    std::cout<<"\nOutput file for MVA input trees: "<<f_savename<<"\n";
    
    // Produce MVA input TTree and store it in output
    TSelectorList* output = new TSelectorList();
    this->writeTrees(output);
    
    // Write file and cleanup
    TIterator* it = output->MakeIterator();
    while(TObject* obj = it->Next()){
        obj->Write();
    }
    outputFile.Close();
    output->SetOwner();
    output->Clear();
}



void MvaTreeHandlerBase::writeTrees(TSelectorList* output)
{
    std::cout<<"--- Beginning production of MVA input trees\n";
    
    // Set pointer to output, so that TTrees are owned by it
    selectorList_ = output;
    
    std::map<TString, TTree*> m_stepTree;
    for(const auto& stepMvaVariables : m_stepMvaVariables_){
        const TString& step = stepMvaVariables.first;
        const std::vector<MvaVariablesBase*>& v_mvaVariables = stepMvaVariables.second;
        TTree* tree = m_stepTree[step];
        tree = this->store(new TTree(prefix_+"mvaVariables"+step, prefix_+"mvaVariables"));
        this->createAndFillBranches(tree, v_mvaVariables);
    }
    
    std::cout<<"MVA variables multiplicity per step (step, no. of entries):\n";
    for(auto vars : m_stepMvaVariables_){
        std::cout<<"\t"<<vars.first<<" , "<<vars.second.size()<<"\n";
    }
    
    std::cout<<"=== Finishing production of MVA input trees\n\n";
}



void MvaTreeHandlerBase::createAndFillBranches(TTree*, const std::vector<MvaVariablesBase*>&)const
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method createAndFillBranches() in MvaTreeHandlerBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



void MvaTreeHandlerBase::createBranch(TTree* tree, const MvaVariableInt& variable)const
{
    std::string name(variable.name());
    std::string nameType(name);
    nameType.append("/").append(variable.type());
    tree->Branch(name.data(), (Long_t)&variable.value_, nameType.data());
}



void MvaTreeHandlerBase::createBranch(TTree* tree, const MvaVariableFloat& variable)const
{
    std::string name(variable.name());
    std::string nameType(name);
    nameType.append("/").append(variable.type());
    tree->Branch(name.c_str(), (Long_t)&variable.value_, nameType.data());
}



void MvaTreeHandlerBase::importTrees(const TString& f_savename, const TString& prefix)
{
    std::cout<<"--- Beginning import of TTrees with MVA variables\n";
    
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
        tth::nameStepPairs(f_savename, treeName);
    
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
        this->importBranches(tree, m_stepMvaVariables_[nameStepPair.second]);
        tree = 0;
        
        std::cout<<"Step, Number of entries: "<<nameStepPair.second<<" , "<<m_stepMvaVariables_[nameStepPair.second].size()<<"\n";
    }
    inputFile->Close();
    
    std::cout<<"=== Finishing import of TTrees with MVA variables\n\n";
}



void MvaTreeHandlerBase::importBranches(TTree*, std::vector<MvaVariablesBase*>&)const
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method importBranches() in MvaTreeHandlerBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



void MvaTreeHandlerBase::importBranch(TTree* tree, MvaVariableInt& variable)const
{
    tree->SetBranchAddress(variable.name().data(), &variable.value_, &variable.branch_);
}



void MvaTreeHandlerBase::importBranch(TTree* tree, MvaVariableFloat& variable)const
{
    tree->SetBranchAddress(variable.name().data(), &variable.value_, &variable.branch_);
}



const std::map<TString, std::vector<MvaVariablesBase*> >& MvaTreeHandlerBase::stepMvaVariablesMap()const
{
    return m_stepMvaVariables_;
}



MvaTreePlotterBase* MvaTreeHandlerBase::setPlotter(const std::map<TString, std::vector<MvaVariablesBase*> >&,
                                                   const bool)const
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method setPlotter() in MvaTreeHandlerBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



void MvaTreeHandlerBase::clear()
{
    selectorList_ = 0;
    
    for(auto& stepMvaVariables : m_stepMvaVariables_){
        MvaVariablesBase::clearVariables(stepMvaVariables.second);
    }
    m_stepMvaVariables_.clear();
}



void MvaTreeHandlerBase::fillVariables(const EventMetadata&,
                                       const RecoObjects&, const CommonGenObjects&,
                                       const TopGenObjects&, const HiggsGenObjects&,
                                       const KinematicReconstructionSolutions&,
                                       const tth::RecoObjectIndices&, const tth::GenObjectIndices&,
                                       const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                       const double&, const TString&,
                                       std::vector<MvaVariablesBase*>&)const
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method fillVariables() in MvaTreeHandlerBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}











