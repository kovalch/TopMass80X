#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <cstdlib>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <Rtypes.h>
#include <TCut.h>
#include <TMVA/Tools.h>
#include <TMVA/Config.h>
#include <TMVA/Factory.h>
#include <TMVA/Types.h>

#include "MvaFactory.h"
#include "MvaVariablesTopJets.h"
#include "mvaSetup.h"
#include "higgsUtils.h"




MvaFactory::MvaFactory(const TString& mvaOutputDir, const char* weightFileDir,
                       const std::vector<mvaSetup::MvaSet>& v_mvaSet,
                       const TString& mergedTreesFileName):
mvaOutputDir_(mvaOutputDir),
weightFileDir_(weightFileDir),
mergedTreesFile_(0)
{
    std::cout<<"--- Beginning MVA training\n";
    
    // Open input file
    mergedTreesFile_ = TFile::Open(mergedTreesFileName);
    if(!mergedTreesFile_){
        std::cerr<<"ERROR in constructor of MvaFactory! File containing the merged trees not found\n...break\n"<<std::endl;
        exit(123);
    }
    
    // Run all MVA trainings for all given sets
    this->train(v_mvaSet);
    
    // Cleanup
    mergedTreesFile_->Close();
    this->clear();
    
    std::cout<<"=== Finishing MVA training\n\n";
}



void MvaFactory::clear()
{
    weightFileDir_ = 0;
    mergedTreesFile_ = 0;
}



void MvaFactory::train(const std::vector<mvaSetup::MvaSet>& v_mvaSet)
{
    // MVA for correct dijet combinations
    constexpr const char* methodPrefixCorrect = "correct";
    const TCut cutSignalCorrect = "correctCombination != 0";
    const TCut cutBackgroundCorrect = "correctCombination == 0 && swappedCombination == 0";
    
    // MVA for swapped dijet combinations
    constexpr const char* methodPrefixSwapped = "swapped";
    const TCut cutSignalSwapped = "swappedCombination != 0";
    const TCut cutBackgroundSwapped = "correctCombination == 0 && swappedCombination == 0";
    
    // Loop over the steps/categories and train MVAs
    for(const auto& mvaSet : v_mvaSet){
        const TString mergedStepName = tth::stepName(mvaSet.step_, mvaSet.v_category_);
        
        TTree* treeTraining = (TTree*)mergedTreesFile_->Get("training"+mergedStepName);
        TTree* treeTesting = (TTree*)mergedTreesFile_->Get("testing"+mergedStepName);
        
        if(mvaSet.v_mvaConfigCorrect_.size())
            this->runMva(methodPrefixCorrect, cutSignalCorrect, cutBackgroundCorrect, treeTraining, treeTesting, mvaSet.v_mvaConfigCorrect_, mergedStepName);
        if(mvaSet.v_mvaConfigSwapped_.size())
            this->runMva(methodPrefixSwapped, cutSignalSwapped, cutBackgroundSwapped, treeTraining, treeTesting, mvaSet.v_mvaConfigSwapped_, mergedStepName);
        
        
    }
    
}



void MvaFactory::runMva(const char* methodPrefix, const TCut& cutSignal, const TCut& cutBackground,
                        TTree* treeTraining, TTree* treeTesting,
                        const std::vector<mvaSetup::MvaConfig>& v_mvaConfig,
                        const TString& stepName)
{
    
    // Get a TMVA instance
    TMVA::Tools::Instance();
    
    // Create a ROOT output file for TMVA
    TString mvaOutputFilename(mvaOutputDir_);
    mvaOutputFilename.Append(methodPrefix);
    mvaOutputFilename.Append(stepName);
    mvaOutputFilename.Append(".root");
    TFile* outputFile = TFile::Open(mvaOutputFilename, "RECREATE");
    
    // Set the output directory for the weights (if not specified, default is "weights")
    TString mvaOutputWeightsFilename(mvaOutputDir_);
    mvaOutputWeightsFilename.Append(weightFileDir_);
    (TMVA::gConfig().GetIONames()).fWeightFileDir = mvaOutputWeightsFilename;
    
    // Create the factory
    TMVA::Factory* factory(0);
    TString factoryName(methodPrefix);
    factoryName.Append(stepName);
    factory = new TMVA::Factory(factoryName, outputFile, "!V:!Silent");
    
    // Set all branches of MVA input which should be used for training
    MvaVariablesTopJets mvaTopJetsVariables;
    
    this->addVariable(factory, mvaTopJetsVariables.jetChargeDiff_);
    this->addVariable(factory, mvaTopJetsVariables.meanDeltaPhi_b_met_);
    this->addVariable(factory, mvaTopJetsVariables.massDiff_recoil_bbbar_);
    this->addVariable(factory, mvaTopJetsVariables.pt_b_antiLepton_);
    this->addVariable(factory, mvaTopJetsVariables.pt_antiB_lepton_);
    this->addVariable(factory, mvaTopJetsVariables.deltaR_b_antiLepton_);
    this->addVariable(factory, mvaTopJetsVariables.deltaR_antiB_lepton_);
    this->addVariable(factory, mvaTopJetsVariables.btagDiscriminatorSum_);
    this->addVariable(factory, mvaTopJetsVariables.deltaPhi_antiBLepton_bAntiLepton_);
    this->addVariable(factory, mvaTopJetsVariables.massDiff_fullBLepton_bbbar_);
    this->addVariable(factory, mvaTopJetsVariables.meanMt_b_met_);
    this->addVariable(factory, mvaTopJetsVariables.massSum_antiBLepton_bAntiLepton_);
    this->addVariable(factory, mvaTopJetsVariables.massDiff_antiBLepton_bAntiLepton_);
    
    // Set all branches of MVA input which should NOT be used for training,
    // but are needed otherwise (e.g. for defining separation cuts)
    this->addSpectator(factory, mvaTopJetsVariables.bQuarkRecoJetMatched_);
    this->addSpectator(factory, mvaTopJetsVariables.correctCombination_);
    this->addSpectator(factory, mvaTopJetsVariables.swappedCombination_);
    
    // Set global weights for individual input
    Double_t signalWeight = 1.;
    Double_t backgroundWeight = 1.;
    
    // Register the training trees
    factory->AddSignalTree(treeTraining, signalWeight, TMVA::Types::kTraining);
    factory->AddBackgroundTree(treeTraining, backgroundWeight, TMVA::Types::kTraining);
    
    // Register the testing trees
    factory->AddSignalTree(treeTesting, signalWeight, TMVA::Types::kTesting);
    factory->AddBackgroundTree(treeTesting, backgroundWeight, TMVA::Types::kTesting);
    
    // Set the branch from which the event weight is taken
    factory->SetSignalWeightExpression(mvaTopJetsVariables.eventWeight_.name());
    factory->SetBackgroundWeightExpression(mvaTopJetsVariables.eventWeight_.name());
    
    // Prepare the training and test trees
    factory->PrepareTrainingAndTestTree(cutSignal, cutBackground,
                                        "SplitMode=Block:SplitSeed=0:NormMode=NumEvents:!V" );
    
    // Book the MVA methods (e.g. boosted decision tree with specific setup)
    for(const auto& mvaConfig : v_mvaConfig){
        const TString methodTitle(mvaConfig.methodAppendix_);
        factory->BookMethod(mvaConfig.mvaType_,
                            methodTitle,
                            mvaConfig.options_);
    }
    
    // Run factory
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();
    
    // Cleanup
    outputFile->Close();
    delete factory;
}



void MvaFactory::addVariable(TMVA::Factory* factory, MvaVariableInt& variable)
{
    factory->AddVariable(variable.name().data(), *variable.type());
}



void MvaFactory::addVariable(TMVA::Factory* factory, MvaVariableFloat& variable)
{
    factory->AddVariable(variable.name().data(), *variable.type());
}



void MvaFactory::addSpectator(TMVA::Factory* factory, MvaVariableInt& variable)
{
    factory->AddSpectator(variable.name().data(), *variable.type());
}



void MvaFactory::addSpectator(TMVA::Factory* factory, MvaVariableFloat& variable)
{
    factory->AddSpectator(variable.name().data(), *variable.type());
}










