#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCut.h>
#include <Rtypes.h>
#include <TMVA/Factory.h>
#include <TMVA/Types.h>

#include "MvaFactoryTopJets.h"
#include "MvaVariablesTopJets.h"
#include "mvaSetup.h"





MvaFactoryTopJets::MvaFactoryTopJets(const TString& mvaOutputDir, const TString& weightFileDir,
                                     const TString& mergedTreesFileName):
MvaFactoryBase(mvaOutputDir, weightFileDir, mergedTreesFileName)
{
    std::cout<<"--- Beginning setting up MVA factory for training\n";
    std::cout<<"=== Finishing setting up MVA factory for training\n\n";
}



void MvaFactoryTopJets::trainMva(TFile* const treeFile, const mvaSetup::MvaSet& mvaSet, const TString& step)const
{
    // MVA for correct dijet combinations
    constexpr const char* methodPrefixCorrect = "correct";
    const TCut cutSignalCorrect = "correctCombination != 0";
    const TCut cutBackgroundCorrect = "correctCombination == 0 && swappedCombination == 0";
    
    // MVA for swapped dijet combinations
    constexpr const char* methodPrefixSwapped = "swapped";
    const TCut cutSignalSwapped = "swappedCombination != 0";
    const TCut cutBackgroundSwapped = "correctCombination == 0 && swappedCombination == 0";
    
    TTree* const treeTraining = (TTree*)treeFile->Get("training"+step);
    TTree* const treeTesting = (TTree*)treeFile->Get("testing"+step);
    
    if(mvaSet.v_mvaConfigCorrect_.size())
        this->runMva(methodPrefixCorrect, cutSignalCorrect, cutBackgroundCorrect, treeTraining, treeTesting, mvaSet.v_mvaConfigCorrect_, step);
    if(mvaSet.v_mvaConfigSwapped_.size())
        this->runMva(methodPrefixSwapped, cutSignalSwapped, cutBackgroundSwapped, treeTraining, treeTesting, mvaSet.v_mvaConfigSwapped_, step);
}



void MvaFactoryTopJets::configureFactory(TMVA::Factory* const factory,
                                         const TCut& cutSignal, const TCut& cutBackground,
                                         TTree* const treeTraining, TTree* const treeTesting)const
{
    MvaVariablesTopJets mvaVariablesTopJets;
    
    // Set all branches of MVA input which should be used for training
    this->addVariable(factory, mvaVariablesTopJets.jetChargeDiff_);
    this->addVariable(factory, mvaVariablesTopJets.meanDeltaPhi_b_met_);
    this->addVariable(factory, mvaVariablesTopJets.massDiff_recoil_bbbar_);
    this->addVariable(factory, mvaVariablesTopJets.pt_b_antiLepton_);
    this->addVariable(factory, mvaVariablesTopJets.pt_antiB_lepton_);
    this->addVariable(factory, mvaVariablesTopJets.deltaR_b_antiLepton_);
    this->addVariable(factory, mvaVariablesTopJets.deltaR_antiB_lepton_);
    this->addVariable(factory, mvaVariablesTopJets.deltaPhi_antiBLepton_bAntiLepton_);
    this->addVariable(factory, mvaVariablesTopJets.massDiff_fullBLepton_bbbar_);
    this->addVariable(factory, mvaVariablesTopJets.meanMt_b_met_);
    this->addVariable(factory, mvaVariablesTopJets.massSum_antiBLepton_bAntiLepton_);
    this->addVariable(factory, mvaVariablesTopJets.massDiff_antiBLepton_bAntiLepton_);
    
    // Set all branches of MVA input which should NOT be used for training,
    // but are needed otherwise (e.g. for defining separation cuts), or are under study
    this->addSpectator(factory, mvaVariablesTopJets.bQuarkRecoJetMatched_);
    this->addSpectator(factory, mvaVariablesTopJets.correctCombination_);
    this->addSpectator(factory, mvaVariablesTopJets.swappedCombination_);
    this->addSpectator(factory, mvaVariablesTopJets.btagDiscriminatorSum_);
    
    // Set global weights for individual input
    constexpr Double_t signalWeight = 1.;
    constexpr Double_t backgroundWeight = 1.;
    
    // Register the training trees
    factory->AddSignalTree(treeTraining, signalWeight, TMVA::Types::kTraining);
    factory->AddBackgroundTree(treeTraining, backgroundWeight, TMVA::Types::kTraining);
    
    // Register the testing trees
    factory->AddSignalTree(treeTesting, signalWeight, TMVA::Types::kTesting);
    factory->AddBackgroundTree(treeTesting, backgroundWeight, TMVA::Types::kTesting);
    
    // Set the branch from which the event weight is taken
    factory->SetSignalWeightExpression(mvaVariablesTopJets.eventWeight_.name());
    factory->SetBackgroundWeightExpression(mvaVariablesTopJets.eventWeight_.name());
    
    // Prepare the training and test trees
    factory->PrepareTrainingAndTestTree(cutSignal, cutBackground,
                                        "SplitMode=Block:SplitSeed=0:NormMode=NumEvents:!V" );
}








