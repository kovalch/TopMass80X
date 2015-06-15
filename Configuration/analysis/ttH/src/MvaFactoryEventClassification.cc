#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCut.h>
#include <Rtypes.h>
#include <TMVA/Factory.h>
#include <TMVA/Types.h>

#include "MvaFactoryEventClassification.h"
#include "MvaVariablesEventClassification.h"
#include "Sample.h"
#include "Samples.h"
#include "mvaSetup.h"
#include "higgsUtils.h"





MvaFactoryEventClassification::MvaFactoryEventClassification(const TString& mvaOutputDir, const TString& weightFileDir,
                                                             const Samples& samples,
                                                             const bool inOneFactory):
MvaFactoryBase(mvaOutputDir, weightFileDir, samples, inOneFactory)
{
    std::cout<<"--- Beginning setting up MVA factory for training\n";
    std::cout<<"=== Finishing setting up MVA factory for training\n\n";
}



void MvaFactoryEventClassification::trainMva2(const mvaSetup::MvaSet& mvaSet, const TString& outputFolder,
                                              const std::vector<Sample>& v_sample, const std::vector<double>& v_weight,
                                              std::ofstream& trainingNameFile, const TString& step)const
{
    if(!mvaSet.v_mvaConfigCorrect_.size()){
        std::cerr<<"ERROR in MvaFactoryEventClassification::trainMva()! No training configs defined\n...break\n"<<std::endl;
        exit(381);
    }
    if(mvaSet.v_mvaConfigSwapped_.size()){
        std::cerr<<"ERROR in MvaFactoryEventClassification::trainMva()! Trainings for 'swapped' defined, but option is here invalid\n...break\n"<<std::endl;
        exit(381);
    }
    this->runMva2(mvaSet.v_mvaConfigCorrect_, outputFolder, v_sample, v_weight, trainingNameFile, step);
}



void MvaFactoryEventClassification::configureFactory2(TMVA::Factory* const factory,
                                                     const TCut& cutSignal, const TCut& cutBackground,
                                                     const std::vector<Sample>& v_sample, const std::vector<double>& v_weight,
                                                     const TString& stepName)const
{
    MvaVariablesEventClassification mvaVariablesEventClassification;
    
    const TString category = tth::extractJetCategory(stepName);
    
    // Set all branches of MVA input which should be used for training (as variable)
    // Set all branches of MVA input which should NOT be used for training (as spectator),
    // but are needed otherwise (e.g. for defining separation cuts)
    if(category == tth::categoryName(0)){
        this->addVariable(factory, mvaVariablesEventClassification.btagDiscriminatorAverage_tagged_);
        this->addVariable(factory, mvaVariablesEventClassification.btagDiscriminatorAverage_untagged_);
        this->addVariable(factory, mvaVariablesEventClassification.ptSum_jets_leptons_);
        this->addVariable(factory, mvaVariablesEventClassification.minDeltaR_jet_jet_);
    }
    else if(category == tth::categoryName(1)){
        this->addVariable(factory, mvaVariablesEventClassification.btagDiscriminatorAverage_tagged_);
        this->addVariable(factory, mvaVariablesEventClassification.minDeltaR_jet_jet_);
        this->addVariable(factory, mvaVariablesEventClassification.ptSum_jets_leptons_);
        this->addVariable(factory, mvaVariablesEventClassification.multiplicity_higgsLikeDijet15_);
        this->addVariable(factory, mvaVariablesEventClassification.mass_higgsLikeDijet_);
    }
    else if(category == tth::categoryName(2)){
        this->addVariable(factory, mvaVariablesEventClassification.btagDiscriminatorAverage_untagged_);
        this->addVariable(factory, mvaVariablesEventClassification.minDeltaR_jet_jet_);
        this->addVariable(factory, mvaVariablesEventClassification.ptSum_jets_leptons_);
        this->addVariable(factory, mvaVariablesEventClassification.multiplicity_jets_);
        this->addVariable(factory, mvaVariablesEventClassification.mass_higgsLikeDijet_);
        this->addVariable(factory, mvaVariablesEventClassification.mass_higgsLikeDijet2_);
    }
    else if(category == tth::categoryName(3)){
        this->addVariable(factory, mvaVariablesEventClassification.multiplicity_jets_);
        this->addVariable(factory, mvaVariablesEventClassification.btagDiscriminatorAverage_tagged_);
        this->addVariable(factory, mvaVariablesEventClassification.btagDiscriminatorAverage_untagged_);
        this->addVariable(factory, mvaVariablesEventClassification.minDeltaR_jet_jet_);
        this->addVariable(factory, mvaVariablesEventClassification.ptSum_jets_leptons_);
        this->addVariable(factory, mvaVariablesEventClassification.multiplicity_higgsLikeDijet15_);
        this->addVariable(factory, mvaVariablesEventClassification.mass_higgsLikeDijet_);
        this->addVariable(factory, mvaVariablesEventClassification.mass_higgsLikeDijet2_);
    }
    else{
        std::cerr<<"Error in MvaFactoryEventClassification::configureFactory2()! No input variables defined for category: "
                 <<category<<"\n...break\n"<<std::endl;
        exit(281);
    }
    
    // Add samples to factory
    for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
        const Sample& sample = v_sample.at(iSample);
        if(sample.sampleType() == Sample::data) continue;
        const double& weight = v_weight.at(iSample);
        TString treename("event_mvaVariables");
        treename.Append(stepName);
        // FIXME: Need to check if file exists, or already done when setting up samples?
        TFile* treeFile = TFile::Open(sample.inputFile());
        TTree* const tree = (TTree*)treeFile->Get(treename);
        if(!tree){
            std::cerr<<"ERROR. Tree not found\n"<<std::endl;
            exit(12);
        }
        //std::cout<<"\nsize: "<<tree->GetEntriesFast()<<"\n\n";
        if(!tree->GetEntriesFast()) continue;
        if(sample.sampleType() == Sample::ttHbb) factory->AddSignalTree(tree, weight);
        else factory->AddBackgroundTree(tree, weight);
        //std::cout<<"\n\nweight: "<<weight<<" , "<<sample.legendEntry()<<"\n\n";
    }
    
    // Set the branch from which the event weight is taken
    factory->SetSignalWeightExpression(mvaVariablesEventClassification.eventWeight_.name());
    factory->SetBackgroundWeightExpression(mvaVariablesEventClassification.eventWeight_.name());
    
    // Prepare the training and test trees
    factory->PrepareTrainingAndTestTree(cutSignal, cutBackground,
                                        "SplitMode=Block:SplitSeed=0:NormMode=NumEvents:!V" );
}








