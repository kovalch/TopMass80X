#include <iostream>

#include <TString.h>
#include <TMVA/Reader.h>

#include "MvaReaderEventClassification.h"
#include "MvaVariablesBase.h"
#include "MvaVariablesEventClassification.h"
#include "higgsUtils.h"





MvaReaderEventClassification::MvaReaderEventClassification(const TString& mvaMethod, const TString& stepInTraining):
MvaReaderBase(mvaMethod),
stepInTraining_(stepInTraining)
{
    std::cout<<"--- Beginning setting up MVA reader for event classification\n";
    std::cout<<"=== Finishing setting up MVA reader for event classification\n\n";
}



void MvaReaderEventClassification::bookVariables(TMVA::Reader* const mvaWeightsReader, MvaVariablesBase*& mvaVariables)const
{
    mvaVariables = new MvaVariablesEventClassification();
    MvaVariablesEventClassification* const mvaVariablesEventClassification = dynamic_cast<MvaVariablesEventClassification* const>(mvaVariables);
    if(!mvaVariablesEventClassification){
        std::cerr<<"ERROR in MvaReaderEventClassification::variablesForReader()! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
        exit(395);
    }
    
    const TString category = tth::extractJetCategory(stepInTraining_);
    
    if(category == tth::categoryName(0)){
        //this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->btagDiscriminatorAverage_tagged_); // Used in Run1, but not good here
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->btagDiscriminatorAverage_untagged_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->ptSum_jets_leptons_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->minDeltaR_jet_jet_);
        //this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->multiplicity_higgsLikeDijet15_); // Might be included
    }
    else if(category == tth::categoryName(1)){
        //this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->mass_higgsLikeDijet2_);
        //this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->minDeltaR_jet_jet_);
        //this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->ptSum_jets_leptons_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->btagDiscriminatorAverage_tagged_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->multiplicity_higgsLikeDijet15_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->mass_higgsLikeDijet_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->multiplicity_jets_);
    }
    else if(category == tth::categoryName(2)){
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->btagDiscriminatorAverage_untagged_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->minDeltaR_jet_jet_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->ptSum_jets_leptons_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->multiplicity_jets_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->mass_higgsLikeDijet_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->mass_higgsLikeDijet2_);
    }
    else if(category == tth::categoryName(3)){
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->multiplicity_jets_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->btagDiscriminatorAverage_tagged_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->multiplicity_higgsLikeDijet15_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->mass_higgsLikeDijet2_);
        //this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->ptSum_jets_leptons_);
    }
    else if(category == tth::categoryName({0, 1, 2, 3})){
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->multiplicity_jets_);
        //this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->btagDiscriminatorAverage_tagged_);
        //this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->btagDiscriminatorAverage_untagged_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->minDeltaR_jet_jet_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->ptSum_jets_leptons_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->multiplicity_higgsLikeDijet15_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->mass_higgsLikeDijet_);
        this->addVariable(mvaWeightsReader, mvaVariablesEventClassification->mass_higgsLikeDijet2_);
    }
    else{
        std::cerr<<"Error in MvaFactoryEventClassification::configureFactory2()! No input variables defined for category: "
                 <<category<<"\n...break\n"<<std::endl;
        exit(281);
    }
}



void MvaReaderEventClassification::readVariables(MvaVariablesBase* const mvaVariables, const MvaVariablesBase* const mvaVariablesTmp)const
{
    const MvaVariablesEventClassification* const mvaVariablesEventClassificationTmp = dynamic_cast<const MvaVariablesEventClassification* const>(mvaVariablesTmp);
    if(!mvaVariablesEventClassificationTmp){
        std::cerr<<"ERROR in MvaReaderEventClassification::readVariables()! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
        exit(395);
    }
    
    MvaVariablesEventClassification* const mvaVariablesEventClassification = dynamic_cast<MvaVariablesEventClassification* const>(mvaVariables);
    if(!mvaVariablesEventClassification){
        std::cerr<<"ERROR in MvaReaderEventClassification::readVariables()! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
        exit(395);
    }
    
    *mvaVariablesEventClassification = *mvaVariablesEventClassificationTmp;
}










