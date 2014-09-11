#include <iostream>

#include <TString.h>
#include <TMVA/Reader.h>

#include "MvaReaderTopJets.h"
#include "MvaVariablesBase.h"
#include "MvaVariablesTopJets.h"





MvaReaderTopJets::MvaReaderTopJets(const TString& mvaMethod):
MvaReaderBase(mvaMethod)
{
    std::cout<<"--- Beginning setting up MVA reader for top jets\n";
    std::cout<<"=== Finishing setting up MVA reader for top jets\n\n";
}



void MvaReaderTopJets::bookVariables(TMVA::Reader* const mvaWeightsReader, MvaVariablesBase*& mvaVariables)const
{
    mvaVariables = new MvaVariablesTopJets();
    MvaVariablesTopJets* const mvaVariablesTopJets = dynamic_cast<MvaVariablesTopJets* const>(mvaVariables);
    if(!mvaVariablesTopJets){
        std::cerr<<"ERROR in MvaReaderTopJets::variablesForReader()! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
        exit(395);
    }
    
    this->addVariable(mvaWeightsReader, mvaVariablesTopJets->jetChargeDiff_);
    this->addVariable(mvaWeightsReader, mvaVariablesTopJets->meanDeltaPhi_b_met_);
    this->addVariable(mvaWeightsReader, mvaVariablesTopJets->massDiff_recoil_bbbar_);
    this->addVariable(mvaWeightsReader, mvaVariablesTopJets->pt_b_antiLepton_);
    this->addVariable(mvaWeightsReader, mvaVariablesTopJets->pt_antiB_lepton_);
    this->addVariable(mvaWeightsReader, mvaVariablesTopJets->deltaR_b_antiLepton_);
    this->addVariable(mvaWeightsReader, mvaVariablesTopJets->deltaR_antiB_lepton_);
    this->addVariable(mvaWeightsReader, mvaVariablesTopJets->deltaPhi_antiBLepton_bAntiLepton_);
    this->addVariable(mvaWeightsReader, mvaVariablesTopJets->massDiff_fullBLepton_bbbar_);
    this->addVariable(mvaWeightsReader, mvaVariablesTopJets->meanMt_b_met_);
    this->addVariable(mvaWeightsReader, mvaVariablesTopJets->massSum_antiBLepton_bAntiLepton_);
    this->addVariable(mvaWeightsReader, mvaVariablesTopJets->massDiff_antiBLepton_bAntiLepton_);
    
    this->addSpectator(mvaWeightsReader, mvaVariablesTopJets->bQuarkRecoJetMatched_);
    this->addSpectator(mvaWeightsReader, mvaVariablesTopJets->correctCombination_);
    this->addSpectator(mvaWeightsReader, mvaVariablesTopJets->swappedCombination_);
    this->addSpectator(mvaWeightsReader, mvaVariablesTopJets->btagDiscriminatorSum_);
}



void MvaReaderTopJets::readVariables(MvaVariablesBase* const mvaVariables, const MvaVariablesBase* const mvaVariablesTmp)const
{
    const MvaVariablesTopJets* const mvaVariablesTopJetsTmp = dynamic_cast<const MvaVariablesTopJets* const>(mvaVariablesTmp);
    if(!mvaVariablesTopJetsTmp){
        std::cerr<<"ERROR in MvaReaderTopJets::readVariables()! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
        exit(395);
    }
    
    MvaVariablesTopJets* const mvaVariablesTopJets = dynamic_cast<MvaVariablesTopJets* const>(mvaVariables);
    if(!mvaVariablesTopJets){
        std::cerr<<"ERROR in MvaReaderTopJets::readVariables()! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
        exit(395);
    }
    
    *mvaVariablesTopJets = *mvaVariablesTopJetsTmp;
}










