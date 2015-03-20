#include <iostream>

#include <TString.h>
#include <TMVA/Reader.h>

#include "MvaReaderJetCharge.h"
#include "MvaVariablesBase.h"
#include "MvaVariablesJetCharge.h"





MvaReaderJetCharge::MvaReaderJetCharge(const TString& mvaMethod):
MvaReaderBase(mvaMethod)
{
    std::cout<<"--- Beginning setting up MVA reader for jet charge\n";
    std::cout<<"=== Finishing setting up MVA reader for jet charge\n\n";
}



void MvaReaderJetCharge::bookVariables(TMVA::Reader* const mvaWeightsReader, MvaVariablesBase*& mvaVariables)const
{
    mvaVariables = new MvaVariablesJetCharge();
    MvaVariablesJetCharge* const mvaVariablesJetCharge = dynamic_cast<MvaVariablesJetCharge* const>(mvaVariables);
    if(!mvaVariablesJetCharge){
        std::cerr<<"ERROR in MvaReaderJetCharge::bookVariables()! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
        exit(395);
    }
    
    this->addVariable(mvaWeightsReader, mvaVariablesJetCharge->trueBJetId_);
    this->addVariable(mvaWeightsReader, mvaVariablesJetCharge->thereIsALeadingLepton_);
    this->addVariable(mvaWeightsReader, mvaVariablesJetCharge->thereIsALeadingMuon_);
    this->addVariable(mvaWeightsReader, mvaVariablesJetCharge->thereIsASecondaryVertex_);
    
    this->addVariable(mvaWeightsReader, mvaVariablesJetCharge->longChargeJet_);
    this->addVariable(mvaWeightsReader, mvaVariablesJetCharge->relChargeJet_);
    this->addVariable(mvaWeightsReader, mvaVariablesJetCharge->leadingTrackPtWeightedCharge_);
    this->addVariable(mvaWeightsReader, mvaVariablesJetCharge->subleadingTrackPtWeightedCharge_);
    this->addVariable(mvaWeightsReader, mvaVariablesJetCharge->thirdleadingTrackPtWeightedCharge_);
    this->addVariable(mvaWeightsReader, mvaVariablesJetCharge->leadingMuonPtWeightedCharge_);
    this->addVariable(mvaWeightsReader, mvaVariablesJetCharge->leadingElectronPtWeightedCharge_);
    this->addVariable(mvaWeightsReader, mvaVariablesJetCharge->trackNumberWeightedJetPt_);
    this->addVariable(mvaWeightsReader, mvaVariablesJetCharge->chargeWeightedTrackId_);
    this->addVariable(mvaWeightsReader, mvaVariablesJetCharge->svChargeWeightedFlightDistance_);
    this->addVariable(mvaWeightsReader, mvaVariablesJetCharge->secondaryVertexCharge_);
    this->addVariable(mvaWeightsReader, mvaVariablesJetCharge->ipSignificanceLeadingTrack_);
}



void MvaReaderJetCharge::readVariables(MvaVariablesBase* const mvaVariables, const MvaVariablesBase* const mvaVariablesTmp)const
{
    const MvaVariablesJetCharge* const mvaVariablesJetChargeTmp = dynamic_cast<const MvaVariablesJetCharge* const>(mvaVariablesTmp);
    if(!mvaVariablesJetChargeTmp){
        std::cerr<<"ERROR in MvaReaderJetCharge::readVariables()! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
        exit(395);
    }
    
    MvaVariablesJetCharge* const mvaVariablesJetCharge = dynamic_cast<MvaVariablesJetCharge* const>(mvaVariables);
    if(!mvaVariablesJetCharge){
        std::cerr<<"ERROR in MvaReaderJetCharge::readVariables()! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
        exit(395);
    }
    
    *mvaVariablesJetCharge = *mvaVariablesJetChargeTmp;
}










