#include <iostream>
#include <fstream>

#include <TString.h>
#include <TH2.h>
#include <TMVA/Reader.h>

#include "MvaReader.h"
#include "MvaVariablesTopJets.h"





MvaReader::MvaReader(const char* mvaWeightsFile):
mvaWeightsReader_(0),
mvaVariables_(0)
{
    std::cout<<"--- Beginning setting up MVA weights from file\n";
    
    ifstream inputFile(mvaWeightsFile);
    if(inputFile.fail()){
        std::cout<<"Input file containing MVA weights not found: "<<mvaWeightsFile
                 <<"\n...running without MVA weights, i.e. setting them all to 1.\n";
    }
    else{
        mvaWeightsReader_ = new TMVA::Reader("Color");
        
        mvaVariables_ = new MvaVariablesTopJets();
        MvaVariablesTopJets* mvaVariablesTopJets = dynamic_cast<MvaVariablesTopJets*>(mvaVariables_);
        if(!mvaVariablesTopJets){
            std::cerr<<"ERROR in constructor of MvaReader! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
            exit(395);
        }
        
        this->addVariable(mvaVariablesTopJets->jetChargeDiff_);
        this->addVariable(mvaVariablesTopJets->meanDeltaPhi_b_met_);
        this->addVariable(mvaVariablesTopJets->massDiff_recoil_bbbar_);
        this->addVariable(mvaVariablesTopJets->pt_b_antiLepton_);
        this->addVariable(mvaVariablesTopJets->pt_antiB_lepton_);
        this->addVariable(mvaVariablesTopJets->deltaR_b_antiLepton_);
        this->addVariable(mvaVariablesTopJets->deltaR_antiB_lepton_);
        this->addVariable(mvaVariablesTopJets->btagDiscriminatorSum_);
        this->addVariable(mvaVariablesTopJets->deltaPhi_antiBLepton_bAntiLepton_);
        this->addVariable(mvaVariablesTopJets->massDiff_fullBLepton_bbbar_);
        this->addVariable(mvaVariablesTopJets->meanMt_b_met_);
        this->addVariable(mvaVariablesTopJets->massSum_antiBLepton_bAntiLepton_);
        this->addVariable(mvaVariablesTopJets->massDiff_antiBLepton_bAntiLepton_);
        
        this->addSpectator(mvaVariablesTopJets->bQuarkRecoJetMatched_);
        this->addSpectator(mvaVariablesTopJets->correctCombination_);
        this->addSpectator(mvaVariablesTopJets->swappedCombination_);
        
        // FIXME: what is first argument, should it be "BDTG" or "BDT method" ???
        mvaWeightsReader_->BookMVA("BDT method", mvaWeightsFile);
    }
    
    std::cout<<"=== Finishing setting up MVA weights from file\n\n";
}



void MvaReader::addVariable(MvaVariableInt& variable)
{
    mvaWeightsReader_->AddVariable(variable.name().data(), &variable.value_);
}



void MvaReader::addVariable(MvaVariableFloat& variable)
{
    mvaWeightsReader_->AddVariable(variable.name().data(), &variable.value_);
}



void MvaReader::addSpectator(MvaVariableInt& variable)
{
    mvaWeightsReader_->AddSpectator(variable.name().data(), &variable.value_);
}



void MvaReader::addSpectator(MvaVariableFloat& variable)
{
    mvaWeightsReader_->AddSpectator(variable.name().data(), &variable.value_);
}



void MvaReader::clear()
{
    v_mvaTopJetsVariablesPerEvent_.clear();
    if(mvaWeightsReader_){
        mvaWeightsReader_->Clear();
        mvaWeightsReader_->Delete();
    }
}



std::vector<float> MvaReader::mvaWeights(const std::vector<MvaVariablesBase*>& v_mvaVariables)
{
    std::vector<float> result;
    
    for(MvaVariablesBase* mvaVariablesTmp : v_mvaVariables){
        if(!mvaWeightsReader_){
            result.push_back(1.);
            continue;
        }
        
        MvaVariablesTopJets* mvaVariablesTopJetsTmp = dynamic_cast<MvaVariablesTopJets*>(mvaVariablesTmp);
        if(!mvaVariablesTopJetsTmp){
            std::cerr<<"ERROR in MvaReader::mvaWeights()! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
            exit(395);
        }
        
        MvaVariablesTopJets* const mvaVariablesTopJets = dynamic_cast<MvaVariablesTopJets* const>(mvaVariables_);
        if(!mvaVariablesTopJets){
            std::cerr<<"ERROR in MvaReader::mvaWeights()! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
            exit(395);
        }
        
        *mvaVariablesTopJets = *mvaVariablesTopJetsTmp;
        const float weight = mvaWeightsReader_->EvaluateMVA("BDT method");
        result.push_back(weight);
    }
    
    return result;
}



std::vector<float> MvaReader::combinedWeight(const TH2* weights2d,
                                             const std::vector<float>& weightCorrect,
                                             const std::vector<float>& weightSwapped)
{
    if(weightCorrect.size() != weightSwapped.size()){
        std::cerr<<"ERROR in combinedWeight()! The two weight vectors are of different size\n...break\n"<<std::endl;
        exit(873);
    }
    if(!weights2d){
        std::cerr<<"ERROR in combinedWeight()! No valid histogram is passed\n...break\n"<<std::endl;
        exit(874);
    }
    
    std::vector<float> result;
    for(size_t iWeight = 0; iWeight < weightCorrect.size(); ++iWeight){
        const float mvaWeightCorrect = weightCorrect.at(iWeight);
        const float mvaWeightSwapped = weightSwapped.at(iWeight);
        
        Int_t binX = weights2d->GetXaxis()->FindBin(mvaWeightCorrect);
        Int_t binY = weights2d->GetYaxis()->FindBin(mvaWeightSwapped);
        const double integral2d = weights2d->Integral(0, binX, 0, binY);
        
        result.push_back(static_cast<float>(integral2d));
    }
    
    return result;
}










