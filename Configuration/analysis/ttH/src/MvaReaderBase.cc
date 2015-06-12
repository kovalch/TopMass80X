#include <iostream>
#include <fstream>

#include <TString.h>
#include <TH2.h>
#include <TMVA/Reader.h>
#include <Rtypes.h>

#include "MvaReaderBase.h"
#include "MvaVariablesBase.h"





MvaReaderBase::MvaReaderBase(const TString& mvaMethod):
mvaMethod_(mvaMethod),
mvaWeightsReader_(0),
mvaVariables_(0)
{}



void MvaReaderBase::book(const TString& mvaWeightsFilename){
    std::cout<<"--- Beginning setting up MVA weights from file\n";
    
    ifstream inputFile(mvaWeightsFilename);
    if(inputFile.fail()){
        std::cout<<"Input file containing MVA weights not found: "<<mvaWeightsFilename
                 <<"\n...running without MVA weights, i.e. setting them all to 1.\n";
    }
    else{
        mvaWeightsReader_ = new TMVA::Reader("Color");
        this->bookVariables(mvaWeightsReader_, mvaVariables_);
        mvaWeightsReader_->BookMVA(mvaMethod_, mvaWeightsFilename);
    }
    
    std::cout<<"=== Finishing setting up MVA weights from file\n\n";
}



void MvaReaderBase::bookVariables(TMVA::Reader* const, MvaVariablesBase*&)const
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method bookVariables() in MvaReaderBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



void MvaReaderBase::addVariable(TMVA::Reader* const mvaWeightsReader, MvaVariableInt& variable)const
{
    mvaWeightsReader->AddVariable(variable.name().data(), &variable.valueFloat_);
}



void MvaReaderBase::addVariable(TMVA::Reader* const mvaWeightsReader, MvaVariableFloat& variable)const
{
    mvaWeightsReader->AddVariable(variable.name().data(), &variable.value_);
}



void MvaReaderBase::addSpectator(TMVA::Reader* const mvaWeightsReader, MvaVariableInt& variable)const
{
    mvaWeightsReader->AddSpectator(variable.name().data(), &variable.valueFloat_);
}



void MvaReaderBase::addSpectator(TMVA::Reader* const mvaWeightsReader, MvaVariableFloat& variable)const
{
    mvaWeightsReader->AddSpectator(variable.name().data(), &variable.value_);
}



void MvaReaderBase::clear()
{
    if(mvaWeightsReader_){
        mvaWeightsReader_->Clear();
        mvaWeightsReader_->Delete();
        mvaWeightsReader_ = 0;
    }
    if(mvaVariables_){
        delete mvaVariables_;
        mvaVariables_ = 0;
    }
}



std::vector<float> MvaReaderBase::mvaWeights(const std::vector<MvaVariablesBase*>& v_mvaVariables)const
{
    std::vector<float> result;
    
    for(MvaVariablesBase* mvaVariables : v_mvaVariables){
        const float weight = this->mvaWeight(mvaVariables);
        result.push_back(weight);
    }
    
    return result;
}



float MvaReaderBase::mvaWeight(const MvaVariablesBase* const mvaVariablesTmp)const
{
    if(!mvaWeightsReader_) return 1.f;
    this->readVariables(mvaVariables_, mvaVariablesTmp);
    const float weight = mvaWeightsReader_->EvaluateMVA(mvaMethod_);
    return weight;
}



void MvaReaderBase::readVariables(MvaVariablesBase* const, const MvaVariablesBase* const)const
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method readVariables() in MvaReaderBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



float MvaReaderBase::combinedWeight(const TH2*const weights2d, const float weight1, const float weight2)
{
    Int_t binX = weights2d->GetXaxis()->FindBin(weight1);
    Int_t binY = weights2d->GetYaxis()->FindBin(weight2);
    const double integral2d = weights2d->Integral(0, binX, 0, binY);
    
    return (static_cast<float>(integral2d));
}



std::vector<float> MvaReaderBase::combinedWeights(const TH2* const weights2d,
                                                  const std::vector<float>& v_weight1,
                                                  const std::vector<float>& v_weight2)
{
    if(v_weight1.size() != v_weight2.size()){
        std::cerr<<"ERROR in MvaReaderBase::combinedWeights()! The two weight vectors are of different size\n...break\n"<<std::endl;
        exit(873);
    }
    if(!weights2d){
        std::cerr<<"ERROR in MvaReaderBase::combinedWeights()! No valid histogram is passed\n...break\n"<<std::endl;
        exit(874);
    }
    
    std::vector<float> result;
    for(size_t iWeight = 0; iWeight < v_weight1.size(); ++iWeight){
        const float combinedWeight = MvaReaderBase::combinedWeight(weights2d, v_weight1.at(iWeight), v_weight2.at(iWeight));
        result.push_back(combinedWeight);
    }
    return result;
}










