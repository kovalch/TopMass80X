#include <iostream>

#include <TString.h>
#include <TH1.h>
#include <TH1D.h>

#include "MvaTreePlotterEventClassification.h"
#include "MvaVariablesBase.h"
#include "MvaVariablesEventClassification.h"





MvaTreePlotterEventClassification::MvaTreePlotterEventClassification(const std::map<TString, std::vector<MvaVariablesBase*> >& m_stepMvaVariables,
                                                                     const bool separationPowerPlots):
MvaTreePlotterBase("mvaEventP_", m_stepMvaVariables, separationPowerPlots)
{
    std::cout<<"--- Beginning setting up top jets MVA variables plotter\n";
    
    // Check that MVA variables are of correct type for plotter
    for(const auto& stepMvaVariables : m_stepMvaVariables){
        const std::vector<MvaVariablesBase*>& v_mvaVariables = stepMvaVariables.second;
        if(v_mvaVariables.size()){
            const MvaVariablesEventClassification* const mvaVariablesEventClassification = dynamic_cast<const MvaVariablesEventClassification* const>(v_mvaVariables.at(0));
            if(!mvaVariablesEventClassification){
                std::cerr<<"ERROR in constructor of MvaTreePlotterEventClassification! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
                exit(395);
            }
            break;
        }
    }
    
    std::cout<<"=== Finishing setting up top jets MVA variables plotter\n\n";
}



void MvaTreePlotterEventClassification::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    const MvaVariablesEventClassification nameDummy;
    
    name = nameDummy.multiplicity_jets_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+";# jets;# events",12,0,12);
    
    name = nameDummy.btagDiscriminatorAverage_tagged_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+";<d>_{tagged};# events",20,0,1);
    
    name = nameDummy.btagDiscriminatorAverage_untagged_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+";<d>_{untagged};# events",20,0,1);
    
    name = nameDummy.minDeltaR_jet_jet_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+";#DeltaR_{min}^{jj};# events",20,0,4);
    
    name = nameDummy.ptSum_jets_leptons_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+"; #Sump_{T}^{l,j}  [GeV];# jet pairs",20,0,800);
}



void MvaTreePlotterEventClassification::fillHistos(const TString&,
                                                   const std::vector<MvaVariablesBase*>& v_mvaVariables,
                                                   std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    for(const MvaVariablesBase* mvaVariables : v_mvaVariables){
        
        const MvaVariablesEventClassification* mvaVariablesEventClassification = dynamic_cast<const MvaVariablesEventClassification*>(mvaVariables);
        if(!mvaVariablesEventClassification){
            std::cerr<<"ERROR in MvaTreePlotterEventClassification::fillHistos()! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
            exit(395);
        }
        
        const double weight(mvaVariablesEventClassification->eventWeight_.value_);
        
        double value(-999.);
        
        name = mvaVariablesEventClassification->multiplicity_jets_.name();
        value = mvaVariablesEventClassification->multiplicity_jets_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesEventClassification, weight);
        
        name = mvaVariablesEventClassification->btagDiscriminatorAverage_tagged_.name();
        value = mvaVariablesEventClassification->btagDiscriminatorAverage_tagged_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesEventClassification, weight);
        
        name = mvaVariablesEventClassification->btagDiscriminatorAverage_untagged_.name();
        value = mvaVariablesEventClassification->btagDiscriminatorAverage_untagged_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesEventClassification, weight);
        
        name = mvaVariablesEventClassification->minDeltaR_jet_jet_.name();
        value = mvaVariablesEventClassification->minDeltaR_jet_jet_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesEventClassification, weight);
        
        name = mvaVariablesEventClassification->ptSum_jets_leptons_.name();
        value = mvaVariablesEventClassification->ptSum_jets_leptons_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesEventClassification, weight);
    }
}



void MvaTreePlotterEventClassification::bookHistosInclExcl(std::map<TString, TH1*>& m_histogram, const TString& prefix, const TString& step,
                                                           const TString& name, const TString& title,
                                                           const int& nBinX, const double& xMin, const double& xMax)
{
    // The target truth is per sample, no need for separation plots
    m_histogram[name] = this->store(new TH1D(prefix+name+step, title, nBinX, xMin, xMax));
}



void MvaTreePlotterEventClassification::fillHistosInclExcl(std::map<TString, TH1*>& m_histogram, const TString& name,
                                                           const double& variable,
                                                           const MvaVariablesBase* mvaVariables, const double& weight)
{
    const MvaVariablesEventClassification* mvaVariablesEventClassification = dynamic_cast<const MvaVariablesEventClassification*>(mvaVariables);
    if(!mvaVariablesEventClassification){
        std::cerr<<"ERROR in MvaTreePlotterEventClassification::fillHistosInclExcl()! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
        exit(395);
    }
    
    // The target truth is per sample, no need for separation plots
    m_histogram[name]->Fill(variable, weight);
}















