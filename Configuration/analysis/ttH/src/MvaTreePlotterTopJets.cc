#include <iostream>

#include <TString.h>
#include <TH1.h>
#include <TH1D.h>

#include "MvaTreePlotterTopJets.h"
#include "MvaVariablesBase.h"
#include "MvaVariablesTopJets.h"





MvaTreePlotterTopJets::MvaTreePlotterTopJets(const std::map<TString, std::vector<MvaVariablesBase*> >& m_stepMvaVariables,
                                             const bool separationPowerPlots):
MvaTreePlotterBase("mvaTopP_", m_stepMvaVariables, separationPowerPlots)
{
    std::cout<<"--- Beginning setting up top jets MVA variables plotter\n";
    
    // Check that MVA variables are of correct type for plotter
    for(const auto& stepMvaVariables : m_stepMvaVariables){
        const std::vector<MvaVariablesBase*>& v_mvaVariables = stepMvaVariables.second;
        if(v_mvaVariables.size()){
            const MvaVariablesTopJets* const mvaVariablesTopJets = dynamic_cast<const MvaVariablesTopJets* const>(v_mvaVariables.at(0));
            if(!mvaVariablesTopJets){
                std::cerr<<"ERROR in constructor of MvaTreePlotterTopJets! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
                exit(395);
            }
            break;
        }
    }
    
    std::cout<<"=== Finishing setting up top jets MVA variables plotter\n\n";
}



void MvaTreePlotterTopJets::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    const MvaVariablesTopJets nameDummy;
    
    name = "trueStatus";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "True status of matched jets;Status;# jet pairs",2,0,2));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "swapped");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "correct");
    
    name = nameDummy.jetChargeDiff_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+";c_{rel}^{#bar{b}}-c_{rel}^{b};# jet pairs",50,0,2);
    
    name = nameDummy.meanDeltaPhi_b_met_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+";0.5(|#Delta#phi(b,MET)|+|#Delta#phi(#bar{b},MET)|)  [rad];# jet pairs",20,0,3.2);
    
    name = nameDummy.massDiff_recoil_bbbar_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+"; m_{recoil}^{jets}-m^{b#bar{b}}  [GeV];# jet pairs",16,-600,600);
    
    name = nameDummy.pt_b_antiLepton_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+"; p_{T}^{bl^{+}}  [GeV];# jet pairs",20,0,500);
    
    name = nameDummy.pt_antiB_lepton_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+"; p_{T}^{#bar{b}l^{-}}  [GeV];# jet pairs",20,0,500);
    
    name = nameDummy.deltaR_b_antiLepton_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+"; #DeltaR(b,l^{+});# jet pairs",25,0,5);
    
    name = nameDummy.deltaR_antiB_lepton_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+"; #DeltaR(#bar{b},l^{-});# jet pairs",25,0,5);
    
    name = nameDummy.btagDiscriminatorSum_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+"; d^{b}+d^{#bar{b}};# jet pairs",20,0,2);
    
    name = nameDummy.deltaPhi_antiBLepton_bAntiLepton_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+"; |#Delta#phi(bl^{+},#bar{b}l^{-})|  [rad];# jet pairs",10,0,3.2);
    
    name = nameDummy.massDiff_fullBLepton_bbbar_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+"; m^{b#bar{b}l^{+}l^{-}}-m^{b#bar{b}}  [GeV];# jet pairs",13,0,1050);
    
    name = nameDummy.meanMt_b_met_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+"; 0.5(m_{T}^{b,MET}+m_{T}^{#bar{b},MET)}  [GeV];# jet pairs",21,0,630);
    
    name = nameDummy.massSum_antiBLepton_bAntiLepton_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+"; m^{#bar{b}l^{-}}+m^{bl^{+}}  [GeV];# jet pairs",21,0,840);
    
    name = nameDummy.massDiff_antiBLepton_bAntiLepton_.name();
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, name+"; m^{#bar{b}l^{-}}-m^{bl^{+}}  [GeV];# jet pairs",41,-400,420);
}



void MvaTreePlotterTopJets::fillHistos(const TString&,
                                       const std::vector<MvaVariablesBase*>& v_mvaVariables,
                                       std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    for(const MvaVariablesBase* mvaVariables : v_mvaVariables){
        
        const MvaVariablesTopJets* mvaVariablesTopJets = dynamic_cast<const MvaVariablesTopJets*>(mvaVariables);
        if(!mvaVariablesTopJets){
            std::cerr<<"ERROR in MvaTreePlotterTopJets::fillHistos()! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
            exit(395);
        }
        
        const double weight(mvaVariablesTopJets->eventWeight_.value_);
        
        name = "trueStatus";
        if(mvaVariablesTopJets->swappedCombination_.value_) m_histogram[name]->Fill(0., weight);
        if(mvaVariablesTopJets->correctCombination_.value_) m_histogram[name]->Fill(1., weight);
        
        double value(-999.);
        
        name = mvaVariablesTopJets->jetChargeDiff_.name();
        value = mvaVariablesTopJets->jetChargeDiff_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesTopJets, weight);
        
        name = mvaVariablesTopJets->meanDeltaPhi_b_met_.name();
        value = mvaVariablesTopJets->meanDeltaPhi_b_met_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesTopJets, weight);
        
        name = mvaVariablesTopJets->massDiff_recoil_bbbar_.name();
        value = mvaVariablesTopJets->massDiff_recoil_bbbar_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesTopJets, weight);
        
        name = mvaVariablesTopJets->pt_b_antiLepton_.name();
        value = mvaVariablesTopJets->pt_b_antiLepton_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesTopJets, weight);
        
        name = mvaVariablesTopJets->pt_antiB_lepton_.name();
        value = mvaVariablesTopJets->pt_antiB_lepton_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesTopJets, weight);
        
        name = mvaVariablesTopJets->deltaR_b_antiLepton_.name();
        value = mvaVariablesTopJets->deltaR_b_antiLepton_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesTopJets, weight);
        
        name = mvaVariablesTopJets->deltaR_antiB_lepton_.name();
        value = mvaVariablesTopJets->deltaR_antiB_lepton_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesTopJets, weight);
        
        name = mvaVariablesTopJets->btagDiscriminatorSum_.name();
        value = mvaVariablesTopJets->btagDiscriminatorSum_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesTopJets, weight);
        
        name = mvaVariablesTopJets->deltaPhi_antiBLepton_bAntiLepton_.name();
        value = mvaVariablesTopJets->deltaPhi_antiBLepton_bAntiLepton_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesTopJets, weight);
        
        name = mvaVariablesTopJets->massDiff_fullBLepton_bbbar_.name();
        value = mvaVariablesTopJets->massDiff_fullBLepton_bbbar_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesTopJets, weight);
        
        name = mvaVariablesTopJets->meanMt_b_met_.name();
        value = mvaVariablesTopJets->meanMt_b_met_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesTopJets, weight);
        
        name = mvaVariablesTopJets->massSum_antiBLepton_bAntiLepton_.name();
        value = mvaVariablesTopJets->massSum_antiBLepton_bAntiLepton_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesTopJets, weight);
        
        name = mvaVariablesTopJets->massDiff_antiBLepton_bAntiLepton_.name();
        value = mvaVariablesTopJets->massDiff_antiBLepton_bAntiLepton_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaVariablesTopJets, weight);
    }
}



void MvaTreePlotterTopJets::bookHistosInclExcl(std::map<TString, TH1*>& m_histogram, const TString& prefix, const TString& step,
                                               const TString& name, const TString& title,
                                               const int& nBinX, const double& xMin, const double& xMax)
{
    const TString correct("correct_");
    const TString swapped("swapped_");
    const TString wrong("wrong_");
    
    if(!plotExclusively_){
        m_histogram[name] = this->store(new TH1D(prefix+name+step, title, nBinX, xMin, xMax));
    }
    else{
        m_histogram[correct+name] = this->store(new TH1D(correct+prefix+name+step, title, nBinX, xMin, xMax));
        m_histogram[swapped+name] = this->store(new TH1D(swapped+prefix+name+step, title, nBinX, xMin, xMax));
        m_histogram[wrong+name] = this->store(new TH1D(wrong+prefix+name+step, title, nBinX, xMin, xMax));
    }
}



void MvaTreePlotterTopJets::fillHistosInclExcl(std::map<TString, TH1*>& m_histogram, const TString& name,
                                               const double& variable,
                                               const MvaVariablesBase* mvaVariables, const double& weight)
{
    const MvaVariablesTopJets* mvaVariablesTopJets = dynamic_cast<const MvaVariablesTopJets*>(mvaVariables);
    if(!mvaVariablesTopJets){
        std::cerr<<"ERROR in MvaTreePlotterTopJets::fillHistosInclExcl()! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
        exit(395);
    }
    
    const TString correct("correct_");
    const TString swapped("swapped_");
    const TString wrong("wrong_");
    
    if(!plotExclusively_){
        m_histogram[name]->Fill(variable, weight);
    }
    else{
        if(mvaVariablesTopJets->correctCombination_.value_) m_histogram[correct+name]->Fill(variable, weight);
        else if(mvaVariablesTopJets->swappedCombination_.value_) m_histogram[swapped+name]->Fill(variable, weight);
        else m_histogram[wrong+name]->Fill(variable, weight);
    }
}















