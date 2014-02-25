#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TSelectorList.h>
#include <TIterator.h>
#include <TObject.h>
#include <TH1.h>
#include <TH1D.h>

#include "MvaTreeAnalyzer.h"
#include "MvaVariablesTopJets.h"





MvaTreeAnalyzer::MvaTreeAnalyzer(const std::map<TString, std::vector<MvaVariablesTopJets> >& m_stepMvaVariables,
                                 const bool separationPowerPlots):
selectorList_(0),
m_stepMvaVariables_(m_stepMvaVariables),
plotExclusively_(separationPowerPlots)
{
    std::cout<<"--- Beginning setting up MVA variables histograms\n";
    std::cout<<"=== Finishing setting up MVA variables histograms\n\n";
}



void MvaTreeAnalyzer::clear()
{
    selectorList_ = 0;
    
    for(auto stepHistograms : m_stepHistograms_){
        stepHistograms.second.m_histogram_.clear();
    }
    m_stepHistograms_.clear();
}



void MvaTreeAnalyzer::plotVariables(const std::string& f_savename)
{
    // Output file
    TFile outputFile(f_savename.c_str(),"RECREATE");
    std::cout<<"\nOutput file for MVA input control plots: "<<f_savename<<"\n";
    
    // Produce MVA input control plots and store them in output
    TSelectorList* output = new TSelectorList();
    this->plotVariables(output);
    
    // Write file and cleanup
    TIterator* it = output->MakeIterator();
    while(TObject* obj = it->Next()){
        obj->Write();
    }
    outputFile.Close();
    output->SetOwner();
    output->Clear();
}



void MvaTreeAnalyzer::plotVariables(TSelectorList* output)
{
    std::cout<<"--- Beginning control plots for MVA variables\n";
    
    // Set pointer to output, so that histograms are owned by it
    selectorList_ = output;
    
    // Loop over steps and plot all histograms
    for(const auto& stepMvaVariables : m_stepMvaVariables_){
        const TString& step(stepMvaVariables.first);
        const std::vector<MvaVariablesTopJets>& v_mvaTopJetsVariables(stepMvaVariables.second);
        this->plotStep(step, v_mvaTopJetsVariables);
    }
    
    std::cout<<"=== Finishing control plots for MVA variables\n\n";
}



void MvaTreeAnalyzer::plotStep(const TString& step, const std::vector<MvaVariablesTopJets>& v_mvaTopJetsVariables)
{    
    const MvaVariablesTopJets nameDummy;
    constexpr const char* prefix = "mvaP_";
    std::map<TString, TH1*>& m_histogram = m_stepHistograms_[step].m_histogram_;
    TString name;
    
    // Book histograms
    name = "trueStatus";
    m_histogram[name] = store(new TH1D(prefix+name+step, "True status of matched jets;Status;# jet pairs",2,0,2));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "swapped");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "correct");
    
    name = nameDummy.jetChargeDiff_.name();
    this->bookHistosInclExcl(m_histogram, prefix, step, name, name+";c_{rel}^{#bar{b}}-c_{rel}^{b};# jet pairs",50,0,2);
    
    name = nameDummy.meanDeltaPhi_b_met_.name();
    this->bookHistosInclExcl(m_histogram, prefix, step, name, name+";0.5(|#Delta#phi(b,MET)|+|#Delta#phi(#bar{b},MET)|)  [rad];# jet pairs",20,0,3.2);
    
    name = nameDummy.massDiff_recoil_bbbar_.name();
    this->bookHistosInclExcl(m_histogram, prefix, step, name, name+"; m_{recoil}^{jets}-m^{b#bar{b}}  [GeV];# jet pairs",16,-600,600);
    
    name = nameDummy.pt_b_antiLepton_.name();
    this->bookHistosInclExcl(m_histogram, prefix, step, name, name+"; p_{T}^{bl^{+}}  [GeV];# jet pairs",20,0,500);
    
    name = nameDummy.pt_antiB_lepton_.name();
    this->bookHistosInclExcl(m_histogram, prefix, step, name, name+"; p_{T}^{#bar{b}l^{-}}  [GeV];# jet pairs",20,0,500);
    
    name = nameDummy.deltaR_b_antiLepton_.name();
    this->bookHistosInclExcl(m_histogram, prefix, step, name, name+"; #DeltaR(b,l^{+});# jet pairs",25,0,5);
    
    name = nameDummy.deltaR_antiB_lepton_.name();
    this->bookHistosInclExcl(m_histogram, prefix, step, name, name+"; #DeltaR(#bar{b},l^{-});# jet pairs",25,0,5);
    
    name = nameDummy.btagDiscriminatorSum_.name();
    this->bookHistosInclExcl(m_histogram, prefix, step, name, name+"; d^{b}+d^{#bar{b}};# jet pairs",20,0,2);
    
    name = nameDummy.deltaPhi_antiBLepton_bAntiLepton_.name();
    this->bookHistosInclExcl(m_histogram, prefix, step, name, name+"; |#Delta#phi(bl^{+},#bar{b}l^{-})|  [rad];# jet pairs",10,0,3.2);
    
    name = nameDummy.massDiff_fullBLepton_bbbar_.name();
    this->bookHistosInclExcl(m_histogram, prefix, step, name, name+"; m^{b#bar{b}l^{+}l^{-}}-m^{b#bar{b}}  [GeV];# jet pairs",13,0,1050);
    
    name = nameDummy.meanMt_b_met_.name();
    this->bookHistosInclExcl(m_histogram, prefix, step, name, name+"; 0.5(m_{T}^{b,MET}+m_{T}^{#bar{b},MET)}  [GeV];# jet pairs",21,0,630);
    
    name = nameDummy.massSum_antiBLepton_bAntiLepton_.name();
    this->bookHistosInclExcl(m_histogram, prefix, step, name, name+"; m^{#bar{b}l^{-}}+m^{bl^{+}}  [GeV];# jet pairs",21,0,840);
    
    name = nameDummy.massDiff_antiBLepton_bAntiLepton_.name();
    this->bookHistosInclExcl(m_histogram, prefix, step, name, name+"; m^{#bar{b}l^{-}}-m^{bl^{+}}  [GeV];# jet pairs",41,-400,420);
    
    
    
    // Fill histograms
    for(const MvaVariablesTopJets& mvaTopJetsVariables : v_mvaTopJetsVariables){
        
        const double weight(mvaTopJetsVariables.eventWeight_.value_);
        
        name = "trueStatus";
        if(mvaTopJetsVariables.swappedCombination_.value_) m_histogram[name]->Fill(0., weight);
        if(mvaTopJetsVariables.correctCombination_.value_) m_histogram[name]->Fill(1., weight);
        
        double value;
        
        name = mvaTopJetsVariables.jetChargeDiff_.name();
        value = mvaTopJetsVariables.jetChargeDiff_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaTopJetsVariables, weight);
        
        name = mvaTopJetsVariables.meanDeltaPhi_b_met_.name();
        value = mvaTopJetsVariables.meanDeltaPhi_b_met_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaTopJetsVariables, weight);
        
        name = mvaTopJetsVariables.massDiff_recoil_bbbar_.name();
        value = mvaTopJetsVariables.massDiff_recoil_bbbar_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaTopJetsVariables, weight);
        
        name = mvaTopJetsVariables.pt_b_antiLepton_.name();
        value = mvaTopJetsVariables.pt_b_antiLepton_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaTopJetsVariables, weight);
        
        name = mvaTopJetsVariables.pt_antiB_lepton_.name();
        value = mvaTopJetsVariables.pt_antiB_lepton_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaTopJetsVariables, weight);
        
        name = mvaTopJetsVariables.deltaR_b_antiLepton_.name();
        value = mvaTopJetsVariables.deltaR_b_antiLepton_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaTopJetsVariables, weight);
        
        name = mvaTopJetsVariables.deltaR_antiB_lepton_.name();
        value = mvaTopJetsVariables.deltaR_antiB_lepton_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaTopJetsVariables, weight);
        
        name = mvaTopJetsVariables.btagDiscriminatorSum_.name();
        value = mvaTopJetsVariables.btagDiscriminatorSum_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaTopJetsVariables, weight);
        
        name = mvaTopJetsVariables.deltaPhi_antiBLepton_bAntiLepton_.name();
        value = mvaTopJetsVariables.deltaPhi_antiBLepton_bAntiLepton_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaTopJetsVariables, weight);
        
        name = mvaTopJetsVariables.massDiff_fullBLepton_bbbar_.name();
        value = mvaTopJetsVariables.massDiff_fullBLepton_bbbar_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaTopJetsVariables, weight);
        
        name = mvaTopJetsVariables.meanMt_b_met_.name();
        value = mvaTopJetsVariables.meanMt_b_met_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaTopJetsVariables, weight);
        
        name = mvaTopJetsVariables.massSum_antiBLepton_bAntiLepton_.name();
        value = mvaTopJetsVariables.massSum_antiBLepton_bAntiLepton_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaTopJetsVariables, weight);
        
        name = mvaTopJetsVariables.massDiff_antiBLepton_bAntiLepton_.name();
        value = mvaTopJetsVariables.massDiff_antiBLepton_bAntiLepton_.value_;
        this->fillHistosInclExcl(m_histogram, name, value, mvaTopJetsVariables, weight);
    }
}



void MvaTreeAnalyzer::bookHistosInclExcl(std::map<TString, TH1*>& m_histogram, const TString& prefix, const TString& step,
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



void MvaTreeAnalyzer::fillHistosInclExcl(std::map<TString, TH1*>& m_histogram, const TString& name,
                                        const double& variable,
                                        const MvaVariablesTopJets& mvaTopJetsVariables, const double& weight)
{
    const TString correct("correct_");
    const TString swapped("swapped_");
    const TString wrong("wrong_");
    
    if(!plotExclusively_){
        m_histogram[name]->Fill(variable, weight);
    }
    else{
        if(mvaTopJetsVariables.correctCombination_.value_) m_histogram[correct+name]->Fill(variable, weight);
        else if(mvaTopJetsVariables.swappedCombination_.value_) m_histogram[swapped+name]->Fill(variable, weight);
        else m_histogram[wrong+name]->Fill(variable, weight);
    }
}















