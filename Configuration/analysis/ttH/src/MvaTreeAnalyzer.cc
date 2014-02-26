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





MvaTreeAnalyzer::MvaTreeAnalyzer(const std::map<TString, std::vector<MvaVariablesBase*> >& m_stepMvaVariables,
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
        const std::vector<MvaVariablesBase*>& v_mvaVariables(stepMvaVariables.second);
        this->plotStep(step, v_mvaVariables);
    }
    
    std::cout<<"=== Finishing control plots for MVA variables\n\n";
}



void MvaTreeAnalyzer::plotStep(const TString& step, const std::vector<MvaVariablesBase*>& v_mvaVariables)
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
    for(const MvaVariablesBase* mvaVariables : v_mvaVariables){
        
        const MvaVariablesTopJets* mvaVariablesTopJets = dynamic_cast<const MvaVariablesTopJets*>(mvaVariables);
        if(!mvaVariablesTopJets){
            std::cerr<<"ERROR in MvaTreeAnalyzer::plotStep()! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
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
                                        const MvaVariablesTopJets* mvaTopJetsVariables, const double& weight)
{
    const TString correct("correct_");
    const TString swapped("swapped_");
    const TString wrong("wrong_");
    
    if(!plotExclusively_){
        m_histogram[name]->Fill(variable, weight);
    }
    else{
        if(mvaTopJetsVariables->correctCombination_.value_) m_histogram[correct+name]->Fill(variable, weight);
        else if(mvaTopJetsVariables->swappedCombination_.value_) m_histogram[swapped+name]->Fill(variable, weight);
        else m_histogram[wrong+name]->Fill(variable, weight);
    }
}















