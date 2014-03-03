#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TSelectorList.h>
#include <TIterator.h>
#include <TObject.h>
#include <TH1.h>

#include "MvaTreePlotterBase.h"
#include "MvaVariablesBase.h"





MvaTreePlotterBase::MvaTreePlotterBase(const TString& prefix,
                                       const std::map<TString, std::vector<MvaVariablesBase*> >& m_stepMvaVariables,
                                       const bool separationPowerPlots):
prefix_(prefix),
selectorList_(0),
m_stepMvaVariables_(m_stepMvaVariables),
plotExclusively_(separationPowerPlots)
{
    if(!m_stepMvaVariables.size()){
        std::cerr<<"ERROR in constructor of MvaTreePlotterBase! No single step containing MVA variables defined\n...break\n"<<std::endl;
        exit(392);
    }
}



void MvaTreePlotterBase::clear()
{
    selectorList_ = 0;
    
    for(auto stepHistograms : m_stepHistograms_){
        stepHistograms.second.m_histogram_.clear();
    }
    m_stepHistograms_.clear();
}



void MvaTreePlotterBase::plotVariables(const std::string& f_savename)
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



void MvaTreePlotterBase::plotVariables(TSelectorList* output)
{
    std::cout<<"--- Beginning control plots for MVA variables\n";
    
    // Set pointer to output, so that histograms are owned by it
    selectorList_ = output;
    
    // Loop over steps and plot all histograms
    for(const auto& stepMvaVariables : m_stepMvaVariables_){
        const TString& step(stepMvaVariables.first);
        const std::vector<MvaVariablesBase*>& v_mvaVariables(stepMvaVariables.second);
        std::map<TString, TH1*>& m_histogram = m_stepHistograms_[step].m_histogram_;
        this->bookHistos(step, m_histogram);
        this->fillHistos(step, v_mvaVariables, m_histogram);
    }
    
    std::cout<<"=== Finishing control plots for MVA variables\n\n";
}



void MvaTreePlotterBase::bookHistos(const TString&, std::map<TString, TH1*>&)
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method bookHistos() in MvaTreeAnalyzerBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



void MvaTreePlotterBase::fillHistos(const TString&,
                                    const std::vector<MvaVariablesBase*>&,
                                    std::map<TString, TH1*>&)
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method fillHistos() in MvaTreeAnalyzerBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



void MvaTreePlotterBase::bookHistosInclExcl(std::map<TString, TH1*>&, const TString&, const TString&,
                                            const TString&, const TString&,
                                            const int&, const double&, const double&)
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method bookHistosInclExcl() in MvaTreeAnalyzerBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



void MvaTreePlotterBase::fillHistosInclExcl(std::map<TString, TH1*>&, const TString&,
                                            const double&,
                                            const MvaVariablesBase*, const double&)
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method fillHistosInclExcl() in MvaTreeAnalyzerBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}















