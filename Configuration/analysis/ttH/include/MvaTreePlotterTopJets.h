#ifndef MvaTreePlotterTopJets_h
#define MvaTreePlotterTopJets_h

#include <vector>
#include <map>

class TString;
class TH1;

#include "MvaTreePlotterBase.h"

class MvaVariablesBase;





/// Class for plotting input variables for MVA,
/// trying to identify the jets coming from (anti)b's from (anti)tops
class MvaTreePlotterTopJets : public MvaTreePlotterBase{
    
public:
    
    /// Constructor for selection steps
    MvaTreePlotterTopJets(const std::map<TString, std::vector<MvaVariablesBase*> >& m_stepMvaVariables,
                          const bool separationPowerPlots =false);
    
    /// Destructor
    ~MvaTreePlotterTopJets(){}
    
    
    
private:
    
    /// Book histograms for given selection step
    virtual void bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram);
    
    /// Fill histograms for given selection step
    virtual void fillHistos(const TString& step,
                            const std::vector<MvaVariablesBase*>& v_mvaVariables,
                            std::map<TString, TH1*>& m_histogram);
    
    /// Book 1-D histograms exclusively for correct, swapped and wrong combinations, and inclusively
    virtual void bookHistosInclExcl(std::map<TString, TH1*>& m_histogram, const TString& prefix, const TString& step,
                                    const TString& name, const TString& title,
                                    const int& nBinX, const double& xMin, const double& xMax);
    
    /// Fill 1-D histograms exclusively for correct, swapped and wrong combinations, and inclusively
    virtual void fillHistosInclExcl(std::map<TString, TH1*>& m_histogram, const TString& name,
                                    const double& variable,
                                    const MvaVariablesBase* mvaTopJetsVariables, const double& weight =1.);
};







#endif







