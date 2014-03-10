#ifndef MvaTreePlotterBase_h
#define MvaTreePlotterBase_h

#include <vector>
#include <string>
#include <map>

#include <TString.h>

class TSelectorList;
class TH1;

#include "../../common/include/storeTemplate.h"

class MvaVariablesBase;





/// Base class for plotting input variables for MVA
class MvaTreePlotterBase{
    
public:
    
    /// Constructor for selection steps
    MvaTreePlotterBase(const TString& prefix,
                       const std::map<TString, std::vector<MvaVariablesBase*> >& m_stepMvaVariables,
                       const bool separationPowerPlots =false);
    
    /// Destructor
    ~MvaTreePlotterBase(){}
    
    
    
    /// Clear the class instance
    void clear();
    
    
    
    /// Plot the variables and write them to the specified folder
    /// If separationPowerPlots==true: plot them exclusively for signal- and background-like cases to see separation power
    /// If separationPowerPlots==false: plot them inclusively (especially for data/MC comparisons)
    void plotVariables(const std::string& f_savename);
    
    /// Plot the variables and store the histograms in the specified TSelectorList
    /// If separationPowerPlots==true: plot them exclusively for signal- and background-like cases to see separation power
    /// If separationPowerPlots==false: plot them inclusively (especially for data/MC comparisons)
    void plotVariables(TSelectorList* output);
    
    
    
protected:
    
    /// Struct holding the histograms for one selection step
    struct StepHistograms{
        /// Constructor
        StepHistograms(){};
        /// Destructor
        ~StepHistograms(){};

        /// The map with all the histograms for one selection step
        std::map<TString, TH1*> m_histogram_;
    };
    
    
    
    /// Book histograms for given selection step (dummy method, override in inherited MvaTreePlotter)
    virtual void bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram);
    
    /// Fill histograms for given selection step (dummy method, override in inherited MvaTreePlotter)
    virtual void fillHistos(const TString& step,
                            const std::vector<MvaVariablesBase*>& v_mvaVariables,
                            std::map<TString, TH1*>& m_histogram);
    
    /// Book 1-D histograms inclusively or exclusively in case of MVA variables which can be signal- and background-like within one sample
    /// (dummy method, override in inherited MvaTreePlotter)
    virtual void bookHistosInclExcl(std::map<TString, TH1*>& m_histogram, const TString& prefix, const TString& step,
                                    const TString& name, const TString& title,
                                    const int& nBinX, const double& xMin, const double& xMax);
    
    /// Fill 1-D histograms inclusively or exclusively in case of MVA variables which can be signal- and background-like within one sample
    /// (dummy method, override in inherited MvaTreePlotter)
    virtual void fillHistosInclExcl(std::map<TString, TH1*>& m_histogram, const TString& name,
                                    const double& variable,
                                    const MvaVariablesBase* mvaVariables, const double& weight);
    
    
    
    /// Store the object in the output list and return it
    template<class T> T* store(T* obj){return common::store(obj, selectorList_);}
    
    
    
    /// The prefix which all histograms of the specific tree plotter should have
    const TString prefix_;
    
    /// Pointer for bookkeeping of histograms, trees, etc.
    TSelectorList* selectorList_;
    
    /// The map containing all the vectors of MVA variables for all selection steps
    const std::map<TString, std::vector<MvaVariablesBase*> >& m_stepMvaVariables_;
    
    /// The map containing all the step histograms for all selection steps
    std::map<TString, StepHistograms> m_stepHistograms_;
    
    /// Whether to produce plots inclusively or exclusively in case of MVA variables which can be signal- and background-like within one sample
    const bool plotExclusively_;
};







#endif







