#ifndef MvaTreeAnalyzer_h
#define MvaTreeAnalyzer_h

#include <vector>
#include <string>
#include <map>

class TSelectorList;
class TString;
class TH1;

#include "MvaVariablesTopJets.h"
#include "../../common/include/storeTemplate.h"





/// Class for plotting input variables for MVA,
/// trying to identify the jets coming from (anti)b's from (anti)tops
class MvaTreeAnalyzer{
    
public:
    
    /// Constructor for selection steps
    MvaTreeAnalyzer(const std::map<TString, std::vector<MvaVariablesTopJets*> >& m_stepMvaVariables,
                    const bool separationPowerPlots =false);
    
    /// Destructor
    ~MvaTreeAnalyzer(){}
    
    
    
    /// Clear the class instance
    void clear();
    
    
    
    /// Plot the variables and write them to the specified folder
    /// If separationPowerPlots==true: plot them exclusively for the cases of correct, swapped and wrong combinations to see separation power
    /// If separationPowerPlots==false: plot them inclusively for data/MC comparisons
    void plotVariables(const std::string& f_savename);
    
    /// Plot the variables and store the histograms in the specified TSelectorList
    /// If separationPowerPlots==true: plot them exclusively for the cases of correct, swapped and wrong combinations to see separation power
    /// If separationPowerPlots==false: plot them inclusively for data/MC comparisons
    void plotVariables(TSelectorList* output);
    
    
    
private:
    
    /// Struct holding the histograms for one selection step
    struct StepHistograms{
        /// Constructor
        StepHistograms(){};
        /// Destructor
        ~StepHistograms(){};

        /// The map with all the histograms for one selection step
        std::map<TString, TH1*> m_histogram_;
    };
    
    
    
    /// Book and fill histograms for given selection step
    void plotStep(const TString& step, const std::vector<MvaVariablesTopJets*>& v_mvaTopJetsVariables);
    
    /// Book 1-D histograms exclusively for correct, swapped and wrong combinations, and inclusively
    void bookHistosInclExcl(std::map<TString, TH1*>& m_histogram, const TString& prefix, const TString& step,
                            const TString& name, const TString& title,
                            const int& nBinX, const double& xMin, const double& xMax);
    
    /// Fill 1-D histograms exclusively for correct, swapped and wrong combinations, and inclusively
    void fillHistosInclExcl(std::map<TString, TH1*>& m_histogram, const TString& name,
                            const double& variable,
                            const MvaVariablesTopJets* mvaTopJetsVariables, const double& weight =1.);
    
    
    
    /// Store the object in the output list and return it
    template<class T> T* store(T* obj){return common::store(obj, selectorList_);}
    
    
    
    /// Pointer for bookkeeping of histograms, trees, etc.
    TSelectorList* selectorList_;
    
    /// The map containing all the vectors of MVA variables for all selection steps
    const std::map<TString, std::vector<MvaVariablesTopJets*> >& m_stepMvaVariables_;
    
    /// The map containing all the step histograms for all selection steps
    std::map<TString, StepHistograms> m_stepHistograms_;
    
    /// Whether to produce plots inclusively, or exclusively for correct, swapped and wrong combinations
    const bool plotExclusively_;
};







#endif







