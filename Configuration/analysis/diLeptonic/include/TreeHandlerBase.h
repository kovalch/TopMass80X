#ifndef TreeHandlerBase_h
#define TreeHandlerBase_h

#include <vector>
#include <map>
#include <utility>

#include <TString.h>

class TTree;
class TSelectorList;

#include "../../common/include/storeTemplate.h"
#include "../../common/include/sampleHelpers.h"

class VariableInt;
class VariableFloat;
class VariablesBase;
class RecoObjects;
class CommonGenObjects;
class TopGenObjects;
class KinRecoObjects;
namespace ttbar{
    class RecoLevelWeights;
    class GenLevelWeights;
    class GenObjectIndices;
    class RecoObjectIndices;
}





/// Base class for handling the trees of input variables
class TreeHandlerBase{
    
public:
    
    /// Constructor for selection steps
    TreeHandlerBase(const TString& prefix,
                       const char* inputDir,
                       const std::vector<TString>& selectionStepsNoCategories);
    
    /// Destructor
    ~TreeHandlerBase(){}
    
    
    
    /// Book the vector of tree variables for each requested selection step
    void book();
    
    void fill(const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
              const TopGenObjects& topGenObjects,
              const KinRecoObjects& kinRecoObjects,
              const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
              const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
              const double& weight, const TString& stepShort);
    
    /// Write trees with MVA input variables in own file
    void writeTrees(const TString& outputFilename,
                    const Channel::Channel& channel, const Systematic::Systematic& systematic);
    
    /// Write trees with MVA input variables owned by a given selectorList
    void writeTrees(TSelectorList* output);
    
    /// Return a constant reference to the map containing all the vectors of MVA variables for all selection steps
    const std::map<TString, std::vector<VariablesBase*> >& stepVariablesMap()const;
    
    /// Import written TTrees
    void importTrees(const TString& f_savename, const TString& prefix ="");
    
    
    
    
    
    /// Clear the class instance
    void clear();
    
    
    
protected:
    
    /// Fill all variables for given selection step (dummy method, override in inherited MvaTreeHandler)
    virtual void fillVariables(const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                               const TopGenObjects& topGenObjects,
                               const KinRecoObjects& kinRecoObjects,
                               const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                               const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                               const double& weight, const TString& step,
                               std::vector<VariablesBase*>& Variables)const;
    
    /// Import all branches from TTree (dummy method, override in inherited MvaTreeHandler)
    virtual void importBranches(TTree* tree, std::vector<VariablesBase*>& Variables)const;
    
    /// Create and fill branches for TTree holding the input variables for MVA (dummy method, override in inherited MvaTreeHandler)
    virtual void createAndFillBranches(TTree* tree, const std::vector<VariablesBase*>& v_Variables)const;
    
    
    
    /// Create single branch for TTree based on MvaInputVariables of type Int_t
    void createBranch(TTree* tree, const VariableInt& variable)const;
    
    /// Create single branch for TTree based on MvaInputVariables of type Float_t
    void createBranch(TTree* tree, const VariableFloat& variable)const;
    
    /// Import branch of type Int_t
    void importBranch(TTree* tree, VariableInt& variable)const;
    
    /// Import branch of type Float_t
    void importBranch(TTree* tree, VariableFloat& variable)const;
    
    
    
private:
    
    /// Store the object in the output list and return it
    template<class T> T* store(T* obj){return common::store(obj, selectorList_);}
    
    
    
    /// Add a new selection step
    void addStep(const TString& step);
    
    /// Check if selection step exists
    bool checkExistence(const TString& step)const;
    
    
    
    /// The prefix which all trees of the specific handler should have
    const TString prefix_;
    
    /// Pointer for bookkeeping of histograms, trees, etc.
    TSelectorList* selectorList_;
    
    /// The map containing all the vectors of MVA variables for all selection steps
    std::map<TString, std::vector<VariablesBase*> > m_stepVariables_;
    
    /// Selection steps for which to run the MVA tool
    const std::vector<TString> selectionSteps_;
    
    /// The folder where to store the input for MVA
    const char* inputDir_;
};






#endif






