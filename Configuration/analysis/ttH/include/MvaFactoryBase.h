#ifndef MvaFactoryBase_h
#define MvaFactoryBase_h

#include <vector>

#include <TString.h>

class TFile;
class TTree;
class TCut;
namespace TMVA{
    class Factory;
}

class MvaVariableInt;
class MvaVariableFloat;
namespace mvaSetup{
    class MvaSet;
    class MvaConfig;
}





/// Base class for handling of TMVA factory
class MvaFactoryBase{
    
public:
    
    /// Constructor
    MvaFactoryBase(const TString& mvaOutputDir, const TString& weightFileDir,
                   const TString& treeFileName);
    
    /// Destructor
    ~MvaFactoryBase(){};
    
    
    
    /// Train MVAs for given MVA sets
    void train(const std::vector<mvaSetup::MvaSet>& v_mvaSet)const;
    
    
    
protected:
    
    /// Train MVAs for given MVA sets with trees from file (dummy method, override in inherited MvaFactory)
    virtual void trainMva(TFile* const treeFile, const mvaSetup::MvaSet& mvaSet, const TString& stepName)const;
    
    /// Configure the factory (dummy method, override in inherited MvaFactory)
    virtual void configureFactory(TMVA::Factory* const factory,
                                  const TCut& cutSignal, const TCut& cutBackground,
                                  TTree* const treeTraining, TTree* const treeTesting)const;
    
    
    
    /// Run the MVA for given parameters
    void runMva(const char* const methodPrefix, const TCut& cutSignal, const TCut& cutBackground,
                TTree* const treeTraining, TTree* const treeTesting,
                const std::vector<mvaSetup::MvaConfig>& v_mvaSet,
                const TString& stepName)const;
    
    /// Add a variable to the factory of type Int_t
    void addVariable(TMVA::Factory* const factory, const MvaVariableInt& variable)const;
    
    /// Add a variable to the factory of type Float_t
    void addVariable(TMVA::Factory* const factory, const MvaVariableFloat& variable)const;
    
    /// Add a spectator to the factory of type Int_t
    void addSpectator(TMVA::Factory* const factory, const MvaVariableInt& variable)const;
    
    /// Add a spectator to the factory of type Float_t
    void addSpectator(TMVA::Factory* const factory, const MvaVariableFloat& variable)const;
    
    
    
private:
    
    /// The folder where to store the input for MVA
    const TString mvaOutputDir_;
    
    /// The sub-folder where to store MVA weights
    const TString weightFileDir_;
    
    /// The file containing the input trees for MVA training
    const TString treeFileName_;
};







#endif







