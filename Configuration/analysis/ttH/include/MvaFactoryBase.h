#ifndef MvaFactoryBase_h
#define MvaFactoryBase_h

#include <vector>
#include <fstream>

#include <TString.h>

class TFile;
class TTree;
class TCut;
namespace TMVA{
    class Factory;
}

#include "SamplesFwd.h"

class MvaVariableInt;
class MvaVariableFloat;
namespace mvaSetup{
    class MvaSet;
    class MvaConfig;
}
class RootFileReader;




/// Base class for handling of TMVA factory
class MvaFactoryBase{
    
public:
    
    MvaFactoryBase(const TString& mvaOutputDir, const TString& weightFileDir,
                   const Samples& samples,
                   const bool inOneFactory =true);
    
    
    
    /// Constructor
    /// If inOneFactory true: run all defined training configurations in one factory (all correlations between trainings available)
    /// If inOneFactory false: run each training in own factory
    MvaFactoryBase(const TString& mvaOutputDir, const TString& weightFileDir,
                   const TString& treeFileName,
                   const bool inOneFactory =true);
    
    /// Destructor
    ~MvaFactoryBase(){};
    
    
    
    /// Train MVAs for given MVA sets
    void train(const std::vector<mvaSetup::MvaSet>& v_mvaSet)const;
    
    /// Train MVAs for given MVA sets
    void train2(const std::vector<mvaSetup::MvaSet>& v_mvaSet)const;
    
    
    
protected:
    
    /// Train MVAs for given MVA sets with trees from file (dummy method, override in inherited MvaFactory)
    virtual void trainMva(TFile* const treeFile, const mvaSetup::MvaSet& mvaSet, const TString& stepName)const;
    
    /// Train MVAs for given MVA sets with trees from file (dummy method, override in inherited MvaFactory)
    virtual void trainMva2(const mvaSetup::MvaSet& mvaSet, const TString& outputFolder,
                           const std::vector<Sample>& v_sample, const std::vector<double>& v_weight,
                           std::ofstream& trainingNameFile, const TString& stepName)const;
    
    /// Configure the factory (dummy method, override in inherited MvaFactory)
    virtual void configureFactory(TMVA::Factory* const factory,
                                  const TCut& cutSignal, const TCut& cutBackground,
                                  TTree* const treeTraining, TTree* const treeTesting)const;
    
    /// Configure the factory (dummy method, override in inherited MvaFactory)
    virtual void configureFactory2(TMVA::Factory* const factory,
                                  const TCut& cutSignal, const TCut& cutBackground,
                                  const std::vector<Sample>& v_sample, const std::vector<double>& v_weight,
                                  const TString& stepName)const;
    
    
    
    /// Run the MVA for given parameters
    void runMva(const char* const methodPrefix, const TCut& cutSignal, const TCut& cutBackground,
                TTree* const treeTraining, TTree* const treeTesting,
                const std::vector<mvaSetup::MvaConfig>& v_mvaSet,
                const TString& stepName)const;
    
    /// Run the MVA for given parameters
    void runMva2(const std::vector<mvaSetup::MvaConfig>& v_mvaSet, const TString& outputFolder,
                const std::vector<Sample>& v_sample, const std::vector<double>& v_weight,
                std::ofstream& trainingNameFile, const TString& stepName)const;
    
    /// Add a variable to the factory of type Int_t
    void addVariable(TMVA::Factory* const factory, const MvaVariableInt& variable)const;
    
    /// Add a variable to the factory of type Float_t
    void addVariable(TMVA::Factory* const factory, const MvaVariableFloat& variable)const;
    
    /// Add a spectator to the factory of type Int_t
    void addSpectator(TMVA::Factory* const factory, const MvaVariableInt& variable)const;
    
    /// Add a spectator to the factory of type Float_t
    void addSpectator(TMVA::Factory* const factory, const MvaVariableFloat& variable)const;
    
    /// Read TTree from sample
    TTree* readTree(const TString& filename, const TString& treename);
    
    
    
private:
    
    /// The folder where to store the input for MVA
    const TString mvaOutputDir_;
    
    /// The sub-folder where to store MVA weights
    const TString weightFileDir_;
    
    /// Samples to be analysed
    const Samples& samples_;
    
    /// The file containing the input trees for MVA training
    const TString treeFileName_;
    
    /// Boolean whether defined training configurations are run in one factory (all correlations between trainings available),
    /// or to run each training in own factory (however memory leak in TMVA leads to huge memory consumtion also in this mode)
    const bool inOneFactory_;
    
    
    
    /// File reader for accessing specific tree from given file
    RootFileReader* fileReader_;
};







#endif







