#ifndef MvaFactoryTopJets_h
#define MvaFactoryTopJets_h

class TFile;
class TTree;
class TString;
class TCut;
namespace TMVA{
    class Factory;
}

#include "../../common/include/sampleHelpers.h"
#include "MvaFactoryBase.h"

class MvaVariableInt;
class MvaVariableFloat;
namespace mvaSetup{
    class MvaSet;
}





/// Class for handling of TMVA factory,
/// trying to identify the jets coming from (anti)b's from (anti)tops
class MvaFactoryTopJets : public MvaFactoryBase{
    
public:
    
    /// Constructor
    MvaFactoryTopJets(const TString& mvaOutputDir, const TString& weightFileDir,
                      const TString& mergedTreesFileName);
    
    /// Destructor
    ~MvaFactoryTopJets(){};
    
    
    
private:
    
    /// Train MVAs for given MVA sets with trees from file
    virtual void trainMva(TFile* const treeFile, const mvaSetup::MvaSet& mvaSet, const TString& step)const;
    
    /// Configure the factory
    virtual void configureFactory(TMVA::Factory* const factory,
                                  const TCut& cutSignal, const TCut& cutBackground,
                                  TTree* const treeTraining, TTree* const treeTesting)const;
};







#endif







