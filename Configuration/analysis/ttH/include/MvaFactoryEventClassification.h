#ifndef MvaFactoryEventClassification_h
#define MvaFactoryEventClassification_h

class TFile;
class TTree;
class TString;
class TCut;
namespace TMVA{
    class Factory;
}

#include "MvaFactoryBase.h"
#include "SamplesFwd.h"
#include "../../common/include/sampleHelpers.h"

class MvaVariableInt;
class MvaVariableFloat;
namespace mvaSetup{
    class MvaSet;
}





/// Class for handling of TMVA factory,
/// classifying events for signal-likeliness
class MvaFactoryEventClassification : public MvaFactoryBase{
    
public:
    
    /// Constructor
    MvaFactoryEventClassification(const TString& mvaOutputDir, const TString& weightFileDir,
                                  const Samples& samples,
                                  const bool inOneFactory =true);
    
    /// Destructor
    ~MvaFactoryEventClassification(){};
    
    
    
private:
    
    /// Train MVAs for given MVA sets with trees from file
    virtual void trainMva2(const mvaSetup::MvaSet& mvaSet, const TString& outputFolder,
                           const std::vector<Sample>& v_sample, const std::vector<double>& v_weight,
                           std::ofstream& trainingNameFile, const TString& step)const;
    
    /// Configure the factory
    virtual void configureFactory2(TMVA::Factory* const factory,
                                  const TCut& cutSignal, const TCut& cutBackground,
                                  const std::vector<Sample>& v_sample, const std::vector<double>& v_weight,
                                  const TString& stepName)const;
};







#endif







