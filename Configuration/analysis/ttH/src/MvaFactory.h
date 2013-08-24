#ifndef MvaFactory_h
#define MvaFactory_h

#include <vector>
#include <string>

class TTree;
class TString;
class TCut;
namespace TMVA{
    class Factory;
}

class MvaVariableInt;
class MvaVariableFloat;





/// Class holding the input variables for MVA,
/// trying to identify the jets coming from (anti)b's from (anti)tops
class MvaFactory{
    
public:
    
    /// Constructor which can optionally set MVA weights and creating TMVA Reader
    MvaFactory(const char* mvaOutputDir, const char* weightFileDir,
               const std::vector<TString>& selectionSteps);
    
    /// Destructor
    ~MvaFactory(){};
    
    
    
    
    /// Clear the class instance
    void clear();
    
    
    
    /// Run the MVA for given parameters
    void runMva(const char* outputFileName,
                const char* methodName, const TCut& cutSignal, const TCut& cutBackground,
                TTree* treeTraining, TTree* treeTesting);
    
    
    
private:
    
    /// Add a variable to the factory of type Int_t
    void addVariable(TMVA::Factory* factory, MvaVariableInt& variable);
    
    /// Add a variable to the factory of type Float_t
    void addVariable(TMVA::Factory* factory, MvaVariableFloat& variable);
    
    /// Add a spectator to the factory of type Int_t
    void addSpectator(TMVA::Factory* factory, MvaVariableInt& variable);
    
    /// Add a spectator to the factory of type Float_t
    void addSpectator(TMVA::Factory* factory, MvaVariableFloat& variable);
    
    
    
    /// Selection steps for which to run the MVA tool
    const std::vector<TString> selectionSteps_;
    
    /// The folder where to store the input for MVA
    const char* mvaOutputDir_;
    
    /// The sub-folder where to store MVA weights
    const char* weightFileDir_;
};







#endif







