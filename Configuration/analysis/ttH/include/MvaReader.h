#ifndef MvaReader_h
#define MvaReader_h

#include <vector>

class TString;
class TH2;
namespace TMVA{
    class Reader;
}

#include "MvaVariablesTopJets.h"





/// Class for reading MVA weights from file,
/// based on set of MVA input variables per event
class MvaReader{
    
public:
    
    /// Constructor which sets MVA weights and creating TMVA Reader
    MvaReader(const char* mvaWeightsFile);
    
    /// Destructor
    ~MvaReader(){};
    
    
    
    /// Clear the class instance
    void clear();
    
    
    
    /// Get the MVA weights for all jet pairs from weights file
    std::vector<float> mvaWeights(const std::vector<MvaVariablesTopJets*>& v_mvaTopJetsVariables);
    
    
    
    /// Get the combined 1D weight from two MVA weights, and the 2D histogram of these weights from training
    static std::vector<float> combinedWeight(const TH2* weights2d,
                                             const std::vector<float>& weightCorrect,
                                             const std::vector<float>& weightSwapped);
    
    
    
private:
    
    /// Add a variable to MVA reader of type Int_t
    void addVariable(MvaVariableInt& variable);
    
    /// Add a variable to MVA reader of type Float_t
    void addVariable(MvaVariableFloat& variable);
    
    /// Add a spectator to MVA reader of type Int_t
    void addSpectator(MvaVariableInt& variable);
    
    /// Add a spectator to MVA reader of type Float_t
    void addSpectator(MvaVariableFloat& variable);
    
    
    
    /// Storage of all entries for the MVA
    std::vector<MvaVariablesTopJetsPerEvent> v_mvaTopJetsVariablesPerEvent_;
    
    /// Pointer to TMVA Reader, i.e. tool for reading in MVA weights
    TMVA::Reader* mvaWeightsReader_;
    
    /// Struct for setting addresses of variables for mvaWeightsReader_
    MvaVariablesTopJets* mvaVariablesTopJets_;
};







#endif







