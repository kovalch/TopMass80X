#ifndef MvaReaderBase_h
#define MvaReaderBase_h

#include <vector>

#include <TString.h>

class TH2;
namespace TMVA{
    class Reader;
}

class MvaVariableInt;
class MvaVariableFloat;
class MvaVariablesBase;





/// Base class for reading MVA weights from file,
/// based on set of MVA input variables per event
class MvaReaderBase{
    
public:
    
    /// Constructor which sets MVA weights and creating TMVA Reader
    MvaReaderBase(const TString& mvaMethod);
    
    /// Destructor
    ~MvaReaderBase(){};
    
    
    
    /// Clear the class instance
    void clear();
    
    
    
    /// Book the MVA reader with corresponding variables
    void book(const TString& mvaWeightsFilename);
    
    /// Get the MVA weight for one set of variables from weights file
    float mvaWeight(const MvaVariablesBase* const mvaVariables)const;
    
    /// Get the MVA weights for several sets of variables per event from weights file
    std::vector<float> mvaWeights(const std::vector<MvaVariablesBase*>& v_mvaVariables)const;
    
    /// Get the combined 1D weights from two MVA weights, and the 2D histogram of these weights from training,
    /// for one set of variables
    static float combinedWeight(const TH2* const weights2d, const float weight1, const float weight2);
    
    /// Get the combined 1D weights from two MVA weights, and the 2D histogram of these weights from training,
    /// for several sets of variables per event
    static std::vector<float> combinedWeights(const TH2* const weights2d,
                                              const std::vector<float>& v_weight1,
                                              const std::vector<float>& v_weight2);
    
    
    
protected:
    
    /// Set the MVA input variables for the MVA reader (dummy method, override in inherited MvaReader)
    virtual void bookVariables(TMVA::Reader* const mvaWeightsReader, MvaVariablesBase*& mvaVariables)const;
    
    /// Read the MVA input variables into the pointer associated to the MVA reader (dummy method, override in inherited MvaReader)
    virtual void readVariables(MvaVariablesBase* const mvaVariables, const MvaVariablesBase* const mvaVariablesTmp)const;
    
    
    
    /// Add a variable to MVA reader of type Int_t
    void addVariable(TMVA::Reader* const mvaWeightsReader, MvaVariableInt& variable)const;
    
    /// Add a variable to MVA reader of type Float_t
    void addVariable(TMVA::Reader* const mvaWeightsReader, MvaVariableFloat& variable)const;
    
    /// Add a spectator to MVA reader of type Int_t
    void addSpectator(TMVA::Reader* const mvaWeightsReader, MvaVariableInt& variable)const;
    
    /// Add a spectator to MVA reader of type Float_t
    void addSpectator(TMVA::Reader* const mvaWeightsReader, MvaVariableFloat& variable)const;
    
    
    
private:
    
    /// Which MVA method
    const TString mvaMethod_;
    
    /// Pointer to TMVA Reader, i.e. tool for reading in MVA weights
    TMVA::Reader* mvaWeightsReader_;
    
    /// Struct for setting addresses of variables for mvaWeightsReader_
    MvaVariablesBase* mvaVariables_;
};







#endif







