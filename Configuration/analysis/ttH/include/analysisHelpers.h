#ifndef analysisHelpers_h
#define analysisHelpers_h

#include <string>
#include <vector>





namespace AnalysisMode{
    
    /// Analysis modes for analysis
    enum AnalysisMode{
        cp,         // Basic control plots analyser
        dijet,      // Dijet analyser
        charge,     // Jet charge analyser
        match,      // Jet matching analyser
        playg,      // Playground analyser
        weight,     // Event weight analyser
        gen,        // Generator event analyser
        mvaTopP,    // Produce MVA input for top system jet assignment
        mvaTopA,    // Apply MVA weights for top system jet assignment
        mvaEventP,  // Produce MVA input for event classification
        mvaEventA,  // Apply MVA weights for event classification
        mvaChargeP, // Produce MVA input for jet charge
        mvaChargeA, // Apply MVA weights for jet charge
        undefined   // Undefined mode
    };
    
    
    
    /// All analysis modes allowed for analysis
    const std::vector<AnalysisMode> allowedAnalysisModes
        {cp, dijet, charge, match, playg, weight, gen, mvaTopP, mvaTopA, mvaEventP, mvaEventA, mvaChargeP, mvaChargeA};
    
    
    
    /// Convert an AnalysisMode from string to typedef
    AnalysisMode convertAnalysisMode(const std::string& analysisMode);
    
    /// Convert an AnalysisMode from typedef to string
    std::string convertAnalysisMode(const AnalysisMode& analysisMode);
    
    /// Convert a vector of AnalysisModes from string to typedef
    std::vector<AnalysisMode> convertAnalysisModes(const std::vector<std::string>& analysisModes);
    
    /// Convert a vector of AnalysisModes from typedef to string
    std::vector<std::string> convertAnalysisModes(const std::vector<AnalysisMode>& analysisModes);
}










#endif // analysisHelpers_h








