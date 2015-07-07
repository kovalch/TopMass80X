#ifndef analysisHelpers_h
#define analysisHelpers_h

#include <vector>
#include <string>

class TString;




namespace AnalysisMode{
    
    /// Analysis modes for analysis
    enum AnalysisMode{
        cp,         // Basic control plots analyser
        dda,        // Double differential analyser
        btop,       // Boosted top analyser
        ddaTree,    // Produce Tree for double differential analysis
        btopTree,   // Produce Tree for boosted top analysis
        kinReco,    // kin. Reco. efficiency analyser
        dijet,      // Dijet analyser
        charge,     // Jet charge analyser
        match,      // Jet matching analyser
        playg,      // Playground analyser
        weight,     // Event weight analyser
        genEvent,   // Generator event analyser
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
        {cp, dda, btop, ddaTree, btopTree, kinReco, dijet, charge, match, playg, weight, genEvent, mvaTopP, mvaTopA, mvaEventP, mvaEventA, mvaChargeP, mvaChargeA};
    
    
    
    /// Convert an AnalysisMode from string to enum
    AnalysisMode convert(const TString& analysisMode);
    
    /// Convert an AnalysisMode from enum to string
    TString convert(const AnalysisMode& analysisMode);
    
    /// Convert a vector of AnalysisModes from string to enum
    std::vector<AnalysisMode> convert(const std::vector<TString>& analysisModes);
    
    /// Convert a vector of AnalysisModes from string to enum
    std::vector<AnalysisMode> convert(const std::vector<std::string>& analysisModes);
    
    /// Convert a vector of AnalysisModes from enum to string
    std::vector<TString> convert(const std::vector<AnalysisMode>& analysisModes);
}










#endif








