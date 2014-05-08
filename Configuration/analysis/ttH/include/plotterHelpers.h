#ifndef plotterHelpers_h
#define plotterHelpers_h

#include <vector>
#include <string>

class TString;




namespace DrawMode{
    
    /// Draw modes for Higgs samples for analysis
    enum DrawMode{stacked, overlaid, scaledoverlaid, scaledoverlaidfixed, undefined};
    
    
    
    /// All draw modes allowed for Higgs samples for analysis
    const std::vector<DrawMode> allowedDrawModes
        {stacked, overlaid, scaledoverlaid, scaledoverlaidfixed};
    
    
    
    /// Convert a DrawMode from string to enum
    DrawMode convert(const TString& drawMode);
    
    /// Convert a DrawMode from enum to string
    TString convert(const DrawMode& drawMode);
    
    /// Convert a vector of DrawModes from string to enum
    std::vector<DrawMode> convert(const std::vector<TString>& drawModes);
    
    /// Convert a vector of DrawModes from string to enum
    std::vector<DrawMode> convert(const std::vector<std::string>& drawModes);
    
    /// Convert a vector of DrawModes from enum to string
    std::vector<TString> convert(const std::vector<DrawMode>& drawModes);
}






namespace GlobalCorrection{
    
    /// Global corrections for analysis, i.e. scale factors which are applied to whole samples
    enum GlobalCorrection{dy, ttbb, undefined};
    
    /// All global corrections allowed for analysis
    const std::vector<GlobalCorrection> allowedGlobalCorrections
        {dy, ttbb};
    
    
    
    /// Convert a GlobalCorrection from string to enum
    GlobalCorrection convert(const TString& globalCorrection);
    
    /// Convert a GlobalCorrection from enum to string
    TString convert(const GlobalCorrection& globalCorrection);
    
    /// Convert a vector of GlobalCorrections from string to enum
    std::vector<GlobalCorrection> convert(const std::vector<TString>& globalCorrections);
    
    /// Convert a vector of GlobalCorrections from string to enum
    std::vector<GlobalCorrection> convert(const std::vector<std::string>& globalCorrections);
    
    /// Convert a vector of GlobalCorrections from enum to string
    std::vector<TString> convert(const std::vector<GlobalCorrection>& globalCorrections);
}






#endif








