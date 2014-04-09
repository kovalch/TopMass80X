#ifndef plotterHelpers_h
#define plotterHelpers_h

#include <string>
#include <vector>





namespace DrawMode{
    
    /// Draw modes for Higgs samples for analysis
    enum DrawMode{stacked, overlaid, scaledoverlaid, scaledoverlaidfixed, undefined};
    
    
    
    /// All draw modes allowed for Higgs samples for analysis
    const std::vector<DrawMode> allowedDrawModes
        {stacked, overlaid, scaledoverlaid, scaledoverlaidfixed};
    
    
    
    /// Convert a DrawMode from string to typedef
    DrawMode convertDrawMode(const std::string& drawMode);
    
    /// Convert a DrawMode from typedef to string
    std::string convertDrawMode(const DrawMode& drawMode);
    
    /// Convert a vector of DrawModes from string to typedef
    std::vector<DrawMode> convertDrawModes(const std::vector<std::string>& drawModes);
    
    /// Convert a vector of DrawModes from typedef to string
    std::vector<std::string> convertDrawModes(const std::vector<DrawMode>& drawModes);
}






namespace GlobalCorrection{
    
    /// Global corrections for analysis, i.e. scale factors which are applied to whole samples
    enum GlobalCorrection{dy, ttbb, undefined};
    
    /// All global corrections allowed for analysis
    const std::vector<GlobalCorrection> allowedGlobalCorrections
        {dy, ttbb};
    
    
    
    /// Convert a GlobalCorrection from string to typedef
    GlobalCorrection convertGlobalCorrection(const std::string& globalCorrection);
    
    /// Convert a GlobalCorrection from typedef to string
    std::string convertGlobalCorrection(const GlobalCorrection& globalCorrection);
    
    /// Convert a vector of GlobalCorrections from string to typedef
    std::vector<GlobalCorrection> convertGlobalCorrections(const std::vector<std::string>& globalCorrections);
    
    /// Convert a vector of GlobalCorrections from typedef to string
    std::vector<std::string> convertGlobalCorrections(const std::vector<GlobalCorrection>& globalCorrections);
}






#endif








