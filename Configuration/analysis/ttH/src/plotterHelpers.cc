#include <iostream>
#include <cstdlib>

#include "plotterHelpers.h"






DrawMode::DrawMode DrawMode::convertDrawMode(const std::string& drawMode)
{
    if(drawMode == "stacked") return stacked;
    if(drawMode == "overlaid") return overlaid;
    if(drawMode == "scaledoverlaid") return scaledoverlaid;
    if(drawMode == "scaledoverlaidfixed") return scaledoverlaidfixed;
    std::cout<<"Warning! The following draw mode conversion is not implemented: "<<drawMode<<std::endl;
    return undefined;
}



std::string DrawMode::convertDrawMode(const DrawMode& drawMode)
{
    if(drawMode == stacked) return "stacked";
    if(drawMode == overlaid) return "overlaid";
    if(drawMode == scaledoverlaid) return "scaledoverlaid";
    if(drawMode == scaledoverlaidfixed) return "scaledoverlaidfixed";
    if(drawMode == undefined) return "";
    std::cout<<"Error! Draw mode conversion is not implemented, break...\n"<<std::endl;
    exit(97);
}



std::vector<DrawMode::DrawMode> DrawMode::convertDrawModes(const std::vector<std::string>& drawModes)
{
    std::vector<DrawMode> v_drawMode;
    for(auto drawMode : drawModes) v_drawMode.push_back(convertDrawMode(drawMode));
    return v_drawMode;
}



std::vector<std::string> DrawMode::convertDrawModes(const std::vector<DrawMode>& drawModes)
{
    std::vector<std::string> v_drawMode;
    for(auto drawMode : drawModes) v_drawMode.push_back(convertDrawMode(drawMode));
    return v_drawMode;
}







GlobalCorrection::GlobalCorrection GlobalCorrection::convertGlobalCorrection(const std::string& globalCorrection)
{
    if(globalCorrection == "dy") return dy;
    if(globalCorrection == "ttbb") return ttbb;
    std::cout<<"Warning! The following global correction conversion is not implemented: "<<globalCorrection<<std::endl;
    return undefined;
}



std::string GlobalCorrection::convertGlobalCorrection(const GlobalCorrection& globalCorrection)
{
    if(globalCorrection == dy) return "dy";
    if(globalCorrection == ttbb) return "ttbb";
    if(globalCorrection == undefined) return "";
    std::cout<<"Error! Global correction conversion is not implemented, break...\n"<<std::endl;
    exit(97);
}



std::vector<GlobalCorrection::GlobalCorrection> GlobalCorrection::convertGlobalCorrections(const std::vector<std::string>& globalCorrections)
{
    std::vector<GlobalCorrection> v_globalCorrection;
    for(auto globalCorrection : globalCorrections) v_globalCorrection.push_back(convertGlobalCorrection(globalCorrection));
    return v_globalCorrection;
}



std::vector<std::string> GlobalCorrection::convertGlobalCorrections(const std::vector<GlobalCorrection>& globalCorrections)
{
    std::vector<std::string> v_globalCorrection;
    for(auto globalCorrection : globalCorrections) v_globalCorrection.push_back(convertGlobalCorrection(globalCorrection));
    return v_globalCorrection;
}




