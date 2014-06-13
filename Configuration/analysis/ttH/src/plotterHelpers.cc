#include <iostream>
#include <cstdlib>

#include <TString.h>

#include "plotterHelpers.h"






DrawMode::DrawMode DrawMode::convert(const TString& drawMode)
{
    if(drawMode == "stacked") return stacked;
    if(drawMode == "overlaid") return overlaid;
    if(drawMode == "scaledoverlaid") return scaledoverlaid;
    if(drawMode == "scaledoverlaidfixed") return scaledoverlaidfixed;
    std::cout<<"Warning! The following draw mode conversion is not implemented: "<<drawMode<<std::endl;
    return undefined;
}



TString DrawMode::convert(const DrawMode& drawMode)
{
    if(drawMode == stacked) return "stacked";
    if(drawMode == overlaid) return "overlaid";
    if(drawMode == scaledoverlaid) return "scaledoverlaid";
    if(drawMode == scaledoverlaidfixed) return "scaledoverlaidfixed";
    if(drawMode == undefined) return "";
    std::cout<<"Error! Draw mode conversion is not implemented, break...\n"<<std::endl;
    exit(97);
}



std::vector<DrawMode::DrawMode> DrawMode::convert(const std::vector<TString>& drawModes)
{
    std::vector<DrawMode> v_drawMode;
    for(auto drawMode : drawModes) v_drawMode.push_back(convert(drawMode));
    return v_drawMode;
}



std::vector<DrawMode::DrawMode> DrawMode::convert(const std::vector<std::string>& drawModes)
{
    std::vector<DrawMode> v_drawMode;
    for(auto drawMode : drawModes) v_drawMode.push_back(convert(drawMode));
    return v_drawMode;
}



std::vector<TString> DrawMode::convert(const std::vector<DrawMode>& drawModes)
{
    std::vector<TString> v_drawMode;
    for(auto drawMode : drawModes) v_drawMode.push_back(convert(drawMode));
    return v_drawMode;
}







GlobalCorrection::GlobalCorrection GlobalCorrection::convert(const TString& globalCorrection)
{
    if(globalCorrection == "dy") return dy;
    if(globalCorrection == "ttbb") return ttbb;
    std::cout<<"Warning! The following global correction conversion is not implemented: "<<globalCorrection<<std::endl;
    return undefined;
}



TString GlobalCorrection::convert(const GlobalCorrection& globalCorrection)
{
    if(globalCorrection == dy) return "dy";
    if(globalCorrection == ttbb) return "ttbb";
    if(globalCorrection == undefined) return "";
    std::cout<<"Error! Global correction conversion is not implemented, break...\n"<<std::endl;
    exit(97);
}



std::vector<GlobalCorrection::GlobalCorrection> GlobalCorrection::convert(const std::vector<TString>& globalCorrections)
{
    std::vector<GlobalCorrection> v_globalCorrection;
    for(auto globalCorrection : globalCorrections) v_globalCorrection.push_back(convert(globalCorrection));
    return v_globalCorrection;
}



std::vector<GlobalCorrection::GlobalCorrection> GlobalCorrection::convert(const std::vector<std::string>& globalCorrections)
{
    std::vector<GlobalCorrection> v_globalCorrection;
    for(auto globalCorrection : globalCorrections) v_globalCorrection.push_back(convert(globalCorrection));
    return v_globalCorrection;
}



std::vector<TString> GlobalCorrection::convert(const std::vector<GlobalCorrection>& globalCorrections)
{
    std::vector<TString> v_globalCorrection;
    for(auto globalCorrection : globalCorrections) v_globalCorrection.push_back(convert(globalCorrection));
    return v_globalCorrection;
}




