#include <iostream>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <sstream>

#include <TString.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <Rtypes.h>

#include "sampleHelpers.h"






// --------------------- Functions defined in namespace Systematic for Type -------------------------



Systematic::Type Systematic::convertType(const TString& type)
{
    // Attention: the order here is important, since the first line where the BeginsWith is true is returned
    if(type.BeginsWith("Nominal")) return nominal;
    if(type.BeginsWith("mH110")) return mH110;
    if(type.BeginsWith("mH115")) return mH115;
    if(type.BeginsWith("mH120")) return mH120;
    if(type.BeginsWith("mH1225")) return mH1225;
    if(type.BeginsWith("mH1275")) return mH1275;
    if(type.BeginsWith("mH130")) return mH130;
    if(type.BeginsWith("mH135")) return mH135;
    if(type.BeginsWith("mH140")) return mH140;
    if(type.BeginsWith("LEPT")) return lept;
    if(type.BeginsWith("TRIG")) return trig;
    if(type.BeginsWith("PU")) return pu;
    if(type.BeginsWith("DY")) return dy;
    if(type.BeginsWith("BG")) return bg;
    if(type.BeginsWith("KIN")) return kin;
    if(type.BeginsWith("BTAGDISCR_BSTAT1")) return btagDiscrBstat1;
    if(type.BeginsWith("BTAGDISCR_BSTAT2")) return btagDiscrBstat2;
    if(type.BeginsWith("BTAGDISCR_LSTAT1")) return btagDiscrLstat1;
    if(type.BeginsWith("BTAGDISCR_LSTAT2")) return btagDiscrLstat2;
    if(type.BeginsWith("BTAGDISCR_PURITY")) return btagDiscrPurity;
    if(type.BeginsWith("BTAG_LJET_PT")) return btagLjetPt;
    if(type.BeginsWith("BTAG_LJET_ETA")) return btagLjetEta;
    if(type.BeginsWith("BTAG_LJET")) return btagLjet;
    if(type.BeginsWith("BTAG_BEFF")) return btagBeff;
    if(type.BeginsWith("BTAG_CEFF")) return btagCeff;
    if(type.BeginsWith("BTAG_LEFF")) return btagLeff;
    if(type.BeginsWith("BTAG_PT")) return btagPt;
    if(type.BeginsWith("BTAG_ETA")) return btagEta;
    if(type.BeginsWith("BTAG")) return btag;
    if(type.BeginsWith("JER")) return jer;
    if(type.BeginsWith("JES")) return jes;
    if(type.BeginsWith("TOP_PT")) return topPt;
    if(type.BeginsWith("MASS")) return mass;
    if(type.BeginsWith("MATCH")) return match;
    if(type.BeginsWith("SCALE")) return scale;
    if(type.BeginsWith("POWHEGHERWIG")) return powhegHerwig;
    if(type.BeginsWith("POWHEG")) return powheg;
    if(type.BeginsWith("MCATNLO")) return mcatnlo;
    if(type.BeginsWith("PERUGIA11NoCR")) return perugia11NoCR;
    if(type.BeginsWith("PERUGIA11")) return perugia11;
    if(type.BeginsWith("PDF")) return pdf;
    if(type.BeginsWith("closure")) return closure;
    if(type.BeginsWith("allAvailable")) return allAvailable;
    if(type.BeginsWith("all")) return all;
    std::cout<<"Warning! The following systematic type conversion is not implemented: "<<type<<std::endl<<std::endl;
    return undefinedType;
}



TString Systematic::convertType(const Type& type)
{
    if(type == nominal) return "Nominal";
    if(type == mH110) return "mH110";
    if(type == mH115) return "mH115";
    if(type == mH120) return "mH120";
    if(type == mH1225) return "mH1225";
    if(type == mH1275) return "mH1275";
    if(type == mH130) return "mH130";
    if(type == mH135) return "mH135";
    if(type == mH140) return "mH140";
    if(type == lept) return "LEPT";
    if(type == trig) return "TRIG";
    if(type == pu) return "PU";
    if(type == dy) return "DY";
    if(type == bg) return "BG";
    if(type == kin) return "KIN";
    if(type == btagDiscrBstat1) return "BTAGDISCR_BSTAT1";
    if(type == btagDiscrBstat2) return "BTAGDISCR_BSTAT2";
    if(type == btagDiscrLstat1) return "BTAGDISCR_LSTAT1";
    if(type == btagDiscrLstat2) return "BTAGDISCR_LSTAT2";
    if(type == btagDiscrPurity) return "BTAGDISCR_PURITY";
    if(type == btagLjetPt) return "BTAG_LJET_PT";
    if(type == btagLjetEta) return "BTAG_LJET_ETA";
    if(type == btagLjet) return "BTAG_LJET";
    if(type == btagBeff) return "BTAG_BEFF";
    if(type == btagCeff) return "BTAG_CEFF";
    if(type == btagLeff) return "BTAG_LEFF";
    if(type == btagPt) return "BTAG_PT";
    if(type == btagEta) return "BTAG_ETA";
    if(type == btag) return "BTAG";
    if(type == jer) return "JER";
    if(type == jes) return "JES";
    if(type == topPt) return "TOP_PT";
    if(type == mass) return "MASS";
    if(type == match) return "MATCH";
    if(type == scale) return "SCALE";
    if(type == powhegHerwig) return "POWHEGHERWIG";
    if(type == powheg) return "POWHEG";
    if(type == mcatnlo) return "MCATNLO";
    if(type == perugia11NoCR) return "PERUGIA11NoCR";
    if(type == perugia11) return "PERUGIA11";
    if(type == pdf) return "PDF";
    if(type == closure) return "closure";
    if(type == allAvailable) return "allAvailable";
    if(type == all) return "all";
    if(type == undefinedType) return "";
    std::cerr<<"Error! Type conversion is not implemented,\n...break\n"<<std::endl;
    exit(99);
}



std::vector<Systematic::Type> Systematic::convertType(const std::vector<TString>& types)
{
    std::vector<Type> v_type;
    for(const auto& type : types) v_type.push_back(convertType(type));
    return v_type;
}



std::vector<Systematic::Type> Systematic::convertType(const std::vector<std::string>& types)
{
    std::vector<Type> v_type;
    for(const auto& type : types) v_type.push_back(convertType(type));
    return v_type;
}



std::vector<TString> Systematic::convertType(const std::vector<Type>& types)
{
    std::vector<TString> v_type;
    for(const auto& type : types) v_type.push_back(convertType(type));
    return v_type;
}










// --------------------- Functions defined in namespace Systematic for Variation -------------------------





Systematic::Variation Systematic::convertVariation(const TString& variation)
{
    if(variation.EndsWith("_UP")) return up;
    if(variation.EndsWith("_DOWN")) return down;
    if(variation.EndsWith("_CENTRAL")) return central;
    //std::cout<<"Warning! The following variation conversion is not implemented: "<<variation<<std::endl<<std::endl;
    return undefinedVariation;
}



TString Systematic::convertVariation(const Variation& variation)
{
    if(variation == up) return "_UP";
    if(variation == down) return "_DOWN";
    if(variation == central) return "_CENTRAL";
    if(variation == undefinedVariation) return "";
    std::cerr<<"Error! Variation conversion is not implemented,\n...break\n"<<std::endl;
    exit(99);
}



std::vector<Systematic::Variation> Systematic::convertVariation(const std::vector<TString>& variations)
{
    std::vector<Variation> v_variation;
    for(const auto& variation : variations) v_variation.push_back(convertVariation(variation));
    return v_variation;
}



std::vector<Systematic::Variation> Systematic::convertVariation(const std::vector<std::string>& variations)
{
    std::vector<Variation> v_variation;
    for(const auto& variation : variations) v_variation.push_back(convertVariation(variation));
    return v_variation;
}



std::vector<TString> Systematic::convertVariation(const std::vector<Variation>& variations)
{
    std::vector<TString> v_variation;
    for(const auto& variation : variations) v_variation.push_back(convertVariation(variation));
    return v_variation;
}










// --------------------- Further functions defined in namespace Systematic -------------------------



void Systematic::isValid(const Type& type, const Variation& variation, const int variationNumber)
{
    // Check validity of variationNumber
    if(variationNumber >= 0){
        if(std::find(centralTypes.begin(), centralTypes.end(), type) == centralTypes.end()){
            std::cerr<<"ERROR in Systematic::isValid()! Given type does not allow variation numbers (type, variationNumber): "
                     <<convertType(type)<<" , "<<variationNumber<<"\n...break\n"<<std::endl;
            exit(7);
        }
    }
    
    // Check validity of variation
    if(variation == undefinedVariation) return;
    else if(variation==up || variation==down){
        if(std::find(upDownTypes.begin(), upDownTypes.end(), type) == upDownTypes.end()){
            std::cerr<<"ERROR in Systematic::isValid()! Given type does not allow variation (type, variation): "
                     <<convertType(type)<<" , "<<convertVariation(variation)<<"\n...break\n"<<std::endl;
            exit(7);
        }
    }
    else if(variation == central){
        if(std::find(centralTypes.begin(), centralTypes.end(), type) == centralTypes.end()){
            std::cerr<<"ERROR in Systematic::isValid()! Given type does not allow variation (type, variation): "
                     <<convertType(type)<<" , "<<convertVariation(variation)<<"\n...break\n"<<std::endl;
            exit(7);
        }
    }
    else{
        std::cerr<<"ERROR in Systematic::isValid()! Variation is not defined for validity check: "
                 <<convertVariation(variation)<<"\n...break\n"<<std::endl;
        exit(7);
    }
}




std::vector<Systematic::Systematic> Systematic::allowedSystematicsAnalysis(const std::vector<Type>& allowedTypes)
{
    std::vector<Systematic> result;
    
    for(const Type& type : allowedTypes){
        // Exclude non-real types
        if(type==all || type==allAvailable) continue;
        
        if(std::find(centralTypes.begin(), centralTypes.end(), type) != centralTypes.end()){
            // Central types need specific treatment using variation numbers, e.g. PDF variations
            // They require detailed specifications at the place where they are used
            result.push_back(Systematic(type, undefinedVariation));
        }
        else if(std::find(upDownTypes.begin(), upDownTypes.end(), type) != upDownTypes.end()){
            // Up/down types need the two variations
            result.push_back(Systematic(type, up));
            result.push_back(Systematic(type, down));
        }
        else{
            // All others have no variations
            result.push_back(Systematic(type, undefinedVariation));
        }
    }
    
    return result;
}



std::vector<Systematic::Systematic> Systematic::setSystematics(const std::vector<std::string>& systematicNames)
{
    std::vector<Systematic> result;
    
    for(const auto& name : systematicNames) result.push_back(Systematic(name));
    
    return result;
}


Systematic::Systematic Systematic::nominalSystematic()
{
    return Systematic(nominal, undefinedVariation);
}



Systematic::Systematic Systematic::undefinedSystematic()
{
    return Systematic(undefinedType, undefinedVariation);
}




// --------------------- Methods of class Systematic in namespace Systematic -------------------------



Systematic::Systematic::Systematic():
type_(undefinedType),
variation_(undefinedVariation),
variationNumber_(-1)
{}



Systematic::Systematic::Systematic(const Type& type, const Variation& variation, const int variationNumber):
type_(type),
variation_(variation),
variationNumber_(variationNumber)
{}



Systematic::Systematic::Systematic(const TString& systematicName):
type_(undefinedType),
variation_(undefinedVariation),
variationNumber_(-1)
{
    TString fragment(systematicName);
    type_ = convertType(systematicName);
    fragment.ReplaceAll(convertType(type_), "");
    variation_ = convertVariation(fragment);
    fragment.ReplaceAll(convertVariation(variation_), "");
    if(fragment != ""){
        fragment.ReplaceAll("_", "");
        int variationNumber = -1;
        std::stringstream stream(fragment.Data());
        if(!(stream>>variationNumber)){
            std::cerr<<"ERROR in constructor of Systematic! Could not fragment systematic name (name --- type, variation, variationNumber): "
                     <<systematicName<<" --- "<<convertType(type_)<<" , "<<convertVariation(variation_)<<" , "<<variationNumber_<<"\n...break\n"<<std::endl;
            exit(8);
        }
        variationNumber_ = variationNumber;
        if(variationNumber_ < 0){
            std::cerr<<"ERROR in constructor of Systematic! Variation numbers must be >=0, but is (systematicName, extracted variationNumber): "
                     <<systematicName<<" , "<<variationNumber_<<"\n...break\n"<<std::endl;
            exit(8);
        }
    }
    isValid(type_, variation_, variationNumber_);
}



TString Systematic::Systematic::name()const
{
    TString result = convertType(type_);
    if(variationNumber_ >= 0){
        std::stringstream stream;
        stream<<"_"<<variationNumber_;
        result.Append(stream.str());
    }
    result.Append(convertVariation(variation_));
    return result;
}










// --------------------- Functions defined in namespace Channel -------------------------



Channel::Channel Channel::convert(const TString& channel)
{
    if(channel == "ee") return ee;
    if(channel == "emu") return emu;
    if(channel == "mumu") return mumu;
    if(channel == "combined") return combined;
    if(channel == "tautau") return tautau;
    if(channel == "") return undefined;
    std::cerr<<"Error! The following channel conversion is not implemented: "<<channel<<"\n...break\n"<<std::endl;
    exit(98);
}



TString Channel::convert(const Channel& channel)
{
    if(channel == ee) return "ee";
    if(channel == emu) return "emu";
    if(channel == mumu) return "mumu";
    if(channel == combined) return "combined";
    if(channel == tautau) return "tautau";
    if(channel == undefined) return "";
    std::cerr<<"Error! Channel conversion is not implemented,\n...break\n"<<std::endl;
    exit(98);
}



TString Channel::label(const Channel& channel)
{
    if(channel == ee) return "ee";
    if(channel == emu) return "e#mu";
    if(channel == mumu) return "#mu#mu";
    if(channel == combined) return "Dilepton Combined";
    if(channel == tautau) return "#tau#tau";
    if(channel == undefined) return "";
    std::cerr<<"Error! Channel label is not implemented,\n...break\n"<<std::endl;
    exit(98);
}



std::vector<Channel::Channel> Channel::convert(const std::vector<TString>& channels)
{
    std::vector<Channel> v_channel;
    for(const auto& channel : channels) v_channel.push_back(convert(channel));
    return v_channel;
}



std::vector<Channel::Channel> Channel::convert(const std::vector<std::string>& channels)
{
    std::vector<Channel> v_channel;
    for(const auto& channel : channels) v_channel.push_back(convert(channel));
    return v_channel;
}



std::vector<TString> Channel::convert(const std::vector<Channel>& channels)
{
    std::vector<TString> v_channel;
    for(const auto& channel : channels) v_channel.push_back(convert(channel));
    return v_channel;
}









// --------------------- Functions defined in namespace common -------------------------



TString common::assignFolder(const char* baseDir, const Channel::Channel& channel, const Systematic::Systematic& systematic)
{
    TString path("");
    
    // Create all subdirectories contained in baseDir
    TObjArray* a_subDir = TString(baseDir).Tokenize("/");
    for(Int_t iSubDir = 0; iSubDir < a_subDir->GetEntriesFast(); ++iSubDir){
        const TString& subDir = a_subDir->At(iSubDir)->GetName();
        path.Append(subDir);
        path.Append("/");
        gSystem->MakeDirectory(path);
    }
    
    // Create subdirectories for systematic and channel
    path.Append(systematic.name());
    path.Append("/");
    gSystem->MakeDirectory(path);
    path.Append(Channel::convert(channel));
    path.Append("/");
    gSystem->MakeDirectory(path);
    
    // FIXME: why not using directly gSystem->mkdir("...", true);   ???
    // Should recursively create all needed directories
    
    return path;
}



TString common::accessFolder(const char* baseDir, const Channel::Channel& channel,
                             const Systematic::Systematic& systematic, const bool allowNonexisting)
{
    // Build directory path
    TString path(baseDir);
    path.Append("/");
    path.Append(systematic.name());
    path.Append("/");
    path.Append(Channel::convert(channel));
    path.Append("/");
    
    // Check if directory really exists
    if(!gSystem->OpenDirectory(path)){
        if(allowNonexisting){
            // It is allowed to request a folder which does not exist, so return empty string silently
            return "";
        }
        else{
            std::cerr<<"ERROR! Request to access directory is not possible, because it does not exist. Directory name: "
                     <<path<<"\n...break\n"<<std::endl;
            exit(237);
        }
    }
    
    return path;
}



Channel::Channel common::finalState(const TString& filename)
{
    std::vector<Channel::Channel> v_channel {Channel::ee, Channel::emu, Channel::mumu};
    for(auto channel : v_channel){
        TString finalState(Channel::convert(channel));
        finalState.Prepend("/");
        finalState.Append("/");
        if(filename.Contains(finalState)){
            return channel;
        }
    }
    return Channel::undefined;
}



TString common::findFilelist(const TString& filelistDirectory,
                             const Channel::Channel& channel,
                             const Systematic::Systematic& systematic)
{
    TString result("");
    
    const TString filelistName(filelistDirectory + "/HistoFileList_" + systematic.name() + "_" + Channel::convert(channel) + ".txt");
    std::ifstream fileList(filelistName);
    if(!fileList.fail()){
        result = filelistName;
        fileList.close();
    }
    
    return result;
}



std::vector<Systematic::Systematic> common::findSystematicsFromFilelists(const TString& filelistDirectory,
                                                                         const std::vector<Channel::Channel>& v_channel,
                                                                         const std::vector<Systematic::Systematic>& v_systematic)
{
    std::vector<Systematic::Systematic> result;
    
    for(const auto& systematic : v_systematic){
        bool foundSystematic(true);
        for(const auto& channel : v_channel){
            const TString filelistName = findFilelist(filelistDirectory, channel, systematic);
            if(filelistName == ""){
                foundSystematic = false;
                break;
            }
        }
        if(foundSystematic) result.push_back(systematic);
    }
    
    return result;
}



std::vector<TString> common::readFilelist(const TString& filelistDirectory,
                                          const Channel::Channel& channel,
                                          const Systematic::Systematic& systematic,
                                          const std::vector<TString>& v_pattern)
{
    // Access name of file list containing list of input root files
    const TString filelistName = findFilelist(filelistDirectory, channel, systematic);
    if(filelistName == ""){
        std::cerr<<"Error in common::readFilelist! Cannot find file (folder, channel, systematic): "
                 <<filelistDirectory<<" , "<<Channel::convert(channel)<<" , "<<systematic.name()<<"\n...break\n"<<std::endl;
        exit(1);
    }
    
    // Read in file list to a vector
    std::vector<TString> v_filename;
    std::cout<<"Reading file: "<<filelistName<<std::endl;
    std::ifstream filelist(filelistName);
    while(!filelist.eof()){
        TString filename;
        filelist>>filename;
        // Skip empty lines
        if(filename == "") continue;
        // Comment lines in FileList with '#'
        if(filename.BeginsWith("#")) continue;
        // Check that all patterns are contained in the filename
        for(const auto& pattern : v_pattern) if(!filename.Contains(pattern)) continue;
        v_filename.push_back(filename);
    }
    filelist.close();
    
    return v_filename;
}








