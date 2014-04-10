#include <iostream>
#include <cstdlib>

#include <TString.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <Rtypes.h>

#include "sampleHelpers.h"






// --------------------- Functions defined in namespace Systematic -------------------------



Systematic::Systematic Systematic::convertSystematic(const std::string& systematic)
{
    if(systematic == "Nominal") return nominal;
    if(systematic == "mH110") return mH110;
    if(systematic == "mH115") return mH115;
    if(systematic == "mH120") return mH120;
    if(systematic == "mH1225") return mH1225;
    if(systematic == "mH1275") return mH1275;
    if(systematic == "mH130") return mH130;
    if(systematic == "mH135") return mH135;
    if(systematic == "mH140") return mH140;
    if(systematic == "LEPT_UP") return lept_up;
    if(systematic == "LEPT_DOWN") return lept_down;
    if(systematic == "TRIG_UP") return trig_up;
    if(systematic == "TRIG_DOWN") return trig_down;
    if(systematic == "PU_UP") return pu_up;
    if(systematic == "PU_DOWN") return pu_down;
    if(systematic == "DY_UP") return dy_up;
    if(systematic == "DY_DOWN") return dy_down;
    if(systematic == "BG_UP") return bg_up;
    if(systematic == "BG_DOWN") return bg_down;
    if(systematic == "KIN_UP") return kin_up;
    if(systematic == "KIN_DOWN") return kin_down;
    if(systematic == "BTAG_UP") return btag_up;
    if(systematic == "BTAG_DOWN") return btag_down;
    if(systematic == "BTAG_PT_UP") return btagPt_up;
    if(systematic == "BTAG_PT_DOWN") return btagPt_down;
    if(systematic == "BTAG_ETA_UP") return btagEta_up;
    if(systematic == "BTAG_ETA_DOWN") return btagEta_down;
    if(systematic == "BTAG_LJET_UP") return btagLjet_up;
    if(systematic == "BTAG_LJET_DOWN") return btagLjet_down;
    if(systematic == "BTAG_LJET_PT_UP") return btagLjetPt_up;
    if(systematic == "BTAG_LJET_PT_DOWN") return btagLjetPt_down;
    if(systematic == "BTAG_LJET_ETA_UP") return btagLjetEta_up;
    if(systematic == "BTAG_LJET_ETA_DOWN") return btagLjetEta_down;
    if(systematic == "BTAG_BEFF_UP") return btagBeff_up;
    if(systematic == "BTAG_BEFF_DOWN") return btagBeff_down;
    if(systematic == "BTAG_CEFF_UP") return btagCeff_up;
    if(systematic == "BTAG_CEFF_DOWN") return btagCeff_down;
    if(systematic == "BTAG_LEFF_UP") return btagLeff_up;
    if(systematic == "BTAG_LEFF_DOWN") return btagLeff_down;
    if(systematic == "JER_UP") return jer_up;
    if(systematic == "JER_DOWN") return jer_down;
    if(systematic == "JES_UP") return jes_up;
    if(systematic == "JES_DOWN") return jes_down;
    if(systematic == "MASS_UP") return mass_up;
    if(systematic == "MASS_DOWN") return mass_down;
    if(systematic == "MATCH_UP") return match_up;
    if(systematic == "MATCH_DOWN") return match_down;
    if(systematic == "SCALE_UP") return scale_up;
    if(systematic == "SCALE_DOWN") return scale_down;
    if(systematic == "POWHEG") return powheg;
    if(systematic == "MCATNLO") return mcatnlo;
    if(systematic == "PDF") return pdf;
    if(systematic == "all") return all;
    std::cout<<"Warning! The following systematic conversion is not implemented: "<<systematic<<std::endl<<std::endl;
    return undefined;
}



std::string Systematic::convertSystematic(const Systematic& systematic)
{
    if(systematic == nominal) return "Nominal";
    if(systematic == mH110) return "mH110";
    if(systematic == mH115) return "mH115";
    if(systematic == mH120) return "mH120";
    if(systematic == mH1225) return "mH1225";
    if(systematic == mH1275) return "mH1275";
    if(systematic == mH130) return "mH130";
    if(systematic == mH135) return "mH135";
    if(systematic == mH140) return "mH140";
    if(systematic == lept_up) return "LEPT_UP";
    if(systematic == lept_down) return "LEPT_DOWN";
    if(systematic == trig_up) return "TRIG_UP";
    if(systematic == trig_down) return "TRIG_DOWN";
    if(systematic == pu_up) return "PU_UP";
    if(systematic == pu_down) return "PU_DOWN";
    if(systematic == dy_up) return "DY_UP";
    if(systematic == dy_down) return "DY_DOWN";
    if(systematic == bg_up) return "BG_UP";
    if(systematic == bg_down) return "BG_DOWN";
    if(systematic == kin_up) return "KIN_UP";
    if(systematic == kin_down) return "KIN_DOWN";
    if(systematic == btag_up) return "BTAG_UP";
    if(systematic == btag_down) return "BTAG_DOWN";
    if(systematic == btagPt_up) return "BTAG_PT_UP";
    if(systematic == btagPt_down) return "BTAG_PT_DOWN";
    if(systematic == btagEta_up) return "BTAG_ETA_UP";
    if(systematic == btagEta_down) return "BTAG_ETA_DOWN";
    if(systematic == btagLjet_up) return "BTAG_LJET_UP";
    if(systematic == btagLjet_down) return "BTAG_LJET_DOWN";
    if(systematic == btagLjetPt_up) return "BTAG_LJET_PT_UP";
    if(systematic == btagLjetPt_down) return "BTAG_LJET_PT_DOWN";
    if(systematic == btagLjetEta_up) return "BTAG_LJET_ETA_UP";
    if(systematic == btagLjetEta_down) return "BTAG_LJET_ETA_DOWN";
    if(systematic == btagBeff_up) return "BTAG_BEFF_UP";
    if(systematic == btagBeff_down) return "BTAG_BEFF_DOWN";
    if(systematic == btagCeff_up) return "BTAG_CEFF_UP";
    if(systematic == btagCeff_down) return "BTAG_CEFF_DOWN";
    if(systematic == btagLeff_up) return "BTAG_LEFF_UP";
    if(systematic == btagLeff_down) return "BTAG_LEFF_DOWN";
    if(systematic == jer_up) return "JER_UP";
    if(systematic == jer_down) return "JER_DOWN";
    if(systematic == jes_up) return "JES_UP";
    if(systematic == jes_down) return "JES_DOWN";
    if(systematic == mass_up) return "MASS_UP";
    if(systematic == mass_down) return "MASS_DOWN";
    if(systematic == match_up) return "MATCH_UP";
    if(systematic == match_down) return "MATCH_DOWN";
    if(systematic == scale_up) return "SCALE_UP";
    if(systematic == scale_down) return "SCALE_DOWN";
    if(systematic == powheg) return "POWHEG";
    if(systematic == mcatnlo) return "MCATNLO";
    if(systematic == pdf) return "PDF";
    if(systematic == all) return "all";
    if(systematic == undefined) return "";
    std::cout<<"Error! Systematic conversion is not implemented, break...\n"<<std::endl;
    exit(99);
}



std::vector<Systematic::Systematic> Systematic::convertSystematics(const std::vector<std::string>& systematics)
{
    std::vector<Systematic> v_systematic;
    for(const auto& systematic : systematics) v_systematic.push_back(convertSystematic(systematic));
    return v_systematic;
}



std::vector<std::string> Systematic::convertSystematics(const std::vector<Systematic>& systematics)
{
    std::vector<std::string> v_systematic;
    for(const auto& systematic : systematics) v_systematic.push_back(convertSystematic(systematic));
    return v_systematic;
}








// --------------------- Functions defined in namespace Channel -------------------------



Channel::Channel Channel::convertChannel(const std::string& channel)
{
    if(channel == "ee") return ee;
    if(channel == "emu") return emu;
    if(channel == "mumu") return mumu;
    if(channel == "combined") return combined;
    if(channel == "tautau") return tautau;
    std::cout<<"Warning! The following channel conversion is not implemented: "<<channel<<std::endl<<std::endl;
    return undefined;
}



std::string Channel::convertChannel(const Channel& channel)
{
    if(channel == ee) return "ee";
    if(channel == emu) return "emu";
    if(channel == mumu) return "mumu";
    if(channel == combined) return "combined";
    if(channel == tautau) return "tautau";
    if(channel == undefined) return "";
    std::cout<<"Error! Channel conversion is not implemented, break...\n"<<std::endl;
    exit(98);
}



std::string Channel::channelLabel(const Channel& channel)
{
    if(channel == ee) return "ee";
    if(channel == emu) return "e#mu";
    if(channel == mumu) return "#mu#mu";
    if(channel == combined) return "Dilepton Combined";
    if(channel == tautau) return "#tau#tau";
    if(channel == undefined) return "";
    std::cout<<"Error! Channel label is not implemented, break...\n"<<std::endl;
    exit(97);
}



std::vector<Channel::Channel> Channel::convertChannels(const std::vector<std::string>& channels)
{
    std::vector<Channel> v_channel;
    for(const auto& channel : channels) v_channel.push_back(convertChannel(channel));
    return v_channel;
}



std::vector<std::string> Channel::convertChannels(const std::vector<Channel>& channels)
{
    std::vector<std::string> v_channel;
    for(const auto& channel : channels) v_channel.push_back(convertChannel(channel));
    return v_channel;
}









// --------------------- Functions defined in namespace common -------------------------



TString common::assignFolder(const char* baseDir, const Channel::Channel& channel, const Systematic::Systematic& systematic)
{
    std::string path("");
    
    // Create all subdirectories contained in baseDir
    TObjArray* a_subDir = TString(baseDir).Tokenize("/");
    for(Int_t iSubDir = 0; iSubDir < a_subDir->GetEntriesFast(); ++iSubDir){
        const TString& subDir = a_subDir->At(iSubDir)->GetName();
        path.append(subDir);
        path.append("/");
        gSystem->MakeDirectory(path.c_str());
    }
    
    // Create subdirectories for systematic and channel
    path.append(Systematic::convertSystematic(systematic));
    path.append("/");
    gSystem->MakeDirectory(path.c_str());
    path.append(Channel::convertChannel(channel));
    path.append("/");
    gSystem->MakeDirectory(path.c_str());
    
    // FIXME: why not using directly gSystem->mkdir("...", true);   ???
    // Should recursively create all needed directories
    
    return path;
}



TString common::accessFolder(const char* baseDir, const Channel::Channel& channel,
                            const Systematic::Systematic& systematic, const bool allowNonexisting)
{
    // Build directory path
    std::string path(baseDir);
    path.append("/");
    path.append(Systematic::convertSystematic(systematic));
    path.append("/");
    path.append(Channel::convertChannel(channel));
    path.append("/");
    
    // Check if directory really exists
    if(!gSystem->OpenDirectory(path.c_str())){
        if(allowNonexisting){
            // It is allowed to request a folder which does not exist, so return empty string silently
            return "";
        }
        else{
            std::cerr<<"ERROR! Request to access directory is not possible, because it does not exist. Directory name: "<<path
                     <<"\n...break\n"<<std::endl;
            exit(237);
        }
    }
    
    return path;
}



Channel::Channel common::finalState(const TString& filename)
{
    std::vector<Channel::Channel> v_channel {Channel::ee, Channel::emu, Channel::mumu};
    for(auto channel : v_channel){
        TString finalState(Channel::convertChannel(channel));
        finalState.Prepend("/");
        finalState.Append("/");
        if(filename.Contains(finalState)){
            return channel;
        }
    }
    return Channel::undefined;
}









