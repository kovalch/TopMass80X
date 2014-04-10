#ifndef sampleHelpers_h
#define sampleHelpers_h

#include <string>
#include <vector>

class TString;







/// Namespace to treat systematics as enum types
namespace Systematic{
    
    /// All systematics as needed in any part of the framework
    enum Systematic{
        nominal,    // nominal, i.e. no systematic variation applied
        mH110,      // Higgs mass of 110 GeV
        mH115,      // Higgs mass of 115 GeV
        mH120,      // Higgs mass of 120 GeV
        mH1225,     // Higgs mass of 122.5 GeV
        mH1275,     // Higgs mass of 127.5 GeV
        mH130,      // Higgs mass of 130 GeV
        mH135,      // Higgs mass of 135 GeV
        mH140,      // Higgs mass of 140 GeV
        lept_up,    // scale up lepton ID/ISO data-to-MC scale factors
        lept_down,  // scale down lepton ID/ISO data-to-MC scale factors
        trig_up,    // scale up trigger data-to-MC scale factors
        trig_down,  // scale down trigger ID/ISO data-to-MC scale factors
        pu_up,
        pu_down,
        dy_up,
        dy_down,
        bg_up,
        bg_down,
        kin_up,
        kin_down,
        btag_up,    // scale up b-tagging data-to-MC scale factors of the b-/c-jets
        btag_down,  // scale down b-tagging data-to-MC scale factors of the b-/c-jets
        btagPt_up,
        btagPt_down,
        btagEta_up,
        btagEta_down,
        btagLjet_up,   // scale up b-tagging data-to-MC scale factors of the l-jets
        btagLjet_down, // scale down b-tagging data-to-MC scale factors of the l-jets
        btagLjetPt_up,
        btagLjetPt_down,
        btagLjetEta_up,
        btagLjetEta_down,
        btagBeff_up,
        btagBeff_down,
        btagCeff_up,
        btagCeff_down,
        btagLeff_up,
        btagLeff_down,
        jer_up,     // scale up jet energy resolution scale factors
        jer_down,   // scale down jet energy resolution scale factors
        jes_up,     // scale up jet energy scale scale factors
        jes_down,   // scale down jet energy scale scale factors
        mass_up,
        mass_down,
        match_up,
        match_down,
        scale_up,
        scale_down,
        powheg,
        mcatnlo,
        pdf,        // PDF variations
        all,        // All allowed systematics
        undefined   // No systematic defined (also not nominal)
    };
    
    
    
    /// All systematics allowed for analysis step in Top analysis
    /// Only systematics which run on the nominal ntuples, e.g. pileup reweighting
    /// (allow undefined to set default behaviour if no option is set, i.e. option is empty)
    const std::vector<Systematic> allowedSystematicsTopAnalysis
        {nominal, undefined};
    
    /// All systematics allowed for plotting step in Top analysis
    const std::vector<Systematic> allowedSystematicsTopPlotting
        {nominal};
    
    
    
    /// All systematics allowed for analysis step in Higgs analysis
    /// Only systematics which run on the nominal ntuples, e.g. pileup reweighting
    /// (allow undefined to set default behaviour if no option is set, i.e. option is empty)
    const std::vector<Systematic> allowedSystematicsHiggsAnalysis
        {nominal, lept_up, lept_down, trig_up, trig_down, undefined};
    
    /// All systematics allowed for plotting step in Higgs analysis
    const std::vector<Systematic> allowedSystematicsHiggsPlotting
        {nominal, mH110, mH115, mH120, mH1225, mH1275, mH130, mH135, mH140};
    
    
    
    /// Convert a systematic from string to typedef
    Systematic convertSystematic(const std::string& systematic);
    
    /// Convert a systematic from typedef to string
    std::string convertSystematic(const Systematic& systematic);
    
    /// Convert a vector of systematics from string to typedef
    std::vector<Systematic> convertSystematics(const std::vector<std::string>& systematics);
    
    /// Convert a vector of systematics from string to typedef
    std::vector<std::string> convertSystematics(const std::vector<Systematic>& systematics);
}







/// Namespace to treat decay channels as enum types
namespace Channel{
    
    /// All dileptonic decay channels as needed in any part of the framework
    enum Channel{ee, emu, mumu, combined, tautau, undefined};
    
    
    
    /// All dileptonic decay channels allowed for analysis step
    /// (allow undefined to select all channels if no option is set, i.e. option is empty)
    const std::vector<Channel> allowedChannelsAnalysis
        {ee, emu, mumu, undefined};
    
    /// All dileptonic decay channels allowed for plotting step
    const std::vector<Channel> allowedChannelsPlotting
        {ee, emu, mumu, combined};
    
    /// Real analysis channels, i.e. all channels which describe a real final state
    const std::vector<Channel> realChannels
        {ee, emu, mumu};
    
    /// Possible Drell-Yan decay channels
    const std::vector<Channel> dyDecayChannels
        {ee, mumu, tautau};
    
    
    
    /// Convert a channel from string to typedef
    Channel convertChannel(const std::string& channel);
    
    /// Convert a channel from typedef to string
    std::string convertChannel(const Channel& channel);
    
    /// Return the label of a channel as used for drawing
    std::string channelLabel(const Channel& channel);
    
    /// Convert a vector of channels from string to typedef
    std::vector<Channel> convertChannels(const std::vector<std::string>& channels);
    
    /// Convert a vector of channels from string to typedef
    std::vector<std::string> convertChannels(const std::vector<Channel>& channels);
}







namespace common{
    
    /// Create and assign an output folder depending on the channel and systematic
    TString assignFolder(const char* baseDir, const Channel::Channel& channel, const Systematic::Systematic& systematic);
    
    /// Access an already existing input folder
    TString accessFolder(const char* baseDir, const Channel::Channel& channel,
                         const Systematic::Systematic& systematic, const bool allowNonexisting =false);
    
    /// Access the real final state from a filename, ie. only "ee", "emu", "mumu", but not "combined"
    Channel::Channel finalState(const TString& filename);
}







#endif






