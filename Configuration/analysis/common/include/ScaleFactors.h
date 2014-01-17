#ifndef ScaleFactors_h
#define ScaleFactors_h

#include <string>
#include <vector>
#include <map>

class TH1;
class TH2;
class TString;
class TSelectorList;

namespace ztop{
    class JetCorrectionUncertainty;
}

#include "classesFwd.h"
#include "storeTemplate.h"




// FIXME: replace these functions with enum type ones
namespace common{

    /// Assign a folder depending on channel and systematic
    std::string assignFolder(const char* baseDir, const TString& channel, const TString& systematic);

    /// Access an already existing input folder
    std::string accessFolder(const char* baseDir, const TString& channel,
                             const TString& systematic, const bool allowNonexisting =false);
}



/// Namespace for functions needed only by ScaleFactors classes
namespace ScaleFactorHelpers{
    
    /// Get 2-dimensional scale factor from histogram
    double get2DSF(TH2* histo, const double x, const double y);
    
    /** Calculate the median of a histogram
     * 
     */
    double median(TH1* h1);
}







class LeptonScaleFactors{
    
private:
    
    /// Enumeration for lepton types
    enum Lepton{electron, muon};
    
    
    
public:
    
    /// Enumeration for possible systematics
    enum Systematic{nominal, vary_up, vary_down};
    
    /// Constructor
    LeptonScaleFactors(const char* electronSFInputFileName,
                       const char* muonSFInputFileName,
                       const Systematic& systematic);
    
    /// Destructor
    ~LeptonScaleFactors(){}
    
    /// Get lepton per-event scale factor for exactly two leptons
    double getLeptonIDSF(const int leadingLeptonIndex, const int nLeadingLeptonIndex,
                         const VLV& leptons, const std::vector<int>& lepPdgIds)const;
    
    /// Get lepton per-event scale factor for all leptons in the event
    double scaleFactorAllLeptons(const std::vector<int>& allLeptonIndices,
                                 const VLV& leptons, const std::vector<int>& lepPdgIds)const;
    
    
    
private:
    
    /// Return the scale factor histogram
    TH2* prepareLeptonIDSF(const std::string& inputFileName,
                           const std::string& histogramName,
                           const Systematic& systematic)const;
    
    /// Electron scale factor histogram differential in eta, pt    
    TH2* h2_ElectronIDSFpteta;
    /// Muon scale factor histogram differential in eta, pt
    TH2* h2_MuonIDSFpteta;
};






class TriggerScaleFactors{
    
public:
    
    /// Enumeration for possible systematics
    enum Systematic{nominal, vary_up, vary_down};
    
    /// Constructor
    TriggerScaleFactors(const char* inputFileSuffix,
                        const std::vector<std::string>& channels,
                        const Systematic& systematic);
    
    /// Destructor
    ~TriggerScaleFactors(){}
    
    
    
    /// Get trigger per-event scale factor
    double getTriggerSF(const int leptonXIndex, const int leptonYIndex,
                        const VLV& leptons, const TString& channel)const;
    
    
    
private:
    
    /// Return the trigger scale factor histogram
    TH2* prepareTriggerSF(const TString& fileName, const Systematic& systematic)const;
    
    /// Trigger scale factor histogram for ee trigger, differential in eta
    TH2* h2_eeTrigSFeta;
    
    /// Trigger scale factor histogram for ee trigger, differential in eta
    TH2* h2_emuTrigSFeta;
    
    /// Trigger scale factor histogram for ee trigger, differential in eta
    TH2* h2_mumuTrigSFeta;
};









class BtagScaleFactors{
    
public:
    
    /// Enumeration for possible systematics
    enum Systematic{nominal,
                    btag_up, btag_down,
                    btagPt_up, btagPt_down,
                    btagEta_up, btagEta_down,
                    btagLjet_up, btagLjet_down,
                    btagLjetPt_up, btagLjetPt_down,
                    btagLjetEta_up, btagLjetEta_down,
                    btagBeff_up, btagBeff_down,
                    btagCeff_up, btagCeff_down,
                    btagLeff_up, btagLeff_down};
    
    
    
    
    
    /// Constructor
    BtagScaleFactors(const char* btagEfficiencyInputDir,
                     const char* btagEfficiencyOutputDir,
                     const std::vector<std::string>& channels,
                     TString systematic);
    
    /// Destructor
    ~BtagScaleFactors(){}
    
    
    
    /// Whether to produce b-tag efficiencies during Analysis
    /// If efficiencies can be found, they are used in the Analysis, but not produced
    /// If efficiencies cannot be found, they are not used in the Analysis, but produced
    bool makeEfficiencies()const;
    
    
    
    /// Book histograms needed for b-tag efficiencies
    void bookBtagHistograms(TSelectorList* output, const std::string& channel);
    
    /// Fill histograms needed for b-tag efficiencies
    void fillBtagHistograms(const std::vector<int>& jetIndices,
                            const std::vector<int>& bjetIndices,
                            const VLV& jets,
                            const std::vector<int>& jetPartonFlavours,
                            const double weight,
                            const std::string& channel);
    
    /// Produce b-tag efficiencies
    void produceBtagEfficiencies(const std::string& channel);
    
    
    
    /// Get b-tag per-event scale factor
    double calculateBtagSF(const std::vector<int>& jetIndices,
                           const VLV& jets,
                           const std::vector<int>& jetPartonFlavours,
                           const std::string& channel)const;
    
    /// Method takes the indices of b-tagged jets,
    /// and overwrites them with the b-tagged jets after randomised tag flipping
    void indexOfBtags(std::vector<int>& bjetIndices,
                      const std::vector<int>& jetIndices,
                      const VLV& jets,
                      const std::vector<int>& jetPartonFlavours,
                      const std::vector<double>& btagDiscriminants,
                      const double TagCut,
                      const std::string& channel)const;
    
    
    
    
private:
    
    struct ChannelStruct{
        
        ChannelStruct();
        ChannelStruct(const std::string& inputFileName, const std::string& outputFileName);
        ~ChannelStruct(){}
        
        /// Store the object in the output list and return it
        template<class T> T* store(T* obj){return common::store(obj, selectorList_);}
        
        /// Input file name for files holding histograms of b-tagging efficiencies
        std::string inputFileName_;
        
        /// Output file name in case b-tagging efficiencies are produced in Analysis
        std::string outputFileName_;
        
        /// Medians in eta, pt for b-tag SF
        double btag_ptmedian_;
        double btag_etamedian_;
        
        /// Histograms of per-jet b-tagging efficiencies
        /// Need to be written out to file previously
        /// Used e.g. to calculate per-event SF
        TH2* h2_bEff;
        TH2* h2_cEff;
        TH2* h2_lEff;
        
        /// Histograms for writing out b-tagging efficiencies
        TH2* h2_bjets;
        TH2* h2_cjets;
        TH2* h2_ljets;
        TH2* h2_btaggedjets;
        TH2* h2_ctaggedjets;
        TH2* h2_ltaggedjets;
        
        /// Pointer for bookkeeping of histograms
        TSelectorList* selectorList_;    
    };
    
    /// Prepare b-tagging scale factors (efficiency histograms and medians of jet eta, pt)
    void prepareBtagSF(ChannelStruct& channelStruct);
    
    /// Methods needed for b-tag scale factor calculation
    double BJetSF(const double&, const double&)const;
    double CJetSF(const double&, const double&)const;
    double LJetSF(const double&, const double&, const TString&)const;
    double BJetSFAbsErr(const double&)const;
    double CJetSFAbsErr(const double&)const;
    
    /// 'Random' decision to tag or not tag a jet.
    /// Method explained in: https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFUtil
    /// and in: https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#2a_Jet_by_jet_updating_of_the_b
    bool isTagged(const LV& jet, const double tagValue, const int flavour,
                  const double tagCut, const std::string& channel)const;
    
    double getEfficiency(const LV& jet, const int partonFlavour,
                         const ChannelStruct& channelStruct)const;
    
    double getSF(const double pt, const double abs_eta,
                 const int flavour, const ChannelStruct& channelStruct)const;
    
    /// Decide wich type of BTag variation is going to be done according to 
    ///   systematics name 
    ///   median value (if applicable)
    double varySF(const double pt, const double eta,
                  const int flavour,
                  const double ptmedian, const double etamedian)const;
    
    /// Check if systematic belongs to a specific category
    bool ljetSystematic()const;
    bool ptUpSystematic()const;
    bool ptDownSystematic()const;
    bool etaUpSystematic()const;
    bool etaDownSystematic()const;
    bool absoluteUpSystematic()const;
    bool absoluteDownSystematic()const;
    
    /// Access the ChannelStruct from the map for a given channel
    ChannelStruct& getChannelStruct(const std::string& channel);
    
    
    
    /// Whether to produce btag efficiencies or whether to read in already produced ones
    bool makeEfficiencies_;
    
    /// Systematic to be used for btag efficiencies
    Systematic systematic_;
    
    /// Map holding the ChannelStruct for each relevant channel
    std::map<std::string, ChannelStruct> m_channelChannelStruct_;
};





class JetEnergyResolutionScaleFactors{
    
public:
    
    /// Enumeration for possible systematics
    enum Systematic{vary_up, vary_down};
    
    /// Constructor
    JetEnergyResolutionScaleFactors(const Systematic& systematic);
    
    /// Destructor
    ~JetEnergyResolutionScaleFactors(){}
    
    /// Scale the jet and MET collections
    void applySystematic(VLV* jets, VLV* jetsForMET, LV* met,
                         const std::vector<double>* jetJERSF, const std::vector<double>* jetForMETJERSF,
                         const VLV* associatedGenJet, const VLV* associatedGenJetForMET)const;
    
    
    
private:
    
    /// The intervals in eta for granularity of scale factor
    std::vector<double> v_etaRange_;
    
    /// The scale factors corresponding to the eta ranges defined in v_etaRange_
    std::vector<double> v_etaScaleFactor_;
};






class JetEnergyScaleScaleFactors{
    
public:
    
    /// Enumeration for possible systematics
    enum Systematic{vary_up, vary_down};
    
    /// Constructor
    JetEnergyScaleScaleFactors(const char* jesUncertaintySourceFile,
                               const Systematic& systematic);
    
    /// Destructor
    ~JetEnergyScaleScaleFactors();
    
    /// Scale the jet and MET collections
    void applySystematic(VLV* jets, VLV* jetsForMET, LV* met)const;
    
    
    
private:
    
    /// Object for retrieving uncertainty values
    ztop::JetCorrectionUncertainty* jetCorrectionUncertainty_;
    
    /// Variation upwards?
    bool varyUp_;
};








#endif





