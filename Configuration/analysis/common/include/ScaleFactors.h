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
    class PUReweighter;
    class JetCorrectionUncertainty;
}

#include "classesFwd.h"
#include "storeTemplate.h"
#include "sampleHelpers.h"
#include "TopAnalysis/ZTopUtils/interface/bTagBase.h"







/// Namespace for functions needed only by ScaleFactors classes
namespace ScaleFactorHelpers{

    /// Get 2-dimensional scale factor from histogram
    double get2DSF(const TH2* const histo, const double x, const double y);
}







class PileupScaleFactors{
    
public:
    
    /// Constructor setting up data and MC input
    PileupScaleFactors(const std::string& inputFilename,
                       const std::string& mcEra, const std::string& pileupScenario,
                       const Systematic::Systematic& systematic);
    
    /// Destructor
    ~PileupScaleFactors(){}
    
    
    
    /// Get the scale factor for a given true vertex multiplicity
    double getSF(const size_t trueVertexMultiplicity)const;
    
    
    
private:
    
    /// Pointer to the pileup reweighter instance
    ztop::PUReweighter* const puReweighter_;
};






class LeptonScaleFactors{
    
public:
    
    /// Constructor
    LeptonScaleFactors(const char* electronSFInputFileName,
                       const char* muonSFInputFileName,
                       const Systematic::Systematic& systematic);
    
    /// Destructor
    ~LeptonScaleFactors(){}
    
    /// Get lepton per-event scale factor for exactly two leptons
    double getSFDilepton(const int leadingLeptonIndex, const int nLeadingLeptonIndex,
                         const VLV& leptons, const std::vector<int>& lepPdgIds)const;
    
    /// Get lepton per-event scale factor for all leptons in the event
    double getSFAllLeptons(const std::vector<int>& allLeptonIndices,
                           const VLV& leptons, const std::vector<int>& lepPdgIds)const;
    
    
    
private:
    
    /// Enumeration for lepton types
    enum Lepton{electron, muon};
    
    /// Enumeration for possible systematics
    enum SystematicInternal{nominal, vary_up, vary_down};
    
    /// Return the scale factor histogram
    const TH2* prepareSF(const std::string& inputFileName,
                         const std::string& histogramName,
                         const SystematicInternal& systematic)const;
    
    
    
    /// Electron scale factor histogram differential in eta, pt    
    const TH2* h2_electronSFpteta_;
    
    /// Muon scale factor histogram differential in eta, pt
    const TH2* h2_muonSFpteta_;
};






class TriggerScaleFactors{
    
public:
    
    /// Constructor
    TriggerScaleFactors(const char* inputFileSuffix,
                        const std::vector<Channel::Channel>& channels,
                        const Systematic::Systematic& systematic);
    
    /// Destructor
    ~TriggerScaleFactors(){}
    
    
    
    /// Get trigger per-event scale factor
    double getSF(const int leptonXIndex, const int leptonYIndex,
                 const VLV& leptons, const Channel::Channel& channel)const;
    
    
    
private:
    
    /// Enumeration for possible systematics
    enum SystematicInternal{nominal, vary_up, vary_down};
    
    /// Return the trigger scale factor histogram
    const TH2* prepareSF(const TString& fileName, const SystematicInternal& systematic)const;
    
    
    
    /// Trigger scale factor histogram for ee trigger, differential in eta
    const TH2* h2_eeSFeta_;
    
    /// Trigger scale factor histogram for ee trigger, differential in eta
    const TH2* h2_emuSFeta_;
    
    /// Trigger scale factor histogram for ee trigger, differential in eta
    const TH2* h2_mumuSFeta_;
};









class BtagScaleFactors : public ztop::bTagBase{
    
public:
    
    /// Enum for the implemented modes of btag corrections
    enum CorrectionMode{
        none,                       // Do not apply any corrections, i.e. scale factors event SF=1
        greaterEqualOneTagReweight, // Correct selection efficiency for given working point via event SF for >=1 b-tag
        randomNumberRetag,          // Random-number based tag flipping for b-/c-/l-jets to correct for selection efficiency
        discriminatorReweight,      // Reweight with event-wise SF to describe b-tag discriminator distribution
        undefinedMode               // Undefined
    };
    
    /// Constructor
    BtagScaleFactors(const char* btagEfficiencyInputDir,
                     const char* btagEfficiencyOutputDir,
                     const std::vector<Channel::Channel>& channels,
                     const Systematic::Systematic& systematic,
                     const CorrectionMode& correctionMode);
    
    /// Destructor
    ~BtagScaleFactors(){}
    
    
    
    /// Whether to produce b-tag efficiencies during Analysis
    /// If efficiencies can be found, they are used in the Analysis, but not produced
    /// If efficiencies cannot be found, they are not used in the Analysis, but produced
    bool makeEfficiencies();



    /// Prepare b-tagging scale factors (efficiency histograms and medians of jet eta, pt)
    void prepareSF(const Channel::Channel& channel);
    
    /// Book histograms needed for b-tag efficiencies
    void bookHistograms(TSelectorList* output);
    
    /// Fill histograms needed for b-tag efficiencies
    void fillHistograms(const std::vector<int>& jetIndices,
                        const std::vector<double>& bTagDiscriminant,
                        const VLV& jets,
                        const std::vector<int>& jetPartonFlavours,
                        const double& weight);

    /// Produce b-tag efficiencies
    void produceEfficiencies();



    /// Get b-tag per-event scale factor for selections >=1 b-tag
    double getSFGreaterEqualOneTag(const std::vector<int>& jetIndices,
                                   const VLV& jets,
                                   const std::vector<int>& jetPartonFlavours);

    /// Method takes the indices of b-tagged jets,
    /// and overwrites them with the b-tagged jets after randomised tag flipping
    /// Method explained in: https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFUtil
    /// and in: https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#2a_Jet_by_jet_updating_of_the_b
    void indexOfBtags(std::vector<int>& bjetIndices,
                      const std::vector<int>& jetIndices,
                      const VLV& jets,
                      const std::vector<int>& jetPartonFlavours,
                      const std::vector<double>& btagDiscriminants)const;

    
    
private:

    /// Pointer for bookkeeping of histograms
    // this is NOT used right now!
    // but it may be good for parallelization in the future
    TSelectorList* selectorList_;
    
    /// Map of the file names for each channel
    std::map<Channel::Channel, TString> m_channelFilename_;
    
    /// Map of the sample names for each channel
    std::map<Channel::Channel, TString> m_channelSamplename_;
    
    /// The channel which is processed
    Channel::Channel channel_;
    
    /// Which correction mode is used
    const CorrectionMode correctionMode_;
};





class JetEnergyResolutionScaleFactors{

public:

    /// Constructor
    JetEnergyResolutionScaleFactors(const Systematic::Systematic& systematic);

    /// Destructor
    ~JetEnergyResolutionScaleFactors(){}

    /// Scale the jet and MET collections
    void applySystematic(VLV* jets, VLV* jetsForMET, LV* met,
                         const std::vector<double>* jetJERSF, const std::vector<double>* jetForMETJERSF,
                         const VLV* associatedGenJet, const VLV* associatedGenJetForMET)const;



private:

    /// Enumeration for possible systematics
    enum SystematicInternal{vary_up, vary_down, undefined};

    /// The intervals in eta for granularity of scale factor
    std::vector<double> v_etaRange_;

    /// The scale factors corresponding to the eta ranges defined in v_etaRange_
    std::vector<double> v_etaScaleFactor_;
};






class JetEnergyScaleScaleFactors{

public:

    /// Constructor
    JetEnergyScaleScaleFactors(const char* jesUncertaintySourceFile,
                               const Systematic::Systematic& systematic);

    /// Destructor
    ~JetEnergyScaleScaleFactors();

    /// Scale the jet and MET collections
    void applySystematic(VLV* jets, VLV* jetsForMET, LV* met)const;



private:

    /// Enumeration for possible systematics
    enum SystematicInternal{vary_up, vary_down, undefined};

    /// Object for retrieving uncertainty values
    ztop::JetCorrectionUncertainty* jetCorrectionUncertainty_;

    /// Variation upwards?
    bool varyUp_;
};








#endif





