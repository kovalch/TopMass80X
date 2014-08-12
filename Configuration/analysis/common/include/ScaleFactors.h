#ifndef ScaleFactors_h
#define ScaleFactors_h

#include <string>
#include <vector>
#include <map>
#include <functional>

class TH1;
class TH2;
class TString;
class TSelectorList;

namespace ztop{
    class PUReweighter;
    class JetCorrectionUncertainty;
}

#include "classesFwd.h"
#include "ScaleFactorsFwd.h"
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
    
    
    
    /// Prepare trigger scale factors per channel
    void prepareChannel(const Channel::Channel& channel);
    
    /// Get trigger per-event scale factor
    double getSF(const int leptonXIndex, const int leptonYIndex,
                 const VLV& leptons)const;
    
    
    
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
    
    
    
    /// Histogram of the channel as set in prepareChannel()
    const TH2* h2_channelSF_;
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
    BtagScaleFactors(const char* efficiencyInputDir,
                     const char* efficiencyOutputDir,
                     const char* inputFileHeavyFlavour,
                     const char* inputFileLightFlavour,
                     const std::vector<Channel::Channel>& channels,
                     const Systematic::Systematic& systematic,
                     const CorrectionMode& correctionMode);
    
    /// Destructor
    ~BtagScaleFactors(){}
    
    
    
    /// Prepare b-tagging scale factors per channel
    /// if a correction mode is applied which needs channel specific settings and can thus not be fully configured in constructor
    void prepareChannel(const Channel::Channel& channel);
    
    /// Set the b-tag algorithm and working point, in order to take associated values centrally from bTagBase
    void algorithmAndWorkingPoint(const Btag::Algorithm& algorithm,
                                  const Btag::WorkingPoint& workingPoint);
    
    
    
    /// Whether to produce b-tag efficiencies during Analysis
    /// If efficiencies can be found, they are used in the Analysis, but not produced
    /// If efficiencies cannot be found, they are not used in the Analysis, but produced
    bool makeEfficiencies()const;
    
    /// Book histograms needed for b-tag efficiencies
    /// if a correction mode is applied which needs b-tag efficiencies,
    /// and if they need to be produced (i.e. no file containing them found)
    void bookEfficiencyHistograms(TSelectorList* output);
    
    /// Fill histograms needed for b-tag efficiencies
    /// if a correction mode is applied which needs b-tag efficiencies,
    /// and if they need to be produced (i.e. no file containing them found)
    void fillEfficiencyHistograms(const std::vector<int>& jetIndices,
                                  const std::vector<double>& btagDiscriminators,
                                  const VLV& jets,
                                  const std::vector<int>& jetPartonFlavours,
                                  const double& weight);
    
    /// Produce b-tag efficiencies
    /// if a correction mode is applied which needs b-tag efficiencies,
    /// and if they need to be produced (i.e. no file containing them found)
    void produceEfficiencies();
    
    
    
    
    /// Get b-tag per-event scale factor
    /// if a correction mode is applied which requires scale factors
    /// else return 1.
    double getSF(const std::vector<int>& jetIndices,
                 const VLV& jets,
                 const std::vector<int>& jetPartonFlavours,
                 const std::vector<double>& btagDiscriminators);
    
    /// Method takes the indices of b-tagged jets,
    /// and overwrites them with the b-tagged jets after randomised tag flipping
    /// if the corresponding correction mode is applied
    /// Method explained in: https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFUtil
    /// and in: https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#2a_Jet_by_jet_updating_of_the_b
    void indexOfBtags(std::vector<int>& bjetIndices,
                      const std::vector<int>& jetIndices,
                      const VLV& jets,
                      const std::vector<int>& jetPartonFlavours,
                      const std::vector<double>& btagDiscriminants)const;
    
    
    
private:
    
    /// Define the internal systematic of the tool
    /// In case the systematic is a variation of b-tagging related effects set it, else use nominal
    void btagSystematic(const Systematic::Systematic& systematic);
    
    /// Prepare the b-tag efficiencies in MC, i.e. check whether they exist already or whether they need to be produced
    /// and set up the tool accordingly
    void prepareEfficiencies(const char* efficiencyInputDir,
                             const char* efficiencyOutputDir,
                             const std::vector<Channel::Channel>& channels,
                             const Systematic::Systematic& systematic);
    
    /// Prepare the b-tag efficiencies in MC per channel, i.e. read the proper histograms for application of the scale factors
    void loadEfficiencies();
    
    /// Prepare the discriminator reweighting, i.e. check whether the file with the official scale factors can be found
    void prepareDiscriminatorReweighting(const char* inputFileHeavyFlavour, const char* inputFileLightFlavour);
    
    
    /// Get b-tag per-event scale factor for selections >=1 b-tag
    double scaleFactorGreaterEqualOneTag(const std::vector<int>& jetIndices,
                                         const VLV& jets,
                                         const std::vector<int>& jetPartonFlavours);
    
    /// Get b-tag per-event scale factor for reweighting the b-tag discriminator distribution
    double scaleFactorDiscriminatorReweight(const std::vector<int>& jetIndices,
                                            const VLV& jets,
                                            const std::vector<int>& jetPartonFlavours,
                                            const std::vector<double>& btagDiscriminators)const;
    
    
    
    /// Pointer for bookkeeping of histograms
    // this is NOT used right now!
    // but it may be good for parallelization in the future
    TSelectorList* selectorList_;
    
    /// Map of the file names for each channel
    std::map<Channel::Channel, TString> m_channelFilename_;
    
    /// Map of the sample names for each channel
    std::map<Channel::Channel, std::string> m_channelSamplename_;
    
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
    
    
    
    /// Scale the jet collection
    void applyJetSystematic(VLV* const v_jet,
                            const std::vector<double>* const v_jetJerSF,
                            const VLV* const v_associatedGenJet)const;
    
    /// Scale the MET collection
    void applyMetSystematic(VLV* const v_jetForMet, LV* const met,
                            const std::vector<double>* const v_jetForMetJerSF,
                            const VLV* const v_associatedGenJetForMet)const;
    
    
    
private:
    
    /// Scale jet according to systematic variation
    void scaleJet(LV& jet, const double& jetJerSF, const LV& associatedGenJet)const;
    
    /// Bin number of eta interval for given jet
    int jetEtaBin(const LV& jet)const;
    
    
    
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
    
    
    
    /// Scale the jet collection
    void applyJetSystematic(VLV* v_jet)const;
    
    /// Scale the MET collection
    void applyMetSystematic(VLV* v_jetForMet, LV* met)const;
    
    
    
private:
    
    /// Scale jet according to systematic variation
    void scaleJet(LV& jet)const;
    
    
    
    /// Enumeration for possible systematics
    enum SystematicInternal{vary_up, vary_down, undefined};
    
    
    
    /// Object for retrieving uncertainty values
    ztop::JetCorrectionUncertainty* jetCorrectionUncertainty_;
    
    /// Variation upwards?
    bool varyUp_;
};




/// Class for reweighting the top pt to the measured spectrum
class TopPtScaleFactors{
    
public:
    
    /// Constructor
    TopPtScaleFactors(const Systematic::Systematic& systematic);
    
    /// Destructor
    ~TopPtScaleFactors();
    
    
    
    /// Get top-pt per-event scale factor
    double getSF(const double& topPt, const double& antiTopPt)const;
    
private:
    
    /// Enumeration for possible systematics
    enum SystematicInternal{nominal, vary_up, vary_down, undefined};
    
    /// Define the reweighting function for given systematic
    void reweightFunction(const SystematicInternal& systematic);
    
    /// Get top-pt scale factor per (anti-)top
    std::function<double(const double&)> reweightValue_;
};







#endif





