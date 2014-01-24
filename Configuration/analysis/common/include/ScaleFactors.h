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
#include "TopAnalysis/ZTopUtils/interface/bTagBase.h"




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









class BtagScaleFactors : public ztop::bTagBase{

public:

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
    bool makeEfficiencies();



    /// Fill histograms needed for b-tag efficiencies
    void fillBtagHistograms(const std::vector<int>& jetIndices,
            const std::vector<double>& bTagDiscriminant,
            const VLV& jets,
            const std::vector<int>& jetPartonFlavours,
            const double weight);

    /// Produce b-tag efficiencies
    void produceBtagEfficiencies(const std::string& channel);



    /// Get b-tag per-event scale factor
    double calculateBtagSF(const std::vector<int>& jetIndices,
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

    /// Prepare b-tagging scale factors (efficiency histograms)
    void prepareBTags(TSelectorList* output, const std::string& channel);



private:

    /// Input directory for files holding histograms of b-tagging efficiencies
    std::string inputDirName_;

    /// Output directory for files holding histograms, in case b-tagging efficiencies are produced in Analysis
    std::string outputDirName_;

    /// Name of the file with btag histograms
    std::string fileName_;



    /// Whether to produce btag efficiencies or whether to read in already produced ones
    //Jan  bool makeEfficiencies_;

    /// Pointer for bookkeeping of histograms
    // this is NOT used right now!
    // but it may be good for parallelization in the future
    TSelectorList* selectorList_;

    /// Map of the file names for each channel
    std::map<std::string, std::string> channelFileNames_;

    /// Map of the sample names for each channel
    std::map<std::string, std::string> channelSampleNames_;

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





