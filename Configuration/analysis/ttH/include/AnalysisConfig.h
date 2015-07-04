#ifndef AnalysisConfig_h
#define AnalysisConfig_h

#include <string>

#include "../../common/include/sampleHelpers.h"





class AnalysisConfig{
    
public:
    
    /// Struct for general info
    struct General{
        /// Constructor
        General();
        
        /// Print variables, return as string, and optionally print to screen
        std::string print(const bool screen =false)const;
        
        /// Analysis era
        Era::Era era_;
        
        /// Data luminosity in pb-1
        double luminosity_;
        
        /// Relative luminosity uncertainty
        double luminosityUncertainty_;
    };
    
    
    
    /// Struct for object/event corrections
    struct Corrections{
        /// Constructor
        Corrections();
        
        /// Print variables, return as string, and optionally print to screen
        std::string print(const bool screen =false)const;
        
        /// Pileup distribution file corresponding to data sample in use
        /// The file ending is automatically adjusted for different systematics
        std::string pileupInputFile_;
        
        /// Pileup MC era
        std::string pileupMcEra_;
        
        /// Pileup scenario
        std::string pileupScenario_;
        
        /// File ending of dilepton trigger scale factors input file
        std::string triggerSFInputSuffix_;
        
        /// Input file for electron ID scale factor
        std::string electronSFInputFile_;
        
        /// Input file for muon ID scale factor
        std::string muonSFInputFile_;
        
        /// Source for the uncertainties associated to JER
        std::string jerUncertaintySourceName_;
        
        /// File containing the uncertainties associated to JES
        std::string jesUncertaintySourceFile_;
        
        /// Correction mode for the b-tagging
        Btag::CorrectionMode btagCorrectionMode_;
        
        /// File for the official heavy flavour scale factors for b-tag discriminator reweighting
        std::string btagHeavyFlavourFile_;
        
        /// File for the official light flavour scale factors for b-tag discriminator reweighting
        std::string btagLightFlavourFile_;
        
        /// File containing the fits for the MVA MET recoil corrections in data
        std::string mvaMetRecoilDataFile_;
        
        /// File containing the fits for the MVA MET recoil corrections in MC
        std::string mvaMetRecoilMcFile_;
    };
    
    
    
    /// Struct for object selections
    struct Selections{
        /// Constructor
        Selections();
        
        /// Print variables, return as string, and optionally print to screen
        std::string print(const bool screen =false)const;
        
        /// Lepton eta selection (absolute value)
        double leptonEtaCut_;
        
        /// Lepton pt selection in GeV
        double leptonPtCut_;
        
        /// Jet eta selection (absolute value)
        double jetEtaCut_;
        
        /// Jet pt selection in GeV
        double jetPtCut_;
        
        /// Minimal deltaR for removal of jets that are close to leptons (if negative, no cut applied)
        double deltaRLeptonJetCut_;
        
        /// Leading 2 jet pt selection in GeV
        double lead2JetPtCut_;
        
        /// B-tag algorithm
        Btag::Algorithm btagAlgorithm_;
        
        /// B-tag working point
        Btag::WorkingPoint btagWorkingPoint_;
        
        /// PF MET or MVA MET
        bool mvaMet_;
        
        /// MET selection for same-flavour channels (ee, mumu)
        double metCut_;
        
        /// Generated jet eta selection (absolute value)
        double genJetEtaCut_;
        
        /// Generated jet pt selection in GeV
        double genJetPtCut_;
        
        /// Minimal deltaR for removal of generated jets that are close to leptons (if negative, no cut applied)
        double genDeltaRLeptonJetCut_;
    };
    
    
    
    /// Struct for sample compostion
    struct SampleComposition{
        /// Constructor
        SampleComposition();
        
        /// Print variables, return as string, and optionally print to screen
        std::string print(const bool screen =false)const;
        
        /// Whether to use pseudodata and how
        /// 0: use real data
        /// 1: use pseudodata as stacksum
        int pseudodata_;
        
        /// Level of merging for different samples
        /// 0: no merging, use all defined processes individually
        /// 1-x: different merging schemes
        int mergeLevel_;
    };
    
    
    
    /// Struct for plot style
    struct PlotStyle{
        /// Constructor
        PlotStyle();
        
        /// Print variables, return as string, and optionally print to screen
        std::string print(const bool screen =false)const;
        
        /// CMS label
        // FIXME: put here definitions of values
        int cmsLabel_;
    };
    
    
    
    /// Constructor reading in config file
    AnalysisConfig(const std::string& configfilename ="config.txt");
    
    /// Destructor
    ~AnalysisConfig(){}
    
    
    
    /// Return constant reference of struct general
    const General& general()const{return general_;}
    
    /// Return constant reference of struct object/event corrections
    const Corrections& corrections()const{return corrections_;}
    
    /// Return constant reference of struct object selections
    const Selections& selections()const{return selections_;}
    
    /// Return constant reference of struct sample compostion
    const SampleComposition& sampleComposition()const{return sampleComposition_;}
    
    /// Return constant reference of struct plot style
    const PlotStyle& plotStyle()const{return plotStyle_;}
    
    /// Print variables for all structs, return as string, and optionally print to screen
    std::string print(const bool screen =false)const;
    
    
    
private:
    
    /// Struct holding general info
    General general_;
    
    /// Struct holding object/event corrections
    Corrections corrections_;
    
    /// Struct holding object selections
    Selections selections_;
    
    /// Struct holding sample compostion
    SampleComposition sampleComposition_;
    
    /// Struct holding plot style
    PlotStyle plotStyle_;
};



#endif





