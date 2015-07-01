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
        
        /// Print variables to screen
        void print()const;
        
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
        
        /// Print variables to screen
        void print()const;
        
        /// Set pileup distribution file corresponding to data sample in use
        /// The file ending is automatically adjusted for different systematics
        std::string pileupInputFile_;
        
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
        
        /// Print variables to screen
        void print()const;
        
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
    
    
    
    /// Constructor reading in config file
    AnalysisConfig(const std::string& configfilename ="config.txt");
    
    /// Destructor
    ~AnalysisConfig(){}
    
    
    
    /// Return constant reference of struct general
    const General& general()const{return general_;}
    
    /// Return constant reference of struct general
    const Corrections& corrections()const{return corrections_;}
    
    /// Return constant reference of struct general
    const Selections& selections()const{return selections_;}
    
    /// Print variables to screen for all structs
    void print()const;
    
    
    
private:
    
    /// Struct holding general info
    General general_;
    
    /// Struct holding object/event corrections
    Corrections corrections_;
    
    /// Struct holding object selections
    Selections selections_;
};



#endif





