#ifndef TopAnalysis_h
#define TopAnalysis_h

#include <map>

class TH1;
class TH2;

#include "../../common/include/AnalysisBase.h"
#include "../../common/include/classesFwd.h"
#include "analysisStructsFwd.h"

class AnalyzerBaseClass;
class TreeHandlerBase;
class EventMetadata;
class RecoObjects;
class CommonGenObjects;
class TopGenObjects;
class KinematicReconstructionSolutions;
namespace ttbar{
    class GenLevelWeights;
    class RecoLevelWeights;
    class GenObjectIndices;
    class RecoObjectIndices;
}




class TopAnalysis : public AnalysisBase
{
    
    /// Histograms
    TH2 *h_GenRecoLeptonpT,*h_GenRecoAntiLeptonpT,*h_GenRecoLeptonEta,*h_GenRecoAntiLeptonEta, *h_GenRecoLLBarMass, *h_GenRecoLLBarpT;
    TH2 *h_GenRecoBJetpT,*h_GenRecoAntiBJetpT, *h_GenRecoBJetEta,*h_GenRecoAntiBJetEta, *h_GenRecoBJetRapidity, *h_GenRecoAntiBJetRapidity;
    TH2 *h_GenRecoToppT,*h_GenRecoAntiToppT,*h_GenRecoTopRapidity,*h_GenRecoAntiTopRapidity, *h_GenRecoTTBarMass, *h_GenRecoTTBarpT, *h_GenRecoTTBarRapidity;
    TH2 *h_GenRecoMet;
    
    TH1 *h_NJetMatching;
    
    TH1 *h_diLepMassFull, *h_diLepMassFull_fullSel,
        *h_jetMultiXSec,*h_jetMultiAll, *h_jetMultiNoPU, *h_jetMultiVisTop,
        *h_jetMulti, *h_jetMulti_noBTag, *h_jetMulti_diLep, *h_BjetMulti, *h_BjetMulti_noBTag;

    TH1 *h_GenAll_RecoCuts, *h_GenAll_RecoCuts_noweight, *h_GenAll, *h_GenAll_noweight, *h_VisGenAll, *h_VisGenAll_noweight;

    TH1 *h_HypTTBarMass, *h_HypTTBarRapidity, *h_HypTTBarpT;
    TH1 *h_HypLLBarMass, *h_HypLLBarpT;
    TH1 *h_HypMet;

    TH1 *h_VisGenTTBarMass,*h_VisGenTTBarRapidity,*h_VisGenTTBarpT;
    TH1 *h_VisGenTopRapidity,*h_VisGenAntiTopRapidity;
    TH1 *h_VisGenLLBarMass,*h_VisGenLLBarpT;
    TH1 *h_VisGenMet;

    TH1 *h_RecoTTBarMass, *h_RecoTTBarRapidity,*h_RecoTTBarpT;
    TH1 *h_RecoToppT,*h_RecoAntiToppT,*h_RecoTopRapidity,*h_RecoAntiTopRapidity;
    TH1 *h_RecoLLBarMass, *h_RecoLLBarpT;
    TH1 *h_RecoLeptonpT,*h_RecoAntiLeptonpT,*h_RecoLeptonEta,*h_RecoAntiLeptonEta;
    TH1 *h_RecoBJetpT,*h_RecoAntiBJetpT, *h_RecoBJetRapidity,*h_RecoAntiBJetRapidity,*h_RecoBJetEta,*h_RecoAntiBJetEta;
    TH1 *h_RecoMet;
    
    TH2 *h_GenRecoHT;
    TH1 *h_VisGenHT, *h_HypHT, *h_RecoHT;

    TH1 *h_MET;

    TH1 *h_jetpT,*h_jetHT;
    TH1 *h_MuonpT, *h_MuonEta;
    TH1 *h_ElectronpT, *h_ElectronEta;
    TH1 *h_LeptonpT, *h_LeptonEta;
    TH1 *h_AntiLeptonpT, *h_AntiLeptonEta;
    TH1 *h_LeptonpT_diLep, *h_LeptonEta_diLep;
    TH1 *h_AntiLeptonpT_diLep, *h_AntiLeptonEta_diLep;

    TH1 *h_MuonpT_postMETcut, *h_MuonEta_postMETcut;
    TH1 *h_ElectronpT_postMETcut, *h_ElectronEta_postMETcut;
    TH1 *h_LeptonpT_postMETcut, *h_LeptonEta_postMETcut;
    TH1 *h_AntiLeptonpT_postMETcut, *h_AntiLeptonEta_postMETcut;
    
    TH1 *h_leptonPtBeforeKinReco, *h_leptonEtaBeforeKinReco;
    TH1 *h_leptonPtAfterKinReco, *h_leptonEtaAfterKinReco;
    TH1 *h_METBeforeKinReco, *h_METAfterKinReco;
    TH1 *h_bjetetaBeforeKinReco, *h_bjetetaAfterKinReco;
    
    TH1 *h_HypAntiToppT, *h_HypAntiTopEta, *h_HypAntiTopMass,*h_HypAntiTopRapidity;
    TH1 *h_HypToppT, *h_HypTopEta,*h_HypTopMass, *h_HypTopRapidity;
    
    TH1 *h_HypNeutrinopT, *h_HypAntiNeutrinopT;
    TH1 *h_RecoNeutrinopT, *h_RecoAntiNeutrinopT;
    TH1 *h_VisGenNeutrinopT, *h_VisGenAntiNeutrinopT;
    TH2 *h_GenRecoNeutrinopT, *h_GenRecoAntiNeutrinopT;
    
    TH1 *h_HypAntiBJetpT, *h_HypAntiBJetEta, *h_HypAntiBJetRapidity;
    TH1 *h_HypBJetpT, *h_HypBJetEta, *h_HypBJetRapidity;

    TH1 *h_HypAntiLeptonpT, *h_HypAntiLeptonEta;
    TH1 *h_HypLeptonpT, *h_HypLeptonEta;

    TH1 *h_VisGenAntiToppT, *h_VisGenAntiTopEta;
    TH1 *h_VisGenToppT, *h_VisGenTopEta;

    TH1 *h_VisGenAntiBJetpT, *h_VisGenAntiBJetEta, *h_VisGenAntiBJetRapidity;
    TH1 *h_VisGenBJetpT, *h_VisGenBJetEta, *h_VisGenBJetRapidity;

    TH1 *h_VisGenAntiBQuarkpT, *h_VisGenAntiBQuarkEta, *h_VisGenAntiBQuarkRapidity;
    TH1 *h_VisGenBQuarkpT, *h_VisGenBQuarkEta, *h_VisGenBQuarkRapidity;

    TH1 *h_VisGenAntiLeptonpT, *h_VisGenAntiLeptonEta;
    TH1 *h_VisGenLeptonpT, *h_VisGenLeptonEta;

    TH2 *h_GenRecoTTBarDeltaPhi, *h_GenRecoTTBarDeltaRapidity;
    TH1 *h_RecoTTBarDeltaPhi, *h_RecoTTBarDeltaRapidity;
    TH1 *h_HypTTBarDeltaPhi, *h_HypTTBarDeltaRapidity;
    TH1 *h_VisGenTTBarDeltaPhi, *h_VisGenTTBarDeltaRapidity;
    
    TH2 *h_GenRecoBBBarpT, *h_GenRecoBBBarMass;
    TH1 *h_RecoBBBarpT, *h_RecoBBBarMass;
    TH1 *h_HypBBBarpT, *h_HypBBBarMass;
    TH1 *h_VisGenBBBarpT, *h_VisGenBBBarMass;
    
    TH1 *h_HypToppTTTRestFrame, *h_HypAntiToppTTTRestFrame;
    TH1 *h_RecoToppTTTRestFrame, *h_RecoAntiToppTTTRestFrame;
    TH1 *h_VisGenToppTTTRestFrame, *h_VisGenAntiToppTTTRestFrame;
    TH2 *h_GenRecoToppTTTRestFrame, *h_GenRecoAntiToppTTTRestFrame;
    
    TH2 *h_GenRecoLLBarDPhi, *h_GenRecoLeptonantiBjetMass, *h_GenRecoAntiLeptonBjetMass, *h_GenRecoJetMult;
    TH1 *h_VisGenLLBarDPhi,  *h_VisGenLeptonantiBjetMass,  *h_VisGenAntiLeptonBjetMass,  *h_VisGenJetMult;
    TH1 *h_HypLLBarDPhi,     *h_HypLeptonantiBjetMass,     *h_HypAntiLeptonBjetMass,     *h_HypJetMult;
    TH1 *h_RecoLLBarDPhi,    *h_RecoLeptonantiBjetMass,    *h_RecoAntiLeptonBjetMass,    *h_RecoJetMult;

    TH1 *h_HypToppTLead,    *h_HypToppTNLead,    *h_HypTopRapidityLead, *h_HypTopRapidityNLead, *h_HypTopMassLead, *h_HypTopMassNLead;
    TH1 *h_HypLeptonpTLead, *h_HypLeptonpTNLead, *h_HypLeptonEtaLead,   *h_HypLeptonEtaNLead;
    TH1 *h_HypBJetpTLead,   *h_HypBJetpTNLead,   *h_HypBJetEtaLead,     *h_HypBJetEtaNLead;

    TH1 *h_RecoToppTLead,    *h_RecoToppTNLead,    *h_RecoTopRapidityLead, *h_RecoTopRapidityNLead, *h_RecoTopMassLead, *h_RecoTopMassNLead;
    TH1 *h_RecoLeptonpTLead, *h_RecoLeptonpTNLead, *h_RecoLeptonEtaLead,   *h_RecoLeptonEtaNLead;
    TH1 *h_RecoBJetpTLead,   *h_RecoBJetpTNLead,   *h_RecoBJetEtaLead,     *h_RecoBJetEtaNLead;

    TH1 *h_VisGenToppTLead,    *h_VisGenToppTNLead,    *h_VisGenTopRapidityLead, *h_VisGenTopRapidityNLead, *h_VisGenTopMassLead, *h_VisGenTopMassNLead;
    TH1 *h_VisGenLeptonpTLead, *h_VisGenLeptonpTNLead, *h_VisGenLeptonEtaLead,   *h_VisGenLeptonEtaNLead;
    TH1 *h_VisGenBJetpTLead,   *h_VisGenBJetpTNLead,   *h_VisGenBJetEtaLead,     *h_VisGenBJetEtaNLead;

    TH2 *h_GenRecoToppTLead,    *h_GenRecoToppTNLead,    *h_GenRecoTopRapidityLead, *h_GenRecoTopRapidityNLead, *h_GenRecoTopMassLead, *h_GenRecoTopMassNLead;
    TH2 *h_GenRecoLeptonpTLead, *h_GenRecoLeptonpTNLead, *h_GenRecoLeptonEtaLead,   *h_GenRecoLeptonEtaNLead;
    TH2 *h_GenRecoBJetpTLead,   *h_GenRecoBJetpTNLead,   *h_GenRecoBJetEtaLead,     *h_GenRecoBJetEtaNLead;

     //Begin: Plots for Carmen
    TH1 *h_RecoLeadingJetpT,    *h_RecoNLeadingJetpT,    *h_RecoLeadingJetEta,    *h_RecoNLeadingJetEta;
    TH1 *h_HypLeadingJetpT,     *h_HypNLeadingJetpT,     *h_HypLeadingJetEta,     *h_HypNLeadingJetEta;
    TH2 *h_GenRecoLeadingJetpT, *h_GenRecoLeadingJetEta, *h_GenRecoNLeadingJetpT, *h_GenRecoNLeadingJetEta;
    TH1 *h_VisGenLeadingJetpT,  *h_VisGenLeadingJetEta,  *h_VisGenNLeadingJetpT,  *h_VisGenNLeadingJetEta;

    TH1 *h_RecoExtraJetpT,  *h_HypExtraJetpT, *h_VisGenExtraJetpT, *h_RecoExtraJetEta, *h_HypExtraJetEta, *h_VisGenExtraJetEta;
    TH1 *h_RecoExtraJetpT2, *h_HypExtraJetpT2, *h_VisGenExtraJetpT2, *h_RecoExtraJetEta2, *h_HypExtraJetEta2, *h_VisGenExtraJetEta2;
    TH1 *h_RecoExtraJetpT3, *h_HypExtraJetpT3, *h_VisGenExtraJetpT3, *h_RecoExtraJetEta3, *h_HypExtraJetEta3, *h_VisGenExtraJetEta3;
    TH1 *h_RecoExtraJetpT4, *h_HypExtraJetpT4, *h_VisGenExtraJetpT4, *h_RecoExtraJetEta4, *h_HypExtraJetEta4, *h_VisGenExtraJetEta4;
    TH2 *h_GenRecoExtraJetpT, *h_GenRecoExtraJetEta, *h_GenRecoExtraJetpT2, *h_GenRecoExtraJetEta2, *h_GenRecoExtraJetpT3, *h_GenRecoExtraJetEta3, *h_GenRecoExtraJetpT4, *h_GenRecoExtraJetEta4;

    TH1 *h_RecoJetMultpt30, *h_RecoJetMultpt40, *h_HypJetMultpt30, *h_VisGenJetMultpt30, *h_HypJetMultpt40, *h_VisGenJetMultpt40, *h_RecoJetMultpt60, *h_HypJetMultpt60, *h_VisGenJetMultpt60;
    TH1 *h_RecoJetMultpt100, *h_HypJetMultpt100, *h_VisGenJetMultpt100;
    TH2 *h_GenRecoJetMultpt30, *h_GenRecoJetMultpt40, *h_GenRecoJetMultpt60, *h_GenRecoJetMultpt100;

    TH1 *h_HypJetMultQ0, *h_RecoJetMultQ0, *h_VisGenJetMultQ0;
    TH1 *h_RecoJetMultTotal, *h_HypJetMultTotal, *h_VisGenJetMultTotal;
    TH1 *h_HypJetMultQsum, *h_RecoJetMultQsum, *h_VisGenJetMultQsum;
    TH1 *h_HypJetExtra2Q0, *h_RecoJetExtra2Q0, *h_VisGenJetExtra2Q0;
    TH2 *h_GenRecoJetMultQ0, *h_GenRecoJetExtra2Q0, *h_GenRecoJetMultQsum, *h_GenRecoJetMultTotal;

    TH2 *h_GenRecoDeltaRExtraJet12;
    TH1 *h_VisGenDeltaRExtraJet12, *h_RecoDeltaRExtraJet12, *h_HypDeltaRExtraJet12;
    TH2 *h_GenRecoDeltaPhiExtraJet12, *h_GenRecoPhiExtraJet12,*h_GenRecoTTBar1stJetMass, *h_GenRecoTTBar0Mass, *h_GenRecoMassExtraJet12; 
    TH1 *h_VisGenDeltaPhiExtraJet12, *h_RecoDeltaPhiExtraJet12, *h_HypDeltaPhiExtraJet12, *h_VisGenPhiExtraJet12, *h_RecoPhiExtraJet12, *h_HypPhiExtraJet12;
    TH1 *h_VisGenMassExtraJet12, *h_RecoMassExtraJet12, *h_HypMassExtraJet12;

    TH1 *h_VisGenTTBar1stJetMass, *h_RecoTTBar1stJetMass, *h_HypTTBar1stJetMass;
    TH1 *h_VisGenTTBar0Mass, *h_RecoTTBar0Mass, *h_HypTTBar0Mass;
    //End: Plots for Carmen

    /// Plots for the parton momentum fraction defined by Olaf
    TH1 *h_HypTopPartonFraction, *h_HypAntiTopPartonFraction;
    TH1 *h_VisGenTopPartonFraction, *h_VisGenAntiTopPartonFraction;
    TH1 *h_RecoTopPartonFraction, *h_RecoAntiTopPartonFraction;
    TH2 *h_GenRecoTopPartonFraction, *h_GenRecoAntiTopPartonFraction;

    /// Histograms for event weights due to specific scale factor
    TH1 *h_PUSF, *h_TrigSF, *h_LepSF, *h_BTagSF, *h_KinRecoSF, *h_EventWeight;
    
    ///Ievgen
       
       TH1 *h_RMSvsGenToppT;
       TH1 *h_RMSvsGenTopRapidity;
       TH1 *h_RMSvsGenTTBarMass; 
       
       TH2 *h_HypTTBarRapidityvsTTBarpT;
       
       TH2 *h_VisGenTTBarRapidityvsTTBarpT;
    /// ... 
    
    /// Do kinematic reconstruction on nTuple level
    bool kinRecoOnTheFly_;
    
    
    /// Histogram for total weight of closure test
    TH1 *h_ClosureTotalWeight;
    
    /// Histogram for total weight of PDF variations
    TH1 *h_PDFTotalWeight;
    
    /// Whether to apply closure test
    bool doClosureTest_;
    
    /// Data for closure test
#ifndef __CINT__
    std::function<double(Long64_t)> closureFunction_;
#endif
    int closureMaxEvents_;
    
    /// Whether it is leptonic decays via tau in ttbar dilepton samples
    bool runViaTau_;
    
    
    
public:
    
    /// Constructor
    TopAnalysis();
    
    /// Inherited from AnalysisBase and overwritten for needs of TopAnalysis
    virtual void Begin(TTree*);
    virtual void SlaveBegin(TTree*);
    virtual void SlaveTerminate();
    virtual void Terminate();
    virtual Bool_t Process(Long64_t entry);
    
    /// Bool for separating direct dileptonic ttbar decays and decays via intermediate taus
    void SetRunViaTau(const bool runViaTau);
    
    /// Set up closure test
    void SetClosureTest(TString closure, double slope);
    
    /// Class definition
    ClassDef(TopAnalysis, 0);    
    
    /// Set up all analysers of type AnalyzerBaseClass
    void SetAllAnalyzers(std::vector<AnalyzerBaseClass*> v_analyzer);
    
    /// Set up all tree handlers of type TreeHandlerBase
    void SetAllTreeHandlers(std::vector<TreeHandlerBase*> v_treeHandler);
    
private:
    
    /// Create binned control plots
    // Create Nbins control plots for the differential distribution h_differential
    // Use h_control for the control plot name and binning
    void CreateBinnedControlPlots(TH1* h_differential, TH1* h_control, const bool fromHistoList =true);
    
    /// Fill binned control plots
    // h: differential distribution histogram
    // binvalue: the value of the quantity in the differential distribution histogram
    // the control plot histogram
    // the value for the control plot
    // weight: event weight
    void FillBinnedControlPlot(TH1* h_differential, double binvalue, 
                               TH1 *h_control, double value, double weight);
    
    
    
    /// Get weight of closure test
    double calculateClosureTestWeight(const Long64_t& entry);
    
    /// Get indices of B hadron and anti-B hadron
    void bHadronIndices(int& bHadronIndex, int& antiBHadronIndex, const TopGenObjects& topGenObjects);
    
    /// Select events from Top signal samples which need to be removed due to generator selection
    bool failsTopGeneratorSelection(const Long64_t& entry)const;
    
    void generatorTopEvent(LV& leadGenTop, LV& nLeadGenTop,
                           LV& leadGenLepton, LV& nLeadGenLepton,
                           LV& leadGenBJet, LV& nLeadGenBJet,
                           double& genHT,
                           const int bHadronIndex, const int antiBHadronIndex,
                           const double trueLevelWeightNoPileup, const double trueLevelWeight,
                           const TopGenObjects& topGenObjects);

    void generatorTTbarjetsEvent(double& jetHTGen,
                                 const int bHadronIndex, const int antiBHadronIndex,
                                 const double trueLevelWeight,
                                 int& GenJets_cut, int& GenJets_cut40, int& GenJets_cut60, int& GenJets_cut100, int& jetnum,
                                 double extragenjet[4],
                                 const TopGenObjects& topGenObjects);
    
    void generatorVisJets(const TopGenObjects& topGenObjects,std::vector<int>& genVisJetIndices);
    

    
    
    /// Map holding binned control plots
    //binnedControlPlots contains:
    //map of name of differential distribution
    // -> pair( histogram with the binning of the differential distribution,
    //          vector(bin) -> map( control plot name -> TH1*))
    std::map<std::string, std::pair<TH1*, std::vector<std::map<std::string, TH1*> > > >* binnedControlPlots_;
    
    /// Fill all analysers and histograms in one method
    void fillAll(const std::string& selectionStep,
                 const EventMetadata& eventMetadata,
                 const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                 const TopGenObjects& topGenObjects,
                 const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                 const ttbar::GenObjectIndices& genObjectIndices, const ttbar::RecoObjectIndices& recoObjectIndices,
                 const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                 const double& defaultWeight)const;
    
    /// Book all histograms of all analysers for all steps in one method
    void bookAll();
    
    /// Clear all analysers in one method
    void clearAll();
    
    /// All analysers of type AnalyzerBaseClass
    std::vector<AnalyzerBaseClass*> v_analyzer_;
    
        /// All tree handlers of type TreeHandlerBase
    std::vector<TreeHandlerBase*> v_treeHandler_;
    
};


#endif


