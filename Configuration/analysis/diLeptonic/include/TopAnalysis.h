#ifndef TopAnalysis_h
#define TopAnalysis_h

#include <map>

class TH1;
class TH2;

#include "../../common/include/AnalysisBase.h"
#include "../../common/include/classesFwd.h"





class TopAnalysis : public AnalysisBase
{
    
    /// Histograms
    TH2 *h_GenRecoLeptonpT,*h_GenRecoAntiLeptonpT,*h_GenRecoLeptonEta,*h_GenRecoAntiLeptonEta, *h_GenRecoLLBarMass, *h_GenRecoLLBarpT;
    TH2 *h_GenRecoBJetpT,*h_GenRecoAntiBJetpT, *h_GenRecoBJetEta,*h_GenRecoAntiBJetEta, *h_GenRecoBJetRapidity, *h_GenRecoAntiBJetRapidity;
    TH2 *h_GenRecoToppT,*h_GenRecoAntiToppT,*h_GenRecoTopRapidity,*h_GenRecoAntiTopRapidity, *h_GenRecoTTBarMass, *h_GenRecoTTBarpT, *h_GenRecoTTBarRapidity;
    TH2 *h_GenRecoMet;
    
    TH1 *h_NJetMatching;

    TH1 *TTh1_postZcut, *Allh1_postZcut, *Zh1_postZcut;
    TH1 *TTh1_post2jets, *Allh1_post2jets, *Zh1_post2jets;
    TH1 *TTh1_postMET, *Allh1_postMET, *Zh1_postMET;
    TH1 *TTh1_post1btag, *Allh1_post1btag, *Zh1_post1btag;
    TH1 *TTh1_postKinReco, *Allh1_postKinReco, *Zh1_postKinReco;
    
    TH1 *Looseh1, *Allh1, *Zh1, *TTh1, *h_diLepMassFull, *h_diLepMassFull_fullSel,
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

    TH1 *h_vertMulti, *h_vertMulti_noPU, *h_MET;

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

    TH1 *h_step3, *h_step4, *h_step5,*h_step6,*h_step7,*h_step8,*h_step9;

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
    TH2 *h_GenRecoDeltaPhiExtraJet12, *h_GenRecoPhiExtraJet12,*h_GenRecoTTBar1stJetMass, *h_GenRecoTTBar0Mass;
    TH1 *h_VisGenDeltaPhiExtraJet12, *h_RecoDeltaPhiExtraJet12, *h_HypDeltaPhiExtraJet12, *h_VisGenPhiExtraJet12, *h_RecoPhiExtraJet12, *h_HypPhiExtraJet12;
    TH1 *h_VisGenTTBar1stJetMass, *h_RecoTTBar1stJetMass, *h_HypTTBar1stJetMass;
    TH1 *h_VisGenTTBar0Mass, *h_RecoTTBar0Mass, *h_HypTTBar0Mass;
    //End: Plots for Carmen

    /// Plots for the parton momentum fraction defined by Olaf
    TH1 *h_HypPartonFraction, *h_HypAntiPartonFraction;
    TH1 *h_VisGenPartonFraction, *h_VisGenAntiPartonFraction;
    TH1 *h_RecoPartonFraction, *h_RecoAntiPartonFraction;
    TH2 *h_GenRecoPartonFraction, *h_GenRecoAntiPartonFraction;

    /// Histograms for event weights due to specific scale factor
    TH1 *h_PUSF, *h_TrigSF, *h_LepSF, *h_BTagSF, *h_KinRecoSF, *h_EventWeight;

    /// Control Plots
    TH1 *h_AllLeptonEta_step2, *h_AllLeptonpT_step2, *h_AllJetsEta_step2, *h_AllJetspT_step2;
    TH1 *h_AllLeptonEta_step3, *h_AllLeptonpT_step3, *h_AllJetsEta_step3, *h_AllJetspT_step3;
    TH1 *h_AllLeptonEta_step4, *h_AllLeptonpT_step4, *h_AllJetsEta_step4, *h_AllJetspT_step4;
    TH1 *h_AllLeptonEta_step5, *h_AllLeptonpT_step5, *h_AllJetsEta_step5, *h_AllJetspT_step5;
    TH1 *h_AllLeptonEta_step6, *h_AllLeptonpT_step6, *h_AllJetsEta_step6, *h_AllJetspT_step6;
    TH1 *h_AllLeptonEta_step7, *h_AllLeptonpT_step7, *h_AllJetsEta_step7, *h_AllJetspT_step7;
    TH1 *h_AllLeptonEta_step8, *h_AllLeptonpT_step8, *h_AllJetsEta_step8, *h_AllJetspT_step8;
    TH1 *h_AllLeptonEta_step9, *h_AllLeptonpT_step9, *h_AllJetsEta_step9, *h_AllJetspT_step9;
    TH1 *h_LeptonEta_step2, *h_LeptonpT_step2, *h_JetsEta_step2, *h_JetspT_step2;
    TH1 *h_LeptonEta_step3, *h_LeptonpT_step3, *h_JetsEta_step3, *h_JetspT_step3;
    TH1 *h_LeptonEta_step4, *h_LeptonpT_step4, *h_JetsEta_step4, *h_JetspT_step4;
    TH1 *h_LeptonEta_step5, *h_LeptonpT_step5, *h_JetsEta_step5, *h_JetspT_step5;
    TH1 *h_LeptonEta_step6, *h_LeptonpT_step6, *h_JetsEta_step6, *h_JetspT_step6;
    TH1 *h_LeptonEta_step7, *h_LeptonpT_step7, *h_JetsEta_step7, *h_JetspT_step7;
    TH1 *h_LeptonEta_step8, *h_LeptonpT_step8, *h_JetsEta_step8, *h_JetspT_step8;
    TH1 *h_LeptonEta_step9, *h_LeptonpT_step9, *h_JetsEta_step9, *h_JetspT_step9;
    TH1 *h_LeptonMult_step2, *h_JetsMult_step2, *h_BJetsMult_step2;
    TH1 *h_LeptonMult_step3, *h_JetsMult_step3, *h_BJetsMult_step3;
    TH1 *h_LeptonMult_step4, *h_JetsMult_step4, *h_BJetsMult_step4;
    TH1 *h_LeptonMult_step5, *h_JetsMult_step5, *h_BJetsMult_step5;
    TH1 *h_LeptonMult_step6, *h_JetsMult_step6, *h_BJetsMult_step6;
    TH1 *h_LeptonMult_step7, *h_JetsMult_step7, *h_BJetsMult_step7;
    TH1 *h_LeptonMult_step8, *h_JetsMult_step8, *h_BJetsMult_step8;
    TH1 *h_LeptonMult_step9, *h_JetsMult_step9, *h_BJetsMult_step9;
    
    /// Plots fot KinReco, Ievgen-11.11.2013
       TH1 *h_signalTopEvents_vs_JetMult; 
       TH1 *h_MatchedJets_vs_JetMult; 
       TH1 *h_nSolTtJets_vs_JetMult;
       TH1 *h_nSolCorrSignJets_vs_JetMult;
       
       TH1 *h_nRecoEvt_vs_JetMult;
       TH1 *h_nKinRecoSol_vs_JetMult;
       
       TH1 *h_nRecoEvt_vs_LepEta;
       TH1 *h_nKinRecoSol_vs_LepEta;

       TH1 *h_nRecoEvt_vs_JetEta;
       TH1 *h_nKinRecoSol_vs_JetEta;

       TH1 *h_nRecoEvt_vs_LeppT;  
       TH1 *h_nKinRecoSol_vs_LeppT;

       TH1 *h_nRecoEvt_vs_MET;
       TH1 *h_nKinRecoSol_vs_MET;
       
       TH1 *h_nRecoEvt_Eff;
       TH1 *h_nKinRecoSol_Eff;
       
       
       TH1 *h_RMSvsGenToppT;
       TH1 *h_RMSvsGenTopRapidity;
       TH1 *h_RMSvsGenTTBarMass; 
       
       TH1 *h_true_vE;
       TH1 *h_reco_vE;
    /// ... 
    
    
    /// Do kinematic reconstruction on nTuple level
    bool kinRecoOnTheFly_;
    
    
    /// Histogram for total weight of closure test
    TH1 *h_ClosureTotalWeight;
    
    /// Histogram for total weight of PDF variations
    TH1 *h_PDFTotalWeight;
    
    /// Whether to apply closure test
    bool doClosureTest_;
    
    /// Variables added from the outside
    int pdf_no_;
    
    
    /// Data for closure test
#ifndef __CINT__
    std::function<double(Long64_t)> closureFunction_;
#endif
    int closureMaxEvents_;
    
    
    
public:
    
    /// Constructor
    TopAnalysis():
        kinRecoOnTheFly_(false),
        doClosureTest_(false),
        pdf_no_(-1),
        closureFunction_(nullptr),
        binnedControlPlots_(0)
        {}
    
    /// Inherited from AnalysisBase and overwritten for needs of TopAnalysis
    virtual void Begin(TTree*);
    virtual void SlaveBegin(TTree*);
    virtual void SlaveTerminate();
    virtual void Terminate();
    virtual Bool_t Process(Long64_t entry);
    
    /// Set up closure test
    void SetClosureTest(TString closure, double slope);
    
    /// Set PDF variation for PDF systematics
    void SetPDF(int pdf_no);
    
    /// Class definition
    ClassDef(TopAnalysis, 0);    
    
    
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
    
    
    
    /// Set global normalisation factors
    double overallGlobalNormalisationFactor();
    
    /// Set global normalisation factor due to closure test
    double globalNormalisationFactorClosureTest();
    
    /// Set global normalisation factor due to PDF systematics
    double globalNormalisationFactorPDF();

    /// Get weight of closure test
    double calculateClosureTestWeight(const Long64_t& entry);
    
    /// Get indices of B hadron and anti-B hadron
    void bHadronIndices(int& bHadronIndex, int& antiBHadronIndex,const CommonGenObjects& commonGenObjects,const TopGenObjects& topGenObjects);
    
    void generatorTopEvent(LV& leadGenTop, LV& nLeadGenTop,
                           LV& leadGenLepton, LV& nLeadGenLepton,
                           LV& leadGenBJet, LV& nLeadGenBJet,
                           double& genHT,
                           const int bHadronIndex, const int antiBHadronIndex,
                           const double trueLevelWeightNoPileup, const double trueLevelWeight,
                           const CommonGenObjects& commonGenObjects,
                           const TopGenObjects& topGenObjects);

    void generatorTTbarjetsEvent(double& jetHTGen,
                                 const int bHadronIndex, const int antiBHadronIndex,
                                 const double trueLevelWeight,
                                 int& GenJets_cut, int& GenJets_cut40, int& GenJets_cut60, int& GenJets_cut100,
                                 double extragenjet[4],
                                 const CommonGenObjects& commonGenObjects,
                                 const TopGenObjects& topGenObjects);
    
    
    /// Map holding binned control plots
    //binnedControlPlots contains:
    //map of name of differential distribution
    // -> pair( histogram with the binning of the differential distribution,
    //          vector(bin) -> map( control plot name -> TH1*))
    std::map<std::string, std::pair<TH1*, std::vector<std::map<std::string, TH1*> > > >* binnedControlPlots_;
};


#endif


