// -*- C++ -*-
//
// Package:    NTupleWriter_two
// Class:      NTupleWriter_two
//
/* *\class NTupleWriter_two NTupleWriter_two.cc TopAnalysis/NTupleWriter_two/src/NTupleWriter_two.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jan Kieseler,,,DESY
//         Created:  Thu Aug 11 16:37:05 CEST 2011
// $Id: NTupleWriter.cc,v 1.35 2013/04/02 14:26:31 hauk Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <map>
#include <boost/lexical_cast.hpp>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "AnalysisDataFormats/TopObjects/interface/TtFullLeptonicEvent.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h" //###############
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TopAnalysis/HiggsUtils/interface/HiggsGenEvent.h"
#include "TopAnalysis/HiggsUtils/interface/JetProperties.h"
#include "TopAnalysis/HiggsUtils/interface/GenZDecayProperties.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "Math/Vector3D.h"

#include <TTree.h>
#include <TLorentzVector.h>

//
// class declaration
//
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LV;

class NTupleWriter : public edm::EDAnalyzer
{
public:
    explicit NTupleWriter(const edm::ParameterSet&);
    ~NTupleWriter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    void clearVariables();
    int getTriggerBits(const edm::Event& iEvent, const edm::Handle< edm::TriggerResults >& trigResults);
    int getTriggerBits (const std::vector< std::string > &trigName);
    int getTriggerBitsTau(const edm::Event& iEvent, const edm::Handle< edm::TriggerResults >& trigResults);

    void AssignLeptonAndTau(const reco::GenParticle* lepton, LV &GenLepton, int &pdgid, LV &GenTau);
    bool isTau( const reco::GenParticle* lepton);
    const reco::GenParticle* getTauDaughter(const reco::GenParticle* tau);
    // ----------member data ---------------------------

    std::map<std::string, int> triggerMap_;
    std::map<std::string, int> triggerMapTau_;
    edm::InputTag inTag_PUSource;
    edm::InputTag elecs_, muons_;
    edm::InputTag jets_;
    edm::InputTag jetProperties_;
    edm::InputTag jetsForMET_;
    edm::InputTag jetsForMETuncorr_;
    edm::InputTag met_;
    edm::InputTag vertices_, genEvent_ ;
    edm::InputTag FullLepEvt_, hypoKey_;
    edm::InputTag genEventHiggs_;
    edm::InputTag genZDecay_;
    edm::InputTag genParticles_;
    edm::InputTag genJets_;
    edm::InputTag bIndex_;
    edm::InputTag antiBIndex_;
    edm::InputTag genMET_;
    edm::InputTag pdfWeightTag_;

    edm::InputTag bHadJetIdx_, antibHadJetIdx_;
    edm::InputTag BHadrons_, AntiBHadrons_;
    edm::InputTag BHadronFromTopB_, AntiBHadronFromTopB_;
    edm::InputTag BHadronVsJet_, AntiBHadronVsJet_;

    edm::InputTag genBHadPlusMothers_, genBHadPlusMothersIndices_;
    edm::InputTag genBHadIndex_, genBHadFlavour_, genBHadJetIndex_;
    edm::InputTag genBHadFromTopWeakDecay_, genBHadFromBHadron_;
    edm::InputTag genBHadLeptonIndex_, genBHadLeptonHadronIndex_;
    
    edm::InputTag genCHadPlusMothers_, genCHadPlusMothersIndices_;
    edm::InputTag genCHadIndex_, genCHadFlavour_, genCHadJetIndex_;
    edm::InputTag genCHadFromTopWeakDecay_, genCHadFromBHadron_;
    edm::InputTag genCHadLeptonIndex_, genCHadLeptonHadronIndex_;
    
    bool saveHadronMothers;

    bool includeTrig_;
    bool isTtBarSample_, isHiggsSample_, isZSample_;
    bool includePDFWeights_;
    bool isMC_;
    std::string sampleName_;
    std::string channelName_;
    std::string systematicsName_;
    double sampleCrossSection_;
    double dataLumi_;
    edm::InputTag trigResults_, decayMode_, higgsDecayMode_;

    bool includeZdecay_;
    edm::InputTag zDecayTag_;

    TTree* Ntuple;

    unsigned int runno;
    unsigned int lumibl;
    unsigned int eventno;
    unsigned int triggerBits;
    unsigned int triggerBitsTau;

    ////////dileptons and leptons/////
    std::vector<LV>     Vlep;
    std::vector<int>    VlepPdgId;
    std::vector<double> VlepID ; //mvaID for electrons (-1 for muon)
    std::vector<double> VlepChargedHadronIso;
    std::vector<double> VlepNeutralHadronIso;
    std::vector<double> VlepPhotonIso;
    std::vector<double> VlepPuChargedHadronIso;
    std::vector<double> VlepPfIso;
    std::vector<double> VlepCombIso;
    std::vector<double> VlepDxyVertex0;
    std::vector<double> VlepDzVertex0;
    std::vector<int>    VlepTrigger;

    //True level info from FullLepGenEvent
    LV GenTop, GenAntiTop;
    LV GenLepton, GenAntiLepton;
    int GenLeptonPdgId, GenAntiLeptonPdgId;
    LV GenTau, GenAntiTau;
    LV GenNeutrino, GenAntiNeutrino;
    LV GenB, GenAntiB;
    //  LV HadronGenB, HadronGenAntiB;
    LV GenWPlus, GenWMinus;
    LV GenMET;

    //True level info for Higgs and bbbar decay
    LV GenH;
    LV GenBFromH, GenAntiBFromH;

    std::vector<int>      VBHadJetIdx;
    std::vector<int>      VAntiBHadJetIdx;
    std::vector<LV>       VBHadron, VAntiBHadron;
    std::vector<bool>     VBHadFromTop, VAntiBHadFromTop;
    std::vector<int>      VBHadVsJet, VAntiBHadVsJet;

    int genExtraTopJetNumberId;
    std::vector<LV>              VgenBHadPlusMothers;
    std::vector<int>             VgenBHadPlusMothersPdg, VgenBHadPlusMothersStatus;
    std::vector<std::vector<int> >  VgenBHadPlusMothersIndices;
    std::vector<int>             VgenBHadIndex, VgenBHadFlavour, VgenBHadJetIndex;
    std::vector<int>             VgenBHadFromTopWeakDecay;
    std::vector<int>             VgenBHadLeptonHadronIndex, VgenBHadLeptonIndex, VgenBHadLeptonViaTau;
    
    std::vector<int>             VgenCHadFromBHadron, VgenCHadJetIndex;
    std::vector<int>             VgenCHadFromTopWeakDecay;
    std::vector<int>             VgenCHadLeptonHadronIndex, VgenCHadLeptonViaTau;

    // True level info from Zs and their decays
    std::vector<LV> VGenZ;
    std::vector<LV> VGenZMeDaughterParticle;
    std::vector<LV> VGenZMeDaughterAntiParticle;
    std::vector<LV> VGenZStableLepton;
    std::vector<LV> VGenZStableAntiLepton;
    std::vector<int> VGenZDecayMode;

    //Complete true level info
    std::vector<LV> GenParticleP4;
    std::vector<int> GenParticlePdgId;
    std::vector<int> GenParticleStatus;

    std::vector<LV> HypTop;
    std::vector<LV> HypAntiTop;
    std::vector<LV> HypLepton;
    std::vector<LV> HypAntiLepton;
    std::vector<LV> HypNeutrino;
    std::vector<LV> HypAntiNeutrino;
    std::vector<LV> HypB;
    std::vector<LV> HypAntiB;
    std::vector<LV> HypWPlus;
    std::vector<LV> HypWMinus;
    std::vector<int> HypJet0index;
    std::vector<int> HypJet1index;


    int TopProductionMode;
    int TopDecayMode;
    int HiggsDecayMode;
    std::vector<int> ZDecayMode;
    int recoInChannel;

    /////////jets///////////

    std::vector<LV> VallGenJets;
    std::vector<LV> VassociatedGenJet;
    std::vector<LV> VassociatedGenJetForMET;
    std::vector<LV> Vjet;
    std::vector<LV> VjetForMET;
    std::vector<double> VjetJERSF;
    std::vector<double> VjetForMETJERSF;
    std::vector<int> VjetPartonFlavourForMET;
    std::vector<int> VjetPartonFlavour;
    std::vector<double> VjetBTagTCHE;
    std::vector<double> VjetBTagTCHP;
    std::vector<double> VjetBTagJetProbability;
    std::vector<double> VjetBTagJetBProbability;
    std::vector<double> VjetBTagSSVHE;
    std::vector<double> VjetBTagSSVHP;
    std::vector<double> VjetBTagCSV;
    std::vector<double> VjetBTagCSVMVA;

    std::vector<double> VjetChargeGlobalPtWeighted;
    std::vector<double> VjetChargeRelativePtWeighted;
    std::vector<int> VjetAssociatedPartonPdgId;
    std::vector<LV> VjetAssociatedParton;
    std::vector<LV> VjetPfCandidateTrack;
    std::vector<int> VjetPfCandidateTrackCharge;
    std::vector<int> VjetPfCandidateTrackId;
    std::vector<int> VjetPfCandidateTrackIndex;
    std::vector<LV> VjetSelectedTrack;
    std::vector<double> VjetSelectedTrackIPValue;
    std::vector<double> VjetSelectedTrackIPSignificance;
    std::vector<int> VjetSelectedTrackCharge;
    std::vector<int> VjetSelectedTrackIndex;

    std::vector<double> VPdfWeights;

    /////////met///////////
    LV met;

    ////////triggers/////////
    std::vector<std::string> VfiredTriggers;

    ///////////weight//////////
    double weightGenerator;

    ////Vertices////
    int vertMulti;
    int vertMultiTrue;

    const LV nullP4;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//



//
// constructors and destructor
//
NTupleWriter::NTupleWriter(const edm::ParameterSet& iConfig):
    inTag_PUSource(iConfig.getParameter<edm::InputTag>("PUSource")),
    elecs_(iConfig.getParameter<edm::InputTag>("elecs")),
    muons_(iConfig.getParameter<edm::InputTag>("muons")),
    jets_(iConfig.getParameter<edm::InputTag>("jets")),
    jetProperties_(iConfig.getParameter<edm::InputTag>("jetProperties")),
    jetsForMET_(iConfig.getParameter<edm::InputTag>("jetsForMET")),
    jetsForMETuncorr_(iConfig.getParameter<edm::InputTag>("jetsForMETuncorr")),
    met_(iConfig.getParameter<edm::InputTag>("met")),

    vertices_(iConfig.getParameter<edm::InputTag>("vertices")),
    genEvent_(iConfig.getParameter<edm::InputTag>("src")),
    FullLepEvt_(iConfig.getParameter<edm::InputTag>("FullLepEvent")),
    hypoKey_(iConfig.getParameter<edm::InputTag>("hypoKey")),
    genEventHiggs_(iConfig.getParameter<edm::InputTag>("genEventHiggs")),
    genZDecay_(iConfig.getParameter<edm::InputTag>("genZDecay")),
    genParticles_( iConfig.getParameter<edm::InputTag>("genParticles")),

    genJets_   (iConfig.getParameter<edm::InputTag>("genJets")),
    bIndex_    (iConfig.getParameter<edm::InputTag>("BJetIndex")),
    antiBIndex_(iConfig.getParameter<edm::InputTag>("AntiBJetIndex")),

    genMET_(iConfig.getParameter<edm::InputTag>("genMET")),
    pdfWeightTag_ (iConfig.getParameter<edm::InputTag>("pdfWeights")),

    //Input from GenLevelBJetProducer
    bHadJetIdx_(iConfig.getParameter<edm::InputTag> ("BHadJetIndex")),
    antibHadJetIdx_(iConfig.getParameter<edm::InputTag> ("AntiBHadJetIndex")),
    BHadrons_(iConfig.getParameter<edm::InputTag> ("BHadrons")),
    AntiBHadrons_(iConfig.getParameter<edm::InputTag> ("AntiBHadrons")),
    BHadronFromTopB_ (iConfig.getParameter<edm::InputTag> ("BHadronFromTopB")),
    AntiBHadronFromTopB_(iConfig.getParameter<edm::InputTag> ("AntiBHadronFromTopB")),
    BHadronVsJet_(iConfig.getParameter<edm::InputTag> ("BHadronVsJet")),
    AntiBHadronVsJet_(iConfig.getParameter<edm::InputTag> ("AntiBHadronVsJet")),

    genBHadPlusMothers_(iConfig.getParameter<edm::InputTag> ("genBHadPlusMothers")),
    genBHadPlusMothersIndices_(iConfig.getParameter<edm::InputTag> ("genBHadPlusMothersIndices")),
    genBHadIndex_(iConfig.getParameter<edm::InputTag> ("genBHadIndex")),
    genBHadFlavour_(iConfig.getParameter<edm::InputTag> ("genBHadFlavour")),
    genBHadJetIndex_(iConfig.getParameter<edm::InputTag> ("genBHadJetIndex")),
    genBHadFromTopWeakDecay_(iConfig.getParameter<edm::InputTag> ("genBHadFromTopWeakDecay")),
    genBHadFromBHadron_(iConfig.getParameter<edm::InputTag> ("genBHadFromBHadron")),
    genBHadLeptonIndex_(iConfig.getParameter<edm::InputTag> ("genBHadLeptonIndex")),
    genBHadLeptonHadronIndex_(iConfig.getParameter<edm::InputTag> ("genBHadLeptonHadronIndex")),
    
    genCHadPlusMothers_(iConfig.getParameter<edm::InputTag> ("genCHadPlusMothers")),
    genCHadPlusMothersIndices_(iConfig.getParameter<edm::InputTag> ("genCHadPlusMothersIndices")),
    genCHadIndex_(iConfig.getParameter<edm::InputTag> ("genCHadIndex")),
    genCHadFlavour_(iConfig.getParameter<edm::InputTag> ("genCHadFlavour")),
    genCHadJetIndex_(iConfig.getParameter<edm::InputTag> ("genCHadJetIndex")),
    genCHadFromTopWeakDecay_(iConfig.getParameter<edm::InputTag> ("genCHadFromTopWeakDecay")),
    genCHadFromBHadron_(iConfig.getParameter<edm::InputTag> ("genCHadFromBHadron")),
    genCHadLeptonIndex_(iConfig.getParameter<edm::InputTag> ("genCHadLeptonIndex")),
    genCHadLeptonHadronIndex_(iConfig.getParameter<edm::InputTag> ("genCHadLeptonHadronIndex")),

    saveHadronMothers(iConfig.getParameter<bool>("saveHadronMothers")),

    includeTrig_(iConfig.getParameter<bool>("includeTrigger")),
    isTtBarSample_(iConfig.getParameter<bool>("isTtBarSample")),
    isHiggsSample_(iConfig.getParameter<bool>("isHiggsSample")),
    isZSample_(iConfig.getParameter<bool>("isZSample")),
    includePDFWeights_ (iConfig.getParameter<bool>("includePDFWeights")),

    //used for header:
    isMC_(iConfig.getParameter<bool>("isMC")),
    sampleName_(iConfig.getParameter<std::string>("sampleName")),
    channelName_(iConfig.getParameter<std::string>("channelName")),
    systematicsName_(iConfig.getParameter<std::string>("systematicsName")),

    trigResults_(iConfig.getParameter<edm::InputTag>("triggerResults")),
    decayMode_(iConfig.getParameter<edm::InputTag>("decayMode")),
    higgsDecayMode_(iConfig.getParameter<edm::InputTag>("higgsDecayMode")),

    includeZdecay_(iConfig.getParameter<bool>("includeZdecay")),
    zDecayTag_(iConfig.getParameter<edm::InputTag>("Zdecay")),


    nullP4(0, 0, 0, 0)
{

    //WARNING: The trigger map can be used either for a specific version, e.g. Trig_v6
    //or for any version, Trig_v*. NOT supported: Trig_v1* - the star captures ALL digits!
    //use first 8 bits for mumu
    triggerMap_["HLT_DoubleMu6_v*"] = 1;
    triggerMap_["HLT_DoubleMu7_v*"] = 2;
    triggerMap_["HLT_Mu13_Mu8_v*"] = 4;
    triggerMap_["HLT_Mu17_Mu8_v*"] = 8;
    triggerMap_["HLT_DoubleMu45_v*"] = 0x10;
    triggerMap_["HLT_Mu17_TkMu8_v*"] = 0x20;

    //use bits 9 to 16 for mu e
    triggerMap_["HLT_Mu8_Ele17_CaloIdL_v*"] = 0x100;
    triggerMap_["HLT_Mu17_Ele8_CaloIdL_v*"] = 0x200;
    triggerMap_["HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*"] = 0x400;
    triggerMap_["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*"] = 0x800;
    triggerMap_["HLT_Mu10_Ele10_CaloIdL_v*"] = 0x1000;
    triggerMap_["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"] = 0x2000;
    triggerMap_["HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"] = 0x4000;


    //use bits 17-24 for ee
    triggerMap_["HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*"] = 0x10000;
    triggerMap_["HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*"] = 0x20000;
    triggerMap_["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"] = 0x40000;
    triggerMap_["HLT_DoubleEle45_CaloIdL_v*"] = 0x80000;


    //use bit 32 for general trigger to avoid false for unmatched triggers
    triggerMap_["general"] = 0x80000000;



    // trigger map for hadronic tau (cross) triggers, analog to trigger map above
    //use first 8 bits for mu + tau
    triggerMapTau_["HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v*"] = 1;
    triggerMapTau_["HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v*"] = 2;

    //use bits 9 to 16 for e + tau
    triggerMapTau_["HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v*"] = 0x100;
    triggerMapTau_["HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*"] = 0x200;

    //use bits 17-24 for tau + tau + jet
    triggerMapTau_["HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30_v*"] = 0x10000;
    triggerMapTau_["HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30_v*"] = 0x20000;
    triggerMapTau_["HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30_v*"] = 0x40000;

    //use bit 32 for general trigger to avoid false for unmatched triggers - is this needed?
    //triggerMapTau_["general"] = 0x80000000;
}



NTupleWriter::~NTupleWriter()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}



//
// member functions
//

// ------------ method called for each event  ------------
void
NTupleWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
    clearVariables();
    //##################### MC weights for MCatNLO ###############
    if ( iEvent.isRealData() ) {
        weightGenerator = 1;
    } else {
        edm::Handle<GenEventInfoProduct> evt_info;
        iEvent.getByType(evt_info);
        weightGenerator = evt_info->weight();
    }

    //#############################################################

    if (includeZdecay_) {
        edm::Handle<std::vector<int> > zDecayModeHandle;
        iEvent.getByLabel(zDecayTag_, zDecayModeHandle);
        ZDecayMode = *zDecayModeHandle;
    }

    edm::Handle<int> TopDecayModeHandle;
    iEvent.getByLabel(decayMode_, TopDecayModeHandle);
    TopDecayMode = TopDecayModeHandle.failedToGet() ? 0 : *TopDecayModeHandle;

    //////////////////Trigger Stuff///////////////hltPaths_[i].c_str()
    edm::Handle<edm::TriggerResults> trigResults;
    iEvent.getByLabel(trigResults_, trigResults);
    triggerBits = getTriggerBits(iEvent, trigResults);
    triggerBitsTau = this->getTriggerBitsTau(iEvent, trigResults);

    //vertices and pu weights
    edm::Handle<std::vector<reco::Vertex> > vertices;
    iEvent.getByLabel(vertices_, vertices);
    vertMulti = vertices->size();

    // default values to allow for tracing errors
    if (! iEvent.isRealData()) {
        edm::Handle<edm::View<PileupSummaryInfo> > pPUInfo;
        iEvent.getByLabel(inTag_PUSource, pPUInfo);
        edm::View<PileupSummaryInfo>::const_iterator iterPU;
        for(iterPU = pPUInfo->begin(); iterPU != pPUInfo->end(); ++iterPU)  // vector size is 3
        {
            if (iterPU->getBunchCrossing() == 0) { // -1: previous BX, 0: current BX,  1: next BX
                vertMultiTrue = iterPU->getTrueNumInteractions();
            }
        }
    } else {
        vertMultiTrue = vertices->size();
    }

    //top event
    edm::Handle<TtFullLeptonicEvent> FullLepEvt;
    iEvent.getByLabel(FullLepEvt_, FullLepEvt);

    edm::Handle<int> hypoKeyHandle;
    iEvent.getByLabel(hypoKey_, hypoKeyHandle);

    edm::Handle<edm::View< pat::Jet > > jets;
    iEvent.getByLabel(jets_, jets);

    edm::Handle<std::vector<JetProperties> > jetPropertiesHandle;
    iEvent.getByLabel(jetProperties_, jetPropertiesHandle);

    edm::Handle<edm::View< pat::Jet > > jetsForMET;
    iEvent.getByLabel(jetsForMET_, jetsForMET);

    edm::Handle<edm::View< pat::Jet > > jetsForMETuncorr;
    iEvent.getByLabel(jetsForMETuncorr_, jetsForMETuncorr);


    if (! hypoKeyHandle.failedToGet())
    {
        TtEvent::HypoClassKey& hypoKey = ( TtEvent::HypoClassKey& ) *hypoKeyHandle;

        //////////////////////////////dilepton and lepton properties/////////////////////


        //the ntuple now only stores maximum 1 solution: the best one
        // - if a two btag solution is found, take it
        // - else if one tag solution is fount, take the best one
        if (FullLepEvt->isHypoAvailable(hypoKey) && FullLepEvt->isHypoValid(hypoKey))
        {
            int best = -1;
            for ( size_t i=0; i<FullLepEvt->numberOfAvailableHypos (hypoKey); ++i )
            {
                const pat::Jet &jet1 = jets->at(FullLepEvt->jetLeptonCombination(hypoKey,i)[0]);
                const pat::Jet &jet2 = jets->at(FullLepEvt->jetLeptonCombination(hypoKey,i)[1]);
                if (jet1.pt() <= 30 || jet2.pt() <= 30 || abs(jet1.eta()) > 2.4 || abs(jet2.eta()) > 2.4 ) continue;
                bool jet1tagged = jet1.bDiscriminator("combinedSecondaryVertexBJetTags")>0.244;
                bool jet2tagged = jet2.bDiscriminator("combinedSecondaryVertexBJetTags")>0.244;
                if (jet1tagged && jet2tagged) { best = i; break; }
                if ((jet1tagged || jet2tagged) && best < 0) { best = i; }
            }

            //std::cout << iEvent.eventAuxiliary().id() <<  ": choose combination #" << best << "\n";
            
            //for ( size_t i=0; i<FullLepEvt->numberOfAvailableHypos (hypoKey); ++i )
            if (best >= 0)
            {
                size_t i = best;

                const reco::Candidate* Top    = FullLepEvt->top(hypoKey, i);
                const reco::Candidate* TopBar = FullLepEvt->topBar(hypoKey, i);
                const reco::Candidate* Lep    = FullLepEvt->lepton(hypoKey, i);
                const reco::Candidate* LepBar = FullLepEvt->leptonBar(hypoKey, i);
                const reco::Candidate* Nu     = FullLepEvt->neutrino(hypoKey, i);
                const reco::Candidate* NuBar  = FullLepEvt->neutrinoBar(hypoKey, i);
                const reco::Candidate* B      = FullLepEvt->b(hypoKey, i);
                const reco::Candidate* BBar   = FullLepEvt->bBar(hypoKey, i);
                const reco::Candidate* Wplus  = FullLepEvt->wPlus(hypoKey, i);
                const reco::Candidate* Wminus = FullLepEvt->wMinus(hypoKey, i);

                HypTop.push_back(Top->polarP4());
                HypAntiTop.push_back( TopBar->polarP4());
                HypLepton.push_back( Lep->polarP4());
                HypAntiLepton.push_back( LepBar->polarP4());
                HypNeutrino.push_back( Nu->polarP4());
                HypAntiNeutrino.push_back( NuBar->polarP4());
                HypB.push_back( B->polarP4());
                HypAntiB.push_back( BBar->polarP4());
                HypWPlus.push_back( Wplus->polarP4());
                HypWMinus.push_back(Wminus->polarP4());

                HypJet0index.push_back(FullLepEvt->jetLeptonCombination(hypoKey,i)[0]);
                HypJet1index.push_back(FullLepEvt->jetLeptonCombination(hypoKey,i)[1]);
            }
        }
    }


    if ( isTtBarSample_ )
    {

        //Generator info
        edm::Handle<TtGenEvent> genEvt;
        iEvent.getByLabel(genEvent_, genEvt);
        if (! genEvt.failedToGet())
        {
	    // Decide which process generates the ttbar event: gluon-gluon-fusion (0), quark-quark-annihilation (1), all other processes (2)
            if (genEvt->fromGluonFusion()) TopProductionMode = 0;
	    else if (genEvt->fromQuarkAnnihilation()) TopProductionMode = 1;
	    else TopProductionMode = 2;

            GenTop = genEvt->top()->polarP4();
            GenAntiTop = genEvt->topBar()->polarP4();
            //for Higgs analysis, also non-leptonic modes might be written
            if (genEvt->lepton()) {
                AssignLeptonAndTau(genEvt->lepton(), GenLepton, GenLeptonPdgId, GenTau);
            } else {
                GenLepton = nullP4; GenLeptonPdgId = 0; GenTau = nullP4;
            }
            if (genEvt->leptonBar()) {
                AssignLeptonAndTau(genEvt->leptonBar(), GenAntiLepton, GenAntiLeptonPdgId, GenAntiTau);
            } else {
                GenAntiLepton = nullP4; GenAntiLeptonPdgId = 0; GenAntiTau = nullP4;
            }
            if (genEvt->b()) GenB = genEvt->b()->polarP4(); else GenB = nullP4;
            if (genEvt->bBar()) GenAntiB = genEvt->bBar()->polarP4(); else GenAntiB = nullP4;
            if (genEvt->neutrino()) GenNeutrino = genEvt->neutrino()->polarP4(); else GenNeutrino = nullP4;
            if (genEvt->neutrinoBar()) GenAntiNeutrino = genEvt->neutrinoBar()->polarP4(); else GenAntiNeutrino = nullP4;
            if (genEvt->wPlus()) GenWPlus = genEvt->wPlus()->polarP4(); else GenWPlus = nullP4;
            if (genEvt->wMinus()) GenWMinus = genEvt->wMinus()->polarP4(); else GenWMinus = nullP4;

            edm::Handle<std::vector<reco::GenMET> > h_genMET;
            iEvent.getByLabel(genMET_, h_genMET);
            GenMET = h_genMET->at(0).polarP4();

            edm::Handle<reco::GenJetCollection> genJets;
            iEvent.getByLabel(genJets_, genJets);

            for ( std::vector< reco::GenJet >::const_iterator aJet = genJets->begin(); aJet != genJets->end(); ++aJet) {
                VallGenJets.push_back(aJet->polarP4());
            }

            //find the b jets corresponding to b QUARKS (DeltaR selection)
            //         edm::Handle<int> BJetIndex;
            //         iEvent.getByLabel(bIndex_, BJetIndex);
            //         edm::Handle<int> AntiBJetIndex;
            //         iEvent.getByLabel(antiBIndex_, AntiBJetIndex);
            //    HadronGenAntiB = dummy;
            //    HadronGenB = dummy;
            //         if (*BJetIndex >= 0) (HadronGenB = genJets->at(*BJetIndex).polarP4());
            //         if (*AntiBJetIndex >= 0) (HadronGenAntiB = genJets->at(*AntiBJetIndex).polarP4());


            edm::Handle<std::vector<int> > BHadJetIndex;
            iEvent.getByLabel(bHadJetIdx_, BHadJetIndex);
            edm::Handle<std::vector<int> > AntiBHadJetIndex;
            iEvent.getByLabel(antibHadJetIdx_, AntiBHadJetIndex);

            for (size_t i = 0; i < BHadJetIndex->size(); ++i) {
                VBHadJetIdx.push_back(BHadJetIndex->at(i));
            }
            for (size_t i = 0; i < AntiBHadJetIndex->size(); ++i) {
                VAntiBHadJetIdx.push_back(AntiBHadJetIndex->at(i));
            }

            edm::Handle<std::vector<reco::GenParticle> > BHadrons;
            iEvent.getByLabel(BHadrons_, BHadrons);
            edm::Handle<std::vector<reco::GenParticle> > AntiBHadrons;
            iEvent.getByLabel(AntiBHadrons_, AntiBHadrons);


            for (std::vector<reco::GenParticle>::const_iterator it=BHadrons->begin(); it!=BHadrons->end(); ++it){
                VBHadron.push_back(it->polarP4());
            }
            for (std::vector<reco::GenParticle>::const_iterator it=AntiBHadrons->begin(); it!=AntiBHadrons->end(); ++it){
                VAntiBHadron.push_back(it->polarP4());
            }
            //    for (int i=0; i< (int) BHadrons->size(); ++i){  VBHadron.push_back(BHadrons->at(i).polarP4());}
            //    for (int i=0; i< (int) AntiBHadrons->size(); ++i){      VAntiBHadron.push_back(AntiBHadrons->at(i).polarP4());}


            edm::Handle<std::vector<bool> > BHadronFromTopB;
            iEvent.getByLabel(BHadronFromTopB_, BHadronFromTopB);
            edm::Handle<std::vector<bool> > AntiBHadronFromTopB;
            iEvent.getByLabel(AntiBHadronFromTopB_, AntiBHadronFromTopB);
            for (size_t i=0; i<BHadronFromTopB->size(); ++i) {
                VBHadFromTop.push_back(BHadronFromTopB->at(i));
            }
            for (size_t i=0; i<AntiBHadronFromTopB->size(); ++i) {
                VAntiBHadFromTop.push_back(AntiBHadronFromTopB->at(i));
            }

            edm::Handle<std::vector<int> > BHadronVsJet;
            iEvent.getByLabel(BHadronVsJet_, BHadronVsJet);
            edm::Handle<std::vector<int> > AntiBHadronVsJet;
            iEvent.getByLabel(AntiBHadronVsJet_, AntiBHadronVsJet);

            for (size_t i=0; i<BHadronVsJet->size(); ++i) {
                VBHadVsJet.push_back(BHadronVsJet->at(i));
            }
            for (size_t i=0; i<AntiBHadronVsJet->size(); ++i) {
                VAntiBHadVsJet.push_back(AntiBHadronVsJet->at(i));
            }


            edm::Handle<std::vector<int> > genBHadIndex;
            iEvent.getByLabel(genBHadIndex_, genBHadIndex);
            if(!genBHadIndex.failedToGet()) {
                for (std::vector<int>::const_iterator it=genBHadIndex->begin(); it!=genBHadIndex->end(); ++it) {
                    VgenBHadIndex.push_back(*it);
                }
            }
            
            edm::Handle<std::vector<std::vector<int> > > genBHadPlusMothersIndices;
            iEvent.getByLabel(genBHadPlusMothersIndices_, genBHadPlusMothersIndices);
            if(!genBHadPlusMothersIndices.failedToGet() && saveHadronMothers) {         // Only if all hadron mothers have to be stored
                for (std::vector<std::vector<int> >::const_iterator it=genBHadPlusMothersIndices->begin(); it!=genBHadPlusMothersIndices->end(); ++it) {
                    VgenBHadPlusMothersIndices.push_back(*it);
                }
            }
            
            edm::Handle<std::vector<int> > genBHadLeptonIndex;
            iEvent.getByLabel(genBHadLeptonIndex_, genBHadLeptonIndex);
            if(!genBHadLeptonIndex.failedToGet()) {
                for(std::vector<int>::const_iterator it=genBHadLeptonIndex->begin(); it!=genBHadLeptonIndex->end(); ++it) {
                    VgenBHadLeptonIndex.push_back(*it);
                }
            }

            edm::Handle<std::vector<reco::GenParticle> > genBHadPlusMothers;
            iEvent.getByLabel(genBHadPlusMothers_, genBHadPlusMothers);
            if(!genBHadPlusMothers.failedToGet()) {
                if(saveHadronMothers) {         // If all particles have to be stored
                    for (std::vector<reco::GenParticle>::const_iterator it=genBHadPlusMothers->begin(); it!=genBHadPlusMothers->end(); ++it){
                        VgenBHadPlusMothers.push_back(it->polarP4());
                        VgenBHadPlusMothersPdg.push_back(it->pdgId());
                        VgenBHadPlusMothersStatus.push_back(it->status());
                    }       // End of loop over all particles
                }   else {      // If only hadrons/leptons have to be stored
                    for (unsigned int i=0; i<genBHadIndex->size(); ++i) {
                        VgenBHadPlusMothers.push_back(genBHadPlusMothers->at(genBHadIndex->at(i)).polarP4());
                        VgenBHadPlusMothersPdg.push_back(genBHadPlusMothers->at(genBHadIndex->at(i)).pdgId());
                        VgenBHadIndex.at(i)=VgenBHadPlusMothers.size()-1;
                    }       // End of loop over hadrons
                    for (unsigned int i=0; i<genBHadLeptonIndex->size(); ++i) {
                        VgenBHadPlusMothers.push_back(genBHadPlusMothers->at(genBHadLeptonIndex->at(i)).polarP4());
                        VgenBHadPlusMothersPdg.push_back(genBHadPlusMothers->at(genBHadLeptonIndex->at(i)).pdgId());
                        VgenBHadLeptonIndex.at(i)=VgenBHadPlusMothers.size()-1;
                        int viaTau = -1;
                        if(!genBHadPlusMothersIndices.failedToGet()) {
                            int leptonMotherIndex = genBHadPlusMothersIndices->at(genBHadLeptonIndex->at(i)).at(0);
                            if(leptonMotherIndex>=0) {
                                viaTau = std::abs(genBHadPlusMothers->at(leptonMotherIndex).pdgId()) == 15?1:0;
                            }
                        }
                        VgenBHadLeptonViaTau.push_back(viaTau);
                        
                    }       // End of loop over leptons
                }
            }       // If genBHadPlusMothers is not empty
            
            edm::Handle<std::vector<int> > genBHadLeptonHadronIndex;
            iEvent.getByLabel(genBHadLeptonHadronIndex_, genBHadLeptonHadronIndex);
            if(!genBHadLeptonHadronIndex.failedToGet()) {
                for(std::vector<int>::const_iterator it=genBHadLeptonHadronIndex->begin(); it!=genBHadLeptonHadronIndex->end(); ++it) {
                    VgenBHadLeptonHadronIndex.push_back(*it);
                }
            }
            
            edm::Handle<std::vector<int> > genBHadFlavour;
            iEvent.getByLabel(genBHadFlavour_, genBHadFlavour);
            if(!genBHadFlavour.failedToGet()) {
                for (std::vector<int>::const_iterator it=genBHadFlavour->begin(); it!=genBHadFlavour->end(); ++it) {
                    VgenBHadFlavour.push_back(*it);
                }
            }
            std::vector<int> genBJetFromTopIds;
            edm::Handle<std::vector<int> > genBHadJetIndex;
            iEvent.getByLabel(genBHadJetIndex_, genBHadJetIndex);
            if(!genBHadJetIndex.failedToGet()) {
                unsigned int iHad=0;
                for (std::vector<int>::const_iterator it=genBHadJetIndex->begin(); it!=genBHadJetIndex->end(); ++it) {
                    VgenBHadJetIndex.push_back(*it);
                    // Building a vector of b-jets from top
                    if(VgenBHadFlavour.size()>iHad) {
                        if(std::abs(VgenBHadFlavour.at(iHad))==6) genBJetFromTopIds.push_back(*it);
                    }
                    iHad++;
                }
            }
            
            edm::Handle<std::vector<int> > genBHadFromTopWeakDecay;
            iEvent.getByLabel(genBHadFromTopWeakDecay_, genBHadFromTopWeakDecay);
            if(!genBHadFromTopWeakDecay.failedToGet()) {
                for (std::vector<int>::const_iterator it=genBHadFromTopWeakDecay->begin(); it!=genBHadFromTopWeakDecay->end(); ++it) {
                    VgenBHadFromTopWeakDecay.push_back(*it);
                }
            }
                        
            edm::Handle<std::vector<int> > genCHadLeptonHadronIndex;
            iEvent.getByLabel(genCHadLeptonHadronIndex_, genCHadLeptonHadronIndex);
            if(!genCHadLeptonHadronIndex.failedToGet()) {
                for(std::vector<int>::const_iterator it=genCHadLeptonHadronIndex->begin(); it!=genCHadLeptonHadronIndex->end(); ++it) {
                    VgenCHadLeptonHadronIndex.push_back(*it);
                }
            }
            
            edm::Handle<std::vector<reco::GenParticle> > genCHadPlusMothers;
            iEvent.getByLabel(genCHadPlusMothers_, genCHadPlusMothers);
            edm::Handle<std::vector<std::vector<int> > > genCHadPlusMothersIndices;
            iEvent.getByLabel(genCHadPlusMothersIndices_, genCHadPlusMothersIndices);
            edm::Handle<std::vector<int> > genCHadLeptonIndex;
            iEvent.getByLabel(genCHadLeptonIndex_, genCHadLeptonIndex);
            if(!genCHadPlusMothers.failedToGet() && !genCHadPlusMothersIndices.failedToGet() && !genCHadLeptonIndex.failedToGet()) {
                    for (unsigned int i=0; i<genCHadLeptonIndex->size(); ++i) {
                        int viaTau = -1;
                        if(!genCHadPlusMothersIndices.failedToGet()) {
                            int leptonMotherIndex = genCHadPlusMothersIndices->at(genCHadLeptonIndex->at(i)).at(0);
                            if(leptonMotherIndex>=0) {
                                viaTau = std::abs(genCHadPlusMothers->at(leptonMotherIndex).pdgId()) == 15?1:0;
                            }
                        }
                        VgenCHadLeptonViaTau.push_back(viaTau);
                    }       // End of loop over leptons
            }       // If genCHadPlusMothers is not empty

            
            edm::Handle<std::vector<int> > genCHadFromBHadron;
            iEvent.getByLabel(genCHadFromBHadron_, genCHadFromBHadron);
            if(!genCHadFromBHadron.failedToGet()) {
                for(std::vector<int>::const_iterator it=genCHadFromBHadron->begin(); it!=genCHadFromBHadron->end(); ++it) {
                    VgenCHadFromBHadron.push_back(*it);
                }
            }
            
            edm::Handle<std::vector<int> > genCHadJetIndex;
            iEvent.getByLabel(genCHadJetIndex_, genCHadJetIndex);
            if(!genCHadJetIndex.failedToGet()) {
                for (std::vector<int>::const_iterator it=genCHadJetIndex->begin(); it!=genCHadJetIndex->end(); ++it) {
                    VgenCHadJetIndex.push_back(*it);
                }
            }
            
            edm::Handle<std::vector<int> > genCHadFromTopWeakDecay;
            iEvent.getByLabel(genCHadFromTopWeakDecay_, genCHadFromTopWeakDecay);
            if(!genCHadFromTopWeakDecay.failedToGet()) {
                for (std::vector<int>::const_iterator it=genCHadFromTopWeakDecay->begin(); it!=genCHadFromTopWeakDecay->end(); ++it) {
                    VgenCHadFromTopWeakDecay.push_back(*it);
                }
            }

            // Id of the event depending on the number of additional b(c)-jets in the event
            // -1  : Something went wrong
            // 2   : 2 or more extra b-jets coming not from the top weak decay chain (and any number of c/light jets)
            // 1   : Exactly 1 extra b-jet coming not from the top weak decay chain (and any number of c/light jets)
            // 12  : 2 or more extra b-jets coming from the top weak decay chain but not directly from top (and any number of c/light jets)
            // 11  : Exactly 1 extra b-jet coming from the top weak decay chain but not directly from top (and any number of c/light jets)
            // 22  : No extra b-jets and 2 or more c-jets coming not from the top weak decay chain (and any number of light jets)
            // 21  : No extra b-jets and exactly 1 c-jet coming not from the top weak decay chain (and any number of light jets)
            // 32  : No extra b-jets and 2 or more c-jets coming from the top weak decay chain (and any number of light jets)
            // 31  : No extra b-jets and exactly 1 c-jet coming from the top weak decay chain (and any number of light jets)
            // 0   : No extra b-jets that come not directly from top and no c-jets that come not from b-hadrons
            genExtraTopJetNumberId = -1;
                
            // Signal definition cuts
            double genSignalJetPt_min = 20.0;
            double genSignalJetEta_max = 2.5;
            
            // Making the list of jets that contain b-hadrons not coming from the top weak decay
            std::vector<int> genBJetNotFromTopDecayChainIndex;
            std::vector<int> genBJetFromTopDecayChainIndex;
            if(!genBHadJetIndex.failedToGet()) {
                for(size_t iHad=0; iHad<genBHadJetIndex->size(); ++iHad) {
                    if(VgenBHadFromTopWeakDecay.size() <= iHad) break;
                    // Skipping b-hadrons from top
                    if(genBHadFlavour->size()>iHad && std::abs(genBHadFlavour->at(iHad))==6) continue;
                    int jetId = genBHadJetIndex->at(iHad);
                    if(jetId<0) continue;
                    // Checking that the jet satisfies the signal definition
                    if(VallGenJets.at(jetId).Pt()<genSignalJetPt_min) continue;
                    if(std::fabs(VallGenJets.at(jetId).Eta())>genSignalJetEta_max) continue;
                    // Skipping the jet if it is already identified as b-jet from top
                    if(std::find(genBJetFromTopIds.begin(), genBJetFromTopIds.end(), jetId) != genBJetFromTopIds.end()) continue;
                    // Checking that the jet is counted only once
                    if(std::find(genBJetNotFromTopDecayChainIndex.begin(), genBJetNotFromTopDecayChainIndex.end(), jetId) != genBJetNotFromTopDecayChainIndex.end()) continue;
                    if(VgenBHadFromTopWeakDecay.at(iHad)==0) {
                        genBJetNotFromTopDecayChainIndex.push_back(jetId);
                    } else if(VgenBHadFromTopWeakDecay.at(iHad)==1) {
                        genBJetFromTopDecayChainIndex.push_back(jetId);
                    }
                }       // End of loop over b-hadrons
            }       // If jet index for b-hadrons exists
            
            // Making the list of jets that contain c-hadrons not coming from the top weak decay
            std::vector<int> genCJetNotFromTopDecayChainIndex;
            std::vector<int> genCJetFromTopDecayChainIndex;
            if(!genCHadJetIndex.failedToGet()) {
                for(size_t iHad=0; iHad<genCHadJetIndex->size(); ++iHad) {
                    if(VgenCHadFromTopWeakDecay.size() <= iHad) break;
                    // Skipping c-hadrons from top
                    if(genBHadFlavour->size()>iHad && std::abs(genBHadFlavour->at(iHad))==6) continue;
                    // Skipping c-hadrons coming from b-hadrons
                    if(genCHadFromBHadron->size()>iHad && genCHadFromBHadron->at(iHad)) continue;
                    int jetId = genCHadJetIndex->at(iHad);
                    if(jetId<0) continue;
                    // Checking that the jet satisfies the signal definition
                    if(VallGenJets.at(jetId).Pt()<genSignalJetPt_min) continue;
                    if(std::fabs(VallGenJets.at(jetId).Eta())>genSignalJetEta_max) continue;
                    // Skipping the jet if it is already identified as b-jet from top
                    if(std::find(genBJetFromTopIds.begin(), genBJetFromTopIds.end(), jetId) != genBJetFromTopIds.end()) continue;
                    // Checking that the jet is counted only once
                    if(std::find(genCJetNotFromTopDecayChainIndex.begin(), genCJetNotFromTopDecayChainIndex.end(), jetId) != genCJetNotFromTopDecayChainIndex.end()) continue;
                    if(VgenCHadFromTopWeakDecay.at(iHad)==0) {
                        genCJetNotFromTopDecayChainIndex.push_back(jetId);
                    } else if(VgenCHadFromTopWeakDecay.at(iHad)==1) {
                        genCJetFromTopDecayChainIndex.push_back(jetId);
                    }
                }       // End of loop over c-hadrons
            }       // If jet index for c-hadrons exists
                
            int nGenBJetNotFromTopDecayChain = genBJetNotFromTopDecayChainIndex.size();
            int nGenBJetFromTopDecayChain = genBJetFromTopDecayChainIndex.size();
            int nGenCJetNotFromTopDecayChain = genCJetNotFromTopDecayChainIndex.size();
            int nGenCJetFromTopDecayChain = genCJetFromTopDecayChainIndex.size();
            if(nGenBJetNotFromTopDecayChain>=2) genExtraTopJetNumberId=2;
            else if(nGenBJetNotFromTopDecayChain==1) genExtraTopJetNumberId=1;
            else if(nGenBJetFromTopDecayChain>=2) genExtraTopJetNumberId = 12;
            else if(nGenBJetFromTopDecayChain==1) genExtraTopJetNumberId=11;
            else if(nGenCJetNotFromTopDecayChain>=2) genExtraTopJetNumberId=22;
            else if(nGenCJetNotFromTopDecayChain==1) genExtraTopJetNumberId=21;
            else if(nGenCJetFromTopDecayChain>=2) genExtraTopJetNumberId = 32;
            else if(nGenCJetFromTopDecayChain==1) genExtraTopJetNumberId=31;
            else genExtraTopJetNumberId = 0;
            
        }       // If GenEvent is available
        else
        {
            std::cerr << "Error: no gen event?!\n";
            TopProductionMode = -1;
            GenTop = nullP4; GenAntiTop = nullP4;
            GenLepton = nullP4; GenAntiLepton = nullP4;
            GenTau = nullP4; GenAntiTau = nullP4;
            GenLeptonPdgId = 0; GenAntiLeptonPdgId = 0;
            GenB = nullP4; GenAntiB = nullP4;
            GenNeutrino = nullP4; GenAntiNeutrino = nullP4;
            GenWMinus = nullP4; GenWPlus = nullP4;
            GenMET = nullP4;
        }

        //put more true info
        //         edm::Handle<std::vector<reco::GenParticle> > genParticles;
        //         iEvent.getByLabel(genParticles_, genParticles);
        //         for ( vector< reco::GenParticle >::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p) {
        //             GenParticleP4.push_back(p->polarP4());
        //             GenParticlePdgId.push_back(p->pdgId());
        //             GenParticleStatus.push_back(p->status());
        //         }
    }
    
        
    if ( isHiggsSample_ ) {
        //Generator info
        edm::Handle<HiggsGenEvent> genEvtHiggs;
        iEvent.getByLabel ( genEventHiggs_, genEvtHiggs );
        if (! genEvtHiggs.failedToGet())
        {
            GenH = genEvtHiggs->higgs()->polarP4();
            if(genEvtHiggs->b()) GenBFromH = genEvtHiggs->b()->polarP4(); else GenBFromH = nullP4;
            if(genEvtHiggs->bBar()) GenAntiBFromH = genEvtHiggs->bBar()->polarP4(); else GenAntiBFromH = nullP4;
        }
        else
        {
            std::cerr << "Error: no Higgs gen event?!\n";
            GenH = nullP4;
            GenBFromH = nullP4; GenAntiBFromH = nullP4;
        }

        edm::Handle<int> HiggsDecayModeHandle;
        iEvent.getByLabel(higgsDecayMode_, HiggsDecayModeHandle);
        HiggsDecayMode = HiggsDecayModeHandle.failedToGet() ? 0 : *HiggsDecayModeHandle;
    }

    if(isZSample_){
        //Generator info
        edm::Handle<std::vector<GenZDecayProperties> > v_genZDecayProperties;
        iEvent.getByLabel(genZDecay_, v_genZDecayProperties);
        if(!v_genZDecayProperties.failedToGet()){
            for(size_t i = 0; i < v_genZDecayProperties->size(); ++i){
                VGenZ.push_back(v_genZDecayProperties->at(i).z()->polarP4());
                VGenZMeDaughterParticle.push_back(v_genZDecayProperties->at(i).meDaughterParticle()->polarP4());
                VGenZMeDaughterAntiParticle.push_back(v_genZDecayProperties->at(i).meDaughterAntiParticle()->polarP4());
                if(v_genZDecayProperties->at(i).stableLepton()) VGenZStableLepton.push_back(v_genZDecayProperties->at(i).stableLepton()->polarP4());
                else VGenZStableLepton.push_back(nullP4);
                if(v_genZDecayProperties->at(i).stableAntiLepton()) VGenZStableAntiLepton.push_back(v_genZDecayProperties->at(i).stableAntiLepton()->polarP4());
                else VGenZStableAntiLepton.push_back(nullP4);
                VGenZDecayMode.push_back(v_genZDecayProperties->at(i).decayMode());
            }
        }
        else{
            std::cerr << "Error: no GenZDecayProperties ?!\n";
        }
    }


    //////fill pfiso///maybe other iso??
    edm::Handle<std::vector<pat::Muon> > muons;
    iEvent.getByLabel(muons_, muons);

    edm::Handle<std::vector<pat::Electron> > electrons;
    iEvent.getByLabel(elecs_, electrons);

    std::vector<pat::Muon>::const_iterator amuon  = muons->begin();
    std::vector<pat::Electron>::const_iterator anelectron  = electrons->begin();


    while ( amuon < muons->end() || anelectron < electrons->end() ) //sort input leptons by pt
    {
        bool writemuon=false;
        bool writeelec=false;

        if ( anelectron == electrons->end() )  writemuon=true;
        else if ( amuon == muons->end() )      writeelec=true;
        else if ( amuon->pt() > anelectron->pt() ) writemuon=true;
        else writeelec=true;


        if ( writemuon )
        {
            //Fill muonstuff
            Vlep.push_back( amuon->polarP4());
            VlepPdgId.push_back(amuon->pdgId());
            VlepID.push_back(-1);

            if (amuon->globalTrack().isAvailable()) {
                VlepDxyVertex0.push_back(amuon->globalTrack()->dxy(vertices->at(0).position()));
                VlepDzVertex0.push_back(amuon->globalTrack()->dz(vertices->at(0).position()));
            } else {
                VlepDxyVertex0.push_back(-1000); //no such value available
                VlepDzVertex0.push_back(-1000); //no such value available
            }

            VlepChargedHadronIso.push_back(amuon->chargedHadronIso());
            VlepNeutralHadronIso.push_back(amuon->neutralHadronIso());
            VlepPhotonIso.push_back(amuon->photonIso());
            VlepPuChargedHadronIso.push_back(amuon->puChargedHadronIso());
            VlepPfIso.push_back((amuon->chargedHadronIso()
                                 + std::max(0., amuon->neutralHadronIso()
                                               + amuon->photonIso()
                                               - 0.5*amuon->puChargedHadronIso())
                                ) / amuon->pt());

            VlepCombIso.push_back(( amuon->trackIso() +amuon->caloIso())/amuon->pt());


            int triggerResult = 0;

            //          std::cout << "processing muon, is lepton object: " << Vlep.size() << std::endl;

            pat::TriggerObjectStandAloneCollection triggers = amuon->triggerObjectMatches();
            if (triggers.size() > 0) triggerResult = triggerMap_["general"];
            for (unsigned int i = 0; i < triggers.size(); ++i)
                triggerResult |= getTriggerBits(triggers[i].pathNames());
            VlepTrigger.push_back(triggerResult);

            //          std::cout << "trigger result: " << std::hex << triggerResult << std::dec << std::endl;
            ++amuon;


        }
        if ( writeelec )
        {
            //Fill elestuff

            //Electron MVAID values

            std::vector<std::pair<std::string,float> > electronMVAIDs = anelectron->electronIDs();

            double idtemp = -9999.0;
            for(unsigned int id = 0; id < electronMVAIDs.size(); id++){
                if(electronMVAIDs[id].first == "mvaTrigV0"){
                    idtemp = electronMVAIDs[id].second;
                    break;
                }
            }

            VlepID.push_back(idtemp);
            const reco::GsfTrack &track = *(anelectron->gsfTrack());
            VlepDxyVertex0.push_back(track.dxy(vertices->at(0).position()));
            VlepDzVertex0.push_back(track.dz(vertices->at(0).position()));
            
            VlepChargedHadronIso.push_back(anelectron->chargedHadronIso());
            VlepNeutralHadronIso.push_back(anelectron->neutralHadronIso());
            VlepPhotonIso.push_back(anelectron->photonIso());
            VlepPuChargedHadronIso.push_back(anelectron->puChargedHadronIso());
            VlepPfIso.push_back((anelectron->chargedHadronIso()
                                 + std::max(0., anelectron->neutralHadronIso()
                                               + anelectron->photonIso()
                                               - 0.5*anelectron->puChargedHadronIso())
                                ) / anelectron->pt());

            Vlep.push_back(anelectron->polarP4());
            VlepPdgId.push_back(anelectron->pdgId());

            if ( fabs(anelectron->superCluster()->eta()) <= 1.479 ) //barrel region
            {
                VlepCombIso.push_back(( anelectron->dr03TkSumPt() +TMath::Max(0.,anelectron->dr03EcalRecHitSumEt()-1.)+anelectron->dr03HcalTowerSumEt())/TMath::Max(20.,anelectron->pt()));
            }
            else   //endcap
            {
                VlepCombIso.push_back(( anelectron->dr03TkSumPt() +anelectron->dr03EcalRecHitSumEt() +anelectron->dr03HcalTowerSumEt())/TMath::Max(20.,anelectron->pt()));
            }

            int triggerResult = 0;

            //          std::cout << "processing electron, is lepton object: " << Vlep.size() << std::endl;

            pat::TriggerObjectStandAloneCollection triggers = anelectron->triggerObjectMatches();
            if (triggers.size() > 0) triggerResult = triggerMap_["general"];
            for (unsigned int i = 0; i < triggers.size(); ++i)
                triggerResult |= getTriggerBits(triggers[i].pathNames());
            VlepTrigger.push_back(triggerResult);

            //          std::cout << "trigger result: " << std::hex << triggerResult << std::dec << std::endl;

            ++anelectron;
        }
    }

    //store the information in which channel the lepton is reconstructed
    for (size_t i = 1; i < VlepPdgId.size(); ++i) {
        int product = VlepPdgId.at(0) * VlepPdgId.at(i);
        if (product < 0) {
            recoInChannel = product;
            break;
        }
    }

    ///////////////////////////////////Jet properties/////////////////////////

    double jetPTThresholdForMET_ =10.;
    double jetEMLimitForMET_ = 0.9;
//     bool ptJetTooLow;    
    for ( edm::View<pat::Jet>::const_iterator ajet  = jets->begin() ; ajet != jets->end(); ++ajet )
    {
        Vjet.push_back(ajet->polarP4());
        VjetJERSF.push_back(ajet->userFloat("jerSF"));
        if (! iEvent.isRealData()) {
            VjetPartonFlavour.push_back( ajet->partonFlavour());

            if (ajet->genJet()) {
                VassociatedGenJet.push_back(ajet->genJet()->polarP4());
            }
            else {
                VassociatedGenJet.push_back(nullP4);
            }
        }
        VjetBTagTCHE.push_back(ajet->bDiscriminator("trackCountingHighEffBJetTags"));
        VjetBTagTCHP.push_back(ajet->bDiscriminator("trackCountingHighPurBJetTags"));
        VjetBTagJetProbability.push_back(ajet->bDiscriminator("jetProbabilityBJetTags"));
        VjetBTagJetBProbability.push_back(ajet->bDiscriminator("jetBProbabilityBJetTags"));
        VjetBTagSSVHE.push_back(ajet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));
        VjetBTagSSVHP.push_back(ajet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
        VjetBTagCSV.push_back(ajet->bDiscriminator("combinedSecondaryVertexBJetTags"));
        VjetBTagCSVMVA.push_back(ajet->bDiscriminator("combinedSecondaryVertexMVABJetTags"));
    }
    
    
    
    for(std::vector<JetProperties>::const_iterator i_jetProperties = jetPropertiesHandle->begin(); i_jetProperties != jetPropertiesHandle->end(); ++i_jetProperties){
        VjetChargeGlobalPtWeighted.push_back(i_jetProperties->jetChargeGlobalPtWeighted());
        VjetChargeRelativePtWeighted.push_back(i_jetProperties->jetChargeRelativePtWeighted());
        VjetAssociatedPartonPdgId.push_back(i_jetProperties->jetAssociatedPartonPdgId());
        VjetAssociatedParton.push_back(i_jetProperties->jetAssociatedParton());
        
        for (size_t i_pfTrack=0;i_pfTrack != i_jetProperties->jetPfCandidateTrack().size();i_pfTrack++)
        {
            VjetPfCandidateTrack.push_back(i_jetProperties->jetPfCandidateTrack().at(i_pfTrack));
            VjetPfCandidateTrackCharge.push_back(i_jetProperties->jetPfCandidateTrackCharge().at(i_pfTrack));
            VjetPfCandidateTrackId.push_back(i_jetProperties->jetPfCandidateTrackId().at(i_pfTrack));
            VjetPfCandidateTrackIndex.push_back(i_jetProperties-jetPropertiesHandle->begin());
        }
        
        for (size_t i_selTrack=0;i_selTrack != i_jetProperties->jetSelectedTrack().size();i_selTrack++)
        {
            VjetSelectedTrack.push_back(i_jetProperties->jetSelectedTrack().at(i_selTrack));
            VjetSelectedTrackIPValue.push_back(i_jetProperties->jetSelectedTrackIPValue().at(i_selTrack));
            VjetSelectedTrackIPSignificance.push_back(i_jetProperties->jetSelectedTrackIPSignificance().at(i_selTrack));
            VjetSelectedTrackCharge.push_back(i_jetProperties->jetSelectedTrackCharge().at(i_selTrack));
            VjetSelectedTrackIndex.push_back(i_jetProperties-jetPropertiesHandle->begin());
        }
    }


    //Here I create a separate jet collection needed for the on-the-fly calculation of jet uncertainties
    //because even bad-id jets are used for the MET

    //    for ( edm::View<pat::Jet>::const_iterator ajet  = jetsForMETuncorr->begin() ; ajet != jetsForMETuncorr->end(); ++ajet )
    for ( size_t jet_it =  0 ; jet_it < jetsForMETuncorr->size(); ++jet_it )
    {
        if(jetsForMETuncorr->at(jet_it).correctedJet("Uncorrected").pt() > jetPTThresholdForMET_ &&
            ((!jetsForMETuncorr->at(jet_it).isPFJet() && jetsForMETuncorr->at(jet_it).emEnergyFraction() < jetEMLimitForMET_) ||
            (jetsForMETuncorr->at(jet_it).isPFJet() && jetsForMETuncorr->at(jet_it).neutralEmEnergyFraction() + jetsForMETuncorr->at(jet_it).chargedEmEnergyFraction() < jetEMLimitForMET_)))
        {
            VjetForMET.push_back(jetsForMET->at(jet_it).polarP4());
            VjetForMETJERSF.push_back(jetsForMET->at(jet_it).userFloat("jerSF"));
            if (! iEvent.isRealData()) {
                VjetPartonFlavourForMET.push_back( jetsForMET->at(jet_it).partonFlavour());

                if (jetsForMET->at(jet_it).genJet()) {
                    VassociatedGenJetForMET.push_back(jetsForMET->at(jet_it).genJet()->polarP4());
                }
                else {
                    VassociatedGenJetForMET.push_back(nullP4);
                }
            }
        }
    }

    ///////////////////////////////////Met properties///////////////////////////

    edm::Handle<edm::View<pat::MET> > h_met;
    iEvent.getByLabel(met_, h_met);
    met = h_met->at(0).polarP4();

    //////////////////////////////Event Info/////////////////////
    runno = iEvent.id().run();
    lumibl = iEvent.id().luminosityBlock();
    eventno = iEvent.id().event();

    //create the PDF weights
    if (includePDFWeights_ && !iEvent.isRealData()) {
        edm::Handle<std::vector<double> > weightHandle;
        iEvent.getByLabel(pdfWeightTag_, weightHandle);
        VPdfWeights = *weightHandle;
    }

    Ntuple->Fill();
}



int NTupleWriter::getTriggerBits (const edm::Event &iEvent, const edm::Handle< edm::TriggerResults > &trigResults )
{
    int n_Triggers = trigResults->size();
    edm::TriggerNames trigName = iEvent.triggerNames(*trigResults);
    int result = 0;

    for ( int i_Trig = 0; i_Trig<n_Triggers; ++i_Trig )
    {
        if ( trigResults.product()->accept(i_Trig))
        {
            if (includeTrig_) VfiredTriggers.push_back(trigName.triggerName(i_Trig));
            const std::string &triggerName = trigName.triggerName(i_Trig);
            std::string triggerNameWithoutVersion(triggerName);
            while (triggerNameWithoutVersion.length() > 0
                    && triggerNameWithoutVersion[triggerNameWithoutVersion.length()-1] >= '0'
                    && triggerNameWithoutVersion[triggerNameWithoutVersion.length()-1] <= '9')
            {
                triggerNameWithoutVersion.replace(triggerNameWithoutVersion.length()-1, 1, "");
            }
            result |= triggerMap_[triggerNameWithoutVersion + "*"];
            result |= triggerMap_[triggerName];
        }
    }
    return result;
}



int NTupleWriter::getTriggerBits(const std::vector< std::string > &trigName)
{
    int result = 0;

    for ( unsigned int i_Trig = 0; i_Trig < trigName.size(); ++i_Trig ) {
        const std::string &triggerName = trigName.at( i_Trig);
        std::string triggerNameWithoutVersion(triggerName);
        while (triggerNameWithoutVersion.length() > 0
                && triggerNameWithoutVersion[triggerNameWithoutVersion.length()-1] >= '0'
                && triggerNameWithoutVersion[triggerNameWithoutVersion.length()-1] <= '9')
        {
            triggerNameWithoutVersion.replace(triggerNameWithoutVersion.length()-1, 1, "");
        }
        result |= triggerMap_[triggerNameWithoutVersion + "*"];
        result |= triggerMap_[triggerName];
    }
    return result;
}



int NTupleWriter::getTriggerBitsTau(const edm::Event &iEvent, const edm::Handle< edm::TriggerResults > &trigResults)
{
    int n_Triggers = trigResults->size();
    edm::TriggerNames trigName = iEvent.triggerNames(*trigResults);
    int result = 0;

    for ( int i_Trig = 0; i_Trig<n_Triggers; ++i_Trig )
    {
        if ( trigResults.product()->accept(i_Trig))
        {
            const std::string &triggerName = trigName.triggerName(i_Trig);
            std::string triggerNameWithoutVersion(triggerName);
            while (triggerNameWithoutVersion.length() > 0
                    && triggerNameWithoutVersion[triggerNameWithoutVersion.length()-1] >= '0'
                    && triggerNameWithoutVersion[triggerNameWithoutVersion.length()-1] <= '9')
            {
                triggerNameWithoutVersion.replace(triggerNameWithoutVersion.length()-1, 1, "");
            }
            result |= triggerMapTau_[triggerNameWithoutVersion + "*"];
            result |= triggerMapTau_[triggerName];
        }
    }
    return result;
}



void NTupleWriter::AssignLeptonAndTau(const reco::GenParticle* lepton, LV& GenLepton, int& pdgid, LV& GenTau )
{
    const reco::GenParticle *finalLepton;
    if (isTau(lepton)) {
        GenTau = lepton->polarP4();
        finalLepton = getTauDaughter(lepton);
    } else {
        GenTau = nullP4;
        finalLepton = lepton;
    }
    if (!isTau(finalLepton)) {
        GenLepton = finalLepton->polarP4();
        pdgid = finalLepton->pdgId();
    } else {
        GenLepton = nullP4;
        pdgid = 0;
    }
}



bool NTupleWriter::isTau(const reco::GenParticle *lepton) {
    return std::abs(lepton->pdgId()) == 15;
}



const reco::GenParticle*
NTupleWriter::getTauDaughter(const reco::GenParticle* tau)
{
    for (size_t i = 0; i < tau->numberOfDaughters(); ++i) {
        const reco::GenParticle* daughter = dynamic_cast<const reco::GenParticle*>(tau->daughter(i));
        if (std::abs(daughter->pdgId()) == 11 || std::abs(daughter->pdgId()) == 13) return daughter;
        else if (isTau(daughter)) return getTauDaughter(daughter);
    }
    return tau;
}





// ------------ method called once each job just before starting event loop  ------------
void
NTupleWriter::beginJob()
{
    edm::Service<TFileService> fs;

    if ( !fs ) throw edm::Exception(edm::errors::Configuration, "TFile Service is not registered in cfg file");

    Ntuple = fs->make<TTree>("NTuple","NTuple");

    //TFileDirectory header = fs->mkdir("Header");
    //header.cd();
    TObjString(sampleName_.c_str()).Write("sampleName");
    TObjString(channelName_.c_str()).Write("channelName");
    TObjString(systematicsName_.c_str()).Write("systematicsName");
    TObjString(boost::lexical_cast<std::string>(isMC_).c_str()).Write("isMC");
    TObjString(boost::lexical_cast<std::string>(isTtBarSample_).c_str()).Write("isSignal");
    TObjString(boost::lexical_cast<std::string>(isHiggsSample_).c_str()).Write("isHiggsSignal");

    ///////////////dilepton and lepton properties//////////
    Ntuple->Branch("leptons", &Vlep);
    Ntuple->Branch("lepPdgId", &VlepPdgId);
    Ntuple->Branch("lepID", &VlepID);
    Ntuple->Branch("lepPfIso", &VlepPfIso);
    Ntuple->Branch("lepChargedHadronIso", &VlepChargedHadronIso);
    Ntuple->Branch("lepNeutralHadronIso", &VlepNeutralHadronIso);
    Ntuple->Branch("lepPhotonIso", &VlepPhotonIso);
    Ntuple->Branch("lepPuChargedHadronIso", &VlepPuChargedHadronIso);
    Ntuple->Branch("lepCombIso", &VlepCombIso);
    Ntuple->Branch("lepDxyVertex0", &VlepDxyVertex0);
    Ntuple->Branch("lepDzVertex0", &VlepDzVertex0);
    Ntuple->Branch("lepTrigger", &VlepTrigger);


    /////////////jet properties////////////
    Ntuple->Branch("jets", &Vjet);
    Ntuple->Branch("jetJERSF", &VjetJERSF);
    Ntuple->Branch("jetsForMET", &VjetForMET);
    Ntuple->Branch("jetForMETJERSF", &VjetForMETJERSF);
    Ntuple->Branch("jetBTagTCHE", &VjetBTagTCHE);
    Ntuple->Branch("jetBTagTCHP", &VjetBTagTCHP);
    Ntuple->Branch("jetBTagSSVHE", &VjetBTagSSVHE);
    Ntuple->Branch("jetBTagSSVHP", &VjetBTagSSVHP);
    Ntuple->Branch("jetBTagJetProbability", &VjetBTagJetProbability);
    Ntuple->Branch("jetBTagJetBProbability", &VjetBTagJetBProbability);
    Ntuple->Branch("jetBTagCSV", &VjetBTagCSV);
    Ntuple->Branch("jetBTagCSVMVA", &VjetBTagCSVMVA);
    Ntuple->Branch("pdfWeights", &VPdfWeights);
    Ntuple->Branch("jetPartonFlavour", &VjetPartonFlavour);
    Ntuple->Branch("jetPartonFlavourForMET", &VjetPartonFlavourForMET);
    Ntuple->Branch("allGenJets", &VallGenJets);
    Ntuple->Branch("associatedGenJet", &VassociatedGenJet);
    Ntuple->Branch("associatedGenJetForMET", &VassociatedGenJetForMET);

    Ntuple->Branch("jetChargeGlobalPtWeighted", &VjetChargeGlobalPtWeighted);
    Ntuple->Branch("jetChargeRelativePtWeighted", &VjetChargeRelativePtWeighted);
    Ntuple->Branch("jetPfCandidateTrack", &VjetPfCandidateTrack);
    Ntuple->Branch("jetPfCandidateTrackCharge", &VjetPfCandidateTrackCharge);
    Ntuple->Branch("jetPfCandidateTrackIndex", &VjetPfCandidateTrackIndex);
    Ntuple->Branch("jetSelectedTrack", &VjetSelectedTrack);
    Ntuple->Branch("jetSelectedTrackIPValue", &VjetSelectedTrackIPValue);
    Ntuple->Branch("jetSelectedTrackIPSignificance", &VjetSelectedTrackIPSignificance);
    Ntuple->Branch("jetSelectedTrackCharge", &VjetSelectedTrackCharge);
    Ntuple->Branch("jetSelectedTrackIndex", &VjetSelectedTrackIndex);
    
    if (isTtBarSample_) Ntuple->Branch("jetAssociatedPartonPdgId", &VjetAssociatedPartonPdgId);

    /////////////met properties///////////
    Ntuple->Branch("met", &met);

    ///////////event info///////////
    Ntuple->Branch("runNumber", &runno, "runNumber/i");
    Ntuple->Branch("lumiBlock", &lumibl,"lumiBlock/i");
    Ntuple->Branch("eventNumber", &eventno, "eventNumber/i");
    Ntuple->Branch("triggerBits", &triggerBits, "triggerBits/i");
    Ntuple->Branch("triggerBitsTau", &triggerBitsTau, "triggerBitsTau/i");


    ////////////triggers//////////////////
    if ( includeTrig_ ) Ntuple->Branch("firedTriggers", &VfiredTriggers);

    ////////////weight////////////
    Ntuple->Branch("weightGenerator", &weightGenerator, "weightGenerator/D");

    /////////vertices
    Ntuple->Branch("vertMulti", &vertMulti, "vertMulti/I");
    Ntuple->Branch("vertMultiTrue", &vertMultiTrue, "vertMultiTrue/I");

    ////////Gen Info
    if (isTtBarSample_) {
        Ntuple->Branch("TopProductionMode", &TopProductionMode, "TopProductionMode/I");
        Ntuple->Branch("GenMET", &GenMET);
        Ntuple->Branch("GenTop", &GenTop);
        Ntuple->Branch("GenAntiTop", &GenAntiTop);
        Ntuple->Branch("GenLepton", &GenLepton);
        Ntuple->Branch("GenAntiLepton", &GenAntiLepton);
        Ntuple->Branch("GenLeptonPdgId", &GenLeptonPdgId, "GenLeptonPdgId/I");
        Ntuple->Branch("GenAntiLeptonPdgId", &GenAntiLeptonPdgId, "GenAntiLeptonPdgId/I");
        Ntuple->Branch("GenTau", &GenTau);
        Ntuple->Branch("GenAntiTau", &GenAntiTau);
        Ntuple->Branch("GenNeutrino", &GenNeutrino);
        Ntuple->Branch("GenAntiNeutrino", &GenAntiNeutrino);
        Ntuple->Branch("GenB", &GenB);
        Ntuple->Branch("GenAntiB", &GenAntiB);
        Ntuple->Branch("GenWPlus", &GenWPlus);
        Ntuple->Branch("GenWMinus", &GenWMinus);
        //     Ntuple->Branch("GenParticleP4", &GenParticleP4);
        //     Ntuple->Branch("GenParticlePdgId", &GenParticlePdgId);
        //     Ntuple->Branch("GenParticleStatus", &GenParticleStatus);
        
        //     Ntuple->Branch("GenJetHadronB", &HadronGenB);
        //     Ntuple->Branch("GenJetHadronAntiB", &HadronGenAntiB);

        Ntuple->Branch("BHadJetIndex", &VBHadJetIdx);
        Ntuple->Branch("AntiBHadJetIndex", &VAntiBHadJetIdx);
        Ntuple->Branch("BHadrons", &VBHadron);
        Ntuple->Branch("AntiBHadrons", &VAntiBHadron);
        Ntuple->Branch("BHadronFromTop", &VBHadFromTop);
        Ntuple->Branch("AntiBHadronFromTopB", &VAntiBHadFromTop);
        Ntuple->Branch("BHadronVsJet", &VBHadVsJet);
        Ntuple->Branch("AntiBHadronVsJet", &VAntiBHadVsJet);

        Ntuple->Branch("genExtraTopJetNumberId", &genExtraTopJetNumberId);
        Ntuple->Branch("genBHadPlusMothers", &VgenBHadPlusMothers);
        Ntuple->Branch("genBHadPlusMothersPdgId", &VgenBHadPlusMothersPdg);
        if(saveHadronMothers) {
            Ntuple->Branch("genBHadPlusMothersStatus", &VgenBHadPlusMothersStatus);
            Ntuple->Branch("genBHadPlusMothersIndices", &VgenBHadPlusMothersIndices);
        }
        Ntuple->Branch("genBHadIndex", &VgenBHadIndex);
        Ntuple->Branch("genBHadFlavour", &VgenBHadFlavour);
        Ntuple->Branch("genBHadJetIndex", &VgenBHadJetIndex);
        Ntuple->Branch("genBHadLeptonIndex", &VgenBHadLeptonIndex);
        Ntuple->Branch("genBHadLeptonHadronIndex", &VgenBHadLeptonHadronIndex);
        Ntuple->Branch("genBHadLeptonViaTau", &VgenBHadLeptonViaTau);
        Ntuple->Branch("genBHadFromTopWeakDecay", &VgenBHadFromTopWeakDecay);
        
        Ntuple->Branch("genCHadJetIndex", &VgenCHadJetIndex);
        Ntuple->Branch("genCHadLeptonHadronIndex", &VgenCHadLeptonHadronIndex);
        Ntuple->Branch("genCHadLeptonViaTau", &VgenCHadLeptonViaTau);
        Ntuple->Branch("genCHadFromBHadron", &VgenCHadFromBHadron);

        Ntuple->Branch("jetAssociatedParton", &VjetAssociatedParton);
    }

    //Gen Info for Higgs and b quarks of decay
    if (isHiggsSample_) {
        Ntuple->Branch("GenH", &GenH);
        Ntuple->Branch("GenBFromH", &GenBFromH);
        Ntuple->Branch("GenAntiBFromH", &GenAntiBFromH);
        Ntuple->Branch("HiggsDecayMode", &HiggsDecayMode, "HiggsDecayMode/I");
    }

    // Gen Info for Zs and decay products
    if(isZSample_){
        Ntuple->Branch("GenZ", &VGenZ);
        Ntuple->Branch("GenZMeDaughterParticle", &VGenZMeDaughterParticle);
        Ntuple->Branch("GenZMeDaughterAntiParticle", &VGenZMeDaughterAntiParticle);
        Ntuple->Branch("GenZStableLepton", &VGenZStableLepton);
        Ntuple->Branch("GenZStableAntiLepton", &VGenZStableAntiLepton);
        Ntuple->Branch("GenZDecayMode", &VGenZDecayMode);
    }

    //Hypothesis Info
    Ntuple->Branch("HypTop", &HypTop);
    Ntuple->Branch("HypAntiTop", &HypAntiTop);
    Ntuple->Branch("HypLepton", &HypLepton);
    Ntuple->Branch("HypAntiLepton", &HypAntiLepton);
    Ntuple->Branch("HypNeutrino", &HypNeutrino);
    Ntuple->Branch("HypAntiNeutrino", &HypAntiNeutrino);
    Ntuple->Branch("HypB", &HypB);
    Ntuple->Branch("HypAntiB", &HypAntiB);
    Ntuple->Branch("HypWPlus", &HypWPlus);
    Ntuple->Branch("HypWMinus", &HypWMinus);
    Ntuple->Branch("HypJet0index", &HypJet0index);
    Ntuple->Branch("HypJet1index", &HypJet1index);

    Ntuple->Branch("TopDecayMode", &TopDecayMode, "TopDecayMode/I");
    Ntuple->Branch("ZDecayMode", &ZDecayMode);
    Ntuple->Branch("recoInChannel", &recoInChannel, "recoInChannel/I");
}



void NTupleWriter::clearVariables()
{
    runno = 0;
    lumibl = 0;
    eventno = 0;
    triggerBits = 0;
    triggerBitsTau = 0;
    recoInChannel = 0;

    ////////dileptons and leptons/////
    Vlep.clear() ;
    VlepID.clear() ;
    VlepPdgId.clear();
    VlepPfIso.clear();
    VlepChargedHadronIso.clear();
    VlepNeutralHadronIso.clear();
    VlepPhotonIso.clear();
    VlepPuChargedHadronIso.clear();
    VlepCombIso.clear();
    VlepDxyVertex0.clear();
    VlepDzVertex0.clear();
    VlepTrigger.clear();

    GenParticleP4.clear();
    GenParticlePdgId.clear();
    GenParticleStatus.clear();

    /////////jets///////////
    Vjet.clear();
    VjetJERSF.clear();
    VjetForMET.clear();
    VjetForMETJERSF.clear();
    VjetPartonFlavour.clear();
    VjetPartonFlavourForMET.clear();
    VjetBTagTCHE.clear();
    VjetBTagTCHP.clear();
    VjetBTagSSVHE.clear();
    VjetBTagSSVHP.clear();
    VjetBTagJetProbability.clear();
    VjetBTagJetBProbability.clear();
    VjetBTagCSV.clear();
    VjetBTagCSVMVA.clear();
    VPdfWeights.clear();
    VallGenJets.clear();
    VassociatedGenJet.clear();
    VassociatedGenJetForMET.clear();

    VjetChargeGlobalPtWeighted.clear();
    VjetChargeRelativePtWeighted.clear();
    VjetAssociatedPartonPdgId.clear();
    VjetAssociatedParton.clear();
    VjetPfCandidateTrack.clear();
    VjetPfCandidateTrackCharge.clear();
    VjetPfCandidateTrackIndex.clear();
    VjetSelectedTrack.clear();
    VjetSelectedTrackIPValue.clear();
    VjetSelectedTrackIPSignificance.clear();
    VjetSelectedTrackCharge.clear();
    VjetSelectedTrackIndex.clear();

    VBHadJetIdx.clear();
    VAntiBHadJetIdx.clear();
    VBHadron.clear();
    VAntiBHadron.clear();
    VBHadFromTop.clear();
    VAntiBHadFromTop.clear();
    VBHadVsJet.clear();
    VAntiBHadVsJet.clear();

    genExtraTopJetNumberId=-1;
    VgenBHadPlusMothers.clear();
    VgenBHadPlusMothersIndices.clear();
    VgenBHadPlusMothersStatus.clear();
    VgenBHadPlusMothersPdg.clear();
    VgenBHadIndex.clear();
    VgenBHadFlavour.clear();
    VgenBHadJetIndex.clear();
    VgenBHadLeptonIndex.clear();
    VgenBHadLeptonHadronIndex.clear();
    VgenBHadLeptonViaTau.clear();
    VgenBHadFromTopWeakDecay.clear();
    
    VgenCHadJetIndex.clear();
    VgenCHadLeptonHadronIndex.clear();
    VgenCHadLeptonViaTau.clear();
    VgenCHadFromTopWeakDecay.clear();
    VgenCHadFromBHadron.clear();

    /////////Triggers/////////
    VfiredTriggers.clear();

    //////vertices
    vertMulti=0;
    vertMultiTrue=0;

    HypTop.clear();
    HypAntiTop.clear();
    HypLepton.clear();
    HypAntiLepton.clear();
    HypNeutrino.clear();
    HypAntiNeutrino.clear();
    HypB.clear();
    HypAntiB.clear();
    HypWPlus.clear();
    HypWMinus.clear();

    HypJet0index.clear();
    HypJet1index.clear();

    VGenZ.clear();
    VGenZMeDaughterParticle.clear();
    VGenZMeDaughterAntiParticle.clear();
    VGenZStableLepton.clear();
    VGenZStableAntiLepton.clear();
    VGenZDecayMode.clear();

    ZDecayMode.clear();
    TopDecayMode = 0;
    HiggsDecayMode = 0;
}



// ------------ method called once each job just after ending the event loop  ------------
void
NTupleWriter::endJob()
{
}



// ------------ method called when starting to processes a run  ------------
void
NTupleWriter::beginRun(edm::Run const&, edm::EventSetup const& )
{
}



// ------------ method called when ending the processing of a run  ------------
void
NTupleWriter::endRun(edm::Run const&, edm::EventSetup const& )
{
}



// ------------ method called when starting to processes a luminosity block  ------------
void
NTupleWriter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const& )
{
}



// ------------ method called when ending the processing of a luminosity block  ------------
void
NTupleWriter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const& )
{
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NTupleWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions )
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}





//define this as a plug-in
DEFINE_FWK_MODULE(NTupleWriter);
