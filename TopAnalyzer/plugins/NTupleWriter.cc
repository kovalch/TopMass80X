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
#include <iostream>
#include <string>
#include <map>
#include <algorithm>
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
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h" //###############
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TopAnalysis/HiggsUtils/interface/HiggsGenEvent.h"
#include "TopAnalysis/HiggsUtils/interface/JetProperties.h"
#include "TopAnalysis/HiggsUtils/interface/GenZDecayProperties.h"

#include <Math/Vector3D.h>
#include <TTree.h>
#include <TLorentzVector.h>





// Typedef for Lorentz Vector as used in ntuples
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LV;



//
// class declaration
//
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
    int getTriggerBits(const edm::Event& iEvent, const edm::Handle<edm::TriggerResults>& triggerResults);
    int getTriggerBits(const std::vector<std::string>& triggerNames);
    int getTriggerBitsTau(const edm::Event& iEvent, const edm::Handle<edm::TriggerResults>& triggerResults);
    
    void assignLeptonAndTau(const reco::GenParticle* lepton, LV& genLepton, int& pdgId, LV& genTau);
    bool isTau(const reco::GenParticle* lepton);
    const reco::GenParticle* tauDaughter(const reco::GenParticle* tau);
    int nUniqueElementsInVector(const std::vector<int> vector, const bool countOnlyNonRepeating = false);
    
    // ----------member data ---------------------------
    
    
    
    // Sample settings
    const std::string sampleName_;
    const std::string channelName_;
    const std::string systematicsName_;
    const bool isMC_;
    const bool isTtbarSample_;
    const bool isHiggsSample_;
    const bool isZSample_;
    const bool isMadgraphSample_;
    const bool includeTrigger_;
    const bool includePdfWeights_;
    const bool saveHadronMothers_;
    const bool saveCHadronParticles_;
    
    
    
    // Input tags...
    
    // ... for reco-level
    const edm::InputTag electronsTag_;
    const edm::InputTag muonsTag_;
    const edm::InputTag jetsTag_;
    const edm::InputTag jetsForMetTag_;
    const edm::InputTag jetsForMetUncorrectedTag_;
    const edm::InputTag jetPropertiesTag_;
    const edm::InputTag metTag_;
    const edm::InputTag mvaMetTag_;
    const edm::InputTag verticesTag_;
    const edm::InputTag triggerResultsTag_;
    
    // ... for gen-level
    const edm::InputTag pileupInfoTag_;
    const edm::InputTag genParticlesTag_;
    const edm::InputTag genJetsTag_;
    const edm::InputTag pdfWeightTag_;
    const edm::InputTag genEventTtbarTag_;
    const edm::InputTag genEventHiggsTag_;
    const edm::InputTag genZDecayTag_;
    const edm::InputTag bHadJetIndexTag_;
    const edm::InputTag antiBHadJetIndexTag_;
    const edm::InputTag BHadronsTag_;
    const edm::InputTag AntiBHadronsTag_;
    const edm::InputTag BHadronFromTopBTag_;
    const edm::InputTag AntiBHadronFromTopBTag_;
    const edm::InputTag BHadronVsJetTag_;
    const edm::InputTag AntiBHadronVsJetTag_;
    const edm::InputTag genBHadPlusMothersTag_;
    const edm::InputTag genBHadPlusMothersIndicesTag_;
    const edm::InputTag genBHadIndexTag_;
    const edm::InputTag genBHadFlavourTag_;
    const edm::InputTag genBHadJetIndexTag_;
    const edm::InputTag genBHadLeptonIndexTag_;
    const edm::InputTag genBHadLeptonHadronIndexTag_;
    const edm::InputTag genBHadLeptonViaTauTag_;
    const edm::InputTag genBHadFromTopWeakDecayTag_;
    const edm::InputTag genCHadPlusMothersTag_;
    const edm::InputTag genCHadPlusMothersIndicesTag_;
    const edm::InputTag genCHadIndexTag_;
    const edm::InputTag genCHadJetIndexTag_;
    const edm::InputTag genCHadLeptonIndexTag_;
    const edm::InputTag genCHadLeptonHadronIndexTag_;
    const edm::InputTag genCHadLeptonViaTauTag_;
    const edm::InputTag genCHadFromTopWeakDecayTag_;
    const edm::InputTag genCHadBHadronIdTag_;
    const edm::InputTag ttbarDecayModeTag_;
    const edm::InputTag higgsDecayModeTag_;
    const edm::InputTag madgraphWDecayTag_;
    
    
    
    // Ntuple and branch variables ...
    TTree* ntuple_;
    
    // ... for event info
    unsigned int runNumber_;
    unsigned int lumiBlock_;
    unsigned int eventNumber_;
    
    // ... for triggers
    unsigned int triggerBits_;
    unsigned int triggerBitsTau_;
    std::vector<std::string> v_firedTriggers_;
    
    // ... for vertices
    int vertexMultiplicity_;
    
    // ... for leptons
    std::vector<LV> v_lepton_;
    std::vector<int> v_leptonPdgId_;
    std::vector<double> v_leptonId_; // mvaID for electrons (-1 for muon)
    std::vector<double> v_leptonChargedHadronIso_;
    std::vector<double> v_leptonNeutralHadronIso_;
    std::vector<double> v_leptonPhotonIso_;
    std::vector<double> v_leptonPuChargedHadronIso_;
    std::vector<double> v_leptonPfIso_;
    std::vector<double> v_leptonCombIso_;
    std::vector<double> v_leptonDxyVertex0_;
    std::vector<double> v_leptonDzVertex0_;
    std::vector<int> v_leptonTrigger_;
    
    // ... for jets
    std::vector<LV> v_jet_;
    std::vector<double> v_jetBtagTCHE_;
    std::vector<double> v_jetBtagTCHP_;
    std::vector<double> v_jetBtagJetProbability_;
    std::vector<double> v_jetBtagJetBProbability_;
    std::vector<double> v_jetBtagSSVHE_;
    std::vector<double> v_jetBtagSSVHP_;
    std::vector<double> v_jetBtagCSV_;
    std::vector<double> v_jetBtagCSVMVA_;
    
    // ... for jet properties
    std::vector<double> v_jetChargeGlobalPtWeighted_;
    std::vector<double> v_jetChargeRelativePtWeighted_;
    std::vector<LV> v_jetPfCandidateTrack_;
    std::vector<int> v_jetPfCandidateTrackCharge_;
    std::vector<int> v_jetPfCandidateTrackId_;
    std::vector<int> v_jetPfCandidatePrimaryVertexId_;
    std::vector<int> v_jetPfCandidateTrackIndex_;
    std::vector<int> v_jetSelectedTrackMatchToPfCandidateIndex_;
    std::vector<LV> v_jetSelectedTrack_;
    std::vector<double> v_jetSelectedTrackIPValue_;
    std::vector<double> v_jetSelectedTrackIPSignificance_;
    std::vector<int> v_jetSelectedTrackCharge_;
    std::vector<int> v_jetSelectedTrackIndex_;
    std::vector<LV> v_jetSecondaryVertex_;
    std::vector<double> v_jetSecondaryVertexPtCorrectedMass_;
    std::vector<int> v_jetSecondaryVertexJetIndex_;
    std::vector<double> v_jetSecondaryVertexFlightDistanceValue_;
    std::vector<double> v_jetSecondaryVertexFlightDistanceSignificance_;
    std::vector<int> v_jetSecondaryVertexTrackMatchToSelectedTrackIndex_;
    std::vector<int> v_jetSecondaryVertexTrackVertexIndex_;
    
    // ... for MET
    LV met_;
    LV mvaMet_;
    
    // ... for jets, but used only for MC corrections (not needed in data)
    std::vector<LV> v_jetForMet_;
    std::vector<double> v_jetJerSF_;
    std::vector<double> v_jetForMetJerSF_;
    
    // ... for general true level info
    int vertexMultiplicityTrue_;
    std::vector<LV> v_associatedGenJet_;
    std::vector<LV> v_associatedGenJetForMet_;
    std::vector<int> v_jetPartonFlavour_;
    std::vector<int> v_jetPartonFlavourForMet_;
    std::vector<double> v_pdfWeight_;
    double weightGenerator_;
    
    // ... for ttbar (+decay) generator info
    int ttbarProductionMode_;
    int ttbarDecayMode_;
    LV genTop_;
    LV genAntiTop_;
    LV genLepton_;
    LV genAntiLepton_;
    int genLeptonPdgId_;
    int genAntiLeptonPdgId_;
    LV genTau_;
    LV genAntiTau_;
    LV genNeutrino_;
    LV genAntiNeutrino_;
    LV genB_;
    LV genAntiB_;
    LV genWPlus_;
    LV genWMinus_;
    LV genMet_;
    std::vector<LV> v_allGenJet_;
    std::vector<LV> v_genParticleP4_;
    std::vector<int> v_genParticlePdgId_;
    std::vector<int> v_genParticleStatus_;
    std::vector<int> v_bHadJetIndex_;
    std::vector<int> v_antiBHadJetIndex_;
    std::vector<LV> v_bHadron_;
    std::vector<LV> v_antiBHadron_;
    std::vector<bool> v_bHadFromTop_;
    std::vector<bool> v_antiBHadFromTop_;
    std::vector<int> v_bHadVsJet_;
    std::vector<int> v_antiBHadVsJet_;
    std::vector<int> v_jetAssociatedPartonPdgId_;
    std::vector<LV> v_jetAssociatedParton_;
    
    // ... for Higgs (+decay) generator info
    int higgsDecayMode_;
    LV genH_;
    LV genBFromH_;
    LV genAntiBFromH_;
    
    // ... for Z boson (+decay) generator info
    std::vector<LV> v_genZ_;
    std::vector<LV> v_genZMeDaughterParticle_;
    std::vector<LV> v_genZMeDaughterAntiParticle_;
    std::vector<LV> v_genZStableLepton_;
    std::vector<LV> v_genZStableAntiLepton_;
    std::vector<int> v_genZDecayMode_;
    
    // ... for full parton-to-genJet matching of b/c quarks
    int genExtraTopJetNumberId_;
    std::vector<LV> v_genBHadPlusMothers_;
    std::vector<int> v_genBHadPlusMothersPdgId_;
    std::vector<int> v_genBHadPlusMothersStatus_;
    std::vector<std::vector<int> > v_genBHadPlusMothersIndices_;
    std::vector<int> v_genBHadIndex_;
    std::vector<int> v_genBHadFlavour_;
    std::vector<int> v_genBHadJetIndex_;
    std::vector<int> v_genBHadFromTopWeakDecay_;
    std::vector<int> v_genBHadLeptonHadronIndex_;
    std::vector<int> v_genBHadLeptonIndex_;
    std::vector<int> v_genBHadLeptonViaTau_;
    std::vector<LV> v_genCHadPlusMothers_;
    std::vector<int> v_genCHadPlusMothersPdgId_;
    std::vector<int> v_genCHadPlusMothersStatus_;
    std::vector<int> v_genCHadIndex_;
    std::vector<int> v_genCHadBHadronId_;
    std::vector<int> v_genCHadJetIndex_;
    std::vector<int> v_genCHadLeptonHadronIndex_;
    std::vector<int> v_genCHadLeptonIndex_;
    std::vector<int> v_genCHadLeptonViaTau_;
    
    // ... for W decays of MadGraph samples
    std::vector<int> v_madgraphWDecay_;
    
    
    
    // Helper variables
    std::map<std::string, int> triggerMap_;
    std::map<std::string, int> triggerMapTau_;
    const LV nullP4_;
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
sampleName_(iConfig.getParameter<std::string>("sampleName")),
channelName_(iConfig.getParameter<std::string>("channelName")),
systematicsName_(iConfig.getParameter<std::string>("systematicsName")),
isMC_(iConfig.getParameter<bool>("isMC")),
isTtbarSample_(iConfig.getParameter<bool>("isTtbarSample")),
isHiggsSample_(iConfig.getParameter<bool>("isHiggsSample")),
isZSample_(iConfig.getParameter<bool>("isZSample")),
isMadgraphSample_(iConfig.getParameter<bool>("isMadgraphSample")),
includeTrigger_(iConfig.getParameter<bool>("includeTrigger")),
includePdfWeights_ (iConfig.getParameter<bool>("includePdfWeights")),
saveHadronMothers_(iConfig.getParameter<bool>("saveHadronMothers")),
saveCHadronParticles_(iConfig.getParameter<bool>("saveCHadronParticles")),

electronsTag_(iConfig.getParameter<edm::InputTag>("electrons")),
muonsTag_(iConfig.getParameter<edm::InputTag>("muons")),
jetsTag_(iConfig.getParameter<edm::InputTag>("jets")),
jetsForMetTag_(iConfig.getParameter<edm::InputTag>("jetsForMet")),
jetsForMetUncorrectedTag_(iConfig.getParameter<edm::InputTag>("jetsForMetUncorrected")),
jetPropertiesTag_(iConfig.getParameter<edm::InputTag>("jetProperties")),
metTag_(iConfig.getParameter<edm::InputTag>("met")),
mvaMetTag_(iConfig.getParameter<edm::InputTag>("mvaMet")),
verticesTag_(iConfig.getParameter<edm::InputTag>("vertices")),
triggerResultsTag_(iConfig.getParameter<edm::InputTag>("triggerResults")),

pileupInfoTag_(iConfig.getParameter<edm::InputTag>("pileupInfo")),
genParticlesTag_(iConfig.getParameter<edm::InputTag>("genParticles")),
genJetsTag_(iConfig.getParameter<edm::InputTag>("genJets")),
pdfWeightTag_(iConfig.getParameter<edm::InputTag>("pdfWeights")),
genEventTtbarTag_(iConfig.getParameter<edm::InputTag>("genEventTtbar")),
genEventHiggsTag_(iConfig.getParameter<edm::InputTag>("genEventHiggs")),
genZDecayTag_(iConfig.getParameter<edm::InputTag>("genZDecay")),
bHadJetIndexTag_(iConfig.getParameter<edm::InputTag>("BHadJetIndex")),
antiBHadJetIndexTag_(iConfig.getParameter<edm::InputTag>("AntiBHadJetIndex")),
BHadronsTag_(iConfig.getParameter<edm::InputTag>("BHadrons")),
AntiBHadronsTag_(iConfig.getParameter<edm::InputTag>("AntiBHadrons")),
BHadronFromTopBTag_(iConfig.getParameter<edm::InputTag>("BHadronFromTopB")),
AntiBHadronFromTopBTag_(iConfig.getParameter<edm::InputTag>("AntiBHadronFromTopB")),
BHadronVsJetTag_(iConfig.getParameter<edm::InputTag>("BHadronVsJet")),
AntiBHadronVsJetTag_(iConfig.getParameter<edm::InputTag>("AntiBHadronVsJet")),
genBHadPlusMothersTag_(iConfig.getParameter<edm::InputTag>("genBHadPlusMothers")),
genBHadPlusMothersIndicesTag_(iConfig.getParameter<edm::InputTag>("genBHadPlusMothersIndices")),
genBHadIndexTag_(iConfig.getParameter<edm::InputTag>("genBHadIndex")),
genBHadFlavourTag_(iConfig.getParameter<edm::InputTag>("genBHadFlavour")),
genBHadJetIndexTag_(iConfig.getParameter<edm::InputTag>("genBHadJetIndex")),
genBHadLeptonIndexTag_(iConfig.getParameter<edm::InputTag>("genBHadLeptonIndex")),
genBHadLeptonHadronIndexTag_(iConfig.getParameter<edm::InputTag>("genBHadLeptonHadronIndex")),
genBHadLeptonViaTauTag_(iConfig.getParameter<edm::InputTag>("genBHadLeptonViaTau")),
genBHadFromTopWeakDecayTag_(iConfig.getParameter<edm::InputTag>("genBHadFromTopWeakDecay")),
genCHadPlusMothersTag_(iConfig.getParameter<edm::InputTag>("genCHadPlusMothers")),
genCHadPlusMothersIndicesTag_(iConfig.getParameter<edm::InputTag>("genCHadPlusMothersIndices")),
genCHadIndexTag_(iConfig.getParameter<edm::InputTag>("genCHadIndex")),
genCHadJetIndexTag_(iConfig.getParameter<edm::InputTag>("genCHadJetIndex")),
genCHadLeptonIndexTag_(iConfig.getParameter<edm::InputTag>("genCHadLeptonIndex")),
genCHadLeptonHadronIndexTag_(iConfig.getParameter<edm::InputTag>("genCHadLeptonHadronIndex")),
genCHadLeptonViaTauTag_(iConfig.getParameter<edm::InputTag>("genCHadLeptonViaTau")),
genCHadFromTopWeakDecayTag_(iConfig.getParameter<edm::InputTag>("genCHadFromTopWeakDecay")),
genCHadBHadronIdTag_(iConfig.getParameter<edm::InputTag>("genCHadBHadronId")),
ttbarDecayModeTag_(iConfig.getParameter<edm::InputTag>("ttbarDecayMode")),
higgsDecayModeTag_(iConfig.getParameter<edm::InputTag>("higgsDecayMode")),
madgraphWDecayTag_(iConfig.getParameter<edm::InputTag>("madgraphWDecay")),

nullP4_(0., 0., 0., 0.)
{
    // WARNING: The trigger map can be used either for a specific version, e.g. Trig_v6
    // or for any version, Trig_v*. NOT supported: Trig_v1* - the star captures ALL digits!
    
    // Use first 8 bits for mumu
    triggerMap_["HLT_DoubleMu6_v*"] = 1;
    triggerMap_["HLT_DoubleMu7_v*"] = 2;
    triggerMap_["HLT_Mu13_Mu8_v*"] = 4;
    triggerMap_["HLT_Mu17_Mu8_v*"] = 8;
    triggerMap_["HLT_DoubleMu45_v*"] = 0x10;
    triggerMap_["HLT_Mu17_TkMu8_v*"] = 0x20;
    
    // Use bits 9 to 16 for mu e
    triggerMap_["HLT_Mu8_Ele17_CaloIdL_v*"] = 0x100;
    triggerMap_["HLT_Mu17_Ele8_CaloIdL_v*"] = 0x200;
    triggerMap_["HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*"] = 0x400;
    triggerMap_["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*"] = 0x800;
    triggerMap_["HLT_Mu10_Ele10_CaloIdL_v*"] = 0x1000;
    triggerMap_["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"] = 0x2000;
    triggerMap_["HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"] = 0x4000;
    
    // Use bits 17-24 for ee
    triggerMap_["HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*"] = 0x10000;
    triggerMap_["HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*"] = 0x20000;
    triggerMap_["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"] = 0x40000;
    triggerMap_["HLT_DoubleEle45_CaloIdL_v*"] = 0x80000;
    
    // Use bit 32 for general trigger to avoid false for unmatched triggers
    triggerMap_["general"] = 0x80000000;
    
    
    
    // Trigger map for hadronic tau (cross) triggers, analog to trigger map above
    
    // Use first 8 bits for mu + tau
    triggerMapTau_["HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v*"] = 1;
    triggerMapTau_["HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v*"] = 2;
    
    // Use bits 9 to 16 for e + tau
    triggerMapTau_["HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v*"] = 0x100;
    triggerMapTau_["HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*"] = 0x200;
    
    // Use bits 17-24 for tau + tau + jet
    triggerMapTau_["HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30_v*"] = 0x10000;
    triggerMapTau_["HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30_v*"] = 0x20000;
    triggerMapTau_["HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30_v*"] = 0x40000;

    // Use bit 32 for general trigger to avoid false for unmatched triggers - is this needed?
    //triggerMapTau_["general"] = 0x80000000;
    
    
    
    // Set ntuple branch variables to default values
    this->clearVariables();
}



NTupleWriter::~NTupleWriter()
{}



//
// member functions
//

// ------------ method called for each event  ------------
void
NTupleWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    this->clearVariables();
    
    // Event info
    runNumber_ = iEvent.id().run();
    lumiBlock_ = iEvent.id().luminosityBlock();
    eventNumber_ = iEvent.id().event();
    
    // MC weights for generators providing weighted events
    if(!iEvent.isRealData()){
        const edm::InputTag genEventInfoTag("generator");
        edm::Handle<GenEventInfoProduct> genEventInfo;
        iEvent.getByLabel(genEventInfoTag, genEventInfo);
        weightGenerator_ = genEventInfo->weight();
    }
    
    // Create PDF weights
    if(includePdfWeights_ && !iEvent.isRealData()){
        edm::Handle<std::vector<double> > pdfWeights;
        iEvent.getByLabel(pdfWeightTag_, pdfWeights);
        v_pdfWeight_ = *pdfWeights;
    }
    
    // True-level vertices for pileup weighting
    if(!iEvent.isRealData()){
        edm::Handle<edm::View<PileupSummaryInfo> > pileupInfo;
        iEvent.getByLabel(pileupInfoTag_, pileupInfo);
        // Loop over pileup vector with size 3
        for(edm::View<PileupSummaryInfo>::const_iterator i_pileup = pileupInfo->begin(); i_pileup != pileupInfo->end(); ++i_pileup){
            // -1: previous BX, 0: current BX,  1: next BX
            if(i_pileup->getBunchCrossing() == 0) vertexMultiplicityTrue_ = i_pileup->getTrueNumInteractions();
        }
    }
    
    // Decay mode of ttbar system
    // FIXME: only for ttbar samples?
    edm::Handle<int> ttbarDecayMode;
    iEvent.getByLabel(ttbarDecayModeTag_, ttbarDecayMode);
    ttbarDecayMode_ = ttbarDecayMode.failedToGet() ? 0 : *ttbarDecayMode;
    
    // Generator info for ttbar samples
    if(isTtbarSample_){
        // Gen-level particles of ttbar system
        edm::Handle<TtGenEvent> ttbarGenEvent;
        iEvent.getByLabel(genEventTtbarTag_, ttbarGenEvent);
        if(!ttbarGenEvent.failedToGet()){
            // Which process generates the ttbar event: gluon-gluon-fusion (0), quark-quark-annihilation (1), all other processes (2)
            if(ttbarGenEvent->fromGluonFusion()) ttbarProductionMode_ = 0;
            else if(ttbarGenEvent->fromQuarkAnnihilation()) ttbarProductionMode_ = 1;
            else ttbarProductionMode_ = 2;
            
            // Top quarks, and for dileptonic decays also decay products
            if(ttbarGenEvent->top()) genTop_ = ttbarGenEvent->top()->polarP4();
            else genTop_ = nullP4_;
            if(ttbarGenEvent->topBar()) genAntiTop_ = ttbarGenEvent->topBar()->polarP4();
            else genAntiTop_ = nullP4_;
            if(ttbarGenEvent->lepton()){
                this->assignLeptonAndTau(ttbarGenEvent->lepton(), genLepton_, genLeptonPdgId_, genTau_);
            }
            else{
                genLepton_ = nullP4_;
                genLeptonPdgId_ = 0;
                genTau_ = nullP4_;
            }
            if(ttbarGenEvent->leptonBar()){
                this->assignLeptonAndTau(ttbarGenEvent->leptonBar(), genAntiLepton_, genAntiLeptonPdgId_, genAntiTau_);
            }
            else{
                genAntiLepton_ = nullP4_;
                genAntiLeptonPdgId_ = 0;
                genAntiTau_ = nullP4_;
            }
            if(ttbarGenEvent->b()) genB_ = ttbarGenEvent->b()->polarP4();
            else genB_ = nullP4_;
            if(ttbarGenEvent->bBar()) genAntiB_ = ttbarGenEvent->bBar()->polarP4();
            else genAntiB_ = nullP4_;
            if(ttbarGenEvent->neutrino()) genNeutrino_ = ttbarGenEvent->neutrino()->polarP4();
            else genNeutrino_ = nullP4_;
            if(ttbarGenEvent->neutrinoBar()) genAntiNeutrino_ = ttbarGenEvent->neutrinoBar()->polarP4();
            else genAntiNeutrino_ = nullP4_;
            if(ttbarGenEvent->wPlus()) genWPlus_ = ttbarGenEvent->wPlus()->polarP4();
            else genWPlus_ = nullP4_;
            if(ttbarGenEvent->wMinus()) genWMinus_ = ttbarGenEvent->wMinus()->polarP4();
            else genWMinus_ = nullP4_;
        }
        else{
            std::cerr<<"\nError: no ttbar gen event?!\n\n";
            ttbarProductionMode_ = -1;
            genTop_ = nullP4_;
            genAntiTop_ = nullP4_;
            genLepton_ = nullP4_;
            genAntiLepton_ = nullP4_;
            genTau_ = nullP4_;
            genAntiTau_ = nullP4_;
            genLeptonPdgId_ = 0;
            genAntiLeptonPdgId_ = 0;
            genB_ = nullP4_;
            genAntiB_ = nullP4_;
            genNeutrino_ = nullP4_;
            genAntiNeutrino_ = nullP4_;
            genWMinus_ = nullP4_;
            genWPlus_ = nullP4_;
        }
        
        // All gen jets
        edm::Handle<reco::GenJetCollection> genJets;
        iEvent.getByLabel(genJetsTag_, genJets);
        for(std::vector<reco::GenJet>::const_iterator i_jet = genJets->begin(); i_jet != genJets->end(); ++i_jet)
            v_allGenJet_.push_back(i_jet->polarP4());
        
        // Put full true info about genParticles
        //edm::Handle<std::vector<reco::GenParticle> > genParticles;
        //iEvent.getByLabel(genParticlesTag_, genParticles);
        //for(std::vector<reco::GenParticle>::const_iterator i_particle = genParticles->begin(); i_particle != genParticles->end(); ++i_particle){
        //    v_genParticleP4_.push_back(i_particle->polarP4());
        //    v_genParticlePdgId_.push_back(i_particle->pdgId());
        //    v_genParticleStatus_.push_back(i_particle->status());
        //}
        
        
        // Old quark-to-genJet association scheme
        // FIXME: should be removed one day from ntuple
        // CSA14: Workaround since this tool cannot work in miniAOD by setting label to "NONE"
        if(bHadJetIndexTag_.label() != "NONE"){
            // Find b jets corresponding to B hadrons
            edm::Handle<std::vector<int> > BHadJetIndex;
            iEvent.getByLabel(bHadJetIndexTag_, BHadJetIndex);
            for(size_t iJet = 0; iJet < BHadJetIndex->size(); ++iJet)
                v_bHadJetIndex_.push_back(BHadJetIndex->at(iJet));
            edm::Handle<std::vector<int> > AntiBHadJetIndex;
            iEvent.getByLabel(antiBHadJetIndexTag_, AntiBHadJetIndex);
            for(size_t iJet = 0; iJet < AntiBHadJetIndex->size(); ++iJet)
                v_antiBHadJetIndex_.push_back(AntiBHadJetIndex->at(iJet));
            
            // All B hadrons
            edm::Handle<std::vector<reco::GenParticle> > BHadrons;
            iEvent.getByLabel(BHadronsTag_, BHadrons);
            for(std::vector<reco::GenParticle>::const_iterator i_hadron = BHadrons->begin(); i_hadron != BHadrons->end(); ++i_hadron)
                v_bHadron_.push_back(i_hadron->polarP4());
            edm::Handle<std::vector<reco::GenParticle> > AntiBHadrons;
            iEvent.getByLabel(AntiBHadronsTag_, AntiBHadrons);
            for(std::vector<reco::GenParticle>::const_iterator i_hadron = AntiBHadrons->begin(); i_hadron != AntiBHadrons->end(); ++i_hadron)
                v_antiBHadron_.push_back(i_hadron->polarP4());
            
            // Which B hadrons stem from ttbar
            edm::Handle<std::vector<bool> > BHadronFromTopB;
            iEvent.getByLabel(BHadronFromTopBTag_, BHadronFromTopB);
            for(size_t iHadron = 0; iHadron < BHadronFromTopB->size(); ++iHadron)
                v_bHadFromTop_.push_back(BHadronFromTopB->at(iHadron));
            edm::Handle<std::vector<bool> > AntiBHadronFromTopB;
            iEvent.getByLabel(AntiBHadronFromTopBTag_, AntiBHadronFromTopB);
            for(size_t iHadron = 0; iHadron < AntiBHadronFromTopB->size(); ++iHadron)
                v_antiBHadFromTop_.push_back(AntiBHadronFromTopB->at(iHadron));
            
            // Another B hadron to jet association ?
            edm::Handle<std::vector<int> > BHadronVsJet;
            iEvent.getByLabel(BHadronVsJetTag_, BHadronVsJet);
            for(size_t iHadron = 0; iHadron < BHadronVsJet->size(); ++iHadron)
                v_bHadVsJet_.push_back(BHadronVsJet->at(iHadron));
            edm::Handle<std::vector<int> > AntiBHadronVsJet;
            iEvent.getByLabel(AntiBHadronVsJetTag_, AntiBHadronVsJet);
            for(size_t iHadron = 0; iHadron < AntiBHadronVsJet->size(); ++iHadron)
                v_antiBHadVsJet_.push_back(AntiBHadronVsJet->at(iHadron));
        }
    }
    
    // Generator info for Higgs samples
    if(isHiggsSample_){
        
        // Gen-level particles of Higgs decay
        edm::Handle<HiggsGenEvent> genEventHiggs;
        iEvent.getByLabel(genEventHiggsTag_, genEventHiggs);
        if(!genEventHiggs.failedToGet()){
            if(genEventHiggs->higgs()) genH_ = genEventHiggs->higgs()->polarP4();
            else genH_ = nullP4_;
            if(genEventHiggs->b()) genBFromH_ = genEventHiggs->b()->polarP4();
            else genBFromH_ = nullP4_;
            if(genEventHiggs->bBar()) genAntiBFromH_ = genEventHiggs->bBar()->polarP4();
            else genAntiBFromH_ = nullP4_;
        }
        else{
            std::cerr<<"\nError: no Higgs gen event?!\n\n";
            genH_ = nullP4_;
            genBFromH_ = nullP4_;
            genAntiBFromH_ = nullP4_;
        }
        
        // Decay mode
        edm::Handle<int> higgsDecayMode;
        iEvent.getByLabel(higgsDecayModeTag_, higgsDecayMode);
        higgsDecayMode_ = higgsDecayMode.failedToGet() ? 0 : *higgsDecayMode;
    }
    
    // Generator info for Z samples
    if(isZSample_){
        // Gen-level particles and decay mode
        edm::Handle<std::vector<GenZDecayProperties> > v_genZDecayProperties;
        iEvent.getByLabel(genZDecayTag_, v_genZDecayProperties);
        if(!v_genZDecayProperties.failedToGet()){
            for(size_t i = 0; i < v_genZDecayProperties->size(); ++i){
                v_genZ_.push_back(v_genZDecayProperties->at(i).z()->polarP4());
                v_genZMeDaughterParticle_.push_back(v_genZDecayProperties->at(i).meDaughterParticle()->polarP4());
                v_genZMeDaughterAntiParticle_.push_back(v_genZDecayProperties->at(i).meDaughterAntiParticle()->polarP4());
                if(v_genZDecayProperties->at(i).stableLepton()) v_genZStableLepton_.push_back(v_genZDecayProperties->at(i).stableLepton()->polarP4());
                else v_genZStableLepton_.push_back(nullP4_);
                if(v_genZDecayProperties->at(i).stableAntiLepton()) v_genZStableAntiLepton_.push_back(v_genZDecayProperties->at(i).stableAntiLepton()->polarP4());
                else v_genZStableAntiLepton_.push_back(nullP4_);
                v_genZDecayMode_.push_back(v_genZDecayProperties->at(i).decayMode());
            }
        }
        else{
            std::cerr<<"\nError: no GenZDecayProperties ?!\n\n";
        }
    }
    
    // W decay modes for MadGraph samples
    if(isMadgraphSample_){
        edm::Handle<std::vector<int> > madgraphWDecays;
        iEvent.getByLabel(madgraphWDecayTag_, madgraphWDecays);
        v_madgraphWDecay_ = *(madgraphWDecays.product());
    }
    
    // Improved quark-to-genJet association scheme
    // Only stored for samples containing ttbar generator info
    if(isTtbarSample_){
        // Signal definition cuts for gen b jets
        constexpr double genSignalJetPt_min = 20.0;
        constexpr double genSignalJetEta_max = 2.4;
        
        // List of jets containing B hadrons (not) coming from the top weak decay for event categorisation
        std::vector<int> genBJetNotFromTopDecayChainIndex;
        std::vector<int> genBJetFromTopDecayChainIndex;
        std::set<int> genBJetFromTopIds_unique;
        
        edm::Handle<std::vector<int> > genBHadIndex;
        iEvent.getByLabel(genBHadIndexTag_, genBHadIndex);
        if(!genBHadIndex.failedToGet()){
            v_genBHadIndex_ = *(genBHadIndex.product());
            
            edm::Handle<std::vector<std::vector<int> > > genBHadPlusMothersIndices;
            iEvent.getByLabel(genBHadPlusMothersIndicesTag_, genBHadPlusMothersIndices);
            // Only if all hadron mothers have to be stored
            if(saveHadronMothers_)
                v_genBHadPlusMothersIndices_ = *(genBHadPlusMothersIndices.product());
            
            edm::Handle<std::vector<int> > genBHadLeptonIndex;
            iEvent.getByLabel(genBHadLeptonIndexTag_, genBHadLeptonIndex);
            v_genBHadLeptonIndex_ = *(genBHadLeptonIndex.product());
            
            edm::Handle<std::vector<reco::GenParticle> > genBHadPlusMothers;
            iEvent.getByLabel(genBHadPlusMothersTag_, genBHadPlusMothers);
            // If all particles have to be stored
            if(saveHadronMothers_){
                for(std::vector<reco::GenParticle>::const_iterator i_particle = genBHadPlusMothers->begin(); i_particle != genBHadPlusMothers->end(); ++i_particle){
                    v_genBHadPlusMothers_.push_back(i_particle->polarP4());
                    v_genBHadPlusMothersPdgId_.push_back(i_particle->pdgId());
                    v_genBHadPlusMothersStatus_.push_back(i_particle->status());
                }
            }
            // If only hadrons/leptons have to be stored
            else{
                for(size_t iParticle = 0; iParticle<genBHadIndex->size(); ++iParticle){
                    v_genBHadPlusMothers_.push_back(genBHadPlusMothers->at(genBHadIndex->at(iParticle)).polarP4());
                    v_genBHadPlusMothersPdgId_.push_back(genBHadPlusMothers->at(genBHadIndex->at(iParticle)).pdgId());
                    v_genBHadIndex_.at(iParticle) = v_genBHadPlusMothers_.size() - 1;
                }
                for(size_t iParticle = 0; iParticle < genBHadLeptonIndex->size(); ++iParticle){
                    v_genBHadPlusMothers_.push_back(genBHadPlusMothers->at(genBHadLeptonIndex->at(iParticle)).polarP4());
                    v_genBHadPlusMothersPdgId_.push_back(genBHadPlusMothers->at(genBHadLeptonIndex->at(iParticle)).pdgId());
                    v_genBHadLeptonIndex_.at(iParticle) = v_genBHadPlusMothers_.size() - 1;
                }
            }
            
            edm::Handle<std::vector<int> > genBHadLeptonHadronIndex;
            iEvent.getByLabel(genBHadLeptonHadronIndexTag_, genBHadLeptonHadronIndex);
            v_genBHadLeptonHadronIndex_ = *(genBHadLeptonHadronIndex.product());
            
            // Detecting whether leptons come via hadron->tau->lepton
            edm::Handle<std::vector<int> > genBHadLeptonViaTau;
            iEvent.getByLabel(genBHadLeptonViaTauTag_, genBHadLeptonViaTau);
            if(!genBHadLeptonViaTau.failedToGet()) {
                v_genBHadLeptonViaTau_ = *(genBHadLeptonViaTau.product());
            } else {
                for(size_t iParticle = 0; iParticle < genBHadLeptonIndex->size(); ++iParticle){
                    int viaTau = -1;
                    const int leptonMotherIndex = genBHadPlusMothersIndices->at(genBHadLeptonIndex->at(iParticle)).at(0);
                    if(leptonMotherIndex >= 0) 
                        viaTau = std::abs(genBHadPlusMothers->at(leptonMotherIndex).pdgId())==15 ? 1 : 0;
                    v_genBHadLeptonViaTau_.push_back(viaTau);
                }
            }
            
            edm::Handle<std::vector<int> > genBHadFlavour;
            iEvent.getByLabel(genBHadFlavourTag_, genBHadFlavour);
            v_genBHadFlavour_ = *(genBHadFlavour.product());
            
            edm::Handle<std::vector<int> > genBHadJetIndex;
            iEvent.getByLabel(genBHadJetIndexTag_, genBHadJetIndex);
            unsigned int iHadron = 0;
            for(std::vector<int>::const_iterator i_index = genBHadJetIndex->begin(); i_index != genBHadJetIndex->end(); ++i_index){
                ++iHadron;
                v_genBHadJetIndex_.push_back(*i_index);
                // Build a vector of b jets from top
                if(v_genBHadFlavour_.size() < iHadron) continue;
                if(std::abs(v_genBHadFlavour_.at(iHadron-1)) != 6) continue;
                genBJetFromTopIds_unique.insert(*i_index);
            }
            
            edm::Handle<std::vector<int> > genBHadFromTopWeakDecay;
            iEvent.getByLabel(genBHadFromTopWeakDecayTag_, genBHadFromTopWeakDecay);
            v_genBHadFromTopWeakDecay_ = *(genBHadFromTopWeakDecay.product());
            
            
            // Select jets from different sources for event categorisation
            for(size_t index = 0; index < genBHadJetIndex->size(); ++index){
                if(v_genBHadFromTopWeakDecay_.size() <= index) break;
                // Skip B hadrons from top
                if(genBHadFlavour->size()>index && std::abs(genBHadFlavour->at(index))==6) continue;
                const int jetId = genBHadJetIndex->at(index);
                if(jetId < 0) continue;
                // Check that jet satisfies the signal definition
                if(v_allGenJet_.at(jetId).Pt() < genSignalJetPt_min) continue;
                if(std::fabs(v_allGenJet_.at(jetId).Eta()) > genSignalJetEta_max) continue;
                // Skip jet if already identified as b jet from top
                if(genBJetFromTopIds_unique.count(jetId) > 0) continue;

                if(v_genBHadFromTopWeakDecay_.at(index) == 0)
                    genBJetNotFromTopDecayChainIndex.push_back(jetId);
                else if(v_genBHadFromTopWeakDecay_.at(index) == 1)
                    genBJetFromTopDecayChainIndex.push_back(jetId);
            }
        }
        
        // List of jets containing C hadrons (not) coming from the top weak decay for event categorisation
        std::vector<int> genCJetNotFromTopDecayChainIndex;
        std::vector<int> genCJetFromTopDecayChainIndex;
        
        edm::Handle<std::vector<int> > genCHadIndex;
        iEvent.getByLabel(genCHadIndexTag_, genCHadIndex);
        if(!genCHadIndex.failedToGet()){
            v_genCHadIndex_ = *(genCHadIndex.product());
            
            edm::Handle<std::vector<int> > genCHadLeptonIndex;
            iEvent.getByLabel(genCHadLeptonIndexTag_, genCHadLeptonIndex);
            v_genCHadLeptonIndex_ = *(genCHadLeptonIndex.product());
            
            edm::Handle<std::vector<reco::GenParticle> > genCHadPlusMothers;
            iEvent.getByLabel(genCHadPlusMothersTag_, genCHadPlusMothers);
            edm::Handle<std::vector<std::vector<int> > > genCHadPlusMothersIndices;
            iEvent.getByLabel(genCHadPlusMothersIndicesTag_, genCHadPlusMothersIndices);
            // If all particles have to be stored
            if(saveHadronMothers_){
                for(std::vector<reco::GenParticle>::const_iterator i_particle = genCHadPlusMothers->begin(); i_particle != genCHadPlusMothers->end(); ++i_particle){
                    v_genCHadPlusMothers_.push_back(i_particle->polarP4());
                    v_genCHadPlusMothersPdgId_.push_back(i_particle->pdgId());
                    v_genCHadPlusMothersStatus_.push_back(i_particle->status());
                }
            }
            // If only hadrons/leptons have to be stored
            else{
                for(unsigned int iParticle = 0; iParticle < genCHadIndex->size(); ++iParticle){
                    v_genCHadIndex_.at(iParticle) = -1;
                    // Do not store actual c-hadrons but store leptons
                    if(!saveCHadronParticles_) continue;
                    v_genCHadPlusMothers_.push_back(genCHadPlusMothers->at(genCHadIndex->at(iParticle)).polarP4());
                    v_genCHadPlusMothersPdgId_.push_back(genCHadPlusMothers->at(genCHadIndex->at(iParticle)).pdgId());
                    v_genCHadIndex_.at(iParticle) = v_genCHadPlusMothers_.size() - 1;
                }
                for(unsigned int iParticle = 0; iParticle < genCHadLeptonIndex->size(); ++iParticle){
                    v_genCHadPlusMothers_.push_back(genCHadPlusMothers->at(genCHadLeptonIndex->at(iParticle)).polarP4());
                    v_genCHadPlusMothersPdgId_.push_back(genCHadPlusMothers->at(genCHadLeptonIndex->at(iParticle)).pdgId());
                    v_genCHadLeptonIndex_.at(iParticle) = v_genCHadPlusMothers_.size() - 1;
                }
            }
            
            edm::Handle<std::vector<int> > genCHadLeptonHadronIndex;
            iEvent.getByLabel(genCHadLeptonHadronIndexTag_, genCHadLeptonHadronIndex);
            v_genCHadLeptonHadronIndex_ = *(genCHadLeptonHadronIndex.product());
            
            // Detecting whether leptons come via hadron->tau->lepton
            edm::Handle<std::vector<int> > genCHadLeptonViaTau;
            iEvent.getByLabel(genCHadLeptonViaTauTag_, genCHadLeptonViaTau);
            if(!genCHadLeptonViaTau.failedToGet()) {
                v_genCHadLeptonViaTau_ = *(genCHadLeptonViaTau.product());
            } else {
                for(size_t iParticle = 0; iParticle < genCHadLeptonIndex->size(); ++iParticle){
                    int viaTau = -1;
                    const int leptonMotherIndex = genCHadPlusMothersIndices->at(genCHadLeptonIndex->at(iParticle)).at(0);
                    if(leptonMotherIndex >= 0) 
                        viaTau = std::abs(genCHadPlusMothers->at(leptonMotherIndex).pdgId())==15 ? 1 : 0;
                    v_genCHadLeptonViaTau_.push_back(viaTau);
                }
            }
            
            edm::Handle<std::vector<int> > genCHadJetIndex;
            iEvent.getByLabel(genCHadJetIndexTag_, genCHadJetIndex);
            v_genCHadJetIndex_ = *(genCHadJetIndex.product());
            
            edm::Handle<std::vector<int> > genCHadFromTopWeakDecay;
            iEvent.getByLabel(genCHadFromTopWeakDecayTag_, genCHadFromTopWeakDecay);
            std::vector<int> v_genCHadFromTopWeakDecay;
            v_genCHadFromTopWeakDecay = *(genCHadFromTopWeakDecay.product());
            
            // Convert index of parent B hadron to point to the specifically identified collection of B hadrons
            // -1  - C hadron doesn't come from a B hadron'
            // -2  - C hadron comes from a B hadron that was not identified in B hadron specific run
            //       [comes from a non-weakly decaying B hadron, decay products of the B hadron not clustered to any jet]
            edm::Handle<std::vector<int> > genCHadBHadronId;
            iEvent.getByLabel(genCHadBHadronIdTag_, genCHadBHadronId);
            for(std::vector<int>::const_iterator i_id = genCHadBHadronId->begin(); i_id!=genCHadBHadronId->end(); ++i_id){
                const int bHadId = *i_id;
                // ID for B hadrons as mothers of C hadrons is different than: ID directly from the B hadrons, as written to ntuple
                // Needs to be corrected
                int bHadIdConverted = -1;
                if(bHadId >= 0){
                    // Look for the same particle in list of hadrons from B hadron related set of particles
                    for(std::vector<int>::const_iterator i_index = v_genBHadIndex_.begin(); i_index != v_genBHadIndex_.end(); ++i_index){
                        if(genCHadPlusMothers->at(bHadId).pdgId() != v_genBHadPlusMothersPdgId_.at(*i_index)) continue;
                        if(std::fabs(genCHadPlusMothers->at(bHadId).pt() - v_genBHadPlusMothers_.at(*i_index).Pt()) > 1e-9) continue;
                        bHadIdConverted = i_index - v_genBHadIndex_.begin();
                        break;
                    }
                    if(bHadIdConverted < 0) bHadIdConverted = -2;
                }
                v_genCHadBHadronId_.push_back(bHadIdConverted);
            }
            
            // Select jets from different sources for event categorisation
            for(size_t index = 0; index < genCHadJetIndex->size(); ++index){
                if(v_genCHadFromTopWeakDecay.size() <= index) break;
                // Skip C hadrons coming from B hadrons
                if(genCHadBHadronId->size()>index && genCHadBHadronId->at(index)>=0) continue;
                const int jetId = genCHadJetIndex->at(index);
                if(jetId < 0) continue;
                // Check that jet satisfies the signal definition
                if(v_allGenJet_.at(jetId).Pt() < genSignalJetPt_min) continue;
                if(std::fabs(v_allGenJet_.at(jetId).Eta()) > genSignalJetEta_max) continue;
                // Skip jet if already identified as b jet from top
                if(genBJetFromTopIds_unique.count(jetId) > 0) continue;
                if(v_genCHadFromTopWeakDecay.at(index) == 0)
                    genCJetNotFromTopDecayChainIndex.push_back(jetId);
                else if(v_genCHadFromTopWeakDecay.at(index) == 1)
                    genCJetFromTopDecayChainIndex.push_back(jetId);
            }
        }
        
        // Id of the event depending on the number of additional b(c) jets in the event
        // -1  : Something went wrong
        // 204 : 2 b jets from top & >=2 extra b jets not from the top weak decay chain with <2 jets having 1 hadron/jet (+ any number of c/light jets)
        // 203 : 2 b jets from top & >=2 extra b jets not from the top weak decay chain with >=2 jets having 1 hadron/jet (+ any number of c/light jets)
        // 202 : 2 b jets from top & exactly 1 extra b jet not from the top weak decay chain with >1 hadrons/jet (+ any number of c/light jets)
        // 201 : 2 b jets from top & exactly 1 extra b jet not from the top weak decay chain containing 1 b hadron (+ any number of c/light jets)
        // 214 : 2 b jets from top & >=2 extra b jets from the top weak decay chain with <2 jets having 1 hadron/jet (+ any number of c/light jets)
        // 213 : 2 b jets from top & >=2 extra b jets from the top weak decay chain with >=2 jets having 1 hadron/jet (+ any number of c/light jets)
        // 212 : 2 b jets from top & exactly 1 extra b jet from the top weak decay chain with >1 hadrons/jet (+ any number of c/light jets)
        // 211 : 2 b jets from top & exactly 1 extra b jet from the top weak decay chain containing 1 b hadron (+ any number of c/light jets)
        // 224 : 2 b jets from top & no extra b jets & >=2 c jets not from the top weak decay chain with <2 jets having 1 hadron/jet (+ any number of light jets)
        // 223 : 2 b jets from top & no extra b jets & >=2 c jets not from the top weak decay chain with >=2 jets having 1 hadron/jet (+ any number of light jets)
        // 222 : 2 b jets from top & no extra b jets & >=2 c jets not from the top weak decay chain with >1 hadrons/jet (+ any number of light jets)
        // 221 : 2 b jets from top & no extra b jets & exactly 1 c jet not from the top weak decay chain containing 1 b hadron (+ any number of light jets)
        // 234 : 2 b jets from top & no extra b jets & >=2 c jets from the top weak decay chain with <2 jets having 1 hadron/jet (+ any number of light jets)
        // 233 : 2 b jets from top & no extra b jets & >=2 c jets from the top weak decay chain with >=2 jets having 1 hadron/jet (+ any number of light jets)
        // 232 : 2 b jets from top & no extra b jets & >=2 c jets from the top weak decay chain with >1 hadrons/jet (+ any number of light jets)
        // 231 : 2 b jets from top & no extra b jets & exactly 1 c jet from the top weak decay chain containing 1 b hadron (+ any number of light jets)
        // 200 : 2 b jets from top & no extra b jets that come not directly from top & no c jets that come not from B hadrons
        // 104 : 1 b jet from top & >=2 extra b jets not from the top weak decay chain with <2 jets having 1 hadron/jet (+ any number of c/light jets)
        //  - - - - - - - - - - - - - - 
        // 100 : 1 b jet from top & no extra b jets that come not directly from top & no c jets that come not from B hadrons
        genExtraTopJetNumberId_ = -1;
            
        const int nGenBJetNotFromTopDecayChain = nUniqueElementsInVector(genBJetNotFromTopDecayChainIndex);
        const int nGenBJetNotFromTopDecayChain_unique = nUniqueElementsInVector(genBJetNotFromTopDecayChainIndex, true);
        const int nGenBJetFromTopDecayChain = nUniqueElementsInVector(genBJetFromTopDecayChainIndex);
        const int nGenBJetFromTopDecayChain_unique = nUniqueElementsInVector(genBJetFromTopDecayChainIndex, true);
        const int nGenCJetNotFromTopDecayChain = nUniqueElementsInVector(genCJetNotFromTopDecayChainIndex);
        const int nGenCJetNotFromTopDecayChain_unique = nUniqueElementsInVector(genCJetNotFromTopDecayChainIndex, true);
        const int nGenCJetFromTopDecayChain = nUniqueElementsInVector(genCJetFromTopDecayChainIndex);
        const int nGenCJetFromTopDecayChain_unique = nUniqueElementsInVector(genCJetFromTopDecayChainIndex, true);
        const int genTopJetNumberId = 100*genBJetFromTopIds_unique.size();
        if(nGenBJetNotFromTopDecayChain >= 2) {
            if(nGenBJetNotFromTopDecayChain_unique >= 2) genExtraTopJetNumberId_ = 3 + genTopJetNumberId;
            else genExtraTopJetNumberId_ = 4 + genTopJetNumberId;
        }
        else if(nGenBJetNotFromTopDecayChain == 1) {
            if(nGenBJetNotFromTopDecayChain_unique >= 1) genExtraTopJetNumberId_ = 1 + genTopJetNumberId;
            else genExtraTopJetNumberId_ = 2 + genTopJetNumberId;
        }
        else if(nGenBJetFromTopDecayChain >= 2) {
            if(nGenBJetFromTopDecayChain_unique >= 2) genExtraTopJetNumberId_ = 13 + genTopJetNumberId;
            else genExtraTopJetNumberId_ = 14 + genTopJetNumberId;
        }
        else if(nGenBJetFromTopDecayChain == 1) {
            if(nGenBJetFromTopDecayChain_unique >= 1) genExtraTopJetNumberId_ = 11 + genTopJetNumberId;
            else genExtraTopJetNumberId_ = 12 + genTopJetNumberId;
        }
        else if(nGenCJetNotFromTopDecayChain >= 2) {
            if(nGenCJetNotFromTopDecayChain_unique >= 2) genExtraTopJetNumberId_ = 23 + genTopJetNumberId;
            else genExtraTopJetNumberId_ = 24 + genTopJetNumberId;
        }
        else if(nGenCJetNotFromTopDecayChain == 1) {
            if(nGenCJetNotFromTopDecayChain_unique >= 1) genExtraTopJetNumberId_ = 21 + genTopJetNumberId;
            else genExtraTopJetNumberId_ = 22 + genTopJetNumberId;
        }
        else if(nGenCJetFromTopDecayChain >= 2) {
            if(nGenCJetFromTopDecayChain_unique >= 2) genExtraTopJetNumberId_ = 33 + genTopJetNumberId;
            else genExtraTopJetNumberId_ = 34 + genTopJetNumberId;
        }
        else if(nGenCJetFromTopDecayChain == 1) {
            if(nGenCJetFromTopDecayChain_unique >= 1) genExtraTopJetNumberId_ = 31 + genTopJetNumberId;
            else genExtraTopJetNumberId_ = 32 + genTopJetNumberId;
        }
        else genExtraTopJetNumberId_ = 0 + genTopJetNumberId;
    }
    
    
    
    // Triggers
    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByLabel(triggerResultsTag_, triggerResults);
    triggerBits_ = this->getTriggerBits(iEvent, triggerResults);
    triggerBitsTau_ = this->getTriggerBitsTau(iEvent, triggerResults);
    
    // Reco-level vertices
    edm::Handle<std::vector<reco::Vertex> > vertices;
    iEvent.getByLabel(verticesTag_, vertices);
    vertexMultiplicity_ = vertices->size();
    
    // Leptons (order by pt)
    edm::Handle<std::vector<pat::Electron> > electrons;
    iEvent.getByLabel(electronsTag_, electrons);
    std::vector<pat::Electron>::const_iterator i_electron  = electrons->begin();
    edm::Handle<std::vector<pat::Muon> > muons;
    iEvent.getByLabel(muonsTag_, muons);
    std::vector<pat::Muon>::const_iterator i_muon  = muons->begin();
    while(i_muon!=muons->end() || i_electron!=electrons->end()){
        // Check pt order: electrons and muons are ordered each, but for common list common ordering is necessary
        bool writeElectron = false;
        bool writeMuon = false;
        if(i_electron == electrons->end()) writeMuon = true;
        else if(i_muon == muons->end()) writeElectron = true;
        else if(i_muon->pt() > i_electron->pt()) writeMuon = true;
        else writeElectron = true;
        
        // Fill electron stuff
        if(writeElectron){
            v_lepton_.push_back(i_electron->polarP4());
            v_leptonPdgId_.push_back(i_electron->pdgId());
            
            const reco::GsfTrack& track = *(i_electron->gsfTrack());
            if(vertexMultiplicity_){
                v_leptonDxyVertex0_.push_back(track.dxy(vertices->at(0).position()));
                v_leptonDzVertex0_.push_back(track.dz(vertices->at(0).position()));
            }
            else{
                v_leptonDxyVertex0_.push_back(-999.);
                v_leptonDzVertex0_.push_back(-999.);
            }
            
            // Electron MVA ID values
            double idTemp = -9999.;
            std::vector<std::pair<std::string,float> > electronMVAIDs = i_electron->electronIDs();
            for(size_t id = 0; id < electronMVAIDs.size(); ++id){
                if(electronMVAIDs[id].first == "mvaTrigV0"){
                    idTemp = electronMVAIDs[id].second;
                    break;
                }
            }
            v_leptonId_.push_back(idTemp);
            
            // Isolation
            v_leptonChargedHadronIso_.push_back(i_electron->chargedHadronIso());
            v_leptonNeutralHadronIso_.push_back(i_electron->neutralHadronIso());
            v_leptonPhotonIso_.push_back(i_electron->photonIso());
            v_leptonPuChargedHadronIso_.push_back(i_electron->puChargedHadronIso());
            v_leptonPfIso_.push_back(
                (i_electron->chargedHadronIso()
                    + std::max(0., i_electron->neutralHadronIso() + i_electron->photonIso() - 0.5*i_electron->puChargedHadronIso())
                ) / i_electron->pt());
            // Barrel region
            if(std::fabs(i_electron->superCluster()->eta()) <= 1.479)
                v_leptonCombIso_.push_back((i_electron->dr03TkSumPt() + std::max(0., i_electron->dr03EcalRecHitSumEt()-1.) + i_electron->dr03HcalTowerSumEt())/std::max<double>(20., i_electron->pt()));
            // Endcap region
            else
                v_leptonCombIso_.push_back((i_electron->dr03TkSumPt() + i_electron->dr03EcalRecHitSumEt() + i_electron->dr03HcalTowerSumEt())/std::max<double>(20., i_electron->pt()));
            
            // Trigger matches
            int triggerResult = 0;
            const pat::TriggerObjectStandAloneCollection& triggers = i_electron->triggerObjectMatches();
            if(triggers.size() > 0) triggerResult = triggerMap_["general"];
            for(unsigned int iTrigger = 0; iTrigger < triggers.size(); ++iTrigger)
                triggerResult |= this->getTriggerBits(triggers.at(iTrigger).pathNames());
            v_leptonTrigger_.push_back(triggerResult);
            
            ++i_electron;
        }
        
        // Fill muon stuff
        if(writeMuon){
            v_lepton_.push_back(i_muon->polarP4());
            v_leptonPdgId_.push_back(i_muon->pdgId());
            v_leptonId_.push_back(-1);
            
            // Following https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId, these functions should be used for dxy and dz
            v_leptonDxyVertex0_.push_back(-i_muon->dB());
            if(i_muon->muonBestTrack().isAvailable() && vertexMultiplicity_) v_leptonDzVertex0_.push_back(i_muon->muonBestTrack()->dz(vertices->at(0).position()));
            else v_leptonDzVertex0_.push_back(-999.);
            
            // Isolation
            v_leptonChargedHadronIso_.push_back(i_muon->chargedHadronIso());
            v_leptonNeutralHadronIso_.push_back(i_muon->neutralHadronIso());
            v_leptonPhotonIso_.push_back(i_muon->photonIso());
            v_leptonPuChargedHadronIso_.push_back(i_muon->puChargedHadronIso());
            v_leptonPfIso_.push_back(
                (i_muon->chargedHadronIso()
                    + std::max(0., i_muon->neutralHadronIso() + i_muon->photonIso() - 0.5*i_muon->puChargedHadronIso())
                ) / i_muon->pt());
            v_leptonCombIso_.push_back((i_muon->trackIso() + i_muon->caloIso()) / i_muon->pt());
            
            // Trigger matches
            int triggerResult = 0;
            const pat::TriggerObjectStandAloneCollection& triggers = i_muon->triggerObjectMatches();
            if(triggers.size() > 0) triggerResult = triggerMap_["general"];
            for(size_t iTrigger = 0; iTrigger < triggers.size(); ++iTrigger)
                triggerResult |= this->getTriggerBits(triggers.at(iTrigger).pathNames());
            v_leptonTrigger_.push_back(triggerResult);
            
            ++i_muon;
        }
    }
    
    // Jets
    edm::Handle<edm::View<pat::Jet> > jets;
    iEvent.getByLabel(jetsTag_, jets);
    for(edm::View<pat::Jet>::const_iterator i_jet  = jets->begin(); i_jet != jets->end(); ++i_jet){
        v_jet_.push_back(i_jet->polarP4());
        if(!iEvent.isRealData()){
            v_jetJerSF_.push_back(i_jet->userFloat("jerSF"));
            v_jetPartonFlavour_.push_back(i_jet->partonFlavour());
            if(i_jet->genJet()) v_associatedGenJet_.push_back(i_jet->genJet()->polarP4());
            else v_associatedGenJet_.push_back(nullP4_);
        }
        
        // Discriminators for b-tagging
        v_jetBtagTCHE_.push_back(i_jet->bDiscriminator("trackCountingHighEffBJetTags"));
        v_jetBtagTCHP_.push_back(i_jet->bDiscriminator("trackCountingHighPurBJetTags"));
        v_jetBtagJetProbability_.push_back(i_jet->bDiscriminator("jetProbabilityBJetTags"));
        v_jetBtagJetBProbability_.push_back(i_jet->bDiscriminator("jetBProbabilityBJetTags"));
        v_jetBtagSSVHE_.push_back(i_jet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));
        v_jetBtagSSVHP_.push_back(i_jet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
        v_jetBtagCSV_.push_back(i_jet->bDiscriminator("combinedSecondaryVertexBJetTags"));
        v_jetBtagCSVMVA_.push_back(i_jet->bDiscriminator("combinedSecondaryVertexMVABJetTags"));
    }
    
    // Separate jet collection needed for the on-the-fly calculation of jet uncertainties,
    // because even bad-id jets are used for MET
    constexpr double jetPtThresholdForMet = 10.;
    constexpr double jetEmLimitForMet = 0.9;
    if(!iEvent.isRealData()){
        edm::Handle<edm::View<pat::Jet> > jetsForMet;
        iEvent.getByLabel(jetsForMetTag_, jetsForMet);
        edm::Handle<edm::View<pat::Jet> > jetsForMetUncorrected;
        iEvent.getByLabel(jetsForMetUncorrectedTag_, jetsForMetUncorrected);
        for(size_t iJet = 0; iJet < jetsForMetUncorrected->size(); ++iJet){
            if(jetsForMetUncorrected->at(iJet).correctedJet("Uncorrected").pt() > jetPtThresholdForMet &&
                ((!jetsForMetUncorrected->at(iJet).isPFJet() && jetsForMetUncorrected->at(iJet).emEnergyFraction() < jetEmLimitForMet) ||
                 (jetsForMetUncorrected->at(iJet).isPFJet() && jetsForMetUncorrected->at(iJet).neutralEmEnergyFraction() + jetsForMetUncorrected->at(iJet).chargedEmEnergyFraction() < jetEmLimitForMet))){
                v_jetForMet_.push_back(jetsForMet->at(iJet).polarP4());
                v_jetForMetJerSF_.push_back(jetsForMet->at(iJet).userFloat("jerSF"));
                v_jetPartonFlavourForMet_.push_back(jetsForMet->at(iJet).partonFlavour());
                if(jetsForMet->at(iJet).genJet()) v_associatedGenJetForMet_.push_back(jetsForMet->at(iJet).genJet()->polarP4());
                else v_associatedGenJetForMet_.push_back(nullP4_);
            }
        }
    }
    
    // Jet properties concerning constituents and secondary vertices
    // CSA14: Workaround as long as jetProperties are not working by setting label to "NONE"
    if(jetPropertiesTag_.label() != "NONE"){
        edm::Handle<std::vector<JetProperties> > v_jetProperties;
        iEvent.getByLabel(jetPropertiesTag_, v_jetProperties);
        for(std::vector<JetProperties>::const_iterator i_jetProperties = v_jetProperties->begin(); i_jetProperties != v_jetProperties->end(); ++i_jetProperties){
            v_jetChargeGlobalPtWeighted_.push_back(i_jetProperties->jetChargeGlobalPtWeighted());
            v_jetChargeRelativePtWeighted_.push_back(i_jetProperties->jetChargeRelativePtWeighted());
            if(isTtbarSample_){
                v_jetAssociatedPartonPdgId_.push_back(i_jetProperties->jetAssociatedPartonPdgId());
                v_jetAssociatedParton_.push_back(i_jetProperties->jetAssociatedParton());
            }
            v_jetSecondaryVertexPtCorrectedMass_.push_back(i_jetProperties->jetSecondaryVertexPtCorrectedMass());
            
            const int jetIndex = i_jetProperties - v_jetProperties->begin();
            
            v_jetPfCandidateTrack_.insert(v_jetPfCandidateTrack_.end(), i_jetProperties->jetPfCandidateTrack().begin(), i_jetProperties->jetPfCandidateTrack().end());
            v_jetPfCandidateTrackCharge_.insert(v_jetPfCandidateTrackCharge_.end(), i_jetProperties->jetPfCandidateTrackCharge().begin(), i_jetProperties->jetPfCandidateTrackCharge().end());
            v_jetPfCandidateTrackId_.insert(v_jetPfCandidateTrackId_.end(), i_jetProperties->jetPfCandidateTrackId().begin(), i_jetProperties->jetPfCandidateTrackId().end());
            v_jetPfCandidatePrimaryVertexId_.insert(v_jetPfCandidatePrimaryVertexId_.end(),i_jetProperties->jetPfCandidatePrimaryVertexId().begin(), i_jetProperties->jetPfCandidatePrimaryVertexId().end());
            v_jetPfCandidateTrackIndex_.insert(v_jetPfCandidateTrackIndex_.end(), i_jetProperties->jetPfCandidateTrack().size(), jetIndex);
            
            v_jetSelectedTrack_.insert(v_jetSelectedTrack_.end(), i_jetProperties->jetSelectedTrack().begin(), i_jetProperties->jetSelectedTrack().end());
            v_jetSelectedTrackIPValue_.insert(v_jetSelectedTrackIPValue_.end(), i_jetProperties->jetSelectedTrackIPValue().begin(), i_jetProperties->jetSelectedTrackIPValue().end());
            v_jetSelectedTrackIPSignificance_.insert(v_jetSelectedTrackIPSignificance_.end(), i_jetProperties->jetSelectedTrackIPSignificance().begin(), i_jetProperties->jetSelectedTrackIPSignificance().end());
            v_jetSelectedTrackCharge_.insert(v_jetSelectedTrackCharge_.end(), i_jetProperties->jetSelectedTrackCharge().begin(), i_jetProperties->jetSelectedTrackCharge().end());
            v_jetSelectedTrackMatchToPfCandidateIndex_.insert(v_jetSelectedTrackMatchToPfCandidateIndex_.end(), i_jetProperties->jetSelectedTrackMatchToPfCandidateIndex().begin(), i_jetProperties->jetSelectedTrackMatchToPfCandidateIndex().end());
            v_jetSelectedTrackIndex_.insert(v_jetSelectedTrackIndex_.end(), i_jetProperties->jetSelectedTrack().size(), jetIndex);
            
            v_jetSecondaryVertexTrackMatchToSelectedTrackIndex_.insert(v_jetSecondaryVertexTrackMatchToSelectedTrackIndex_.end(), i_jetProperties->jetSecondaryVertexTrackMatchToSelectedTrackIndex().begin(), i_jetProperties->jetSecondaryVertexTrackMatchToSelectedTrackIndex().end());
            v_jetSecondaryVertexTrackVertexIndex_.insert(v_jetSecondaryVertexTrackVertexIndex_.end(), i_jetProperties->jetSecondaryVertexTrackVertexIndex().begin(), i_jetProperties->jetSecondaryVertexTrackVertexIndex().end());
            
            v_jetSecondaryVertex_.insert(v_jetSecondaryVertex_.end(), i_jetProperties->jetSecondaryVertex().begin(), i_jetProperties->jetSecondaryVertex().end());
            v_jetSecondaryVertexFlightDistanceValue_.insert(v_jetSecondaryVertexFlightDistanceValue_.end(), i_jetProperties->jetSecondaryVertexFlightDistanceValue().begin(), i_jetProperties->jetSecondaryVertexFlightDistanceValue().end());
            v_jetSecondaryVertexFlightDistanceSignificance_.insert(v_jetSecondaryVertexFlightDistanceSignificance_.end(), i_jetProperties->jetSecondaryVertexFlightDistanceSignificance().begin(), i_jetProperties->jetSecondaryVertexFlightDistanceSignificance().end());
            v_jetSecondaryVertexJetIndex_.insert(v_jetSecondaryVertexJetIndex_.end(), i_jetProperties->jetSecondaryVertex().size(), jetIndex);
        }
    }
    
    // PF MET and genMET
    edm::Handle<edm::View<pat::MET> > met;
    iEvent.getByLabel(metTag_, met);
    met_ = met->at(0).polarP4();
    if(isTtbarSample_) genMet_ = met->at(0).genMET()->polarP4();
    
    // MVA MET --- for CSA14, need to check whether really produced
    edm::Handle<edm::View<pat::MET> > mvaMet;
    iEvent.getByLabel(mvaMetTag_, mvaMet);
    if(!mvaMet.failedToGet()) mvaMet_ = mvaMet->at(0).polarP4();

    // Fill ntuple
    ntuple_->Fill();
}



int NTupleWriter::getTriggerBits(const edm::Event& iEvent, const edm::Handle<edm::TriggerResults>& triggerResults)
{
    int result = 0;

    const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);
    const unsigned int nTrigger = triggerResults->size();
    for(unsigned int iTrigger = 0; iTrigger < nTrigger; ++iTrigger){
        if(!(triggerResults.product()->accept(iTrigger))) continue;
        if(includeTrigger_) v_firedTriggers_.push_back(triggerNames.triggerName(iTrigger));
        const std::string& triggerName = triggerNames.triggerName(iTrigger);
        std::string triggerNameWithoutVersion(triggerName);
        while(triggerNameWithoutVersion.length() > 0
               && triggerNameWithoutVersion[triggerNameWithoutVersion.length()-1] >= '0'
               && triggerNameWithoutVersion[triggerNameWithoutVersion.length()-1] <= '9')
            triggerNameWithoutVersion.replace(triggerNameWithoutVersion.length()-1, 1, "");
        result |= triggerMap_[triggerNameWithoutVersion + "*"];
        result |= triggerMap_[triggerName];
    }
    
    return result;
}



int NTupleWriter::getTriggerBits(const std::vector<std::string>& triggerNames)
{
    int result = 0;

    for(unsigned int iTrigger = 0; iTrigger < triggerNames.size(); ++iTrigger){
        const std::string& triggerName = triggerNames.at(iTrigger);
        std::string triggerNameWithoutVersion(triggerName);
        while(triggerNameWithoutVersion.length() > 0
               && triggerNameWithoutVersion[triggerNameWithoutVersion.length()-1] >= '0'
               && triggerNameWithoutVersion[triggerNameWithoutVersion.length()-1] <= '9')
            triggerNameWithoutVersion.replace(triggerNameWithoutVersion.length()-1, 1, "");
        result |= triggerMap_[triggerNameWithoutVersion + "*"];
        result |= triggerMap_[triggerName];
    }
    
    return result;
}



int NTupleWriter::getTriggerBitsTau(const edm::Event& iEvent, const edm::Handle<edm::TriggerResults>& triggerResults)
{
    int result = 0;

    const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);
    const unsigned int nTrigger = triggerResults->size();
    for(unsigned int iTrigger = 0; iTrigger < nTrigger; ++iTrigger){
        if(!(triggerResults.product()->accept(iTrigger))) continue;
        const std::string& triggerName = triggerNames.triggerName(iTrigger);
        std::string triggerNameWithoutVersion(triggerName);
        while(triggerNameWithoutVersion.length() > 0
                && triggerNameWithoutVersion[triggerNameWithoutVersion.length()-1] >= '0'
                && triggerNameWithoutVersion[triggerNameWithoutVersion.length()-1] <= '9')
            triggerNameWithoutVersion.replace(triggerNameWithoutVersion.length()-1, 1, "");
        result |= triggerMapTau_[triggerNameWithoutVersion + "*"];
        result |= triggerMapTau_[triggerName];
    }
    
    return result;
}



void NTupleWriter::assignLeptonAndTau(const reco::GenParticle* lepton, LV& genLepton, int& pdgId, LV& genTau)
{
    const reco::GenParticle* finalLepton;
    if(this->isTau(lepton)){
        genTau = lepton->polarP4();
        finalLepton = this->tauDaughter(lepton);
    }
    else{
        genTau = nullP4_;
        finalLepton = lepton;
    }
    
    if(!this->isTau(finalLepton)){
        genLepton = finalLepton->polarP4();
        pdgId = finalLepton->pdgId();
    }
    else{
        genLepton = nullP4_;
        pdgId = 0;
    }
}



bool NTupleWriter::isTau(const reco::GenParticle* lepton)
{
    return std::abs(lepton->pdgId()) == 15;
}



const reco::GenParticle* NTupleWriter::tauDaughter(const reco::GenParticle* tau)
{
    for(size_t iDaughter = 0; iDaughter < tau->numberOfDaughters(); ++iDaughter){
        const reco::GenParticle* daughter = dynamic_cast<const reco::GenParticle*>(tau->daughter(iDaughter));
        if(std::abs(daughter->pdgId())==11 || std::abs(daughter->pdgId())==13) return daughter;
        else if(this->isTau(daughter)) return this->tauDaughter(daughter);
    }
    return tau;
}





// ------------ method called once each job just before starting event loop  ------------
void
NTupleWriter::beginJob()
{
    edm::Service<TFileService> fileService;
    if(!fileService) throw edm::Exception(edm::errors::Configuration, "TFileService is not registered in cfg file");
    ntuple_ = fileService->make<TTree>("NTuple", "NTuple");
    
    // Sample header
    //TFileDirectory header = fs->mkdir("Header");
    //header.cd();
    TObjString(sampleName_.c_str()).Write("sampleName");
    TObjString(channelName_.c_str()).Write("channelName");
    TObjString(systematicsName_.c_str()).Write("systematicsName");
    TObjString(boost::lexical_cast<std::string>(isMC_).c_str()).Write("isMC");
    TObjString(boost::lexical_cast<std::string>(isTtbarSample_).c_str()).Write("isSignal");
    TObjString(boost::lexical_cast<std::string>(isHiggsSample_).c_str()).Write("isHiggsSignal");
    
    // Event info
    ntuple_->Branch("runNumber", &runNumber_, "runNumber/i");
    ntuple_->Branch("lumiBlock", &lumiBlock_,"lumiBlock/i");
    ntuple_->Branch("eventNumber", &eventNumber_, "eventNumber/i");
    
    // Triggers
    ntuple_->Branch("triggerBits", &triggerBits_, "triggerBits/i");
    ntuple_->Branch("triggerBitsTau", &triggerBitsTau_, "triggerBitsTau/i");
    if(includeTrigger_) ntuple_->Branch("firedTriggers", &v_firedTriggers_);
    
    // Vertices
    ntuple_->Branch("vertMulti", &vertexMultiplicity_, "vertMulti/I");
    
    // Leptons
    ntuple_->Branch("leptons", &v_lepton_);
    ntuple_->Branch("lepPdgId", &v_leptonPdgId_);
    ntuple_->Branch("lepID", &v_leptonId_);
    ntuple_->Branch("lepChargedHadronIso", &v_leptonChargedHadronIso_);
    ntuple_->Branch("lepNeutralHadronIso", &v_leptonNeutralHadronIso_);
    ntuple_->Branch("lepPhotonIso", &v_leptonPhotonIso_);
    ntuple_->Branch("lepPuChargedHadronIso", &v_leptonPuChargedHadronIso_);
    ntuple_->Branch("lepPfIso", &v_leptonPfIso_);
    ntuple_->Branch("lepCombIso", &v_leptonCombIso_);
    ntuple_->Branch("lepDxyVertex0", &v_leptonDxyVertex0_);
    ntuple_->Branch("lepDzVertex0", &v_leptonDzVertex0_);
    ntuple_->Branch("lepTrigger", &v_leptonTrigger_);
    
    // Jets
    ntuple_->Branch("jets", &v_jet_);
    ntuple_->Branch("jetBTagTCHE", &v_jetBtagTCHE_);
    ntuple_->Branch("jetBTagTCHP", &v_jetBtagTCHP_);
    ntuple_->Branch("jetBTagJetProbability", &v_jetBtagJetProbability_);
    ntuple_->Branch("jetBTagJetBProbability", &v_jetBtagJetBProbability_);
    ntuple_->Branch("jetBTagSSVHE", &v_jetBtagSSVHE_);
    ntuple_->Branch("jetBTagSSVHP", &v_jetBtagSSVHP_);
    ntuple_->Branch("jetBTagCSV", &v_jetBtagCSV_);
    ntuple_->Branch("jetBTagCSVMVA", &v_jetBtagCSVMVA_);
    
    // Jet properties
    ntuple_->Branch("jetChargeGlobalPtWeighted", &v_jetChargeGlobalPtWeighted_);
    ntuple_->Branch("jetChargeRelativePtWeighted", &v_jetChargeRelativePtWeighted_);
    ntuple_->Branch("jetPfCandidateTrack", &v_jetPfCandidateTrack_);
    ntuple_->Branch("jetPfCandidateTrackCharge", &v_jetPfCandidateTrackCharge_);
    ntuple_->Branch("jetPfCandidateTrackId", &v_jetPfCandidateTrackId_);
    ntuple_->Branch("jetPfCandidatePrimaryVertexId", &v_jetPfCandidatePrimaryVertexId_);
    ntuple_->Branch("jetPfCandidateTrackIndex", &v_jetPfCandidateTrackIndex_);
    ntuple_->Branch("jetSelectedTrackMatchToPfCandidateIndex",&v_jetSelectedTrackMatchToPfCandidateIndex_);
    ntuple_->Branch("jetSelectedTrack", &v_jetSelectedTrack_);
    ntuple_->Branch("jetSelectedTrackIPValue", &v_jetSelectedTrackIPValue_);
    ntuple_->Branch("jetSelectedTrackIPSignificance", &v_jetSelectedTrackIPSignificance_);
    ntuple_->Branch("jetSelectedTrackCharge", &v_jetSelectedTrackCharge_);
    ntuple_->Branch("jetSelectedTrackIndex", &v_jetSelectedTrackIndex_);
    ntuple_->Branch("jetSecondaryVertex", &v_jetSecondaryVertex_);
    ntuple_->Branch("jetSecondaryVertexPtCorrectedMass", &v_jetSecondaryVertexPtCorrectedMass_);
    ntuple_->Branch("jetSecondaryVertexJetIndex", &v_jetSecondaryVertexJetIndex_);
    ntuple_->Branch("jetSecondaryVertexFlightDistanceValue", &v_jetSecondaryVertexFlightDistanceValue_);
    ntuple_->Branch("jetSecondaryVertexFlightDistanceSignificance", &v_jetSecondaryVertexFlightDistanceSignificance_);
    ntuple_->Branch("jetSecondaryVertexTrackMatchToSelectedTrackIndex", &v_jetSecondaryVertexTrackMatchToSelectedTrackIndex_);
    ntuple_->Branch("jetSecondaryVertexTrackVertexIndex", &v_jetSecondaryVertexTrackVertexIndex_);
    
    // Jets, but used only for MC corrections (not needed in data)
    if(isMC_){
        ntuple_->Branch("jetsForMET", &v_jetForMet_);
        ntuple_->Branch("jetJERSF", &v_jetJerSF_);
        ntuple_->Branch("jetForMETJERSF", &v_jetForMetJerSF_);
    }
    
    // MET
    ntuple_->Branch("met", &met_);
    ntuple_->Branch("mvamet", &mvaMet_);
    
    // General true-level info
    if(isMC_){
        ntuple_->Branch("vertMultiTrue", &vertexMultiplicityTrue_, "vertMultiTrue/I");
        ntuple_->Branch("associatedGenJet", &v_associatedGenJet_);
        ntuple_->Branch("associatedGenJetForMET", &v_associatedGenJetForMet_);
        ntuple_->Branch("jetPartonFlavour", &v_jetPartonFlavour_);
        ntuple_->Branch("jetPartonFlavourForMET", &v_jetPartonFlavourForMet_);
        if(includePdfWeights_) ntuple_->Branch("pdfWeights", &v_pdfWeight_);
        ntuple_->Branch("weightGenerator", &weightGenerator_, "weightGenerator/D");
    }
    
    // Ttbar generator info
    ntuple_->Branch("TopDecayMode", &ttbarDecayMode_, "TopDecayMode/I");
    if(isTtbarSample_){
        ntuple_->Branch("TopProductionMode", &ttbarProductionMode_, "TopProductionMode/I");
        ntuple_->Branch("GenTop", &genTop_);
        ntuple_->Branch("GenAntiTop", &genAntiTop_);
        ntuple_->Branch("GenLepton", &genLepton_);
        ntuple_->Branch("GenAntiLepton", &genAntiLepton_);
        ntuple_->Branch("GenLeptonPdgId", &genLeptonPdgId_, "GenLeptonPdgId/I");
        ntuple_->Branch("GenAntiLeptonPdgId", &genAntiLeptonPdgId_, "GenAntiLeptonPdgId/I");
        ntuple_->Branch("GenTau", &genTau_);
        ntuple_->Branch("GenAntiTau", &genAntiTau_);
        ntuple_->Branch("GenNeutrino", &genNeutrino_);
        ntuple_->Branch("GenAntiNeutrino", &genAntiNeutrino_);
        ntuple_->Branch("GenB", &genB_);
        ntuple_->Branch("GenAntiB", &genAntiB_);
        ntuple_->Branch("GenWPlus", &genWPlus_);
        ntuple_->Branch("GenWMinus", &genWMinus_);
        ntuple_->Branch("GenMET", &genMet_);
        ntuple_->Branch("allGenJets", &v_allGenJet_);
        //ntuple_->Branch("GenParticleP4", &v_genParticleP4_);
        //ntuple_->Branch("GenParticlePdgId", &v_genParticlePdgId_);
        //ntuple_->Branch("GenParticleStatus", &v_genParticleStatus_);
        ntuple_->Branch("BHadJetIndex", &v_bHadJetIndex_);
        ntuple_->Branch("AntiBHadJetIndex", &v_antiBHadJetIndex_);
        ntuple_->Branch("BHadrons", &v_bHadron_);
        ntuple_->Branch("AntiBHadrons", &v_antiBHadron_);
        ntuple_->Branch("BHadronFromTop", &v_bHadFromTop_);
        ntuple_->Branch("AntiBHadronFromTopB", &v_antiBHadFromTop_);
        ntuple_->Branch("BHadronVsJet", &v_bHadVsJet_);
        ntuple_->Branch("AntiBHadronVsJet", &v_antiBHadVsJet_);
        // These two branches are not in general use, could probably be deleted
        ntuple_->Branch("jetAssociatedPartonPdgId", &v_jetAssociatedPartonPdgId_);
        ntuple_->Branch("jetAssociatedParton", &v_jetAssociatedParton_);
    }
    
    // Higgs generator info
    if(isHiggsSample_){
        ntuple_->Branch("HiggsDecayMode", &higgsDecayMode_, "HiggsDecayMode/I");
        ntuple_->Branch("GenH", &genH_);
        ntuple_->Branch("GenBFromH", &genBFromH_);
        ntuple_->Branch("GenAntiBFromH", &genAntiBFromH_);
    }
    
    // Z boson generator info
    if(isZSample_){
        ntuple_->Branch("GenZ", &v_genZ_);
        ntuple_->Branch("GenZMeDaughterParticle", &v_genZMeDaughterParticle_);
        ntuple_->Branch("GenZMeDaughterAntiParticle", &v_genZMeDaughterAntiParticle_);
        ntuple_->Branch("GenZStableLepton", &v_genZStableLepton_);
        ntuple_->Branch("GenZStableAntiLepton", &v_genZStableAntiLepton_);
        ntuple_->Branch("GenZDecayMode", &v_genZDecayMode_);
    }
    
    // MadGraph W decay modes
    if(isMadgraphSample_)
        ntuple_->Branch("madgraphWDecay", &v_madgraphWDecay_);
    
    // Full parton-to-genJet matching of b/c quarks
    if(isTtbarSample_){
        ntuple_->Branch("genExtraTopJetNumberId", &genExtraTopJetNumberId_);
        ntuple_->Branch("genBHadPlusMothers", &v_genBHadPlusMothers_);
        ntuple_->Branch("genBHadPlusMothersPdgId", &v_genBHadPlusMothersPdgId_);
        if(saveHadronMothers_) {
            ntuple_->Branch("genBHadPlusMothersStatus", &v_genBHadPlusMothersStatus_);
            ntuple_->Branch("genBHadPlusMothersIndices", &v_genBHadPlusMothersIndices_);
        }
        ntuple_->Branch("genBHadIndex", &v_genBHadIndex_);
        ntuple_->Branch("genBHadFlavour", &v_genBHadFlavour_);
        ntuple_->Branch("genBHadJetIndex", &v_genBHadJetIndex_);
        ntuple_->Branch("genBHadFromTopWeakDecay", &v_genBHadFromTopWeakDecay_);
        ntuple_->Branch("genBHadLeptonHadronIndex", &v_genBHadLeptonHadronIndex_);
        ntuple_->Branch("genBHadLeptonIndex", &v_genBHadLeptonIndex_);
        ntuple_->Branch("genBHadLeptonViaTau", &v_genBHadLeptonViaTau_);
        ntuple_->Branch("genCHadPlusMothers", &v_genCHadPlusMothers_);
        ntuple_->Branch("genCHadPlusMothersPdgId", &v_genCHadPlusMothersPdgId_);
        if(saveHadronMothers_) {
            ntuple_->Branch("genCHadPlusMothersStatus", &v_genCHadPlusMothersStatus_);
        }
        if(saveHadronMothers_ || saveCHadronParticles_){
            ntuple_->Branch("genCHadIndex", &v_genCHadIndex_);
        }
        ntuple_->Branch("genCHadBHadronId", &v_genCHadBHadronId_);
        ntuple_->Branch("genCHadJetIndex", &v_genCHadJetIndex_);
        ntuple_->Branch("genCHadLeptonHadronIndex", &v_genCHadLeptonHadronIndex_);
        ntuple_->Branch("genCHadLeptonIndex", &v_genCHadLeptonIndex_);
        ntuple_->Branch("genCHadLeptonViaTau", &v_genCHadLeptonViaTau_);
    }
}



void NTupleWriter::clearVariables()
{
    // Event info
    runNumber_ = 0;
    lumiBlock_ = 0;
    eventNumber_ = 0;
    
    // Triggers
    triggerBits_ = 0;
    triggerBitsTau_ = 0;
    v_firedTriggers_.clear();
    
    // Vertices
    vertexMultiplicity_ = 0;
    
    // Leptons
    v_lepton_.clear();
    v_leptonPdgId_.clear();
    v_leptonId_.clear() ;
    v_leptonChargedHadronIso_.clear();
    v_leptonNeutralHadronIso_.clear();
    v_leptonPhotonIso_.clear();
    v_leptonPuChargedHadronIso_.clear();
    v_leptonPfIso_.clear();
    v_leptonCombIso_.clear();
    v_leptonDxyVertex0_.clear();
    v_leptonDzVertex0_.clear();
    v_leptonTrigger_.clear();
    
    // Jets
    v_jet_.clear();
    v_jetBtagTCHE_.clear();
    v_jetBtagTCHP_.clear();
    v_jetBtagJetProbability_.clear();
    v_jetBtagJetBProbability_.clear();
    v_jetBtagSSVHE_.clear();
    v_jetBtagSSVHP_.clear();
    v_jetBtagCSV_.clear();
    v_jetBtagCSVMVA_.clear();
    
    // Jet properties
    v_jetChargeGlobalPtWeighted_.clear();
    v_jetChargeRelativePtWeighted_.clear();
    v_jetPfCandidateTrack_.clear();
    v_jetPfCandidateTrackCharge_.clear();
    v_jetPfCandidateTrackId_.clear();
    v_jetPfCandidatePrimaryVertexId_.clear();
    v_jetPfCandidateTrackIndex_.clear();
    v_jetSelectedTrackMatchToPfCandidateIndex_.clear();
    v_jetSelectedTrack_.clear();
    v_jetSelectedTrackIPValue_.clear();
    v_jetSelectedTrackIPSignificance_.clear();
    v_jetSelectedTrackCharge_.clear();
    v_jetSelectedTrackIndex_.clear();
    v_jetSecondaryVertex_.clear();
    v_jetSecondaryVertexPtCorrectedMass_.clear();
    v_jetSecondaryVertexJetIndex_.clear();
    v_jetSecondaryVertexFlightDistanceValue_.clear();
    v_jetSecondaryVertexFlightDistanceSignificance_.clear();
    v_jetSecondaryVertexTrackMatchToSelectedTrackIndex_.clear();
    v_jetSecondaryVertexTrackVertexIndex_.clear();
    
    // MET
    met_ = nullP4_;
    mvaMet_ = nullP4_;
    
    // Jets, but used only for MC corrections (not needed in data)
    v_jetForMet_.clear();
    v_jetJerSF_.clear();
    v_jetForMetJerSF_.clear();
    
    // General true-level info
    vertexMultiplicityTrue_ = 0;
    v_associatedGenJet_.clear();
    v_associatedGenJetForMet_.clear();
    v_jetPartonFlavour_.clear();
    v_jetPartonFlavourForMet_.clear();
    v_pdfWeight_.clear();
    weightGenerator_ = 0.;
    
    // Ttbar generator info
    ttbarProductionMode_ = 0;
    ttbarDecayMode_ = 0;
    genTop_ = nullP4_;
    genAntiTop_ = nullP4_;
    genLepton_ = nullP4_;
    genAntiLepton_ = nullP4_;
    genLeptonPdgId_ = 0;
    genAntiLeptonPdgId_ = 0;
    genTau_ = nullP4_;
    genAntiTau_ = nullP4_;
    genNeutrino_ = nullP4_;
    genAntiNeutrino_ = nullP4_;
    genB_ = nullP4_;
    genAntiB_ = nullP4_;
    genWPlus_ = nullP4_;
    genWMinus_ = nullP4_;
    genMet_ = nullP4_;
    v_allGenJet_.clear();
    v_genParticleP4_.clear();
    v_genParticlePdgId_.clear();
    v_genParticleStatus_.clear();
    v_bHadJetIndex_.clear();
    v_antiBHadJetIndex_.clear();
    v_bHadron_.clear();
    v_antiBHadron_.clear();
    v_bHadFromTop_.clear();
    v_antiBHadFromTop_.clear();
    v_bHadVsJet_.clear();
    v_antiBHadVsJet_.clear();
    v_jetAssociatedPartonPdgId_.clear();
    v_jetAssociatedParton_.clear();
    
    // Higgs generator info
    higgsDecayMode_ = 0;
    genH_ = nullP4_;
    genBFromH_ = nullP4_;
    genAntiBFromH_ = nullP4_;
    
    // Z boson generator info
    v_genZ_.clear();
    v_genZMeDaughterParticle_.clear();
    v_genZMeDaughterAntiParticle_.clear();
    v_genZStableLepton_.clear();
    v_genZStableAntiLepton_.clear();
    v_genZDecayMode_.clear();
    
    // MadGraph W decay modes
    v_madgraphWDecay_.clear();
    
    // Full parton-to-genJet matching of b/c quarks
    genExtraTopJetNumberId_ = -1;
    v_genBHadPlusMothers_.clear();
    v_genBHadPlusMothersPdgId_.clear();
    v_genBHadPlusMothersStatus_.clear();
    v_genBHadPlusMothersIndices_.clear();
    v_genBHadIndex_.clear();
    v_genBHadFlavour_.clear();
    v_genBHadJetIndex_.clear();
    v_genBHadFromTopWeakDecay_.clear();
    v_genBHadLeptonHadronIndex_.clear();
    v_genBHadLeptonIndex_.clear();
    v_genBHadLeptonViaTau_.clear();
    v_genCHadPlusMothers_.clear();
    v_genCHadPlusMothersPdgId_.clear();
    v_genCHadPlusMothersStatus_.clear();
    v_genCHadIndex_.clear();
    v_genCHadBHadronId_.clear();
    v_genCHadJetIndex_.clear();
    v_genCHadLeptonHadronIndex_.clear();
    v_genCHadLeptonIndex_.clear();
    v_genCHadLeptonViaTau_.clear();
}


int NTupleWriter::nUniqueElementsInVector(const std::vector<int> vector, const bool countOnlyNonRepeating)
{
    // Producing a set removing duplicates
    std::set<int> set(vector.begin(), vector.end());
    if(!countOnlyNonRepeating) return set.size();
    
    int nUniqueElements = 0;
    // Looping over unique elements of the vector skipping those that appear more than once
    for(std::set<int>::iterator it = set.begin(); it!=set.end(); ++it) {
        const int nOccurances = std::count(vector.begin(), vector.end(), *it);
        if(nOccurances > 1) continue;
        nUniqueElements++;
    }
    
    return nUniqueElements;
}



// ------------ method called once each job just after ending the event loop  ------------
void
NTupleWriter::endJob()
{
}



// ------------ method called when starting to processes a run  ------------
void
NTupleWriter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}



// ------------ method called when ending the processing of a run  ------------
void
NTupleWriter::endRun(edm::Run const&, edm::EventSetup const&)
{
}



// ------------ method called when starting to processes a luminosity block  ------------
void
NTupleWriter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}



// ------------ method called when ending the processing of a luminosity block  ------------
void
NTupleWriter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NTupleWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}





//define this as a plug-in
DEFINE_FWK_MODULE(NTupleWriter);
