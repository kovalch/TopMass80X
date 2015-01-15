#ifndef analysisObjectStructs_h
#define analysisObjectStructs_h

#include <vector>

#include <Rtypes.h>

#include "classesFwd.h"





/// Struct for holding variables associated to nTuple branches holding event meta data
struct EventMetadata{
    EventMetadata();
    ~EventMetadata(){}
    void clear();
    
    bool valuesSet_;
    
    UInt_t runNumber_;
    UInt_t lumiBlock_;
    UInt_t eventNumber_;
};


    
/// Struct for holding variables associated to nTuple branches relevant for reconstruction level
struct RecoObjects{
    RecoObjects();
    ~RecoObjects(){}
    void clear();
    
    bool valuesSet_;
    
    VLV* allLeptons_;
    std::vector<int>* lepPdgId_;
    //std::vector<double>* lepID_;
    std::vector<double>* lepPfIso_;
    //std::vector<double>* lepChargedHadronIso_;
    //std::vector<double>* lepNeutralHadronIso_;
    //std::vector<double>* lepPhotonIso_;
    //std::vector<double>* lepPuChargedHadronIso_;
    std::vector<double>* lepCombIso_;
    std::vector<double>* lepDxyVertex0_;
    std::vector<double>* lepDzVertex0_;
    //std::vector<int>* lepTrigger_;
    VLV* jets_;
    std::vector<double>* jetBTagTCHE_;
    //std::vector<double>* jetBTagTCHP_;
    std::vector<double>* jetBTagSSVHE_;
    //std::vector<double>* jetBTagSSVHP_;
    //std::vector<double>* jetBTagJetProbability_;
    //std::vector<double>* jetBTagJetBProbability_;
    std::vector<double>* jetBTagCSV_;
    //std::vector<double>* jetBTagCSVMVA_;
    std::vector<double>* jetChargeGlobalPtWeighted_;
    std::vector<double>* jetChargeRelativePtWeighted_;
    std::vector<LV>* jetPfCandidateTrack_;
    std::vector<int>* jetPfCandidateTrackCharge_;
    std::vector<int>* jetPfCandidateTrackId_;
    std::vector<int>* jetPfCandidateTrackIndex_;
    std::vector<LV>* jetSelectedTrack_;
    std::vector<double>* jetSelectedTrackIPValue_;
    std::vector<double>* jetSelectedTrackIPSignificance_;
    std::vector<int>* jetSelectedTrackCharge_;
    std::vector<int>* jetSelectedTrackIndex_;
    std::vector<int>* jetSelectedTrackMatchToPfCandidateIndex_;
    std::vector<LV>* jetSecondaryVertex_;
    std::vector<double>* jetSecondaryVertexPtCorrectedMass_;
    std::vector<int>* jetSecondaryVertexJetIndex_;
    std::vector<double>* jetSecondaryVertexFlightDistanceValue_;
    std::vector<double>* jetSecondaryVertexFlightDistanceSignificance_;
    std::vector<int>* jetSecondaryVertexTrackVertexIndex_;
    std::vector<int>* jetSecondaryVertexTrackMatchToSelectedTrackIndex_;
    std::vector<int>* eventVerticesIndices_;
    std::vector<float>* eventVerticesWeights_;
    std::vector<double>* jetPfCandidateZDistanceToVertices_;
    std::vector<int>* jetPfCandidateMatchToVerticesIndex_;
    std::vector<int>* jetPfCandidateVerticesJetIndex_;
    LV* met_;
    Int_t vertMulti_;
};



/// Struct for holding variables associated to nTuple branches holding generator information for all MC samples
struct CommonGenObjects{
    CommonGenObjects();
    ~CommonGenObjects(){}
    void clear();
    
    bool valuesSet_;
    
    std::vector<double>* jetJERSF_;
    VLV* jetsForMET_;
    std::vector<double>* jetForMETJERSF_;
    VLV* associatedGenJet_;
    VLV* associatedGenJetForMET_;
    std::vector<int>* jetPartonFlavour_;
    //std::vector<int>* jetPartonFlavourForMET_;
};



/// Struct for holding variables associated to nTuple branches for Top signal samples on generator level
struct TopGenObjects{
    TopGenObjects();
    ~TopGenObjects(){}
    void clear();
    
    bool valuesSet_;
    
    LV* GenTop_;
    LV* GenAntiTop_;
    LV* GenLepton_;
    LV* GenAntiLepton_;
    //int GenLeptonPdgId_;
    //int GenAntiLeptonPdgId_;
    //LV* GenTau_;
    //LV* GenAntiTau_;
    LV* GenNeutrino_;
    LV* GenAntiNeutrino_;
    LV* GenB_;
    LV* GenAntiB_;
    LV* GenWPlus_;
    LV* GenWMinus_;
    LV* GenMet_;
    VLV* allGenJets_;
    //std::vector<LV>* GenParticleP4_;
    //std::vector<int>* GenParticlePdgId_;
    //std::vector<int>* GenParticleStatus_;
    std::vector<int>* BHadJetIndex_;
    std::vector<int>* AntiBHadJetIndex_;
    VLV* BHadrons_;
    VLV* AntiBHadrons_;
    std::vector<bool>* BHadronFromTopB_;
    std::vector<bool>* AntiBHadronFromTopB_;
    std::vector<int>* BHadronVsJet_;
    std::vector<int>* AntiBHadronVsJet_;
    //std::vector<int>* jetAssociatedPartonPdgId_;
    //std::vector<LV>* jetAssociatedParton_;
    
    std::vector<int>* genBHadPlusMothersPdgId_;
    //std::vector<int>* genBHadPlusMothersStatus_;
    //std::vector<std::vector<int> >* genBHadPlusMothersIndices_;
    std::vector<LV>* genBHadPlusMothers_;
    std::vector<int>* genBHadIndex_;
    std::vector<int>* genBHadFlavour_;
    std::vector<int>* genBHadJetIndex_;
    std::vector<int>* genBHadLeptonIndex_;
    std::vector<int>* genBHadLeptonHadronIndex_;
    std::vector<int>* genBHadLeptonViaTau_;
    std::vector<int>* genBHadFromTopWeakDecay_;
    std::vector<int>* genCHadPlusMothersPdgId_;
    std::vector<LV>* genCHadPlusMothers_;
    std::vector<int>* genCHadJetIndex_;
    std::vector<int>* genCHadLeptonIndex_;
    std::vector<int>* genCHadLeptonHadronIndex_;
    std::vector<int>* genCHadLeptonViaTau_;
    int genExtraTopJetNumberId_;
};



/// Struct for holding variables associated to nTuple branches for Higgs signal samples on generator level
struct HiggsGenObjects{
    HiggsGenObjects();
    ~HiggsGenObjects(){}
    void clear();
    
    bool valuesSet_;
    
    LV* GenH_;
    LV* GenBFromH_;
    LV* GenAntiBFromH_;
};



/// Struct for holding variables associated to nTuple branches for Z signal samples on generator level
struct ZGenObjects{
    ZGenObjects();
    ~ZGenObjects(){}
    void clear();
    
    bool valuesSet_;
    
    VLV* GenZ_;
    VLV* GenZMeDaughterParticle_;
    VLV* GenZMeDaughterAntiParticle_;
    VLV* GenZStableLepton_;
    VLV* GenZStableAntiLepton_;
};






#endif




