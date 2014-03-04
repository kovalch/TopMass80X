#ifndef analysisObjectStructs_h
#define analysisObjectStructs_h

#include <vector>

#include <Rtypes.h>

#include "classesFwd.h"




    
/// Struct for holding variables associated to nTuple branches relevant for reconstruction level
struct RecoObjects{
    RecoObjects();
    ~RecoObjects(){}
    void clear();
    
    bool valuesSet_;
    
    // Concerning physics objects
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
    std::vector<int>* jetTrackIndex_;
    std::vector<int>* jetTrackCharge_;
    std::vector<LV>* jetTrack_;
    std::vector<int>* trackPdgId_;
    LV* met_;
    std::vector<double>* jetJERSF_;
    VLV* jetsForMET_;
    std::vector<double>* jetForMETJERSF_;
    
    // Concerning event
    UInt_t runNumber_;
    UInt_t lumiBlock_;
    UInt_t eventNumber_;
    //int recoInChannel_;
    Int_t vertMulti_;
};



/// Struct for holding variables associated to nTuple branches holding generator information for all MC samples
struct CommonGenObjects{
    CommonGenObjects();
    ~CommonGenObjects(){}
    void clear();
    
    bool valuesSet_;
    
    // Concerning physics objects
    VLV* allGenJets_;
    std::vector<int>* jetPartonFlavour_;
    VLV* associatedGenJet_;
    VLV* associatedGenJetForMET_;
    //std::vector<int>* jetAssociatedPartonPdgId_;
    //std::vector<LV>* jetAssociatedParton_;
};



/// Struct for holding variables associated to nTuple branches for Top signal samples on generator level
struct TopGenObjects{
    TopGenObjects();
    ~TopGenObjects(){}
    void clear();
    
    bool valuesSet_;
    
    LV* GenMet_;
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
    std::vector<int>* genBHadPlusMothersPdgId_;
    //std::vector<int>* genBHadPlusMothersStatus_;
    //std::vector<std::vector<int> >* genBHadPlusMothersIndices_;
    std::vector<LV>* genBHadPlusMothers_;
    std::vector<int>* genBHadIndex_;
    std::vector<int>* genBHadFlavour_;
    std::vector<int>* genBHadJetIndex_;
    //std::vector<int>* genBHadLeptonIndex_;
    //std::vector<int>* genBHadLeptonHadronIndex_;
    //std::vector<int>* genBHadLeptonViaTau_;
    std::vector<int>* genBHadFromTopWeakDecay_;
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



/// Struct for holding variables associated to nTuple branches of kinematic reconstruction
struct KinRecoObjects{
    KinRecoObjects();
    ~KinRecoObjects(){}
    void clear();
    
    bool valuesSet_;
    
    VLV* HypTop_;
    VLV* HypAntiTop_;
    VLV* HypLepton_;
    VLV* HypAntiLepton_;
    VLV* HypNeutrino_;
    VLV* HypAntiNeutrino_;
    VLV* HypBJet_;
    VLV* HypAntiBJet_;
    //VLV* HypWPlus_;
    //VLV* HypWMinus_;
    std::vector<int>* HypJet0index_;
    std::vector<int>* HypJet1index_;
};






#endif




