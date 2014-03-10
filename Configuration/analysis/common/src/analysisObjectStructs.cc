#include "analysisObjectStructs.h"





RecoObjects::RecoObjects()
{
    this->clear();
}



void RecoObjects::clear()
{
    valuesSet_ = false;
    
    // Concerning physics objects
    allLeptons_ = 0;
    lepPdgId_ = 0;
    //lepID_ = 0;
    lepPfIso_ = 0;
    //lepChargedHadronIso_ = 0;
    //lepNeutralHadronIso_ = 0;
    //lepPhotonIso_ = 0;
    //lepPuChargedHadronIso_ = 0;
    lepCombIso_ = 0;
    lepDxyVertex0_ = 0;
    //lepDzVertex0_ = 0;
    //lepTrigger_ = 0;
    jets_ = 0;
    jetBTagTCHE_ = 0;
    //jetBTagTCHP_ = 0;
    jetBTagSSVHE_ = 0;
    //jetBTagSSVHP_ = 0;
    //jetBTagJetProbability_ = 0;
    //jetBTagJetBProbability_ = 0;
    jetBTagCSV_ = 0;
    //jetBTagCSVMVA_ = 0;
    jetChargeGlobalPtWeighted_ = 0;
    jetChargeRelativePtWeighted_ = 0;
    jetChargeRelativePtWeighted_ = 0;
    jetTrackIndex_ = 0;
    jetTrackCharge_ = 0;
    jetTrack_ = 0;
    //jetPfCandidateTrack = 0;
    //jetPfCandidateTrackCharge = 0;
    //jetPfCandidateTrackId = 0;
    //jetPfCandidateTrackIndex = 0;
    //jetSelectedTrack = 0;
    //jetSelectedTrackIPValue = 0;
    //jetSelectedTrackIPSignificance = 0;
    //jetSelectedTrackCharge = 0;
    //jetSelectedTrackIndex = 0;
    met_ = 0;
    jetJERSF_ = 0;
    jetsForMET_ = 0;
    jetForMETJERSF_ = 0;
    
    // Concerning event
    runNumber_ = 0;
    lumiBlock_ = 0;
    eventNumber_ = 0;
    //recoInChannel_ = 0;
    vertMulti_ = 0;
}



CommonGenObjects::CommonGenObjects()
{
    this->clear();
}



void CommonGenObjects::clear()
{
    valuesSet_ = false;
    
    allGenJets_ = 0;
    jetPartonFlavour_ = 0;
    associatedGenJet_ = 0;
    associatedGenJetForMET_ = 0;
    //jetAssociatedPartonPdgId_ = 0;
    //jetAssociatedParton_ = 0;
}



TopGenObjects::TopGenObjects()
{
    this->clear();
}



void TopGenObjects::clear()
{
    valuesSet_ = false;
    
    GenMet_ = 0;
    GenTop_ = 0;
    GenAntiTop_ = 0;
    GenLepton_ = 0;
    GenAntiLepton_ = 0;
    //GenLeptonPdgId_ = 0;
    //GenAntiLeptonPdgId_ = 0;
    //GenTau_ = 0;
    //GenAntiTau_ = 0;
    GenNeutrino_ = 0;
    GenAntiNeutrino_ = 0;
    GenB_ = 0;
    GenAntiB_ = 0;
    //GenWPlus_ = 0;
    //GenWMinus_ = 0;
    //GenParticleP4_= 0;
    //GenParticlePdgId_= 0;
    //GenParticleStatus_= 0;
    BHadJetIndex_ = 0;
    AntiBHadJetIndex_ = 0;
    BHadrons_ = 0;
    AntiBHadrons_ = 0;
    BHadronFromTopB_ = 0;
    AntiBHadronFromTopB_ = 0;
    BHadronVsJet_ = 0;
    AntiBHadronVsJet_ = 0;
    genBHadPlusMothersPdgId_ = 0;
    //genBHadPlusMothersStatus_ = 0;
    //genBHadPlusMothersIndices_ = 0;
    genBHadPlusMothers_ = 0;
    genBHadIndex_ = 0;
    genBHadFlavour_ = 0;
    genBHadJetIndex_ = 0;
    //genBHadLeptonIndex_ = 0;
    //genBHadLeptonHadronIndex_ = 0;
    //genBHadLeptonViaTau_ = 0;
    genBHadFromTopWeakDecay_ = 0;
}



HiggsGenObjects::HiggsGenObjects()
{
    this->clear();
}



void HiggsGenObjects::clear()
{
    valuesSet_ = false;
    
    GenH_ = 0;
    GenBFromH_ = 0;
    GenAntiBFromH_ = 0;
}



KinRecoObjects::KinRecoObjects()
{
    this->clear();
}



void KinRecoObjects::clear()
{
    valuesSet_ = false;
    
    HypTop_ = 0;
    HypAntiTop_ = 0;
    HypLepton_ = 0;
    HypAntiLepton_ = 0;
    HypNeutrino_ = 0;
    HypAntiNeutrino_ = 0;
    HypBJet_ = 0;
    HypAntiBJet_ = 0;
    //HypWPlus_ = 0;
    //HypWMinus_ = 0;
    HypJet0index_ = 0;
    HypJet1index_ = 0;
}







