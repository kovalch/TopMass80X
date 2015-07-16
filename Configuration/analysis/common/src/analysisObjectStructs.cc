#include "analysisObjectStructs.h"





EventMetadata::EventMetadata()
{
    this->clear();
}



void EventMetadata::clear()
{
    runNumber_ = 0;
    lumiBlock_ = 0;
    eventNumber_ = 0;
}



RecoObjects::RecoObjects()
{
    this->clear();
}



void RecoObjects::clear()
{
    valuesSet_ = false;
    
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
    lepDzVertex0_ = 0;
    //lepTrigger_ = 0;
    jets_ = 0;
    jetBtags_ = 0;
    jetChargeGlobalPtWeighted_ = 0;
    jetChargeRelativePtWeighted_ = 0;
    jetPfCandidateTrack_ = 0;
    jetPfCandidateTrackCharge_ = 0;
    jetPfCandidateTrackId_ = 0;
    jetPfCandidateTrackIndex_ = 0;
    jetSelectedTrack_ = 0;
    jetSelectedTrackIPValue_ = 0;
    jetSelectedTrackIPSignificance_ = 0;
    jetSelectedTrackCharge_ = 0;
    jetSelectedTrackIndex_ = 0;
    jetSelectedTrackMatchToPfCandidateIndex_ = 0;
    jetSecondaryVertex_ = 0;
    jetSecondaryVertexPtCorrectedMass_ = 0;
    jetSecondaryVertexJetIndex_ = 0;
    jetSecondaryVertexFlightDistanceValue_ = 0;
    jetSecondaryVertexFlightDistanceSignificance_ = 0;
    jetSecondaryVertexTrackVertexIndex_ = 0;
    jetSecondaryVertexTrackMatchToSelectedTrackIndex_ = 0;
    jetPfCandidatePrimaryVertexId_ = 0;
    met_ = 0;
    vertMulti_ = 0;
    
    m_userInts_.clear();
    m_userDoubles_.clear();
}



CommonGenObjects::CommonGenObjects()
{
    this->clear();
}



void CommonGenObjects::clear()
{
    valuesSet_ = false;
    
    associatedGenJet_ = 0;
    jetPartonFlavour_ = 0;
}



JecObjects::JecObjects()
{
    this->clear();
}



void JecObjects::clear()
{
    valuesSet_ = false;
    
    rho_ = -999.;
    jetArea_ = 0;
    jetJERSF_ = 0;
    jetsForMET_ = 0;
    jetForMETJERSF_ = 0;
    associatedGenJetForMET_ = 0;
    //jetPartonFlavourForMET_ = 0;
}


TopGenObjects::TopGenObjects()
{
    this->clear();
}



void TopGenObjects::clear()
{
    valuesSet_ = false;
    
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
    GenMet_ = 0;
    allGenJets_ = 0;
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
    //jetAssociatedPartonPdgId_ = 0;
    //jetAssociatedParton_ = 0;
    
    genBHadPlusMothersPdgId_ = 0;
    //genBHadPlusMothersStatus_ = 0;
    //genBHadPlusMothersIndices_ = 0;
    genBHadPlusMothers_ = 0;
    genBHadIndex_ = 0;
    genBHadFlavour_ = 0;
    genBHadJetIndex_ = 0;
    genBHadLeptonIndex_ = 0;
    genBHadLeptonHadronIndex_ = 0;
    genBHadLeptonViaTau_ = 0;
    genBHadFromTopWeakDecay_ = 0;
    genCHadPlusMothersPdgId_ = 0;
    genCHadPlusMothers_ = 0;
    genCHadJetIndex_ = 0;
    genCHadLeptonIndex_ = 0;
    genCHadLeptonHadronIndex_ = 0;
    genCHadLeptonViaTau_ = 0;
    genExtraTopJetNumberId_ = -2;
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



ZGenObjects::ZGenObjects()
{
    this->clear();
}



void ZGenObjects::clear()
{
    valuesSet_ = false;
    
    GenZ_ = 0;
    GenZMeDaughterParticle_ = 0;
    GenZMeDaughterAntiParticle_ = 0;
    GenZStableLepton_ = 0;
    GenZStableAntiLepton_ = 0;
}







