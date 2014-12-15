#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iomanip>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TString.h>
#include <TIterator.h>
#include <TObject.h>
#include <TObjString.h>
#include <Math/VectorUtil.h>

#include "AnalysisBase.h"
#include "analysisUtils.h"
#include "classes.h"
#include "ScaleFactors.h"
#include "Correctors.h"
#include "KinematicReconstruction.h"
#include "KinematicReconstructionSolution.h"
#include "analysisObjectStructs.h"








// ------------------------------------- Basic TSelector methods for the processing of a sample -------------------------------------


AnalysisBase::AnalysisBase(TTree*):
chain_(0),
eventMetadata_(0),
recoObjects_(0),
commonGenObjects_(0),
topGenObjects_(0),
higgsGenObjects_(0),
zGenObjects_(0),
h_weightedEvents(0),
samplename_(""),
channel_(Channel::undefined),
systematic_(),
isMC_(false),
isTopSignal_(false),
isHiggsSignal_(false),
isDrellYan_(false),
isTtbarPlusTauSample_(false),
correctMadgraphBR_(false),
channelPdgIdProduct_(0),
checkZDecayMode_(0),
outputfilename_(""),
runViaTau_(false),
isTtbarSample_(false),
eventCounter_(0),
analysisOutputBase_(0),
kinematicReconstruction_(0),
pileupScaleFactors_(0),
leptonScaleFactors_(0),
triggerScaleFactors_(0),
jetEnergyResolutionScaleFactors_(0),
jetEnergyScaleScaleFactors_(0),
btagScaleFactors_(0),
topPtScaleFactors_(0),
isSampleForBtagEfficiencies_(false),
mvaMet_(false),
metRecoilCorrector_(0)
{
    this->clearBranches();
    this->clearBranchVariables();
}



void AnalysisBase::Begin(TTree*)
{
    // WARNING! In general do not make changes here, but in your analysis' Begin function
    
    TSelector::Begin(0);
    
    eventCounter_ = 0;
}



void AnalysisBase::SlaveBegin(TTree*)
{
    // WARNING! In general do not make changes here, but in your analysis' SlaveBegin function
    
    TSelector::SlaveBegin(0);
    
    eventMetadata_ = new EventMetadata();
    recoObjects_ = new RecoObjects();
    commonGenObjects_ = new CommonGenObjects();
    topGenObjects_ = new TopGenObjects();
    higgsGenObjects_ = new HiggsGenObjects();
    zGenObjects_ = new ZGenObjects();
}



void AnalysisBase::SlaveTerminate()
{
    // WARNING! In general do not make changes here, but in your analysis' SlaveTerminate function

    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.
    
    if(eventMetadata_) delete eventMetadata_;
    if(recoObjects_) delete recoObjects_;
    if(commonGenObjects_) delete commonGenObjects_;
    if(topGenObjects_) delete topGenObjects_;
    if(higgsGenObjects_) delete higgsGenObjects_;
    if(zGenObjects_) delete zGenObjects_;
}



Bool_t AnalysisBase::Process(Long64_t)
{
    // WARNING! In general do not make changes here, but in your analysis' Process function

    if(++eventCounter_ % 100000 == 0)
        std::cout<<"Event Counter: "<<eventCounter_<<"\t--  Channel, Systematic, Sample:"
                 <<std::setw(5)<<Channel::convert(channel_)<<" ,"
                 <<std::setw(10)<<systematic_.name()<<" , "
                 <<samplename_<<std::endl;
    return kTRUE;
}



void AnalysisBase::Terminate()
{
    // WARNING! In general do not make changes here, but in your analysis' Terminate function
    
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.
    
    // Write analysis output to file
    this->writeOutput();
    
    // Cleanup
    fOutput->SetOwner();
    fOutput->Clear();
}



void AnalysisBase::writeOutput()
{
    // Do not write files in case b-tagging efficiencies are produced
    if(this->makeBtagEfficiencies()) return;
    
    // Open output file for writing
    if(!analysisOutputBase_){
        std::cerr<<"ERROR! No base directory for analysis output specified\n...break\n"<<std::endl;
        exit(3);
    }
    TString f_savename = common::assignFolder(analysisOutputBase_, channel_, systematic_);
    f_savename.Append(outputfilename_);
    std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!Finishing: "<<samplename_<<"!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    TFile outputFile(f_savename, "RECREATE");
    if (outputFile.IsZombie()) {
        std::cerr<<"ERROR in AnalysisBase::Terminate()! Cannot open file for writing: "
                 <<f_savename<<"\n...break\n"<<std::endl;
        exit(2);
    }
    
    // Write everything held by fOutput
    TIterator* iterator = fOutput->MakeIterator();
    while(TObject* obj = iterator->Next()) obj->Write();
    
    // Write additional information into file
    h_weightedEvents->Write();
    TObjString(Channel::convert(channel_)).Write("channelName");
    TObjString(systematic_.name()).Write("systematicsName");
    TObjString(samplename_).Write("sampleName");
    TObjString(isTopSignal_ ? "1" : "0").Write("isSignal");
    TObjString(isHiggsSignal_ ? "1" : "0").Write("isHiggsSignal");
    TObjString(isMC_ ? "1" : "0").Write("isMC");
    
    // Cleanup
    outputFile.Close();
    std::cout<<"Created: \033[1;1m"<<f_savename<<"\033[1;m\n\n";
}



void AnalysisBase::Init(TTree *tree)
{
    // WARNING! In general do not make changes here, but in your analysis' Init function

    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Reset all branches and their associated variables
    this->clearBranches();
    this->clearBranchVariables();

    // Set branch addresses and branch pointers
    if(!tree) return;
    chain_ = tree;
    chain_->SetMakeClass(0);
    this->SetEventMetadataBranchAddresses();
    this->SetRecoBranchAddresses();
    this->SetTriggerBranchAddresses();
    //this->SetKinRecoBranchAddresses();
    if(isMC_) this->SetCommonGenBranchAddresses();
    if(isMC_) this->SetVertMultiTrueBranchAddress();
    if(isMC_) this->SetWeightGeneratorBranchAddress();
    this->SetPdfBranchAddress();
    this->SetTopDecayBranchAddress();
    if(isTtbarPlusTauSample_ && topPtScaleFactors_) this->SetGenTopBranchAddresses();
    if(isTopSignal_) this->SetTopSignalBranchAddresses();
    if(isHiggsSignal_) this->SetHiggsDecayBranchAddress();
    if(isHiggsSignal_) this->SetHiggsSignalBranchAddresses();
    if(isDrellYan_) this->SetZDecayBranchAddress();
    if(isDrellYan_) this->SetZSignalBranchAddresses();
}





// ------------------------------------- Public methods for steering the class from outside -------------------------------------


void AnalysisBase::SetChannel(const Channel::Channel& channel)
{
    channel_ = channel;
    channelPdgIdProduct_ =
        channel == Channel::ee ? -11*11
        : channel == Channel::emu ? -11*13
        : -13*13;
}



void AnalysisBase::SetTopSignal(const bool isTopSignal)
{
    isTopSignal_ = isTopSignal;
}



void AnalysisBase::SetHiggsSignal(const bool higgsSignal)
{
    isHiggsSignal_ = higgsSignal;
}



void AnalysisBase::SetDrellYan(const bool isDrellYan)
{
    isDrellYan_ = isDrellYan;
}



void AnalysisBase::SetSystematic(const Systematic::Systematic& systematic)
{
    systematic_ = systematic;
}



void AnalysisBase::SetGeneratorBools(const TString& samplename, const Systematic::Systematic& systematic)
{
    isTtbarSample_ = samplename.BeginsWith("ttbar") && !samplename.BeginsWith("ttbarhiggs") &&
                        !(samplename=="ttbarw") && !(samplename=="ttbarz");
    isTtbarPlusTauSample_ = isTtbarSample_ && !samplename.BeginsWith("ttbarbg");
    
    const TString systematicName = systematic.name();
    correctMadgraphBR_ = samplename.BeginsWith("ttbar") && !samplename.BeginsWith("ttbarhiggs") && !systematicName.Contains("SPIN") &&
                            !systematicName.Contains("POWHEG") && !systematicName.Contains("MCATNLO");
}



void AnalysisBase::SetSamplename(const TString& samplename)
{
    samplename_ = samplename;
}



void AnalysisBase::SetMC(const bool isMC)
{
    isMC_ = isMC;
}



void AnalysisBase::SetOutputfilename(const TString& outputfilename)
{
    outputfilename_ = outputfilename;
}



void AnalysisBase::SetWeightedEvents(TH1* weightedEvents)
{
    h_weightedEvents = weightedEvents;
}



void AnalysisBase::SetRunViaTau(const bool runViaTau)
{
    runViaTau_ = runViaTau;
    if (runViaTau) isTopSignal_ = false;
}



void AnalysisBase::SetAnalysisOutputBase(const char* analysisOutputBase)
{
    analysisOutputBase_ = analysisOutputBase;
}


void AnalysisBase::SetMetRecoilCorrector(const MetRecoilCorrector* const metRecoilCorrector)
{
    metRecoilCorrector_ = metRecoilCorrector;
}


void AnalysisBase::SetKinematicReconstruction(const KinematicReconstruction* const kinematicReconstruction,
                                              const KinematicReconstructionScaleFactors* const kinematicReconstructionScaleFactors)
{
    kinematicReconstruction_ = kinematicReconstruction;
    kinematicReconstructionScaleFactors_ = kinematicReconstructionScaleFactors;
}



void AnalysisBase::SetPileupScaleFactors(const PileupScaleFactors* pileupScaleFactors)
{
    pileupScaleFactors_ = pileupScaleFactors;
}



void AnalysisBase::SetLeptonScaleFactors(const LeptonScaleFactors& scaleFactors)
{
    leptonScaleFactors_ = &scaleFactors;
}



void AnalysisBase::SetTriggerScaleFactors(const TriggerScaleFactors& scaleFactors)
{
    triggerScaleFactors_ = &scaleFactors;
}



void AnalysisBase::SetJetEnergyResolutionScaleFactors(const JetEnergyResolutionScaleFactors* jetEnergyResolutionScaleFactors)
{
    jetEnergyResolutionScaleFactors_ = jetEnergyResolutionScaleFactors;
}



void AnalysisBase::SetJetEnergyScaleScaleFactors(const JetEnergyScaleScaleFactors* jetEnergyScaleScaleFactors)
{
    jetEnergyScaleScaleFactors_ = jetEnergyScaleScaleFactors;
}



void AnalysisBase::SetBtagScaleFactors(BtagScaleFactors& scaleFactors)
{
    btagScaleFactors_ = &scaleFactors;
}



void AnalysisBase::SetSampleForBtagEfficiencies(const bool isSampleForBtagEfficiencies)
{
    isSampleForBtagEfficiencies_ = isSampleForBtagEfficiencies;
}



void AnalysisBase::SetTopPtScaleFactors(const TopPtScaleFactors* topPtScaleFactors)
{
    topPtScaleFactors_ = topPtScaleFactors;
}






// ------------------------------------- Methods for handling of ntuple branches -------------------------------------


void AnalysisBase::clearBranches()
{
    // nTuple branches holding event meta data
    b_runNumber = 0;
    b_lumiBlock = 0;
    b_eventNumber = 0;
    
    
    // nTuple branches relevant for reconstruction level
    b_lepton = 0;
    b_lepPdgId = 0;
    b_lepID = 0;
    b_lepPfIso = 0;
    b_lepChargedHadronIso = 0;
    b_lepNeutralHadronIso = 0;
    b_lepPhotonIso = 0;
    b_lepPuChargedHadronIso = 0;
    b_lepCombIso = 0;
    b_lepDxyVertex0 = 0;
    b_lepDzVertex0 = 0;
    b_lepTrigger = 0;
    b_jet = 0;
    b_jetBTagTCHE = 0;
    b_jetBTagTCHP = 0;
    b_jetBTagSSVHE = 0;
    b_jetBTagSSVHP = 0;
    b_jetBTagJetProbability = 0;
    b_jetBTagJetBProbability = 0;
    b_jetBTagCSV = 0;
    b_jetBTagCSVMVA = 0;
    b_jetChargeGlobalPtWeighted = 0;
    b_jetChargeRelativePtWeighted = 0;
    b_jetPfCandidateTrack = 0;
    b_jetPfCandidateTrackCharge = 0;
    b_jetPfCandidateTrackId = 0;
    b_jetPfCandidateTrackIndex = 0;
    b_jetSelectedTrack = 0;
    b_jetSelectedTrackIPValue = 0;
    b_jetSelectedTrackIPSignificance = 0;
    b_jetSelectedTrackCharge = 0;
    b_jetSelectedTrackIndex = 0;
    b_jetSelectedTrackMatchToPfCandidateIndex = 0;
    b_jetSecondaryVertex = 0;
    b_jetSecondaryVertexPtCorrectedMass = 0;
    b_jetSecondaryVertexJetIndex = 0;
    b_jetSecondaryVertexFlightDistanceValue = 0;
    b_jetSecondaryVertexFlightDistanceSignificance = 0;
    b_jetSecondaryVertexTrackVertexIndex = 0;
    b_jetSecondaryVertexTrackMatchToSelectedTrackIndex = 0;
    b_met = 0;
    b_vertMulti = 0;
    
    
    // nTuple branches holding trigger bits
    b_triggerBits = 0;
    b_triggerBitsTau = 0;
    b_firedTriggers = 0;
    
    
    // nTuple branches holding generator information for all MC samples
    b_jetJERSF = 0;
    b_jetForMET = 0;
    b_jetForMETJERSF = 0;
    b_associatedGenJet = 0;
    b_associatedGenJetForMET = 0;
    b_jetPartonFlavour = 0;
    b_jetPartonFlavourForMET = 0;
    
    
    // nTuple branch for true vertex multiplicity
    b_vertMultiTrue = 0;
    
    
    // nTuple branch for generator event weight
    b_weightGenerator = 0;
    
    
    // nTuple branch for PDF weights
    b_weightPDF = 0;
    
    
    // nTuple branch for Top decay mode
    b_TopDecayMode = 0;
    
    
    // nTuple branch for Higgs decay mode
    b_HiggsDecayMode = 0;
    
    
    // nTuple branch for Drell-Yan decay mode
    b_GenZDecayMode = 0;
    
    
    // nTuple branches for Top signal samples on generator level
    b_GenTop = 0;
    b_GenAntiTop = 0;
    b_GenLepton = 0;
    b_GenAntiLepton = 0;
    b_GenLeptonPdgId = 0;
    b_GenAntiLeptonPdgId = 0;
    b_GenTau = 0;
    b_GenAntiTau = 0;
    b_GenNeutrino = 0;
    b_GenAntiNeutrino = 0;
    b_GenB = 0;
    b_GenAntiB = 0;
    b_GenWPlus = 0;
    b_GenWMinus = 0;
    b_GenMet = 0;
    b_allGenJets = 0;
    b_GenParticleP4 = 0;
    b_GenParticlePdgId = 0;
    b_GenParticleStatus = 0;
    b_BHadJetIndex = 0;
    b_AntiBHadJetIndex = 0;
    b_BHadrons = 0;
    b_AntiBHadrons = 0;
    b_BHadronFromTopB = 0;
    b_AntiBHadronFromTopB = 0;
    b_BHadronVsJet = 0;
    b_AntiBHadronVsJet = 0;
    b_jetAssociatedPartonPdgId = 0;
    b_jetAssociatedParton = 0;
    b_genBHadPlusMothersPdgId = 0;
    b_genBHadPlusMothersStatus = 0;
    b_genBHadPlusMothersIndices = 0;
    b_genBHadPlusMothers = 0;
    b_genBHadIndex = 0;
    b_genBHadFlavour = 0;
    b_genBHadJetIndex = 0;
    b_genBHadLeptonIndex = 0;
    b_genBHadLeptonHadronIndex = 0;
    b_genBHadLeptonViaTau = 0;
    b_genBHadFromTopWeakDecay = 0;
    b_genCHadPlusMothersPdgId = 0;
    b_genCHadPlusMothers = 0;
    b_genCHadJetIndex = 0;
    b_genCHadLeptonIndex = 0;
    b_genCHadLeptonHadronIndex = 0;
    b_genCHadLeptonViaTau = 0;
    b_genExtraTopJetNumberId = 0;


    // nTuple branches for Higgs signal samples on generator level
    b_GenH = 0;
    b_GenBFromH = 0;
    b_GenAntiBFromH = 0;
    
    
    // nTuple branches for Z signal samples on generator level
    b_GenZ = 0;
    b_GenZMeDaughterParticle = 0;
    b_GenZMeDaughterAntiParticle = 0;
    b_GenZStableLepton = 0;
    b_GenZStableAntiLepton = 0;
}



void AnalysisBase::clearBranchVariables()
{
    // Set values to null for all variables arranged in structs
    if(eventMetadata_) eventMetadata_->clear();
    if(recoObjects_) recoObjects_->clear();
    if(commonGenObjects_) commonGenObjects_->clear();
    if(topGenObjects_) topGenObjects_->clear();
    if(higgsGenObjects_) higgsGenObjects_->clear();
    if(zGenObjects_) zGenObjects_->clear();
    
    // Set values to null for trigger bits
    triggerBits_ = 0;
    //triggerBitsTau_ = 0;
    //firedTriggers_ = 0;
    
    // Set values to null for true vertex multiplicity
    vertMultiTrue_ = 0;
    
    // Set values to null for generator event weight
    weightGenerator_ = 0;
    
    // Set values to null for PDF weight
    weightPDF_ = 0;
    
    // Set values to null for Top decay branch
    topDecayMode_ = 0;
    
    // Set values to null for Higgs decay branch
    higgsDecayMode_ = 0;
    
    // Set values to null for Z decay branch
    v_genZDecayMode_ = 0;
}



void AnalysisBase::SetEventMetadataBranchAddresses()
{
    chain_->SetBranchAddress("runNumber", &eventMetadata_->runNumber_, &b_runNumber);
    chain_->SetBranchAddress("lumiBlock", &eventMetadata_->lumiBlock_, &b_lumiBlock);
    chain_->SetBranchAddress("eventNumber", &eventMetadata_->eventNumber_, &b_eventNumber);
}



void AnalysisBase::SetRecoBranchAddresses()
{
    chain_->SetBranchAddress("leptons", &recoObjects_->allLeptons_, &b_lepton);
    chain_->SetBranchAddress("lepPdgId", &recoObjects_->lepPdgId_, &b_lepPdgId);
    //chain_->SetBranchAddress("lepID", &recoObjects_->lepID_, &b_lepID);
    chain_->SetBranchAddress("lepPfIso", &recoObjects_->lepPfIso_, &b_lepPfIso);
    //chain_->SetBranchAddress("lepChargedHadronIso", &recoObjects_->lepChargedHadronIso_, &b_lepChargedHadronIso);
    //chain_->SetBranchAddress("lepNeutralHadronIso", &recoObjects_->lepNeutralHadronIso_, &b_lepNeutralHadronIso);
    //chain_->SetBranchAddress("lepPhotonIso", &recoObjects_->lepPhotonIso_, &b_lepPhotonIso);
    //chain_->SetBranchAddress("lepPuChargedHadronIso", &recoObjects_->lepPuChargedHadronIso_, &b_lepPuChargedHadronIso);
    chain_->SetBranchAddress("lepCombIso", &recoObjects_->lepCombIso_, &b_lepCombIso);
    chain_->SetBranchAddress("lepDxyVertex0", &recoObjects_->lepDxyVertex0_, &b_lepDxyVertex0);
    if(chain_->GetBranch("lepDzVertex0")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("lepDzVertex0", &recoObjects_->lepDzVertex0_, &b_lepDzVertex0);
    //chain_->SetBranchAddress("lepTrigger", &recoObjects_->lepTrigger_, &b_lepTrigger);
    chain_->SetBranchAddress("jets", &recoObjects_->jets_, &b_jet);
    chain_->SetBranchAddress("jetBTagTCHE", &recoObjects_->jetBTagTCHE_, &b_jetBTagTCHE);
    //chain_->SetBranchAddress("jetBTagTCHP", &recoObjects_->jetBTagTCHP_, &b_jetBTagTCHP);
    chain_->SetBranchAddress("jetBTagSSVHE", &recoObjects_->jetBTagSSVHE_, &b_jetBTagSSVHE);
    //chain_->SetBranchAddress("jetBTagSSVHP", &recoObjects_->jetBTagSSVHP_, &b_jetBTagSSVHP);
    //chain_->SetBranchAddress("jetBTagJetProbability", &recoObjects_->jetBTagJetProbability_, &b_jetBTagJetProbability);
    //chain_->SetBranchAddress("jetBTagJetBProbability", &recoObjects_->jetBTagJetBProbability_, &b_jetBTagJetBProbability);
    chain_->SetBranchAddress("jetBTagCSV", &recoObjects_->jetBTagCSV_, &b_jetBTagCSV);
    //chain_->SetBranchAddress("jetBTagCSVMVA", &recoObjects_->jetBTagCSVMVA_, &b_jetBTagCSVMVA);
    if(chain_->GetBranch("jetChargeGlobalPtWeighted")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetChargeGlobalPtWeighted", &recoObjects_->jetChargeGlobalPtWeighted_, &b_jetChargeGlobalPtWeighted);
    if(chain_->GetBranch("jetChargeRelativePtWeighted")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetChargeRelativePtWeighted", &recoObjects_->jetChargeRelativePtWeighted_, &b_jetChargeRelativePtWeighted);
    if(chain_->GetBranch("jetPfCandidateTrack")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetPfCandidateTrack", &recoObjects_->jetPfCandidateTrack_, &b_jetPfCandidateTrack);
    if(chain_->GetBranch("jetPfCandidateTrackCharge")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetPfCandidateTrackCharge", &recoObjects_->jetPfCandidateTrackCharge_, &b_jetPfCandidateTrackCharge);
    if(chain_->GetBranch("jetPfCandidateTrackId")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetPfCandidateTrackId", &recoObjects_->jetPfCandidateTrackId_, &b_jetPfCandidateTrackId);
    if(chain_->GetBranch("jetPfCandidateTrackIndex")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetPfCandidateTrackIndex", &recoObjects_->jetPfCandidateTrackIndex_, &b_jetPfCandidateTrackIndex);
    if(chain_->GetBranch("jetSelectedTrack")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetSelectedTrack", &recoObjects_->jetSelectedTrack_, &b_jetSelectedTrack);
    if(chain_->GetBranch("jetSelectedTrackIPValue")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetSelectedTrackIPValue", &recoObjects_->jetSelectedTrackIPValue_, &b_jetSelectedTrackIPValue);
    if(chain_->GetBranch("jetSelectedTrackIPSignificance")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetSelectedTrackIPSignificance", &recoObjects_->jetSelectedTrackIPSignificance_, &b_jetSelectedTrackIPSignificance);
    if(chain_->GetBranch("jetSelectedTrackCharge")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetSelectedTrackCharge", &recoObjects_->jetSelectedTrackCharge_, &b_jetSelectedTrackCharge);
    if(chain_->GetBranch("jetSelectedTrackIndex")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetSelectedTrackIndex", &recoObjects_->jetSelectedTrackIndex_, &b_jetSelectedTrackIndex);
    if(chain_->GetBranch("jetSelectedTrackMatchToPfCandidateIndex")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetSelectedTrackMatchToPfCandidateIndex", &recoObjects_->jetSelectedTrackMatchToPfCandidateIndex_, &b_jetSelectedTrackMatchToPfCandidateIndex);
    if(chain_->GetBranch("jetSecondaryVertex")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetSecondaryVertex", &recoObjects_->jetSecondaryVertex_, &b_jetSecondaryVertex);
    if(chain_->GetBranch("jetSecondaryVertexPtCorrectedMass")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetSecondaryVertexPtCorrectedMass", &recoObjects_->jetSecondaryVertexPtCorrectedMass_, &b_jetSecondaryVertexPtCorrectedMass);
    if(chain_->GetBranch("jetSecondaryVertexJetIndex")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetSecondaryVertexJetIndex", &recoObjects_->jetSecondaryVertexJetIndex_, &b_jetSecondaryVertexJetIndex);
    if(chain_->GetBranch("jetSecondaryVertexFlightDistanceValue")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetSecondaryVertexFlightDistanceValue", &recoObjects_->jetSecondaryVertexFlightDistanceValue_, &b_jetSecondaryVertexFlightDistanceValue);
    if(chain_->GetBranch("jetSecondaryVertexFlightDistanceSignificance")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetSecondaryVertexFlightDistanceSignificance", &recoObjects_->jetSecondaryVertexFlightDistanceSignificance_, &b_jetSecondaryVertexFlightDistanceSignificance);
    if(chain_->GetBranch("jetSecondaryVertexTrackVertexIndex")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetSecondaryVertexTrackVertexIndex", &recoObjects_->jetSecondaryVertexTrackVertexIndex_, &b_jetSecondaryVertexTrackVertexIndex);
    if(chain_->GetBranch("jetSecondaryVertexTrackMatchToSelectedTrackIndex")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("jetSecondaryVertexTrackMatchToSelectedTrackIndex", &recoObjects_->jetSecondaryVertexTrackMatchToSelectedTrackIndex_, &b_jetSecondaryVertexTrackMatchToSelectedTrackIndex);
    if(mvaMet_){
        if(chain_->GetBranch("mvamet")) // new variable, keep check a while for compatibility
            chain_->SetBranchAddress("mvamet", &recoObjects_->met_, &b_met);
        else{
            std::cerr<<"ERROR in AnalysisBase::SetRecoBranchAddresses()! MVA MET is requested, but not found in nTuple\n...break\n"<<std::endl;
            exit(237);
        }
    }
    else{
        chain_->SetBranchAddress("met", &recoObjects_->met_, &b_met);
    }
    chain_->SetBranchAddress("vertMulti", &recoObjects_->vertMulti_, &b_vertMulti);
}



void AnalysisBase::SetTriggerBranchAddresses()
{
    chain_->SetBranchAddress("triggerBits", &triggerBits_, &b_triggerBits);
    //chain_->SetBranchAddress("triggerBitsTau", &triggerBitsTau_, &b_triggerBitsTau);
    //chain_->SetBranchAddress("firedTriggers", &firedTriggers_, &b_firedTriggers);
}



void AnalysisBase::SetCommonGenBranchAddresses()
{
    if(jetEnergyResolutionScaleFactors_ || jetEnergyScaleScaleFactors_){
        chain_->SetBranchAddress("jetsForMET", &commonGenObjects_->jetsForMET_, &b_jetForMET);
    }
    if(jetEnergyResolutionScaleFactors_){
        chain_->SetBranchAddress("jetJERSF", &commonGenObjects_->jetJERSF_, &b_jetJERSF);
        chain_->SetBranchAddress("jetForMETJERSF", &commonGenObjects_->jetForMETJERSF_, &b_jetForMETJERSF);
    }
    chain_->SetBranchAddress("associatedGenJet", &commonGenObjects_->associatedGenJet_, &b_associatedGenJet);
    chain_->SetBranchAddress("jetPartonFlavour", &commonGenObjects_->jetPartonFlavour_, &b_jetPartonFlavour);
    if(jetEnergyResolutionScaleFactors_){
        chain_->SetBranchAddress("associatedGenJetForMET", &commonGenObjects_->associatedGenJetForMET_, &b_associatedGenJetForMET);
        //chain_->SetBranchAddress("jetPartonFlavourForMET", &commonGenObjects_->jetPartonFlavourForMET_, &b_jetPartonFlavourForMET);
    }
}



void AnalysisBase::SetVertMultiTrueBranchAddress()
{
    chain_->SetBranchAddress("vertMultiTrue", &vertMultiTrue_, &b_vertMultiTrue);
}



void AnalysisBase::SetWeightGeneratorBranchAddress()
{
    chain_->SetBranchAddress("weightGenerator", &weightGenerator_, &b_weightGenerator);
}



void AnalysisBase::SetPdfBranchAddress()
{
    if(chain_->GetBranch("pdfWeights")) chain_->SetBranchAddress("pdfWeights", &weightPDF_, &b_weightPDF);
}



void AnalysisBase::SetTopDecayBranchAddress()
{
    chain_->SetBranchAddress("TopDecayMode", &topDecayMode_, &b_TopDecayMode);
}



void AnalysisBase::SetHiggsDecayBranchAddress()
{
    chain_->SetBranchAddress("HiggsDecayMode", &higgsDecayMode_, &b_HiggsDecayMode);
}



void AnalysisBase::SetZDecayBranchAddress()
{
    chain_->SetBranchAddress("GenZDecayMode", &v_genZDecayMode_, &b_GenZDecayMode);
}



void AnalysisBase::SetGenTopBranchAddresses()
{
    chain_->SetBranchAddress("GenTop", &topGenObjects_->GenTop_, &b_GenTop);
    chain_->SetBranchAddress("GenAntiTop", &topGenObjects_->GenAntiTop_, &b_GenAntiTop);
}



void AnalysisBase::SetTopSignalBranchAddresses()
{
    this->SetGenTopBranchAddresses();
    chain_->SetBranchAddress("GenLepton", &topGenObjects_->GenLepton_, &b_GenLepton);
    chain_->SetBranchAddress("GenAntiLepton", &topGenObjects_->GenAntiLepton_, &b_GenAntiLepton);
    //chain_->SetBranchAddress("GenLeptonPdgId", &topGenObjects_->GenLeptonPdgId_, &b_GenLeptonPdgId);
    //chain_->SetBranchAddress("GenAntiLeptonPdgId", &topGenObjects_->GenAntiLeptonPdgId_, &b_GenAntiLeptonPdgId);
    //chain_->SetBranchAddress("GenTau", &topGenObjects_->GenTau_, &b_GenTau);
    //chain_->SetBranchAddress("GenAntiTau", &topGenObjects_->GenAntiTau_, &b_GenAntiTau);
    chain_->SetBranchAddress("GenNeutrino", &topGenObjects_->GenNeutrino_, &b_GenNeutrino);
    chain_->SetBranchAddress("GenAntiNeutrino", &topGenObjects_->GenAntiNeutrino_, &b_GenAntiNeutrino);
    chain_->SetBranchAddress("GenB", &topGenObjects_->GenB_, &b_GenB);
    chain_->SetBranchAddress("GenAntiB", &topGenObjects_->GenAntiB_, &b_GenAntiB);
    //chain_->SetBranchAddress("GenWPlus", &topGenObjects_->GenWPlus_, &b_GenWPlus);
    //chain_->SetBranchAddress("GenWMinus", &topGenObjects_->GenWMinus_, &b_GenWMinus);
    //chain_->SetBranchAddress("GenWPlus.fCoordinates.fX", &topGenObjects_->GenWPluspX_, &b_GenWPluspX);
    //chain_->SetBranchAddress("GenWMinus.fCoordinates.fX", &topGenObjects_->GenWMinuspX_, &b_GenWMinuspX);
    chain_->SetBranchAddress("GenMET", &topGenObjects_->GenMet_, &b_GenMet);
    chain_->SetBranchAddress("allGenJets", &topGenObjects_->allGenJets_, &b_allGenJets);
    //chain_->SetBranchAddress("GenParticleP4", &topGenObjects_->GenParticleP4_, &b_GenParticleP4);
    //chain_->SetBranchAddress("GenParticlePdgId", &topGenObjects_->GenParticlePdgId_, &b_GenParticlePdgId);
    //chain_->SetBranchAddress("GenParticleStatus", &topGenObjects_->GenParticleStatus_, &b_GenParticleStatus);
    chain_->SetBranchAddress("BHadJetIndex", &topGenObjects_->BHadJetIndex_, &b_BHadJetIndex);
    chain_->SetBranchAddress("AntiBHadJetIndex", &topGenObjects_->AntiBHadJetIndex_, &b_AntiBHadJetIndex);
    chain_->SetBranchAddress("BHadrons", &topGenObjects_->BHadrons_, &b_BHadrons);
    chain_->SetBranchAddress("AntiBHadrons", &topGenObjects_->AntiBHadrons_, &b_AntiBHadrons);
    chain_->SetBranchAddress("BHadronFromTop", &topGenObjects_->BHadronFromTopB_, &b_BHadronFromTopB);
    chain_->SetBranchAddress("AntiBHadronFromTopB", &topGenObjects_->AntiBHadronFromTopB_, &b_AntiBHadronFromTopB);
    chain_->SetBranchAddress("BHadronVsJet", &topGenObjects_->BHadronVsJet_, &b_BHadronVsJet);
    chain_->SetBranchAddress("AntiBHadronVsJet", &topGenObjects_->AntiBHadronVsJet_, &b_AntiBHadronVsJet);
    //chain_->SetBranchAddress("jetAssociatedPartonPdgId", &topGenObjects_->jetAssociatedPartonPdgId_, &b_jetAssociatedPartonPdgId);
    //chain_->SetBranchAddress("jetAssociatedParton", &topGenObjects_->jetAssociatedParton_, &b_jetAssociatedParton);
    
    if(chain_->GetBranch("genBHadPlusMothersPdgId")) // need to check whether branch exists
        chain_->SetBranchAddress("genBHadPlusMothersPdgId", &topGenObjects_->genBHadPlusMothersPdgId_, &b_genBHadPlusMothersPdgId);
    //if(chain_->GetBranch("genBHadPlusMothersStatus")) // need to check whether branch exists
    //    chain_->SetBranchAddress("genBHadPlusMothersStatus", &topGenObjects_->genBHadPlusMothersStatus_, &b_genBHadPlusMothersStatus);
    //if(chain_->GetBranch("genBHadPlusMothersIndices")) // need to check whether branch exists
    //    chain_->SetBranchAddress("genBHadPlusMothersIndices", &topGenObjects_->genBHadPlusMothersIndices_, &b_genBHadPlusMothersIndices);
    if(chain_->GetBranch("genBHadPlusMothers")) // need to check whether branch exists
        chain_->SetBranchAddress("genBHadPlusMothers", &topGenObjects_->genBHadPlusMothers_, &b_genBHadPlusMothers);
    if(chain_->GetBranch("genBHadIndex")) // need to check whether branch exists
        chain_->SetBranchAddress("genBHadIndex", &topGenObjects_->genBHadIndex_, &b_genBHadIndex);
    if(chain_->GetBranch("genBHadFlavour")) // need to check whether branch exists
        chain_->SetBranchAddress("genBHadFlavour", &topGenObjects_->genBHadFlavour_, &b_genBHadFlavour);
    if(chain_->GetBranch("genBHadJetIndex")) // need to check whether branch exists
        chain_->SetBranchAddress("genBHadJetIndex", &topGenObjects_->genBHadJetIndex_, &b_genBHadJetIndex);
    if(chain_->GetBranch("genBHadLeptonIndex")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("genBHadLeptonIndex", &topGenObjects_->genBHadLeptonIndex_, &b_genBHadLeptonIndex);
    if(chain_->GetBranch("genBHadLeptonHadronIndex")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("genBHadLeptonHadronIndex", &topGenObjects_->genBHadLeptonHadronIndex_, &b_genBHadLeptonHadronIndex);
    if(chain_->GetBranch("genBHadLeptonViaTau")) // new variable, keep check a while for compatibility
        chain_->SetBranchAddress("genBHadLeptonViaTau", &topGenObjects_->genBHadLeptonViaTau_, &b_genBHadLeptonViaTau);
    if(chain_->GetBranch("genBHadFromTopWeakDecay")) // need to check whether branch exists
        chain_->SetBranchAddress("genBHadFromTopWeakDecay", &topGenObjects_->genBHadFromTopWeakDecay_, &b_genBHadFromTopWeakDecay);
    
    //if(chain_->GetBranch("genCHadPlusMothersPdgId")) // need to check whether branch exists
    //   chain_->SetBranchAddress("genCHadPlusMothersPdgId", &topGenObjects_->genCHadPlusMothersPdgId_, &b_genCHadPlusMothersPdgId);
    //if(chain_->GetBranch("genCHadPlusMothers")) // need to check whether branch exists
    //    chain_->SetBranchAddress("genCHadPlusMothers", &topGenObjects_->genCHadPlusMothers_, &b_genCHadPlusMothers);
    if(chain_->GetBranch("genCHadJetIndex")) // need to check whether branch exists
        chain_->SetBranchAddress("genCHadJetIndex", &topGenObjects_->genCHadJetIndex_, &b_genCHadJetIndex);
    //if(chain_->GetBranch("genCHadLeptonIndex")) // new variable, keep check a while for compatibility
    //    chain_->SetBranchAddress("genCHadLeptonIndex", &topGenObjects_->genCHadLeptonIndex_, &b_genCHadLeptonIndex);
    //if(chain_->GetBranch("genCHadLeptonHadronIndex")) // new variable, keep check a while for compatibility
    //    chain_->SetBranchAddress("genCHadLeptonHadronIndex", &topGenObjects_->genCHadLeptonHadronIndex_, &b_genCHadLeptonHadronIndex);
    //if(chain_->GetBranch("genCHadLeptonViaTau")) // new variable, keep check a while for compatibility
    //    chain_->SetBranchAddress("genCHadLeptonViaTau", &topGenObjects_->genCHadLeptonViaTau_, &b_genCHadLeptonViaTau);
    
    if(chain_->GetBranch("genExtraTopJetNumberId")) // need to check whether branch exists
        chain_->SetBranchAddress("genExtraTopJetNumberId", &topGenObjects_->genExtraTopJetNumberId_, &b_genExtraTopJetNumberId);
}



void AnalysisBase::SetHiggsSignalBranchAddresses()
{
    chain_->SetBranchAddress("GenH", &higgsGenObjects_->GenH_, &b_GenH);
    chain_->SetBranchAddress("GenBFromH", &higgsGenObjects_->GenBFromH_, &b_GenBFromH);
    chain_->SetBranchAddress("GenAntiBFromH", &higgsGenObjects_->GenAntiBFromH_, &b_GenAntiBFromH);
}



void AnalysisBase::SetZSignalBranchAddresses()
{
    chain_->SetBranchAddress("GenZ", &zGenObjects_->GenZ_, &b_GenZ);
    chain_->SetBranchAddress("GenZMeDaughterParticle", &zGenObjects_->GenZMeDaughterParticle_, &b_GenZMeDaughterParticle);
    chain_->SetBranchAddress("GenZMeDaughterAntiParticle", &zGenObjects_->GenZMeDaughterAntiParticle_, &b_GenZMeDaughterAntiParticle);
    chain_->SetBranchAddress("GenZStableLepton", &zGenObjects_->GenZStableLepton_, &b_GenZStableLepton);
    chain_->SetBranchAddress("GenZStableAntiLepton", &zGenObjects_->GenZStableAntiLepton_, &b_GenZStableAntiLepton);
}



void AnalysisBase::GetEventMetadataBranchesEntry(const Long64_t& entry)const
{
    // Check if branches' entry is already read
    if(eventMetadata_->valuesSet_) return;
    eventMetadata_->valuesSet_ = true;

    b_runNumber->GetEntry(entry);
    b_lumiBlock->GetEntry(entry);
    b_eventNumber->GetEntry(entry);
}



void AnalysisBase::GetRecoBranchesEntry(const Long64_t& entry)const
{
    // Check if branches' entry is already read
    if(recoObjects_->valuesSet_) return;
    recoObjects_->valuesSet_ = true;

    // Concerning physics objects
    b_lepton->GetEntry(entry);
    b_lepPdgId->GetEntry(entry);
    //b_lepID->GetEntry(entry);
    //b_lepPfIso->GetEntry(entry);
    //b_lepChargedHadronIso->GetEntry(entry);
    //b_lepNeutralHadronIso->GetEntry(entry);
    //b_lepPhotonIso->GetEntry(entry);
    //b_lepPuChargedHadronIso->GetEntry(entry);
    //b_lepCombIso->GetEntry(entry);
    b_lepDxyVertex0->GetEntry(entry);
    if(b_lepDzVertex0) b_lepDzVertex0->GetEntry(entry);
    //b_lepTrigger->GetEntry(entry);
    b_jet->GetEntry(entry);
    //b_jetBTagTCHE->GetEntry(entry);
    //b_jetBTagTCHP->GetEntry(entry);
    //b_jetBTagSSVHE->GetEntry(entry);
    //b_jetBTagSSVHP->GetEntry(entry);
    //b_jetBTagJetProbability->GetEntry(entry);
    //b_jetBTagJetBProbability->GetEntry(entry);
    b_jetBTagCSV->GetEntry(entry);
    //b_jetBTagCSVMVA->GetEntry(entry);
    if(b_jetChargeGlobalPtWeighted) b_jetChargeGlobalPtWeighted->GetEntry(entry);
    if(b_jetChargeRelativePtWeighted) b_jetChargeRelativePtWeighted->GetEntry(entry);
//     if(b_jetPfCandidateTrack) b_jetPfCandidateTrack->GetEntry(entry);
//     if(b_jetPfCandidateTrackCharge) b_jetPfCandidateTrackCharge->GetEntry(entry);
//     if(b_jetPfCandidateTrackId) b_jetPfCandidateTrackId->GetEntry(entry);
//     if(b_jetPfCandidateTrackIndex) b_jetPfCandidateTrackIndex->GetEntry(entry);
//     if(b_jetSelectedTrack) b_jetSelectedTrack->GetEntry(entry);
//     if(b_jetSelectedTrackMatchToPfCandidateIndex) b_jetSelectedTrackMatchToPfCandidateIndex->GetEntry(entry);
//     if(b_jetSelectedTrackIPValue) b_jetSelectedTrackIPValue->GetEntry(entry);
//     if(b_jetSelectedTrackIPSignificance) b_jetSelectedTrackIPSignificance->GetEntry(entry);
//     if(b_jetSelectedTrackCharge) b_jetSelectedTrackCharge->GetEntry(entry);
//     if(b_jetSelectedTrackIndex) b_jetSelectedTrackIndex->GetEntry(entry);
//     if(b_jetSecondaryVertex) b_jetSecondaryVertex->GetEntry(entry);
//     if(b_jetSecondaryVertexPtCorrectedMass) b_jetSecondaryVertexPtCorrectedMass->GetEntry(entry);
//     if(b_jetSecondaryVertexJetIndex) b_jetSecondaryVertexJetIndex->GetEntry(entry);
//     if(b_jetSecondaryVertexFlightDistanceValue) b_jetSecondaryVertexFlightDistanceValue->GetEntry(entry);
//     if(b_jetSecondaryVertexFlightDistanceSignificance) b_jetSecondaryVertexFlightDistanceSignificance->GetEntry(entry);
//     if(b_jetSecondaryVertexTrackVertexIndex) b_jetSecondaryVertexTrackVertexIndex->GetEntry(entry);
//     if(b_jetSecondaryVertexTrackMatchToSelectedTrackIndex) b_jetSecondaryVertexTrackMatchToSelectedTrackIndex->GetEntry(entry);
    b_met->GetEntry(entry);
    b_vertMulti->GetEntry(entry);

    // Concerning event
    b_runNumber->GetEntry(entry);
    b_lumiBlock->GetEntry(entry);
    b_eventNumber->GetEntry(entry);
}



void AnalysisBase::GetTriggerBranchesEntry(const Long64_t& entry)const
{
    b_triggerBits->GetEntry(entry);
    //b_triggerBitsTau->GetEntry(entry);
    //b_firedTriggers->GetEntry(entry);
}



void AnalysisBase::GetCommonGenBranchesEntry(const Long64_t& entry)const
{
    // Check if branches' entry is already read
    if(commonGenObjects_->valuesSet_) return;
    commonGenObjects_->valuesSet_ = true;

    // Concerning physics objects
    if(b_jetForMET) b_jetForMET->GetEntry(entry);
    if(b_jetJERSF) b_jetJERSF->GetEntry(entry);
    if(b_jetForMETJERSF) b_jetForMETJERSF->GetEntry(entry);
    b_associatedGenJet->GetEntry(entry);
    b_jetPartonFlavour->GetEntry(entry);
    if(b_associatedGenJetForMET) b_associatedGenJetForMET->GetEntry(entry);
    //if(b_jetPartonFlavourForMET) b_jetPartonFlavourForMET->GetEntry(entry);
}



void AnalysisBase::GetVertMultiTrueEntry(const Long64_t& entry)const
{
    b_vertMultiTrue->GetEntry(entry);
}



void AnalysisBase::GetWeightGeneratorEntry(const Long64_t& entry)const
{
    b_weightGenerator->GetEntry(entry);
}



void AnalysisBase::GetPDFEntry(const Long64_t& entry)const
{
    b_weightPDF->GetEntry(entry);
}



void AnalysisBase::SetTrueLevelDYChannel(const int dy)
{
    if(dy){
        std::cout<<"Include true-level filter for Z decay to PDG ID: "<<dy<<"\n";
        
        // Create function to check the Z decay channel
        checkZDecayMode_ = [&, dy](Long64_t entry) -> bool {
            this->GetZDecayModeEntry(entry);
            bool pass = false;
            // Loop over all Zs
            for(const auto decayMode : *v_genZDecayMode_){
                if((dy == 15 && decayMode >= 150000) ||
                   (dy == 13 && decayMode == 13) ||
                   (dy == 11 && decayMode == 11)){
                    pass = true;
                    break;
                }
            }
            return pass;
        };

    }
    else{
        checkZDecayMode_ = nullptr;
    }
}



void AnalysisBase::GetTopDecayModeEntry(const Long64_t& entry)const
{
    b_TopDecayMode->GetEntry(entry);
}



void AnalysisBase::GetHiggsDecayModeEntry(const Long64_t& entry)const
{
    b_HiggsDecayMode->GetEntry(entry);
}



void AnalysisBase::GetZDecayModeEntry(const Long64_t& entry)const
{
    b_GenZDecayMode->GetEntry(entry);
}


void AnalysisBase::GetGenTopBranchesEntry(const Long64_t& entry)const
{
    b_GenTop->GetEntry(entry);
    b_GenAntiTop->GetEntry(entry);
}



void AnalysisBase::GetTopSignalBranchesEntry(const Long64_t& entry)const
{
    // Check if branches' entry is already read
    if(topGenObjects_->valuesSet_) return;
    topGenObjects_->valuesSet_ = true;

    this->GetGenTopBranchesEntry(entry);
    b_GenLepton->GetEntry(entry);
    b_GenAntiLepton->GetEntry(entry);
    //b_GenLeptonPdgId->GetEntry(entry);
    //b_GenAntiLeptonPdgId->GetEntry(entry);
    //b_GenTau->GetEntry(entry);
    //b_GenAntiTau->GetEntry(entry);
    b_GenNeutrino->GetEntry(entry);
    b_GenAntiNeutrino->GetEntry(entry);
    b_GenB->GetEntry(entry);
    b_GenAntiB->GetEntry(entry);
    //b_GenWPlus->GetEntry(entry);
    //b_GenWMinus->GetEntry(entry);
    b_GenMet->GetEntry(entry);
    b_allGenJets->GetEntry(entry);
    //b_GenParticleP4->GetEntry(entry);
    //b_GenParticlePdgId->GetEntry(entry);
    //b_GenParticleStatus->GetEntry(entry);
    b_BHadJetIndex->GetEntry(entry);
    b_AntiBHadJetIndex->GetEntry(entry);
    b_BHadrons->GetEntry(entry);
    b_AntiBHadrons->GetEntry(entry);
    b_BHadronFromTopB->GetEntry(entry);
    b_AntiBHadronFromTopB->GetEntry(entry);
    b_BHadronVsJet->GetEntry(entry);
    b_AntiBHadronVsJet->GetEntry(entry);
    //if(b_jetAssociatedPartonPdgId) b_jetAssociatedPartonPdgId->GetEntry(entry);
    //if(b_jetAssociatedParton) b_jetAssociatedParton->GetEntry(entry);
    if(b_genBHadPlusMothersPdgId) b_genBHadPlusMothersPdgId->GetEntry(entry);
    //if(b_genBHadPlusMothersStatus) b_genBHadPlusMothersStatus->GetEntry(entry);
    //if(b_genBHadPlusMothersIndices) b_genBHadPlusMothersIndices->GetEntry(entry);
    if(b_genBHadPlusMothers) b_genBHadPlusMothers->GetEntry(entry);
    if(b_genBHadIndex) b_genBHadIndex->GetEntry(entry);
    if(b_genBHadFlavour) b_genBHadFlavour->GetEntry(entry);
    if(b_genBHadJetIndex) b_genBHadJetIndex->GetEntry(entry);
    if(b_genBHadLeptonIndex) b_genBHadLeptonIndex->GetEntry(entry);
    if(b_genBHadLeptonHadronIndex) b_genBHadLeptonHadronIndex->GetEntry(entry);
    if(b_genBHadLeptonViaTau) b_genBHadLeptonViaTau->GetEntry(entry);
    if(b_genBHadFromTopWeakDecay) b_genBHadFromTopWeakDecay->GetEntry(entry);
    //if(b_genCHadPlusMothersPdgId) b_genCHadPlusMothersPdgId->GetEntry(entry);
    //if(b_genCHadPlusMothers) b_genCHadPlusMothers->GetEntry(entry);
    if(b_genCHadJetIndex) b_genCHadJetIndex->GetEntry(entry);
    //if(b_genCHadLeptonIndex) b_genCHadLeptonIndex->GetEntry(entry);
    //if(b_genCHadLeptonHadronIndex) b_genCHadLeptonHadronIndex->GetEntry(entry);
    //if(b_genCHadLeptonViaTau) b_genCHadLeptonViaTau->GetEntry(entry);
    if(b_genExtraTopJetNumberId) b_genExtraTopJetNumberId->GetEntry(entry);
}



void AnalysisBase::GetHiggsSignalBranchesEntry(const Long64_t& entry)const
{
    // Check if branches' entry is already read
    if(higgsGenObjects_->valuesSet_) return;
    higgsGenObjects_->valuesSet_ = true;

    b_GenH->GetEntry(entry);
    b_GenBFromH->GetEntry(entry);
    b_GenAntiBFromH->GetEntry(entry);
}



void AnalysisBase::GetZSignalBranchesEntry(const Long64_t& entry)const
{
    // Check if branches' entry is already read
    if(zGenObjects_->valuesSet_) return;
    zGenObjects_->valuesSet_ = true;
    
    b_GenZ->GetEntry(entry);
    b_GenZMeDaughterParticle->GetEntry(entry);
    b_GenZMeDaughterAntiParticle->GetEntry(entry);
    b_GenZStableLepton->GetEntry(entry);
    b_GenZStableAntiLepton->GetEntry(entry);
}





// ------------------------------------- Methods for accessing the object structs of ntuple branches -------------------------------------


const EventMetadata& AnalysisBase::getEventMetadata(const Long64_t& entry)const
{
    if(eventMetadata_->valuesSet_) return *eventMetadata_;
    
    this->GetEventMetadataBranchesEntry(entry);
    return *eventMetadata_;
}



const RecoObjects& AnalysisBase::getRecoObjects(const Long64_t& entry)const
{
    if(recoObjects_->valuesSet_) return *recoObjects_;
    
    this->GetRecoBranchesEntry(entry);
    
    // Apply systematic variations in the following, of course only in MC
    if(!isMC_) return *recoObjects_;
    
    if(jetEnergyResolutionScaleFactors_){
        // Get references for all relevant objects for jet variation
        VLV* v_jet = recoObjects_->jets_;
        if(!commonGenObjects_->valuesSet_) this->GetCommonGenBranchesEntry(entry);
        const std::vector<double>* v_jetJerSF = commonGenObjects_->jetJERSF_;
        const VLV* v_associatedGenJet = commonGenObjects_->associatedGenJet_;
        // Apply systematic variation to jets
        jetEnergyResolutionScaleFactors_->applyJetSystematic(v_jet, v_jetJerSF, v_associatedGenJet);
        
        if(!mvaMet_){
            // Get references for all relevant objects for MET variation
            LV* met = recoObjects_->met_;
            VLV* v_jetForMet = commonGenObjects_->jetsForMET_;
            const std::vector<double>* v_jetForMetJerSF = commonGenObjects_->jetForMETJERSF_;
            const VLV* v_associatedGenJetForMet = commonGenObjects_->associatedGenJetForMET_;
            // Apply systematic variation to MET
            jetEnergyResolutionScaleFactors_->applyMetSystematic(v_jetForMet, met, v_jetForMetJerSF, v_associatedGenJetForMet);
        }
    }
    
    if(jetEnergyScaleScaleFactors_){
        // Get references for all relevant objects for jet variation
        VLV* v_jet = recoObjects_->jets_;
        // Apply systematic variation to jets
        jetEnergyScaleScaleFactors_->applyJetSystematic(v_jet);
        
        if(!mvaMet_){
            // Get references for all relevant objects for MET variation
            LV* met = recoObjects_->met_;
            VLV* v_jetForMet = commonGenObjects_->jetsForMET_;
            // Apply systematic variation to MET
            jetEnergyScaleScaleFactors_->applyMetSystematic(v_jetForMet, met);
        }
    }
    
    return *recoObjects_;
}



const CommonGenObjects& AnalysisBase::getCommonGenObjects(const Long64_t& entry)const
{
    if(!isMC_) return *commonGenObjects_;
    if(commonGenObjects_->valuesSet_) return *commonGenObjects_;

    this->GetCommonGenBranchesEntry(entry);
    return *commonGenObjects_;
}



const TopGenObjects& AnalysisBase::getTopGenObjects(const Long64_t& entry)const
{
    if(!isTopSignal_) return *topGenObjects_;
    if(topGenObjects_->valuesSet_) return *topGenObjects_;

    this->GetTopSignalBranchesEntry(entry);
    return *topGenObjects_;
}



const HiggsGenObjects& AnalysisBase::getHiggsGenObjects(const Long64_t& entry)const
{
    if(!isHiggsSignal_) return *higgsGenObjects_;
    if(higgsGenObjects_->valuesSet_) return *higgsGenObjects_;

    this->GetHiggsSignalBranchesEntry(entry);
    return *higgsGenObjects_;
}



const ZGenObjects& AnalysisBase::getZGenObjects(const Long64_t& entry)const
{
    if(!isDrellYan_) return *zGenObjects_;
    if(zGenObjects_->valuesSet_) return *zGenObjects_;

    this->GetZSignalBranchesEntry(entry);
    return *zGenObjects_;
}



void AnalysisBase::resetObjectStructEntry()const
{
    eventMetadata_->valuesSet_ = false;
    recoObjects_->valuesSet_ = false;
    commonGenObjects_->valuesSet_ = false;
    topGenObjects_->valuesSet_ = false;
    higgsGenObjects_->valuesSet_ = false;
    zGenObjects_->valuesSet_ = false;
}





// ------------------------------------- Methods for event and object selection -------------------------------------


void AnalysisBase::setBtagAlgorithmAndWorkingPoint(const Btag::Algorithm& algorithm,
                                                   const Btag::WorkingPoint& workingPoint)
{
    btagScaleFactors_->algorithmAndWorkingPoint(algorithm, workingPoint);
}



double AnalysisBase::btagCutValue()const
{
    return static_cast<double>(btagScaleFactors_->getWPDiscrValue());
}



bool AnalysisBase::failsDrellYanGeneratorSelection(const Long64_t& entry)const
{
    if(checkZDecayMode_ && !checkZDecayMode_(entry)) return true;
    return false;
}



bool AnalysisBase::failsTopGeneratorSelection(const Long64_t& entry)const
{
    if(!isTtbarPlusTauSample_) return false;
    this->GetTopDecayModeEntry(entry);
    //decayMode contains the decay of the top (*10) + the decay of the antitop
    //1=hadron, 2=e, 3=mu, 4=tau->hadron, 5=tau->e, 6=tau->mu
    //i.e. 23 == top decays to e, tbar decays to mu
    bool isViaTau = topDecayMode_ > 40 || (topDecayMode_ % 10 > 4);
    bool isCorrectChannel = false;
    switch (channelPdgIdProduct_) {
        case -11*13: isCorrectChannel = topDecayMode_ == 23 || topDecayMode_ == 32 //emu prompt
                        || topDecayMode_ == 53 || topDecayMode_ == 35 //e via tau, mu prompt
                        || topDecayMode_ == 26 || topDecayMode_ == 62 //e prompt, mu via tau
                        || topDecayMode_ == 56 || topDecayMode_ == 65; //both via tau
                        break;
        case -11*11: isCorrectChannel = topDecayMode_ == 22  //ee prompt
                        || topDecayMode_ == 52 || topDecayMode_ == 25 //e prompt, e via tau
                        || topDecayMode_ == 55; break; //both via tau
        case -13*13: isCorrectChannel = topDecayMode_ == 33  //mumu prompt
                        || topDecayMode_ == 36 || topDecayMode_ == 63 //mu prompt, mu via tau
                        || topDecayMode_ == 66; break; //both via tau
        default: std::cerr << "Invalid channel! Product = " << channelPdgIdProduct_ << "\n";
    };
    bool isBackgroundInSignalSample = !isCorrectChannel || isViaTau;
    if(runViaTau_ != isBackgroundInSignalSample) return true;
    return false;
}



bool AnalysisBase::hasLeptonPair(const int leadingLeptonIndex, const int nLeadingLeptonIndex,
                                 const std::vector<int>& lepPdgId)const
{
    bool hasLeptonPair(false);
    if(leadingLeptonIndex!=-1 && nLeadingLeptonIndex!=-1){
        // Check if lepton pair is correct flavour combination for the specified analysis channel (ee, emu, mumu)
        const int pdgIdProduct = lepPdgId.at(leadingLeptonIndex) * lepPdgId.at(nLeadingLeptonIndex);
        if(pdgIdProduct == channelPdgIdProduct_) hasLeptonPair = true;
    }
    return hasLeptonPair;
}



bool AnalysisBase::failsDileptonTrigger(const Long64_t& entry)const
{
    this->GetTriggerBranchesEntry(entry);

    //our triggers (bits: see the ntuplewriter!)
    constexpr int mumuTriggers = 0x8 + 0x20; //17/8 + 17Tr8
    constexpr int emuTriggers = 0x2000 + 0x4000;
    constexpr int eeTriggers = 0x40000;

    if (((triggerBits_ & mumuTriggers) && channelPdgIdProduct_ == -13*13)     // mumu triggers in rightmost byte
        || ((triggerBits_ & emuTriggers) && channelPdgIdProduct_ == -11*13)   // emu in 2nd byte
        || ((triggerBits_ & eeTriggers) && channelPdgIdProduct_ == -11*11))  // ee in 3rd byte
    {
        return false;
    }
    return true;
}



void AnalysisBase::mvaMet()
{
    mvaMet_ = true;
}

int AnalysisBase::topDecayMode(const Long64_t& entry)const
{
    this->GetTopDecayModeEntry(entry);
    
    return topDecayMode_;
}





// ----------------------- Methods for genJet selection, gen b/c jet identification and gen-reco jet matching -----------------------


void AnalysisBase::leptonCleanedJetIndices(std::vector<int>& jetIndices, const VLV& allJets, 
                                           const VLV& allLeptons, const double& minDeltaR)const
{
    // Not apply any cleaning if no leptons are given or minDeltaR set to negative value
    if(minDeltaR < 0.) return;
    if(!allLeptons.size()) return;
    
    std::vector<int> selectedIndices;
    for(int index : jetIndices){
        bool leptonNearby = false;
        const LV& jet = allJets.at(index);
        for(const LV& lepton : allLeptons){
            const double deltaR = ROOT::Math::VectorUtil::DeltaR(jet, lepton);
            if(deltaR >= minDeltaR) continue;
            leptonNearby = true;
            break;
        }
        if(leptonNearby) continue;
        selectedIndices.push_back(index);
    }
    
    jetIndices.clear();
    jetIndices = selectedIndices;
}



std::vector<std::vector<int> > AnalysisBase::matchHadronsToGenJets(const std::vector<int>& genJetIndices, const VLV& allGenJets, 
                                                                   const std::vector<int>& genHadJetIndices)const
{
    std::vector<std::vector<int> > result = std::vector<std::vector<int> >(allGenJets.size());
    
    for(size_t iHadron = 0; iHadron < genHadJetIndices.size(); ++iHadron){
        const int& jetIndex = genHadJetIndices.at(iHadron);
        // Protect against hadrons not clustered to any jet
        if(jetIndex < 0) continue;
        if(std::find(genJetIndices.begin(), genJetIndices.end(), jetIndex) == genJetIndices.end()) continue;
        result.at(jetIndex).push_back(iHadron);
    }
    
    return result;
}



std::vector<int> AnalysisBase::genBjetIndices(const std::vector<std::vector<int> >& genJetBhadronIndices)const
{
    std::vector<int> result;
    
    for(size_t iJet = 0; iJet < genJetBhadronIndices.size(); ++iJet){
        if(genJetBhadronIndices.at(iJet).size()) result.push_back(iJet);
    }
    
    return result;
}



std::vector<int> AnalysisBase::genCjetIndices(const std::vector<std::vector<int> >& genJetBhadronIndices,
                                              const std::vector<std::vector<int> >& genJetChadronIndices)const
{
    if(genJetBhadronIndices.size() != genJetChadronIndices.size()){
        std::cerr<<"ERROR in AnalysisBase::genCjetIndices()! Input vectors are of different size: "
                 <<genJetBhadronIndices.size()<<" , "<<genJetChadronIndices.size()
                 <<"\n...break\n"<<std::endl;
        exit(489);
    }
    
    std::vector<int> result;
    
    for(size_t iJet = 0; iJet < genJetChadronIndices.size(); ++iJet){
        // First exclude b jets, for remaining ones assign c jets
        if(genJetBhadronIndices.at(iJet).size()) continue;
        if(genJetChadronIndices.at(iJet).size()) result.push_back(iJet);
    }
    
    return result;
}



int AnalysisBase::genBjetIndex(const TopGenObjects& topGenObjects, const int pdgId)const
{
    int result(-1);
    
    const std::vector<int>& genBHadFlavour(*topGenObjects.genBHadFlavour_);
    const std::vector<int>& genBHadJetIndex(*topGenObjects.genBHadJetIndex_);
    
    bool alreadyFound(false);
    for(size_t iBHadron = 0; iBHadron < genBHadFlavour.size(); ++iBHadron){
        const int flavour = genBHadFlavour.at(iBHadron);
        if(flavour != pdgId) continue;
        
        // Assign jet index of corresponding hadron. Set to -2 if >1 hadrons found for the same flavour
        if(alreadyFound) return -2;
        result = genBHadJetIndex.at(iBHadron);
        alreadyFound = true;
    }
    
    return result;
}



int AnalysisBase::matchRecoToGenJet(const std::vector<int>& jetIndices, const VLV& jets,
                                    const int genJetIndex, const VLV& genJets)const
{
    if(genJetIndex < 0) return -1;
    const LV& genJet = genJets.at(genJetIndex);
    
    int result(-999);
    
    // Find closest jet and its distance in deltaR
    double deltaRJet(999.);
    for(const auto& index : jetIndices){
        double deltaR = ROOT::Math::VectorUtil::DeltaR(genJet, jets.at(index));
        if(deltaR < deltaRJet){
            deltaRJet = deltaR;
            result = index;
        }
    }
    
    // Call a jet matched if it is close enough, and has similar pt
    if(deltaRJet > 0.4) return -2;
    if(result >= 0){
        const double ptRecoJet = jets.at(result).pt();
        const double ptJet = genJet.pt();
        const double deltaPtRel = (ptJet - ptRecoJet)/ptJet;
        if(deltaPtRel<-0.5 || deltaPtRel>0.6) return -3;
    }
    
    return result;
}



std::vector<int> AnalysisBase::matchRecoToGenJets(const std::vector<int>& jetIndices, const VLV& jets,
                                                  const std::vector<int>& genJetIndices, const VLV& allGenJets)const
{
    // Set all values to -1
    std::vector<int> result = std::vector<int>(allGenJets.size(), -1);
    
    for(const int index : genJetIndices)
        result.at(index) = this->matchRecoToGenJet(jetIndices, jets, index, allGenJets);
    
    return result;
}





// ------------------------------------- Methods for application of corrections (e.g. scale factors) stored in ntuple -------------------------------------


double AnalysisBase::madgraphWDecayCorrection(const Long64_t& entry)const
{
    if(!correctMadgraphBR_) return 1.;
    this->GetTopDecayModeEntry(entry);

    // We must correct for the madGraph branching fraction being 1/9 for dileptons (PDG average is .108)
    if(topDecayMode_ == 11){ //all hadronic decay
        return (0.676*1.5) * (0.676*1.5);
    }
    else if(topDecayMode_< 20 || ( topDecayMode_ % 10 == 1)){ //semileptonic Decay
        return (0.108*9.) * (0.676*1.5);
    }
    else{ //dileptonic decay (including taus!)
        return (0.108*9.) * (0.108*9.);
    }
}



double AnalysisBase::weightPileup(const Long64_t& entry)const
{
    if(!isMC_ || !pileupScaleFactors_) return 1.;
    this->GetVertMultiTrueEntry(entry);
    return pileupScaleFactors_->getSF(vertMultiTrue_);
}



double AnalysisBase::weightGenerator(const Long64_t& entry)const
{
    if(!isMC_) return 1.;
    this->GetWeightGeneratorEntry(entry);
    return weightGenerator_;
}



double AnalysisBase::weightPdf(const Long64_t& entry, const int pdfNo)const
{
    if(pdfNo < 0) return 1.;
    this->GetPDFEntry(entry);
    const double pdfWeight = weightPDF_->at(pdfNo);
    return pdfWeight;
}






// ------------------------------------- Methods for application of corrections (e.g. scale factors) NOT stored in ntuple -------------------------------------


double AnalysisBase::weightTopPtReweighting(const Long64_t& entry)const
{
    // Apply only for ttbar dilepton sample
    if(!isTtbarPlusTauSample_) return 1.;
    // Do not apply if topPtScaleFactors_ is not set up
    if(!topPtScaleFactors_) return 1.;
    
    if(!topGenObjects_->valuesSet_) this->GetGenTopBranchesEntry(entry);
    
    return topPtScaleFactors_->getSF(topGenObjects_->GenTop_->pt(), topGenObjects_->GenAntiTop_->pt());
}



double AnalysisBase::weightLeptonSF(const int leadingLeptonIndex, const int nLeadingLeptonIndex,
                                    const VLV& allLeptons, const std::vector<int>& lepPdgId)const
{
    if(!isMC_) return 1.;
    if(leadingLeptonIndex<0 || nLeadingLeptonIndex<0) return 1.;
    return leptonScaleFactors_->getSFDilepton(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId);
}



double AnalysisBase::weightTriggerSF(const int leptonXIndex, const int leptonYIndex,
                                     const VLV& allLeptons)const
{
    if(!isMC_) return 1.;
    if(leptonXIndex<0 || leptonYIndex<0) return 1.;
    return triggerScaleFactors_->getSF(leptonXIndex, leptonYIndex, allLeptons);
}



double AnalysisBase::weightBtagSF(const std::vector<int>& jetIndices,
                                  const VLV& jets, const std::vector<int>& jetPartonFlavour,
                                  const std::vector<double>& btagDiscriminators)const
{
    if(!isMC_) return 1.;
    return btagScaleFactors_->getSF(jetIndices, jets, jetPartonFlavour, btagDiscriminators);
}



void AnalysisBase::retagJets(std::vector<int>& bjetIndices, const std::vector<int>& jetIndices,
                             const VLV& jets, const std::vector<int>& jetPartonFlavours,
                             const std::vector<double>& btagDiscriminants)const
{
    if(!isMC_) return;
    btagScaleFactors_->indexOfBtags(bjetIndices, jetIndices, jets, jetPartonFlavours, btagDiscriminants);
}



double AnalysisBase::weightKinReco()const{
    if(!isMC_) return 1.;
    if(!kinematicReconstructionScaleFactors_) return 1.;
    return kinematicReconstructionScaleFactors_->getSF();
}



bool AnalysisBase::makeBtagEfficiencies()const
{
    return btagScaleFactors_->makeEfficiencies();
}



void AnalysisBase::bookBtagEfficiencyHistos()
{
    if(!isMC_) return;
    if(isSampleForBtagEfficiencies_) btagScaleFactors_->bookEfficiencyHistograms(fOutput);
}



void AnalysisBase::fillBtagEfficiencyHistos(const std::vector<int>& jetIndices,
                                            const std::vector<double>& btagDiscriminators,
                                            const VLV& jets,
                                            const std::vector<int>& jetPartonFlavours,
                                            const double& weight)
{
    if(!isMC_) return;
    if(isSampleForBtagEfficiencies_)
        btagScaleFactors_->fillEfficiencyHistograms(jetIndices, btagDiscriminators, jets, jetPartonFlavours, weight);
}
    


void AnalysisBase::produceBtagEfficiencies()
{
    if(!isMC_) return;
    if(isSampleForBtagEfficiencies_) btagScaleFactors_->produceEfficiencies();
}



void AnalysisBase::correctMvaMet(const LV& dilepton, const int nJet, const Long64_t& entry)const
{
    if(!isDrellYan_ || !mvaMet_) return;
    
    if(!metRecoilCorrector_) return;
    
    const ZGenObjects& zGenObjects = this->getZGenObjects(entry);
    if(zGenObjects.GenZ_->size() != 1){
        std::cerr<<"ERROR in AnalysisBase::correctMvaMet()! Not exactly one Z stored, but: "<<zGenObjects.GenZ_->size()
                 <<"\n...break\n"<<std::endl;
        exit(292);
    }
    const LV& genZ = zGenObjects.GenZ_->at(0);
    
    LV& met = *recoObjects_->met_;
    
    metRecoilCorrector_->applyCorrection(met, genZ, dilepton, nJet);
}







// ------------------------------------- Various helper methods -------------------------------------


double AnalysisBase::getJetHT(const std::vector<int>& jetIndices, const VLV& jets)const
{
    double result = 0.;
    for(const int index : jetIndices){
        const double pt = jets.at(index).pt();
        result += pt;
    }
    return result;
}



KinematicReconstructionSolutions AnalysisBase::kinematicReconstructionSolutions(const int leptonIndex, const int antiLeptonIndex,
                                                                                const std::vector<int>& jetIndices, const std::vector<int>& bjetIndices,
                                                                                const VLV& allLeptons, const VLV& jets,
                                                                                const std::vector<double>& jetBTagCSV, const LV& met)const
{
    // If kinematic reconstruction is not initialised, do not run it but return dummy
    if(!kinematicReconstruction_) return KinematicReconstructionSolutions();
    
    return kinematicReconstruction_->solutions({leptonIndex}, {antiLeptonIndex}, jetIndices, bjetIndices,
                                               allLeptons, jets, jetBTagCSV, met);
}



const std::string AnalysisBase::topDecayModeString()const
{
    const std::vector<std::string> WMode {"unknown", "hadronic", "e", "mu", "tau->hadron", "tau->e", "tau->mu"};
    const int top = topDecayMode_ / 10;
    const int antitop = topDecayMode_ % 10;
    const std::string result = WMode[top] + "/" + WMode[antitop];
    return result;
}



int AnalysisBase::higgsDecayMode(const Long64_t& entry)const
{
    if(!isHiggsSignal_) return -1;
    this->GetHiggsDecayModeEntry(entry);
    
    return higgsDecayMode_;
}










