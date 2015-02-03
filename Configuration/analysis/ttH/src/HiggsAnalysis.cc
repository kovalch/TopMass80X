#define HiggsAnalysis_cxx

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <algorithm>

#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TSystem.h>
#include <Math/VectorUtil.h>

#include "HiggsAnalysis.h"
#include "higgsUtils.h"
#include "analysisStructs.h"
#include "AnalyzerBase.h"
#include "MvaTreeHandlerBase.h"
#include "MvaTreePlotterBase.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/classes.h"
#include "../../common/include/KinematicReconstructionSolution.h"
#include "../../common/include/ScaleFactors.h"





/// Lepton eta selection (absolute value)
constexpr double LeptonEtaCUT = 2.4;

/// Lepton pt selection in GeV
constexpr double LeptonPtCut = 20.;

/// Jet eta selection (absolute value)
constexpr double JetEtaCUT = 2.4;

/// Jet pt selection in GeV
constexpr double JetPtCUT = 20.;

/// Leading 2 jet pt selection in GeV
constexpr double Lead2JetPtCUT = 30.;

/// Minimal deltaR for removal of jets that are close to leptons (if negative, no cut applied)
constexpr double DeltaRLeptonJetCUT = -1.;

/// B-tag algorithm and working point
constexpr Btag::Algorithm BtagALGO = Btag::csv;
constexpr Btag::WorkingPoint BtagWP = Btag::M;

/// PF MET or MVA MET
constexpr bool MvaMET = true;

/// MET selection for same-flavour channels (ee, mumu)
constexpr double MetCUT = 40.;

/// Generated jet eta selection (absolute value)
constexpr double GenJetEtaCUT = 2.4;

/// Generated jet pt selection in GeV
constexpr double GenJetPtCUT = JetPtCUT;

/// Minimal deltaR for removal of generated jets that are close to leptons (if negative, no cut applied)
constexpr double GenDeltaRLeptonJetCUT = -1.;






HiggsAnalysis::HiggsAnalysis(TTree*):
inclusiveHiggsDecayMode_(-999),
additionalBjetMode_(-999),
reweightingName_(""),
reweightingSlope_(0.0)
{
    if(MvaMET) this->mvaMet();
}



HiggsAnalysis::~HiggsAnalysis()
{}



void HiggsAnalysis::Begin(TTree*)
{
    // Defaults from AnalysisBase
    AnalysisBase::Begin(0);
    
    // Set b-tagging working point
    this->setBtagAlgorithmAndWorkingPoint(BtagALGO, BtagWP);
    
    // Set up selection steps of MVA tree handlers
    for(MvaTreeHandlerBase* mvaTreeHandler : v_mvaTreeHandler_){
        if(mvaTreeHandler) mvaTreeHandler->book();
    }
}




void HiggsAnalysis::Terminate()
{
    // Produce b-tag efficiencies if required for given correction mode
    this->produceBtagEfficiencies();
    
    // Do everything needed for MVA
    for(MvaTreeHandlerBase* mvaTreeHandler : v_mvaTreeHandler_){
        if(mvaTreeHandler){
            // Produce and write tree
            mvaTreeHandler->writeTrees(this->outputFilename(), this->channel(), this->systematic());
            //mvaTreeHandler->writeTrees(fOutput);

            // Create and store control plots in fOutput
            MvaTreePlotterBase* mvaTreePlotter = mvaTreeHandler->setPlotter(mvaTreeHandler->stepMvaVariablesMap());
            mvaTreePlotter->plotVariables(fOutput);

            // Cleanup
            mvaTreePlotter->clear();
            delete mvaTreePlotter;
            mvaTreeHandler->clear();
        }
    }
    
    // Defaults from AnalysisBase
    AnalysisBase::Terminate();
}



void HiggsAnalysis::SlaveBegin(TTree *)
{
    // Defaults from AnalysisBase
    AnalysisBase::SlaveBegin(0);
    
    // Book histograms for b-tagging efficiencies if required for given correction mode
    this->bookBtagEfficiencyHistos();
    
    // Book histograms of all analyzers
    this->bookAll();
}



void HiggsAnalysis::SlaveTerminate()
{
    this->clearAll();

    // Defaults from AnalysisBase
    AnalysisBase::SlaveTerminate();
}



Bool_t HiggsAnalysis::Process(Long64_t entry)
{
    // Defaults from AnalysisBase
    if(!AnalysisBase::Process(entry)) return kFALSE;


    // Use utilities without namespaces
    using namespace common;


    // Entry for object structs are not yet read, so reset
    this->resetObjectStructEntry();


    // Define the selection steps as strings
    std::string selectionStep("");


    //===CUT===
    // this is step0a, no cut application
    selectionStep = "0a";
    
    // Create dummies for objects, non-dummies are created only as soon as needed
    const EventMetadata eventMetadataDummy;
    const RecoObjects recoObjectsDummy;
    const CommonGenObjects commonGenObjectsDummy;
    const TopGenObjects topGenObjectsDummy;
    const HiggsGenObjects higgsGenObjectsDummy;
    const KinematicReconstructionSolutions kinematicReconstructionSolutionsDummy;

    // Set up dummies for weights and indices, as needed for generic functions
    const tth::GenObjectIndices genObjectIndicesDummy({}, {}, {}, {}, {}, {}, {}, {}, {}, -1, -1, -1, -1, -1, -1, -1, -1);
    const tth::RecoObjectIndices recoObjectIndicesDummy({}, {}, {}, -1, -1, -1, -1, -1, -1, {}, {}, {});
    const tth::GenLevelWeights genLevelWeightsDummy(0., 0., 0., 0., 0., 0.);
    const tth::RecoLevelWeights recoLevelWeightsDummy(0., 0., 0., 0., 0., 0.);
    
    // ++++ Control Plots ++++
    
    this->fillAll(selectionStep,
                  eventMetadataDummy,
                  recoObjectsDummy, commonGenObjectsDummy,
                  topGenObjectsDummy, higgsGenObjectsDummy,
                  kinematicReconstructionSolutionsDummy,
                  genObjectIndicesDummy, recoObjectIndicesDummy,
                  genLevelWeightsDummy, recoLevelWeightsDummy,
                  1.);
    
    
    
    //===CUT===
    // this is step0b, select events on generator level and access true level weights
    selectionStep = "0b";
    
    // Separate DY dilepton decays in lepton flavours
    if(this->failsDrellYanGeneratorSelection(entry)) return kTRUE;
    
    // Separate dileptonic ttbar decays via tau
    //if(this->failsTopGeneratorSelection(entry)) return kTRUE;
    
    // Separate inclusive ttH sample in decays H->bbbar and others
    const int higgsDecayMode = this->higgsDecayMode(entry);
    if(this->failsHiggsGeneratorSelection(higgsDecayMode)) return kTRUE;
    
    // Separate tt+bb from tt+other
    const int topDecayMode = this->topDecayMode(entry);
    const int additionalJetFlavourId = this->additionalJetFlavourId(entry);
    if(this->failsAdditionalJetFlavourSelection(topDecayMode, additionalJetFlavourId)) return kTRUE;
    
    // Correct for the MadGraph branching fraction being 1/9 for dileptons (PDG average is .108)
    const double weightMadgraphCorrection = this->madgraphWDecayCorrection(entry);
    
    // Get weight due to pileup reweighting
    const double weightPU = this->weightPileup(entry);
    
    // Get weight due to generator weights
    const double weightGenerator = this->weightGenerator(entry);
    
    // Get weight due to top-pt reweighting
    const double weightTopPt = this->weightTopPtReweighting(entry);
    
    // Get true level weights
    const double trueLevelWeightNoPileup = weightTopPt*weightGenerator*weightMadgraphCorrection;
    const double trueLevelWeight = trueLevelWeightNoPileup*weightPU;
    
    const tth::GenLevelWeights genLevelWeights(weightMadgraphCorrection, weightPU,
                                               weightGenerator, weightTopPt,
                                               trueLevelWeightNoPileup, trueLevelWeight);
    
    // ++++ Control Plots ++++
    
    this->fillAll(selectionStep,
                  eventMetadataDummy,
                  recoObjectsDummy, commonGenObjectsDummy,
                  topGenObjectsDummy, higgsGenObjectsDummy,
                  kinematicReconstructionSolutionsDummy,
                  genObjectIndicesDummy, recoObjectIndicesDummy,
                  genLevelWeights, recoLevelWeightsDummy,
                  1.);
    
    
    
    //===CUT===
    selectionStep = "1";
    
    // Check if event was triggered with the same dilepton trigger as the specified analysis channel
    if(this->failsDileptonTrigger(entry)) return kTRUE;
    
    
    
    // === FULL RECO OBJECT SELECTION === (can thus be used at each selection step)
    
    // Access reco objects, and common generator objects
    const RecoObjects& recoObjects = this->getRecoObjects(entry);
    const CommonGenObjects& commonGenObjects = this->getCommonGenObjects(entry);
    
    // Get allLepton indices, apply selection cuts and order them by pt (beginning with the highest value)
    const VLV& allLeptons = *recoObjects.allLeptons_;
    const std::vector<int>& lepPdgId = *recoObjects.lepPdgId_;
    std::vector<int> allLeptonIndices = initialiseIndices(allLeptons);
    selectIndices(allLeptonIndices, allLeptons, LVeta, LeptonEtaCUT, false);
    selectIndices(allLeptonIndices, allLeptons, LVeta, -LeptonEtaCUT);
    selectIndices(allLeptonIndices, allLeptons, LVpt, LeptonPtCut);
    orderIndices(allLeptonIndices, allLeptons, LVpt);
    //const int numberOfAllLeptons = allLeptonIndices.size();
    
    // Get indices of leptons and antiLeptons separated by charge, and get the leading ones if they exist
    std::vector<int> leptonIndices = allLeptonIndices;
    std::vector<int> antiLeptonIndices = allLeptonIndices;
    selectIndices(leptonIndices, lepPdgId, 0);
    selectIndices(antiLeptonIndices, lepPdgId, 0, false);
    const int numberOfLeptons = leptonIndices.size();
    const int numberOfAntiLeptons = antiLeptonIndices.size();
    const int leptonIndex = numberOfLeptons>0 ? leptonIndices.at(0) : -1;
    const int antiLeptonIndex = numberOfAntiLeptons>0 ? antiLeptonIndices.at(0) : -1;
    
    // In case of an existing opposite-charge dilepton system,
    // get their indices for leading and next-to-leading lepton
    int leadingLeptonIndex(-1);
    int nLeadingLeptonIndex(-1);
    if(numberOfLeptons>0 && numberOfAntiLeptons>0){
        leadingLeptonIndex = leptonIndex;
        nLeadingLeptonIndex = antiLeptonIndex;
        orderIndices(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, LVpt);
    }
    const bool hasLeptonPair = this->hasLeptonPair(leadingLeptonIndex, nLeadingLeptonIndex, lepPdgId);
    
    // Get two indices of the two leptons in the right order for trigger scale factor, if existing
    int leptonXIndex(leadingLeptonIndex);
    int leptonYIndex(nLeadingLeptonIndex);
    if(hasLeptonPair){
        //in ee and mumu channel leptonX must be the highest pt lepton, i.e. this is already correct
        // in emu channel leptonX must be electron
        if(std::abs(lepPdgId.at(leptonXIndex)) != std::abs(lepPdgId.at(leptonYIndex))){
            orderIndices(leptonYIndex, leptonXIndex, lepPdgId, true);
        }
    }
    
    // Get dilepton system, if existing
    const LV dummyLV(0.,0.,0.,0.);
    const LV dilepton(hasLeptonPair ? allLeptons.at(leadingLeptonIndex)+allLeptons.at(nLeadingLeptonIndex) : dummyLV);
    
    // Get jet indices, apply selection cuts and order them by pt (beginning with the highest value)
    const VLV& jets = *recoObjects.jets_;
    std::vector<int> jetIndices = initialiseIndices(jets);
    selectIndices(jetIndices, jets, LVeta, JetEtaCUT, false);
    selectIndices(jetIndices, jets, LVeta, -JetEtaCUT);
    selectIndices(jetIndices, jets, LVpt, JetPtCUT);
    if(DeltaRLeptonJetCUT > 0.){
        // Vector of leptons from which jets need to be separated in deltaR
        VLV leptonsForJetCleaning;
        for(const int index : allLeptonIndices) leptonsForJetCleaning.push_back(allLeptons.at(index));
        this->leptonCleanedJetIndices(jetIndices, jets, leptonsForJetCleaning, DeltaRLeptonJetCUT);
    }
    orderIndices(jetIndices, jets, LVpt);
    const int numberOfJets = jetIndices.size();
    const bool has2Jets = numberOfJets > 1 && jets.at(jetIndices.at(1)).pt() >= Lead2JetPtCUT;
    
    // Fill a vector with all jet pair indices, while sorting each pair by the jet charge:
    // first entry is antiBIndex i.e. with higher jet charge, second entry is bIndex
    //const std::vector<double>& jetChargeGlobalPtWeighted = *recoObjects.jetChargeGlobalPtWeighted_;
    const std::vector<double>& jetChargeRelativePtWeighted = *recoObjects.jetChargeRelativePtWeighted_;
    const tth::IndexPairs& jetIndexPairs = this->chargeOrderedJetPairIndices(jetIndices, jetChargeRelativePtWeighted);
    
    // Get b-jet indices, apply selection cuts
    // and apply b-tag efficiency MC correction using random number based tag flipping (if requested correction mode is applied)
    // and order b-jets by btag discriminator (beginning with the highest value)
    const std::vector<double>& jetBTagCSV = *recoObjects.jetBTagCSV_;
    const std::vector<int>& jetPartonFlavour = *commonGenObjects.jetPartonFlavour_;
    std::vector<int> bjetIndices = jetIndices;
    selectIndices(bjetIndices, jetBTagCSV, this->btagCutValue());
    this->retagJets(bjetIndices, jetIndices, jets, jetPartonFlavour, jetBTagCSV);
    orderIndices(bjetIndices, jetBTagCSV);
    const int numberOfBjets = bjetIndices.size();
    const bool hasBtag = numberOfBjets > 0;
    
    // Get MET, and in case of MVA MET apply recoil correction for Drell-Yan sample
    this->correctMvaMet(dilepton, numberOfJets, entry);
    const LV& met = *recoObjects.met_;
    const bool hasMet = met.pt() > MetCUT;
    
    const tth::RecoObjectIndices recoObjectIndices(allLeptonIndices,
                                                   leptonIndices, antiLeptonIndices,
                                                   leptonIndex, antiLeptonIndex,
                                                   leadingLeptonIndex, nLeadingLeptonIndex,
                                                   leptonXIndex, leptonYIndex,
                                                   jetIndices, jetIndexPairs,
                                                   bjetIndices);
    
    
    // Determine all reco level weights
    const double weightLeptonSF = this->weightLeptonSF(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId);
    const double weightTriggerSF = this->weightTriggerSF(leptonXIndex, leptonYIndex, allLeptons);
    const double weightNoPileup = trueLevelWeightNoPileup*weightTriggerSF*weightLeptonSF;
    const double weightBtagSF = this->weightBtagSF(jetIndices, jets, jetPartonFlavour, jetBTagCSV);
    const double weightKinReco = this->weightKinReco();
    
    // The weight to be used for filling the histograms
    double weight = weightNoPileup*weightPU;
    
    
    tth::RecoLevelWeights recoLevelWeights(weightLeptonSF, weightTriggerSF,
                                           weightBtagSF, weightKinReco,
                                           weightNoPileup, weight);
    
    
    // ++++ Control Plots ++++
    
    this->fillAll(selectionStep,
                  eventMetadataDummy,
                  recoObjects, commonGenObjects,
                  topGenObjectsDummy, higgsGenObjectsDummy,
                  kinematicReconstructionSolutionsDummy,
                  genObjectIndicesDummy, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  1.);
    
    
    
    //===CUT===
    selectionStep = "2";
    
    // we need an OS lepton pair matching the trigger selection...
    if (!hasLeptonPair) return kTRUE;
    
    // ++++ Control Plots ++++
    
    this->fillAll(selectionStep,
                  eventMetadataDummy,
                  recoObjects, commonGenObjects,
                  topGenObjectsDummy, higgsGenObjectsDummy,
                  kinematicReconstructionSolutionsDummy,
                  genObjectIndicesDummy, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);
    
    
    
    //===CUT===
    selectionStep = "3";
    
    // ...with at least 20 GeV invariant mass
    if(dilepton.M() < 20.) return kTRUE;
    
    // ++++ Control Plots ++++

    this->fillAll(selectionStep,
                  eventMetadataDummy,
                  recoObjects, commonGenObjects,
                  topGenObjectsDummy, higgsGenObjectsDummy,
                  kinematicReconstructionSolutionsDummy,
                  genObjectIndicesDummy, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);
    
    
    
    //=== CUT ===
    selectionStep = "4";
    
    // Exclude the Z window in analysis cutflow, but keep these events for Drell-Yan corrections
    const bool isZregion = dilepton.M() > 76. && dilepton.M() < 106.;
    const bool isEmu = this->channel() == Channel::emu;
    
    // ++++ Z-window plots ++++
    
    if(isZregion){
        this->fillAll("4zWindow",
                      eventMetadataDummy,
                      recoObjects, commonGenObjects,
                      topGenObjectsDummy, higgsGenObjectsDummy,
                      kinematicReconstructionSolutionsDummy,
                      genObjectIndicesDummy, recoObjectIndices,
                      genLevelWeights, recoLevelWeights,
                      weight);
    }
    
    // ++++ Control Plots ++++
    
    if(isEmu || !isZregion){
        this->fillAll(selectionStep,
                      eventMetadataDummy,
                      recoObjects, commonGenObjects,
                      topGenObjectsDummy, higgsGenObjectsDummy,
                      kinematicReconstructionSolutionsDummy,
                      genObjectIndicesDummy, recoObjectIndices,
                      genLevelWeights, recoLevelWeights,
                      weight);
    }
    
    
    
    //=== CUT ===
    selectionStep = "5";
    
    //Require at least two jets
    if(!has2Jets) return kTRUE;
    
    // ++++ Z-window plots ++++
    
    if(isZregion){
        this->fillAll("5zWindow",
                      eventMetadataDummy,
                      recoObjects, commonGenObjects,
                      topGenObjectsDummy, higgsGenObjectsDummy,
                      kinematicReconstructionSolutionsDummy,
                      genObjectIndicesDummy, recoObjectIndices,
                      genLevelWeights, recoLevelWeights,
                      weight);
    }
    
    // ++++ Control Plots ++++
    
    if(isEmu || !isZregion){
        this->fillAll(selectionStep,
                      eventMetadataDummy,
                      recoObjects, commonGenObjects,
                      topGenObjectsDummy, higgsGenObjectsDummy,
                      kinematicReconstructionSolutionsDummy,
                      genObjectIndicesDummy, recoObjectIndices,
                      genLevelWeights, recoLevelWeights,
                      weight);
    }
    
    
    
    //=== CUT ===
    selectionStep = "6";
    
    //Require MET > 40 GeV in non-emu channels
    if(!(hasMet || isEmu)) return kTRUE;
    
    // ++++ Z-window plots ++++
    
    if(isZregion){
        this->fillAll("6zWindow",
                      eventMetadataDummy,
                      recoObjects, commonGenObjects,
                      topGenObjectsDummy, higgsGenObjectsDummy,
                      kinematicReconstructionSolutionsDummy,
                      genObjectIndicesDummy, recoObjectIndices,
                      genLevelWeights, recoLevelWeights,
                      weight);
    }
    
    // ++++ Control Plots ++++
    
    if(isEmu || !isZregion){
        this->fillAll(selectionStep,
                      eventMetadataDummy,
                      recoObjects, commonGenObjects,
                      topGenObjectsDummy, higgsGenObjectsDummy,
                      kinematicReconstructionSolutionsDummy,
                      genObjectIndicesDummy, recoObjectIndices,
                      genLevelWeights, recoLevelWeights,
                      weight);
        
        // Fill b-tagging efficiencies if required for given correction mode
        this->fillBtagEfficiencyHistos(jetIndices, jetBTagCSV, jets, jetPartonFlavour, weight);
    }
    
    // In case of filling b-tagging efficiencies, do not process further steps
    if(this->makeBtagEfficiencies()) return kTRUE;
    
    
    
    //=== CUT ===
    selectionStep = "7";
    
    //Require at least one b tagged jet
    if(!hasBtag) return kTRUE;
    
    weight *= weightBtagSF;
    
    // ++++ Z-window plots ++++
    
    if(isZregion){
        this->fillAll("7zWindow",
                      eventMetadataDummy,
                      recoObjects, commonGenObjects,
                      topGenObjectsDummy, higgsGenObjectsDummy,
                      kinematicReconstructionSolutionsDummy,
                      genObjectIndicesDummy, recoObjectIndices,
                      genLevelWeights, recoLevelWeights,
                      weight);
    }
    
    // Gen-level info and kinematic reconstruction not needed for Z region
    if(!isEmu && isZregion) return kTRUE;
    
    
    
    // Access event meta data
    //const EventMetadata eventMetadata = eventMetadataDummy;
    const EventMetadata eventMetadata = this->getEventMetadata(entry);
    
    // Access kinematic reconstruction info
    //const KinematicReconstructionSolutions kinematicReconstructionSolutions = kinematicReconstructionSolutionsDummy;
    const KinematicReconstructionSolutions kinematicReconstructionSolutions = this->kinematicReconstructionSolutions(leptonIndex, antiLeptonIndex, jetIndices, bjetIndices, allLeptons, jets, jetBTagCSV, met);
    //const bool hasSolution = kinematicReconstructionSolutions.numberOfSolutions();
    //std::cout<<"\n\n\nNew event - solutions(): "<<kinematicReconstructionSolutions.numberOfSolutions()<<"\n";
    //for(size_t iSolution = 0; iSolution < kinematicReconstructionSolutions.numberOfSolutions(); ++iSolution){
    //    kinematicReconstructionSolutions.solution(KinematicReconstructionSolution::averagedSumSmearings_mlb, iSolution).print();
    //    std::cout<<"\n";
    //}
    
    
    
    // === FULL GEN OBJECT SELECTION ===
    
    // Access top generator object struct, and higgs generator object struct
    const TopGenObjects& topGenObjects = this->getTopGenObjects(entry);
    const HiggsGenObjects& higgsGenObjects = this->getHiggsGenObjects(entry);
    
    // Generated jets
    const VLV& allGenJets =  topGenObjects.valuesSet_ ? *topGenObjects.allGenJets_ : VLV();
    std::vector<int> allGenJetIndices = initialiseIndices(allGenJets);
    std::vector<int> genJetIndices = allGenJetIndices;
    selectIndices(genJetIndices, allGenJets, LVeta, GenJetEtaCUT, false);
    selectIndices(genJetIndices, allGenJets, LVeta, -GenJetEtaCUT);
    selectIndices(genJetIndices, allGenJets, LVpt, GenJetPtCUT);
    if(GenDeltaRLeptonJetCUT > 0.){
        // Vector of genLeptons from which genJets need to be separated in deltaR
        VLV allGenLeptons;
        if(topGenObjects.valuesSet_){
            if(topGenObjects.GenLepton_) allGenLeptons.push_back(*topGenObjects.GenLepton_);
            if(topGenObjects.GenAntiLepton_) allGenLeptons.push_back(*topGenObjects.GenAntiLepton_);
        }
        this->leptonCleanedJetIndices(genJetIndices, allGenJets, allGenLeptons, GenDeltaRLeptonJetCUT);
    }
    
    // Match for all genJets all B hadrons
    std::vector<std::vector<int> > allGenJetBhadronIndices;
    std::vector<std::vector<int> > genJetBhadronIndices;
    std::vector<int> allGenBjetIndices;
    std::vector<int> genBjetIndices;
    std::vector<int> genJetMatchedRecoBjetIndices;
    if(topGenObjects.valuesSet_){
        allGenJetBhadronIndices = this->matchHadronsToGenJets(allGenJetIndices, allGenJets, *topGenObjects.genBHadJetIndex_);
        genJetBhadronIndices = this->matchHadronsToGenJets(genJetIndices, allGenJets, *topGenObjects.genBHadJetIndex_);
        allGenBjetIndices = this->genBjetIndices(allGenJetBhadronIndices);
        genBjetIndices = this->genBjetIndices(genJetBhadronIndices);
        genJetMatchedRecoBjetIndices = this->matchRecoToGenJets(jetIndices, jets, allGenBjetIndices, allGenJets);
    }
    
    // Match for all genJets all C hadrons
    std::vector<std::vector<int> > allGenJetChadronIndices;
    std::vector<std::vector<int> > genJetChadronIndices;
    std::vector<int> allGenCjetIndices;
    std::vector<int> genCjetIndices;
    std::vector<int> genJetMatchedRecoCjetIndices;
    if(topGenObjects.valuesSet_){
        allGenJetChadronIndices = this->matchHadronsToGenJets(allGenJetIndices, allGenJets, *topGenObjects.genCHadJetIndex_);
        genJetChadronIndices = this->matchHadronsToGenJets(genJetIndices, allGenJets, *topGenObjects.genCHadJetIndex_);
        allGenCjetIndices = this->genCjetIndices(allGenJetBhadronIndices, allGenJetChadronIndices);
        genCjetIndices = this->genCjetIndices(genJetBhadronIndices, genJetChadronIndices);
        genJetMatchedRecoCjetIndices = this->matchRecoToGenJets(jetIndices, jets, allGenCjetIndices, allGenJets);
    }
    
    // Jet matchings for ttbar system
    int genBjetFromTopIndex(-1);
    int genAntiBjetFromTopIndex(-1);
    int matchedBjetFromTopIndex(-1);
    int matchedAntiBjetFromTopIndex(-1);
    if(topGenObjects.valuesSet_){
        genBjetFromTopIndex = this->genBjetIndex(topGenObjects, 6);
        genAntiBjetFromTopIndex = this->genBjetIndex(topGenObjects, -6);
        matchedBjetFromTopIndex = genBjetFromTopIndex>=0 ? genJetMatchedRecoBjetIndices.at(genBjetFromTopIndex) : -1;
        matchedAntiBjetFromTopIndex = genAntiBjetFromTopIndex>=0 ? genJetMatchedRecoBjetIndices.at(genAntiBjetFromTopIndex) : -1;
    }
    
    // Jet matchings for Higgs system
    int genBjetFromHiggsIndex(-1);
    int genAntiBjetFromHiggsIndex(-1);
    int matchedBjetFromHiggsIndex(-1);
    int matchedAntiBjetFromHiggsIndex(-1);
    if(topGenObjects.valuesSet_ && higgsDecayMode == 5){
        genBjetFromHiggsIndex = this->genBjetIndex(topGenObjects, 25);
        genAntiBjetFromHiggsIndex = this->genBjetIndex(topGenObjects, -25);
        matchedBjetFromHiggsIndex = genBjetFromHiggsIndex>=0 ? genJetMatchedRecoBjetIndices.at(genBjetFromHiggsIndex) : -1;
        matchedAntiBjetFromHiggsIndex = genAntiBjetFromHiggsIndex>=0 ? genJetMatchedRecoBjetIndices.at(genAntiBjetFromHiggsIndex) : -1;
    }
    
    
    const tth::GenObjectIndices genObjectIndices(genJetIndices,
                                                 allGenBjetIndices,
                                                 genBjetIndices,
                                                 genJetBhadronIndices,
                                                 genJetMatchedRecoBjetIndices,
                                                 allGenCjetIndices,
                                                 genCjetIndices,
                                                 genJetChadronIndices,
                                                 genJetMatchedRecoCjetIndices,
                                                 genBjetFromTopIndex, genAntiBjetFromTopIndex,
                                                 matchedBjetFromTopIndex, matchedAntiBjetFromTopIndex,
                                                 genBjetFromHiggsIndex, genAntiBjetFromHiggsIndex,
                                                 matchedBjetFromHiggsIndex, matchedAntiBjetFromHiggsIndex);
    
        
    // Applying the reweighting
    weight *= this->reweightingWeight(topGenObjects, genObjectIndices);
    
    
    // ++++ Control Plots ++++
    
    this->fillAll(selectionStep,
                  eventMetadata,
                  recoObjects, commonGenObjects,
                  topGenObjects, higgsGenObjects,
                  kinematicReconstructionSolutions,
                  genObjectIndices, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);
    
    
    
    return kTRUE;
}



tth::IndexPairs HiggsAnalysis::chargeOrderedJetPairIndices(const std::vector<int>& jetIndices,
                                                           const std::vector<double>& jetCharges)
{
    tth::IndexPairs result;
    if(jetIndices.size() < 2) return result;

    // Loop over all jet combinations
    for(std::vector<int>::const_iterator i_jetIndex = jetIndices.begin(); i_jetIndex != --(jetIndices.end()); ++i_jetIndex){
        std::vector<int>::const_iterator incrementIterator(i_jetIndex);
        ++incrementIterator;
        for(std::vector<int>::const_iterator j_jetIndex = incrementIterator; j_jetIndex != jetIndices.end(); ++j_jetIndex){
            // Get the indices of b and anti-b jet defined by jet charge
            int bIndex = *i_jetIndex;
            int antiBIndex = *j_jetIndex;
            common::orderIndices(antiBIndex, bIndex, jetCharges);

            result.push_back(std::make_pair(antiBIndex, bIndex));
        }
    }

    return result;
}



void HiggsAnalysis::SetInclusiveHiggsDecayMode(const int inclusiveHiggsDecayMode)
{
    inclusiveHiggsDecayMode_ = inclusiveHiggsDecayMode;
}



void HiggsAnalysis::SetAdditionalBjetMode(const int additionalBjetMode)
{
    additionalBjetMode_ = additionalBjetMode;
}



void HiggsAnalysis::SetReweightingName(const TString reweightingName)
{
    reweightingName_ = reweightingName;
}



void HiggsAnalysis::SetReweightingSlope(const double reweightingSlope)
{
    reweightingSlope_ = reweightingSlope;
}



void HiggsAnalysis::SetAllAnalyzers(std::vector<AnalyzerBase*> v_analyzer)
{
    v_analyzer_ = v_analyzer;
}



void HiggsAnalysis::SetAllTreeHandlers(std::vector<MvaTreeHandlerBase*> v_mvaTreeHandler)
{
    v_mvaTreeHandler_ = v_mvaTreeHandler;
}



bool HiggsAnalysis::failsHiggsGeneratorSelection(const int higgsDecayMode)const
{
    if(inclusiveHiggsDecayMode_ == -999) return false;
    
    // Separate ttH events from inclusve decay into H->bbbar and other decays
    if(inclusiveHiggsDecayMode_==0 && higgsDecayMode==5) return true;
    if(inclusiveHiggsDecayMode_==5 && higgsDecayMode!=5) return true;
    return false;
}



bool HiggsAnalysis::failsAdditionalJetFlavourSelection(const int topDecayMode, const int additionalJetFlavourId)const
{
    // Use the full ttbar dilepton sample for creating btag efficiencies
    if(this->makeBtagEfficiencies()) return false;
    
    // Check steering parameter if any separation is requested
    int additionalBjetMode(additionalBjetMode_);
    if(additionalBjetMode == -999) return false;
    
    // topDecayMode contains the decay of the top (*10) + the decay of the antitop (plus 100 or 200 for non-b decays of tops)
    // 1=hadron, 2=e, 3=mu, 4=tau->hadron, 5=tau->e, 6=tau->mu
    // i.e. 23 == top decays to e, tbar decays to mu
    const int decay = topDecayMode % 100;
    
    // Extract part telling how leptonic tau decays should be treated
    const int viaTauMode = additionalBjetMode / 100;
    if(viaTauMode == 0){
        // any dileptonic final state allowed
    }
    else if(viaTauMode == 1){
        // both leptons from W->e/mu
        if(!(decay/10 < 4 && decay%10 < 4)) return true;
    }
    else if(viaTauMode == 2){
        // at least 1 lepton from W->tau->e/mu
        if(!(decay/10 > 4 || decay%10 > 4)) return true;
    }
    else{
        std::cerr<<"ERROR in HiggsAnalysis::failsAdditionalJetFlavourSelection()! Undefined tau mode requested: "
                 <<viaTauMode<<"\n...break\n"<<std::endl;
        exit(421);
    }
    
    // additionalJetFlavourId contains number of jets from top (*100)
    // + the second digit encodes the flavour of the additional jets and whether they stem from before or after the top weak decay
    // + the last digit encodes the number of jets of this flavour, and the number of hadrons contained in them
    const int jetId = additionalJetFlavourId % 100;
    
    // Leave only part telling which additional jets should be present
    additionalBjetMode %= 100;
    if(additionalBjetMode == 4){
        // tt+c (tt+cc)
        if(jetId>20 && jetId<30) return false;
    }
    else if(additionalBjetMode == 3){
        // tt+bb
        if(jetId==3 || jetId==4) return false;
    }
    else if(additionalBjetMode == 2){
        // tt+2b
        if(jetId==2) return false;
    }
    else if(additionalBjetMode == 1){
        // tt+b
        if(jetId==1) return false;
    }
    else if(additionalBjetMode == 0){
        // tt+other
        if(jetId==0 || (jetId>4 && jetId<20) ||  jetId>30 ) return false;
    }
    else{
        std::cerr<<"ERROR in HiggsAnalysis::failsAdditionalJetFlavourSelection()! Undefined additional jet mode requested: "
                 <<additionalBjetMode<<"\n...break\n"<<std::endl;
        exit(422);
    }
    
    return true;
}



double HiggsAnalysis::reweightingWeight(const TopGenObjects& topGenObjects, const tth::GenObjectIndices& genObjectIndices)const
{
    if(reweightingName_ == "") return 1.0;
    if(reweightingSlope_ == 0.0) return 1.0;
    if(reweightingName_ == "nominal") return 1.0;
    
    // Getting generator level information
    const VLV& allGenJets = (topGenObjects.valuesSet_) ? *topGenObjects.allGenJets_ : VLV(0);
    std::vector<int> genJetsId = genObjectIndices.genJetIndices_;
    std::vector<int> genBJetsId = genObjectIndices.genBjetIndices_;
    std::vector<int> topBJetsId_gen;
    std::vector<int> addBJetsId_gen;
    
    if(genJetsId.size() > allGenJets.size() || genJetsId.size() < 1) {
        std::cerr << "ERROR! LorentzVectors of genJets not available for reweighting. Stopping..." << std::endl;
        exit(1);
    }
    
    // Sorting bjets by Pt
    common::orderIndices(addBJetsId_gen, allGenJets, common::LVpt);
    
    // Selecting additional b jets
    for(int jetId : genBJetsId) {
        if(jetId == genObjectIndices.genBjetFromTopIndex_) {
            topBJetsId_gen.push_back(jetId);
            continue;
        }
        if(jetId == genObjectIndices.genAntiBjetFromTopIndex_) {
            topBJetsId_gen.push_back(jetId);
            continue;
        }
        addBJetsId_gen.push_back(jetId);
    }
    
    
    if(reweightingName_ == "1st_add_bjet_pt") {
        if(addBJetsId_gen.size() < 1) return 1.0;
        const int jetId = addBJetsId_gen.at(0);
        return 1 + reweightingSlope_*(allGenJets.at(jetId).Pt() - 100.)/100.;
    }
    
    if(reweightingName_ == "1st_add_bjet_eta") {
        if(addBJetsId_gen.size() < 1) return 1.0;
        const int jetId = addBJetsId_gen.at(0);
        return 1. + reweightingSlope_*(std::fabs(allGenJets.at(jetId).Eta()) - 1.);
    }
    
    if(reweightingName_ == "2nd_add_bjet_pt") {
        if(addBJetsId_gen.size() < 2) return 1.0;
        const int jetId = addBJetsId_gen.at(1);
        return 1. + reweightingSlope_*(allGenJets.at(jetId).Pt() - 60.)/60.;
    }
    
    if(reweightingName_ == "2nd_add_bjet_eta") {
        if(addBJetsId_gen.size() < 2) return 1.0;
        const int jetId = addBJetsId_gen.at(1);
        return 1. + reweightingSlope_*(std::fabs(allGenJets.at(jetId).Eta()) - 1.);
    }
    
    if(reweightingName_ == "add_bjet_dR") {
        if(addBJetsId_gen.size() < 2) return 1.0;
        const int jetId_1 = addBJetsId_gen.at(0);
        const int jetId_2 = addBJetsId_gen.at(1);
        return 1. + reweightingSlope_*(ROOT::Math::VectorUtil::DeltaR(allGenJets.at(jetId_1), allGenJets.at(jetId_2)) - 1.5)/1.5;
    }
    
    if(reweightingName_ == "add_bjet_Mjj") {
        if(addBJetsId_gen.size() < 2) return 1.0;
        const int jetId_1 = addBJetsId_gen.at(0);
        const int jetId_2 = addBJetsId_gen.at(1);
        return 1. + reweightingSlope_*((allGenJets.at(jetId_1) + allGenJets.at(jetId_2)).M() - 100.)/100.;
    }
    
    std::cerr<<"ERROR in HiggsAnalysis::reweightingWeight()! Provided reweighting name is not supported: "<<reweightingName_<<"\n...break\n"<<std::endl;
    exit(1);

    return 1.0;
}



void HiggsAnalysis::fillAll(const std::string& selectionStep,
                            const EventMetadata& eventMetadata,
                            const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                            const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                            const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                            const tth::GenObjectIndices& genObjectIndices, const tth::RecoObjectIndices& recoObjectIndices,
                            const tth::GenLevelWeights& genLevelWeights, const tth::RecoLevelWeights& recoLevelWeights,
                            const double& defaultWeight)const
{
    // In case b-tag efficiencies are produced, analysis output is not
    if(this->makeBtagEfficiencies()) return;
    
    for(AnalyzerBase* analyzer : v_analyzer_){
        if(analyzer) analyzer->fill(eventMetadata,
                                    recoObjects, commonGenObjects,
                                    topGenObjects, higgsGenObjects,
                                    kinematicReconstructionSolutions,
                                    recoObjectIndices, genObjectIndices,
                                    genLevelWeights, recoLevelWeights,
                                    defaultWeight, selectionStep);
    }
    
    for(MvaTreeHandlerBase* mvaTreeHandler : v_mvaTreeHandler_){
        if(mvaTreeHandler) mvaTreeHandler->fill(eventMetadata,
                                                recoObjects, commonGenObjects,
                                                topGenObjects, higgsGenObjects,
                                                kinematicReconstructionSolutions,
                                                recoObjectIndices, genObjectIndices,
                                                genLevelWeights, recoLevelWeights,
                                                defaultWeight, selectionStep);
    }
}



void HiggsAnalysis::bookAll()
{
    for(AnalyzerBase* analyzer : v_analyzer_){
        if(analyzer) analyzer->book(fOutput);
    }
}



void HiggsAnalysis::clearAll()
{
    for(AnalyzerBase* analyzer : v_analyzer_){
        if(analyzer) analyzer->clear();
    }
}








