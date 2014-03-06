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
#include "AnalyzerBaseClass.h"
#include "MvaTreeHandlerBase.h"
#include "MvaTreePlotterBase.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/classes.h"
#include "../../common/include/ScaleFactors.h"





/// Lepton eta selection (absolute value)
constexpr double LeptonEtaCUT = 2.4;

/// Lepton pt selection in GeV
constexpr double LeptonPtCut = 20.;

/// Jet eta selection (absolute value)
constexpr double JetEtaCUT = 2.4;

/// Jet pt selection in GeV
constexpr double JetPtCUT = 30.;

/// Leading 2 jet pt selection in GeV (For cut based approach)
constexpr double Lead2JetPtCUT = JetPtCUT;


/// B-tag working point
/// Available options: 
///  csvl_wp (0.244) 
///  csvm_wp (0.679)
///  csvt_wp (0.898)
constexpr BtagScaleFactors::workingPoints BtagWP = BtagScaleFactors::csvl_wp;


/// MET selection for same-flavour channels (ee, mumu)
constexpr double MetCUT = 40.;










HiggsAnalysis::HiggsAnalysis(TTree*):
inclusiveHiggsDecayMode_(-999),
additionalBjetMode_(-999)
{}



HiggsAnalysis::~HiggsAnalysis()
{}



void HiggsAnalysis::Begin(TTree*)
{
    // Defaults from AnalysisBase
    AnalysisBase::Begin(0);

    // Set up selection steps of MVA tree handlers
    for(MvaTreeHandlerBase* mvaTreeHandler : v_mvaTreeHandler_){
        if(mvaTreeHandler) mvaTreeHandler->book();
    }
}




void HiggsAnalysis::Terminate()
{
    // Produce b-tag efficiencies
    // FIXME: Shouldn't we also clear b-tagging efficiency histograms if they are produced ?
    if(this->makeBtagEfficiencies()) btagScaleFactors_->produceBtagEfficiencies(static_cast<std::string>(this->channel()));

    // Do everything needed for MVA
    for(MvaTreeHandlerBase* mvaTreeHandler : v_mvaTreeHandler_){
        if(mvaTreeHandler){
            // Produce and write tree
            mvaTreeHandler->writeTrees(static_cast<std::string>(this->outputFilename()),
                                        Channel::convertChannel(static_cast<std::string>(this->channel())),
                                        Systematic::convertSystematic(static_cast<std::string>(this->systematic())));
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
    
    // Set b-tagging working point
    btagScaleFactors_->setWorkingPoint(BtagWP);
    
    // Book histograms for b-tagging efficiencies
    if(this->makeBtagEfficiencies()) btagScaleFactors_->bookBtagHistograms(fOutput);
    
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

    // ++++ Control Plots ++++



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
    if(this->failsAdditionalJetFlavourSelection(entry)) return kTRUE;

    // Correct for the MadGraph branching fraction being 1/9 for dileptons (PDG average is .108)
    const double weightMadgraphCorrection = this->madgraphWDecayCorrection(entry);

    // Get weight due to pileup reweighting
    const double weightPU = this->weightPileup(entry);

    // Get weight due to generator weights
    const double weightGenerator = this->weightGenerator(entry);

    // Get true level weights
    const double trueLevelWeightNoPileup = weightGenerator*weightMadgraphCorrection;
    const double trueLevelWeight = trueLevelWeightNoPileup*weightPU;

    const tth::GenLevelWeights genLevelWeights(weightMadgraphCorrection, weightPU, weightGenerator,
                                               trueLevelWeightNoPileup, trueLevelWeight);

    // ++++ Control Plots ++++



    //===CUT===
    selectionStep = "1";

    // Check if event was triggered with the same dilepton trigger as the specified analysis channel
    if(this->failsDileptonTrigger(entry)) return kTRUE;



    // === FULL RECO OBJECT SELECTION === (can thus be used at each selection step)

    // Access reco objects, and common generator objects
    const RecoObjects& recoObjects = this->getRecoObjects(entry);
    const CommonGenObjects& commonGenObjects = this->getCommonGenObjects(entry);

    // Create dummies for other objects, non-dummies are created only as soon as needed
    const KinRecoObjects kinRecoObjectsDummy;
    const TopGenObjects topGenObjectsDummy;
    const HiggsGenObjects higgsGenObjectsDummy;

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
    orderIndices(jetIndices, jets, LVpt);
    const int numberOfJets = jetIndices.size();
    const bool has2Jets = numberOfJets > 1 && jets.at(jetIndices.at(0)).pt() >= Lead2JetPtCUT && jets.at(jetIndices.at(1)).pt() >= Lead2JetPtCUT;

    // Fill a vector with all jet pair indices, while sorting each pair by the jet charge:
    // first entry is antiBIndex i.e. with higher jet charge, second entry is bIndex
    //const std::vector<double>& jetChargeGlobalPtWeighted = *recoObjects.jetChargeGlobalPtWeighted_;
    const std::vector<double>& jetChargeRelativePtWeighted = *recoObjects.jetChargeRelativePtWeighted_;
    const tth::IndexPairs& jetIndexPairs = this->chargeOrderedJetPairIndices(jetIndices, jetChargeRelativePtWeighted);

    // Get b-jet indices, apply selection cuts
    // and apply b-tag efficiency MC correction using random number based tag flipping
    // and order b-jets by btag discriminator (beginning with the highest value)
    const std::vector<double>& jetBTagCSV = *recoObjects.jetBTagCSV_;
    const std::vector<int>& jetPartonFlavour = *commonGenObjects.jetPartonFlavour_;
    std::vector<int> bjetIndices = jetIndices;
    selectIndices(bjetIndices, jetBTagCSV, (double)btagScaleFactors_->getWPDiscrValue());
    this->retagJets(bjetIndices, jetIndices, jets, jetPartonFlavour, jetBTagCSV);
    orderIndices(bjetIndices, jetBTagCSV);
    const int numberOfBjets = bjetIndices.size();
    const bool hasBtag = numberOfBjets > 0;

    // Get MET
    const LV& met = *recoObjects.met_;
    const bool hasMetOrEmu = this->channel()=="emu" || met.pt()>MetCUT;

    const tth::RecoObjectIndices recoObjectIndices(allLeptonIndices,
                                                   leptonIndices, antiLeptonIndices,
                                                   leptonIndex, antiLeptonIndex,
                                                   leadingLeptonIndex, nLeadingLeptonIndex,
                                                   leptonXIndex, leptonYIndex,
                                                   jetIndices, jetIndexPairs,
                                                   bjetIndices);

    const tth::GenObjectIndices genObjectIndicesDummy(-1, -1, -1, -1, -1, -1, -1, -1);


    // Determine all reco level weights
    const double weightLeptonSF = this->weightLeptonSF(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId);
    const double weightTriggerSF = this->weightTriggerSF(leptonXIndex, leptonYIndex, allLeptons);
    const double weightNoPileup = trueLevelWeightNoPileup*weightTriggerSF*weightLeptonSF;
    // We do not apply a b-tag scale factor
    //const double weightBtagSF = ReTagJet ? 1. : this->weightBtagSF(jetIndices, jets, jetPartonFlavour);
    constexpr double weightBtagSF = 1.;

    // The weight to be used for filling the histograms
    double weight = weightNoPileup*weightPU;


    tth::RecoLevelWeights recoLevelWeights(weightLeptonSF, weightTriggerSF, weightBtagSF,
                                           weightNoPileup, weight);


    // ++++ Control Plots ++++

    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjectsDummy, higgsGenObjectsDummy,
                  kinRecoObjectsDummy,
                  genObjectIndicesDummy, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  1.);



    //===CUT===
    selectionStep = "2";

    // we need an OS lepton pair matching the trigger selection...
    if (!hasLeptonPair) return kTRUE;

    // ++++ Control Plots ++++

    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjectsDummy, higgsGenObjectsDummy,
                  kinRecoObjectsDummy,
                  genObjectIndicesDummy, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);



    //===CUT===
    selectionStep = "3";

    // ...with at least 20 GeV invariant mass
    if(dilepton.M() < 20.) return kTRUE;

    const bool isZregion = dilepton.M() > 76 && dilepton.M() < 106;
    //const KinRecoObjects& kinRecoObjects = this->getKinRecoObjects(entry);
    //const KinRecoObjects& kinRecoObjects = this->getKinRecoObjectsOnTheFly(leptonIndex, antiLeptonIndex, jetIndices,
    //                                                                       allLeptons, jets, jetBTagCSV, met);
    const KinRecoObjects& kinRecoObjects = kinRecoObjectsDummy;
    //const bool hasSolution = kinRecoObjects.valuesSet_;



    // ++++ Control Plots ++++

    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjectsDummy, higgsGenObjectsDummy,
                  kinRecoObjects,
                  genObjectIndicesDummy, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);



    // ****************************************
    //handle inverted Z cut
    // Z window plots need to be filled here, in order to rescale the contribution to data
    if(isZregion){
        selectionStep = "4zWindow";

        this->fillAll(selectionStep,
                      recoObjects, commonGenObjects,
                      topGenObjectsDummy, higgsGenObjectsDummy,
                      kinRecoObjects,
                      genObjectIndicesDummy, recoObjectIndices,
                      genLevelWeights, recoLevelWeights,
                      weight);

        if(has2Jets){
            selectionStep = "5zWindow";

            this->fillAll(selectionStep,
                          recoObjects, commonGenObjects,
                          topGenObjectsDummy, higgsGenObjectsDummy,
                          kinRecoObjects,
                          genObjectIndicesDummy, recoObjectIndices,
                          genLevelWeights, recoLevelWeights,
                          weight);

            if(hasMetOrEmu){
                selectionStep = "6zWindow";

                this->fillAll(selectionStep,
                              recoObjects, commonGenObjects,
                              topGenObjectsDummy, higgsGenObjectsDummy,
                              kinRecoObjects,
                              genObjectIndicesDummy, recoObjectIndices,
                              genLevelWeights, recoLevelWeights,
                              weight);

                if(hasBtag){
                    selectionStep = "7zWindow";

                    // FIXME: do not use b-tag scale factor
                    //fullWeights *= weightBtagSF;

                    this->fillAll(selectionStep,
                                  recoObjects, commonGenObjects,
                                  topGenObjectsDummy, higgsGenObjectsDummy,
                                  kinRecoObjects,
                                  genObjectIndicesDummy, recoObjectIndices,
                                  genLevelWeights, recoLevelWeights,
                                  weight);
                }
            }
        }
    }



    //=== CUT ===
    selectionStep = "4";

    //Exclude the Z window
    if(this->channel()!="emu" && isZregion) return kTRUE;

    // ++++ Control Plots ++++

    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjectsDummy, higgsGenObjectsDummy,
                  kinRecoObjects,
                  genObjectIndicesDummy, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);



    //=== CUT ===
    selectionStep = "5";

    //Require at least two jets
    if(!has2Jets) return kTRUE;

    // ++++ Control Plots ++++

    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjectsDummy, higgsGenObjectsDummy,
                  kinRecoObjects,
                  genObjectIndicesDummy, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);



    //=== CUT ===
    selectionStep = "6";

    //Require MET > 40 GeV in non-emu channels
    if(!hasMetOrEmu) return kTRUE;

    // ++++ Control Plots ++++

    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjectsDummy, higgsGenObjectsDummy,
                  kinRecoObjects,
                  genObjectIndicesDummy, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);

    // Fill the b-tagging efficiency plots
    if(this->makeBtagEfficiencies()){
        btagScaleFactors_->fillBtagHistograms(jetIndices, jetBTagCSV,
                                              jets, jetPartonFlavour,
                                              weight);
    }

    

    //=== CUT ===
    selectionStep = "7";

    //Require at least one b tagged jet
    if(!hasBtag) return kTRUE;

    weight *= weightBtagSF;



    // === FULL GEN OBJECT SELECTION ===

    // Access top generator object struct, and higgs generator object struct
    const TopGenObjects& topGenObjects = this->getTopGenObjects(entry);
    const HiggsGenObjects& higgsGenObjects = this->getHiggsGenObjects(entry);

    // Do jet matchings for ttbar system
    int genBjetFromTopIndex(-1);
    int genAntiBjetFromTopIndex(-1);
    int matchedBjetFromTopIndex(-1);
    int matchedAntiBjetFromTopIndex(-1);
    if(topGenObjects.valuesSet_){
        const VLV& allGenJets = *commonGenObjects.allGenJets_;
        // Find gen-level b jet and anti-b jet corresponding to (anti)b from (anti)top
        // FIXME: should one clean the genJetCollection to remove low-pt (or high-eta) jets?
        if(this->getGenBjetIndices(genBjetFromTopIndex, genAntiBjetFromTopIndex, topGenObjects, 6)){
            // Match recoJets to the two selected genJets from (anti)top
            this->matchRecoToGenJets(matchedBjetFromTopIndex, matchedAntiBjetFromTopIndex,
                                     jetIndices,
                                     jets,
                                     &allGenJets.at(genBjetFromTopIndex), &allGenJets.at(genAntiBjetFromTopIndex));
        }
    }

    // Do jet matchings for Higgs system
    int genBjetFromHiggsIndex(-1);
    int genAntiBjetFromHiggsIndex(-1);
    int matchedBjetFromHiggsIndex(-1);
    int matchedAntiBjetFromHiggsIndex(-1);
    if(higgsDecayMode == 5){
        const VLV& allGenJets = *commonGenObjects.allGenJets_;
        // Find gen-level b jet and anti-b jet corresponding to (anti)b from Higgs
        // FIXME: should one clean the genJetCollection to remove low-pt (or high-eta) jets?
        if(this->getGenBjetIndices(genBjetFromHiggsIndex, genAntiBjetFromHiggsIndex, topGenObjects, 25)){
            // Match recoJets to the two selected genJets from Higgs
            this->matchRecoToGenJets(matchedBjetFromHiggsIndex, matchedAntiBjetFromHiggsIndex,
                                     jetIndices,
                                     jets,
                                     &allGenJets.at(genBjetFromHiggsIndex), &allGenJets.at(genAntiBjetFromHiggsIndex));
        }
    }

    const tth::GenObjectIndices genObjectIndices(genBjetFromTopIndex, genAntiBjetFromTopIndex,
                                                 matchedBjetFromTopIndex, matchedAntiBjetFromTopIndex,
                                                 genBjetFromHiggsIndex, genAntiBjetFromHiggsIndex,
                                                 matchedBjetFromHiggsIndex, matchedAntiBjetFromHiggsIndex);



    // ++++ Control Plots ++++

    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjects, higgsGenObjects,
                  kinRecoObjects,
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



bool HiggsAnalysis::getGenBjetIndices(int& genBjetIndex, int& genAntiBjetIndex,
                                      const TopGenObjects& topGenObjects, const int pdgId)
{
    if(!pdgId>0){
        std::cerr<<"ERROR! Method getGenBJetIndices needs a pdgId>0, but used is: "<<pdgId
                 <<"\n...break\n\n";
        exit(71);
    }

    const std::vector<int>& genBHadFlavour(*topGenObjects.genBHadFlavour_);
    const std::vector<int>& genBHadJetIndex(*topGenObjects.genBHadJetIndex_);

    for(size_t iBHadron=0; iBHadron<genBHadFlavour.size(); ++iBHadron){
        const int flavour = genBHadFlavour.at(iBHadron);
        if(std::abs(flavour) != std::abs(pdgId)) continue;     // Skipping hadrons with the wrong flavour
        // Assigning jet index of corresponding hadron. Set to -2 if >1 hadrons found for the same flavour
        if(flavour>0) genBjetIndex = (genBjetIndex==-1) ? genBHadJetIndex.at(iBHadron) : -2;
        else if(flavour<0) genAntiBjetIndex = (genAntiBjetIndex==-1) ? genBHadJetIndex.at(iBHadron) : -2;
    }

    // If no unique match of jets from (anti)b from (anti)top is found, return false
    if(genBjetIndex<0 || genAntiBjetIndex<0 || genBjetIndex==genAntiBjetIndex){
        return false;
    }
    return true;
}



bool HiggsAnalysis::matchRecoToGenJets(int& matchedBjetIndex, int& matchedAntiBjetIndex,
                                       const std::vector<int>& jetIndices,
                                       const VLV& jets,
                                       const LV* genBjet, const LV* genAntiBjet)
{
    using ROOT::Math::VectorUtil::DeltaR;
    
    // Find closest jet and its distance in deltaR
    double deltaRBjet(999.);
    double deltaRAntiBjet(999.);
    for(const auto& index : jetIndices){
        double deltaR = DeltaR(*genBjet, jets.at(index));
        if(deltaR < deltaRBjet){
            deltaRBjet = deltaR;
            matchedBjetIndex = index;
        }
        deltaR = DeltaR(*genAntiBjet, jets.at(index));
        if(deltaR < deltaRAntiBjet){
            deltaRAntiBjet = deltaR;
            matchedAntiBjetIndex = index;
        }
    }
    
    // Call a jet matched if it is close enough, and has similar pt
    if(deltaRBjet>0.4){
        matchedBjetIndex = -2;
    }
    else if(matchedBjetIndex >= 0){
        const double ptRecoJet = jets.at(matchedBjetIndex).pt();
        const double ptBjet = genBjet->pt();
        const double deltaPtRel = (ptBjet - ptRecoJet)/ptBjet;
        if(deltaPtRel<-0.4 || deltaPtRel>0.6) matchedBjetIndex = -3;
    }
    
    if(deltaRAntiBjet>0.4){
        matchedAntiBjetIndex = -2;
    }
    else if(matchedAntiBjetIndex >= 0){
        const double ptRecoJet = jets.at(matchedAntiBjetIndex).pt();
        const double ptAntiBjet = genAntiBjet->pt();
        const double deltaPtRel = (ptAntiBjet - ptRecoJet)/ptAntiBjet;
        if(deltaPtRel<-0.4 || deltaPtRel>0.6) matchedAntiBjetIndex = -3;
    }
    
    // Check if both gen jets are successfully matched to different reco jets
    if(matchedBjetIndex<0 || matchedAntiBjetIndex<0 || matchedBjetIndex==matchedAntiBjetIndex) return false;
    
    return true;
}



void HiggsAnalysis::SetInclusiveHiggsDecayMode(const int inclusiveHiggsDecayMode)
{
    inclusiveHiggsDecayMode_ = inclusiveHiggsDecayMode;
}



void HiggsAnalysis::SetAdditionalBjetMode(const int additionalBjetMode)
{
    additionalBjetMode_ = additionalBjetMode;
}




void HiggsAnalysis::SetAllAnalyzers(std::vector<AnalyzerBaseClass*> v_analyzer)
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



bool HiggsAnalysis::failsAdditionalJetFlavourSelection(const Long64_t& entry)const
{
    if(additionalBjetMode_ == -999) return false;
    
    // Use the full sample for creating btag efficiencies
    if(this->makeBtagEfficiencies()) return false;
    
    // FIXME: this is a workaround as long as there is no specific additional jet flavour info written to nTuple
    const TopGenObjects& topGenObjects = this->getTopGenObjects(entry);
    const CommonGenObjects& commonGenObjects = this->getCommonGenObjects(entry);
    
    std::vector<int> genAddBJetIdNotFromTop;

    float signalJetPt_min = 20.;
    float signalJetEta_max = 2.5;

    for(size_t iHad=0; iHad<topGenObjects.genBHadJetIndex_->size(); iHad++) {
//         printf("hadId: %d\n", iHad);
        if(topGenObjects.genBHadFromTopWeakDecay_->at(iHad)==0) {
            int genJetId = topGenObjects.genBHadJetIndex_->at(iHad);
//             printf(" genJetId: %d\n", genJetId);
            if(genJetId>=0) {
                if(commonGenObjects.allGenJets_->at(genJetId).Pt()>signalJetPt_min && std::fabs(commonGenObjects.allGenJets_->at(genJetId).Eta())<signalJetEta_max) {
                    if(std::find(genAddBJetIdNotFromTop.begin(), genAddBJetIdNotFromTop.end(), genJetId) == genAddBJetIdNotFromTop.end()) {
                        genAddBJetIdNotFromTop.push_back(genJetId);
                    }
                }   // If the jet suits the signal selection criterea
            }   // If the hadron is clustered to any jet
        }   // If the hadron is additional to b-hadrons from Top
    }   // End of loop over all b-hadrons

    const unsigned int nExtraBjets = genAddBJetIdNotFromTop.size();
//     printf("Mode: %d  nAddJets: %d\n", additionalBJetMode_, nExtraBjets);

    if(additionalBjetMode_==2 && nExtraBjets>=2) return false;
    if(additionalBjetMode_==1 && nExtraBjets==1) return false;
    if(additionalBjetMode_==0 && nExtraBjets==0) return false;

    return true;
}



void HiggsAnalysis::fillAll(const std::string& selectionStep,
                            const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                            const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                            const KinRecoObjects& kinRecoObjects,
                            const tth::GenObjectIndices& genObjectIndices, const tth::RecoObjectIndices& recoObjectIndices,
                            const tth::GenLevelWeights& genLevelWeights, const tth::RecoLevelWeights& recoLevelWeights,
                            const double& defaultWeight)const
{
    for(AnalyzerBaseClass* analyzer : v_analyzer_){
        if(analyzer) analyzer->fill(recoObjects, commonGenObjects,
                                    topGenObjects, higgsGenObjects,
                                    kinRecoObjects,
                                    recoObjectIndices, genObjectIndices,
                                    genLevelWeights, recoLevelWeights,
                                    defaultWeight, selectionStep);
    }
    
    for(MvaTreeHandlerBase* mvaTreeHandler : v_mvaTreeHandler_){
        if(mvaTreeHandler) mvaTreeHandler->fill(recoObjects, commonGenObjects,
                                                topGenObjects, higgsGenObjects,
                                                kinRecoObjects,
                                                recoObjectIndices, genObjectIndices,
                                                genLevelWeights, recoLevelWeights,
                                                defaultWeight, selectionStep);
    }
}



void HiggsAnalysis::bookAll()
{
    for(AnalyzerBaseClass* analyzer : v_analyzer_){
        if(analyzer) analyzer->book(fOutput);
    }
}



void HiggsAnalysis::clearAll()
{
    for(AnalyzerBaseClass* analyzer : v_analyzer_){
        if(analyzer) analyzer->clear();
    }
}








