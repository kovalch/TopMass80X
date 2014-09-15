#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <exception>

#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TString.h>
#include <Math/VectorUtil.h>
#include <TSelectorList.h>
#include <Rtypes.h>

#include "AnalyzerDijet.h"
#include "analysisStructs.h"
#include "higgsUtils.h"
#include "MvaReaderTopJets.h"
#include "MvaVariablesTopJets.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"
#include "../../common/include/KinematicReconstructionSolution.h"





AnalyzerDijet::AnalyzerDijet(const char* mva2dWeightsFile, const std::string& corName, const std::string& swpName,
                             const std::vector<TString>& selectionStepsNoCategories,
                             const std::vector<TString>& stepsForCategories,
                             const JetCategories* jetCategories, bool doHadronMatchingComparison, bool doLeadingJetsAnalysis):
AnalyzerBase("dijet_", selectionStepsNoCategories, stepsForCategories, jetCategories),
weightsCorrect_(0),
weightsSwapped_(0),
doHadronMatchingComparison_(doHadronMatchingComparison),
doLeadingJetsAnalysis_(doLeadingJetsAnalysis)
{
    std::cout<<"--- Beginning setting up dijet analyzer\n";
    
    // Setting up the MVA weights if available
    std::string mvaWeightsFolder(mva2dWeightsFile);
    mvaWeightsFolder.erase(mvaWeightsFolder.rfind('/'));
    if(!gSystem->OpenDirectory(mvaWeightsFolder.c_str()) && corName.length()>0 && swpName.length()>0) {
        throw std::logic_error("AnalyzerDijet::AnalyzerDijet     WARNING! Folder with MVA weights doesn't exist\n");
    }
    
    std::string tempStr;
    // Setting up the correct training
    if(corName.length()>0) {
        tempStr = mvaWeightsFolder;
        tempStr.append("/").append(corName).append(".weights.xml");
        weightsCorrect_ = new MvaReaderTopJets("BDT method");
        weightsCorrect_->book(tempStr);
    }
    // Setting up the swapped training
    if(swpName.length()>0) {
        tempStr = mvaWeightsFolder;
        tempStr.append("/").append(swpName).append(".weights.xml");
        weightsSwapped_ = new MvaReaderTopJets("BDT method");
        weightsSwapped_->book(tempStr);
    }
    

    std::cout<<"=== Finishing setting up dijet analyzer\n\n";
}



void AnalyzerDijet::fillHistos(const EventMetadata& eventMetadata,
                               const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                               const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                               const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                               const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                               const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                               const double& weight, const TString&,
                               std::map<TString, TH1*>& m_histogram)
{
    // Extracting input data to more comfortable variables
    const VLV& allJets = (recoObjects.valuesSet_) ? *recoObjects.jets_ : VLV();
    const std::vector<int>& jetsId = recoObjectIndices.jetIndices_;           // Selected jets (point to jets from allJets)
    const std::vector<int>& bJetsId = recoObjectIndices.bjetIndices_;         // B-tagged jets (point to jets from allJets)
//     printf("NJETS: %d  NBJETS: %d ###############################\n", (int)jetsId.size(), (int)bJetsId.size());
    std::vector<int> topJetsId;                                               // Jets from ttbar by KinReco (Point to jets from allJets)
    if(kinematicReconstructionSolutions.numberOfSolutions()) {
        topJetsId.push_back(kinematicReconstructionSolutions.solution().bjetIndex());
        topJetsId.push_back(kinematicReconstructionSolutions.solution().antiBjetIndex());
    }
    const std::vector<double>& allJetsBtagDiscriminant = (recoObjects.valuesSet_) ? *recoObjects.jetBTagCSV_ : std::vector<double>(0);
    

    // Setting variables of gen. level if available
    const VLV& genAllJets = (topGenObjects.valuesSet_) ? *topGenObjects.allGenJets_ : VLV();
    const std::vector<int>& jetPartonFlavour = (commonGenObjects.valuesSet_) ? *commonGenObjects.jetPartonFlavour_ : std::vector<int>(0);
    std::vector<int> genAllJetsId = common::initialiseIndices(genAllJets);
    std::vector<int> genJetsId = genObjectIndices.genJetIndices_;
    std::vector<int> genBJetsId = genObjectIndices.genBjetIndices_;
    const std::vector<int>& bHadJetIndex = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadJetIndex_ : std::vector<int>(0);
    const std::vector<int>& bHadFlavour = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadFlavour_ : std::vector<int>(0);
    const std::vector<int>& bHadIndex = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadIndex_ : std::vector<int>(0);
    const std::vector<LV>&  bHadPlusMothers = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadPlusMothers_ : std::vector<LV>(0);
    const std::vector<int>& bHadFromTopWeakDecay = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadFromTopWeakDecay_ : std::vector<int>(0);
    std::vector<int> trueTopAllJetsId;          // Reco jets coming from top (Point to jets from allJets) [Can be any jet]
    std::vector<int> trueTopJetsId;             // Reco jets coming from top (Point to jets from allJets) [Only selected jets]
    std::vector<int> trueHiggsAllJetsId;        // Reco jets coming from higgs (Point to jets from allJets) [Can be any jet]
    std::vector<int> trueHiggsJetsId;           // Reco jets coming from higgs (Point to jets from allJets) [Only selected jets]
    std::vector<int> trueTopBJetsId;            // Reco jets coming from top (Point to b-tagged jets from allJets) [Only selected jets]
    std::vector<int> trueHiggsBJetsId;          // Reco jets coming from higgs (Point to b-tagged jets from allJets) [Only selected jets]
    std::vector<int> topBJetsId;                // Reco jets coming from top according to KinReco (Point to b-tagged jets from allJets) [Only selected jets]
    std::vector<int> genTopAllJetsId;           // Gen jets coming from top
    std::vector<int> genHiggsAllJetsId;         // Gen jets coming from higgs
    std::vector<int> genTopJetsId;              // Gen jets coming from top
    std::vector<int> genHiggsJetsId;            // Gen jets coming from higgs

    // Checking how we loose additional b-jets from tt+bb events causing tt+b, tt+other events
    checkAdditionalGenBJetAcceptance(topGenObjects, genObjectIndices, m_histogram, weight);
    
    // Checking how many additional b-hadrons there are in b-jets
    std::vector<std::vector<int> > genJetBhadronIndices;
    if(topGenObjects.valuesSet_) {
        genJetBhadronIndices = this->matchBhadronsToGenJets(genAllJets, topGenObjects);
        for(unsigned int jetId : genBJetsId) {
            std::vector<int> hadIds = genJetBhadronIndices.at(jetId);
            int nHads_add = 0;
            int nHads_top = 0;
            for(unsigned int hadId : hadIds) {
                if(std::abs(bHadFlavour.at(hadId)) == 6) nHads_top++;
                if(std::abs(bHadFromTopWeakDecay.at(hadId)) == 0) nHads_add++;
            }
            if(nHads_top>0) m_histogram["nBHadAddInGenBJet_top"]->Fill(nHads_add, weight);
            else if(nHads_add>0) m_histogram["nBHadAddInGenBJet_add"]->Fill(nHads_add, weight);
        }
    }

    // Creating a list of LorentzVectors of the bHadrons
    VLV bHadLVs;
    for(const size_t index : bHadIndex) {
        LV bHadLV = bHadPlusMothers.at(index);
        bHadLVs.push_back(bHadLV);
    }

    // Filling vectors of true reco top jets
    const int trueTopBJetId = genObjectIndices.recoBjetFromTopIndex_;
    
    if(trueTopBJetId >= 0) {
        trueTopAllJetsId.push_back(trueTopBJetId);
        if(isInVector(jetsId, trueTopBJetId)) trueTopJetsId.push_back(trueTopBJetId);
        if(isInVector(bJetsId, trueTopBJetId)) trueTopBJetsId.push_back(trueTopBJetId);
    }
    const int trueTopAntiBJetId = genObjectIndices.recoAntiBjetFromTopIndex_;
    if(trueTopAntiBJetId >= 0) {
        trueTopAllJetsId.push_back(trueTopAntiBJetId);
        if(isInVector(jetsId, trueTopAntiBJetId)) trueTopJetsId.push_back(trueTopAntiBJetId);
        if(isInVector(bJetsId, trueTopAntiBJetId)) trueTopBJetsId.push_back(trueTopAntiBJetId);
    }
    
//     printf("trueTopJetId: %d | %d  jets: ", trueTopBJetId, trueTopAntiBJetId, (int)jetsId.at(jetsId.size()-1), (int)bJetsId.at(bJetsId.size()-1));
//     for(int i: jetsId) printf("%d ", i);
//     printf(" bjets: ");
//     for(int i: bJetsId) printf("%d ", i);
//     printf("\n");

    // Filling vectors of true reco Higgs jets
    const int trueHiggsBJetId = genObjectIndices.recoBjetFromHiggsIndex_;
    if(trueHiggsBJetId >= 0) {
        trueHiggsAllJetsId.push_back(trueHiggsBJetId);
        if(isInVector(jetsId, trueHiggsBJetId)) trueHiggsJetsId.push_back(trueHiggsBJetId);
        if(isInVector(bJetsId, trueHiggsBJetId)) trueHiggsBJetsId.push_back(trueHiggsBJetId);
    }
    const int trueHiggsAntiBJetId = genObjectIndices.recoAntiBjetFromHiggsIndex_;
    if(trueHiggsAntiBJetId >= 0) {
        trueHiggsAllJetsId.push_back(trueHiggsAntiBJetId);
        if(isInVector(jetsId, trueHiggsAntiBJetId)) trueHiggsJetsId.push_back(trueHiggsAntiBJetId);
        if(isInVector(bJetsId, trueHiggsAntiBJetId)) trueHiggsBJetsId.push_back(trueHiggsAntiBJetId);
    }
    // Getting indices of true gen Higgs jets
//     const int trueHiggsGenBJetId = genObjectIndices.genBjetFromHiggsIndex_;
//     const int trueHiggsGenAntiBJetId = genObjectIndices.genAntiBjetFromHiggsIndex_;

    // Filling the dijet mass of true reco jets from Higgs
    if(trueHiggsBJetId >= 0 && trueHiggsAntiBJetId >= 0) {
        LV dijet_trueH = allJets.at(trueHiggsBJetId) + allJets.at(trueHiggsAntiBJetId);
        m_histogram["dijet_mass_trueHiggsJets"]->Fill(dijet_trueH.M(), weight);
    }

    for(size_t hadId =0; hadId<bHadFlavour.size(); ++hadId) {
        int flavour = std::abs(bHadFlavour.at(hadId));
        if(flavour!=6 && flavour!=25) continue;
        int jetId = bHadJetIndex.at(hadId);
        if(jetId<0) continue;
        if(flavour==6) putUniquelyInVector(genTopAllJetsId, jetId);
        else if(flavour==25) putUniquelyInVector(genHiggsAllJetsId, jetId);
        if(!isInVector(genBJetsId, jetId)) continue;
        if(flavour==6) putUniquelyInVector(genTopJetsId, jetId);
        else if(flavour==25) putUniquelyInVector(genHiggsJetsId, jetId);
    }
    
    if(trueTopAllJetsId.size() == 2 && genTopAllJetsId.size()==2) {
        for(int jetId : trueTopAllJetsId) {
            m_histogram["topBJet_Pt_true"]->Fill(allJets.at(jetId).Pt(), weight);
            m_histogram["topBJet_Btag_true"]->Fill(allJetsBtagDiscriminant.at(jetId), weight);
        }
        for(int jetId : genTopAllJetsId) {
            m_histogram["topBJet_Pt_gen"]->Fill(genAllJets.at(jetId).Pt(), weight);
        }
    }
    if(trueHiggsAllJetsId.size() == 2 && genHiggsAllJetsId.size()==2) {
        for(int jetId : trueHiggsAllJetsId) {
            m_histogram["higgsBJet_Pt_true"]->Fill(allJets.at(jetId).Pt(), weight);
            m_histogram["higgsBJet_Btag_true"]->Fill(allJetsBtagDiscriminant.at(jetId), weight);
        }
        for(int jetId : genHiggsAllJetsId) {
            m_histogram["higgsBJet_Pt_gen"]->Fill(genAllJets.at(jetId).Pt(), weight);
        }
    }
    
    m_histogram["topJet_multiplicity_allGen"]->Fill(genTopAllJetsId.size(), weight);
    m_histogram["higgsJet_multiplicity_allGen"]->Fill(genHiggsAllJetsId.size(), weight);
    m_histogram["topJet_multiplicity_gen"]->Fill(genTopJetsId.size(), weight);
    m_histogram["higgsJet_multiplicity_gen"]->Fill(genHiggsJetsId.size(), weight);
    
    
    // Filling the dijet mass of true gen jets from top
    if(genTopJetsId.size()==2) {
        LV dijet_genH = genAllJets.at(genTopJetsId.at(0)) + genAllJets.at(genTopJetsId.at(0));
        m_histogram["dijet_mass_genTopJets"]->Fill(dijet_genH.M(), weight);
    }
    // Filling the dijet mass of true gen jets from Higgs
    if(genHiggsJetsId.size() == 2) {
        LV dijet_genH = genAllJets.at(genHiggsJetsId.at(0)) + genAllJets.at(genHiggsJetsId.at(1));
        m_histogram["dijet_mass_genHiggsJets"]->Fill(dijet_genH.M(), weight);
    }


    int nTopJetsInAllTrue=0;
    int nTopJetsInTrue=0;
    // Selecting b-tagged jets among top jets from KinReco
    for(size_t iJet=0, iJetN=topJetsId.size(); iJet<iJetN; iJet++) {
        int iAllJet = topJetsId.at(iJet);
        if(isInVector(trueTopAllJetsId,iAllJet)) nTopJetsInAllTrue++;
        if(isInVector(trueTopJetsId,iAllJet)) nTopJetsInTrue++;
        if(!isInVector(bJetsId,iAllJet)) continue;
        topBJetsId.push_back(iAllJet);
    }


    // Identifying numbers of jets and b-jets
    unsigned int nAllJets = allJets.size();
    unsigned int nJets = jetsId.size();
    unsigned int nBJets = bJetsId.size();
    unsigned int nGenJets = genAllJets.size();
    

    // Analyzing selected jets
    int nJetsPtLt30 = 0;
    int nBJetsPtLt30 = 0;
    float bJet_dRmin = 999.9;
    for(size_t iJet = 0, iJetN=jetsId.size(); iJet<iJetN; iJet++) {
        unsigned int iAllJet = jetsId.at(iJet);
        LV jet = allJets.at(iAllJet);

        // Filling information about the jets
        m_histogram["jet_pt"]->Fill(jet.Pt(), weight);
        m_histogram["jet_eta"]->Fill(jet.Eta(), weight);
        m_histogram["jet_btagDiscriminator"]->Fill(allJetsBtagDiscriminant.at(iAllJet), weight);
        if(jetPartonFlavour.size()>0) m_histogram["jet_partonFlavour"]->Fill(jetPartonFlavour.at(iAllJet), weight);

        // Filling information about b-tagged jets
        if(isInVector(bJetsId, iAllJet)) {
            m_histogram["bJet_pt"]->Fill(jet.Pt(), weight);
            m_histogram["bJet_eta"]->Fill(jet.Eta(), weight);
            m_histogram["bJet_btagDiscriminator"]->Fill(allJetsBtagDiscriminant.at(iAllJet), weight);
            if(jetPartonFlavour.size()>0) m_histogram["bJet_partonFlavour"]->Fill(jetPartonFlavour.at(iAllJet), weight);

            // Checking dR between closest pair of reco b-jets
            for(size_t iJet2 = 0; iJet2<iJet; iJet2++) {
                unsigned int iAllJet2 = jetsId.at(iJet2);
                if(!isInVector(bJetsId, iAllJet2)) continue;
                LV jet2 = allJets.at(iAllJet2);
                float dR = ROOT::Math::VectorUtil::DeltaR(jet, jet2);
                if(dR<bJet_dRmin) bJet_dRmin = dR;
            }       // End of loop over all other reco b-jets
        }       // End of loop over all reco b-jets


        ////////////////////////////////////////////////////////////////// Analyzing only jets with Pt<30
        if(jet.Pt()>=30) continue;
        nJetsPtLt30++;
        m_histogram["jet_PtLt30_btagDiscriminator"]->Fill(allJetsBtagDiscriminant.at(iAllJet), weight);
        if(isInVector(bJetsId, iAllJet)) {
            nBJetsPtLt30++;
            m_histogram["bJet_PtLt30_btagDiscriminator"]->Fill(allJetsBtagDiscriminant.at(iAllJet), weight);
        }

        // Finding the corresponding genJetId and its flavours
        int genJetId = genJetIdOfRecoJet(jet, genAllJets);
        std::vector<int> flavJet = bHadFlavoursInGenJet(genJetId, bHadJetIndex, bHadFlavour, true);
        if(!isInVector(bJetsId,iAllJet)) continue;            // Skip jet if it is not b-tagged
        for(size_t iFlav = 0, iFlavN = flavJet.size(); iFlav<iFlavN; iFlav++) {
            m_histogram["jet_PtLt30_flavour"]->Fill(flavJet.at(iFlav), weight);
            if(!isInVector(bJetsId, iAllJet)) continue;          // Skip jet if it was assigned to the Top
            m_histogram["bJet_PtLt30_flavour"]->Fill(flavJet.at(iFlav), weight);
        }       // End of loop over flavours of the jet

    }       // End of loop over all reco jets in the event
    m_histogram["bJet_dRmin"]->Fill(bJet_dRmin, weight);


    m_histogram["allJet_multiplicity"]->Fill(nAllJets, weight);
    m_histogram["jet_multiplicity"]->Fill(nJets, weight);
    m_histogram["bJet_multiplicity"]->Fill(nBJets, weight);
    m_histogram["genJet_multiplicity"]->Fill(nGenJets, weight);
    m_histogram["topJet_multiplicity_allTrue"]->Fill(trueTopAllJetsId.size(), weight);
    m_histogram["topJet_multiplicity_true"]->Fill(trueTopJetsId.size(), weight);
    m_histogram["topJet_multiplicity_reco"]->Fill(topJetsId.size(), weight);
    m_histogram["topJet_multiplicity_recoInTrue"]->Fill(nTopJetsInTrue, weight);
    m_histogram["topJet_multiplicity_recoInAllTrue"]->Fill(nTopJetsInAllTrue, weight);
    m_histogram["topBJet_multiplicity_true"]->Fill(trueTopBJetsId.size(), weight);
    m_histogram["topBJet_multiplicity_reco"]->Fill(topBJetsId.size(), weight);
    m_histogram["higgsJet_multiplicity_allTrue"]->Fill(trueHiggsAllJetsId.size(), weight);
    m_histogram["higgsJet_multiplicity_true"]->Fill(trueHiggsJetsId.size(), weight);
    m_histogram["higgsBJet_multiplicity_true"]->Fill(trueHiggsBJetsId.size(), weight);
    m_histogram["jet_PtLt30_multiplicity"]->Fill(nJetsPtLt30, weight);
    m_histogram["bJet_PtLt30_multiplicity"]->Fill(nBJetsPtLt30, weight);
    m_histogram["weight"]->Fill(weight, weight);
    
    // Counting multipllicity of jets with different Pt
    int nBJets_pt[3] = {0,0,0};
    for(int iJet : bJetsId) {
        float pt = allJets.at(iJet).Pt();
        if(pt<30.f) nBJets_pt[0]++;
        else if(pt<70) nBJets_pt[1]++;
        else if(pt<500) nBJets_pt[2]++;
    }
    m_histogram["bJet_multiplicity_pt_0_30"]->Fill(nBJets_pt[0], weight);
    m_histogram["bJet_multiplicity_pt_30_70"]->Fill(nBJets_pt[1], weight);
    m_histogram["bJet_multiplicity_pt_70_500"]->Fill(nBJets_pt[2], weight);

    std::vector<int> emptyVector;

    /////////////////////////////////////////////////////////////////// FILLING DIJET MASS FOR DIFFERENT JET SELECTIONS
    const bool fillAllCombinations = true;
    ///////////////////////////////////////////////////////////////////////////////////////////////// ALL SELECTED JETS
    // Analyzing jet pairs for all jet combinations
//     float correctPairFraction_all = correctPairFraction(allJets, jetsId, bJetsId, allJetsBtagDiscriminant, emptyVector, trueHiggsJetsId, weight, m_histogram, "all", fillAllCombinations);
//     m_histogram["dijet_correctPairFraction_all"]->Fill(correctPairFraction_all, weight);

    // Analyzing jet pairs for all combinations of jets tagged with CSVM
    std::vector<int> bMJetsId;
    for(size_t jetId : bJetsId) { if(allJetsBtagDiscriminant.at(jetId)>=0.679) bMJetsId.push_back(jetId); };
    float correctPairFraction_allM = correctPairFraction(allJets, jetsId, bMJetsId, allJetsBtagDiscriminant, emptyVector, trueHiggsJetsId, weight, m_histogram, "allM", fillAllCombinations);
    m_histogram["dijet_correctPairFraction_allM"]->Fill(correctPairFraction_allM, weight);

    // Analyzing jet pairs for all jet combinations except true b-jets from top
    float correctPairFraction_trueTopJets = correctPairFraction(allJets, jetsId, bJetsId, allJetsBtagDiscriminant, trueTopJetsId, trueHiggsJetsId, weight, m_histogram, "trueTopJets", fillAllCombinations);
    m_histogram["dijet_correctPairFraction_trueTopJets"]->Fill(correctPairFraction_trueTopJets, weight);

    // Analyzing jet pairs for all jet combinations except reco b-jets from top found by kinematic reconstruction
    float correctPairFraction_recoTopJets = correctPairFraction(allJets, jetsId, bJetsId, allJetsBtagDiscriminant, topJetsId, trueHiggsJetsId, weight, m_histogram, "recoTopJets", fillAllCombinations);
    m_histogram["dijet_correctPairFraction_recoTopJets"]->Fill(correctPairFraction_recoTopJets, weight);
    

    ///////////////////////////////////////////////////////////////////////////////////////////////// Checking min dR and flavours of Jets, Hadrons
    // Counting number of gen b-jets
    int nGenBjets = 0;
    int nGenBhads = bHadJetIndex.size();
    int nGenBhadsNoG = 0;
    int nGenBhadsNoT = 0;
    int nGenBhadsT = 0;
    int nGenBhadsH = 0;
    float genBJet_dRmin = 999.9;
    int nGenAddBhadsFromTop = 0;
    int nGenAddBhadsNotFromTop = 0;
    std::vector<int> genAddBJetIdFromTop;
    std::vector<int> genAddBJetIdNotFromTop;
    
    float signalJetPt_min = 20.;
    float signalJetEta_max = 2.5;
    std::vector<int> genAddBJetIdFromTop_signal;
    std::vector<int> genAddBJetIdNotFromTop_signal;
    // Finding jets that are matched to b-hadrons
    for(size_t iJet=0, nJets=genAllJets.size(); iJet<nJets; iJet++) {
        if(!isInVector(bHadJetIndex, iJet)) continue;       // Skipping if jet is not a true b-jet
        nGenBjets++;
        m_histogram["genBJet_pt"]->Fill(genAllJets.at(iJet).Pt(), weight);
        m_histogram["genBJet_eta"]->Fill(genAllJets.at(iJet).Eta(), weight);

        // Checking dR between the closest pair of gen b-jets
        for(size_t iJet2=0; iJet2<iJet; iJet2++) {
            if(!isInVector(bHadJetIndex, iJet2)) continue;
            float dR = ROOT::Math::VectorUtil::DeltaR(genAllJets.at(iJet), genAllJets.at(iJet2));
            if(dR<genBJet_dRmin) genBJet_dRmin = dR;
        }       // End of loop over all other gen b-jets

        float genBJet_recoJet_dRmin = 999.9;
        float genBJet_recoBJet_dRmin = 999.9;
        // Checking dR between each gen b-jet and closest reco jet
        for(size_t iJet2 = 0; iJet2<jetsId.size(); iJet2++) {
            unsigned int iAllJet2 = jetsId.at(iJet2);
            LV jet2 = allJets.at(iAllJet2);
            float dR = ROOT::Math::VectorUtil::DeltaR(genAllJets.at(iJet), jet2);
            if(dR<genBJet_recoJet_dRmin) genBJet_recoJet_dRmin = dR;
            if(!isInVector(bJetsId, iAllJet2)) continue;
            float dR_b = ROOT::Math::VectorUtil::DeltaR(genAllJets.at(iJet), jet2);
            if(dR_b<genBJet_recoBJet_dRmin) genBJet_recoBJet_dRmin = dR_b;
        }       // End of loop over all reco b-jets
        m_histogram["genBJet_recoJet_dRmin"]->Fill(genBJet_recoJet_dRmin, weight);
        m_histogram["genBJet_recoBJet_dRmin"]->Fill(genBJet_recoBJet_dRmin, weight);

    }       // End of loop over all gen b-jets
    m_histogram["genBJet_dRmin"]->Fill(genBJet_dRmin, weight);

    float bHad_dRmin = 999.9;
    // Getting number of hadrons of different flavours
    for(size_t iHad=0, nHads=bHadFlavour.size(); iHad<nHads; iHad++) {
        if(std::abs(bHadFlavour.at(iHad))!=21) nGenBhadsNoG++;
        if(std::abs(bHadFlavour.at(iHad))!=6) nGenBhadsNoT++;
        if(std::abs(bHadFlavour.at(iHad))==6) nGenBhadsT++;
        if(std::abs(bHadFlavour.at(iHad))==25) nGenBhadsH++;
        m_histogram["bHad_flavour"]->Fill(bHadFlavour.at(iHad), weight);
        
        if(bHadFromTopWeakDecay.size()>iHad && bHadFromTopWeakDecay.at(iHad)==1 && std::abs(bHadFlavour.at(iHad))!=6) {
            nGenAddBhadsFromTop++;
            int genJetId = bHadJetIndex.at(iHad);
            if(genJetId>=0 && !isInVector(genTopAllJetsId,genJetId)) {
                putUniquelyInVector(genAddBJetIdFromTop, genJetId);
                if(genAllJets.at(genJetId).Pt()>signalJetPt_min && std::fabs(genAllJets.at(genJetId).Eta())<signalJetEta_max && !isInVector(genTopJetsId, genJetId)) {
                    putUniquelyInVector(genAddBJetIdFromTop_signal, genJetId);
                }
            }
        }
        if(bHadFromTopWeakDecay.size()>iHad && bHadFromTopWeakDecay.at(iHad)==0) {
            nGenAddBhadsNotFromTop++;
            int genJetId = bHadJetIndex.at(iHad);
            if(genJetId>=0 && !isInVector(genTopAllJetsId,genJetId)) {
                putUniquelyInVector(genAddBJetIdNotFromTop, genJetId);
                if(genAllJets.at(genJetId).Pt()>signalJetPt_min && std::fabs(genAllJets.at(genJetId).Eta())<signalJetEta_max && !isInVector(genTopJetsId, genJetId)) {
                    putUniquelyInVector(genAddBJetIdNotFromTop_signal, genJetId);
                }
            }
        }
        

        LV had1 = bHadLVs.at(iHad);
        for(size_t iHad2=0; iHad2<iHad; iHad2++) {
            LV had2 = bHadLVs.at(iHad2);
            float dR = ROOT::Math::VectorUtil::DeltaR(had1, had2);
            if(dR<bHad_dRmin) bHad_dRmin = dR;
        }       // End of loop over all remaining hadrons in event
    }       // End of loop over all hadrons in event

    m_histogram["bHad_dRmin"]->Fill(bHad_dRmin, weight);

    m_histogram["genBJet_multiplicity"]->Fill(nGenBjets, weight);
    m_histogram["bHad_multiplicity"]->Fill(nGenBhads, weight);
    m_histogram["bHad_top_multiplicity"]->Fill(nGenBhadsT, weight);
    m_histogram["bHad_higgs_multiplicity"]->Fill(nGenBhadsH, weight);
    m_histogram["bHadNoG_multiplicity"]->Fill(nGenBhadsNoG, weight);
    ((TH2*)m_histogram["bHadVsBJet_multiplicity"])->Fill(nGenBhads, nBJets, weight);
    ((TH2*)m_histogram["bHadVsBHadNoT_multiplicity"])->Fill(nGenBhads, nGenBhadsNoT, weight);

    m_histogram["nHadAddFromTopWeakDecay"]->Fill(nGenAddBhadsFromTop);
    m_histogram["nHadAddNotFromTopWeakDecay"]->Fill(nGenAddBhadsNotFromTop);
    m_histogram["nGenJetAddFromTopWeakDecay"]->Fill( (int)genAddBJetIdFromTop.size() );
    m_histogram["nGenJetAddFromTopWeakDecay_signal"]->Fill( (int)genAddBJetIdFromTop_signal.size() );
    m_histogram["nGenJetAddNotFromTopWeakDecay"]->Fill( (int)genAddBJetIdNotFromTop.size() );
    m_histogram["nGenJetAddNotFromTopWeakDecay_signal"]->Fill( (int)genAddBJetIdNotFromTop_signal.size() );
    

    ///////////////////////////////////////////////////////////////////////// COMPARISON OF GEN MATCHING AND dR MATCHING
    if(doHadronMatchingComparison_) {
        fillGenRecoMatchingComparisonHistos(topGenObjects, higgsGenObjects, bHadLVs, bHadFlavour, bHadJetIndex, genAllJets, m_histogram, weight);
    }

    ///////////////////////////////////////////////////////////////////////// B-TAGGING DISCRIMINANT VS B-HADRON MULTIPLICITY IN JET
    for(int allJetId : jetsId) {
        int genJetId = genJetIdOfRecoJet(allJets.at(allJetId), genAllJets);
        if(genJetId<0) continue;

        std::vector<int> genJetFlavours = bHadFlavoursInGenJet(genJetId, bHadJetIndex, bHadFlavour, true);

        float discrVal = allJetsBtagDiscriminant.at(allJetId);
        if(discrVal<0) discrVal = -0.01;

        int nHads = genJetFlavours.size();
        int nHadsNotG = 0;
        int nHadsT = 0;
        for(size_t flavour : genJetFlavours) {
            if(flavour!=21 && flavour!=0) nHadsNotG++;
            if(std::abs(flavour)==6) nHadsT++;
        }
        ((TH2D*)m_histogram["jetBtagDiscriminantVsHadronMultiplicity"])->Fill(discrVal, nHads, weight);
        ((TH2D*)m_histogram["jetBtagDiscriminantVsHadronTopMultiplicity"])->Fill(discrVal, nHadsT, weight);
        ((TH2D*)m_histogram["jetBtagDiscriminantVsHadronNotGMultiplicity"])->Fill(discrVal, nHadsNotG, weight);

    }
    
    
    ///////////////////////////////////////////////////////////////////////// DIJET MASS DISTRIBUTION FOR ALL AND GOOD PAIRS
    LV* dijet_trueHiggs = 0;
    if(trueHiggsJetsId.size()==2) dijet_trueHiggs = new LV(allJets.at(trueHiggsAllJetsId.at(0)) + allJets.at(trueHiggsAllJetsId.at(1)));
    
//     if(trueTopJetsId.size()>1 && trueHiggsJetsId.size()>1) {
    std::vector<std::pair<int, int> > allJetPairs = recoObjectIndices.jetIndexPairs_;
    std::vector<std::pair<int, int> > allBJetPairs;
    for(std::pair<int, int> pair : allJetPairs) {
        if(!isInVector(bJetsId, pair.first) || !isInVector(bJetsId, pair.second)) continue;
        allBJetPairs.push_back(pair);
    }
    fillDijetMassForPairs(allJets, trueHiggsJetsId, allBJetPairs, weight, m_histogram, "dijet_mass_all" );
    fillDijetMassForPairs(allJets, trueHiggsJetsId, allBJetPairs, weight, m_histogram, "dijet_mass_allNormWeight", true);
    
    
    std::vector<std::pair<int, int> > addBJetPairs;
    for(std::pair<int, int> pair : allBJetPairs) {
        if(isInVector(trueTopAllJetsId, pair.first) || isInVector(trueTopAllJetsId, pair.second)) continue;
        addBJetPairs.push_back(pair);
    }
    fillDijetMassForPairs(allJets, trueHiggsJetsId, addBJetPairs, weight, m_histogram, "dijet_mass_recoAdditionalBJets" );
    
    
    std::vector<std::pair<int, int> > goodJetPairs;
    if(recoObjectIndices.jetIndices_.size()>1) {
        goodJetPairs = jetPairsFromMVA(m_histogram, recoObjectIndices, genObjectIndices, recoObjects, trueTopJetsId, trueHiggsJetsId, weight);
    }
    fillDijetMassForPairs(allJets, trueHiggsJetsId, goodJetPairs, weight, m_histogram, "dijet_mass_goodPairs" );
    fillDijetMassForPairs(allJets, trueHiggsJetsId, goodJetPairs, weight, m_histogram, "dijet_mass_goodPairsNormWeight", true);
    
    
    // Filling true-all correlation plots
    for(std::pair<int, int> pair : allBJetPairs) {
        LV dijet_all = allJets.at(pair.first) + allJets.at(pair.second);
        for(std::pair<int, int> truePair : addBJetPairs) {
            LV dijet_true = allJets.at(truePair.first) + allJets.at(truePair.second);
            ((TH2*)m_histogram["dijet_mass_all_vs_trueAdd"])->Fill(dijet_all.M(), dijet_true.M(), weight);
        }
        if(!dijet_trueHiggs) continue;
        ((TH2*)m_histogram["dijet_mass_all_vs_trueHiggs"])->Fill(dijet_all.M(), dijet_trueHiggs->M(), weight);
    }
    
    // Filling true-good correlation plots
    for(std::pair<int, int> pair : goodJetPairs) {
        LV dijet_good = allJets.at(pair.first) + allJets.at(pair.second);
        for(std::pair<int, int> truePair : addBJetPairs) {
            LV dijet_true = allJets.at(truePair.first) + allJets.at(truePair.second);
            ((TH2*)m_histogram["dijet_mass_goodPairs_vs_trueAdd"])->Fill(dijet_good.M(), dijet_true.M(), weight);
        }
        if(!dijet_trueHiggs) continue;
        ((TH2*)m_histogram["dijet_mass_goodPairs_vs_trueHiggs"])->Fill(dijet_good.M(), dijet_trueHiggs->M(), weight);
    }
    
    std::vector<std::pair<int, int> > genAllBJetPairs;
    for(const int i: genBJetsId){
        if(!isInVector(genJetsId, i)) continue;
        for(const int j: genBJetsId){
            if(j<=i) continue;
            if(!isInVector(genJetsId, j)) continue;
            genAllBJetPairs.push_back(std::pair<int,int>(i,j));
        }
    }
    fillDijetMassForPairs(genAllJets, genHiggsAllJetsId, genAllBJetPairs, weight, m_histogram, "dijet_mass_genAllBJets" );
    
    std::vector<std::pair<int, int> > genAddBJetPairs;
    for(std::pair<int, int> pair : genAllBJetPairs) {
        if(isInVector(genTopAllJetsId, pair.first)) continue;
        if(isInVector(genTopAllJetsId, pair.second)) continue;
        genAddBJetPairs.push_back(pair);
    }
    fillDijetMassForPairs(genAllJets, genHiggsAllJetsId, genAddBJetPairs, weight, m_histogram, "dijet_mass_genAdditionalBJets" );
    
//     std::vector<std::pair<int, int> > correctJetPairs;
//     const tth::IndexPairs& jetIndexPairs = recoObjectIndices.jetIndexPairs_;
//     for(std::pair<int, int> pair : jetIndexPairs){
//         if(areAmongPairs(wrongJetPairs, pair.first, pair.second)) continue;
//         correctJetPairs.push_back(pair);
//     }
// 
//     // Plotting dijet mass for combinations of jets depending on MVA weight
//     correctPairFraction(allJets, bJetsId, allJetsBtagDiscriminant, emptyVector, trueHiggsJetsId, weight, m_histogram["dijet_mass_noWrongMVApairs"], 0, fillAllCombinations, 0.0, 1, wrongJetPairs);
//     correctPairFraction(allJets, bJetsId, allJetsBtagDiscriminant, emptyVector, trueHiggsJetsId, weight, m_histogram["dijet_mass_noCorrectMVApairs"], 0, fillAllCombinations, 0.0, 1, correctJetPairs);
    

    if(doLeadingJetsAnalysis_) fillTopAdditionalJetsHistos(eventMetadata,
                                                           recoObjects, topGenObjects, kinematicReconstructionSolutions,
                                                           recoObjectIndices, genObjectIndices, weight, m_histogram);
    
}


void AnalyzerDijet::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;

    std::stringstream ss_label;
    TString categoryName = tth::extractJetCategory(step);
    if(categoryName != ""){
        categoryName.ReplaceAll("_cate", "");
        int category = atoi(categoryName);
        ss_label<<" ["<<jetCategories_->binLabels().at(category)<<"]";
    }
    const TString label = ss_label.str();
    
    const int nBinsX_dijetM = 26;
    double binsX_dijetM[nBinsX_dijetM+1];
    for(unsigned int i=0; i<nBinsX_dijetM+1; ++i) binsX_dijetM[i]=float(i)*20.0;

    name = "allJet_multiplicity";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "All jet multiplicity;N jets_{all}"+label+";Events",20,0,20));
    name = "jet_multiplicity";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet multiplicity;N jets_{reco}"+label+";Events",20,0,20));
    name = "jet_pt";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet Pt;jet Pt{reco}"+label+";Jets",30,0,300));
    name = "jet_eta";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet Eta;jet #eta{reco}"+label+";Jets",30,-2.6,2.6));
    name = "jet_btagDiscriminator";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet btagDiscriminator d;jet d{reco}"+label+";Jets",36,-0.1,1.1));
    name = "jet_partonFlavour";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet partonFlavour;jet flavour_{reco}^{parton}"+label+";Jets",20,-10,10));

    name = "bJet_multiplicity";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "B-jet multiplicity;N b-jets_{reco}"+label+";Events",20,0,20));
    name = "bJet_multiplicity_pt_0_30";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "B-jet multiplicity (0 #lt p_{T} #lt 30);N b-jets_{reco}^{0 #lt p_{T} #lt 30}"+label+";Events",15,0,15));
    name = "bJet_multiplicity_pt_30_70";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "B-jet multiplicity (30 #lt p_{T} #lt 70);N b-jets_{reco}^{30 #lt p_{T} #lt 70}"+label+";Events",15,0,15));
    name = "bJet_multiplicity_pt_70_500";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "B-jet multiplicity (70 #lt p_{T} #lt 500);N b-jets_{reco}^{70 #lt p_{T} #lt 500}"+label+";Events",15,0,15));
    name = "bJet_pt";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "B-jet Pt;bJet Pt{reco}"+label+";Jets",30,0,300));
    name = "bJet_eta";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "B-jet Eta;bJet #eta{reco}"+label+";Jets",30,-2.6,2.6));
    name = "bJet_btagDiscriminator";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "B-jet btagDiscriminator d;bJet d{reco}"+label+";Jets",36,-0.1,1.1));
    name = "bJet_partonFlavour";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "B-jet partonFlavour;bJet flavour_{reco}^{parton}"+label+";Jets",20,-10,10));
    name = "bJet_dRmin";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "B-jet dRmin;bJet dR_min{reco}"+label+";Events",50,0,3.5));

    name = "genJet_multiplicity";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Gen. jet multiplicity;N jets_{gen}"+label+";Events",50,0,50));
    name = "genBJet_multiplicity";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Gen. b-jet multiplicity;N(b-jet)_{gen}"+label+";Events",15,0,15));
    name = "genBJet_pt";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Gen. b-jet Pt;bJet Pt{gen}"+label+";Jets",30,0,300));
    name = "genBJet_eta";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Gen. b-jet Eta;bJet #eta{gen}"+label+";Jets",30,-2.6,2.6));
    name = "genBJet_dRmin";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Gen. b-jet dRmin;bJet dR_min_{gen}"+label+";Events",50,0,3.5));
    name = "genBJet_recoJet_dRmin";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Gen. b-jet - Reco jet dRmin;bJet_{gen}-Jet_{reco} dR_min"+label+";Events",50,0,3.5));
    name = "genBJet_recoBJet_dRmin";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Gen. b-jet - Reco b-jet dRmin;bJet_{gen}-bJet_{reco} dR_min"+label+";Events",50,0,3.5));
    name = "bHad_multiplicity";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Gen. b-had multiplicity;N(b-had)_{gen}"+label+";Events",15,0,15));
    name = "bHad_top_multiplicity";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Gen. b-had multiplicity;N(b-had)_{gen}^{top}"+label+";Events",6,0,6));
    name = "bHad_higgs_multiplicity";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Gen. b-had multiplicity;N(b-had)_{gen}^{higgs}"+label+";Events",6,0,6));
    name = "bHad_flavour";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Gen. b-had flavour;Flavour(b-had)_{gen}"+label+";B-hadrons",60,-30,30));
    name = "bHad_dRmin";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Gen. b-had dRmin;dR_min(b-had)_{gen}"+label+";Jets",50,0,3.5));
    name = "bHadNoG_multiplicity";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Gen. b-had multiplicity (not from gluons);N(b-had NoG)_{gen}"+label+";Events",15,0,15));
    name = "bHadVsBJet_multiplicity";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "Gen. b-had/b-jet multiplicity;N(b-had)_{gen}"+label+";N(b-jet)_{reco}",15,0,15,15,0,15));
    name = "bHadVsBHadNoT_multiplicity";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "Gen. b-had/b-had^{not t#bar{t}} multiplicity;N(b-had)_{all}"+label+";N(b-had)_{not t#bar{t}}",15,0,15,15,0,15));

    name = "topJet_multiplicity_allTrue";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Top jet multiplicity (true);N jets_{top}^{true} (all)"+label+";Events",5,0,5));
    name = "topJet_multiplicity_true";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Top jet multiplicity (true);N jets_{top}^{true} (sel.)"+label+";Events",5,0,5));
    name = "topJet_multiplicity_reco";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Top jet multiplicity (reco);N jets_{top}^{reco} (sel.)"+label+";Events",5,0,5));
    name = "topJet_multiplicity_recoInAllTrue";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Top jet multiplicity (reco in true);N jets_{top}^{reco} (among all true)"+label+";Events",5,0,5));
    name = "topJet_multiplicity_recoInTrue";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Top jet multiplicity (reco in true);N jets_{top}^{reco} (among sel. true)"+label+";Events",5,0,5));
    name = "topJet_multiplicity_gen";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Top jet multiplicity (gen);N jets_{top}^{gen}"+label+";Events",5,0,5));
    name = "topJet_multiplicity_allGen";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Top jet multiplicity (gen);N jets_{top}^{all gen}"+label+";Events",5,0,5));
    name = "topBJet_multiplicity_true";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Top b-jet multiplicity (true);N b-tagged jets_{top}^{true} (sel. b-tagged)"+label+";Events",5,0,5));
    name = "topBJet_multiplicity_reco";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Top b-jet multiplicity (true);N b-tagged jets_{top}^{reco} (sel. b-tagged)"+label+";Events",5,0,5));
    name = "higgsJet_multiplicity_allTrue";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Higgs jet multiplicity (true);N jets_{higgs}^{true} (all)"+label+";Events",5,0,5));
    name = "higgsJet_multiplicity_true";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Higgs jet multiplicity (true);N jets_{higgs}^{true} (sel.)"+label+";Events",5,0,5));
    name = "higgsBJet_multiplicity_true";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Higgs b-jet multiplicity (true);N b-tagged jets_{higgs}^{true} (sel. b-tagged)"+label+";Events",5,0,5));
    name = "higgsJet_multiplicity_gen";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Higgs jet multiplicity (gen);N jets_{higgs}^{gen}"+label+";Events",5,0,5));
    name = "higgsJet_multiplicity_allGen";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Higgs jet multiplicity (gen);N jets_{higgs}^{all gen}"+label+";Events",5,0,5));
    
    name = "topBJet_Pt_gen";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Top b-jet_{gen} Pt;Pt_{jet}^{gen}"+label+";Top Jets",50,0,100));
    name = "topBJet_Pt_true";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Top b-jet_{reco} Pt;Pt_{jet}^{reco}"+label+";Top Jets",50,0,100));
    name = "topBJet_Btag_true";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Top b-jet_{reco} bTag disriminant;bTag_{jet}^{reco}"+label+";Top Jets",20,0,1));
    
    name = "higgsBJet_Pt_gen";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Higgs b-jet_{gen} Pt;Pt_{jet}^{gen}"+label+";Higgs Jets",50,0,100));
    name = "higgsBJet_Pt_true";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Higgs b-jet_{reco} Pt;Pt_{jet}^{reco}"+label+";Higgs Jets",50,0,100));
    name = "higgsBJet_Btag_true";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Higgs b-jet_{reco} bTag disriminant;bTag_{jet}^{reco}"+label+";Higgs Jets",20,0,1));

    name = "topLeptonAntilepton_dR";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Top lepton - antiLepton dR;Lepton_{top}-Antilepton_{antitop} dR"+label+";Events",50,0,3.5));
    name = "topDiLeptonMet_dR";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Top diLepton - Met dR;diLepton_{top}-Met_{event} dR"+label+";Events",50,0,3.5));

    name = "dijet_correctJet_multiplicity_all";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Number of correct jets;N jets^{correct} (all)"+label+";Jet pairs",4,0,4));
    name = "dijet_correctJet_multiplicity_allM";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Number of correct jets;N jets^{correct} (all M)"+label+";Jet pairs",4,0,4));
    name = "dijet_correctJet_multiplicity_trueTopJets";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Number of correct jets;N jets^{correct} (trueTopJets)"+label+";Jet pairs",4,0,4));
    name = "dijet_correctJet_multiplicity_recoTopJets";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Number of correct jets;N jets^{correct} (recoTopJets)"+label+";Jet pairs",4,0,4));

//     name = "dijet_correctJet_multiplicity_all_jetPtLt30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Number of correct jets;N jets^{correct}_{Pt<30} (all)"+label+";Jet pairs",4,0,4));
//     name = "dijet_correctJet_multiplicity_allM_jetPtLt30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Number of correct jets;N jets^{correct}_{Pt<30} (all M)"+label+";Jet pairs",4,0,4));
//     name = "dijet_correctJet_multiplicity_trueTopJets_jetPtLt30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Number of correct jets;N jets^{correct}_{Pt<30} (trueTopJets)"+label+";Jet pairs",4,0,4));
//     name = "dijet_correctJet_multiplicity_recoTopJets_jetPtLt30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Number of correct jets;N jets^{correct}_{Pt<30} (recoTopJets)"+label+";Jet pairs",4,0,4));
// 
//     name = "dijet_correctJet_multiplicity_all_jetPtGe30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Number of correct jets;N jets^{correct}_{Pt>=30} (all)"+label+";Jet pairs",4,0,4));
//     name = "dijet_correctJet_multiplicity_allM_jetPtGe30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Number of correct jets;N jets^{correct}_{Pt>=30} (all M)"+label+";Jet pairs",4,0,4));
//     name = "dijet_correctJet_multiplicity_trueTopJets_jetPtGe30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Number of correct jets;N jets^{correct}_{Pt>=30} (trueTopJets)"+label+";Jet pairs",4,0,4));
//     name = "dijet_correctJet_multiplicity_recoTopJets_jetPtGe30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Number of correct jets;N jets^{correct}_{Pt>=30} (recoTopJets)"+label+";Jet pairs",4,0,4));

    name = "dijet_correctPairFraction_all";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Correct pair fraction;correct/all (all)"+label+";Events",15,-0.3,1.2));
    name = "dijet_correctPairFraction_allM";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Correct pair fraction;correct/all (all M)"+label+";Events",15,-0.3,1.2));
    name = "dijet_correctPairFraction_trueTopJets";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Correct pair fraction;correct/all (trueTopJets)"+label+";Events",15,-0.3,1.2));
    name = "dijet_correctPairFraction_recoTopJets";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Correct pair fraction;correct/all (recoTopJets)"+label+";Events",15,-0.3,1.2));

//     name = "dijet_correctPairFraction_all_jetPtLt30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Correct pair fraction;correct/all_{Pt<30} (all)"+label+";Events",15,-0.3,1.2));
//     name = "dijet_correctPairFraction_allM_jetPtLt30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Correct pair fraction;correct/all_{Pt<30} (all M)"+label+";Events",15,-0.3,1.2));
//     name = "dijet_correctPairFraction_trueTopJets_jetPtLt30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Correct pair fraction;correct/all_{Pt<30} (trueTopJets)"+label+";Events",15,-0.3,1.2));
//     name = "dijet_correctPairFraction_recoTopJets_jetPtLt30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Correct pair fraction;correct/all_{Pt<30} (recoTopJets)"+label+";Events",15,-0.3,1.2));
// 
//     name = "dijet_correctPairFraction_all_jetPtGe30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Correct pair fraction;correct/all_{Pt>=30} (all)"+label+";Events",15,-0.3,1.2));
//     name = "dijet_correctPairFraction_allM_jetPtGe30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Correct pair fraction;correct/all_{Pt>=30} (all M)"+label+";Events",15,-0.3,1.2));
//     name = "dijet_correctPairFraction_trueTopJets_jetPtGe30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Correct pair fraction;correct/all_{Pt>=30} (trueTopJets)"+label+";Events",15,-0.3,1.2));
//     name = "dijet_correctPairFraction_recoTopJets_jetPtGe30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Correct pair fraction;correct/all_{Pt>=30} (recoTopJets)"+label+";Events",15,-0.3,1.2));

    name = "dijet_mass_all";
    bookPairHistos(new TH1D(prefix_+name+step, "Dijet mass;M(dijet) (all)"+label+";Jet pairs",nBinsX_dijetM,binsX_dijetM), m_histogram, name);
    name = "dijet_mass_all_nEntriesPerEvent";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "N jet pairs per event;N jet pairs/event (all)"+label+";Events",16,0,16));
    name = "dijet_mass_all_diffMin";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass difference_{min};M1_{dijet}-M2_{dijet} (all)"+label+";Events_{N pairs>1}",100,0,500));
    name = "dijet_mass_all_nEntriesPerBin";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "N jet pairs per bin;M(dijet) (all)"+label+";N jet pairs",nBinsX_dijetM,binsX_dijetM, 10, 0, 10));
    name = "dijet_mass_all_vs_trueHiggs";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "Dijet mass;M(dijet)_{all}"+label+";M(dijet)_{true}^{Higgs}",nBinsX_dijetM,binsX_dijetM, nBinsX_dijetM,binsX_dijetM));
    name = "dijet_mass_all_vs_trueAdd";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "Dijet mass;M(dijet)_{all}"+label+";M(dijet)_{true}^{add}",nBinsX_dijetM,binsX_dijetM, nBinsX_dijetM,binsX_dijetM));
    
    name = "dijet_mass_allNormWeight";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass;M(dijet) (all)"+label+";Jet pairs",nBinsX_dijetM,binsX_dijetM));
    
    name = "dijet_mass_allM";
    bookPairHistos(new TH1D(prefix_+name+step, "Dijet mass;M(dijet) (all M)"+label+";Jet pairs",nBinsX_dijetM,binsX_dijetM), m_histogram, name);
    name = "dijet_mass_trueTopJets";
    bookPairHistos(new TH1D(prefix_+name+step, "Dijet mass;M(dijet) (trueTopJets)"+label+";Jet pairs",nBinsX_dijetM,binsX_dijetM), m_histogram, name);
    name = "dijet_mass_recoTopJets";
    bookPairHistos(new TH1D(prefix_+name+step, "Dijet mass;M(dijet) (recoTopJets)"+label+";Jet pairs",nBinsX_dijetM,binsX_dijetM), m_histogram, name);
    name = "dijet_mass_genTopJets";
    bookPairHistos(new TH1D(prefix_+name+step, "Dijet mass;M(dijet) (genTopJets)"+label+";Jet pairs",nBinsX_dijetM,binsX_dijetM), m_histogram, name);

//     name = "dijet_mass_all_jetPtLt30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass;M(dijet)_{Pt<30} (all)"+label+";Jet pairs",25,0,500));
//     name = "dijet_mass_allM_jetPtLt30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass;M(dijet)_{Pt<30} (all M)"+label+";Jet pairs",25,0,500));
//     name = "dijet_mass_trueTopJets_jetPtLt30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass;M(dijet)_{Pt<30} (trueTopJets)"+label+";Jet pairs",25,0,500));
//     name = "dijet_mass_recoTopJets_jetPtLt30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass;M(dijet)_{Pt<30} (recoTopJets)"+label+";Jet pairs",25,0,500));
// 
//     name = "dijet_mass_all_jetPtGe30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass;M(dijet)_{Pt>=30} (all)"+label+";Jet pairs",25,0,500));
//     name = "dijet_mass_allM_jetPtGe30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass;M(dijet)_{Pt>=30} (all M)"+label+";Jet pairs",25,0,500));
//     name = "dijet_mass_trueTopJets_jetPtGe30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass;M(dijet)_{Pt>=30} (trueTopJets)"+label+";Jet pairs",25,0,500));
//     name = "dijet_mass_recoTopJets_jetPtGe30";
//     m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass;M(dijet)_{Pt>=30} (recoTopJets)"+label+";Jet pairs",25,0,500));

    name = "dijet_mass_trueHiggsJets";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass;M(dijet)_{H} (true reco)"+label+";Jet pairs",nBinsX_dijetM,binsX_dijetM));
    name = "dijet_mass_genHiggsJets";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass;M(dijet)_{H} (true gen)"+label+";Jet pairs",nBinsX_dijetM,binsX_dijetM));
    
    name = "dijet_mass_genAdditionalBJets";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass;M(dijet)_{add} (true gen)"+label+";Jet pairs",nBinsX_dijetM,binsX_dijetM));
    name = "dijet_mass_genAllBJets";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass;M(dijet)_{all} (true gen)"+label+";Jet pairs",nBinsX_dijetM,binsX_dijetM));
    name = "dijet_mass_recoAdditionalBJets";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass;M(dijet) (recoAdditionalJets)"+label+";Jet pairs",nBinsX_dijetM,binsX_dijetM));

    name = "jet_PtLt30_multiplicity";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet multiplicity (Pt<30);N jets_{reco}^{Pt<30}"+label+";Events",20,0,20));
    name = "bJet_PtLt30_multiplicity";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "B-jet multiplicity (Pt<30);N b-jets_{reco}^{Pt<30}"+label+";Events",20,0,20));
    name = "jet_PtLt30_flavour";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet flavour (Pt<30);flavour(jet)_{Pt<30}"+label+";Jets",30,0,30));
    name = "bJet_PtLt30_flavour";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "B-jet flavour (Pt<30);flavour(b-jet)_{Pt<30}"+label+";B-tagged jets",30,0,30));
    name = "jet_PtLt30_btagDiscriminator";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet btagDiscriminator d (Pt<30);d_{reco}^{Pt<30}"+label+";Jets",20,0,1));
    name = "bJet_PtLt30_btagDiscriminator";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "B-jet btagDiscriminator d (Pt<30);d_{reco}^{Pt<30}"+label+";B-tagged jets",20,0,1));


    name = "jetBtagDiscriminantVsHadronMultiplicity";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "Jet_{d}^{reco} vs Hadron multiplicity;jet_{reco} d"+label+";N(b-had)_{all}",21,-0.05,1,7,0,7));
    name = "jetBtagDiscriminantVsHadronTopMultiplicity";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "Jet_{d}^{reco} vs Hadron_{tt} multiplicity;jet_{reco} d"+label+";N(b-had)_{tt}",21,-0.05,1,7,0,7));
    name = "jetBtagDiscriminantVsHadronNotGMultiplicity";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "Jet_{d}^{reco} vs Hadron_{not g} multiplicity;jet_{reco} d"+label+";N(b-had)_{not g}",21,-0.05,1,7,0,7));

    name = "mvaWeight_correctTop_correct";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "MVA^{correct} weight for correct tt combinations;weight_{MVA}^{correct}"+label+";Correct jet pairs",20,-1.2,0.2));
    name = "mvaWeight_correctTop_swapped";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "MVA^{swapped} weight for correct tt combinations;weight_{MVA}^{swapped}"+label+";Correct jet pairs",20,-1.2,0.2));
    
    name = "mvaWeight_wrongTop_correct";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "MVA^{correct} weight for wrong tt combinations;weight_{MVA}^{correct}"+label+";Wrong jet pairs",20,-1.2,0.2));

    name = "mvaWeight_correctHiggs_correct";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "MVA^{correct} weight for correct H combinations;weight_{MVA}^{correct}"+label+";Correct jet pairs",20,-1.2,0.2));
    name = "mvaWeight_correctHiggs_swapped";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "MVA^{swapped} weight for correct H combinations;weight_{MVA}^{swapped}"+label+";Correct jet pairs",20,-1.2,0.2));
    
    name="mvaWeight_correctTop_correctHiggs";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "MVA_{top}^{higgs} weight for correct combinations;MVA weight_{correct}^{top}"+label+";MVA weight_{correct}^{higgs}",20,-1.2,0.2,20,-1.2,0.2));
    name = "mvaWeight_wrongTopAndHiggs_correct";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "MVA^{correct} weight for combinations not from Top, Higgs;weight_{MVA}^{correct}"+label+";Wrong jet pairs",20,-1.2,0.2));

    name = "dijet_mass_goodPairs";
    bookPairHistos(new TH1D(prefix_+name+step, "Dijet mass (no wrong pairs from MVA);M(dijet) (MVA)"+label+";Jet pairs",nBinsX_dijetM,binsX_dijetM), m_histogram, name);
    name = "dijet_mass_goodPairs_nEntriesPerEvent";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass (now wrong pairs from MVA);N jet pairs/event (MVA)+"+label+";Events",16,0,16));
    name = "dijet_mass_goodPairs_diffMin";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass difference_{min} (no wrong pairs from MVA);M1_{dijet}-M2_{dijet} (MVA_{correct}^{high})"+label+";Events_{N pairs>1}",100,0,500));
    name = "dijet_mass_goodPairs_nEntriesPerBin";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "N jet pairs per bin;M(dijet) (MVA)"+label+";N jets pairs",nBinsX_dijetM,binsX_dijetM, 10, 0, 10));
    name = "dijet_mass_goodPairs_vs_trueHiggs";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "Dijet mass;M(dijet)_{MVA}^{add}"+label+";M(dijet)_{true}^{Higgs}",nBinsX_dijetM,binsX_dijetM, nBinsX_dijetM,binsX_dijetM));
    name = "dijet_mass_goodPairs_vs_trueAdd";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "Dijet mass;M(dijet)_{MVA}^{add}"+label+";M(dijet)_{true}^{add}",nBinsX_dijetM,binsX_dijetM, nBinsX_dijetM,binsX_dijetM));
    
    name = "dijet_mass_goodPairsNormWeight";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass (no correct pairs from MVA);M(dijet) (MVA)"+label+";Jet pairs",nBinsX_dijetM,binsX_dijetM));

    name = "weight";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Weight of the event;weight"+label+";Events",30,0,3));
    
    name = "nHadAddFromTopWeakDecay";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "N b-hadrons not from top but from top weak decay;N b-hadrons_{extra}^{from top}"+label+";Events",10,0,10));
    name = "nHadAddNotFromTopWeakDecay";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "N b-hadrons not from top and not from top weak decay;N b-hadrons_{extra}^{real}"+label+";Events",10,0,10));
    name = "nGenJetAddFromTopWeakDecay";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "N b-jets not from top but from top weak decay;N b-jets_{extra}^{from top}"+label+";Events",10,0,10));
    name = "nGenJetAddNotFromTopWeakDecay";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "N b-jets not from top and not from top weak decay;N b-jets_{extra}^{real}"+label+";Events",10,0,10));
    name = "nGenJetAddFromTopWeakDecay_signal";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "N b-jets not from top but from top weak decay;(p_{T}>20,|#eta|<2.5) N b-jets_{extra}^{from top}"+label+";Events",10,0,10));
    name = "nGenJetAddNotFromTopWeakDecay_signal";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "N b-jets not from top and not from top weak decay;(p_{T}>20,|#eta|<2.5) N b-jets_{extra}^{real}"+label+";Events",10,0,10));
    
    name = "extraGenBHadronMultiplicity";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "N additional gen. b-hadrons;N hadrons"+label+";Events",10,0,10));
    
    name = "extraGenBJetMultiplicity";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "N additional gen. b-hadrons;N jets"+label+";Events",10,0,10));
    
    name = "extraGenBHadronsJets";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "N additional gen. b-jets missed due to different reasons;Requirement passed"+label+";B-jets_{gen}",10,0,10));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "Additional hadron");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "Clustered to jet");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "Jet in acceptance");
    m_histogram[name]->GetXaxis()->SetBinLabel(4, "No overlap with top jets");
    m_histogram[name]->GetXaxis()->SetBinLabel(5, "No overlap with extra jets");
    
            
    name = "nBHadAddInGenBJet_top";
    m_histogram[name] = store(new TH1D(prefix_+name+step, ";N b-hadron_{add}"+label+";b-jet_{t#bar{t}}",5,0,5));
    name = "nBHadAddInGenBJet_add";
    m_histogram[name] = store(new TH1D(prefix_+name+step, ";N b-hadron_{add}"+label+";b-jet_{add}",5,0,5));


    if(doHadronMatchingComparison_) {
        name = "bQuarkT_bHad_dRmin";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Min. dR b-quark_{tt} ^ bHadron;dR(bQ,bHad)^{tt}"+label+";b-quarks_{tt}",50,0,2));
        name = "bQuarkH_bHad_dRmin";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Min. dR b-quark_{H} ^ bHadron;dR(bQ,bHad)^{H}"+label+";b-quarks_{H}",50,0,2));
        name = "bQuarkT_bHad_dRmatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "dR b-quark_{tt} ^ bHadron_{tt};dR(bQ,bHad)^{tt}"+label+";b-quarks_{tt}",50,0,2));
        name = "bQuarkH_bHad_dRmatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "dR b-quark_{H} ^ bHadron_{H};dR(bQ,bHad)^{H}_{matched}"+label+";b-quarks_{H}",50,0,2));
        name = "bHadT_genJet_dRmin";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Min. dR b-hadron_{tt} ^ jet_{gen};dR(bHad,genJet)^{tt}"+label+";b-hadrons_{tt}",50,0,2));
        name = "bHadH_genJet_dRmin";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Min. dR b-hadron_{H} ^ jet_{gen};dR(bHad,genJet)^{H}"+label+";b-hadrons_{H}",50,0,2));
        name = "bHadT_genJet_dRmin_05";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Min. dR b-hadron_{tt} ^ jet_{gen} (dR<0.5);dR(bHad,genJet)^{tt}"+label+";b-hadrons_{tt}",50,0,2));
        name = "bHadH_genJet_dRmin_05";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Min. dR b-hadron_{H} ^ jet_{gen} (dR<0.5);dR(bHad,genJet)^{H}"+label+";b-hadrons_{H}",50,0,2));
        name = "bHadT_genJet_dRmatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "dR b-hadron_{tt} ^ jet_{gen};dR(bHad,genJet)^{tt}"+label+";b-hadrons_{tt}",50,0,2));
        name = "bHadT_genJet_dRmin_matched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Min. dR b-hadron_{tt} ^ jet_{gen};dR(bHad,genJet)^{tt}"+label+";b-hadrons_{tt}^{matched}",50,0,2));
        name = "bHadH_genJet_dRmatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "dR b-hadron_{H} ^ jet_{gen};dR(bHad,genJet)^{H}"+label+";b-hadrons_{H}",50,0,2));
        name = "bHadH_genJet_dRmin_matched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Min. dR b-hadron_{H} ^ jet_{gen};dR(bHad,genJet)^{H}"+label+";b-hadrons_{tt}^{matched}",50,0,2));
        name = "bHadG_genJet_dRmatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "dR b-hadron_{gluon} ^ jet_{gen};dR(bHad,genJet)^{g}"+label+";b-hadrons_{gluon}",50,0,2));
        name = "bHadG_genJet_dRmin_matched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Min. dR b-hadron_{gluon} ^ jet_{gen};dR(bHad,genJet)^{H}g"+label+";b-hadrons_{gluon}^{matched}",50,0,2));
        name = "bQuarkT_genJet_dRmin";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Min. dR b-quark_{tt} ^ jet_{gen};dR(bQ,genJet)^{tt}"+label+";b-quarks_{tt}",50,0,2));
        name = "bQuarkH_genJet_dRmin";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Min. dR b-quark_{H} ^ jet_{gen};dR(bQ,genJet)^{H}"+label+";b-quarks_{H}",50,0,2));
        name = "bQuarkT_genJet_dRmin_05";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Min. dR b-quark_{tt} ^ jet_{gen} (dR<0.5);dR(bQ,genJet)^{tt}"+label+";b-quarks_{tt}",50,0,2));
        name = "bQuarkH_genJet_dRmin_05";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Min. dR b-quark_{H} ^ jet_{gen} (dR<0.5);dR(bQ,genJet)^{H}"+label+";b-quarks_{H}",50,0,2));
        name = "bQuarkT_genJet_dRmatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "dR b-quark_{tt} ^ jet_{gen};dR(bQ,genJet)^{tt}"+label+";b-quarks_{tt}",50,0,2));
        name = "bQuarkTpt_genJet_dRmatched_low";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Pt b-quark_{tt} ^ jet_{gen};dR(bQ,genJet)^{tt}<1"+label+";b-quarks_{tt}",50,0,500));
        name = "bQuarkTpt_genJet_dRmatched_high";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Pt b-quark_{tt} ^ jet_{gen};dR(bQ,genJet)^{tt}>1"+label+";b-quarks_{tt}",50,0,500));
        name = "bQuarkHpt_genJet_dRmatched_low";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Pt b-quark_{H} ^ jet_{gen};dR(bQ,genJet)^{H}<1"+label+";b-quarks_{H}",50,0,500));
        name = "bQuarkHpt_genJet_dRmatched_high";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Pt b-quark_{H} ^ jet_{gen};dR(bQ,genJet)^{H}>1"+label+";b-quarks_{H}",50,0,500));

        name = "bQuarkT_genJet_directQ_dRmin";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Min. dR b-quark_{tt} ^ jet_{gen};dR(bQ,genJet)^{tt}"+label+";b-quarks_{tt}",50,0,2));
        name = "bQuarkH_genJet_directQ_dRmin";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Min. dR b-quark_{H} ^ jet_{gen};dR(bQ,genJet)^{H}"+label+";b-quarks_{H}",50,0,2));


        name = "bHadT_unique_multiplicity_dR";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "B-hadron multiplicity (tt);N hadrons_{unique}^{tt}"+label+";Events",6,0,6));
        name = "bHadH_unique_multiplicity_dR";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "B-hadron multiplicity (H);N hadrons_{unique}^{H}"+label+";Events",6,0,6));
        name = "bHadTH_unique_multiplicity_dR";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "B-hadron multiplicity (ttH);N hadrons_{unique}^{ttH}"+label+";Events",6,0,6));
        name = "bHadT_unique_multiplicity_dR_05";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "B-hadron multiplicity (tt);N hadrons_{unique}^{tt}"+label+";Events",6,0,6));
        name = "bHadH_unique_multiplicity_dR_05";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "B-hadron multiplicity (H);N hadrons_{unique}^{H}"+label+";Events",6,0,6));
        name = "bHadTH_unique_multiplicity_dR_05";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "B-hadron multiplicity (ttH);N hadrons_{unique}^{ttH}"+label+";Events",6,0,6));
        name = "bHadT_unique_multiplicity_match";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "B-hadron multiplicity (tt);N hadrons_{unique}^{tt}"+label+";Events",6,0,6));
        name = "bHadH_unique_multiplicity_match";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "B-hadron multiplicity (H);N hadrons_{unique}^{H}"+label+";Events",6,0,6));
        name = "bHadTH_unique_multiplicity_match";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "B-hadron multiplicity (ttH);N hadrons_{unique}^{ttH}"+label+";Events",6,0,6));
        name = "genJetT_unique_multiplicity_dR";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet_{gen} multiplicity (tt);N genJets_{unique}^{tt}"+label+";Events",6,0,6));
        name = "genJetH_unique_multiplicity_dR";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet_{gen} multiplicity (H);N genJets_{unique}^{H}"+label+";Events",6,0,6));
        name = "genJetTH_unique_multiplicity_dR";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet_{gen} multiplicity (ttH);N genJets_{unique}^{ttH}"+label+";Events",6,0,6));
        name = "genJetT_unique_multiplicity_dR_05";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet_{gen} multiplicity (tt);N genJets_{unique}^{tt}"+label+";Events",6,0,6));
        name = "genJetH_unique_multiplicity_dR_05";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet_{gen} multiplicity (H);N genJets_{unique}^{H}"+label+";Events",6,0,6));
        name = "genJetTH_unique_multiplicity_dR_05";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet_{gen} multiplicity (ttH);N genJets_{unique}^{ttH}"+label+";Events",6,0,6));
        name = "genJetT_unique_multiplicity_match";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet_{gen} multiplicity (tt);N genJets_{unique}^{tt}"+label+";Events",6,0,6));
        name = "genJetH_unique_multiplicity_match";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet_{gen} multiplicity (H);N genJets_{unique}^{H}"+label+";Events",6,0,6));
        name = "genJetTH_unique_multiplicity_match";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet_{gen} multiplicity (ttH);N genJets_{unique}^{ttH}"+label+";Events",6,0,6));

        name = "genJetT_unique_multiplicity_directQ_dR_05";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet_{gen} multiplicity (tt);N genJets_{unique}^{tt}"+label+";Events",6,0,6));
        name = "genJetH_unique_multiplicity_directQ_dR_05";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet_{gen} multiplicity (H);N genJets_{unique}^{H}"+label+";Events",6,0,6));
        name = "genJetTH_unique_multiplicity_directQ_dR_05";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Jet_{gen} multiplicity (ttH);N genJets_{unique}^{ttH}"+label+";Events",6,0,6));

        name = "bQuarkT_genJet_dRmin_05_notAsMatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "dR b-quark_{tt} ^ jet_{gen} (if not as matched);dR(bQ,genJet)^{tt}"+label+";b-quarks_{tt}",50,0,2));
        name = "bQuarkH_genJet_dRmin_05_notAsMatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "dR b-quark_{H} ^ jet_{gen} (if not as matched);dR(bQ,genJet)^{H}"+label+";b-quarks_{H}",50,0,2));
        name = "bQuarkTH_genJet_dRmin_05_notAsMatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "dR b-quark_{ttH} ^ jet_{gen} (if not as matched);dR(bQ,genJet)^{ttH}"+label+";b-quarks_{ttH}",50,0,2));
        name = "bQuarkT_genJet_dRmin_05_isAsMatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Is b-quark_{tt} - jet_{gen} matched as with QdR;dR(bQ,genJet)^{tt}"+label+";b-quarks_{tt}",3,0,3));
        name = "bQuarkH_genJet_dRmin_05_isAsMatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Is b-quark_{H} - jet_{gen} matched as with QdR;dR(bQ,genJet)^{H}"+label+";b-quarks_{H}",3,0,3));
        name = "bQuarkTH_genJet_dRmin_05_isAsMatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Is b-quark_{ttH} - jet_{gen} matched as with QdR;dR(bQ,genJet)^{ttH}"+label+";b-quarks_{ttH}",3,0,3));

        name = "bQuarkT_genJet_directQ_dRmin_05_notAsMatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "dR b-quark_{tt} ^ jet_{gen} (if not as matched);dR(bQ,genJet)^{tt}"+label+";b-quarks_{tt}",50,0,2));
        name = "bQuarkH_genJet_directQ_dRmin_05_notAsMatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "dR b-quark_{H} ^ jet_{gen} (if not as matched);dR(bQ,genJet)^{H}"+label+";b-quarks_{H}",50,0,2));
        name = "bQuarkTH_genJet_directQ_dRmin_05_notAsMatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "dR b-quark_{ttH} ^ jet_{gen} (if not as matched);dR(bQ,genJet)^{ttH}"+label+";b-quarks_{ttH}",50,0,2));
        name = "bQuarkT_genJet_directQ_dRmin_05_isAsMatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Is b-quark_{tt} - jet_{gen} matched as with direct QdR;dR(bQ,genJet)^{tt}"+label+";b-quarks_{tt}",3,0,3));
        name = "bQuarkH_genJet_directQ_dRmin_05_isAsMatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Is b-quark_{H} - jet_{gen} matched as with direct QdR;dR(bQ,genJet)^{H}"+label+";b-quarks_{H}",3,0,3));
        name = "bQuarkTH_genJet_directQ_dRmin_05_isAsMatched";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "Is b-quark_{ttH} - jet_{gen} matched as with direct QdR;dR(bQ,genJet)^{ttH}"+label+";b-quarks_{ttH}",3,0,3));
    }
    
    if(doLeadingJetsAnalysis_) {
        
        const int nBins_Pt = 20;
//         const double bins_Pt[nBins_Pt+1] = {0., 20., 40., 60., 80., 100., 120., 140., 160., 180., 200., 220., 240., 260., 280., 300., 320., 360., 420.};
        
        // Histograms with selection at gen level and matching to reco level
        // Generic info about top and additional jets
        bookLeadingJetsHistos(m_histogram, "gen", step, label);
        bookLeadingJetsHistos(m_histogram, "true", step, label);
        bookLeadingJetsHistos(m_histogram, "kinReco", step, label);
        bookLeadingJetsHistos(m_histogram, "mva", step, label);
        // Control plots
        bookLeadingJetsCPHistos(m_histogram, "true", step, label);
        bookLeadingJetsCPHistos(m_histogram, "kinReco", step, label);
        bookLeadingJetsCPHistos(m_histogram, "mva", step, label);
        // Correlation with gen level jets
        bookAddGenJetsCorrelationHistos(m_histogram, "true", step, label, true);
        bookAddGenJetsCorrelationHistos(m_histogram, "kinReco", step, label);
        bookAddGenJetsCorrelationHistos(m_histogram, "mva", step, label);
        
        // Histograms with selection at reco level and matching to gen level
        // Generic info about top and additional jets
        bookLeadingJetsHistos(m_histogram, "vis_kinReco", step, label);
        bookLeadingJetsHistos(m_histogram, "vis_mva", step, label);
        // Correlation with gen level jets
        bookAddGenJetsCorrelationHistos(m_histogram, "vis_kinReco", step, label);
        bookAddGenJetsCorrelationHistos(m_histogram, "vis_mva", step, label);

        // True-KinReco-MVA comparison of top jets
        name = "jet_top_mult_genAll";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "; N jet^{top}_{genAll}"+label+"; Events",5,0,5));
        name = "jet_top_mult_gen";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "; N jet^{top}_{gen}"+label+"; Events",5,0,5));
        name = "jet_top_mult_true";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "; N jet^{top}_{true}"+label+"; Events",5,0,5));
        name = "jet_top_mult_kinReco";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "; N jet^{top}_{kinReco}"+label+"; Events",5,0,5));
        name = "jet_top_mult_mva";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "; N jet^{top}_{mva}"+label+"; Events",5,0,5));
        name = "bJet_top_mult_true";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "; N b-jet^{top}_{true}"+label+"; Events",5,0,5));
        name = "bJet_top_mult_kinReco";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "; N b-jet^{top}_{kinReco}"+label+"; Events",5,0,5));
        name = "bJet_top_mult_mva";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "; N b-jet^{top}_{mva}"+label+"; Events",5,0,5));
        
        name = "bJet_add_mult_gen";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "; N b-jet^{add}_{gen}"+label+"; Events",10,0,10));
        
        name = "leadingJet_1st_nHad_addB_gen";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "; N b hadrons/jet"+label+"; B jets^{add}_{gen}", 5,0,5));
        name = "leadingJet_2nd_nHad_addB_gen";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "; N b hadrons/jet"+label+"; B jets^{add}_{gen}", 5,0,5));
        
        
        // For jets assigned to tt
        name = "jet_top_id_trueVsKinReco";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; ID jet^{top}_{true}"+label+"; ID jet^{top}_{kinReco}",16,-1,15,16,-1,15));
        name = "jet_top_id_trueVsMva";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; ID jet^{top}_{true}"+label+"; ID jet^{top}_{mva}",16,-1,15,16,-1,15));
        name = "jet_top_id_kinRecoVsMva";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; ID jet^{top}_{kinReco}"+label+"; ID jet^{top}_{mva}",16,-1,15,16,-1,15));
        
        name = "jet_top_Pt_trueVsKinReco";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; p_{T} jet^{top}_{true}"+label+"; p_{T} jet^{top}_{kinReco}", nBins_Pt, 0., 400., nBins_Pt, 0., 400.));
        name = "jet_top_Pt_trueVsMva";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; p_{T} jet^{top}_{true}"+label+"; p_{T} jet^{top}_{mva}", nBins_Pt, 0., 400., nBins_Pt, 0., 400.));
        name = "jet_top_Pt_kinRecoVsMva";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; p_{T} jet^{top}_{kinReco}"+label+"; p_{T} jet^{top}_{mva}", nBins_Pt, 0., 400., nBins_Pt, 0., 400.));
        
        name = "jet_top_mult_kinReco_corr";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "; N correct jet^{top}_{kinReco}"+label+"; Events",5,0,5));
        name = "jet_top_mult_mva_corr";
        m_histogram[name] = store(new TH1D(prefix_+name+step, "; N correct jet^{top}_{mva}"+label+"; Events",5,0,5));
        
        

        // For jets that where assigned to tt by both KinReco and MVA
        name = "jet_top_id_sameJet_trueVsKinReco";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; ID jet^{top}_{true} same_{kinReco}^{MVA}"+label+"; ID jet^{top}_{kinReco}",16,-1,15,16,-1,15));
        name = "jet_top_id_sameJet_trueVsMva";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; ID jet^{top}_{true} same_{kinReco}^{MVA}"+label+"; ID jet^{top}_{mva}",16,-1,15,16,-1,15));
        name = "jet_top_id_sameJet_kinRecoVsMva";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; ID jet^{top}_{kinReco} same_{kinReco}^{MVA}"+label+"; ID jet^{top}_{mva}",16,-1,15,16,-1,15));
        
        name = "jet_top_Pt_sameJet_trueVsKinReco";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; p_{T} jet^{top}_{true} same_{kinReco}^{MVA}"+label+"; p_{T} jet^{top}_{kinReco}", nBins_Pt, 0., 400., nBins_Pt, 0., 400.));
        name = "jet_top_Pt_sameJet_trueVsMva";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; p_{T} jet^{top}_{true} same_{kinReco}^{MVA}"+label+"; p_{T} jet^{top}_{mva}", nBins_Pt, 0., 400., nBins_Pt, 0., 400.));
        name = "jet_top_Pt_sameJet_kinRecoVsMva";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; p_{T} jet^{top}_{kinReco} same_{kinReco}^{MVA}"+label+"; p_{T} jet^{top}_{mva}", nBins_Pt, 0., 400., nBins_Pt, 0., 400.));
        

        // For events where both jets were assigned to tt by both KinReco and MVA in the same way
        name = "jet_top_id_same2Jets_trueVsKinReco";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; ID jet^{top}_{true} 2*same_{kinReco}^{MVA}"+label+"; ID jet^{top}_{kinReco}",16,-1,15,16,-1,15));
        name = "jet_top_id_same2Jets_trueVsMva";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; ID jet^{top}_{true} 2*same_{kinReco}^{MVA}"+label+"; ID jet^{top}_{mva}",16,-1,15,16,-1,15));
        name = "jet_top_id_same2Jets_kinRecoVsMva";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; ID jet^{top}_{kinReco} 2*same_{kinReco}^{MVA}"+label+"; ID jet^{top}_{mva}",16,-1,15,16,-1,15));
        
        name = "jet_top_Pt_same2Jets_trueVsKinReco";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; p_{T} jet^{top}_{true} 2*same_{kinReco}^{MVA}"+label+"; p_{T} jet^{top}_{kinReco}", nBins_Pt, 0., 400., nBins_Pt, 0., 400.));
        name = "jet_top_Pt_same2Jets_trueVsMva";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; p_{T} jet^{top}_{true} 2*same_{kinReco}^{MVA}"+label+"; p_{T} jet^{top}_{mva}", nBins_Pt, 0., 400., nBins_Pt, 0., 400.));
        name = "jet_top_Pt_same2Jets_kinRecoVsMva";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; p_{T} jet^{top}_{kinReco} 2*same_{kinReco}^{MVA}"+label+"; p_{T} jet^{top}_{mva}", nBins_Pt, 0., 400., nBins_Pt, 0., 400.));
        
        name = "jet_top_mult_same2Jets_corr";                                                                            
        m_histogram[name] = store(new TH1D(prefix_+name+step, "; N correct jet^{top}_{kinReco} (2*same_{kinReco}^{MVA})"+label+"; Events",5,0,5));
    }

}


void AnalyzerDijet::bookAddGenJetsCorrelationHistos (std::map<TString, TH1*>& m_histogram, const TString addName, const TString& step, 
                                                     const TString& label, const bool bookJetwiseHistos)
{
    const int nBins_Pt_j1 = 5;
    const double bins_Pt_j1[nBins_Pt_j1+1] = {20., 50., 80., 140., 200., 400.};
    const int nBins_Pt_j2 = 5;
    const double bins_Pt_j2[nBins_Pt_j2+1] = {20., 50., 80., 140., 200., 400.};
    const int nBins_Eta_j1 = 8;
    const double bins_Eta_j1[nBins_Eta_j1+1] = {-2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4};
    const int nBins_Eta_j2 = 8;
    const double bins_Eta_j2[nBins_Eta_j2+1] = {-2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4};
    const int nBins_Mjj = 4;
    const double bins_Mjj[nBins_Mjj+1] = {0., 80., 140., 240., 400.};
    const int nBins_dR = 5;
    const double bins_dR[nBins_dR+1] = {0., 1., 2., 3., 4., 5.};
    
    
    TString name;
    name = "leadingJet_1st_Pt_add_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; 1st jet^{add}_{reco} p_{T} "+addName+label+"; p_{T} jet^{add}_{gen}", nBins_Pt_j1, bins_Pt_j1, nBins_Pt_j1, bins_Pt_j1));
    name = "leadingJet_2nd_Pt_add_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; 2nd jet^{add}_{reco} p_{T} "+addName+label+"; p_{T} jet^{add}_{gen}", nBins_Pt_j2, bins_Pt_j2, nBins_Pt_j2, bins_Pt_j2));
    name = "leadingJet_3rd_Pt_add_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; 3rd jet^{add}_{reco} p_{T} "+addName+label+"; p_{T} jet^{add}_{gen}", nBins_Pt_j2, bins_Pt_j2, nBins_Pt_j2, bins_Pt_j2));
    name = "leadingJet_4th_Pt_add_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; 4th jet^{add}_{reco} p_{T} "+addName+label+"; p_{T} jet^{add}_{gen}", nBins_Pt_j2, bins_Pt_j2, nBins_Pt_j2, bins_Pt_j2));
    name = "leadingJet_1st_Eta_add_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; 1st jet^{add}_{reco} #eta "+addName+label+"; #eta jet^{add}_{gen}", nBins_Eta_j1, bins_Eta_j1, nBins_Eta_j1, bins_Eta_j1));
    name = "leadingJet_2nd_Eta_add_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; 2nd jet^{add}_{reco} #eta "+addName+label+"; #eta jet^{add}_{gen}", nBins_Eta_j2, bins_Eta_j2, nBins_Eta_j2, bins_Eta_j2));
    name = "leadingJet_3rd_Eta_add_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; 3rd jet^{add}_{reco} #eta "+addName+label+"; #eta jet^{add}_{gen}", nBins_Eta_j2, bins_Eta_j2, nBins_Eta_j2, bins_Eta_j2));
    name = "leadingJet_4th_Eta_add_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; 4th jet^{add}_{reco} #eta "+addName+label+"; #eta jet^{add}_{gen}", nBins_Eta_j2, bins_Eta_j2, nBins_Eta_j2, bins_Eta_j2));
    name = "leadingJet_dR_add_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; #DeltaR jet^{1,add}_{2,add} reco "+addName+label+"; #DeltaR jet^{1,add}_{2,add} gen",nBins_dR, bins_dR, nBins_dR, bins_dR));
    name = "leadingJet_Mjj_add_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; M(jet^{add}_{"+addName+"}, jet^{add}_{"+addName+"}) "+addName+label+"; M(jet^{add}_{"+addName+"}, jet^{add}_{"+addName+"})", nBins_Mjj, bins_Mjj, nBins_Mjj, bins_Mjj));
    name = "leadingJet_dPhi_add_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; #Delta#phi jet^{1,add}_{2,add} reco "+addName+label+"; #Delta#phi jet^{1,add}_{2,add} gen",10,0.,5.,10,0.,5.));
    
    name = "leadingJet_1st_Pt_addB_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; 1st b-jet^{add}_{reco} p_{T} "+addName+label+"; p_{T} b-jet^{add}_{gen}", nBins_Pt_j1, bins_Pt_j1, nBins_Pt_j1, bins_Pt_j1));
    name = "leadingJet_2nd_Pt_addB_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; 2nd b-jet^{add}_{reco} p_{T} "+addName+label+"; p_{T} b-jet^{add}_{gen}", nBins_Pt_j2, bins_Pt_j2, nBins_Pt_j2, bins_Pt_j2));
    name = "leadingJet_3rd_Pt_addB_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; 3rd b-jet^{add}_{reco} p_{T} "+addName+label+"; p_{T} b-jet^{add}_{gen}", nBins_Pt_j2, bins_Pt_j2, nBins_Pt_j2, bins_Pt_j2));
    name = "leadingJet_4th_Pt_addB_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; 4th b-jet^{add}_{reco} p_{T} "+addName+label+"; p_{T} b-jet^{add}_{gen}", nBins_Pt_j2, bins_Pt_j2, nBins_Pt_j2, bins_Pt_j2));
    name = "leadingJet_1st_Eta_addB_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; 1st b-jet^{add}_{reco} #eta "+addName+label+"; #eta b-jet^{add}_{gen}", nBins_Eta_j1, bins_Eta_j1, nBins_Eta_j1, bins_Eta_j1));
    name = "leadingJet_2nd_Eta_addB_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; 2nd b-jet^{add}_{reco} #eta "+addName+label+"; #eta b-jet^{add}_{gen}", nBins_Eta_j2, bins_Eta_j2, nBins_Eta_j2, bins_Eta_j2));
    name = "leadingJet_3rd_Eta_addB_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; 3rd b-jet^{add}_{reco} #eta "+addName+label+"; #eta b-jet^{add}_{gen}", nBins_Eta_j2, bins_Eta_j2, nBins_Eta_j2, bins_Eta_j2));
    name = "leadingJet_4th_Eta_addB_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; 4th b-jet^{add}_{reco} #eta "+addName+label+"; #eta b-jet^{add}_{gen}", nBins_Eta_j2, bins_Eta_j2, nBins_Eta_j2, bins_Eta_j2));
    name = "leadingJet_dR_addB_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; #DeltaR b-jet^{1,add}_{2,add} reco "+addName+label+"; #DeltaR b-jet^{1,add}_{2,add} gen",nBins_dR,bins_dR,nBins_dR,bins_dR));
    name = "leadingJet_Mjj_addB_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; M(b-jet^{add}_{"+addName+"}, b-jet^{add}_{"+addName+"}) "+addName+label+"; M(b-jet^{add}_{"+addName+"}, b-jet^{add}_{"+addName+"})", nBins_Mjj, bins_Mjj, nBins_Mjj, bins_Mjj));
    name = "leadingJet_dPhi_addB_"+addName+"VsGen";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "; #Delta#phi b-jet^{1,add}_{2,add} reco "+addName+label+"; #Delta#phi b-jet^{1,add}_{2,add} gen",10,0.,5.,10,0.,5.));
    
    if(bookJetwiseHistos) {
        // Leading top jets ordering
        name = "leadingJet_index_top_"+addName+"VsGen";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; ID of recoJet_{"+addName+"}^{match}"+label+"; ID of genJet_{top}; Events",22,-2,20,22,-2,20));
        name = "leadingJet_Pt_top_"+addName+"VsGen";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; Pt recoJet_{"+addName+"}^{match}"+label+"; Pt genJet_{top}; Events", nBins_Pt_j2, bins_Pt_j2, nBins_Pt_j2, bins_Pt_j2));
        name = "leadingJet_Eta_top_"+addName+"VsGen";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; Eta recoJet_{"+addName+"}^{match}"+label+"; Eta genJet_{top}; Events", nBins_Eta_j2, bins_Eta_j2, nBins_Eta_j2, bins_Eta_j2));
        // Leading additional jets ordering
        name = "leadingJet_index_add_"+addName+"VsGen";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; ID of recoJet_{"+addName+"}^{match}"+label+"; ID of genJet_{add}; Events",22,-2,20,22,-2,20));
        name = "leadingJet_Pt_add_"+addName+"VsGen";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; Pt recoJet_{"+addName+"}^{match}"+label+"; Pt genJet_{add}; Events", nBins_Pt_j2, bins_Pt_j2, nBins_Pt_j2, bins_Pt_j2));
        name = "leadingJet_Eta_add_"+addName+"VsGen";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; Eta recoJet_{"+addName+"}^{match}"+label+"; Eta genJet_{add}; Events", nBins_Eta_j2, bins_Eta_j2, nBins_Eta_j2, bins_Eta_j2));
        name = "leadingJet_ttdR_add_"+addName+"VsGen";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; dR(recoJet_{"+addName+"}^{match},t#bar{t}"+label+"; dR (genJet_{add},t#bar{t}; Events",20,0.,5.,20,0.,5.));
        name = "leadingJet_ttdEta_add_"+addName+"VsGen";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; #Delta#eta(recoJet_{"+addName+"}^{match},t#bar{t})"+label+"; #Delta#eta(genJet_{add},t#bar{t}); Events",50,-5.,5.,50,-5.,5.));
        // Leading additional b-jets order
        name = "leadingJet_index_addB_"+addName+"VsGen";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; ID of recoJet_{"+addName+"}^{match}"+label+"; ID of genJet_{addB}; Events",22,-2,20,22,-2,20));
        name = "leadingJet_Pt_addB_"+addName+"VsGen";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; Pt recoJet_{"+addName+"}^{match}"+label+"; Pt genJet_{addB}; Events", nBins_Pt_j2, bins_Pt_j2, nBins_Pt_j2, bins_Pt_j2));
        name = "leadingJet_Eta_addB_"+addName+"VsGen";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; Eta recoJet_{"+addName+"}^{match}"+label+"; Eta genJet_{addB}; Events", nBins_Eta_j2, bins_Eta_j2, nBins_Eta_j2, bins_Eta_j2));
        name = "leadingJet_ttdR_addB_"+addName+"VsGen";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; dR(recoJet_{"+addName+"}^{match},t#bar{t}"+label+"; dR (genJet_{addB},t#bar{t}; Events",20,0.,5.,20,0.,5.));
        name = "leadingJet_ttdEta_addB_"+addName+"VsGen";
        m_histogram[name] = store(new TH2D(prefix_+name+step, "; #Delta#eta(recoJet_{"+addName+"}^{match},t#bar{t})"+label+"; #Delta#eta(genJet_{addB},t#bar{t}); Events",50,-5.,5.,50,-5.,5.));
    }
}


void AnalyzerDijet::bookLeadingJetsHistos (std::map<TString, TH1*>& m_histogram, const TString addName, const TString& step, 
                                           const TString& label)
{
    const int nBins_jet_Pt = 20;
//     const double bins_Pt[nBins_jet_Pt+1] = {0., 20., 60., 120., 180., 240., 300., 360., 420.};
    
    const int nBins_Pt_j1 = 5;
    const double bins_Pt_j1[nBins_Pt_j1+1] = {20., 50., 80., 140., 200., 400.};
    const int nBins_Pt_j2 = 3;
    const double bins_Pt_j2[nBins_Pt_j2+1] = {20., 50., 80., 250.};
    const int nBins_Eta_j1 = 8;
    const double bins_Eta_j1[nBins_Eta_j1+1] = {-2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4};
    const int nBins_Eta_j2 = 4;
    const double bins_Eta_j2[nBins_Eta_j2+1] = {-2.4, -1.2, 0., 1.2, 2.4};
    const int nBins_Mjj = 4;
    const double bins_Mjj[nBins_Mjj+1] = {0., 80., 140., 240., 400.};
    const int nBins_dR = 4;
    const double bins_dR[nBins_dR+1] = {0., 1., 2., 3., 5.};
    
    TString name;
    // Top jets
    name = "leadingJet_1st_Pt_top_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 1st jet^{top}_{"+addName+"} p_{T}"+label+"; Events", nBins_jet_Pt, 0., 400.));
    name = "leadingJet_2nd_Pt_top_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 2nd jet^{top}_{"+addName+"} p_{T}"+label+"; Events", nBins_jet_Pt, 0., 400.));
    name = "leadingJet_1st_Eta_top_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 1st jet^{top}_{"+addName+"} #eta"+label+"; Events",25,-2.5,2.5));
    name = "leadingJet_2nd_Eta_top_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 2nd jet^{top}_{"+addName+"} #eta"+label+"; Events",25,-2.5,2.5));
    name = "leadingJet_dR_top_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; #DeltaR(jet^{top}_{"+addName+"}, jet^{top}_{"+addName+"})"+label+"; Events",10,0.,5.));
    name = "leadingJet_dPhi_top_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; #Delta#phi(jet^{top}_{"+addName+"}, jet^{top}_{"+addName+"})"+label+"; Events",10,0.,5.));
    name = "leadingJet_Mjj_top_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; M(jet^{top}_{"+addName+"}, jet^{top}_{"+addName+"})"+label+"; Events", nBins_jet_Pt, 0., 400.));
    name = "leadingJet_1st_Pt_2nd_Pt_top_";
    m_histogram[name+addName] = store(new TH2D(prefix_+name+addName+step, "; 1st jet^{top}_{"+addName+"} p_{T}"+label+"; 2nd jet^{top}_{"+addName+"} p_{T}", nBins_jet_Pt, 0., 400.,  nBins_jet_Pt, 0., 400.));
    if(!addName.Contains("gen")) {
        name = "leadingJets_btagDiscriminator_top_";
        m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; jet^{top}_{"+addName+"} btag_{D}"+label+"; Events",60,-0.1,1.1));
    }
    // Additional jets
    name = "leadingJet_1st_Pt_add_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 1st jet^{add}_{"+addName+"} p_{T}"+label+"; Events", nBins_jet_Pt, 0., 400.));
    name = "leadingJet_2nd_Pt_add_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 2nd jet^{add}_{"+addName+"} p_{T}"+label+"; Events", nBins_jet_Pt, 0., 400.));
    name = "leadingJet_3rd_Pt_add_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 2nd jet^{add}_{"+addName+"} p_{T}"+label+"; Events", nBins_jet_Pt, 0., 400.));
    name = "leadingJet_4th_Pt_add_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 2nd jet^{add}_{"+addName+"} p_{T}"+label+"; Events", nBins_jet_Pt, 0., 400.));
    name = "leadingJet_1st_Eta_add_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 1st jet^{add}_{"+addName+"} #eta"+label+"; Events",25,-2.5,2.5));
    name = "leadingJet_2nd_Eta_add_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 2nd jet^{add}_{"+addName+"} #eta"+label+"; Events",25,-2.5,2.5));
    name = "leadingJet_3rd_Eta_add_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 3rd jet^{add}_{"+addName+"} #eta"+label+"; Events",25,-2.5,2.5));
    name = "leadingJet_4th_Eta_add_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 4th jet^{add}_{"+addName+"} #eta"+label+"; Events",25,-2.5,2.5));
    name = "leadingJet_dR_add_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; #DeltaR(jet^{1,add}_{"+addName+"}, jet^{2,add}_{"+addName+"}))"+label+"; Events",10,0.,5.));
    name = "leadingJet_dPhi_add_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; #Delta#phi(jet^{1,add}_{"+addName+"}, jet^{2,add}_{"+addName+"}))"+label+"; Events",10,0.,5.));
    name = "leadingJet_Mjj_add_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; M(jet^{1,add}_{"+addName+"}, jet^{2,add}_{"+addName+"}))"+label+"; Events", nBins_jet_Pt, 0., 400.));
    name = "leadingJet_1st_Pt_2nd_Pt_add_";
    m_histogram[name+addName] = store(new TH2D(prefix_+name+addName+step, "; 1st jet^{add}_{"+addName+"} p_{T}"+label+"; 2nd jet^{add}_{"+addName+"} p_{T}", nBins_jet_Pt, 0., 400.,  nBins_jet_Pt, 0., 400.));
    if(!addName.Contains("gen")) {
        name = "leadingJets_btagDiscriminator_add_";
        m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; jet^{add}_{"+addName+"} btag_{D}"+label+"; Events",60,-0.1,1.1));
    }
    // Additional b-jets
    name = "leadingJet_1st_Pt_addB_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 1st b-jet^{add}_{"+addName+"} p_{T}"+label+"; Events", nBins_Pt_j1, bins_Pt_j1));
    name = "leadingJet_2nd_Pt_addB_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 2nd b-jet^{add}_{"+addName+"} p_{T}"+label+"; Events", nBins_Pt_j2, bins_Pt_j2));
    name = "leadingJet_1st_Eta_addB_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 1st b-jet^{add}_{"+addName+"} #eta"+label+"; Events", nBins_Eta_j1, bins_Eta_j1));
    name = "leadingJet_2nd_Eta_addB_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 2nd b-jet^{add}_{"+addName+"} #eta"+label+"; Events", nBins_Eta_j2, bins_Eta_j2));
    name = "leadingJet_dR_addB_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; #DeltaR(b-jet^{1,add}_{"+addName+"}, b-jet^{2,add}_{"+addName+"}))"+label+"; Events", nBins_dR, bins_dR));
    name = "leadingJet_dPhi_addB_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; #Delta#phi(b-jet^{1,add}_{"+addName+"}, b-jet^{2,add}_{"+addName+"}))"+label+"; Events", nBins_dR, bins_dR));
    name = "leadingJet_Mjj_addB_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; M(b-jet^{1,add}_{"+addName+"}, b-jet^{2,add}_{"+addName+"}))"+label+"; Events", nBins_Mjj, bins_Mjj));
    if(!addName.Contains("gen")) {
        name = "leadingJets_btagDiscriminator_addB_";
        m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; b-jet^{add}_{"+addName+"} btag_{D}"+label+"; Events",60,-0.1,1.1));
    }
}


void AnalyzerDijet::bookLeadingJetsCPHistos (std::map<TString, TH1*>& m_histogram, const TString addName, const TString& step, 
                                           const TString& label)
{
    const int nBins_jet_Pt = 19;
    const int nBins_jet_Eta = 13;
    const int nBins_jet_Mjj = 8;

    
    TString name;
    // Top jets
    name = "leadingJet_1st_Pt_top_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 1st jet^{top}_{"+addName+"} p_{T}"+label+"; Events", nBins_jet_Pt, 20., 400.));
    name = "leadingJet_2nd_Pt_top_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 2nd jet^{top}_{"+addName+"} p_{T}"+label+"; Events", nBins_jet_Pt, 20., 400.));
    name = "leadingJet_1st_Eta_top_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 1st jet^{top}_{"+addName+"} #eta"+label+"; Events",nBins_jet_Eta,-2.6,2.6));
    name = "leadingJet_2nd_Eta_top_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 2nd jet^{top}_{"+addName+"} #eta"+label+"; Events",nBins_jet_Eta,-2.6,2.6));
    name = "leadingJet_dR_top_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; #DeltaR(jet^{top}_{"+addName+"}, jet^{top}_{"+addName+"})"+label+"; Events",10,0.,5.));
    name = "leadingJet_dPhi_top_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; #Delta#phi(jet^{top}_{"+addName+"}, jet^{top}_{"+addName+"})"+label+"; Events",10,0.,5.));
    name = "leadingJet_Mjj_top_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; M(jet^{top}_{"+addName+"}, jet^{top}_{"+addName+"})"+label+"; Events", nBins_jet_Mjj, 20., 400.));

    // Additional jets
    name = "leadingJet_1st_Pt_add_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 1st jet^{add}_{"+addName+"} p_{T}"+label+"; Events", nBins_jet_Pt, 20., 400.));
    name = "leadingJet_2nd_Pt_add_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 2nd jet^{add}_{"+addName+"} p_{T}"+label+"; Events", nBins_jet_Pt, 20., 400.));
    name = "leadingJet_3rd_Pt_add_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 2nd jet^{add}_{"+addName+"} p_{T}"+label+"; Events", nBins_jet_Pt, 20., 400.));
    name = "leadingJet_4th_Pt_add_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 2nd jet^{add}_{"+addName+"} p_{T}"+label+"; Events", nBins_jet_Pt, 20., 400.));
    name = "leadingJet_1st_Eta_add_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 1st jet^{add}_{"+addName+"} #eta"+label+"; Events",nBins_jet_Eta,-2.6,2.6));
    name = "leadingJet_2nd_Eta_add_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 2nd jet^{add}_{"+addName+"} #eta"+label+"; Events",nBins_jet_Eta,-2.6,2.6));
    name = "leadingJet_3rd_Eta_add_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 3rd jet^{add}_{"+addName+"} #eta"+label+"; Events",nBins_jet_Eta,-2.6,2.6));
    name = "leadingJet_4th_Eta_add_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 4th jet^{add}_{"+addName+"} #eta"+label+"; Events",nBins_jet_Eta,-2.6,2.6));
    name = "leadingJet_dR_add_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; #DeltaR(jet^{1,add}_{"+addName+"}, jet^{2,add}_{"+addName+"}))"+label+"; Events",10,0.,5.));
    name = "leadingJet_dPhi_add_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; #Delta#phi(jet^{1,add}_{"+addName+"}, jet^{2,add}_{"+addName+"}))"+label+"; Events",10,0.,5.));
    name = "leadingJet_Mjj_add_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; M(jet^{1,add}_{"+addName+"}, jet^{2,add}_{"+addName+"}))"+label+"; Events", nBins_jet_Mjj, 20., 400.));

    // Additional b-jets
    name = "leadingJet_1st_Pt_addB_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 1st b-jet^{add}_{"+addName+"} p_{T}"+label+"; Events", nBins_jet_Pt, 20., 400.));
    name = "leadingJet_2nd_Pt_addB_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 2nd b-jet^{add}_{"+addName+"} p_{T}"+label+"; Events", nBins_jet_Pt, 20., 400.));
    name = "leadingJet_1st_Eta_addB_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 1st b-jet^{add}_{"+addName+"} #eta"+label+"; Events", nBins_jet_Eta,-2.6,2.6));
    name = "leadingJet_2nd_Eta_addB_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; 2nd b-jet^{add}_{"+addName+"} #eta"+label+"; Events", nBins_jet_Eta,-2.6,2.6));
    name = "leadingJet_dR_addB_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; #DeltaR(b-jet^{1,add}_{"+addName+"}, b-jet^{2,add}_{"+addName+"}))"+label+"; Events", 10,0.,5.));
    name = "leadingJet_dPhi_addB_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; #Delta#phi(b-jet^{1,add}_{"+addName+"}, b-jet^{2,add}_{"+addName+"}))"+label+"; Events", 10,0.,5.));
    name = "leadingJet_Mjj_addB_cp_";
    m_histogram[name+addName] = store(new TH1D(prefix_+name+addName+step, "; M(b-jet^{1,add}_{"+addName+"}, b-jet^{2,add}_{"+addName+"}))"+label+"; Events", nBins_jet_Mjj, 20., 400.));

}


void AnalyzerDijet::bookPairHistos(TH1* histo, std::map<TString, TH1*>& m_histogram, const TString& name)
{
    m_histogram[name] = store((TH1D*)histo);
    // Adding histogram for the correct pairs
    TH1D* histo_corrrect = new TH1D(*(TH1D*)histo);
    histo_corrrect->SetName(TString(histo_corrrect->GetName()).ReplaceAll(name, name+"_correct"));
    histo_corrrect->SetTitle(TString(histo_corrrect->GetTitle())+" [correct]");
    m_histogram[name+"_correct"] = store(histo_corrrect);
    // Adding histogram for the wrong pairs
    TH1* histo_wrong = new TH1D(*(TH1D*)histo);
    histo_wrong->SetName(TString(histo_wrong->GetName()).ReplaceAll(name, name+"_wrong"));
    histo_wrong->SetTitle(TString(histo_wrong->GetTitle())+" [wrong]");
    m_histogram[name+"_wrong"] = store(histo_wrong);
}


void AnalyzerDijet::fillPairHistos(std::map<TString, TH1*>& m_histogram, const TString& name, const double value, const bool isCorrect, const double weight)
{
    if(m_histogram[name]) m_histogram[name]->Fill(value, weight);
    
    if(isCorrect && m_histogram[name+"_correct"]) m_histogram[name+"_correct"]->Fill(value, weight);
    else if (m_histogram[name+"_wrong"]) m_histogram[name+"_wrong"]->Fill(value, weight);
}


std::vector<std::pair<int,int> > AnalyzerDijet::jetPairsFromMVA(std::map<TString, TH1*>& m_histogram, const tth::RecoObjectIndices& recoObjectIndices,
                                                                const tth::GenObjectIndices& genObjectIndices, const RecoObjects& recoObjects,
                                                                const std::vector<int>& trueTopJetsId, const std::vector<int>& trueHiggsJetsId, 
                                                                const double weight)
{
    std::vector<std::pair<int,int> > goodJetPairs;
    if(!weightsCorrect_) return goodJetPairs;
//     if(!weightsSwapped_) return goodJetPairs;
    bool requireBJetPairs = true;
    
    // Setting up the MVA input #################################
    // Loop over all jet combinations and get MVA input variables
    std::vector<MvaVariablesBase*> v_mvaVariables = MvaVariablesTopJets::fillVariables(recoObjectIndices, genObjectIndices, recoObjects, weight);

    // Getting the MVA weights from weights file as vector, one entry per jet pair
    std::vector<float> v_mvaWeightsCorrect, v_mvaWeightsSwapped;
    if(weightsCorrect_) v_mvaWeightsCorrect = weightsCorrect_->mvaWeights(v_mvaVariables);
//     if(weightsSwapped_) v_mvaWeightsSwapped = weightsSwapped_->mvaWeights(v_mvaVariables);
    MvaVariablesTopJets::clearVariables(v_mvaVariables);
    
    // Get the indices of the jet pairs and order them by MVA weights, biggest value first
    const tth::IndexPairs& jetIndexPairs = recoObjectIndices.jetIndexPairs_;
//     printf("nAllJets: %d  nJets: %d  nBJets: %d  nPairs: %d\n", (int)recoObjects.jets_->size(), (int)recoObjectIndices.jetIndices_.size(), (int)recoObjectIndices.bjetIndices_.size(), (int)jetIndexPairs.size() );
    
    
    std::vector<int> bJetPairIndices;
    
    // Testing the top/higgs pairs MVA weights
    int topPairId = -1;
    int higgsPairId = -1;
    for(size_t i=0; i<jetIndexPairs.size(); ++i) {
        std::pair<int,int> pair = jetIndexPairs.at(i);
        if(isInVector(recoObjectIndices.bjetIndices_, pair.first) && isInVector(recoObjectIndices.bjetIndices_, pair.second)) bJetPairIndices.push_back(i);
        if(requireBJetPairs && !isInVector(bJetPairIndices, i)) continue;
//         printf("Dijet pair: %d  <%d,%d>  weight: %.3f  ", (int)i, pair.first, pair.second, v_mvaWeightsCorrect.at(i));
        if(isInVector(trueTopJetsId, pair.first) && isInVector(trueTopJetsId, pair.second)) {
            topPairId = i;
            m_histogram["mvaWeight_correctTop_correct"]->Fill(v_mvaWeightsCorrect.at(i), weight);
//             printf("top\n");
        }
        else if (isInVector(trueHiggsJetsId, pair.first) && isInVector(trueHiggsJetsId, pair.second)) {
            higgsPairId = i;
            m_histogram["mvaWeight_correctHiggs_correct"]->Fill(v_mvaWeightsCorrect.at(i), weight);
            m_histogram["mvaWeight_wrongTop_correct"]->Fill(v_mvaWeightsCorrect.at(i), weight);
//             printf("higgs\n");
        }
        else {
            m_histogram["mvaWeight_wrongTopAndHiggs_correct"]->Fill(v_mvaWeightsCorrect.at(i), weight);
            m_histogram["mvaWeight_wrongTop_correct"]->Fill(v_mvaWeightsCorrect.at(i), weight);
//             printf("wrong\n");
        }
    }
    
    double bottomWeight = -0.3;
    // Adding the good jet pairs based on the MVA weight of other pairs
    for(size_t i=0; i<jetIndexPairs.size(); ++i) {
        if(!isInVector(bJetPairIndices, i)) continue;
        std::pair<int,int> pair1 = jetIndexPairs.at(i);
        // Add it to the list of good pairs if there is at least one pair of other jets satisfying the cut
        for(size_t j=0; j<jetIndexPairs.size(); ++j) {
            if(j==i) continue;
            if(requireBJetPairs && !isInVector(bJetPairIndices, j)) continue;
            std::pair<int,int> pair2 = jetIndexPairs.at(j);
            if(pair1.first == pair2.first || pair1.second == pair2.first || pair1.first == pair2.second || pair1.second == pair2.second) continue;
            double weight2 = v_mvaWeightsCorrect.at(j);
            if(weight2<bottomWeight) continue;
            goodJetPairs.push_back(pair1);
            break;
        }
    }
    
    if(topPairId>=0 && higgsPairId>=0) {
        ((TH2D*)m_histogram["mvaWeight_correctTop_correctHiggs"])->Fill(v_mvaWeightsCorrect.at(topPairId), v_mvaWeightsCorrect.at(higgsPairId), weight);
    }
    
    return goodJetPairs;
}


void AnalyzerDijet::checkAdditionalGenBJetAcceptance(const TopGenObjects& topGenObjects, const tth::GenObjectIndices& genObjectIndices, 
                                                     std::map<TString, TH1*>& m_histogram, const double weight) 
{
    // Setting variables of gen. level if available
    std::vector<int> genJetIndices = genObjectIndices.genJetIndices_;
    const std::vector<int>& bHadJetIndex = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadJetIndex_ : std::vector<int>(0);
    std::vector<int> genBJetsId = genObjectIndices.genBjetIndices_;
    std::vector<std::vector<int> > genJetBhadronIndices = genObjectIndices.genJetBhadronIndices_;
    const std::vector<int>& bHadFlavour = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadFlavour_ : std::vector<int>(0);
    const std::vector<int>& bHadFromTopWeakDecay = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadFromTopWeakDecay_ : std::vector<int>(0);
    
    std::vector<int> addHadronIds;
    std::vector<int> addJetIds;
    
    std::vector<int> topJetIds;
    // Creating the list of gen b-jets coming from  or not from top
    for(size_t hadId = 0; hadId < bHadFlavour.size(); ++hadId) {
        if(bHadJetIndex.at(hadId)<0) continue;
        if(!isInVector(genJetIndices, bHadJetIndex.at(hadId))) continue;
        if(std::abs(bHadFlavour.at(hadId))==6) {
            putUniquelyInVector(topJetIds, bHadJetIndex.at(hadId));
        } else {
            if(bHadFromTopWeakDecay.at(hadId)==1) continue;
            putUniquelyInVector(addJetIds, bHadJetIndex.at(hadId));
        }
    }
    
    for(size_t hadId = 0; hadId < bHadFlavour.size(); ++hadId) {
        int flavour = std::abs(bHadFlavour.at(hadId));
        if(flavour==6) continue;
        if(bHadFromTopWeakDecay.at(hadId)==1) continue;
        putUniquelyInVector(addHadronIds, hadId);
        m_histogram["extraGenBHadronsJets"]->Fill(0., weight);
        // Hadron not clustered to any jet
        if(bHadJetIndex.at(hadId)<0) continue;
        m_histogram["extraGenBHadronsJets"]->Fill(1., weight);
        int jetId = bHadJetIndex.at(hadId);
        // Jet not in acceptance
        if(!isInVector(genJetIndices, jetId)) continue;
        m_histogram["extraGenBHadronsJets"]->Fill(2., weight);
        // Jet overlaps with b-jets from tt
        if(isInVector(topJetIds, jetId)) continue;
        m_histogram["extraGenBHadronsJets"]->Fill(3., weight);
        // Jet overlaps with extra b-jets
        if(genJetBhadronIndices.at(jetId).size()>1) continue;
        m_histogram["extraGenBHadronsJets"]->Fill(4., weight);
    }       // End of loop over all b-hadrons
    
    m_histogram["extraGenBHadronMultiplicity"]->Fill((int)addHadronIds.size(), weight);
    m_histogram["extraGenBJetMultiplicity"]->Fill((int)addJetIds.size(), weight);
}


void AnalyzerDijet::fillDijetMassForPairs(const VLV& allJets, const std::vector<int>& higgsJetsId,
                                          const std::vector<std::pair<int, int> > &jetPairs, const double weight, 
                                          std::map<TString, TH1*>& m_histogram, std::string histoName, const bool normaliseWeight )
{
    unsigned int nPairs = jetPairs.size();
    int nEntries = 0;
    float diffMin = 999.9;
    TH1D* h_nEntriesPerBin = 0;
    int nBins = m_histogram[histoName]->GetNbinsX();
    h_nEntriesPerBin = new TH1D("h_nEntriesPerBin_temp","",nBins,m_histogram[histoName]->GetXaxis()->GetXbins()->GetArray());
    
    for(unsigned int i = 0; i<nPairs; ++i) {
        const std::pair<int,int> &allJetPair = jetPairs.at(i);
        bool isCorrect = isInVector(higgsJetsId, allJetPair.first) && isInVector(higgsJetsId, allJetPair.second);
        
        LV dijet = allJets.at(allJetPair.first) + allJets.at(allJetPair.second);
        if(!normaliseWeight) fillPairHistos(m_histogram, histoName, dijet.M(), isCorrect, weight);
        if(h_nEntriesPerBin)h_nEntriesPerBin->Fill(dijet.M());
        nEntries++;
        for(unsigned int j = i+1; j<nPairs; ++j) {
            std::pair<int,int> allJetPair = jetPairs.at(j);
            
            LV dijet2 = allJets.at(allJetPair.first) + allJets.at(allJetPair.second);
            
            float diff = std::fabs(dijet.M() - dijet2.M());
            if(diff<diffMin) diffMin = diff;
        }
    }
    
    if(normaliseWeight) {
        // Looping through bins of the dijet.M plot and counting entries to fill 1 entry/bin
        for(int iBin=0; iBin<nBins; ++iBin) {
            int nPairs = h_nEntriesPerBin->GetBinContent(iBin+1);
            if(nPairs<1) continue;
            double fillWeight = weight*double(nPairs);
            fillPairHistos(m_histogram, histoName, h_nEntriesPerBin->GetBinLowEdge(iBin+1), false, fillWeight);
        }
    }
    
    
    
    if(m_histogram[histoName+"_nEntriesPerEvent"]) m_histogram[histoName+"_nEntriesPerEvent"]->Fill(nEntries, weight);
    if(m_histogram[histoName+"_diffMin"] && diffMin < 999.) m_histogram[histoName+"_diffMin"]->Fill(diffMin, weight);
    
    if(m_histogram[histoName+"_nEntriesPerBin"]) {
        for(int iBin=0; iBin<nBins; ++iBin) {
            if(iBin >= m_histogram[histoName+"_nEntriesPerBin"]->GetNbinsX()) break;
            int N = h_nEntriesPerBin->GetBinContent(iBin+1);
            if(N<1) continue;
//             printf("reading bin: %d  N: %d Entries: %d\n", iBin, N, int(h_nEntriesPerBin->GetEntries()));
           ((TH2D*) m_histogram[histoName+"_nEntriesPerBin"])->Fill(h_nEntriesPerBin->GetBinLowEdge(iBin+1), N, weight);
        }
    }
    delete h_nEntriesPerBin;
}


float AnalyzerDijet::correctPairFraction(const VLV& allJets, const std::vector<int>& jetsId,
                                         const std::vector<int>& bJetsId, const std::vector<double>& jetsBtagDiscriminant,
                                         const std::vector<int>& topJetsId, const std::vector<int>& higgsJetsId,
                                         const double weight, std::map<TString, TH1*>& m_histogram, std::string histoName, bool fillAllCombinations,
                                         const double jetPt_threshold, const int lowerHigher, const std::vector<std::pair<int, int> > &pairsToIgnore)
{
    unsigned int nJets = jetsId.size();
    int nCorrectJetPairs = 0;
    int nAllJetPairs = 0;
    
    std::string histoName_dijetM("dijet_mass_");
    histoName_dijetM.append(histoName);
    
    std::string histoName_corJetMult("dijet_correctJet_multiplicity");
    histoName_corJetMult.append(histoName);


    std::vector<int> higgsJetCandidatesId(0);

    for(size_t iJet1 = 0; iJet1<nJets; iJet1++) {
        const int iAllJet1 = jetsId.at(iJet1);
        if(isInVector(topJetsId, iAllJet1)) continue;              // Skip jet if it was assigned to the Top
        if(!isInVector(bJetsId, iAllJet1)) continue;               // Skip jet if it is not b-tagged
        if(lowerHigher==1 && allJets.at(iAllJet1).Pt()<jetPt_threshold) continue;       // Skip jets with low Pt

        for(size_t iJet2 = 0; iJet2<iJet1; iJet2++) {
            const int iAllJet2 = jetsId.at(iJet2);
            if(isInVector(topJetsId, iAllJet2)) continue;          // Skip jet if it was assigned to the Top
            if(!isInVector(bJetsId, iAllJet2)) continue;           // Skip jet if it is not b-tagged

            if(areAmongPairs(pairsToIgnore, iJet1, iJet2)) continue;      // Skip jet pairs that should be ignored


            if(lowerHigher==1 && allJets.at(iAllJet2).Pt()<jetPt_threshold) continue;       // 2 jets are high Pt
            if(lowerHigher==-1 && allJets.at(iAllJet1).Pt()>=jetPt_threshold && allJets.at(iAllJet2).Pt()>=jetPt_threshold) continue;  // >=1 jet is low Pt

            int nCorrectJets = 0;
            if(fillAllCombinations) {
                // Filling the dijet mass
                LV dijet = allJets.at(iAllJet1) + allJets.at(iAllJet2);
                fillPairHistos(m_histogram, histoName_dijetM, dijet.M(), (isInVector(higgsJetsId, iAllJet1) && isInVector(higgsJetsId, iAllJet2)), weight);
//                 if(m_histogram[histoName_dijetM]) m_histogram[histoName_dijetM]->Fill(dijet.M(), weight);


                if(isInVector(higgsJetsId, iAllJet1)) nCorrectJets++;
                if(isInVector(higgsJetsId, iAllJet2)) nCorrectJets++;
                // Filling the number of correct jets in the dijet pair
                if(m_histogram[histoName_corJetMult]) m_histogram[histoName_corJetMult]->Fill(nCorrectJets, weight);
            }       // If fill all combinations of jets

            // Adding the 2 jets to the list of candidates for b-jets from Higgs
            if(!isInVector(higgsJetCandidatesId,iAllJet1)) higgsJetCandidatesId.push_back(iAllJet1);
            if(!isInVector(higgsJetCandidatesId,iAllJet2)) higgsJetCandidatesId.push_back(iAllJet2);

            // Updating the number of correct/wrong pairs
            if(nCorrectJets>=2) nCorrectJetPairs++;
            nAllJetPairs++;

        }       // End of the second loop over jets
    }       // End of the first loop over jets

    float correctPairFraction;
    if(higgsJetCandidatesId.size()<2) {
        correctPairFraction = -0.19;                   // If less than 2 true Higgs jets in acceptance
        return correctPairFraction;
    }

    correctPairFraction = (float)nCorrectJetPairs/nAllJetPairs;
//     if(bJetsId.size()>3 && topJetsId.size()==0) printf("All: %d\tCor: %d\tFrac: %.2f\n",nAllJetPairs,nCorrectJetPairs,correctPairFraction);

    if(higgsJetsId.size()<2) correctPairFraction = -0.09;          // If less than 2 candidates were found

    // Stopping if best combination is not needed
    if(fillAllCombinations) return correctPairFraction;

    // Finding the two jets that are assumed to come from the Higgs
    // Ordering jets by b-tagging discriminant if there are more than 2 of them
    if(higgsJetCandidatesId.size()>2) common::orderIndices(higgsJetCandidatesId, jetsBtagDiscriminant);
    LV dijet = allJets.at(higgsJetCandidatesId.at(0)) + allJets.at(higgsJetCandidatesId.at(1));
    fillPairHistos(m_histogram, histoName_dijetM, dijet.M(), (isInVector(higgsJetsId, higgsJetCandidatesId.at(0)) && isInVector(higgsJetsId, higgsJetCandidatesId.at(1))), weight);


    // Counting the number of correct higgs jets among the candidates
    int nCorrectJets = 0;
    for(size_t iJet=0; iJet<higgsJetsId.size(); iJet++) {
        int allJetId = higgsJetsId.at(iJet);
        if(!isInVector(higgsJetCandidatesId, allJetId)) continue;
        nCorrectJets++;
    }
    if(m_histogram[histoName_corJetMult]) m_histogram[histoName_corJetMult]->Fill(nCorrectJets, weight);

    if(higgsJetsId.size()>=2)  {
        correctPairFraction = (int)(isInVector(higgsJetsId, higgsJetCandidatesId.at(0)) && isInVector(higgsJetsId, higgsJetCandidatesId.at(1)));
    }

    return correctPairFraction;
}


void AnalyzerDijet::fillGenRecoMatchingComparisonHistos(const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                                                        const VLV& bHadLVs, const std::vector<int>& bHadFlavour, const std::vector<int>& bHadJetIndex,
                                                        const VLV& genJets, std::map<TString, TH1*>& m_histogram, const double weight )
{
    // Finding the closest b-hadrons to the b-quarks
    LV* bQt = topGenObjects.GenB_;
    LV* bQat = topGenObjects.GenAntiB_;
    LV* bQh = higgsGenObjects.GenBFromH_;
    LV* bQah = higgsGenObjects.GenAntiBFromH_;
    int bJt_id_matched = -1;
    int bJat_id_matched = -1;
    int bJh_id_matched = -1;
    int bJah_id_matched = -1;
    int bHt_id=-1;
    int bHat_id=-1;
    int bHh_id=-1;
    int bHah_id=-1;
    float dRQHt_min = 999.9;
    float dRQHat_min = 999.9;
    float dRQHh_min = 999.9;
    float dRQHah_min = 999.9;
    int bHt_id_05=-1;
    int bHat_id_05=-1;
    int bHh_id_05=-1;
    int bHah_id_05=-1;
    float dRQHt_min_05 = 0.5;
    float dRQHat_min_05 = 0.5;
    float dRQHh_min_05 = 0.5;
    float dRQHah_min_05 = 0.5;
    std::vector<int> bHt_unique_ids_dR(0);
    std::vector<int> bHh_unique_ids_dR(0);
    std::vector<int> bHth_unique_ids_dR(0);
    std::vector<int> bHt_unique_ids_dR_05(0);
    std::vector<int> bHh_unique_ids_dR_05(0);
    std::vector<int> bHth_unique_ids_dR_05(0);
    std::vector<int> bHt_unique_ids_matched(0);
    std::vector<int> bHh_unique_ids_matched(0);
    std::vector<int> bHth_unique_ids_matched(0);
    std::vector<int> bJt_unique_ids_dR(0);
    std::vector<int> bJh_unique_ids_dR(0);
    std::vector<int> bJth_unique_ids_dR(0);
    std::vector<int> bJt_unique_ids_dR_05(0);
    std::vector<int> bJh_unique_ids_dR_05(0);
    std::vector<int> bJth_unique_ids_dR_05(0);
    std::vector<int> bJt_unique_ids_matched(0);
    std::vector<int> bJh_unique_ids_matched(0);
    std::vector<int> bJth_unique_ids_matched(0);

    // Looping over all b-hadrons in the event to find closest ones
    using ROOT::Math::VectorUtil::DeltaR;
    for(size_t iHad = 0; iHad < bHadLVs.size(); iHad++ ) {
        LV bHadLV = bHadLVs.at(iHad);
        float dRQHt = bQt?DeltaR(*bQt, bHadLV):1000.f;
        float dRQHat = bQat?DeltaR(*bQat, bHadLV):1000.f;
        float dRQHh = bQh?DeltaR(*bQh, bHadLV):1000.f;
        float dRQHah = bQah?DeltaR(*bQah, bHadLV):1000.f;
        if(dRQHt < dRQHt_min)   {dRQHt_min = dRQHt; bHt_id = iHad;}
        if(dRQHat < dRQHat_min) {dRQHat_min = dRQHat; bHat_id = iHad;}
        if(dRQHh < dRQHh_min)   {dRQHh_min = dRQHh; bHh_id = iHad;}
        if(dRQHah < dRQHah_min) {dRQHah_min = dRQHah; bHah_id = iHad;}
        if(dRQHt < dRQHt_min_05)   {dRQHt_min_05 = dRQHt; bHt_id_05 = iHad;}
        if(dRQHat < dRQHat_min_05) {dRQHat_min_05 = dRQHat; bHat_id_05 = iHad;}
        if(dRQHh < dRQHh_min_05)   {dRQHh_min_05 = dRQHh; bHh_id_05 = iHad;}
        if(dRQHah < dRQHah_min_05) {dRQHah_min_05 = dRQHah; bHah_id_05 = iHad;}

        // Filling multiplicity of unique hadrons and dR to closest quark
        switch(bHadFlavour.at(iHad)) {
            case 6:
                bHt_unique_ids_matched.push_back(iHad);
                if(bQt) m_histogram["bQuarkT_bHad_dRmatched"]->Fill(DeltaR(*bQt, bHadLV), weight);
                if(!isInVector(bHth_unique_ids_matched, iHad)) bHth_unique_ids_matched.push_back(iHad);
                break;
            case -6:
                bHt_unique_ids_matched.push_back(iHad);
                if(bQat) m_histogram["bQuarkT_bHad_dRmatched"]->Fill(DeltaR(*bQat, bHadLV), weight);
                if(!isInVector(bHth_unique_ids_matched, iHad)) bHth_unique_ids_matched.push_back(iHad);
                break;
            case 25:
                bHh_unique_ids_matched.push_back(iHad);
                if(bQh) m_histogram["bQuarkH_bHad_dRmatched"]->Fill(DeltaR(*bQh, bHadLV), weight);
                if(!isInVector(bHth_unique_ids_matched, iHad)) bHth_unique_ids_matched.push_back(iHad);
                break;
            case -25:
                bHh_unique_ids_matched.push_back(iHad);
                if(bQah) m_histogram["bQuarkH_bHad_dRmatched"]->Fill(DeltaR(*bQah, bHadLV), weight);
                if(!isInVector(bHth_unique_ids_matched, iHad)) bHth_unique_ids_matched.push_back(iHad);
                break;
            default:
                break;
        }

        int flavour = bHadFlavour.at(iHad);
        int flavourAbs = std::abs(flavour);

        if(flavourAbs!=6 && flavourAbs!=25 && flavourAbs!=21) continue;

        int genJetId_dR = genJetIdOfRecoJet(bHadLV, genJets);
        int genJetId_matched = bHadJetIndex.at(iHad);
        if(genJetId_matched<0) continue;

        float hadJetdR_matched = DeltaR(bHadLV, genJets.at(genJetId_matched));
        float hadJetdR_dR = DeltaR(bHadLV, genJets.at(genJetId_dR));


        // Filling dR to matched and closest jets
        switch(flavourAbs) {
            case 6:
                m_histogram["bHadT_genJet_dRmatched"]->Fill(hadJetdR_matched, weight);
                putUniquelyInVector(bJt_unique_ids_matched, genJetId_matched);
                putUniquelyInVector(bJth_unique_ids_matched, genJetId_matched);
                if(genJetId_dR>=0) m_histogram["bHadT_genJet_dRmin_matched"]->Fill(hadJetdR_dR, weight);
                if(flavour<0 && bQat) {
                    m_histogram["bQuarkT_genJet_dRmatched"]->Fill(DeltaR(*bQat, genJets.at(genJetId_matched)));
                    if(DeltaR(*bQat, genJets.at(genJetId_matched)) < 1.0) m_histogram["bQuarkTpt_genJet_dRmatched_low"]->Fill(bQat->Pt()); else
                        m_histogram["bQuarkTpt_genJet_dRmatched_high"]->Fill(bQat->Pt());
                    bJat_id_matched = genJetId_matched;
                }
                else if(flavour>0 && bQt) {
                    m_histogram["bQuarkT_genJet_dRmatched"]->Fill(DeltaR(*bQt, genJets.at(genJetId_matched)));
                    if(DeltaR(*bQt, genJets.at(genJetId_matched)) < 1.0) m_histogram["bQuarkTpt_genJet_dRmatched_low"]->Fill(bQt->Pt()); else
                        m_histogram["bQuarkTpt_genJet_dRmatched_high"]->Fill(bQt->Pt());
                    bJt_id_matched = genJetId_matched;
                }
                break;
            case 25:
                m_histogram["bHadH_genJet_dRmatched"]->Fill(hadJetdR_matched, weight);
                putUniquelyInVector(bJh_unique_ids_matched, genJetId_matched);
                putUniquelyInVector(bJth_unique_ids_matched, genJetId_matched);
                if(genJetId_dR>=0) m_histogram["bHadH_genJet_dRmin_matched"]->Fill(hadJetdR_dR, weight);
                if(flavour<0 && bQah) {
                    m_histogram["bQuarkH_genJet_dRmatched"]->Fill(DeltaR(*bQah, genJets.at(genJetId_matched)));
                    if(DeltaR(*bQah, genJets.at(genJetId_matched)) < 1.0) m_histogram["bQuarkHpt_genJet_dRmatched_low"]->Fill(bQah->Pt()); else
                        m_histogram["bQuarkHpt_genJet_dRmatched_high"]->Fill(bQah->Pt());
                    bJah_id_matched = genJetId_matched;
                }
                else if(flavour>0 && bQh) {
                    m_histogram["bQuarkH_genJet_dRmatched"]->Fill(DeltaR(*bQh, genJets.at(genJetId_matched)));
                    if(DeltaR(*bQh, genJets.at(genJetId_matched)) < 1.0) m_histogram["bQuarkHpt_genJet_dRmatched_low"]->Fill(bQh->Pt()); else
                        m_histogram["bQuarkHpt_genJet_dRmatched_high"]->Fill(bQh->Pt());
                    bJh_id_matched = genJetId_matched;
                }
                break;
            case 21:
                if(genJetId_dR>=0) m_histogram["bHadG_genJet_dRmin_matched"]->Fill(hadJetdR_dR, weight);
                m_histogram["bHadG_genJet_dRmatched"]->Fill(hadJetdR_matched, weight);
                break;
            default:
                break;
        }
    }       // End of loop over b-hadrons


    if(bHt_id >= 0) {
        int hadId = bHt_id;
        if(bQt) m_histogram["bQuarkT_bHad_dRmin"]->Fill(DeltaR(*bQt, bHadLVs.at(hadId)), weight);
        putUniquelyInVector(bHt_unique_ids_dR, hadId);
        putUniquelyInVector(bHth_unique_ids_dR, hadId);
        int genJetId = genJetIdOfRecoJet(bHadLVs.at(hadId), genJets);
        if(genJetId>=0) {
            putUniquelyInVector(bJt_unique_ids_dR, genJetId);
            putUniquelyInVector(bJth_unique_ids_dR, genJetId);
            m_histogram["bHadT_genJet_dRmin"]->Fill(DeltaR(bHadLVs.at(hadId), genJets.at(genJetId)), weight);
            LV* bQ = bQt;
            if(bQ) m_histogram["bQuarkT_genJet_dRmin"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
        }
    }
    if(bHat_id >= 0) {
        int hadId = bHat_id;
        if(bQat) m_histogram["bQuarkT_bHad_dRmin"]->Fill(DeltaR(*bQat, bHadLVs.at(hadId)), weight);
        putUniquelyInVector(bHt_unique_ids_dR, hadId);
        putUniquelyInVector(bHth_unique_ids_dR, hadId);
        int genJetId = genJetIdOfRecoJet(bHadLVs.at(hadId), genJets);
        if(genJetId>=0) {
            putUniquelyInVector(bJt_unique_ids_dR, genJetId);
            putUniquelyInVector(bJth_unique_ids_dR, genJetId);
            m_histogram["bHadT_genJet_dRmin"]->Fill(DeltaR(bHadLVs.at(hadId), genJets.at(genJetId)), weight);
            LV* bQ = bQat;
            if(bQ) m_histogram["bQuarkT_genJet_dRmin"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
        }
    }
    if(bHh_id >= 0) {
        int hadId = bHh_id;
        if(bQh) m_histogram["bQuarkH_bHad_dRmin"]->Fill(DeltaR(*bQh, bHadLVs.at(hadId)), weight);
        putUniquelyInVector(bHh_unique_ids_dR, hadId);
        putUniquelyInVector(bHth_unique_ids_dR, hadId);
        int genJetId = genJetIdOfRecoJet(bHadLVs.at(hadId), genJets);
        if(genJetId>=0) {
            putUniquelyInVector(bJh_unique_ids_dR, genJetId);
            putUniquelyInVector(bJth_unique_ids_dR, genJetId);
            m_histogram["bHadT_genJet_dRmin"]->Fill(DeltaR(bHadLVs.at(hadId), genJets.at(genJetId)), weight);
            LV* bQ = bQh;
            if(bQ) m_histogram["bQuarkT_genJet_dRmin"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
        }
    }
    if(bHah_id >= 0) {
        int hadId = bHah_id;
        if(bQah) m_histogram["bQuarkH_bHad_dRmin"]->Fill(DeltaR(*bQah, bHadLVs.at(hadId)), weight);
        putUniquelyInVector(bHh_unique_ids_dR, hadId);
        putUniquelyInVector(bHth_unique_ids_dR, hadId);
        int genJetId = genJetIdOfRecoJet(bHadLVs.at(hadId), genJets);
        if(genJetId>=0) {
            putUniquelyInVector(bJh_unique_ids_dR, genJetId);
            putUniquelyInVector(bJth_unique_ids_dR, genJetId);
            m_histogram["bHadT_genJet_dRmin"]->Fill(DeltaR(bHadLVs.at(hadId), genJets.at(genJetId)), weight);
            LV* bQ = bQah;
            if(bQ) m_histogram["bQuarkT_genJet_dRmin"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
        }
    }

    if(bHt_id_05 >= 0) {
        int hadId = bHt_id_05;
        putUniquelyInVector(bHt_unique_ids_dR_05, bHt_id_05);
        putUniquelyInVector(bHth_unique_ids_dR_05, bHt_id_05);
        int genJetId = genJetIdOfRecoJet(bHadLVs.at(bHt_id_05), genJets, 0.5);
        if(genJetId>=0) {
            putUniquelyInVector(bJt_unique_ids_dR_05, genJetId);
            putUniquelyInVector(bJth_unique_ids_dR_05, genJetId);
            m_histogram["bHadT_genJet_dRmin_05"]->Fill(DeltaR(bHadLVs.at(hadId), genJets.at(genJetId)), weight);
            LV* bQ = bQt;
            if(bQ) m_histogram["bQuarkT_genJet_dRmin_05"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
            if(genJetId!=bJt_id_matched && bJt_id_matched>=0) {
                m_histogram["bQuarkT_genJet_dRmin_05_notAsMatched"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
                m_histogram["bQuarkTH_genJet_dRmin_05_notAsMatched"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
                m_histogram["bQuarkT_genJet_dRmin_05_isAsMatched"]->Fill(0.0, weight);
                m_histogram["bQuarkTH_genJet_dRmin_05_isAsMatched"]->Fill(0.0, weight);
            } else {
                m_histogram["bQuarkT_genJet_dRmin_05_isAsMatched"]->Fill(1.0, weight);
                m_histogram["bQuarkTH_genJet_dRmin_05_isAsMatched"]->Fill(1.0, weight);
            }
        }
    }
    if(bHat_id_05 >= 0) {
        int hadId = bHat_id_05;
        putUniquelyInVector(bHt_unique_ids_dR_05, bHat_id_05);
        putUniquelyInVector(bHth_unique_ids_dR_05, bHat_id_05);
        int genJetId = genJetIdOfRecoJet(bHadLVs.at(bHat_id_05), genJets, 0.5);
        if(genJetId>=0) {
            putUniquelyInVector(bJt_unique_ids_dR_05, genJetId);
            putUniquelyInVector(bJth_unique_ids_dR_05, genJetId);
            m_histogram["bHadT_genJet_dRmin_05"]->Fill(DeltaR(bHadLVs.at(hadId), genJets.at(genJetId)), weight);
            LV* bQ = bQat;
            if(bQ) m_histogram["bQuarkT_genJet_dRmin_05"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
            if(genJetId!=bJat_id_matched && bJat_id_matched>=0) {
                m_histogram["bQuarkT_genJet_dRmin_05_notAsMatched"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
                m_histogram["bQuarkTH_genJet_dRmin_05_notAsMatched"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
                m_histogram["bQuarkT_genJet_dRmin_05_isAsMatched"]->Fill(0.0, weight);
                m_histogram["bQuarkTH_genJet_dRmin_05_isAsMatched"]->Fill(0.0, weight);
            } else {
                m_histogram["bQuarkT_genJet_dRmin_05_isAsMatched"]->Fill(1.0, weight);
                m_histogram["bQuarkTH_genJet_dRmin_05_isAsMatched"]->Fill(1.0, weight);
            }
        }
    }
    if(bHh_id_05 >= 0) {
        int hadId = bHh_id_05;
        putUniquelyInVector(bHh_unique_ids_dR_05, bHh_id_05);
        putUniquelyInVector(bHth_unique_ids_dR_05, bHh_id_05);
        int genJetId = genJetIdOfRecoJet(bHadLVs.at(bHh_id_05), genJets, 0.5);
        if(genJetId>=0) {
            putUniquelyInVector(bJt_unique_ids_dR_05, genJetId);
            putUniquelyInVector(bJth_unique_ids_dR_05, genJetId);
            m_histogram["bHadH_genJet_dRmin_05"]->Fill(DeltaR(bHadLVs.at(hadId), genJets.at(genJetId)), weight);
            LV* bQ = bQh;
            if(bQ) m_histogram["bQuarkH_genJet_dRmin_05"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
            if(genJetId!=bJh_id_matched && bJh_id_matched>=0) {
                m_histogram["bQuarkH_genJet_dRmin_05_notAsMatched"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
                m_histogram["bQuarkTH_genJet_dRmin_05_notAsMatched"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
                m_histogram["bQuarkH_genJet_dRmin_05_isAsMatched"]->Fill(0.0, weight);
                m_histogram["bQuarkTH_genJet_dRmin_05_isAsMatched"]->Fill(0.0, weight);
            } else {
                m_histogram["bQuarkH_genJet_dRmin_05_isAsMatched"]->Fill(1.0, weight);
                m_histogram["bQuarkTH_genJet_dRmin_05_isAsMatched"]->Fill(1.0, weight);
            }
        }
    }
    if(bHah_id_05 >= 0) {
        int hadId = bHah_id_05;
        putUniquelyInVector(bHh_unique_ids_dR_05, bHah_id_05);
        putUniquelyInVector(bHth_unique_ids_dR_05, bHah_id_05);
        int genJetId = genJetIdOfRecoJet(bHadLVs.at(bHah_id_05), genJets, 0.5);
        if(genJetId>=0) {
            putUniquelyInVector(bJt_unique_ids_dR_05, genJetId);
            putUniquelyInVector(bJth_unique_ids_dR_05, genJetId);
            m_histogram["bHadH_genJet_dRmin_05"]->Fill(DeltaR(bHadLVs.at(hadId), genJets.at(genJetId)), weight);
            LV* bQ = bQah;
            if(bQ) m_histogram["bQuarkH_genJet_dRmin_05"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
            if(genJetId!=bJah_id_matched && bJah_id_matched>=0) {
                m_histogram["bQuarkH_genJet_dRmin_05_notAsMatched"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
                m_histogram["bQuarkTH_genJet_dRmin_05_notAsMatched"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
                m_histogram["bQuarkH_genJet_dRmin_05_isAsMatched"]->Fill(0.0, weight);
                m_histogram["bQuarkTH_genJet_dRmin_05_isAsMatched"]->Fill(0.0, weight);
            } else {
                m_histogram["bQuarkH_genJet_dRmin_05_isAsMatched"]->Fill(1.0, weight);
                m_histogram["bQuarkTH_genJet_dRmin_05_isAsMatched"]->Fill(1.0, weight);
            }
        }
    }

    // Unique hadrons
    m_histogram["bHadT_unique_multiplicity_dR"]->Fill(bHt_unique_ids_dR.size());
    m_histogram["bHadH_unique_multiplicity_dR"]->Fill(bHh_unique_ids_dR.size());
    m_histogram["bHadTH_unique_multiplicity_dR"]->Fill(bHth_unique_ids_dR.size());

    m_histogram["bHadT_unique_multiplicity_dR_05"]->Fill(bHt_unique_ids_dR_05.size());
    m_histogram["bHadH_unique_multiplicity_dR_05"]->Fill(bHh_unique_ids_dR_05.size());
    m_histogram["bHadTH_unique_multiplicity_dR_05"]->Fill(bHth_unique_ids_dR_05.size());

    m_histogram["bHadT_unique_multiplicity_match"]->Fill(bHt_unique_ids_matched.size());
    m_histogram["bHadH_unique_multiplicity_match"]->Fill(bHh_unique_ids_matched.size());
    m_histogram["bHadTH_unique_multiplicity_match"]->Fill(bHth_unique_ids_matched.size());

    // Unique jets
    m_histogram["genJetT_unique_multiplicity_dR"]->Fill(bJt_unique_ids_dR.size());
    m_histogram["genJetH_unique_multiplicity_dR"]->Fill(bJh_unique_ids_dR.size());
    m_histogram["genJetTH_unique_multiplicity_dR"]->Fill(bJth_unique_ids_dR.size());

    m_histogram["genJetT_unique_multiplicity_dR_05"]->Fill(bJt_unique_ids_dR_05.size());
    m_histogram["genJetH_unique_multiplicity_dR_05"]->Fill(bJh_unique_ids_dR_05.size());
    m_histogram["genJetTH_unique_multiplicity_dR_05"]->Fill(bJth_unique_ids_dR_05.size());

    m_histogram["genJetT_unique_multiplicity_match"]->Fill(bJt_unique_ids_matched.size());
    m_histogram["genJetH_unique_multiplicity_match"]->Fill(bJh_unique_ids_matched.size());
    m_histogram["genJetTH_unique_multiplicity_match"]->Fill(bJth_unique_ids_matched.size());


    std::vector<int> bJt_unique_ids_directQ_dR;
    std::vector<int> bJh_unique_ids_directQ_dR;
    std::vector<int> bJth_unique_ids_directQ_dR;
    std::vector<int> bJt_unique_ids_directQ_dR_05;
    std::vector<int> bJh_unique_ids_directQ_dR_05;
    std::vector<int> bJth_unique_ids_directQ_dR_05;
    // Direct quark-jet matching
    if(bQt) {
        LV* bQ = bQt;
        int genJetId = genJetIdOfRecoJet(*bQ, genJets);
        int genJetId_05 = genJetIdOfRecoJet(*bQ, genJets, 0.5);
        if(genJetId>=0) m_histogram["bQuarkT_genJet_directQ_dRmin"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
        if(genJetId_05>=0) {
            putUniquelyInVector(bJt_unique_ids_directQ_dR_05, genJetId_05);
            putUniquelyInVector(bJth_unique_ids_directQ_dR_05, genJetId_05);
            if(genJetId_05!=bJt_id_matched && bJt_id_matched>=0) {
                m_histogram["bQuarkT_genJet_directQ_dRmin_05_notAsMatched"]->Fill(DeltaR(*bQ, genJets.at(genJetId_05)), weight);
                m_histogram["bQuarkTH_genJet_directQ_dRmin_05_notAsMatched"]->Fill(DeltaR(*bQ, genJets.at(genJetId_05)), weight);
                m_histogram["bQuarkT_genJet_directQ_dRmin_05_isAsMatched"]->Fill(0.0, weight);
                m_histogram["bQuarkTH_genJet_directQ_dRmin_05_isAsMatched"]->Fill(0.0, weight);
            } else {
                m_histogram["bQuarkT_genJet_directQ_dRmin_05_isAsMatched"]->Fill(1.0, weight);
                m_histogram["bQuarkTH_genJet_directQ_dRmin_05_isAsMatched"]->Fill(1.0, weight);
            }
        }
    }
    if(bQat) {
        LV* bQ = bQat;
        int genJetId = genJetIdOfRecoJet(*bQ, genJets);
        int genJetId_05 = genJetIdOfRecoJet(*bQ, genJets, 0.5);
        if(genJetId>=0) m_histogram["bQuarkT_genJet_directQ_dRmin"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
        if(genJetId_05>=0) {
            putUniquelyInVector(bJt_unique_ids_directQ_dR_05, genJetId_05);
            putUniquelyInVector(bJth_unique_ids_directQ_dR_05, genJetId_05);
            if(genJetId_05!=bJat_id_matched && bJat_id_matched>=0) {
                m_histogram["bQuarkT_genJet_directQ_dRmin_05_notAsMatched"]->Fill(DeltaR(*bQ, genJets.at(genJetId_05)), weight);
                m_histogram["bQuarkTH_genJet_directQ_dRmin_05_notAsMatched"]->Fill(DeltaR(*bQ, genJets.at(genJetId_05)), weight);
                m_histogram["bQuarkT_genJet_directQ_dRmin_05_isAsMatched"]->Fill(0.0, weight);
                m_histogram["bQuarkTH_genJet_directQ_dRmin_05_isAsMatched"]->Fill(0.0, weight);
            } else {
                m_histogram["bQuarkT_genJet_directQ_dRmin_05_isAsMatched"]->Fill(1.0, weight);
                m_histogram["bQuarkTH_genJet_directQ_dRmin_05_isAsMatched"]->Fill(1.0, weight);
            }
        }
    }
    if(bQh) {
        LV* bQ = bQh;
        int genJetId = genJetIdOfRecoJet(*bQ, genJets);
        int genJetId_05 = genJetIdOfRecoJet(*bQ, genJets, 0.5);
        if(genJetId>=0) m_histogram["bQuarkH_genJet_directQ_dRmin"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
        if(genJetId_05>=0) {
            putUniquelyInVector(bJh_unique_ids_directQ_dR_05, genJetId_05);
            putUniquelyInVector(bJth_unique_ids_directQ_dR_05, genJetId_05);
            if(genJetId_05!=bJh_id_matched && bJh_id_matched >=0) {
                m_histogram["bQuarkH_genJet_directQ_dRmin_05_notAsMatched"]->Fill(DeltaR(*bQ, genJets.at(genJetId_05)), weight);
                m_histogram["bQuarkTH_genJet_directQ_dRmin_05_notAsMatched"]->Fill(DeltaR(*bQ, genJets.at(genJetId_05)), weight);
                m_histogram["bQuarkH_genJet_directQ_dRmin_05_isAsMatched"]->Fill(0.0, weight);
                m_histogram["bQuarkTH_genJet_directQ_dRmin_05_isAsMatched"]->Fill(0.0, weight);
            } else {
                m_histogram["bQuarkH_genJet_directQ_dRmin_05_isAsMatched"]->Fill(1.0, weight);
                m_histogram["bQuarkTH_genJet_directQ_dRmin_05_isAsMatched"]->Fill(1.0, weight);
            }
        }
    }
    if(bQah) {
        LV* bQ = bQah;
        int genJetId = genJetIdOfRecoJet(*bQ, genJets);
        int genJetId_05 = genJetIdOfRecoJet(*bQ, genJets, 0.5);
        if(genJetId>=0) m_histogram["bQuarkH_genJet_directQ_dRmin"]->Fill(DeltaR(*bQ, genJets.at(genJetId)), weight);
        if(genJetId_05>=0) {
            putUniquelyInVector(bJh_unique_ids_directQ_dR_05, genJetId_05);
            putUniquelyInVector(bJth_unique_ids_directQ_dR_05, genJetId_05);
            if(genJetId_05!=bJah_id_matched && bJah_id_matched>=0) {
                m_histogram["bQuarkH_genJet_directQ_dRmin_05_notAsMatched"]->Fill(DeltaR(*bQ, genJets.at(genJetId_05)), weight);
                m_histogram["bQuarkTH_genJet_directQ_dRmin_05_notAsMatched"]->Fill(DeltaR(*bQ, genJets.at(genJetId_05)), weight);
                m_histogram["bQuarkH_genJet_directQ_dRmin_05_isAsMatched"]->Fill(0.0, weight);
                m_histogram["bQuarkTH_genJet_directQ_dRmin_05_isAsMatched"]->Fill(0.0, weight);
            } else {
                m_histogram["bQuarkH_genJet_directQ_dRmin_05_isAsMatched"]->Fill(1.0, weight);
                m_histogram["bQuarkTH_genJet_directQ_dRmin_05_isAsMatched"]->Fill(1.0, weight);
            }
        }
    }

    // Unique jets
    m_histogram["genJetT_unique_multiplicity_directQ_dR_05"]->Fill(bJt_unique_ids_directQ_dR_05.size());
    m_histogram["genJetH_unique_multiplicity_directQ_dR_05"]->Fill(bJh_unique_ids_directQ_dR_05.size());
    m_histogram["genJetTH_unique_multiplicity_directQ_dR_05"]->Fill(bJth_unique_ids_directQ_dR_05.size());



    // Filling plots for jets that are different between dR and advanced matching


    //     if((int)bHt_unique_ids_dR.size() < 2) printf("nHadrons: %d\tnSelJets: %d\tnColesJets: %d\n", (int)bHt_unique_ids_dR.size(), (int)genJets.size(), (int)bJt_unique_ids_dR.size());
}


int AnalyzerDijet::genJetIdOfRecoJet(const LV& recoJet, const VLV& genJets, const float dR_max)
{
    float dR_min = 999.9;
    int genJetId = -1;
    // Loop over gen jets to find the closest one to reco jet
    for(size_t iJet=0, nJets=genJets.size(); iJet<nJets; iJet++) {
        float dR = ROOT::Math::VectorUtil::DeltaR(recoJet,genJets.at(iJet));
        if(dR>=dR_min) continue;
        if(dR>=dR_max) continue;
        dR_min = dR;
        genJetId = iJet;
    }

    return genJetId;
}


std::vector<int> AnalyzerDijet::bHadIdsInGenJet(const int jetId, const std::vector<int>& hadJetIndices)
{
    std::vector<int> hadIndices;

    if(jetId<0) return hadIndices;

    for(size_t iHad=0, nHads=hadJetIndices.size(); iHad<nHads; iHad++) {
        if(hadJetIndices.at(iHad)!=jetId) continue;
        hadIndices.push_back(iHad);
    }

    return hadIndices;
}


std::vector<int> AnalyzerDijet::bHadFlavoursInGenJet(const int jetId, const std::vector<int>& hadJetIndices,
                                                     const std::vector<int>& hadFlavours, const bool absFlavour)
{
    std::vector<int> flavours;
    std::vector<int> hadIndices = bHadIdsInGenJet(jetId, hadJetIndices);

    for(size_t iHad=0, nHads=hadIndices.size(); iHad<nHads; iHad++) {
        int hadId = hadIndices.at(iHad);
        int flavour = hadFlavours.at(hadId);
        if(absFlavour) flavour = std::abs(flavour);

        putUniquelyInVector(flavours, flavour);
    }

    return flavours;
}


bool AnalyzerDijet::isInVector(const std::vector<int>& vector, const int id)
{
    bool isIn = std::find(vector.begin(), vector.end(), id) != vector.end();

    return isIn;
}

bool AnalyzerDijet::isInVector(std::vector<std::pair<int,int> >& vector, const std::pair<int,int> id)
{
    bool isIn = std::find(vector.begin(), vector.end(), id) != vector.end();

    return isIn;
}

bool AnalyzerDijet::putUniquelyInVector(std::vector<int>& vector, const int id)
{
    if(isInVector(vector, id)) return false;

    vector.push_back(id);
    return true;
}

bool AnalyzerDijet::putUniquelyInVector(std::vector<std::pair<int,int> >& vector, const std::pair<int,int> id)
{
    if(isInVector(vector, id)) return false;

    vector.push_back(id);
    return true;
}

std::vector<std::pair<int,int> > AnalyzerDijet::allPairsFromJets(const std::vector<int>& allJets, const std::vector<int>& jetsToIgnore)
{
    std::vector<std::pair<int,int> > jetPairs;
    
    for(size_t jetId1=0; jetId1<allJets.size(); ++jetId1) {
        if( isInVector(jetsToIgnore, allJets.at(jetId1)) ) continue;
        for(size_t jetId2=jetId1+1; jetId2<allJets.size(); ++jetId2) {
            if( isInVector(jetsToIgnore, allJets.at(jetId2)) ) continue;
            
            jetPairs.push_back( std::pair<int,int>(allJets.at(jetId1), allJets.at(jetId2)) );
        }
    }
    
    return jetPairs;
}


bool AnalyzerDijet::areAmongPairs( const std::vector< std::pair<int,int> >& pairs, const int idx1, const int idx2 )
{
    for(size_t iPair=0; iPair<pairs.size(); iPair++) {
        std::pair<int,int> pair = pairs.at(iPair);
//         printf("Chk <%d,%d> | <%d,%d>\n", pair.first, pair.second, idx1, idx2);
        if( (pair.first==idx1 || pair.first==idx2) && (pair.second==idx1 || pair.second==idx2) ) return true;
//         printf("-   <%d,%d> | <%d,%d>\n", pair.first, pair.second, idx1, idx2);
    }
    return false;
}


void AnalyzerDijet::fillTopAdditionalJetsHistos(const EventMetadata& eventMetadata,
                                                const RecoObjects& recoObjects, const TopGenObjects& topGenObjects,
                                                const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                                                const tth::RecoObjectIndices& recoObjectIndices,
                                                const tth::GenObjectIndices& genObjectIndices,
                                                const double& weight, std::map<TString, TH1*>& m_histogram)
{
    // Setting the ordering for gen and reco jet indices
    common::LVParameter ordering = common::LVpt;
    // Extracting input data to more comfortable variables
    const VLV& allJets = (recoObjects.valuesSet_) ? *recoObjects.jets_ : VLV();
    std::vector<int> jetsId = recoObjectIndices.jetIndices_;           // Selected jets (point to jets from allJets)
    std::vector<int> bJetsId = recoObjectIndices.bjetIndices_;         // B-tagged jets (point to jets from allJets)
    std::vector<int> topJetsId_kinReco;                                       // Jets from ttbar by KinReco (Point to jets from allJets)
    std::vector<int> addJetsId_kinReco;                                       // Jets not from ttbar by KinReco (Point to jets from allJets)
    std::vector<int> addBJetsId_kinReco;                                      // B-jets not from ttbar by KinReco (Point to jets from allJets)
    std::vector<int> topJetsId_mva;                                           // Jets from ttbar by MVA (Point to jets from allJets)
    std::vector<int> addJetsId_mva;                                           // Jets not from ttbar by MVA (Point to jets from allJets)
    std::vector<int> addBJetsId_mva;                                          // B-ets not from ttbar by MVA (Point to jets from allJets)
    // Ordering reco jet indices
    common::orderIndices(jetsId, allJets, ordering);
    common::orderIndices(bJetsId, allJets, ordering);
    // Identifying reco jets from tt by KinReco
    if(kinematicReconstructionSolutions.numberOfSolutions()) {
        topJetsId_kinReco.push_back(kinematicReconstructionSolutions.solution().bjetIndex());
        topJetsId_kinReco.push_back(kinematicReconstructionSolutions.solution().antiBjetIndex());
    }

    // Identifying reco jets from tt by MVA
    std::vector<float> v_mvaWeights;
    if(recoObjectIndices.jetIndexPairs_.size()>0) {
	std::vector<MvaVariablesBase*> v_mvaVariables = MvaVariablesTopJets::fillVariables(recoObjectIndices, genObjectIndices, recoObjects, weight);
        if(weightsCorrect_) {
            v_mvaWeights = weightsCorrect_->mvaWeights(v_mvaVariables);
            const tth::IndexPairs& jetIndexPairs = recoObjectIndices.jetIndexPairs_;
            std::vector<int> jetIndexPairsIndices = common::initialiseIndices(jetIndexPairs);
            common::orderIndices(jetIndexPairsIndices, v_mvaWeights);
            if(jetIndexPairsIndices.size()>0) {
                topJetsId_mva.push_back(jetIndexPairs.at(jetIndexPairsIndices.at(0)).first);
                topJetsId_mva.push_back(jetIndexPairs.at(jetIndexPairsIndices.at(0)).second);
            }
        }
    }
    // Setting variables of gen. level if available
    const VLV& allGenJets = (topGenObjects.valuesSet_) ? *topGenObjects.allGenJets_ : VLV(0);
    std::vector<int> genJetsId = genObjectIndices.genJetIndices_;
    std::vector<int> genBJetsId = genObjectIndices.genBjetIndices_;
    std::vector<int> genJetsRecoId;
//     const std::vector<int>& bHadJetIndex = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadJetIndex_ : std::vector<int>(0);
//     const std::vector<int>& bHadFlavour = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadFlavour_ : std::vector<int>(0);
//     const std::vector<int>& bHadIndex = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadIndex_ : std::vector<int>(0);
//     const std::vector<int>& bHadFromTopWeakDecay = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadFromTopWeakDecay_ : std::vector<int>(0);
    std::vector<int> topJetsId_true;          // Reco jets coming from top (Point to jets from allJets) [Can be only selected jet]
    std::vector<int> addJetsId_true;          // Reco additional jets (Point to jets from allJets) [Can be only selected jet]
    std::vector<int> addBJetsId_true;         // Reco additional b-jets (Point to jets from allJets) [Can be only selected jet]
    std::vector<int> topAllJetsId_gen;        // Gen jets coming from top (Point to jets from allJets)
    std::vector<int> topJetsId_gen;           // Gen jets coming from top (Point to jets from allGenJets) [Can be only selected jet]
    std::vector<int> addJetsId_gen;           // Gen additional jets (Point to jets from allGenJets) [Can be only selected jet]
    std::vector<int> addBJetsId_gen;          // Gen additional b-jets (Point to jets from allGenJets) [Can be only selected jet]
    
    // Selecting true reco jets from top
    if(genObjectIndices.recoBjetFromTopIndex_>=0) putUniquelyInVector(topJetsId_true ,genObjectIndices.recoBjetFromTopIndex_);
    if(genObjectIndices.recoAntiBjetFromTopIndex_>=0) putUniquelyInVector(topJetsId_true, genObjectIndices.recoAntiBjetFromTopIndex_);
    
    // Reordering all jet indices by the appropriate variable
    common::orderIndices(topJetsId_true, allJets, ordering);
    common::orderIndices(topJetsId_kinReco, allJets, ordering);
    common::orderIndices(topJetsId_mva, allJets, ordering);
    
    
    // Selecting gen jets/b-jets in acceptance
    for(int jetId = 0; jetId<(int)allGenJets.size(); ++jetId) {
        if(jetId==genObjectIndices.genBjetFromTopIndex_) topAllJetsId_gen.push_back(jetId);
        if(jetId==genObjectIndices.genAntiBjetFromTopIndex_) topAllJetsId_gen.push_back(jetId);
        if(!isInVector(genJetsId, jetId)) continue;
        // Jets from top
        if(jetId==genObjectIndices.genBjetFromTopIndex_) {
            topJetsId_gen.push_back(jetId);
            continue;
        }
        if(jetId==genObjectIndices.genAntiBjetFromTopIndex_) {
            topJetsId_gen.push_back(jetId);
            continue;
        }
        // Additional jets/b-jets
        addJetsId_gen.push_back(jetId);
        if(isInVector(genBJetsId, jetId)) addBJetsId_gen.push_back(jetId);
    }
    m_histogram["jet_top_mult_genAll"]->Fill(topAllJetsId_gen.size(), weight);
    m_histogram["jet_top_mult_gen"]->Fill(topJetsId_gen.size(), weight);
    
    // Ordering reco jet indices
    common::orderIndices(genJetsId, allGenJets, ordering);
    common::orderIndices(genBJetsId, allGenJets, ordering);
    common::orderIndices(topJetsId_gen, allGenJets, ordering);
    common::orderIndices(topAllJetsId_gen, allGenJets, ordering);
    common::orderIndices(addJetsId_gen, allGenJets, ordering);
    common::orderIndices(addBJetsId_gen, allGenJets, ordering);
    
//    printf("\nAll gen b-jets:  %d\n", (int)genBJetsId.size());
//    for(int jetId : genBJetsId) printf("%d. \tPt: %.3f \tEta: %.3f\n", jetId, allGenJets.at(jetId).Pt(), allGenJets.at(jetId).Eta());
//    printf("All gen jets: %d\n", (int)genJetsId.size());
//    for(int jetId : genJetsId) printf("%d. \tPt: %.3f \tEta: %.3f\n", jetId, allGenJets.at(jetId).Pt(), allGenJets.at(jetId).Eta());
    
    // Matching each reco jet to the gen jet
    genJetsRecoId = std::vector<int>(allGenJets.size(), -1);
    for(size_t index= 0; index<allGenJets.size(); ++index){
        genJetsRecoId.at(index) = matchRecoToGenJet(jetsId, allJets, index, allGenJets);
    }
    
//     for(int id : genJetsId) {
//         printf("Gen: Pt: %.2f Eta: %.2f Phi: %.2f\n", allGenJets.at(id).Pt(), allGenJets.at(id).Eta(), allGenJets.at(id).Phi());
//         int recoId = genJetsRecoId.at(id);
//         if(recoId<0) printf("  -1\n");
//         else printf("  Reco: Pt: %.2f Eta: %.2f Phi: %.2f\n", allJets.at(recoId).Pt(), allJets.at(recoId).Eta(), allJets.at(recoId).Phi());
//     }
    
    // Selecting additional jets/b-jets
    for(int jetId : jetsId) { 
        if(!isInVector(topJetsId_kinReco, jetId)) addJetsId_kinReco.push_back(jetId);
        if(!isInVector(topJetsId_mva, jetId)) addJetsId_mva.push_back(jetId);
        if(!isInVector(topJetsId_true, jetId)) addJetsId_true.push_back(jetId);
        if(!isInVector(bJetsId, jetId)) continue;
        if(!isInVector(topJetsId_kinReco, jetId)) addBJetsId_kinReco.push_back(jetId);
        if(!isInVector(topJetsId_mva, jetId)) addBJetsId_mva.push_back(jetId);
        if(!isInVector(topJetsId_true, jetId)) addBJetsId_true.push_back(jetId);
    }
    
//     printf("\nReco jets:\n");
//     for(int jetId : jetsId) {
//         LV jet = allJets.at(jetId);
//         printf("%d. pT: %.2f eta: %.2f\n", jetId, jet.Pt(), jet.Eta());
//     }
//     printf("Reco add jets true:\n");
//     for(int jetId : addJetsId_true) {
//         LV jet = allJets.at(jetId);
//         printf("%d. pT: %.2f eta: %.2f\n", jetId, jet.Pt(), jet.Eta());
//     }
//     printf("Reco add jets kinReco:\n");
//     for(int jetId : addJetsId_kinReco) {
//         LV jet = allJets.at(jetId);
//         printf("%d. pT: %.2f eta: %.2f\n", jetId, jet.Pt(), jet.Eta());
//     }
//     printf("Gen jets:\n");
//     for(int jetId : genJetsId) {
//         LV jet = allGenJets.at(jetId);
//         printf("%d. pT: %.2f eta: %.2f\n", jetId, jet.Pt(), jet.Eta());
//     }
//     printf("Gen add jets:\n");
//     for(int jetId : addJetsId_gen) {
//         LV jet = allGenJets.at(jetId);
//         printf("%d. pT: %.2f eta: %.2f\n", jetId, jet.Pt(), jet.Eta());
//     }
    
    int nJets_top_true = 0;
    int nBJets_top_true = 0;
    int nJets_top_kinReco = topJetsId_kinReco.size();
    int nJets_top_kinReco_corr = 0;
    int nBJets_top_kinReco = 0;
    int nJets_top_mva = topJetsId_mva.size();
    int nJets_top_mva_corr = 0;
    int nBJets_top_mva = 0;
    int nBJets_add_gen = addBJetsId_gen.size();
    m_histogram["bJet_add_mult_gen"]->Fill(nBJets_add_gen, weight);
    // Counting true top jets/b-jets
    for(int jetId : topJetsId_true) {
        if(isInVector(jetsId, jetId)) nJets_top_true++;
        if(isInVector(bJetsId, jetId)) nBJets_top_true++;
        if(isInVector(topJetsId_kinReco, jetId)) nJets_top_kinReco_corr++;
        if(isInVector(topJetsId_mva, jetId)) nJets_top_mva_corr++;
    }
    if(nJets_top_true==2 && !genObjectIndices.uniqueRecoTopMatching()) { nJets_top_true--; nBJets_top_true--; }
    m_histogram["jet_top_mult_true"]->Fill(nJets_top_true, weight);
    m_histogram["bJet_top_mult_true"]->Fill(nBJets_top_true, weight);
    // Counting KinReco top jets/b-jets
    for(int jetId : topJetsId_kinReco) {
        if(isInVector(bJetsId, jetId)) nBJets_top_kinReco++;
    }
    m_histogram["jet_top_mult_kinReco"]->Fill(nJets_top_kinReco, weight);
    m_histogram["jet_top_mult_kinReco_corr"]->Fill(nJets_top_kinReco_corr, weight);
    m_histogram["bJet_top_mult_kinReco"]->Fill(nBJets_top_kinReco, weight);
    // Counting MVA top jets/b-jets
    for(int jetId : topJetsId_mva) {
        if(isInVector(bJetsId, jetId)) nBJets_top_mva++;
    }
    m_histogram["jet_top_mult_mva"]->Fill(nJets_top_mva, weight);
    m_histogram["jet_top_mult_mva_corr"]->Fill(nJets_top_mva_corr, weight);
    m_histogram["bJet_top_mult_mva"]->Fill(nBJets_top_mva, weight);
    
    // Checking how many proper jets assigned if KinReco and MVA find the same jets
    if(topJetsId_mva.size() == topJetsId_kinReco.size() && topJetsId_kinReco.size() == 2) {
        if(isInVector(topJetsId_kinReco, topJetsId_mva.at(0)) && isInVector(topJetsId_kinReco, topJetsId_mva.at(1))) {
            m_histogram["jet_top_mult_same2Jets_corr"]->Fill(nJets_top_kinReco_corr, weight);
        }
    }
    // Counting N hadrons in each additional b jet
    if(nBJets_add_gen>0) m_histogram["leadingJet_1st_nHad_addB_gen"]->Fill(genObjectIndices.genJetBhadronIndices_.at(addBJetsId_gen.at(0)).size());
    if(nBJets_add_gen>1) m_histogram["leadingJet_2nd_nHad_addB_gen"]->Fill(genObjectIndices.genJetBhadronIndices_.at(addBJetsId_gen.at(1)).size());
    
    
    // Filling histograms about leading jets with tt jets identified at generator level (GEN -> RECO)
//     fillLeadingJetsHistosVsGen("top_gen", eventMetadata, allGenJets, topAllJetsId_gen, allGenJets, topAllJetsId_gen, genJetsRecoId, topAllJetsId_gen, topAllJetsId_gen, weight, m_histogram, recoObjects, false);
//     fillLeadingJetsHistosVsGen("add_gen", eventMetadata, allGenJets, addJetsId_gen, allGenJets, addJetsId_gen, genJetsRecoId, topAllJetsId_gen, topAllJetsId_gen, weight, m_histogram, recoObjects, false);
    fillLeadingJetsHistosVsGen("addB_gen", eventMetadata, allGenJets, addBJetsId_gen, allGenJets, addBJetsId_gen, genJetsRecoId, topAllJetsId_gen, topAllJetsId_gen, weight, m_histogram, recoObjects, false);
    
    // Filling histograms about leading jets with tt jets identified by True
    fillLeadingJetsHistosVsGen("top_true", eventMetadata, allGenJets, topAllJetsId_gen, allJets, topJetsId_true, genJetsRecoId, topAllJetsId_gen, topJetsId_true, weight, m_histogram, recoObjects, false);
    fillLeadingJetsHistosVsGen("add_true", eventMetadata, allGenJets, addJetsId_gen, allJets, addJetsId_true, genJetsRecoId, topAllJetsId_gen, topJetsId_true, weight, m_histogram, recoObjects, false);
    fillLeadingJetsHistosVsGen("addB_true", eventMetadata, allGenJets, addBJetsId_gen, allJets, addBJetsId_true, genJetsRecoId, topAllJetsId_gen, topJetsId_true, weight, m_histogram, recoObjects, false);
    // Control plots
    fillLeadingJetsHistosVsGen("top_cp_true", eventMetadata, allGenJets, topAllJetsId_gen, allJets, topJetsId_true, genJetsRecoId, topAllJetsId_gen, topJetsId_true, weight, m_histogram, recoObjects, false);
    fillLeadingJetsHistosVsGen("add_cp_true", eventMetadata, allGenJets, addJetsId_gen, allJets, addJetsId_true, genJetsRecoId, topAllJetsId_gen, topJetsId_true, weight, m_histogram, recoObjects, false);
    fillLeadingJetsHistosVsGen("addB_cp_true", eventMetadata, allGenJets, addBJetsId_gen, allJets, addBJetsId_true, genJetsRecoId, topAllJetsId_gen, topJetsId_true, weight, m_histogram, recoObjects, false);
    
    // Filling histograms about leading jets with tt jets identified by KinReco
    fillLeadingJetsHistosVsGen("top_kinReco", eventMetadata, allGenJets, topAllJetsId_gen, allJets, topJetsId_kinReco, genJetsRecoId, topAllJetsId_gen, topJetsId_kinReco, weight, m_histogram, recoObjects);
    fillLeadingJetsHistosVsGen("add_kinReco", eventMetadata, allGenJets, addJetsId_gen, allJets, addJetsId_kinReco, genJetsRecoId, topAllJetsId_gen, topJetsId_kinReco, weight, m_histogram, recoObjects);
    fillLeadingJetsHistosVsGen("addB_kinReco", eventMetadata, allGenJets, addBJetsId_gen, allJets, addBJetsId_kinReco, genJetsRecoId, topAllJetsId_gen, topJetsId_kinReco, weight, m_histogram, recoObjects);
    // Control plots
    fillLeadingJetsHistosVsGen("top_cp_kinReco", eventMetadata, allGenJets, topAllJetsId_gen, allJets, topJetsId_kinReco, genJetsRecoId, topAllJetsId_gen, topJetsId_kinReco, weight, m_histogram, recoObjects);
    fillLeadingJetsHistosVsGen("add_cp_kinReco", eventMetadata, allGenJets, addJetsId_gen, allJets, addJetsId_kinReco, genJetsRecoId, topAllJetsId_gen, topJetsId_kinReco, weight, m_histogram, recoObjects);
    fillLeadingJetsHistosVsGen("addB_cp_kinReco", eventMetadata, allGenJets, addBJetsId_gen, allJets, addBJetsId_kinReco, genJetsRecoId, topAllJetsId_gen, topJetsId_kinReco, weight, m_histogram, recoObjects);
    
    // Filling histograms about leading jets with tt jets identified by MVA
    fillLeadingJetsHistosVsGen("top_mva", eventMetadata, allGenJets, topAllJetsId_gen, allJets, topJetsId_mva, genJetsRecoId, topAllJetsId_gen, topJetsId_mva, weight, m_histogram, recoObjects);
    fillLeadingJetsHistosVsGen("add_mva", eventMetadata, allGenJets, addJetsId_gen, allJets, addJetsId_mva, genJetsRecoId, topAllJetsId_gen, topJetsId_mva, weight, m_histogram, recoObjects);
    fillLeadingJetsHistosVsGen("addB_mva", eventMetadata, allGenJets, addBJetsId_gen, allJets, addBJetsId_mva, genJetsRecoId, topAllJetsId_gen, topJetsId_mva, weight, m_histogram, recoObjects);
    // Control plots
    fillLeadingJetsHistosVsGen("top_cp_mva", eventMetadata, allGenJets, topAllJetsId_gen, allJets, topJetsId_mva, genJetsRecoId, topAllJetsId_gen, topJetsId_mva, weight, m_histogram, recoObjects);
    fillLeadingJetsHistosVsGen("add_cp_mva", eventMetadata, allGenJets, addJetsId_gen, allJets, addJetsId_mva, genJetsRecoId, topAllJetsId_gen, topJetsId_mva, weight, m_histogram, recoObjects);
    fillLeadingJetsHistosVsGen("addB_cp_mva", eventMetadata, allGenJets, addBJetsId_gen, allJets, addBJetsId_mva, genJetsRecoId, topAllJetsId_gen, topJetsId_mva, weight, m_histogram, recoObjects);
    
    // Filling histograms about leading jets with tt jets identified at reconstructed level (RECO -> GEN)
    // tt jets from the Kinematic Reconstruction
    if(topJetsId_kinReco.size() == 2) {
        fillLeadingJetsHistosVsGen("top_vis_kinReco", eventMetadata, allGenJets, topAllJetsId_gen, allJets, topJetsId_kinReco, genJetsRecoId, topAllJetsId_gen, topJetsId_kinReco,  weight, m_histogram, recoObjects);
        fillLeadingJetsHistosVsGen("add_vis_kinReco", eventMetadata, allGenJets, addJetsId_gen, allJets, addJetsId_kinReco, genJetsRecoId, topAllJetsId_gen, topJetsId_kinReco,  weight, m_histogram, recoObjects);
        fillLeadingJetsHistosVsGen("addB_vis_kinReco", eventMetadata, allGenJets, addBJetsId_gen, allJets, addBJetsId_kinReco, genJetsRecoId, topAllJetsId_gen, topJetsId_kinReco,  weight, m_histogram, recoObjects);
    }
    // tt jets from the MVA
    if(topJetsId_mva.size() == 2) {
        fillLeadingJetsHistosVsGen("top_vis_mva", eventMetadata, allGenJets, topAllJetsId_gen, allJets, topJetsId_mva, genJetsRecoId, topAllJetsId_gen, topJetsId_mva,  weight, m_histogram, recoObjects);
        fillLeadingJetsHistosVsGen("add_vis_mva", eventMetadata, allGenJets, addJetsId_gen, allJets, addJetsId_mva, genJetsRecoId, topAllJetsId_gen, topJetsId_mva,  weight, m_histogram, recoObjects);
        fillLeadingJetsHistosVsGen("addB_vis_mva", eventMetadata, allGenJets, addBJetsId_gen, allJets, addBJetsId_mva, genJetsRecoId, topAllJetsId_gen, topJetsId_mva,  weight, m_histogram, recoObjects);
    }
    
    
    // Correlation hitograms for true vs Mva/KinReco jets
    std::vector<int> sameJets_kinReco;
    std::vector<int> sameJets_mva;
    for(size_t iJet = 0; iJet<2; ++iJet) {
        int trueAllJetId = (iJet<topJetsId_true.size()) ? topJetsId_true.at(iJet) : -1;
        int trueJetId = std::find(jetsId.begin(), jetsId.end(), trueAllJetId) - jetsId.begin();
        if(trueJetId>=(int)jetsId.size()) trueJetId=-1;
        // Finding true jets from KinReco
        int sameJet_kinReco = -1;
        for(int jetId : topJetsId_kinReco) {
            if(trueJetId != jetId) continue;
            sameJet_kinReco = jetId;
            break;
        }
        sameJets_kinReco.push_back(sameJet_kinReco);
        // Finding true jets from MVA
        int sameJet_mva = -1;
        for(int jetId : topJetsId_mva) {
            if(trueJetId != jetId) continue;
            sameJet_mva = jetId;
            break;
        }
        sameJets_mva.push_back(sameJet_mva);
    }
    // Setting remaining jet to true jets that haven't been selected by KinReco/MVA
    for(size_t iJet = 0; iJet<2; ++iJet) {
        if(sameJets_kinReco.at(iJet) < 0) {
            for(int jetId : topJetsId_kinReco) {
                if(isInVector(sameJets_kinReco, jetId)) continue;
                sameJets_kinReco.at(iJet) = jetId;
            }
        }
        if(sameJets_mva.at(iJet) < 0) {
            for(int jetId : topJetsId_mva) {
                if(isInVector(sameJets_mva, jetId)) continue;
                sameJets_mva.at(iJet) = jetId;
            }
        }
    }
    // Getting indices of jets coming from the tt system as proposed by true/kinReco/MVA
    std::vector<int> jetIds_kinReco, jetIds_mva, jetIds_true;
    for(size_t iJet=0; iJet<2; ++iJet ) {
        int jetAllId_true = (iJet<topJetsId_true.size()) ? topJetsId_true.at(iJet) : -1;
        int jetId_true = (jetAllId_true>=0) ? std::find(jetsId.begin(), jetsId.end(), jetAllId_true) - jetsId.begin() : -1;
        if(jetId_true>=(int)jetsId.size()) jetId_true = -1;
        jetIds_true.push_back(jetId_true);
        
        int jetAllId_kinReco = (iJet<sameJets_kinReco.size()) ? sameJets_kinReco.at(iJet) : -1;
        int jetId_kinReco = (jetAllId_kinReco>=0) ? std::find(jetsId.begin(), jetsId.end(), jetAllId_kinReco) - jetsId.begin() : -1;
        if(jetId_kinReco>=(int)jetsId.size()) jetId_kinReco = -1;
        jetIds_kinReco.push_back(jetId_kinReco);
        
        int jetAllId_mva = (iJet<sameJets_mva.size()) ? sameJets_mva.at(iJet) : -1;
        int jetId_mva = (jetAllId_mva>=0) ? std::find(jetsId.begin(), jetsId.end(), jetAllId_mva) - jetsId.begin() : -1;
        if(jetId_mva>=(int)jetsId.size()) jetId_mva = -1;
        jetIds_mva.push_back(jetId_mva);
    }
    
    // Filling actual histograms for each of the jets from tt as proposed by true/kinReco/MVA
    for(size_t iJet=0; iJet<2; ++iJet ) {
        const int jetId_true = jetIds_true.at(iJet);
        const int jetId_kinReco = jetIds_kinReco.at(iJet);
        const int jetId_mva = jetIds_mva.at(iJet);
        std::string hAdd[3] = {"","_sameJet","_same2Jets"};
        for(int i=0; i<3; ++i) {
            if(i==1) {
                if(jetId_kinReco < 0) continue;
                if(jetId_mva < 0) continue;
                if(jetId_kinReco != jetId_mva) continue;
            }
            else if(i==2) {
                if(jetIds_kinReco.at(1) < 0 || jetIds_kinReco.at(0)  < 0 ) continue;
                if(jetIds_mva.at(1) < 0 || jetIds_mva.at(0) < 0) continue;
                if(jetIds_kinReco.at(1) != jetIds_mva.at(1)) continue;
                if(jetIds_kinReco.at(0) != jetIds_mva.at(0)) continue;
            }
            ((TH2*)m_histogram["jet_top_id"+hAdd[i]+"_trueVsKinReco"])->Fill(jetId_true, sameJets_kinReco.at(iJet), weight);
            if(jetId_true>=0 && jetId_kinReco>=0) {
                ((TH2*)m_histogram["jet_top_Pt"+hAdd[i]+"_trueVsKinReco"])->Fill(allJets.at(jetsId.at(jetId_true)).Pt(), allJets.at(jetsId.at(jetId_kinReco)).Pt(), weight);
            }
            
            ((TH2*)m_histogram["jet_top_id"+hAdd[i]+"_trueVsMva"])->Fill(jetId_true, sameJets_mva.at(iJet), weight);
            if(jetId_true>=0 && jetId_mva>=0) {
                ((TH2*)m_histogram["jet_top_Pt"+hAdd[i]+"_trueVsMva"])->Fill(allJets.at(jetsId.at(jetId_true)).Pt(), allJets.at(jetsId.at(jetId_mva)).Pt(), weight);
            }
            
            ((TH2*)m_histogram["jet_top_id"+hAdd[i]+"_kinRecoVsMva"])->Fill(sameJets_kinReco.at(iJet), sameJets_mva.at(iJet), weight);
            if(jetId_mva>=0 && jetId_kinReco>=0) {
                ((TH2*)m_histogram["jet_top_Pt"+hAdd[i]+"_kinRecoVsMva"])->Fill(allJets.at(jetsId.at(jetId_kinReco)).Pt(), allJets.at(jetsId.at(jetId_mva)).Pt(), weight);
            }
        }
    }
    
}


void AnalyzerDijet::fillLeadingJetsHistosVsGen(const std::string& name,
                                               const EventMetadata&, const VLV& allGenJets,
                                               const std::vector<int>& genJetsId, const VLV& allJets, 
                                               const std::vector<int>& jetsId, const std::vector<int>& genJetsRecoId,
                                               const std::vector<int>& topJetsId_gen, const std::vector<int>& topJetsId_reco,
                                               const double& weight, std::map<TString, TH1*>& m_histogram,
                                               const RecoObjects& recoObjects, const bool require2TopJets) 
{
    
    double no = -100.;
    const int nTopJets_reco = topJetsId_reco.size();
    
    const std::vector<double>& bTagDiscriminator = (recoObjects.valuesSet_) ? *recoObjects.jetBTagCSV_ : std::vector<double>(0);
    
    std::vector<std::string> idStr;
    idStr.push_back("1st");
    idStr.push_back("2nd");
    idStr.push_back("3rd");
    idStr.push_back("4th");
    const size_t nJetsToPlot = idStr.size();
    std::string histoName;
    // Properties of leading jets
    for(size_t iJet=0; iJet<jetsId.size(); ++iJet) {
        if(iJet >= nJetsToPlot) continue;
        int jetId = jetsId.at(iJet);
        LV jet = allJets.at(jetId);
        
        histoName = "leadingJet_"+idStr.at(iJet)+"_Pt_"+name;
        if(m_histogram[histoName]) m_histogram[histoName]->Fill(jet.Pt(), weight);
        histoName = "leadingJet_"+idStr.at(iJet)+"_Eta_"+name;
        if(m_histogram[histoName]) m_histogram[histoName]->Fill(jet.Eta(), weight);
        histoName = "leadingJets_btagDiscriminator_"+name;
        if(m_histogram[histoName]) m_histogram[histoName]->Fill(bTagDiscriminator.at(jetId), weight);
        
//         // Correlation to generator level
//         double genJetPt = no;
//         double genJetEta = no;
//         if(genJetsId.size()>iJet) {
//             int genJetId = genJetsId.at(iJet);
//             LV genJet = allGenJets.at(genJetId);
//             genJetPt = genJet.Pt();
//             genJetEta = genJet.Eta();
//         } 
//         else if(!fillAllGen) continue;
//         
//         histoName = "leadingJet_"+idStr.at(iJet)+"_Pt_"+name+"VsGen";
//         if(m_histogram[histoName]) ((TH2*)m_histogram[histoName])->Fill(jet.Pt(), genJetPt, weight);
//         histoName = "leadingJet_"+idStr.at(iJet)+"_Eta_"+name+"VsGen";
//         if(m_histogram[histoName]) ((TH2*)m_histogram[histoName])->Fill(jet.Eta(), genJetEta, weight);
    }
    // Properties of jet pairs
    double reco_dR = no;
    double reco_dPhi = no;
    double reco_Mjj = no;
    double reco_Pt1 = no;
    double reco_Pt2 = no;
    double reco_Eta1 = no;
    double reco_Eta2 = no;
    LV *jet_1(0), *jet_2(0);
    if(jetsId.size() > 0 && (!require2TopJets || nTopJets_reco > 1)) {
        jet_1 = new LV(allJets.at(jetsId.at(0)));
        reco_Pt1 = jet_1->Pt();
        reco_Eta1 = jet_1->Eta();
    }
    if(jetsId.size() > 1 && (!require2TopJets || nTopJets_reco > 1)) {
        jet_2 = new LV(allJets.at(jetsId.at(1)));
        reco_Pt2 = jet_2->Pt();
        reco_Eta2 = jet_2->Eta();
        reco_dR = ROOT::Math::VectorUtil::DeltaR(*jet_1, *jet_2);
        reco_dPhi = ROOT::Math::VectorUtil::DeltaPhi(*jet_1, *jet_2);
        LV dijet = *jet_1 + *jet_2;
        reco_Mjj = dijet.M();
    }

    m_histogram["leadingJet_dR_"+name]->Fill(reco_dR, weight);
    m_histogram["leadingJet_dPhi_"+name]->Fill(reco_dPhi, weight);
    m_histogram["leadingJet_Mjj_"+name]->Fill(reco_Mjj, weight);
    
    histoName = "leadingJet_1st_Pt_2nd_Pt_"+name;
    if(m_histogram[histoName]) ((TH2*)m_histogram[histoName])->Fill(reco_Pt1, reco_Pt2, weight);
        
    // Correlation to generator level
    double gen_dR = no;
    double gen_dPhi = no;
    double gen_Mjj = no;
    double gen_Pt1 = no;
    double gen_Pt2 = no;
    double gen_Eta1 = no;
    double gen_Eta2 = no;
    
    LV *genJet_1(0), *genJet_2(0);
    if(genJetsId.size() > 0) {
        genJet_1 = new LV(allGenJets.at(genJetsId.at(0)));
        gen_Pt1 = genJet_1->Pt();
        gen_Eta1 = genJet_1->Eta();
    }
    if(genJetsId.size() > 1) {
        genJet_2 = new LV(allGenJets.at(genJetsId.at(1)));
        gen_Pt2 = genJet_2->Pt();
        gen_Eta2 = genJet_2->Eta();
        gen_dR = ROOT::Math::VectorUtil::DeltaR(*genJet_1, *genJet_2);
        gen_dPhi = ROOT::Math::VectorUtil::DeltaPhi(*genJet_1, *genJet_2);
        LV dijet = *genJet_1 + *genJet_2;
        gen_Mjj = dijet.M();
    }
    

//     if(name=="addB_kinReco") {
//             if(reco_dR > 2.5 && reco_dR < 3. && gen_dR < 1. && gen_dR > 0.)
//             {
//                 const int nTopJets_gen = topJetsId_gen.size();
//                 std::cout << "              " << name << std::endl;
//                 std::cout << "  Event: " << eventMetadata.eventNumber_ << " Lumi: " << eventMetadata.lumiBlock_ << " weight: " << weight << std::endl;
//                 std::cout << "  reco_dR: " << reco_dR << " gen_dR: " << gen_dR << std::endl;
//                 std::cout << "  nJetsAdd_reco: " << jetsId.size() << " nJetsAdd_gen: " << genJetsId.size() << " nJetsTop_reco: " << nTopJets_reco << " nJetsTop_gen: " << nTopJets_gen << std::endl;
//                 std::cout << std::endl;
//                 std::cout << "Top jets: GEN" << std::endl;
//                 for(int jetId : topJetsId_gen) printf("%d. Pt: %.2f  Eta: %.2f\n", jetId, allGenJets.at(jetId).Pt(), allGenJets.at(jetId).Eta());
//                 std::cout << "Top jets: RECO" << std::endl;
//                 for(int jetId : topJetsId_reco) printf("%d. Pt: %.2f  Eta: %.2f\n", jetId, allJets.at(jetId).Pt(), allJets.at(jetId).Eta());
//                 std::cout << std::endl;
//                 std::cout << "Add jets: GEN" << std::endl;
//                 for(int jetId : genJetsId) printf("%d. Pt: %.2f  Eta: %.2f\n", jetId, allGenJets.at(jetId).Pt(), allGenJets.at(jetId).Eta());
//                 std::cout << "Add jets: RECO" << std::endl;
//                 for(int jetId : jetsId) printf("%d. Pt: %.2f  Eta: %.2f\n", jetId, allJets.at(jetId).Pt(), allJets.at(jetId).Eta());
//                 std::cout << std::endl;
//                 std::cout << "All jets: RECO" << std::endl;
//                 std::vector<int> jetsIndices = common::initialiseIndices(allJets);
//                 for(int jetId : jetsIndices) printf("%d. Pt: %.2f  Eta: %.2f  bTag: %.2f\n", jetId, allJets.at(jetId).Pt(), allJets.at(jetId).Eta(), recoObjects.jetBTagCSV_->at(jetId));
// 
//                 std::cout << std::endl;
//                 std::cout << std::endl;
//             }
//     }
    histoName = "leadingJet_dR_"+name+"VsGen";
    if(m_histogram[histoName]) ((TH2*)m_histogram[histoName])->Fill(reco_dR, gen_dR, weight);
    histoName = "leadingJet_dPhi_"+name+"VsGen";
    if(m_histogram[histoName]) ((TH2*)m_histogram[histoName])->Fill(reco_dPhi, gen_dPhi, weight);
    histoName = "leadingJet_Mjj_"+name+"VsGen";
    if(m_histogram[histoName]) ((TH2*)m_histogram[histoName])->Fill(reco_Mjj, gen_Mjj, weight);
    
    histoName = "leadingJet_1st_Pt_"+name+"VsGen";
    if(m_histogram[histoName]) ((TH2*)m_histogram[histoName])->Fill(reco_Pt1, gen_Pt1, weight);
    histoName = "leadingJet_1st_Eta_"+name+"VsGen";
    if(m_histogram[histoName]) ((TH2*)m_histogram[histoName])->Fill(reco_Eta1, gen_Eta1, weight);
    histoName = "leadingJet_2nd_Pt_"+name+"VsGen";
    if(m_histogram[histoName]) ((TH2*)m_histogram[histoName])->Fill(reco_Pt2, gen_Pt2, weight);
    histoName = "leadingJet_2nd_Eta_"+name+"VsGen";
    if(m_histogram[histoName]) ((TH2*)m_histogram[histoName])->Fill(reco_Eta2, gen_Eta2, weight);
    
    for(size_t iJet_gen = 0; iJet_gen<genJetsId.size(); ++iJet_gen) {
        int genJetId = genJetsId.at(iJet_gen);
        int recoJetId = genJetsRecoId.at(genJetId);
        int index_reco = std::find(jetsId.begin(), jetsId.end(), recoJetId) - jetsId.begin();
        if(index_reco>=(int)jetsId.size()) index_reco = -1;
        histoName = "leadingJet_index_"+name+"VsGen";
        if(m_histogram[histoName]) ((TH2*)m_histogram[histoName])->Fill(index_reco, iJet_gen, weight);
        if(index_reco < 0) continue;
        histoName = "leadingJet_Pt_"+name+"VsGen";
        if(m_histogram[histoName]) ((TH2*)m_histogram[histoName])->Fill(allJets.at(recoJetId).Pt(), allGenJets.at(genJetId).Pt(), weight);
        histoName = "leadingJet_Eta_"+name+"VsGen";
        if(m_histogram[histoName]) ((TH2*)m_histogram[histoName])->Fill(allJets.at(recoJetId).Eta(), allGenJets.at(genJetId).Eta(), weight);
        if(topJetsId_gen.size()<2 || topJetsId_reco.size()<2) continue;
        LV tt_gen = allGenJets.at(topJetsId_gen.at(0)) + allGenJets.at(topJetsId_gen.at(1));
        LV tt_reco = allJets.at(topJetsId_reco.at(0)) + allJets.at(topJetsId_reco.at(1));
        double dR_gen = ROOT::Math::VectorUtil::DeltaR(tt_gen, allGenJets.at(genJetId));
        double dR_reco = ROOT::Math::VectorUtil::DeltaR(tt_reco, allJets.at(recoJetId));
        histoName = "leadingJet_ttdR_"+name+"VsGen";
        if(m_histogram[histoName]) ((TH2*)m_histogram[histoName])->Fill(dR_reco, dR_gen, weight);
        double dEta_gen = tt_gen.Eta() - allGenJets.at(genJetId).Eta();
        double dEta_reco = tt_reco.Eta() - allJets.at(recoJetId).Eta();
        histoName = "leadingJet_ttdEta_"+name+"VsGen";
        if(m_histogram[histoName]) ((TH2*)m_histogram[histoName])->Fill(dEta_reco, dEta_gen, weight);
    }    
}


void AnalyzerDijet::fillLeadingJetsHistosVsTrue(const std::string& name, const std::vector<int>& trueJetsId,
                                                const std::vector<int>& jetsId, const double& weight, std::map<TString, TH1*>& m_histogram)
{
    if(jetsId.size()<2) return;
    if(trueJetsId.size()<2) return;
    
    unsigned int nJets_same = 0;
    unsigned int nJets_sameOrder = 0;
    for(size_t iJet_true = 0; iJet_true<trueJetsId.size(); ++iJet_true) {
        if(iJet_true>=2) continue;
        int jetId_true = trueJetsId.at(iJet_true);
        if(jetId_true<0) continue;
        if(!isInVector(jetsId, jetId_true)) continue;
        nJets_same++;
        for(size_t iJet = 0; iJet < jetsId.size(); ++iJet) {
            if(iJet>=2) continue;
            int jetId = jetsId.at(iJet);
            if(jetId<0) continue;
            if(iJet != iJet_true) continue;
            if(jetId != jetId_true) continue;
            nJets_sameOrder++;
        }
    }
    
    TString histoName;
    
    histoName = "leadingJet_ordering_"+name;
    if(!m_histogram[histoName]) return;
    
    m_histogram[histoName]->Fill(nJets_same, weight);
    if(nJets_sameOrder == jetsId.size() || nJets_sameOrder==2) m_histogram[histoName]->Fill(3, weight);
//     printf("Same: %d  Matched: %d  All: %d\n", nJets_same, nJets_sameOrder, (int)jetsId.size());
}

int AnalyzerDijet::matchRecoToGenJet(const std::vector<int>& jetIndices, const VLV& jets, const int genJetIndex, const VLV& genJets)const
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
    if(deltaRJet>0.4) return -2;
    if(result >= 0){
        const double ptRecoJet = jets.at(result).pt();
        const double ptJet = genJet.pt();
        const double deltaPtRel = (ptJet - ptRecoJet)/ptJet;
        if(deltaPtRel<-0.4 || deltaPtRel>0.6) return -3;
    }
    
    return result;
}

std::vector<std::vector<int> > AnalyzerDijet::matchBhadronsToGenJets(const VLV& allGenJets, const TopGenObjects& topGenObjects)const
{
    std::vector<std::vector<int> > result = std::vector<std::vector<int> >(allGenJets.size());
    
    const std::vector<int>& genBHadJetIndex(*topGenObjects.genBHadJetIndex_);
    for(size_t iHadron = 0; iHadron < genBHadJetIndex.size(); ++iHadron){
        const int& jetIndex = genBHadJetIndex.at(iHadron);
        // Protect against hadrons not clustered to any jet
        if(jetIndex == -1) continue;
        result.at(jetIndex).push_back(iHadron);
    }
    
    return result;
}

