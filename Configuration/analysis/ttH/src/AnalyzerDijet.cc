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


AnalyzerDijet::AnalyzerDijet(const char* mva2dWeightsFile, const std::string& corName, const std::string& swpName,
                             const std::vector<TString>& selectionStepsNoCategories,
                             const std::vector<TString>& stepsForCategories,
                             const JetCategories* jetCategories, bool doHadronMatchingComparison):
AnalyzerBaseClass("dijet_", selectionStepsNoCategories, stepsForCategories, jetCategories),
weightsCorrect_(0),
weightsSwapped_(0),
doHadronMatchingComparison_(doHadronMatchingComparison),
sigJetPt_min_(20.0),
sigJetEta_max_(2.5)
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



void AnalyzerDijet::fillHistos(const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                               const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                               const KinRecoObjects& kinRecoObjects,
                               const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                               const tth::GenLevelWeights&, const tth::RecoLevelWeights& recoLevelWeights,
                               const double& weight, const TString&,
                               std::map<TString, TH1*>& m_histogram)
{
    // Extracting input data to more comfortable variables
    const VLV& allJets = *recoObjects.jets_;
    const LV& lepton = recoObjects.allLeptons_->at(recoObjectIndices.leptonIndex_);
    const LV& antilepton = recoObjects.allLeptons_->at(recoObjectIndices.antiLeptonIndex_);
    const LV& met = *recoObjects.met_;
    const std::vector<int>& jetsId = recoObjectIndices.jetIndices_;           // Selected jets (point to jets from allJets)
    const std::vector<int>& bJetsId = recoObjectIndices.bjetIndices_;         // B-tagged jets (point to jets from allJets)
//     printf("NJETS: %d  NBJETS: %d ###############################\n", (int)jetsId.size(), (int)bJetsId.size());
    std::vector<int> topJetsId;                                               // Jets from ttbar by KinReco (Point to jets from allJets)
    if(kinRecoObjects.valuesSet_) {
        if(kinRecoObjects.HypJet0index_) topJetsId.push_back(kinRecoObjects.HypJet0index_->at(0));
        if(kinRecoObjects.HypJet1index_) topJetsId.push_back(kinRecoObjects.HypJet1index_->at(0));
    }
    const std::vector<double>& allJetsBtagDiscriminant = *recoObjects.jetBTagCSV_;

    // Setting variables of gen. level if available
    const VLV& genJets = (commonGenObjects.valuesSet_) ? *commonGenObjects.allGenJets_ : VLV(0);
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
    checkAdditionalGenBJetAcceptance(topGenObjects, commonGenObjects, m_histogram, 20.0, 2.5, weight);

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
        if(genJets.at(jetId).Pt()<sigJetPt_min_) continue;
        if(genJets.at(jetId).Eta()>sigJetEta_max_) continue;
        if(flavour==6) putUniquelyInVector(genTopJetsId, jetId);
        else if(flavour==25) putUniquelyInVector(genHiggsJetsId, jetId);
    }
    
    if(trueTopAllJetsId.size() == 2 && genTopAllJetsId.size()==2) {
        for(int jetId : trueTopAllJetsId) {
            m_histogram["topBJet_Pt_true"]->Fill(allJets.at(jetId).Pt(), weight);
            m_histogram["topBJet_Btag_true"]->Fill(allJetsBtagDiscriminant.at(jetId), weight);
        }
        for(int jetId : genTopAllJetsId) {
            m_histogram["topBJet_Pt_gen"]->Fill(genJets.at(jetId).Pt(), weight);
        }
    }
    if(trueHiggsAllJetsId.size() == 2 && genHiggsAllJetsId.size()==2) {
        for(int jetId : trueHiggsAllJetsId) {
            m_histogram["higgsBJet_Pt_true"]->Fill(allJets.at(jetId).Pt(), weight);
            m_histogram["higgsBJet_Btag_true"]->Fill(allJetsBtagDiscriminant.at(jetId), weight);
        }
        for(int jetId : genHiggsAllJetsId) {
            m_histogram["higgsBJet_Pt_gen"]->Fill(genJets.at(jetId).Pt(), weight);
        }
    }
    
    m_histogram["topJet_multiplicity_allGen"]->Fill(genTopAllJetsId.size(), weight);
    m_histogram["higgsJet_multiplicity_allGen"]->Fill(genHiggsAllJetsId.size(), weight);
    m_histogram["topJet_multiplicity_gen"]->Fill(genTopJetsId.size(), weight);
    m_histogram["higgsJet_multiplicity_gen"]->Fill(genHiggsJetsId.size(), weight);
    
    
    // Filling the dijet mass of true gen jets from top
    if(genTopJetsId.size()==2) {
        LV dijet_genH = genJets.at(genTopJetsId.at(0)) + genJets.at(genTopJetsId.at(0));
        m_histogram["dijet_mass_genTopJets"]->Fill(dijet_genH.M(), weight);
    }
    // Filling the dijet mass of true gen jets from Higgs
    if(genHiggsJetsId.size() == 2) {
        LV dijet_genH = genJets.at(genHiggsJetsId.at(0)) + genJets.at(genHiggsJetsId.at(1));
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
    unsigned int nGenJets = genJets.size();


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

        // Filling information about b-tagged jets
        if(isInVector(bJetsId, iAllJet)) {
            m_histogram["bJet_pt"]->Fill(jet.Pt(), weight);
            m_histogram["bJet_eta"]->Fill(jet.Eta(), weight);
            m_histogram["bJet_btagDiscriminator"]->Fill(allJetsBtagDiscriminant.at(iAllJet), weight);

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
        int genJetId = genJetIdOfRecoJet(jet, genJets);
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
    m_histogram["weightBTagSF"]->Fill(recoLevelWeights.weightBtagSF_, weight);
    
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

    // DeltaR between leptons, Met
    m_histogram["topLeptonAntilepton_dR"]->Fill(ROOT::Math::VectorUtil::DeltaR(lepton, antilepton));
    LV dilepton = lepton + antilepton;
    m_histogram["topDiLeptonMet_dR"]->Fill(ROOT::Math::VectorUtil::DeltaR(dilepton, met));


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
    for(size_t iJet=0, nJets=genJets.size(); iJet<nJets; iJet++) {
        if(!isInVector(bHadJetIndex, iJet)) continue;       // Skipping if jet is not a true b-jet
        nGenBjets++;
        m_histogram["genBJet_pt"]->Fill(genJets.at(iJet).Pt(), weight);
        m_histogram["genBJet_eta"]->Fill(genJets.at(iJet).Eta(), weight);

        // Checking dR between the closest pair of gen b-jets
        for(size_t iJet2=0; iJet2<iJet; iJet2++) {
            if(!isInVector(bHadJetIndex, iJet2)) continue;
            float dR = ROOT::Math::VectorUtil::DeltaR(genJets.at(iJet), genJets.at(iJet2));
            if(dR<genBJet_dRmin) genBJet_dRmin = dR;
        }       // End of loop over all other gen b-jets

        float genBJet_recoJet_dRmin = 999.9;
        float genBJet_recoBJet_dRmin = 999.9;
        // Checking dR between each gen b-jet and closest reco jet
        for(size_t iJet2 = 0; iJet2<jetsId.size(); iJet2++) {
            unsigned int iAllJet2 = jetsId.at(iJet2);
            LV jet2 = allJets.at(iAllJet2);
            float dR = ROOT::Math::VectorUtil::DeltaR(genJets.at(iJet), jet2);
            if(dR<genBJet_recoJet_dRmin) genBJet_recoJet_dRmin = dR;
            if(!isInVector(bJetsId, iAllJet2)) continue;
            float dR_b = ROOT::Math::VectorUtil::DeltaR(genJets.at(iJet), jet2);
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
                if(genJets.at(genJetId).Pt()>signalJetPt_min && std::fabs(genJets.at(genJetId).Eta())<signalJetEta_max && !isInVector(genTopJetsId, genJetId)) {
                    putUniquelyInVector(genAddBJetIdFromTop_signal, genJetId);
                }
            }
        }
        if(bHadFromTopWeakDecay.size()>iHad && bHadFromTopWeakDecay.at(iHad)==0) {
            nGenAddBhadsNotFromTop++;
            int genJetId = bHadJetIndex.at(iHad);
            if(genJetId>=0 && !isInVector(genTopAllJetsId,genJetId)) {
                putUniquelyInVector(genAddBJetIdNotFromTop, genJetId);
                if(genJets.at(genJetId).Pt()>signalJetPt_min && std::fabs(genJets.at(genJetId).Eta())<signalJetEta_max && !isInVector(genTopJetsId, genJetId)) {
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
        fillGenRecoMatchingComparisonHistos(topGenObjects, higgsGenObjects, bHadLVs, bHadFlavour, bHadJetIndex, genJets, m_histogram, weight);
    }


    ///////////////////////////////////////////////////////////////////////// B-TAGGING DISCRIMINANT VS B-HADRON MULTIPLICITY IN JET
    for(size_t iJet = 0; iJet < jetsId.size(); iJet++) {
        int allJetId = jetsId.at(iJet);
        int genJetId = genJetIdOfRecoJet(allJets.at(allJetId), genJets);
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
    
//     if(trueTopJetsId.size()>1 && trueHiggsJetsId.size()>1) {
    std::vector<std::pair<int, int> > allJetPairs = recoObjectIndices.jetIndexPairs_;
    std::vector<std::pair<int, int> > allBJetPairs;
    for(std::pair<int, int> pair : allJetPairs) {
        if(!isInVector(bJetsId, pair.first) || !isInVector(bJetsId, pair.second)) continue;
        allBJetPairs.push_back(pair);
    }
    fillDijetMassForPairs(allJets, jetsId, trueHiggsJetsId, allBJetPairs, weight, m_histogram, "dijet_mass_all" );
    
    std::vector<std::pair<int, int> > goodJetPairs = jetPairsFromMVA(m_histogram, recoObjectIndices, genObjectIndices, recoObjects, 
                                                                     trueTopJetsId, trueHiggsJetsId, weight);
    fillDijetMassForPairs(allJets, jetsId, trueHiggsJetsId, goodJetPairs, weight, m_histogram, "dijet_mass_goodPairs" );
    fillDijetMassForPairs(allJets, jetsId, trueHiggsJetsId, goodJetPairs, weight, m_histogram, "dijet_mass_goodPairsNormWeight" , true);
    
//     std::vector<std::pair<int, int> > correctJetPairs;
//     const tth::IndexPairs& jetIndexPairs = recoObjectIndices.jetIndexPairs_;
//     for(std::pair<int, int> pair : jetIndexPairs){
//         if(areAmongPairs(wrongJetPairs, pair.first, pair.second)) continue;
//         correctJetPairs.push_back(pair);
//     }
// 
//     // Plotting dijet mass for combinations of jets depending on MVA weight
//     correctPairFraction(allJets, jetsId, bJetsId, allJetsBtagDiscriminant, emptyVector, trueHiggsJetsId, weight, m_histogram["dijet_mass_noWrongMVApairs"], 0, fillAllCombinations, 0.0, 1, wrongJetPairs);
//     correctPairFraction(allJets, jetsId, bJetsId, allJetsBtagDiscriminant, emptyVector, trueHiggsJetsId, weight, m_histogram["dijet_mass_noCorrectMVApairs"], 0, fillAllCombinations, 0.0, 1, correctJetPairs);
    
    
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
    
    const int nBinsX_dijetM = 25;
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
    m_histogram[name] = store(new TH1D(prefix_+name+step, "B-jet btagDiscriminator d;jet d{reco}"+label+";Jets",36,-0.1,1.1));
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
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Gen. b-jet dRmin;bJet dR_min{gen}"+label+";Events",50,0,3.5));
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
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass;N entries/event(M_{dijet}) (all)"+label+";Events",16,0,16));
    name = "dijet_mass_all_diffMin";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass difference_{min};M1_{dijet}-M2_{dijet} (all)"+label+";Events_{N pairs>1}",100,0,500));
    
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

    name = "mvaWeight_correctHiggs_correct";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "MVA^{correct} weight for correct H combinations;weight_{MVA}^{correct}"+label+";Correct jet pairs",20,-1.2,0.2));
    name = "mvaWeight_correctHiggs_swapped";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "MVA^{swapped} weight for correct H combinations;weight_{MVA}^{swapped}"+label+";Correct jet pairs",20,-1.2,0.2));
    
    name="mvaWeight_correctTop_correctHiggs";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "MVA_{top}^{higgs} weight for correct combinations;MVA weight_{correct}^{top}"+label+";MVA weight_{correct}^{higgs}",20,-1.2,0.2,20,-1.2,0.2));
    name = "mvaWeight_wrongTopAndHiggs_correct";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "MVA^{correct} weight for combinations not from Top, Higgs;weight_{MVA}^{correct}"+label+";Wrong jet pairs",20,-1.2,0.2));

    name = "dijet_mass_goodPairs";
    bookPairHistos(new TH1D(prefix_+name+step, "Dijet mass (no wrong pairs from MVA);M(dijet) (MVA_{correct}^{high})"+label+";Jet pairs",nBinsX_dijetM,binsX_dijetM), m_histogram, name);
    name = "dijet_mass_goodPairs_nEntriesPerEvent";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass (now wrong pairs from MVA);N entries/event(M_{dijet}) (MVA_{correct}^{high})"+label+";Events",16,0,16));
    name = "dijet_mass_goodPairs_diffMin";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Dijet mass difference_{min} (no wrong pairs from MVA);M1_{dijet}-M2_{dijet} (MVA_{correct}^{high})"+label+";Events_{N pairs>1}",100,0,500));
    name = "dijet_mass_goodPairs_nEntriesPerBin";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "Dijet mass difference_{min} (no wrong pairs from MVA);M_{dijet}(MVA_{correct}^{high})"+label+";Entries",nBinsX_dijetM,binsX_dijetM, 10, 0, 10));
    
    name = "dijet_mass_goodPairsNormWeight";
    bookPairHistos(new TH1D(prefix_+name+step, "Dijet mass (no correct pairs from MVA);M(dijet) (MVA_{correct}^{high})"+label+";Jet pairs",nBinsX_dijetM,binsX_dijetM), m_histogram, name);

    name = "weight";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "Weight of the event;weight"+label+";Events",30,0,3));
    name = "weightBTagSF";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "B-tag weight of the event;weight"+label+";Events",60,0.8,1.1));
    
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
//             printf("higgs\n");
        }
        else {
            m_histogram["mvaWeight_wrongTopAndHiggs_correct"]->Fill(v_mvaWeightsCorrect.at(i), weight);
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


void AnalyzerDijet::checkAdditionalGenBJetAcceptance(const TopGenObjects& topGenObjects, const CommonGenObjects& commonGenObjects, 
                                                     std::map<TString, TH1*>& m_histogram, const float jetPt_min, const float jetEta_max, const double weight) 
{
    // Setting variables of gen. level if available
    const VLV& genJets = (commonGenObjects.valuesSet_) ? *commonGenObjects.allGenJets_ : VLV(0);
    const std::vector<int>& bHadJetIndex = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadJetIndex_ : std::vector<int>(0);
    const std::vector<int>& bHadFlavour = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadFlavour_ : std::vector<int>(0);
    const std::vector<int>& bHadFromTopWeakDecay = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadFromTopWeakDecay_ : std::vector<int>(0);
    
    std::vector<int> addHadronIds;
    std::vector<int> addJetIds;
    
    std::vector<int> topJetIds;
    // Creating the list of gen b-jets coming from top
    for(size_t hadId = 0; hadId < bHadFlavour.size(); ++hadId) {
        if(std::abs(bHadFlavour.at(hadId))!=6) continue;
        if(bHadJetIndex.at(hadId)<0) continue;
        putUniquelyInVector(topJetIds, bHadJetIndex.at(hadId));
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
        if(genJets.at(jetId).Pt()<jetPt_min || std::fabs(genJets.at(jetId).Eta())>jetEta_max) continue;
        m_histogram["extraGenBHadronsJets"]->Fill(2., weight);
        // Jet overlaps with b-jets from tt
        if(isInVector(topJetIds, jetId)) continue;
        m_histogram["extraGenBHadronsJets"]->Fill(3., weight);
        // Jet overlaps with extra b-jets
        if(isInVector(addJetIds, jetId)) continue;
        m_histogram["extraGenBHadronsJets"]->Fill(4., weight);
        
        putUniquelyInVector(addJetIds, jetId);
    }       // End of loop over all b-hadrons
    
    m_histogram["extraGenBHadronMultiplicity"]->Fill((int)addHadronIds.size(), weight);
    m_histogram["extraGenBJetMultiplicity"]->Fill((int)addJetIds.size(), weight);
    
    for(int iJet1 : bHadJetIndex) {
        if(iJet1<0) continue;
        if(genJets.at(iJet1).Pt()<sigJetPt_min_) continue;
        if(std::fabs(genJets.at(iJet1).Eta())>sigJetEta_max_) continue;
        
        for(int iJet2 : bHadJetIndex) {
            if(iJet2<0) continue;
            if(iJet2==iJet1) continue;
            if(genJets.at(iJet2).Pt()<sigJetPt_min_) continue;
            if(std::fabs(genJets.at(iJet2).Eta())>sigJetEta_max_) continue;
            LV dijet = genJets.at(iJet1) + genJets.at(iJet2);
            m_histogram["dijet_mass_genAllBJets"]->Fill(dijet.M(), weight);
            if(!isInVector(addJetIds, iJet1)) continue;
            if(!isInVector(addJetIds, iJet2)) continue;
            m_histogram["dijet_mass_genAdditionalBJets"]->Fill(dijet.M(), weight);
        }
    }
    
}


void AnalyzerDijet::fillDijetMassForPairs(const VLV& allJets, const std::vector<int>& jetsId, const std::vector<int>& higgsJetsId,
                                          const std::vector<std::pair<int, int> > &jetPairs, const double weight, 
                                          std::map<TString, TH1*>& m_histogram, std::string histoName, const bool normaliseWeight )
{
    int nPairs = jetPairs.size();
    int nEntries = 0;
    float diffMin = 999.9;
    TH1D* h_nEntriesPerBin = 0;
    if(m_histogram[histoName+"_nEntriesPerBin"]) {
    int nBins = m_histogram[histoName+"_nEntriesPerBin"]->GetNbinsX();
    h_nEntriesPerBin = new TH1D("h_nEntriesPerBin_temp","",nBins,m_histogram[histoName+"_nEntriesPerBin"]->GetXaxis()->GetXbins()->GetArray());
    }
    
    for(unsigned int i = 0; i<jetPairs.size(); ++i) {
        const std::pair<int,int> &jetPair = jetPairs.at(i);
        std::pair<int,int> allJetPair(jetsId.at(jetPair.first), jetsId.at(jetPair.second));
        bool isCorrect = isInVector(higgsJetsId, allJetPair.first) && isInVector(higgsJetsId, allJetPair.second);
        
        LV dijet = allJets.at(allJetPair.first) + allJets.at(allJetPair.second);
        double fillWeight = normaliseWeight ? weight/double(nPairs) : weight;
        fillPairHistos(m_histogram, histoName, dijet.M(), isCorrect, fillWeight);
        if(h_nEntriesPerBin)h_nEntriesPerBin->Fill(dijet.M());
        nEntries++;
        for(unsigned int j = i+1; j<jetPairs.size(); ++j) {
            const std::pair<int,int> &jetPair = jetPairs.at(j);
            std::pair<int,int> allJetPair(jetsId.at(jetPair.first), jetsId.at(jetPair.second));
            
            LV dijet2 = allJets.at(allJetPair.first) + allJets.at(allJetPair.second);
            
            float diff = std::fabs(dijet.M() - dijet2.M());
            if(diff<diffMin) diffMin = diff;
        }
    }
    
    if(m_histogram[histoName+"_nEntriesPerEvent"]) m_histogram[histoName+"_nEntriesPerEvent"]->Fill(nEntries, weight);
    if(m_histogram[histoName+"_diffMin"] && diffMin < 999.) m_histogram[histoName+"_diffMin"]->Fill(diffMin, weight);
    
    if(m_histogram[histoName+"_nEntriesPerBin"]) {
        for(int iBin=0; iBin<h_nEntriesPerBin->GetNbinsX(); ++iBin) {
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
