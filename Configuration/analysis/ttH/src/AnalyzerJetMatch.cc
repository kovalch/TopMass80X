#include <map>
#include <vector>
#include <iostream>
#include <algorithm>

#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <Math/VectorUtil.h>

#include "AnalyzerJetMatch.h"
#include "analysisStructs.h"
#include "JetCategories.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"








AnalyzerJetMatch::AnalyzerJetMatch(const std::vector<TString>& selectionStepsNoCategories,
                                   const std::vector<TString>& stepsForCategories,
                                   const JetCategories* jetCategories):
AnalyzerBase("jetMatch_", selectionStepsNoCategories, stepsForCategories, jetCategories)
{
    std::cout<<"--- Beginning setting up jet match analyzer\n";
    std::cout<<"=== Finishing setting up jet match analyzer\n\n";
}



void AnalyzerJetMatch::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    this->bookJetHistos("initial_", step, m_histogram);
    this->bookJetHistos("initialAmbiguous_", step, m_histogram);
    this->bookJetHistos("initialUnambiguous_", step, m_histogram);
    this->bookJetHistos("mismatchedInR_", step, m_histogram);
    this->bookJetHistos("matchedInR_", step, m_histogram);
    this->bookJetHistos("matchedInRAmbiguous_", step, m_histogram);
    this->bookJetHistos("matchedInRUnambiguous_", step, m_histogram);
    
    this->bookJetHistos("mismatchedInPt_m05_p07_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguousMatchedInR_m05_p07_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguousMatchedInR_m05_p07_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguous_m05_p07_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguous_m05_p07_", step, m_histogram);
    this->bookJetHistos("matchedInPt_m05_p07_", step, m_histogram);
    this->bookJetHistos("matchedInPtAmbiguous_m05_p07_", step, m_histogram);
    this->bookJetHistos("matchedInPtUnambiguous_m05_p07_", step, m_histogram);
    
    this->bookJetHistos("mismatchedInPt_m05_p06_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguousMatchedInR_m05_p06_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguousMatchedInR_m05_p06_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguous_m05_p06_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguous_m05_p06_", step, m_histogram);
    this->bookJetHistos("matchedInPt_m05_p06_", step, m_histogram);
    this->bookJetHistos("matchedInPtAmbiguous_m05_p06_", step, m_histogram);
    this->bookJetHistos("matchedInPtUnambiguous_m05_p06_", step, m_histogram);
    
    this->bookJetHistos("mismatchedInPt_m05_p05_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguousMatchedInR_m05_p05_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguousMatchedInR_m05_p05_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguous_m05_p05_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguous_m05_p05_", step, m_histogram);
    this->bookJetHistos("matchedInPt_m05_p05_", step, m_histogram);
    this->bookJetHistos("matchedInPtAmbiguous_m05_p05_", step, m_histogram);
    this->bookJetHistos("matchedInPtUnambiguous_m05_p05_", step, m_histogram);
    
    this->bookJetHistos("mismatchedInPt_m04_p07_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguousMatchedInR_m04_p07_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguousMatchedInR_m04_p07_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguous_m04_p07_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguous_m04_p07_", step, m_histogram);
    this->bookJetHistos("matchedInPt_m04_p07_", step, m_histogram);
    this->bookJetHistos("matchedInPtAmbiguous_m04_p07_", step, m_histogram);
    this->bookJetHistos("matchedInPtUnambiguous_m04_p07_", step, m_histogram);
    
    this->bookJetHistos("mismatchedInPt_m04_p06_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguousMatchedInR_m04_p06_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguousMatchedInR_m04_p06_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguous_m04_p06_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguous_m04_p06_", step, m_histogram);
    this->bookJetHistos("matchedInPt_m04_p06_", step, m_histogram);
    this->bookJetHistos("matchedInPtAmbiguous_m04_p06_", step, m_histogram);
    this->bookJetHistos("matchedInPtUnambiguous_m04_p06_", step, m_histogram);
    
    this->bookJetHistos("mismatchedInPt_m04_p05_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguousMatchedInR_m04_p05_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguousMatchedInR_m04_p05_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguous_m04_p05_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguous_m04_p05_", step, m_histogram);
    this->bookJetHistos("matchedInPt_m04_p05_", step, m_histogram);
    this->bookJetHistos("matchedInPtAmbiguous_m04_p05_", step, m_histogram);
    this->bookJetHistos("matchedInPtUnambiguous_m04_p05_", step, m_histogram);
    
    this->bookJetHistos("mismatchedInPt_m03_p07_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguousMatchedInR_m03_p07_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguousMatchedInR_m03_p07_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguous_m03_p07_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguous_m03_p07_", step, m_histogram);
    this->bookJetHistos("matchedInPt_m03_p07_", step, m_histogram);
    this->bookJetHistos("matchedInPtAmbiguous_m03_p07_", step, m_histogram);
    this->bookJetHistos("matchedInPtUnambiguous_m03_p07_", step, m_histogram);
    
    this->bookJetHistos("mismatchedInPt_m03_p06_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguousMatchedInR_m03_p06_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguousMatchedInR_m03_p06_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguous_m03_p06_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguous_m03_p06_", step, m_histogram);
    this->bookJetHistos("matchedInPt_m03_p06_", step, m_histogram);
    this->bookJetHistos("matchedInPtAmbiguous_m03_p06_", step, m_histogram);
    this->bookJetHistos("matchedInPtUnambiguous_m03_p06_", step, m_histogram);
    
    this->bookJetHistos("mismatchedInPt_m03_p05_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguousMatchedInR_m03_p05_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguousMatchedInR_m03_p05_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtAmbiguous_m03_p05_", step, m_histogram);
    this->bookJetHistos("mismatchedInPtUnambiguous_m03_p05_", step, m_histogram);
    this->bookJetHistos("matchedInPt_m03_p05_", step, m_histogram);
    this->bookJetHistos("matchedInPtAmbiguous_m03_p05_", step, m_histogram);
    this->bookJetHistos("matchedInPtUnambiguous_m03_p05_", step, m_histogram);
}



void AnalyzerJetMatch::bookJetHistos(const TString& whichSelection, const TString& step, std::map<TString, TH1*>& m_histogram)
{
    this->bookJetHistosInclExcl(whichSelection, "all_", step, m_histogram);
    this->bookJetHistosInclExcl(whichSelection, "topORhiggs_", step, m_histogram);
    this->bookJetHistosInclExcl(whichSelection, "top_", step, m_histogram);
    this->bookJetHistosInclExcl(whichSelection, "higgs_", step, m_histogram);
}



void AnalyzerJetMatch::bookJetHistosInclExcl(const TString& whichSelection, const TString& whichJets,
                                             const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    const TString identifier = whichSelection + whichJets;
    
    name = identifier + "deltaR";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, name+";#DeltaR_{gen,reco};Entries",50,0,5));

    name = identifier + "deltaR_zoom";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, name+";#DeltaR_{gen,reco};Entries",50,0,1));

    name = identifier + "recoPtOverGenPt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, name+";p_{T,reco}/p_{T,gen};Entries",30,0,3));
    
    name = identifier + "deltaPtOverGenPt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, name+";(p_{T,gen}-p_{T,reco})/p_{T,gen};Entries",120,-1.2,1.2));
    
    name = identifier + "genPtVsRecoPt";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step, name+";p_{T,reco} [GeV];p_{T,gen} [GeV]",50,0,300,50,0,300));

    name = identifier + "genEtaVsRecoEta";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step, name+";#eta_{reco};#eta_{gen}",50,-2.6,2.6,50,-2.6,2.6));

    name = identifier + "genPhiVsRecoPhi";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step, name+";#phi_{reco};#phi_{gen}",40,-4,4,40,-4,4));

    name = identifier + "genPVsRecoP";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step, name+";|p_{reco}| [GeV];|p_{gen}| [GeV]",30,0,150,30,0,150));

    name = identifier + "genPtVsDeltaPt";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step, name+";p_{T,gen}-p_{T,reco} [GeV];p_{T,gen} [GeV]",30,-100,200,50,0,300));

    name = identifier + "deltaRVsGenPt";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step, name+";p_{T,gen} [GeV];#DeltaR_{gen,reco}",50,0,300,260,0,2)); 
    
    name = identifier + "recoPtOverGenPtVsGenPt";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step, name+";p_{T,gen} [GeV];p_{T,reco}/p_{T,gen}",50,0,300,30,0,3));

    name = identifier + "recoPtOverGenPtVsGenEta";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step, name+";#eta_{gen};p_{T,reco}/p_{T,gen}",50,-2.6,2.6,30,0,3));

    name = identifier + "recoPtOverGenPtVsGenPhi";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step, name+";#phi_{gen};p_{T,reco}/p_{T,gen}",40,-4,4,30,0,3));

    name = identifier + "recoPtOverGenPtVsDeltaR";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step, name+";#DeltaR_{gen,reco};p_{T,reco}/p_{T,gen}",250,0,5,30,0,3));

    name = identifier + "deltaPtOverGenPtVsDeltaR";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step, name+";#DeltaR_{gen,reco};(p_{T,gen}-p_{T,reco})/p_{T,gen}",250,0,5,120,-1.2,1.2));
}



void AnalyzerJetMatch::fillHistos(const RecoObjects& recoObjects, const CommonGenObjects&,
                                  const TopGenObjects& topGenObjects, const HiggsGenObjects&,
                                  const KinRecoObjects&,
                                  const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                                  const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                  const double& weight, const TString&,
                                  std::map<TString, TH1*>& m_histogram)
{
    if(!topGenObjects.valuesSet_) return;
    
    //For each gen jet find the reco jet with the minimum deltaR, and store values for selection criteria in R and pt
    std::vector<int> selectedGenIndices;
    std::vector<int> closestRecoJetIndices;
    std::vector<double> v_deltaR;
    std::vector<double> v_deltaPtRel;
    for(const int genIndex : common::initialiseIndices(*topGenObjects.allGenJets_)){
        const LV& genJet = topGenObjects.allGenJets_->at(genIndex);
        double minDeltaR(999.);
        double deltaPtRel(999.);
        int recoJetIndex(-1);
        for(const int recoIndex : recoObjectIndices.jetIndices_){
            const LV& recoJet = recoObjects.jets_->at(recoIndex);
            const double deltaR = ROOT::Math::VectorUtil::DeltaR(recoJet, genJet);
            if(deltaR < minDeltaR){
                minDeltaR = deltaR;
                recoJetIndex = recoIndex;
                const double& genPt = genJet.pt();
                const double& recoPt = recoJet.pt();
                const double deltaPt = genPt - recoPt;
                deltaPtRel = deltaPt/genPt;
            }
        }
        // The values always need to be set, such that the vectors have the same size as the genJets
        closestRecoJetIndices.push_back(recoJetIndex);
        v_deltaR.push_back(minDeltaR);
        v_deltaPtRel.push_back(deltaPtRel);
        
        // For selectedGenIndices, store only those where a recoJet is found
        // Should be always true except in events without any selected recoJet
        if(recoJetIndex != -1) selectedGenIndices.push_back(genIndex);
    }
    
    // Build collections for all selections 
    std::vector<int> genIndices_mismatchedInR = selectedGenIndices;
    common::selectIndices(genIndices_mismatchedInR, v_deltaR, 0.5);
    std::vector<int> genIndices_matchedInR = selectedGenIndices;
    common::selectIndices(genIndices_matchedInR, v_deltaR, 0.4, false);
    
    std::vector<int> genIndices_mismatchedInPt_up = genIndices_matchedInR;
    std::vector<int> genIndices_mismatchedInPt_down = genIndices_matchedInR;
    common::selectIndices(genIndices_mismatchedInPt_up, v_deltaPtRel, 0.7);
    common::selectIndices(genIndices_mismatchedInPt_down, v_deltaPtRel, -0.5, false);
    const std::vector<int> genIndices_mismatchedInPt_m05_p07 = common::mergeIndices(genIndices_mismatchedInPt_up, genIndices_mismatchedInPt_down);
    std::vector<int> genIndices_matchedInPt_m05_p07 = genIndices_matchedInR;
    common::selectIndices(genIndices_matchedInPt_m05_p07, v_deltaPtRel, 0.7, false);
    common::selectIndices(genIndices_matchedInPt_m05_p07, v_deltaPtRel, -0.5);
    genIndices_mismatchedInPt_up.clear(); genIndices_mismatchedInPt_down.clear();
    
    genIndices_mismatchedInPt_up = genIndices_matchedInR;
    genIndices_mismatchedInPt_down = genIndices_matchedInR;
    common::selectIndices(genIndices_mismatchedInPt_up, v_deltaPtRel, 0.6);
    common::selectIndices(genIndices_mismatchedInPt_down, v_deltaPtRel, -0.5, false);
    const std::vector<int> genIndices_mismatchedInPt_m05_p06 = common::mergeIndices(genIndices_mismatchedInPt_up, genIndices_mismatchedInPt_down);
    std::vector<int> genIndices_matchedInPt_m05_p06 = genIndices_matchedInR;
    common::selectIndices(genIndices_matchedInPt_m05_p06, v_deltaPtRel, 0.6, false);
    common::selectIndices(genIndices_matchedInPt_m05_p06, v_deltaPtRel, -0.5);
    genIndices_mismatchedInPt_up.clear(); genIndices_mismatchedInPt_down.clear();
    
    genIndices_mismatchedInPt_up = genIndices_matchedInR;
    genIndices_mismatchedInPt_down = genIndices_matchedInR;
    common::selectIndices(genIndices_mismatchedInPt_up, v_deltaPtRel, 0.5);
    common::selectIndices(genIndices_mismatchedInPt_down, v_deltaPtRel, -0.5, false);
    const std::vector<int> genIndices_mismatchedInPt_m05_p05 = common::mergeIndices(genIndices_mismatchedInPt_up, genIndices_mismatchedInPt_down);
    std::vector<int> genIndices_matchedInPt_m05_p05 = genIndices_matchedInR;
    common::selectIndices(genIndices_matchedInPt_m05_p05, v_deltaPtRel, 0.5, false);
    common::selectIndices(genIndices_matchedInPt_m05_p05, v_deltaPtRel, -0.5);
    genIndices_mismatchedInPt_up.clear(); genIndices_mismatchedInPt_down.clear();
    
    genIndices_mismatchedInPt_up = genIndices_matchedInR;
    genIndices_mismatchedInPt_down = genIndices_matchedInR;
    common::selectIndices(genIndices_mismatchedInPt_up, v_deltaPtRel, 0.7);
    common::selectIndices(genIndices_mismatchedInPt_down, v_deltaPtRel, -0.4, false);
    const std::vector<int> genIndices_mismatchedInPt_m04_p07 = common::mergeIndices(genIndices_mismatchedInPt_up, genIndices_mismatchedInPt_down);
    std::vector<int> genIndices_matchedInPt_m04_p07 = genIndices_matchedInR;
    common::selectIndices(genIndices_matchedInPt_m04_p07, v_deltaPtRel, 0.7, false);
    common::selectIndices(genIndices_matchedInPt_m04_p07, v_deltaPtRel, -0.4);
    genIndices_mismatchedInPt_up.clear(); genIndices_mismatchedInPt_down.clear();
    
    genIndices_mismatchedInPt_up = genIndices_matchedInR;
    genIndices_mismatchedInPt_down = genIndices_matchedInR;
    common::selectIndices(genIndices_mismatchedInPt_up, v_deltaPtRel, 0.6);
    common::selectIndices(genIndices_mismatchedInPt_down, v_deltaPtRel, -0.4, false);
    const std::vector<int> genIndices_mismatchedInPt_m04_p06 = common::mergeIndices(genIndices_mismatchedInPt_up, genIndices_mismatchedInPt_down);
    std::vector<int> genIndices_matchedInPt_m04_p06 = genIndices_matchedInR;
    common::selectIndices(genIndices_matchedInPt_m04_p06, v_deltaPtRel, 0.6, false);
    common::selectIndices(genIndices_matchedInPt_m04_p06, v_deltaPtRel, -0.4);
    genIndices_mismatchedInPt_up.clear(); genIndices_mismatchedInPt_down.clear();
    
    genIndices_mismatchedInPt_up = genIndices_matchedInR;
    genIndices_mismatchedInPt_down = genIndices_matchedInR;
    common::selectIndices(genIndices_mismatchedInPt_up, v_deltaPtRel, 0.5);
    common::selectIndices(genIndices_mismatchedInPt_down, v_deltaPtRel, -0.4, false);
    const std::vector<int> genIndices_mismatchedInPt_m04_p05 = common::mergeIndices(genIndices_mismatchedInPt_up, genIndices_mismatchedInPt_down);
    std::vector<int> genIndices_matchedInPt_m04_p05 = genIndices_matchedInR;
    common::selectIndices(genIndices_matchedInPt_m04_p05, v_deltaPtRel, 0.5, false);
    common::selectIndices(genIndices_matchedInPt_m04_p05, v_deltaPtRel, -0.4);
    genIndices_mismatchedInPt_up.clear(); genIndices_mismatchedInPt_down.clear();
    
    genIndices_mismatchedInPt_up = genIndices_matchedInR;
    genIndices_mismatchedInPt_down = genIndices_matchedInR;
    common::selectIndices(genIndices_mismatchedInPt_up, v_deltaPtRel, 0.7);
    common::selectIndices(genIndices_mismatchedInPt_down, v_deltaPtRel, -0.3, false);
    const std::vector<int> genIndices_mismatchedInPt_m03_p07 = common::mergeIndices(genIndices_mismatchedInPt_up, genIndices_mismatchedInPt_down);
    std::vector<int> genIndices_matchedInPt_m03_p07 = genIndices_matchedInR;
    common::selectIndices(genIndices_matchedInPt_m03_p07, v_deltaPtRel, 0.7, false);
    common::selectIndices(genIndices_matchedInPt_m03_p07, v_deltaPtRel, -0.3);
    genIndices_mismatchedInPt_up.clear(); genIndices_mismatchedInPt_down.clear();
    
    genIndices_mismatchedInPt_up = genIndices_matchedInR;
    genIndices_mismatchedInPt_down = genIndices_matchedInR;
    common::selectIndices(genIndices_mismatchedInPt_up, v_deltaPtRel, 0.6);
    common::selectIndices(genIndices_mismatchedInPt_down, v_deltaPtRel, -0.3, false);
    const std::vector<int> genIndices_mismatchedInPt_m03_p06 = common::mergeIndices(genIndices_mismatchedInPt_up, genIndices_mismatchedInPt_down);
    std::vector<int> genIndices_matchedInPt_m03_p06 = genIndices_matchedInR;
    common::selectIndices(genIndices_matchedInPt_m03_p06, v_deltaPtRel, 0.6, false);
    common::selectIndices(genIndices_matchedInPt_m03_p06, v_deltaPtRel, -0.3);
    genIndices_mismatchedInPt_up.clear(); genIndices_mismatchedInPt_down.clear();
    
    genIndices_mismatchedInPt_up = genIndices_matchedInR;
    genIndices_mismatchedInPt_down = genIndices_matchedInR;
    common::selectIndices(genIndices_mismatchedInPt_up, v_deltaPtRel, 0.5);
    common::selectIndices(genIndices_mismatchedInPt_down, v_deltaPtRel, -0.3, false);
    const std::vector<int> genIndices_mismatchedInPt_m03_p05 = common::mergeIndices(genIndices_mismatchedInPt_up, genIndices_mismatchedInPt_down);
    std::vector<int> genIndices_matchedInPt_m03_p05 = genIndices_matchedInR;
    common::selectIndices(genIndices_matchedInPt_m03_p05, v_deltaPtRel, 0.5, false);
    common::selectIndices(genIndices_matchedInPt_m03_p05, v_deltaPtRel, -0.3);    
    genIndices_mismatchedInPt_up.clear(); genIndices_mismatchedInPt_down.clear();
    
    // Loop over jets and fill the histograms for all selections
    for(const int genIndex : selectedGenIndices){
        const int recoIndex = closestRecoJetIndices.at(genIndex);
        
        // Fill histograms for jets without any selection, for
        // 1) inclusively
        // 2) for jets with ambiguous matching (i.e. another genJet is matched to same recoJet - both of the genJets should be compatible with all the applied cuts)
        // 3) for the unambiguous cases
        this->fillJetHistos("initial_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        if(this->isAmbiguous(genIndex, genIndices_matchedInR, closestRecoJetIndices))
            this->fillJetHistos("initialAmbiguous_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        else
            this->fillJetHistos("initialUnambiguous_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);

        // Fill histograms for jets mismatched in R
        if(std::find(genIndices_mismatchedInR.begin(), genIndices_mismatchedInR.end(), genIndex) != genIndices_mismatchedInR.end())
            this->fillJetHistos("mismatchedInR_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        
        // Select jets matched in R
        if(std::find(genIndices_matchedInR.begin(), genIndices_matchedInR.end(), genIndex) == genIndices_matchedInR.end()) continue;
        
        // Fill histograms for jets matched in R, for
        // 1) inclusively
        // 2) for jets with ambiguous matching (i.e. another genJet is matched to same recoJet - both of the genJets should be compatible with all the applied cuts)
        // 3) for the unambiguous cases
        this->fillJetHistos("matchedInR_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        if(this->isAmbiguous(genIndex, genIndices_matchedInR, closestRecoJetIndices))
            this->fillJetHistos("matchedInRAmbiguous_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        else
            this->fillJetHistos("matchedInRUnambiguous_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        
        
        // Select jets mismatched in pt
        // Fill histograms for jets mismatched in pt, same separation of plots as for the matched in R
        // Additional plots for jets with ambiguousMatchedInR matching (i.e another genJet is matched to the same recoJet, if at least one of these genJets is compatible with all the applied cuts)
  
        // -0.5< pt < 0.7    
        if(std::find(genIndices_mismatchedInPt_m05_p07.begin(), genIndices_mismatchedInPt_m05_p07.end(), genIndex) != genIndices_mismatchedInPt_m05_p07.end()){
            this->fillJetHistos("mismatchedInPt_m05_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            if(this->isAmbiguous(genIndex, genIndices_mismatchedInPt_m05_p07, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguous_m05_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguous_m05_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            
            if(this->isAmbiguous(genIndex, genIndices_matchedInR, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguousMatchedInR_m05_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguousMatchedInR_m05_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        }
         
        // -0.5< pt < 0.6
        if(std::find(genIndices_mismatchedInPt_m05_p06.begin(), genIndices_mismatchedInPt_m05_p06.end(), genIndex) != genIndices_mismatchedInPt_m05_p06.end()){
            this->fillJetHistos("mismatchedInPt_m05_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            if(this->isAmbiguous(genIndex, genIndices_mismatchedInPt_m05_p06, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguous_m05_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguous_m05_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            
            if(this->isAmbiguous(genIndex, genIndices_matchedInR, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguousMatchedInR_m05_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguousMatchedInR_m05_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        }
        
        // -0.5< pt < 0.5
        if(std::find(genIndices_mismatchedInPt_m05_p05.begin(), genIndices_mismatchedInPt_m05_p05.end(), genIndex) != genIndices_mismatchedInPt_m05_p05.end()){
            this->fillJetHistos("mismatchedInPt_m05_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            if(this->isAmbiguous(genIndex, genIndices_mismatchedInPt_m05_p05, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguous_m05_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguous_m05_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            
            if(this->isAmbiguous(genIndex, genIndices_matchedInR, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguousMatchedInR_m05_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguousMatchedInR_m05_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        }
        
        // -0.4< pt < 0.7
        if(std::find(genIndices_mismatchedInPt_m04_p07.begin(), genIndices_mismatchedInPt_m04_p07.end(), genIndex) != genIndices_mismatchedInPt_m04_p07.end()){
            this->fillJetHistos("mismatchedInPt_m04_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            if(this->isAmbiguous(genIndex, genIndices_mismatchedInPt_m04_p07, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguous_m04_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguous_m04_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            
            if(this->isAmbiguous(genIndex, genIndices_matchedInR, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguousMatchedInR_m04_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguousMatchedInR_m04_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        }
        
        //-0.4< pt < 0.6
        if(std::find(genIndices_mismatchedInPt_m04_p06.begin(), genIndices_mismatchedInPt_m04_p06.end(), genIndex) != genIndices_mismatchedInPt_m04_p06.end()){
            this->fillJetHistos("mismatchedInPt_m04_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            if(this->isAmbiguous(genIndex, genIndices_mismatchedInPt_m04_p06, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguous_m04_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguous_m04_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            
            if(this->isAmbiguous(genIndex, genIndices_matchedInR, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguousMatchedInR_m04_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguousMatchedInR_m04_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        }
        
        // -0.4< pt < 0.5
        if(std::find(genIndices_mismatchedInPt_m04_p05.begin(), genIndices_mismatchedInPt_m04_p05.end(), genIndex) != genIndices_mismatchedInPt_m04_p05.end()){
            this->fillJetHistos("mismatchedInPt_m04_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            if(this->isAmbiguous(genIndex, genIndices_mismatchedInPt_m04_p05, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguous_m04_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguous_m04_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            
            if(this->isAmbiguous(genIndex, genIndices_matchedInR, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguousMatchedInR_m04_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguousMatchedInR_m04_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        }
        
        // -0.3< pt < 0.7
        if(std::find(genIndices_mismatchedInPt_m03_p07.begin(), genIndices_mismatchedInPt_m03_p07.end(), genIndex) != genIndices_mismatchedInPt_m03_p07.end()){
            this->fillJetHistos("mismatchedInPt_m03_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            if(this->isAmbiguous(genIndex, genIndices_mismatchedInPt_m03_p07, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguous_m03_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguous_m03_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            
            if(this->isAmbiguous(genIndex, genIndices_matchedInR, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguousMatchedInR_m03_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguousMatchedInR_m03_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        }
        
        // -0.3< pt < 0.6
        if(std::find(genIndices_mismatchedInPt_m03_p06.begin(), genIndices_mismatchedInPt_m03_p06.end(), genIndex) != genIndices_mismatchedInPt_m03_p06.end()){
            this->fillJetHistos("mismatchedInPt_m03_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            if(this->isAmbiguous(genIndex, genIndices_mismatchedInPt_m03_p06, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguous_m03_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguous_m03_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            
            if(this->isAmbiguous(genIndex, genIndices_matchedInR, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguousMatchedInR_m03_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguousMatchedInR_m03_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        }
        
        // -0.3< pt < 0.5
        if(std::find(genIndices_mismatchedInPt_m03_p05.begin(), genIndices_mismatchedInPt_m03_p05.end(), genIndex) != genIndices_mismatchedInPt_m03_p05.end()){
            this->fillJetHistos("mismatchedInPt_m03_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            if(this->isAmbiguous(genIndex, genIndices_mismatchedInPt_m03_p05, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguous_m03_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguous_m03_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            
            if(this->isAmbiguous(genIndex, genIndices_matchedInR, closestRecoJetIndices))
                this->fillJetHistos("mismatchedInPtAmbiguousMatchedInR_m03_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("mismatchedInPtUnambiguousMatchedInR_m03_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        }
        
        
        // Select jets matched in pt
        // Fill histograms for jets matched in pt, same separation of plots as for the matched in R
        
        // -0.5< pt < 0.7
        if(std::find(genIndices_matchedInPt_m05_p07.begin(), genIndices_matchedInPt_m05_p07.end(), genIndex) == genIndices_matchedInPt_m05_p07.end()) continue;
        
        // Fill histograms for jets matched in pt, same separation of plots as for the mismatched in pt
        this->fillJetHistos("matchedInPt_m05_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        if(this->isAmbiguous(genIndex, genIndices_matchedInPt_m05_p07, closestRecoJetIndices))
            this->fillJetHistos("matchedInPtAmbiguous_m05_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        else
            this->fillJetHistos("matchedInPtUnambiguous_m05_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        
        
        // -0.5< pt < 0.6
        if(std::find(genIndices_matchedInPt_m05_p06.begin(), genIndices_matchedInPt_m05_p06.end(), genIndex) != genIndices_matchedInPt_m05_p06.end()){
            // Fill histograms for jets matched in pt, same separation of plots as for the matched in R
            this->fillJetHistos("matchedInPt_m05_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            if(this->isAmbiguous(genIndex, genIndices_matchedInPt_m05_p06, closestRecoJetIndices))
                this->fillJetHistos("matchedInPtAmbiguous_m05_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("matchedInPtUnambiguous_m05_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        }
        
        // -0.5< pt < 0.5
        if(std::find(genIndices_matchedInPt_m05_p05.begin(), genIndices_matchedInPt_m05_p05.end(), genIndex) != genIndices_matchedInPt_m05_p05.end()){
            // Fill histograms for jets matched in pt, same separation of plots as for the matched in R
            this->fillJetHistos("matchedInPt_m05_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            if(this->isAmbiguous(genIndex, genIndices_matchedInPt_m05_p05, closestRecoJetIndices))
                this->fillJetHistos("matchedInPtAmbiguous_m05_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("matchedInPtUnambiguous_m05_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        }
        
        // -0.4< pt < 0.7
        if(std::find(genIndices_matchedInPt_m04_p07.begin(), genIndices_matchedInPt_m04_p07.end(), genIndex) != genIndices_matchedInPt_m04_p07.end()){
            // Fill histograms for jets matched in pt, same separation of plots as for the matched in R
            this->fillJetHistos("matchedInPt_m04_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            if(this->isAmbiguous(genIndex, genIndices_matchedInPt_m04_p07, closestRecoJetIndices))
                this->fillJetHistos("matchedInPtAmbiguous_m04_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("matchedInPtUnambiguous_m04_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        }
        
        //-0.4< pt < 0.6
        if(std::find(genIndices_matchedInPt_m04_p06.begin(), genIndices_matchedInPt_m04_p06.end(), genIndex) != genIndices_matchedInPt_m04_p06.end()){
            // Fill histograms for jets matched in pt, same separation of plots as for the matched in R
            this->fillJetHistos("matchedInPt_m04_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            if(this->isAmbiguous(genIndex, genIndices_matchedInPt_m04_p06, closestRecoJetIndices))
                this->fillJetHistos("matchedInPtAmbiguous_m04_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("matchedInPtUnambiguous_m04_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        }
                
        // -0.4< pt < 0.5
        if(std::find(genIndices_matchedInPt_m04_p05.begin(), genIndices_matchedInPt_m04_p05.end(), genIndex) != genIndices_matchedInPt_m04_p05.end()){
            // Fill histograms for jets matched in pt, same separation of plots as for the matched in R
            this->fillJetHistos("matchedInPt_m04_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            if(this->isAmbiguous(genIndex, genIndices_matchedInPt_m04_p05, closestRecoJetIndices))
                this->fillJetHistos("matchedInPtAmbiguous_m04_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("matchedInPtUnambiguous_m04_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        }
        
        // -0.3< pt < 0.7
        if(std::find(genIndices_matchedInPt_m03_p07.begin(), genIndices_matchedInPt_m03_p07.end(), genIndex) != genIndices_matchedInPt_m03_p07.end()){
            // Fill histograms for jets matched in pt, same separation of plots as for the matched in R
            this->fillJetHistos("matchedInPt_m03_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            if(this->isAmbiguous(genIndex, genIndices_matchedInPt_m03_p07, closestRecoJetIndices))
                this->fillJetHistos("matchedInPtAmbiguous_m03_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("matchedInPtUnambiguous_m03_p07_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        }
        
        // -0.3< pt < 0.6
        if(std::find(genIndices_matchedInPt_m03_p06.begin(), genIndices_matchedInPt_m03_p06.end(), genIndex) != genIndices_matchedInPt_m03_p06.end()){ 
            // Fill histograms for jets matched in pt, same separation of plots as for the matched in R
            this->fillJetHistos("matchedInPt_m03_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            if(this->isAmbiguous(genIndex, genIndices_matchedInPt_m03_p06, closestRecoJetIndices))
                this->fillJetHistos("matchedInPtAmbiguous_m03_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("matchedInPtUnambiguous_m03_p06_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        }
        
        // -0.3< pt < 0.5
        if(std::find(genIndices_matchedInPt_m03_p05.begin(), genIndices_matchedInPt_m03_p05.end(), genIndex) != genIndices_matchedInPt_m03_p05.end()){ 
            // Fill histograms for jets matched in pt, same separation of plots as for the matched in R
            this->fillJetHistos("matchedInPt_m03_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            if(this->isAmbiguous(genIndex, genIndices_matchedInPt_m03_p05, closestRecoJetIndices))
                this->fillJetHistos("matchedInPtAmbiguous_m03_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
            else
                this->fillJetHistos("matchedInPtUnambiguous_m03_p05_", recoObjects, topGenObjects, genIndex, recoIndex, genObjectIndices, weight, m_histogram);
        }
 
    }
}



void AnalyzerJetMatch::fillJetHistos(const TString& whichSelection,
                                     const RecoObjects& recoObjects, const TopGenObjects& topGenObjects,
                                     const int genIndex, const int recoIndex,
                                     const tth::GenObjectIndices& genObjectIndices,
                                     const double& weight,
                                     std::map<TString, TH1*>& m_histogram)
{
    const LV& genJet = topGenObjects.allGenJets_->at(genIndex);
    const LV& recoJet = recoObjects.jets_->at(recoIndex);
    
    // Check whether the jet is a top and/or a Higgs jet
    const bool isTopJet = genObjectIndices.uniqueGenTopMatching() &&
                          (genIndex==genObjectIndices.genBjetFromTopIndex_ || genIndex==genObjectIndices.genAntiBjetFromTopIndex_);
    const bool isHiggsJet = genObjectIndices.uniqueGenHiggsMatching() &&
                            (genIndex==genObjectIndices.genBjetFromHiggsIndex_ || genIndex==genObjectIndices.genAntiBjetFromHiggsIndex_);
    
    this->fillJetHistosInclExcl(whichSelection, "all_", genJet, recoJet, weight, m_histogram);
    if((isTopJet || isHiggsJet) && genObjectIndices.uniqueGenMatching()) this->fillJetHistosInclExcl(whichSelection, "topORhiggs_", genJet, recoJet, weight, m_histogram);
    if(isTopJet) this->fillJetHistosInclExcl(whichSelection, "top_", genJet, recoJet, weight, m_histogram);
    if(isHiggsJet) this->fillJetHistosInclExcl(whichSelection, "higgs_", genJet, recoJet, weight, m_histogram);
}



void AnalyzerJetMatch::fillJetHistosInclExcl(const TString& whichSelection, const TString& whichJets,
                                             const LV& genJet, const LV& recoJet,
                                             const double& weight,
                                             std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    const TString identifier = whichSelection + whichJets;
    
    const double deltaR = ROOT::Math::VectorUtil::DeltaR(recoJet, genJet);
    const double deltaPt = genJet.pt() - recoJet.pt();
    const double recoOverGenPt = recoJet.pt()/genJet.pt();
    
    // FIXME: why using the jet mass, what is this meant for ?
    const double genP = genJet.M();
    const double recoP = recoJet.M();          
    
    name = identifier + "genPVsRecoP";
    ((TH2D*)m_histogram.at(name))->Fill(recoP, genP, weight);
    
    name = identifier + "deltaR";
    m_histogram.at(name)->Fill(deltaR, weight);
    
    name = identifier + "deltaR_zoom";
    m_histogram.at(name)->Fill(deltaR, weight);
    
    name = identifier + "deltaRVsGenPt";
    ((TH2D*)m_histogram.at(name))->Fill(genJet.pt(), deltaR, weight);
    
    name = identifier + "recoPtOverGenPt";
    m_histogram.at(name)->Fill(recoOverGenPt, weight);
    
    name = identifier + "deltaPtOverGenPt";
    m_histogram.at(name)->Fill(deltaPt/genJet.pt(), weight);
    
    name = identifier + "genPtVsRecoPt";
    ((TH2D*)m_histogram.at(name))->Fill(recoJet.pt(), genJet.pt(), weight);
    
    name = identifier + "genEtaVsRecoEta";
    ((TH2D*)m_histogram.at(name))->Fill(recoJet.eta(), genJet.eta(), weight);
    
    name = identifier + "genPhiVsRecoPhi";
    ((TH2D*)m_histogram.at(name))->Fill(recoJet.phi(), genJet.phi(), weight);
    
    name = identifier + "genPtVsDeltaPt";
    ((TH2D*)m_histogram.at(name))->Fill(deltaPt, genJet.pt(), weight);
    
    name = identifier + "recoPtOverGenPtVsGenPt";
    ((TH2D*)m_histogram.at(name))->Fill(genJet.pt(), recoOverGenPt, weight);
    
    name = identifier + "recoPtOverGenPtVsGenEta";
    ((TH2D*)m_histogram.at(name))->Fill(genJet.eta(), recoOverGenPt, weight);
    
    name = identifier + "recoPtOverGenPtVsGenPhi";
    ((TH2D*)m_histogram.at(name))->Fill(genJet.phi(), recoOverGenPt, weight);
    
    name = identifier + "recoPtOverGenPtVsDeltaR";
    ((TH2D*)m_histogram.at(name))->Fill(deltaR, recoOverGenPt, weight);
    
    name = identifier + "deltaPtOverGenPtVsDeltaR";
    ((TH2D*)m_histogram.at(name))->Fill(deltaR, deltaPt/genJet.pt(), weight);
}



bool AnalyzerJetMatch::isAmbiguous(const int genIndex, const std::vector<int>& genIndices, const std::vector<int>& closestRecoIndices)const
{
    const int recoIndex = closestRecoIndices.at(genIndex);
    for(const int genIndexProbe : genIndices){
        if(genIndex == genIndexProbe) continue;
        if(recoIndex == closestRecoIndices.at(genIndexProbe)) return true;
    }
    
    return false;
}








