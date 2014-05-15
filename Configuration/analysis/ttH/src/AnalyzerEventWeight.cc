#include <iostream>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <algorithm>

#include <TString.h>
#include <TH1.h>
#include <TH1D.h>
#include <Math/VectorUtil.h>
#include <TProfile.h>

#include "AnalyzerEventWeight.h"
#include "analysisStructs.h"
#include "JetCategories.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/classes.h"





AnalyzerEventWeight::AnalyzerEventWeight(const std::vector<TString>& selectionStepsNoCategories,
                                         const std::vector<TString>& stepsForCategories,
                                         const JetCategories* jetCategories):
AnalyzerBase("weight_", selectionStepsNoCategories, stepsForCategories, jetCategories)
{
    std::cout<<"--- Beginning setting up event weight analyzer\n";
    std::cout<<"=== Finishing setting up event weight analyzer\n\n";
}



void AnalyzerEventWeight::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    // Generator level weights
    name = "madgraphCorrection";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"MadGraph W decay BR;weight;# events",100,0.5,2.));
    name = "pileup";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Pileup;weight;# events",100,0,10));
    name = "generator";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Generator;weight;# events",100,0,10));
    name = "topPt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Top pt;weight;# events",100,0.5,2.));
    name = "trueLevel_noPileup";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"True level (no pileup);weight;# events",100,0,10));
    name = "trueLevel";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"True level;weight;# events",100,0,10));
    
    // Reconstruction level weights
    name = "leptonSF";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Lepton scale factor;weight;# events",100,0.5,2.));
    name = "triggerSF";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Trigger scale factor;weight;# events",100,0.5,2.));
    name = "btagSF";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Btag scale factor;weight;# events",100,0,10));
    name = "kinReco";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Kin reco scale factor;weight;# events",100,0.5,2.));
    name = "recoLevel_noPileup";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Reco level (no pileup);weight;# events",100,0,10));
    name = "recoLevel";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Reco level;weight;# events",100,0,10));
}



void AnalyzerEventWeight::fillHistos(const RecoObjects&, const CommonGenObjects&,
                                     const TopGenObjects&, const HiggsGenObjects&,
                                     const KinRecoObjects&,
                                     const tth::RecoObjectIndices&, const tth::GenObjectIndices&,
                                     const tth::GenLevelWeights& genLevelWeights, const tth::RecoLevelWeights& recoLevelWeights,
                                     const double&, const TString&, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    // Generator level weights
    name = "madgraphCorrection";
    m_histogram.at(name)->Fill(genLevelWeights.weightMadgraphCorrection_);
    name = "pileup";
    m_histogram.at(name)->Fill(genLevelWeights.weightPileup_);
    name = "generator";
    m_histogram.at(name)->Fill(genLevelWeights.weightGenerator_);
    name = "topPt";
    m_histogram.at(name)->Fill(genLevelWeights.weightTopPt_);
    name = "trueLevel_noPileup";
    m_histogram.at(name)->Fill(genLevelWeights.trueLevelWeightNoPileup_);
    name = "trueLevel";
    m_histogram.at(name)->Fill(genLevelWeights.trueLevelWeight_);
    
    // Reconstruction level weights
    name = "leptonSF";
    m_histogram.at(name)->Fill(recoLevelWeights.weightLeptonSF_);
    name = "triggerSF";
    m_histogram.at(name)->Fill(recoLevelWeights.weightTriggerSF_);
    name = "btagSF";
    m_histogram.at(name)->Fill(recoLevelWeights.weightBtagSF_);
    name = "kinReco";
    m_histogram.at(name)->Fill(recoLevelWeights.weightKinReco_);
    name = "recoLevel_noPileup";
    m_histogram.at(name)->Fill(recoLevelWeights.weightNoPileup_);
    name = "recoLevel";
    m_histogram.at(name)->Fill(recoLevelWeights.weight_);
}






