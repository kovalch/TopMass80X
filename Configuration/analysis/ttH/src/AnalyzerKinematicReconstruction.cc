#include <map>
#include <vector>
#include <iostream>

#include <TH1.h>
#include <TH1D.h>
#include <TString.h>
#include <Math/VectorUtil.h>

#include "AnalyzerKinematicReconstruction.h"
#include "analysisStructs.h"
#include "JetCategories.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"








AnalyzerKinematicReconstruction::AnalyzerKinematicReconstruction(const std::vector<TString>& selectionStepsNoCategories,
                                                                 const std::vector<TString>& stepsForCategories,
                                                                 const JetCategories* jetCategories):
AnalyzerBase("kinReco_", selectionStepsNoCategories, stepsForCategories, jetCategories)
{
    std::cout<<"--- Beginning setting up kinematic reconstruction histograms\n";
    std::cout<<"=== Finishing setting up kinematic reconstruction histograms\n\n";
}



void AnalyzerKinematicReconstruction::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    name = "top_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,";pt [GeV];Events",50,0,500));
    
    name = "antiTop_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,";pt [GeV];Events",50,0,500));
    
    name = "lepton_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,";pt [GeV];Events",50,0,500));
    
    name = "antiLepton_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,";pt [GeV];Events",50,0,500));
    
    name = "neutrino_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,";pt [GeV];Events",50,0,500));
    
    name = "antiNeutrino_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,";pt [GeV];Events",50,0,500));
    
    name = "bjet_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,";pt [GeV];Events",50,0,500));
    
    name = "antiBjet_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,";pt [GeV];Events",50,0,500));
    
    name = "bjet_index";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,";index;Events",20,0,20));
    
    name = "antiBjet_index";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,";index;Events",20,0,20));
    
}



void AnalyzerKinematicReconstruction::fillHistos(const RecoObjects&, const CommonGenObjects&,
                                                 const TopGenObjects&, const HiggsGenObjects&,
                                                 const KinRecoObjects& kinRecoObjects,
                                                 const tth::RecoObjectIndices&, const tth::GenObjectIndices&,
                                                 const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                                 const double& weight, const TString&,
                                                 std::map<TString, TH1*>& m_histogram)
{
    if(!kinRecoObjects.valuesSet_) return;
    
    TString name;
    
    name = "top_pt";
    m_histogram.at(name)->Fill(kinRecoObjects.HypTop_->at(0).pt(), weight);
    
    name = "antiTop_pt";
    m_histogram.at(name)->Fill(kinRecoObjects.HypAntiTop_->at(0).pt(), weight);
    
    name = "lepton_pt";
    m_histogram.at(name)->Fill(kinRecoObjects.HypLepton_->at(0).pt(), weight);
    
    name = "antiLepton_pt";
    m_histogram.at(name)->Fill(kinRecoObjects.HypAntiLepton_->at(0).pt(), weight);
    
    name = "neutrino_pt";
    m_histogram.at(name)->Fill(kinRecoObjects.HypNeutrino_->at(0).pt(), weight);
    
    name = "antiNeutrino_pt";
    m_histogram.at(name)->Fill(kinRecoObjects.HypAntiNeutrino_->at(0).pt(), weight);
    
    name = "bjet_pt";
    m_histogram.at(name)->Fill(kinRecoObjects.HypBJet_->at(0).pt(), weight);
    
    name = "antiBjet_pt";
    m_histogram.at(name)->Fill(kinRecoObjects.HypAntiBJet_->at(0).pt(), weight);
    
    name = "bjet_index";
    m_histogram.at(name)->Fill(kinRecoObjects.HypJet0index_->at(0), weight);
    
    name = "antiBjet_index";
    m_histogram.at(name)->Fill(kinRecoObjects.HypJet1index_->at(0), weight);
}








