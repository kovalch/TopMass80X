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
#include "../../common/include/KinematicReconstructionSolution.h"





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
    
    
    name = "efficiency";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,";;Events",6,0,6));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "events");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "solution");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "tt in event");
    m_histogram[name]->GetXaxis()->SetBinLabel(4, "found");
    m_histogram[name]->GetXaxis()->SetBinLabel(5, "correct");
    m_histogram[name]->GetXaxis()->SetBinLabel(6, "swapped");
    
    name = "top_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,";pt [GeV];Events",50,0,500));
    
    name = "antiTop_pt";
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



void AnalyzerKinematicReconstruction::fillHistos(const EventMetadata&,
                                                 const RecoObjects&, const CommonGenObjects&,
                                                 const TopGenObjects&, const HiggsGenObjects&,
                                                 const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                                                 const tth::RecoObjectIndices&, const tth::GenObjectIndices& genObjectIndices,
                                                 const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                                 const double& weight, const TString&,
                                                 std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    name = "efficiency";
    m_histogram.at(name)->Fill(0., weight);
    
    if(!kinematicReconstructionSolutions.numberOfSolutions()) return;
    const KinematicReconstructionSolution& solution = kinematicReconstructionSolutions.solution();
    
    m_histogram.at(name)->Fill(1., weight);
    if(genObjectIndices.uniqueRecoTopMatching()){
        m_histogram.at(name)->Fill(2., weight);
        const int bjetIndex = genObjectIndices.recoBjetFromTopIndex_;
        const int antiBjetIndex = genObjectIndices.recoAntiBjetFromTopIndex_;
        if(bjetIndex==solution.bjetIndex() && antiBjetIndex==solution.antiBjetIndex()){
            m_histogram.at(name)->Fill(3., weight);
            m_histogram.at(name)->Fill(4., weight);
        }
        if(antiBjetIndex==solution.bjetIndex() && bjetIndex==solution.antiBjetIndex()){
            m_histogram.at(name)->Fill(3., weight);
            m_histogram.at(name)->Fill(5., weight);
        }
    }
    
    name = "top_pt";
    m_histogram.at(name)->Fill(solution.top().pt(), weight);
    
    name = "antiTop_pt";
    m_histogram.at(name)->Fill(solution.antiTop().pt(), weight);
    
    name = "neutrino_pt";
    m_histogram.at(name)->Fill(solution.neutrino().pt(), weight);
    
    name = "antiNeutrino_pt";
    m_histogram.at(name)->Fill(solution.antiNeutrino().pt(), weight);
    
    name = "bjet_pt";
    m_histogram.at(name)->Fill(solution.bjet().pt(), weight);
    
    name = "antiBjet_pt";
    m_histogram.at(name)->Fill(solution.antiBjet().pt(), weight);
    
    name = "bjet_index";
    m_histogram.at(name)->Fill(solution.bjetIndex(), weight);
    
    name = "antiBjet_index";
    m_histogram.at(name)->Fill(solution.antiBjetIndex(), weight);
}








