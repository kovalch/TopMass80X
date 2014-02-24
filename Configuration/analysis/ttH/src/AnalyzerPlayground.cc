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

#include "AnalyzerPlayground.h"
#include "analysisStructs.h"
#include "JetCategories.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/classes.h"





AnalyzerPlayground::AnalyzerPlayground(const std::vector<TString>& selectionStepsNoCategories,
                                       const std::vector<TString>& stepsForCategories,
                                       const JetCategories* jetCategories):
AnalyzerBaseClass("test_", selectionStepsNoCategories, stepsForCategories, jetCategories)
{
    std::cout<<"--- Beginning setting up playground\n";
    std::cout<<"=== Finishing setting up playground\n\n";
}



void AnalyzerPlayground::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    
    // Book histograms here
    name = "blah1";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"histo title;x axis;y axis",10,0,100));
    
    name = "blah2_blubb";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"histo title;x axis;y axis",10,0,100));
}



void AnalyzerPlayground::fillHistos(const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                    const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                                    const KinRecoObjects& kinRecoObjects,
                                    const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                                    const tth::GenLevelWeights& genLevelWeights, const tth::RecoLevelWeights& recoLevelWeights,
                                    const double& weight, const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    
    // Do calculations and filling of histograms
    const bool nonsenseBool = recoObjects.valuesSet_ && commonGenObjects.valuesSet_ && topGenObjects.valuesSet_ &&
                              higgsGenObjects.valuesSet_ && kinRecoObjects.valuesSet_;
    
    const int nonsenseInt = recoObjectIndices.antiLeptonIndex_ + genObjectIndices.genAntiBjetFromHiggsIndex_;
    
    const double nonsenseDouble = genLevelWeights.trueLevelWeight_ * recoLevelWeights.weight_ * 10.;
    
    
    
    name = "blah1";
    if(nonsenseBool) m_histogram[name]->Fill(1., weight);
    else if(step == "nonsense") m_histogram[name]->Fill(2., weight);
    else m_histogram[name]->Fill(3., weight);
    
    name = "blah2_blubb";
    m_histogram[name]->Fill(nonsenseInt, weight);
    m_histogram[name]->Fill(nonsenseDouble, weight);
}






