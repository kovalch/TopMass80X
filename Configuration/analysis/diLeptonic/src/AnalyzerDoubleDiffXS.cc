#include <map>
#include <vector>
#include <iostream>

#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TString.h>
#include <Math/VectorUtil.h>

#include "AnalyzerDoubleDiffXS.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"








AnalyzerDoubleDiffXS::AnalyzerDoubleDiffXS(const std::vector<TString>& selectionStepsNoCategories):
AnalyzerBaseClass("dda_", selectionStepsNoCategories)
{
    std::cout<<"--- Beginning setting up basic histograms\n";
    std::cout<<"=== Finishing setting up basic histograms\n\n";
}



void AnalyzerDoubleDiffXS::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    name = "bjet_multiplicity_vs_x1";
    m_histogram[name] = this->store(new TH2D(prefix_+name,"B-Jet Multiplicity vs x1;;N b-jets",400,0,1,21,-0.5,20.5));
    name = "bjet_multiplicity_vs_x2";
    m_histogram[name] = this->store(new TH2D(prefix_+name,"B-Jet Multiplicity vs x2;;N b-jets",400,0,1,21,-0.5,20.5));
    name = "jet_multiplicity_vs_x1";
    m_histogram[name] = this->store(new TH2D(prefix_+name,"Jet Multiplicity vs x1;;N b-jets",400,0,1,21,-0.5,20.5));
    name = "jet_multiplicity_vs_x2";
    m_histogram[name] = this->store(new TH2D(prefix_+name,"Jet Multiplicity vs x2;;N b-jets",400,0,1,21,-0.5,20.5));
    name = "dummy_vs_x1";
    m_histogram[name] = this->store(new TH2D(prefix_+name,"dummy vs x1;;",400,0,1,1,0,10));
    name = "dummy_vs_x2";
    m_histogram[name] = this->store(new TH2D(prefix_+name,"dummy vs x2;;",400,0,1,1,0,10));
    
}



void AnalyzerDoubleDiffXS::fillHistos(const RecoObjects& recoObjects, const CommonGenObjects&,
                                      const TopGenObjects&,
                                      const KinRecoObjects& kinRecoObjects,
                                      const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices&,
                                      const ttbar::GenLevelWeights&, const ttbar::RecoLevelWeights&,
                                      const double& weight, const TString&,
                                      std::map< TString, TH1* >& m_histogram)
{
   //proton Energy [GeV]  
   double protonE = 4000;
   
   double x1 = ((*kinRecoObjects.HypTop_).at(0).E()+(*kinRecoObjects.HypAntiTop_).at(0).E()+(*kinRecoObjects.HypTop_).at(0).Pz()+(*kinRecoObjects.HypAntiTop_).at(0).Pz())/(2*protonE);
   double x2 = ((*kinRecoObjects.HypTop_).at(0).E()+(*kinRecoObjects.HypAntiTop_).at(0).E()-(*kinRecoObjects.HypTop_).at(0).Pz()-(*kinRecoObjects.HypAntiTop_).at(0).Pz())/(2*protonE);
   ((TH2D*)m_histogram["bjet_multiplicity_vs_x1"])->Fill(x1,recoObjectIndices.bjetIndices_.size(),weight);
   ((TH2D*)m_histogram["bjet_multiplicity_vs_x2"])->Fill(x2,recoObjectIndices.bjetIndices_.size(),weight);
   ((TH2D*)m_histogram["jet_multiplicity_vs_x1"])->Fill(x1,recoObjectIndices.jetIndices_.size(),weight);
   ((TH2D*)m_histogram["jet_multiplicity_vs_x2"])->Fill(x2,recoObjectIndices.jetIndices_.size(),weight);
   
   if(recoObjectIndices.bjetIndices_.size()==2&&recoObjectIndices.jetIndices_.size()==2)
   {
      ((TH2D*)m_histogram["dummy_vs_x1"])->Fill(x1,1,weight);
      ((TH2D*)m_histogram["dummy_vs_x2"])->Fill(x2,1,weight);
   }
}








