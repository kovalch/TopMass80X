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
    
     name = "top_rapidity_vs_top_pt";
     m_histogram[name] = this->store(new TH2D (prefix_+name+step,"Top Rapidity vs Top pT; p_{T}^{t} [GeV];y(t)",1200,0,1200,200,-5,5));
    
     name = "gen_top_rapidity_vs_top_pt";
     m_histogram[name] = this->store(new TH2D (prefix_+name+step,"Top Rapidity vs Top pT; p_{T}^{t} [GeV];y(t)",1200,0,1200,200,-5,5));
    
    name = "bjet_multiplicity_vs_x1";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"B-Jet Multiplicity vs x1;;N b-jets",400,0,1,21,-0.5,20.5));
    name = "bjet_multiplicity_vs_x2";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"B-Jet Multiplicity vs x2;;N b-jets",400,0,1,21,-0.5,20.5));
    name = "jet_multiplicity_vs_x1";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs x1;;N b-jets",400,0,1,21,-0.5,20.5));
    name = "jet_multiplicity_vs_x2";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs x2;;N b-jets",400,0,1,21,-0.5,20.5));
    name = "dummy_vs_x1";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"dummy vs x1;;",400,0,1,1,0,10));
    name = "dummy_vs_x2";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"dummy vs x2;;",400,0,1,1,0,10));
    
    name = "jet_multiplicity_vs_top_pt";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs top pt;pt(top), [GeV];N b-jets",1200,0,1200,21,-0.5,20.5));
    name = "jet_multiplicity_vs_top_y";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs top y;y(top) ;N b-jets",200,-5,5,21,-0.5,20.5));
    name = "jet_multiplicity_vs_ttbar_pt";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs ttbar pt;pt(ttbar), [GeV];N b-jets",400,0,400,21,-0.5,20.5));
    name = "jet_multiplicity_vs_ttbar_y";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs ttbar y;y(ttbar) ;N b-jets",200,-5,5,21,-0.5,20.5));
    name = "jet_multiplicity_vs_ttbar_mass";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"Jet Multiplicity vs ttbar mass;mass(ttbar), [GeV];N b-jets",2000,0,2000,21,-0.5,20.5));
    
    name = "ttbar_mass_vs_top_pt";
    m_histogram[name] = this->store(new TH2D(prefix_+name+step,"ttbar mass vs top pt;pt(top), [GeV];N b-jets",1200,0,1200,2000,0,2000));
    
//     name = "ttbar_mass_vs_x";
//     m_histogram[name] = this->store(new TH2D(prefix_+name,"ttbar mass vs x;x;N b-jets",400,0,1,2000,0,2000));
//     name = "ttbar_pt_vs_x";
//     m_histogram[name] = this->store(new TH2D(prefix_+name," vs x;x;N b-jets",400,0,1,400,0,400));
//     name = "ttbar_rapidity_vs_x";
//     m_histogram[name] = this->store(new TH2D(prefix_+name," vs x;x;N b-jets",400,0,1,200,-5,5));
//     name = "top_pt_vs_x";
//     m_histogram[name] = this->store(new TH2D(prefix_+name," vs x;x;N b-jets",400,0,1,1200,0,1200));
//     name = "top_rapidity_vs_x";
//     m_histogram[name] = this->store(new TH2D(prefix_+name," vs x;x;N b-jets",400,0,1,200,-5,5));
    
}



void AnalyzerDoubleDiffXS::fillHistos(const RecoObjects& , const CommonGenObjects&,
                                      const TopGenObjects& topGenObjects,
                                      const KinRecoObjects& kinRecoObjects,
                                      const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices&,
                                      const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights&,
                                      const double& weight, const TString&,
                                      std::map< TString, TH1* >& m_histogram)
{
    
   if(kinRecoObjects.valuesSet_){
    
    ((TH2D*)m_histogram["top_rapidity_vs_top_pt"])->Fill((*kinRecoObjects.HypTop_).at(0).Pt(),(*kinRecoObjects.HypTop_).at(0).Rapidity(),weight);
    ((TH2D*)m_histogram["top_rapidity_vs_top_pt"])->Fill((*kinRecoObjects.HypAntiTop_).at(0).Pt(),(*kinRecoObjects.HypAntiTop_).at(0).Rapidity(),weight);

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
   
   
   ((TH2D*)m_histogram["jet_multiplicity_vs_top_pt"])->Fill((*kinRecoObjects.HypTop_).at(0).Pt(),recoObjectIndices.jetIndices_.size(),weight);
   ((TH2D*)m_histogram["jet_multiplicity_vs_top_pt"])->Fill((*kinRecoObjects.HypAntiTop_).at(0).Pt(),recoObjectIndices.jetIndices_.size(),weight);
   
   ((TH2D*)m_histogram["jet_multiplicity_vs_top_y"])->Fill((*kinRecoObjects.HypTop_).at(0).Rapidity(),recoObjectIndices.jetIndices_.size(),weight);
   ((TH2D*)m_histogram["jet_multiplicity_vs_top_y"])->Fill((*kinRecoObjects.HypAntiTop_).at(0).Rapidity(),recoObjectIndices.jetIndices_.size(),weight);
   
   ((TH2D*)m_histogram["jet_multiplicity_vs_ttbar_pt"])->Fill(((*kinRecoObjects.HypTop_).at(0)+(*kinRecoObjects.HypAntiTop_).at(0)).Pt(),recoObjectIndices.jetIndices_.size(),weight);
   ((TH2D*)m_histogram["jet_multiplicity_vs_ttbar_y"])->Fill(((*kinRecoObjects.HypTop_).at(0)+(*kinRecoObjects.HypAntiTop_).at(0)).Rapidity(),recoObjectIndices.jetIndices_.size(),weight);
   ((TH2D*)m_histogram["jet_multiplicity_vs_ttbar_mass"])->Fill(((*kinRecoObjects.HypTop_).at(0)+(*kinRecoObjects.HypAntiTop_).at(0)).M(),recoObjectIndices.jetIndices_.size(),weight);
    
   ((TH2D*)m_histogram["ttbar_mass_vs_top_pt"])->Fill((*kinRecoObjects.HypTop_).at(0).Pt(),((*kinRecoObjects.HypTop_).at(0)+(*kinRecoObjects.HypAntiTop_).at(0)).M(),weight);
   ((TH2D*)m_histogram["ttbar_mass_vs_top_pt"])->Fill((*kinRecoObjects.HypAntiTop_).at(0).Pt(),((*kinRecoObjects.HypTop_).at(0)+(*kinRecoObjects.HypAntiTop_).at(0)).M(),weight);
   
   }
   
   if(topGenObjects.valuesSet_){
     double trueLevelWeight = genLevelWeights.trueLevelWeight_;
    ((TH2D*)m_histogram["gen_top_rapidity_vs_top_pt"])->Fill((*topGenObjects.GenTop_).Pt(),(*topGenObjects.GenTop_).Rapidity(),trueLevelWeight);
    ((TH2D*)m_histogram["gen_top_rapidity_vs_top_pt"])->Fill((*topGenObjects.GenAntiTop_).Pt(),(*topGenObjects.GenAntiTop_).Rapidity(),trueLevelWeight);
    
   }
}








