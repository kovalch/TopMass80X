#include <map>
#include <vector>
#include <iostream>

#include <TH1.h>
#include <TH1D.h>
#include <TString.h>
#include <Math/VectorUtil.h>

#include "AnalyzerControlPlots.h"
#include "analysisStructs.h"
#include "JetCategories.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"








AnalyzerControlPlots::AnalyzerControlPlots(const std::vector<TString>& selectionStepsNoCategories,
                                           const std::vector<TString>& stepsForCategories,
                                           const JetCategories* jetCategories):
AnalyzerBase("basic_", selectionStepsNoCategories, stepsForCategories, jetCategories)
{
    std::cout<<"--- Beginning setting up basic histograms\n";
    std::cout<<"=== Finishing setting up basic histograms\n\n";
}



void AnalyzerControlPlots::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    // Vertices
    name = "vertex_multiplicity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Primary Vertex Multiplicity;N Vertex;Events",50,0,50));
    
    // Leptons
    name = "lepton_multiplicity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Lepton multiplicity;N leptons;Events",10,0,10));
    name = "lepton_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Lepton p_{t};p_{t}^{l} [GeV];Leptons",50,0,250));
    name = "lepton_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Lepton #eta;#eta^{l};Leptons",50,-2.6,2.6));
    name = "lepton_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Lepton #phi;#phi^{l};Leptons",50,-3.2,3.2));
    name = "lepton_dxy";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Lepton d_{xy};d_{xy}^{l} [cm];Events",60,-0.25,0.25));
    //name = "lepton_dz";
    //m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Lepton d_{z};d_{z} [cm];Events",60,-1,1));
    
    // Leading lepton and antilepton
    name = "lepton1st_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st Lepton p_{t};p_{t}^{l_{1}} [GeV];Leptons",50,0,250));
    name = "lepton1st_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st Lepton #eta;#eta^{l_{1}};Leptons",50,-2.6,2.6));
    name = "lepton1st_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st Lepton #phi;#phi^{l_{1}};Leptons",50,-3.2,3.2));
    name = "lepton2nd_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd Lepton p_{t};p_{t}^{l_{2}} [GeV];Leptons",50,0,250));
    name = "lepton2nd_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd Lepton #eta;#eta^{l_{2}};Leptons",50,-2.6,2.6));
    name = "lepton2nd_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd Lepton #phi;#phi^{l_{2}};Leptons",50,-3.2,3.2));

    // Dilepton
    name = "dilepton_mass";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Dilepton mass;m^{l^{+}l^{-}} [GeV];Events",50,0,350));
    name = "dilepton_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Dilepton p_{t};p_{t}^{l^{+}l^{-}} [GeV];Events",50,0,300));
    name = "dilepton_rapidity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Dilepton rapidity;y^{l^{+}l^{-}};Events",50,-2.6,2.6));
    name = "dilepton_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Dilepton #phi;#phi^{l^{+}l^{-}};Events",50,-3.2,3.2));
    name = "dilepton_deltaEta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Dilepton #Delta#eta;#eta^{l^{+}}-#eta^{l^{-}};Events",50,-5,5));
    name = "dilepton_deltaPhi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Dilepton #Delta#phi;#phi^{l^{+}}-#phi^{l^{-}};Events",50,-3.2,3.2));
    name = "dilepton_deltaDxy";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Dilepton {#Delta}d_{xy};d_{xy}^{l^{+}}-d_{xy}^{l^{-}} [cm];Events",40,-0.3,0.3));
    //name = "dilepton_deltaDz";
    //m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Dilepton {#Delta}d_{z};d_{z}^{l^{+}}-d_{z}^{l^{-}} [cm];Events",40,-1,1));
    name = "dilepton_alpha";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Dilepton opening angle #alpha^{l^{+}l^{-}};#alpha^{l^{+}l^{-}};Events",40,0,3.2));
    
    // Jets
    name = "jet_multiplicity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Jet Multiplicity;N jets;Events",20,0,20));
    name = "jet_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Jet p_{t};p_{t}^{jet} [GeV];Jets",50,0,300));
    name = "jet_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Jet #eta;#eta^{jet};Jets",50,-2.6,2.6));
    name = "jet_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Jet #phi;#phi^{jet};Jets",50,-3.2,3.2));
    name = "jet_btagDiscriminator";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "b-tag Discriminator d;d;Jets",60,-0.1,1.1));
    name = "jet_chargeGlobalPtWeighted";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "jetChargeGlobalPtWeighted c_{glob}^{jet}; c_{glob}^{jet};# jets", 110, -1.1, 1.1));
    name = "jet_chargeRelativePtWeighted";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "jetChargeRelativePtWeighted c_{rel}^{jet}; c_{rel}^{jet};# jets", 110, -1.1, 1.1));
    
    // Bjets
    name = "bjet_multiplicity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-Jet Multiplicity;N b-jets;Events",20,0,20));
    name = "bjet_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-Jet p_{t};p_{t}^{b-jet} [GeV];B-Jets",50,0,300));
    name = "bjet_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-Jet #eta;#eta^{b-jet};B-Jets",50,-2.6,2.6));
    name = "bjet_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-Jet #phi;#phi^{b-jet};B-Jets",50,-3.2,3.2));
    name = "bjet_chargeRelativePtWeighted";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-JetChargeRelativePtWeighted c_{rel}^{jet}; c_{rel}^{jet};# B-Jets", 110, -1.1, 1.1));

    // Met
    name = "met_et";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Met E_{t};E_{t}^{met};Events",50,0,300));
    name = "met_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Met #phi;#phi^{met};Events",50,-3.2,3.2));
}



void AnalyzerControlPlots::fillHistos(const RecoObjects& recoObjects, const CommonGenObjects&,
                                      const TopGenObjects&, const HiggsGenObjects&,
                                      const KinRecoObjects&,
                                      const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices&,
                                      const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                      const double& weight, const TString&,
                                      std::map< TString, TH1* >& m_histogram)
{
    // Vertices
    m_histogram["vertex_multiplicity"]->Fill(recoObjects.vertMulti_, weight);    
    
    
    // Leptons
    m_histogram["lepton_multiplicity"]->Fill(recoObjectIndices.allLeptonIndices_.size(), weight);
    for(const int index : recoObjectIndices.leptonIndices_){
        m_histogram["lepton_pt"]->Fill(recoObjects.allLeptons_->at(index).Pt(), weight);
        m_histogram["lepton_eta"]->Fill(recoObjects.allLeptons_->at(index).Eta(), weight);
        m_histogram["lepton_phi"]->Fill(recoObjects.allLeptons_->at(index).Phi(), weight);
        m_histogram["lepton_dxy"]->Fill(recoObjects.lepDxyVertex0_->at(index), weight);
    }
    for(const int index : recoObjectIndices.antiLeptonIndices_){
        m_histogram["lepton_pt"]->Fill(recoObjects.allLeptons_->at(index).Pt(), weight);
        m_histogram["lepton_eta"]->Fill(recoObjects.allLeptons_->at(index).Eta(), weight);
        m_histogram["lepton_phi"]->Fill(recoObjects.allLeptons_->at(index).Phi(), weight);
        m_histogram["lepton_dxy"]->Fill(recoObjects.lepDxyVertex0_->at(index), weight);
    }
    
    const int leptonIndex = recoObjectIndices.leptonIndices_.size()>0 ? recoObjectIndices.leptonIndices_.at(0) : -1;
    const int antiLeptonIndex = recoObjectIndices.antiLeptonIndices_.size()>0 ? recoObjectIndices.antiLeptonIndices_.at(0) : -1;
    const bool hasLeptonPair = (leptonIndex!=-1 && antiLeptonIndex!=-1);
    
    // Leading lepton and antilepton
    int leadingLeptonIndex(leptonIndex);
    int nLeadingLeptonIndex(antiLeptonIndex);
    if(hasLeptonPair){
        common::orderIndices(leadingLeptonIndex, nLeadingLeptonIndex, *recoObjects.allLeptons_, common::LVpt);
        
        m_histogram["lepton1st_pt"]->Fill(recoObjects.allLeptons_->at(leadingLeptonIndex).Pt(), weight);
        m_histogram["lepton1st_eta"]->Fill(recoObjects.allLeptons_->at(leadingLeptonIndex).Eta(), weight);
        m_histogram["lepton1st_phi"]->Fill(recoObjects.allLeptons_->at(leadingLeptonIndex).Phi(), weight);
        
        m_histogram["lepton2nd_pt"]->Fill(recoObjects.allLeptons_->at(nLeadingLeptonIndex).Pt(), weight);
        m_histogram["lepton2nd_eta"]->Fill(recoObjects.allLeptons_->at(nLeadingLeptonIndex).Eta(), weight);
        m_histogram["lepton2nd_phi"]->Fill(recoObjects.allLeptons_->at(nLeadingLeptonIndex).Phi(), weight);
    }
    
    
    // Dilepton
    if(hasLeptonPair){
        LV dilepton(0.,0.,0.,0.);
        dilepton = recoObjects.allLeptons_->at(leadingLeptonIndex) + recoObjects.allLeptons_->at(nLeadingLeptonIndex);
        
        m_histogram["dilepton_mass"]->Fill(dilepton.M(), weight);
        m_histogram["dilepton_pt"]->Fill(dilepton.Pt(), weight);
        m_histogram["dilepton_rapidity"]->Fill(dilepton.Rapidity(), weight);
        m_histogram["dilepton_phi"]->Fill(dilepton.Phi(), weight);

        m_histogram["dilepton_deltaEta"]->Fill(recoObjects.allLeptons_->at(leptonIndex).Eta() - recoObjects.allLeptons_->at(antiLeptonIndex).Eta(), weight);
        m_histogram["dilepton_deltaPhi"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(recoObjects.allLeptons_->at(antiLeptonIndex), recoObjects.allLeptons_->at(leptonIndex)), weight);
        m_histogram["dilepton_deltaDxy"]->Fill(recoObjects.lepDxyVertex0_->at(leptonIndex) - recoObjects.lepDxyVertex0_->at(antiLeptonIndex), weight);
        
        m_histogram["dilepton_alpha"]->Fill(ROOT::Math::VectorUtil::Angle(recoObjects.allLeptons_->at(antiLeptonIndex), recoObjects.allLeptons_->at(leptonIndex)), weight);
    }


    // Jets
    m_histogram["jet_multiplicity"]->Fill(recoObjectIndices.jetIndices_.size(), weight);
    for(const int index : recoObjectIndices.jetIndices_){
        m_histogram["jet_pt"]->Fill(recoObjects.jets_->at(index).Pt(), weight);
        m_histogram["jet_eta"]->Fill(recoObjects.jets_->at(index).Eta(), weight);
        m_histogram["jet_phi"]->Fill(recoObjects.jets_->at(index).Phi(), weight);
        double btagDiscriminator = recoObjects.jetBTagCSV_->at(index);
        if(btagDiscriminator < -0.1) btagDiscriminator = -0.05;
        m_histogram["jet_btagDiscriminator"]->Fill(btagDiscriminator, weight);
        m_histogram["jet_chargeGlobalPtWeighted"]->Fill(recoObjects.jetChargeGlobalPtWeighted_->at(index), weight);
        m_histogram["jet_chargeRelativePtWeighted"]->Fill(recoObjects.jetChargeRelativePtWeighted_->at(index), weight);
    }
    
    
    // Bjets
     m_histogram["bjet_multiplicity"]->Fill(recoObjectIndices.bjetIndices_.size(), weight);
    for(const int index : recoObjectIndices.bjetIndices_){
        m_histogram["bjet_pt"]->Fill(recoObjects.jets_->at(index).Pt(), weight);
        m_histogram["bjet_eta"]->Fill(recoObjects.jets_->at(index).Eta(), weight);
        m_histogram["bjet_phi"]->Fill(recoObjects.jets_->at(index).Phi(), weight);
        m_histogram["bjet_chargeRelativePtWeighted"]->Fill(recoObjects.jetChargeRelativePtWeighted_->at(index), weight);
    }
    
    
    // Met
    m_histogram["met_et"]->Fill(recoObjects.met_->E(), weight);
    m_histogram["met_phi"]->Fill(recoObjects.met_->Phi(), weight);
}








