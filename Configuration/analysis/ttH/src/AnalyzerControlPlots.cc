#include <map>
#include <vector>
#include <iostream>
#include <algorithm>

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
#include "../../common/include/KinematicReconstructionSolution.h"





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
    
    // Primary vertices
    name = "vertex_multiplicity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Primary Vertex Multiplicity;N Vertex;Events",50,0,50));
    name = "vertex_multiplicity_unweighted";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Unweighted Primary Vertex Multiplicity;N Vertex;Events",50,0,50));
    
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
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Lepton d_{xy};d_{xy}^{l} [cm];Leptons",60,-0.25,0.25));
    name = "lepton_dxy_zoom";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Lepton d_{xy};d_{xy}^{l} [cm];Leptons",40,-0.04,0.04));
    name = "lepton_dz";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Lepton d_{z};d_{z} [cm];Leptons",60,-1,1));
    name = "lepton_dz_zoom";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Lepton d_{z};d_{z} [cm];Leptons",60,-0.06,0.06));
    
    // Leading and next-to-leading lepton
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
    name = "dilepton_deltaDxy_zoom";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Dilepton {#Delta}d_{xy};d_{xy}^{l^{+}}-d_{xy}^{l^{-}} [cm];Events",30,-0.06,0.06));
    name = "dilepton_deltaDz";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Dilepton {#Delta}d_{z};d_{z}^{l^{+}}-d_{z}^{l^{-}} [cm];Events",40,-1,1));
    name = "dilepton_deltaDz_zoom";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Dilepton {#Delta}d_{z};d_{z}^{l^{+}}-d_{z}^{l^{-}} [cm];Events",30,-0.06,0.06));
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
    name = "jet_btagDiscriminator_min";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "b-tag Discriminator d;d;Lowest d jet",60,-0.1,1.1));
    
    // Bjets
    name = "bjet_multiplicity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-Jet Multiplicity;N b-jets;Events",20,0,20));
    name = "bjet_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-Jet p_{t};p_{t}^{b-jet} [GeV];B-Jets",50,0,300));
    name = "bjet_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-Jet #eta;#eta^{b-jet};B-Jets",50,-2.6,2.6));
    name = "bjet_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-Jet #phi;#phi^{b-jet};B-Jets",50,-3.2,3.2));

    // Met
    name = "met_et";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Met E_{t};E_{t}^{met};Events",50,0,300));
    name = "met_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Met #phi;#phi^{met};Events",50,-3.2,3.2));
}



void AnalyzerControlPlots::fillHistos(const EventMetadata&,
                                      const RecoObjects& recoObjects, const CommonGenObjects&,
                                      const TopGenObjects&, const HiggsGenObjects&,
                                      const KinematicReconstructionSolutions&,
                                      const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices&,
                                      const tth::GenLevelWeights&, const tth::RecoLevelWeights& recoLevelWeights,
                                      const double& weight, const TString&,
                                      std::map<TString, TH1*>& m_histogram)
{
    const int leptonIndex = recoObjectIndices.leptonIndex_;
    const int antiLeptonIndex = recoObjectIndices.antiLeptonIndex_;
    const bool hasLeptonPair = leptonIndex>=0 && antiLeptonIndex>=0;
    const bool hasLepton = !hasLeptonPair && leptonIndex>=0;
    const bool hasAntiLepton = !hasLeptonPair && antiLeptonIndex>=0;
    
    
    // Primary vertices
    m_histogram.at("vertex_multiplicity")->Fill(recoObjects.vertMulti_, weight);    
    m_histogram.at("vertex_multiplicity_unweighted")->Fill(recoObjects.vertMulti_, recoLevelWeights.weightNoPileup_);
    
    
    // Leptons
    m_histogram.at("lepton_multiplicity")->Fill(recoObjectIndices.allLeptonIndices_.size(), weight);
    for(const int index : recoObjectIndices.allLeptonIndices_){
        const LV& lepton = recoObjects.allLeptons_->at(index);
        const double& dxy = recoObjects.lepDxyVertex0_->at(index);
        const double& dz = recoObjects.lepDzVertex0_->at(index);
        m_histogram.at("lepton_pt")->Fill(lepton.pt(), weight);
        m_histogram.at("lepton_eta")->Fill(lepton.eta(), weight);
        m_histogram.at("lepton_phi")->Fill(lepton.phi(), weight);
        m_histogram.at("lepton_dxy")->Fill(dxy, weight);
        m_histogram.at("lepton_dxy_zoom")->Fill(dxy, weight);
        m_histogram.at("lepton_dz")->Fill(dz, weight);
        m_histogram.at("lepton_dz_zoom")->Fill(dz, weight);
    }
    
    
    // Leading and next-to-leading lepton
    if(hasLepton){
        const LV& leadingLepton = recoObjects.allLeptons_->at(leptonIndex);
        m_histogram.at("lepton1st_pt")->Fill(leadingLepton.pt(), weight);
        m_histogram.at("lepton1st_eta")->Fill(leadingLepton.eta(), weight);
        m_histogram.at("lepton1st_phi")->Fill(leadingLepton.phi(), weight);
    }
    else if(hasAntiLepton){
        const LV& leadingLepton = recoObjects.allLeptons_->at(antiLeptonIndex);
        m_histogram.at("lepton1st_pt")->Fill(leadingLepton.pt(), weight);
        m_histogram.at("lepton1st_eta")->Fill(leadingLepton.eta(), weight);
        m_histogram.at("lepton1st_phi")->Fill(leadingLepton.phi(), weight);
    }
    else if(hasLeptonPair){
        const LV& leadingLepton = recoObjects.allLeptons_->at(recoObjectIndices.leadingLeptonIndex_);
        const LV& nLeadingLepton = recoObjects.allLeptons_->at(recoObjectIndices.nLeadingLeptonIndex_);
        m_histogram.at("lepton1st_pt")->Fill(leadingLepton.pt(), weight);
        m_histogram.at("lepton1st_eta")->Fill(leadingLepton.eta(), weight);
        m_histogram.at("lepton1st_phi")->Fill(leadingLepton.phi(), weight);
        m_histogram.at("lepton2nd_pt")->Fill(nLeadingLepton.pt(), weight);
        m_histogram.at("lepton2nd_eta")->Fill(nLeadingLepton.eta(), weight);
        m_histogram.at("lepton2nd_phi")->Fill(nLeadingLepton.phi(), weight);
    }
    
    
    // Dilepton
    if(hasLeptonPair){
        const LV& lepton = recoObjects.allLeptons_->at(leptonIndex);
        const LV& antiLepton = recoObjects.allLeptons_->at(antiLeptonIndex);
        const LV dilepton = lepton + antiLepton;
        const double deltaDxy = recoObjects.lepDxyVertex0_->at(leptonIndex) - recoObjects.lepDxyVertex0_->at(antiLeptonIndex);
        const double deltaDz = recoObjects.lepDzVertex0_->at(leptonIndex) - recoObjects.lepDzVertex0_->at(antiLeptonIndex);
        m_histogram.at("dilepton_mass")->Fill(dilepton.M(), weight);
        m_histogram.at("dilepton_pt")->Fill(dilepton.Pt(), weight);
        m_histogram.at("dilepton_rapidity")->Fill(dilepton.Rapidity(), weight);
        m_histogram.at("dilepton_phi")->Fill(dilepton.Phi(), weight);
        m_histogram.at("dilepton_alpha")->Fill(ROOT::Math::VectorUtil::Angle(antiLepton, lepton), weight);
        m_histogram.at("dilepton_deltaEta")->Fill(lepton.Eta() - antiLepton.Eta(), weight);
        m_histogram.at("dilepton_deltaPhi")->Fill(ROOT::Math::VectorUtil::DeltaPhi(antiLepton, lepton), weight);
        m_histogram.at("dilepton_deltaDxy")->Fill(deltaDxy, weight);
        m_histogram.at("dilepton_deltaDxy_zoom")->Fill(deltaDxy, weight);
        m_histogram.at("dilepton_deltaDz")->Fill(deltaDz, weight);
        m_histogram.at("dilepton_deltaDz_zoom")->Fill(deltaDz, weight);
    }
    
    
    // Jets
    m_histogram.at("jet_multiplicity")->Fill(recoObjectIndices.jetIndices_.size(), weight);
    double btagDiscriminator_min = 999.;
    for(const int index : recoObjectIndices.jetIndices_){
        const LV& jet = recoObjects.jets_->at(index);
        m_histogram.at("jet_pt")->Fill(jet.pt(), weight);
        m_histogram.at("jet_eta")->Fill(jet.eta(), weight);
        m_histogram.at("jet_phi")->Fill(jet.phi(), weight);
        double btagDiscriminator = recoObjects.jetBTagCSV_->at(index);
        if(btagDiscriminator < -0.1) btagDiscriminator = -0.05;
        m_histogram.at("jet_btagDiscriminator")->Fill(btagDiscriminator, weight);
        if(btagDiscriminator < btagDiscriminator_min) btagDiscriminator_min = btagDiscriminator;
    }
    m_histogram.at("jet_btagDiscriminator_min")->Fill(btagDiscriminator_min, weight);
    
    
    // Bjets
     m_histogram.at("bjet_multiplicity")->Fill(recoObjectIndices.bjetIndices_.size(), weight);
    for(const int index : recoObjectIndices.bjetIndices_){
        const LV& jet = recoObjects.jets_->at(index);
        m_histogram.at("bjet_pt")->Fill(jet.pt(), weight);
        m_histogram.at("bjet_eta")->Fill(jet.eta(), weight);
        m_histogram.at("bjet_phi")->Fill(jet.phi(), weight);
    }
    
    
    // Met
    const LV& met = *recoObjects.met_;
    m_histogram.at("met_et")->Fill(met.E(), weight);
    m_histogram.at("met_phi")->Fill(met.phi(), weight);
}








