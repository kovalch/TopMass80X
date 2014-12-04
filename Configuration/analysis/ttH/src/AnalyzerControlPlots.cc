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
    
    // Vertices
    name = "vertex_multiplicity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Primary Vertex Multiplicity;N Vertex;Events",50,0,50));
    
    // Secondary Vertex
    name = "secondaryVertex_multiplicityPerEvent";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Secondary Vertex Multiplicity;SecondaryVertex multiplicity;Events",10,0,10));    
    name = "secondaryVertex_multiplicityPerJet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Secondary Vertex Multiplicity;SecondaryVertex multiplicity;Jets",8,0,8));    
    name = "secondaryVertex_trackMultiplicityPerEvent";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Secondary Vertex Track Multiplicity;SecondaryVertex track multiplicity;Events",30,0,30));    
    name = "secondaryVertex_trackMultiplicityPerJet";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Secondary Vertex Track Multiplicity;SecondaryVertex track multiplicity;Jets",20,0,20));
    name = "secondaryVertex_trackMultiplicityPerSecondaryVertex";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Secondary Vertex Track Multiplicity;SecondaryVertex track multiplicity;Secondary Vertices",20,0,20));
    name = "secondaryVertex_mass";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Secondary Vertex Mass;M(Secondary Vertex) [GeV];Secondary Vertices",80,0,8));
    name = "secondaryVertex_correctedMass";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Secondary Vertex Mass (corrected);MC corrected M(Secondary Vertex) [GeV];Secondary Vertices",80,0,8));
    name = "allSecondaryVertexTracks_mass";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Secondary Vertex Mass (p_{t}-corrected);p_{t}-corrected M(Secondary Vertex) [GeV];Jets",80,0,8));
    name = "allSecondaryVertexTracks_correctedMass";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Secondary Vertex Mass (MC p_{t}-corrected);MC p_{t}-corrected M(Secondary Vertex) [GeV];Jets",80,0,8));
    name = "sumSecondaryVertex_mass";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Sum of the Secondary Vertex Mass;#sum M(Secondary Vertex) [GeV];Jets",80,0,8));
    name = "sumSecondaryVertex_correctedMass";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Sum of the Secondary Vertex Mass (MC corrected);#sum MC corrected M(Secondary Vertex) [GeV];Jets",80,0,8));
    name = "secondaryVertex_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Secondary Vertex p_{t};p_{t}^{SecondaryVertex} [GeV];Secondary Vertices",50,0,250));
    name = "secondaryVertex_ptOverJetPt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Secondary Vertex p_{t} Over Jet p_{t};p_{t}^{SecondaryVertex} / p_{t}^{jet};Events",30,0,12));
    name = "secondaryVertex_ptOverJetPt_zoom";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Secondary Vertex p_{t} Over Jet p_{t};p_{t}^{SecondaryVertex} / p_{t}^{jet};Events",50,0,1));
    name = "secondaryVertex_flightDistance";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Secondary Vertex Flight Distance; Secondary Vertex Flight Distance [cm];Secondary Vertices",50,0,10));
    name = "secondaryVertex_flightDistanceSignificance";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Secondary Vertex Flight Distance Significance; Secondary Vertex Flight Distance Significance; Secondary Vertices",50,0,250));
    
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
    m_histogram[name]->Sumw2();
    name = "jet_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Jet p_{t};p_{t}^{jet} [GeV];Jets",50,0,300));
    m_histogram[name]->Sumw2();
    name = "jet_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Jet #eta;#eta^{jet};Jets",50,-2.6,2.6));
    m_histogram[name]->Sumw2();
    name = "jet_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Jet #phi;#phi^{jet};Jets",50,-3.2,3.2));
    m_histogram[name]->Sumw2();
    name = "jet_btagDiscriminator";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "b-tag Discriminator d;d;Jets",60,-0.1,1.1));
    m_histogram[name]->Sumw2();
    name = "jet_btagDiscriminator_min";
    m_histogram[name] = store(new TH1D(prefix_+name+step, "b-tag Discriminator d;d;Lowest d jet",60,-0.1,1.1));
    m_histogram[name]->Sumw2();
    name = "jet_chargeGlobalPtWeighted";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "jetChargeGlobalPtWeighted c_{glob}^{jet}; c_{glob}^{jet};# jets", 110, -1.1, 1.1));
    m_histogram[name]->Sumw2();
    name = "jet_chargeRelativePtWeighted";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "jetChargeRelativePtWeighted c_{rel}^{jet}; c_{rel}^{jet};# jets", 110, -1.1, 1.1));
    m_histogram[name]->Sumw2();
    
    // Bjets
    name = "bjet_multiplicity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-Jet Multiplicity;N b-jets;Events",20,0,20));
    m_histogram[name]->Sumw2();
    name = "bjet_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-Jet p_{t};p_{t}^{b-jet} [GeV];B-Jets",50,0,300));
    m_histogram[name]->Sumw2();
    name = "bjet_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-Jet #eta;#eta^{b-jet};B-Jets",50,-2.6,2.6));
    m_histogram[name]->Sumw2();
    name = "bjet_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-Jet #phi;#phi^{b-jet};B-Jets",50,-3.2,3.2));
    m_histogram[name]->Sumw2();
    name = "bjet_chargeRelativePtWeighted";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-JetChargeRelativePtWeighted c_{rel}^{jet}; c_{rel}^{jet};# B-Jets", 110, -1.1, 1.1));
    m_histogram[name]->Sumw2();

    // Met
    name = "met_et";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Met E_{t};E_{t}^{met};Events",50,0,300));
    name = "met_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Met #phi;#phi^{met};Events",50,-3.2,3.2));
}



void AnalyzerControlPlots::fillHistos(const EventMetadata&,
                                      const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                      const TopGenObjects&, const HiggsGenObjects&,
                                      const KinematicReconstructionSolutions&,
                                      const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices&,
                                      const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                      const double& weight, const TString&,
                                      std::map<TString, TH1*>& m_histogram)
{
    // Vertices
    m_histogram.at("vertex_multiplicity")->Fill(recoObjects.vertMulti_, weight);    
    
    // Secondary Vertex
    const std::vector<int>& jetSecondaryVertexTrackVertexIndex = (recoObjects.valuesSet_) ? *recoObjects.jetSecondaryVertexTrackVertexIndex_ : std::vector<int>(0);
    const std::vector<double>&  jetSecondaryVertexFlightDistanceValue = (recoObjects.valuesSet_) ? *recoObjects.jetSecondaryVertexFlightDistanceValue_ : std::vector<double>(0);
    const std::vector<double>&  jetSecondaryVertexFlightDistanceSignificance = (recoObjects.valuesSet_) ? *recoObjects.jetSecondaryVertexFlightDistanceSignificance_ : std::vector<double>(0);
    const std::vector<double>&  jetSecondaryVertexPtCorrectedMass = (recoObjects.valuesSet_) ? *recoObjects.jetSecondaryVertexPtCorrectedMass_ : std::vector<double>(0);
    const std::vector<int>&  jetSecondaryVertexJetIndex = (recoObjects.valuesSet_) ? *recoObjects.jetSecondaryVertexJetIndex_ : std::vector<int>(0);
    const std::vector<LV>&  jetSecondaryVertex = (recoObjects.valuesSet_) ? *recoObjects.jetSecondaryVertex_ : std::vector<LV>(0);
    const VLV& recoJets = *recoObjects.jets_;
    
    // Find if the currently analyzed file is data or MC 
    bool isMC = true;
    if(!commonGenObjects.valuesSet_) isMC = false;
    
    unsigned int nSecondaryVertex=jetSecondaryVertex.size();
    m_histogram.at("secondaryVertex_multiplicityPerEvent")->Fill(nSecondaryVertex, weight);    
    
    for(size_t iSecondaryVertex=0; iSecondaryVertex<jetSecondaryVertex.size(); ++iSecondaryVertex) {
        m_histogram.at("secondaryVertex_mass")->Fill(jetSecondaryVertex.at(iSecondaryVertex).M(), weight);
        
        // Correct for effects on MC by multiplying the secondary vertex mass of MC with 0.98
        if (isMC)  m_histogram.at("secondaryVertex_correctedMass")->Fill(0.98*jetSecondaryVertex.at(iSecondaryVertex).M(), weight);
        else m_histogram.at("secondaryVertex_correctedMass")->Fill(jetSecondaryVertex.at(iSecondaryVertex).M(), weight);
        
        m_histogram.at("secondaryVertex_pt")->Fill(jetSecondaryVertex.at(iSecondaryVertex).Pt(), weight);        
    }
    
    
    // Loop over all the jets
    for(size_t iJet=0; iJet<jetSecondaryVertexPtCorrectedMass.size(); ++iJet) {
        unsigned int secondaryVertexMultiplicityPerJet = 0;
        
        // Sum of the secondary vertex masses for the secondary vertices that belong to the same jet
        double SumSecondaryVertexMass = 0.0;
        double SumSecondaryVertexMass_correctedForMC = 0.0;
        
        // Fill the histogram with the allSecondaryVertexTracks_mass as it was retrieved by the CSV algorithm: 
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSoftware#FAQ (check pT-corrected secondary vertex mass)
        m_histogram.at("allSecondaryVertexTracks_mass")->Fill(jetSecondaryVertexPtCorrectedMass.at(iJet), weight);
        
        // Apply the correction for MC effects (factor 0.98) on the allSecondaryVertexTracks_mass
        if (isMC)  m_histogram.at("allSecondaryVertexTracks_correctedMass")->Fill(0.98*jetSecondaryVertexPtCorrectedMass.at(iJet), weight);
        else m_histogram.at("allSecondaryVertexTracks_correctedMass")->Fill(jetSecondaryVertexPtCorrectedMass.at(iJet), weight);
        
        // Loop over all the secondary vertices of the current analyzed event
        for(size_t iSecondaryVertex=0; iSecondaryVertex<jetSecondaryVertex.size(); ++iSecondaryVertex) {
            // Find the Secondary Vertices that belong to the same jet
            if(jetSecondaryVertexJetIndex.at(iSecondaryVertex)!=static_cast<int>(iJet)) continue;
            
            // Calculate the sum of the secondary vertex masses for the secondary vertices that belong to the same jet
            for(size_t iSecondaryVertex2=iSecondaryVertex; iSecondaryVertex2<jetSecondaryVertex.size(); ++iSecondaryVertex2){
                if(jetSecondaryVertexJetIndex.at(iSecondaryVertex2)!=jetSecondaryVertexJetIndex.at(iSecondaryVertex))continue; 
                ++secondaryVertexMultiplicityPerJet;
                SumSecondaryVertexMass = SumSecondaryVertexMass + jetSecondaryVertex.at(iSecondaryVertex2).M();
            }
            break;
        }

        m_histogram.at("secondaryVertex_multiplicityPerJet")->Fill(secondaryVertexMultiplicityPerJet, weight);
        if(secondaryVertexMultiplicityPerJet>0) {
            m_histogram.at("sumSecondaryVertex_mass")->Fill(SumSecondaryVertexMass, weight);
            SumSecondaryVertexMass_correctedForMC=0.98*SumSecondaryVertexMass;
            if (isMC)  m_histogram.at("sumSecondaryVertex_correctedMass")->Fill(SumSecondaryVertexMass_correctedForMC, weight);
            else m_histogram.at("sumSecondaryVertex_correctedMass")->Fill(SumSecondaryVertexMass, weight);
        }        
    }

    for(size_t iSecondaryVertex=0; iSecondaryVertex<jetSecondaryVertex.size(); iSecondaryVertex++) {
        m_histogram.at("secondaryVertex_ptOverJetPt")->Fill((jetSecondaryVertex.at(iSecondaryVertex).Pt())/(recoJets.at(jetSecondaryVertexJetIndex.at(iSecondaryVertex)).Pt()), weight);
        m_histogram.at("secondaryVertex_ptOverJetPt_zoom")->Fill((jetSecondaryVertex.at(iSecondaryVertex).Pt())/(recoJets.at(jetSecondaryVertexJetIndex.at(iSecondaryVertex)).Pt()), weight);
        m_histogram.at("secondaryVertex_flightDistance")->Fill(jetSecondaryVertexFlightDistanceValue.at(iSecondaryVertex), weight);
        m_histogram.at("secondaryVertex_flightDistanceSignificance")->Fill(jetSecondaryVertexFlightDistanceSignificance.at(iSecondaryVertex), weight);
    }
    
    // Calculate the secondaryVertex_trackMultiplicityPerEvent
    unsigned int nSecondaryVertexTrack = jetSecondaryVertexTrackVertexIndex.size();
    m_histogram.at("secondaryVertex_trackMultiplicityPerEvent")->Fill(nSecondaryVertexTrack,weight);
    
    // Calculate the secondaryVertexTrackMultiplicityPerJet and the secondaryVertexTrackMultiplicityPerSecondaryVertex
    for(size_t iJet=0; iJet<jetSecondaryVertexPtCorrectedMass.size(); ++iJet) {
        unsigned int secondaryVertexTrackMultiplicityPerJet = 0;
        unsigned int secondaryVertexTrackMultiplicityPerSecondaryVertex = 0;
        for(size_t iSecondaryVertex=0; iSecondaryVertex<jetSecondaryVertex.size(); ++iSecondaryVertex) {
            if(jetSecondaryVertexJetIndex.at(iSecondaryVertex)!=static_cast<int>(iJet)) continue;
            secondaryVertexTrackMultiplicityPerJet = secondaryVertexTrackMultiplicityPerSecondaryVertex;
            secondaryVertexTrackMultiplicityPerSecondaryVertex = 0;
            for(size_t iSecondaryVertexTrack=0; iSecondaryVertexTrack<jetSecondaryVertexTrackVertexIndex.size(); ++iSecondaryVertexTrack) {
                if(jetSecondaryVertexTrackVertexIndex.at(iSecondaryVertexTrack)!=static_cast<int>(iSecondaryVertex)) continue;
                for(size_t iSecondaryVertexTrack2=iSecondaryVertexTrack; iSecondaryVertexTrack2<jetSecondaryVertexTrackVertexIndex.size(); ++iSecondaryVertexTrack2) {
                    if(jetSecondaryVertexTrackVertexIndex.at(iSecondaryVertexTrack2)!=jetSecondaryVertexTrackVertexIndex.at(iSecondaryVertexTrack))continue;
                    ++secondaryVertexTrackMultiplicityPerSecondaryVertex;
                }
                break;
            }
            m_histogram.at("secondaryVertex_trackMultiplicityPerSecondaryVertex")->Fill(secondaryVertexTrackMultiplicityPerSecondaryVertex, weight);
            secondaryVertexTrackMultiplicityPerJet = secondaryVertexTrackMultiplicityPerJet + secondaryVertexTrackMultiplicityPerSecondaryVertex;
        }
        m_histogram.at("secondaryVertex_trackMultiplicityPerJet")->Fill(secondaryVertexTrackMultiplicityPerJet, weight);
    }
        
        
    // Leptons
    m_histogram.at("lepton_multiplicity")->Fill(recoObjectIndices.allLeptonIndices_.size(), weight);
    for(const int index : recoObjectIndices.leptonIndices_){
        m_histogram.at("lepton_pt")->Fill(recoObjects.allLeptons_->at(index).Pt(), weight);
        m_histogram.at("lepton_eta")->Fill(recoObjects.allLeptons_->at(index).Eta(), weight);
        m_histogram.at("lepton_phi")->Fill(recoObjects.allLeptons_->at(index).Phi(), weight);
        m_histogram.at("lepton_dxy")->Fill(recoObjects.lepDxyVertex0_->at(index), weight);
        m_histogram.at("lepton_dxy_zoom")->Fill(recoObjects.lepDxyVertex0_->at(index), weight);
        m_histogram.at("lepton_dz")->Fill(recoObjects.lepDzVertex0_->at(index), weight);
        m_histogram.at("lepton_dz_zoom")->Fill(recoObjects.lepDzVertex0_->at(index), weight);
    }
    for(const int index : recoObjectIndices.antiLeptonIndices_){
        m_histogram.at("lepton_pt")->Fill(recoObjects.allLeptons_->at(index).Pt(), weight);
        m_histogram.at("lepton_eta")->Fill(recoObjects.allLeptons_->at(index).Eta(), weight);
        m_histogram.at("lepton_phi")->Fill(recoObjects.allLeptons_->at(index).Phi(), weight);
        m_histogram.at("lepton_dxy")->Fill(recoObjects.lepDxyVertex0_->at(index), weight);
        m_histogram.at("lepton_dxy_zoom")->Fill(recoObjects.lepDxyVertex0_->at(index), weight);
        m_histogram.at("lepton_dz")->Fill(recoObjects.lepDzVertex0_->at(index), weight);
        m_histogram.at("lepton_dz_zoom")->Fill(recoObjects.lepDzVertex0_->at(index), weight);
    }
    
    const int leptonIndex = recoObjectIndices.leptonIndices_.size()>0 ? recoObjectIndices.leptonIndices_.at(0) : -1;
    const int antiLeptonIndex = recoObjectIndices.antiLeptonIndices_.size()>0 ? recoObjectIndices.antiLeptonIndices_.at(0) : -1;
    const bool hasLeptonPair = (leptonIndex!=-1 && antiLeptonIndex!=-1);
    
    // Leading lepton and antilepton
    int leadingLeptonIndex(leptonIndex);
    int nLeadingLeptonIndex(antiLeptonIndex);
    if(hasLeptonPair){
        common::orderIndices(leadingLeptonIndex, nLeadingLeptonIndex, *recoObjects.allLeptons_, common::LVpt);
        
        m_histogram.at("lepton1st_pt")->Fill(recoObjects.allLeptons_->at(leadingLeptonIndex).Pt(), weight);
        m_histogram.at("lepton1st_eta")->Fill(recoObjects.allLeptons_->at(leadingLeptonIndex).Eta(), weight);
        m_histogram.at("lepton1st_phi")->Fill(recoObjects.allLeptons_->at(leadingLeptonIndex).Phi(), weight);
        
        m_histogram.at("lepton2nd_pt")->Fill(recoObjects.allLeptons_->at(nLeadingLeptonIndex).Pt(), weight);
        m_histogram.at("lepton2nd_eta")->Fill(recoObjects.allLeptons_->at(nLeadingLeptonIndex).Eta(), weight);
        m_histogram.at("lepton2nd_phi")->Fill(recoObjects.allLeptons_->at(nLeadingLeptonIndex).Phi(), weight);
    }
    
    
    // Dilepton
    if(hasLeptonPair){
        LV dilepton(0.,0.,0.,0.);
        dilepton = recoObjects.allLeptons_->at(leadingLeptonIndex) + recoObjects.allLeptons_->at(nLeadingLeptonIndex);
        
        m_histogram.at("dilepton_mass")->Fill(dilepton.M(), weight);
        m_histogram.at("dilepton_pt")->Fill(dilepton.Pt(), weight);
        m_histogram.at("dilepton_rapidity")->Fill(dilepton.Rapidity(), weight);
        m_histogram.at("dilepton_phi")->Fill(dilepton.Phi(), weight);

        m_histogram.at("dilepton_deltaEta")->Fill(recoObjects.allLeptons_->at(leptonIndex).Eta() - recoObjects.allLeptons_->at(antiLeptonIndex).Eta(), weight);
        m_histogram.at("dilepton_deltaPhi")->Fill(ROOT::Math::VectorUtil::DeltaPhi(recoObjects.allLeptons_->at(antiLeptonIndex), recoObjects.allLeptons_->at(leptonIndex)), weight);
        m_histogram.at("dilepton_deltaDxy")->Fill(recoObjects.lepDxyVertex0_->at(leptonIndex) - recoObjects.lepDxyVertex0_->at(antiLeptonIndex), weight);
        m_histogram.at("dilepton_deltaDxy_zoom")->Fill(recoObjects.lepDxyVertex0_->at(leptonIndex) - recoObjects.lepDxyVertex0_->at(antiLeptonIndex), weight);
        m_histogram.at("dilepton_deltaDz")->Fill(recoObjects.lepDzVertex0_->at(leptonIndex) - recoObjects.lepDzVertex0_->at(antiLeptonIndex), weight);
        m_histogram.at("dilepton_deltaDz_zoom")->Fill(recoObjects.lepDzVertex0_->at(leptonIndex) - recoObjects.lepDzVertex0_->at(antiLeptonIndex), weight);
        
        m_histogram.at("dilepton_alpha")->Fill(ROOT::Math::VectorUtil::Angle(recoObjects.allLeptons_->at(antiLeptonIndex), recoObjects.allLeptons_->at(leptonIndex)), weight);
    }


    // Jets
    m_histogram.at("jet_multiplicity")->Fill(recoObjectIndices.jetIndices_.size(), weight);
    double btagDiscriminator_min = 1.1;
    for(const int index : recoObjectIndices.jetIndices_){
        m_histogram.at("jet_pt")->Fill(recoObjects.jets_->at(index).Pt(), weight);
        m_histogram.at("jet_eta")->Fill(recoObjects.jets_->at(index).Eta(), weight);
        m_histogram.at("jet_phi")->Fill(recoObjects.jets_->at(index).Phi(), weight);
        double btagDiscriminator = recoObjects.jetBTagCSV_->at(index);
        if(btagDiscriminator < -0.1) btagDiscriminator = -0.05;
        m_histogram.at("jet_btagDiscriminator")->Fill(btagDiscriminator, weight);
        m_histogram.at("jet_chargeGlobalPtWeighted")->Fill(recoObjects.jetChargeGlobalPtWeighted_->at(index), weight);
        m_histogram.at("jet_chargeRelativePtWeighted")->Fill(recoObjects.jetChargeRelativePtWeighted_->at(index), weight);
        if(btagDiscriminator < btagDiscriminator_min) btagDiscriminator_min = btagDiscriminator;
    }
    m_histogram.at("jet_btagDiscriminator_min")->Fill(btagDiscriminator_min, weight);
    
    
    // Bjets
     m_histogram.at("bjet_multiplicity")->Fill(recoObjectIndices.bjetIndices_.size(), weight);
    for(const int index : recoObjectIndices.bjetIndices_){
        m_histogram.at("bjet_pt")->Fill(recoObjects.jets_->at(index).Pt(), weight);
        m_histogram.at("bjet_eta")->Fill(recoObjects.jets_->at(index).Eta(), weight);
        m_histogram.at("bjet_phi")->Fill(recoObjects.jets_->at(index).Phi(), weight);
        m_histogram.at("bjet_chargeRelativePtWeighted")->Fill(recoObjects.jetChargeRelativePtWeighted_->at(index), weight);
    }
    
    
    // Met
    m_histogram.at("met_et")->Fill(recoObjects.met_->E(), weight);
    m_histogram.at("met_phi")->Fill(recoObjects.met_->Phi(), weight);
}








