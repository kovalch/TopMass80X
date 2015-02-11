#include <map>
#include <vector>
#include <iostream>
#include <algorithm>

#include <TH1.h>
#include <TH1D.h>
#include <TString.h>
#include <Math/VectorUtil.h>

#include "AnalyzerJetProperties.h"
#include "analysisStructs.h"
#include "JetCategories.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"
#include "../../common/include/KinematicReconstructionSolution.h"





AnalyzerJetProperties::AnalyzerJetProperties(const std::vector<TString>& selectionStepsNoCategories,
                                             const std::vector<TString>& stepsForCategories,
                                             const JetCategories* jetCategories):
AnalyzerBase("jetProp_", selectionStepsNoCategories, stepsForCategories, jetCategories)
{
    std::cout<<"--- Beginning setting up basic histograms\n";
    std::cout<<"=== Finishing setting up basic histograms\n\n";
}



void AnalyzerJetProperties::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    // Secondary vertices
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
    name = "allSecondaryVertexTracks_ptCorrectedMass";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Secondary Vertex Mass (p_{t}-corrected);p_{t}-corrected M(Secondary Vertex) [GeV];Jets",80,0,8));
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
    
    // Jet charge
    name = "jet_chargeGlobalPtWeighted";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "jetChargeGlobalPtWeighted c_{glob}^{jet}; c_{glob}^{jet};# jets", 110, -1.1, 1.1));
    name = "jet_chargeRelativePtWeighted";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "jetChargeRelativePtWeighted c_{rel}^{jet}; c_{rel}^{jet};# jets", 110, -1.1, 1.1));
    name = "bjet_chargeRelativePtWeighted";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-JetChargeRelativePtWeighted c_{rel}^{jet}; c_{rel}^{jet};# B-Jets", 110, -1.1, 1.1));
}



void AnalyzerJetProperties::fillHistos(const EventMetadata&,
                                       const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                       const TopGenObjects&, const HiggsGenObjects&,
                                       const KinematicReconstructionSolutions&,
                                       const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices&,
                                       const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                       const double& weight, const TString&,
                                       std::map<TString, TH1*>& m_histogram)
{
    const std::vector<int>& secondaryVertexTrackVertexIndex = *recoObjects.jetSecondaryVertexTrackVertexIndex_;
    const std::vector<double>& secondaryVertexFlightDistance = *recoObjects.jetSecondaryVertexFlightDistanceValue_;
    const std::vector<double>& secondaryVertexFlightDistanceSignificance = *recoObjects.jetSecondaryVertexFlightDistanceSignificance_;
    const std::vector<double>& secondaryVertexPtCorrectedMass = *recoObjects.jetSecondaryVertexPtCorrectedMass_;
    const std::vector<int>& secondaryVertexJetIndex = *recoObjects.jetSecondaryVertexJetIndex_;
    const std::vector<LV>& secondaryVertex = *recoObjects.jetSecondaryVertex_;
    const VLV& jets = *recoObjects.jets_;
    const std::vector<int>& jetIndices = recoObjectIndices.jetIndices_;
    const bool isMc = commonGenObjects.valuesSet_;
    
    
    // Secondary vertices
    // Loop over jets
    unsigned int secondaryVertexMultiplicityPerEvent(0);
    unsigned int secondaryVertexTrackMultiplicityPerEvent(0);
    for(const int jetIndex : jetIndices){
        // Fill the histogram with the mass as it was retrieved by the CSV algorithm: 
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSoftware#FAQ (check pT-corrected secondary vertex mass)
        m_histogram.at("allSecondaryVertexTracks_ptCorrectedMass")->Fill(secondaryVertexPtCorrectedMass.at(jetIndex), weight);
        
        // Loop over all secondary vertices belonging to this jet
        double sumSecondaryVertexMassPerJet(0.);
        unsigned int secondaryVertexMultiplicityPerJet(0);
        unsigned int secondaryVertexTrackMultiplicityPerJet(0);
        for(size_t iSecondaryVertex = 0; iSecondaryVertex < secondaryVertex.size(); ++iSecondaryVertex){
            if(secondaryVertexJetIndex.at(iSecondaryVertex) != jetIndex) continue;
            sumSecondaryVertexMassPerJet += secondaryVertex.at(iSecondaryVertex).M();
            ++secondaryVertexMultiplicityPerEvent;
            ++secondaryVertexMultiplicityPerJet;
            
            m_histogram.at("secondaryVertex_pt")->Fill(secondaryVertex.at(iSecondaryVertex).pt(), weight);        
            m_histogram.at("secondaryVertex_mass")->Fill(secondaryVertex.at(iSecondaryVertex).M(), weight);
            // Correct for effects on MC by multiplying the secondary vertex mass of MC with 0.98
            if(isMc) m_histogram.at("secondaryVertex_correctedMass")->Fill(0.98*secondaryVertex.at(iSecondaryVertex).M(), weight);
            else m_histogram.at("secondaryVertex_correctedMass")->Fill(secondaryVertex.at(iSecondaryVertex).M(), weight);
            m_histogram.at("secondaryVertex_ptOverJetPt")->Fill((secondaryVertex.at(iSecondaryVertex).pt())/(jets.at(jetIndex).pt()), weight);
            m_histogram.at("secondaryVertex_ptOverJetPt_zoom")->Fill((secondaryVertex.at(iSecondaryVertex).pt())/(jets.at(jetIndex).pt()), weight);
            m_histogram.at("secondaryVertex_flightDistance")->Fill(secondaryVertexFlightDistance.at(iSecondaryVertex), weight);
            m_histogram.at("secondaryVertex_flightDistanceSignificance")->Fill(secondaryVertexFlightDistanceSignificance.at(iSecondaryVertex), weight);
            
            // Loop over all tracks belonging to this secondary vertex
            unsigned int secondaryVertexTrackMultiplicityPerSecondaryVertex(0);
            for(size_t iSecondaryVertexTrack = 0; iSecondaryVertexTrack < secondaryVertexTrackVertexIndex.size(); ++iSecondaryVertexTrack){
                if(secondaryVertexTrackVertexIndex.at(iSecondaryVertexTrack) != static_cast<int>(iSecondaryVertex)) continue;
                ++secondaryVertexTrackMultiplicityPerSecondaryVertex;
            }
            m_histogram.at("secondaryVertex_trackMultiplicityPerSecondaryVertex")->Fill(secondaryVertexTrackMultiplicityPerSecondaryVertex, weight);
            secondaryVertexTrackMultiplicityPerEvent += secondaryVertexTrackMultiplicityPerSecondaryVertex;
            secondaryVertexTrackMultiplicityPerJet += secondaryVertexTrackMultiplicityPerSecondaryVertex;
        }
        m_histogram.at("secondaryVertex_multiplicityPerJet")->Fill(secondaryVertexMultiplicityPerJet, weight);
        m_histogram.at("secondaryVertex_trackMultiplicityPerJet")->Fill(secondaryVertexTrackMultiplicityPerJet, weight);
        if(secondaryVertexMultiplicityPerJet > 0){
            m_histogram.at("sumSecondaryVertex_mass")->Fill(sumSecondaryVertexMassPerJet, weight);
            if(isMc) m_histogram.at("sumSecondaryVertex_correctedMass")->Fill(0.98*sumSecondaryVertexMassPerJet, weight);
            else m_histogram.at("sumSecondaryVertex_correctedMass")->Fill(sumSecondaryVertexMassPerJet, weight);
        }
    }
    m_histogram.at("secondaryVertex_multiplicityPerEvent")->Fill(secondaryVertexMultiplicityPerEvent, weight);
    m_histogram.at("secondaryVertex_trackMultiplicityPerEvent")->Fill(secondaryVertexTrackMultiplicityPerEvent, weight);
    
    
    // Jet charge
    for(const int index : recoObjectIndices.jetIndices_){
        m_histogram.at("jet_chargeGlobalPtWeighted")->Fill(recoObjects.jetChargeGlobalPtWeighted_->at(index), weight);
        m_histogram.at("jet_chargeRelativePtWeighted")->Fill(recoObjects.jetChargeRelativePtWeighted_->at(index), weight);
    }
    for(const int index : recoObjectIndices.bjetIndices_){
        m_histogram.at("bjet_chargeRelativePtWeighted")->Fill(recoObjects.jetChargeRelativePtWeighted_->at(index), weight);
    }
}








