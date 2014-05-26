#include "TopAnalysis/HiggsUtils/interface/JetProperties.h"




JetProperties::JetProperties():
jetChargeGlobalPtWeighted_(0), jetChargeRelativePtWeighted_(0),
jetAssociatedPartonPdgId_(0), jetAssociatedParton_(math::PtEtaPhiMLorentzVectorD(0,0,0,0)) 
{}
    
JetProperties::JetProperties(const double& jetChargeGlobalPtWeighted, const double& jetChargeRelativePtWeighted,
			     const int& jetAssociatedPartonPdgId, const math::PtEtaPhiMLorentzVectorD& jetAssociatedParton, 
			     const std::vector<math::PtEtaPhiMLorentzVectorD>& jetPfCandidateTrack, const std::vector<int>& jetPfCandidateTrackCharge,
                 const std::vector<int>& jetPfCandidateTrackId, const std::vector<int>& jetSelectedTrackMatchToPfCandidateIndex,
                 const std::vector<math::PtEtaPhiMLorentzVectorD>& jetSelectedTrack, 
                 const std::vector<double>& jetSelectedTrackIPValue, const std::vector<double>& jetSelectedTrackIPSignificance, 
                 const std::vector<int>& jetSelectedTrackCharge, const std::vector<int>& jetSecondaryVertexTrackSelectedTrackIndex, 
                 const std::vector<int>& jetSecondaryVertexTrackVertexIndex, const std::vector<math::PtEtaPhiMLorentzVectorD>& jetSecondaryVertex, 
                 const std::vector<double>& jetSecondaryVertexFlightDistanceValue, const std::vector<double>& jetSecondaryVertexFlightDistanceSignificance):
			     
jetChargeGlobalPtWeighted_(jetChargeGlobalPtWeighted), jetChargeRelativePtWeighted_(jetChargeRelativePtWeighted),
jetAssociatedPartonPdgId_(jetAssociatedPartonPdgId), jetAssociatedParton_(jetAssociatedParton),
jetPfCandidateTrack_(jetPfCandidateTrack), jetPfCandidateTrackCharge_(jetPfCandidateTrackCharge), jetPfCandidateTrackId_(jetPfCandidateTrackId), jetSelectedTrackMatchToPfCandidateIndex_(jetSelectedTrackMatchToPfCandidateIndex),
jetSelectedTrack_(jetSelectedTrack), jetSelectedTrackIPValue_(jetSelectedTrackIPValue), jetSelectedTrackIPSignificance_(jetSelectedTrackIPSignificance), jetSelectedTrackCharge_(jetSelectedTrackCharge),
jetSecondaryVertexTrackSelectedTrackIndex_(jetSecondaryVertexTrackSelectedTrackIndex), jetSecondaryVertexTrackVertexIndex_(jetSecondaryVertexTrackVertexIndex), jetSecondaryVertex_(jetSecondaryVertex),
jetSecondaryVertexFlightDistanceValue_(jetSecondaryVertexFlightDistanceValue), jetSecondaryVertexFlightDistanceSignificance_(jetSecondaryVertexFlightDistanceSignificance)
{}
    
JetProperties::JetProperties(const JetProperties& jetProperties):
jetChargeGlobalPtWeighted_(jetProperties.jetChargeGlobalPtWeighted_), jetChargeRelativePtWeighted_(jetProperties.jetChargeRelativePtWeighted_),
jetAssociatedPartonPdgId_(jetProperties.jetAssociatedPartonPdgId_), jetAssociatedParton_(jetProperties.jetAssociatedParton_),jetPfCandidateTrack_(jetProperties.jetPfCandidateTrack_), 
jetPfCandidateTrackCharge_(jetProperties.jetPfCandidateTrackCharge_), jetPfCandidateTrackId_(jetProperties.jetPfCandidateTrackId_), jetSelectedTrackMatchToPfCandidateIndex_(jetProperties.jetSelectedTrackMatchToPfCandidateIndex_),
jetSelectedTrack_(jetProperties.jetSelectedTrack_), jetSelectedTrackIPValue_(jetProperties.jetSelectedTrackIPValue_), jetSelectedTrackIPSignificance_(jetProperties.jetSelectedTrackIPSignificance_), jetSelectedTrackCharge_(jetProperties.jetSelectedTrackCharge_),
jetSecondaryVertexTrackSelectedTrackIndex_(jetProperties.jetSecondaryVertexTrackSelectedTrackIndex_), jetSecondaryVertexTrackVertexIndex_(jetProperties.jetSecondaryVertexTrackVertexIndex_), jetSecondaryVertex_(jetProperties.jetSecondaryVertex_),
jetSecondaryVertexFlightDistanceValue_(jetProperties.jetSecondaryVertexFlightDistanceValue_), jetSecondaryVertexFlightDistanceSignificance_(jetProperties.jetSecondaryVertexFlightDistanceSignificance_)
{}
    






double
JetProperties::jetChargeGlobalPtWeighted()const{return jetChargeGlobalPtWeighted_;}



double
JetProperties::jetChargeRelativePtWeighted()const{return jetChargeRelativePtWeighted_;}



int
JetProperties::jetAssociatedPartonPdgId()const{return jetAssociatedPartonPdgId_;}



math::PtEtaPhiMLorentzVectorD
JetProperties::jetAssociatedParton()const{return jetAssociatedParton_;}



std::vector<math::PtEtaPhiMLorentzVectorD>
JetProperties::jetPfCandidateTrack()const{return jetPfCandidateTrack_;}



std::vector<int>
JetProperties::jetPfCandidateTrackCharge()const{return jetPfCandidateTrackCharge_;}



std::vector<int>
JetProperties::jetPfCandidateTrackId()const{return jetPfCandidateTrackId_;}



std::vector<int>
JetProperties::jetSelectedTrackMatchToPfCandidateIndex()const{return jetSelectedTrackMatchToPfCandidateIndex_;}



std::vector<math::PtEtaPhiMLorentzVectorD>
JetProperties::jetSelectedTrack()const{return jetSelectedTrack_;}



std::vector<double>
JetProperties::jetSelectedTrackIPValue()const{return jetSelectedTrackIPValue_;}



std::vector<double>
JetProperties::jetSelectedTrackIPSignificance()const{return jetSelectedTrackIPSignificance_;}



std::vector<int>
JetProperties::jetSelectedTrackCharge()const{return jetSelectedTrackCharge_;}



std::vector<int>
JetProperties::jetSecondaryVertexTrackSelectedTrackIndex()const{return jetSecondaryVertexTrackSelectedTrackIndex_;}



std::vector<int>
JetProperties::jetSecondaryVertexTrackVertexIndex()const{return jetSecondaryVertexTrackVertexIndex_;}



std::vector<math::PtEtaPhiMLorentzVectorD>
JetProperties::jetSecondaryVertex()const{return jetSecondaryVertex_;}



std::vector<double>
JetProperties::jetSecondaryVertexFlightDistanceValue()const{return jetSecondaryVertexFlightDistanceValue_;}



std::vector<double>
JetProperties::jetSecondaryVertexFlightDistanceSignificance()const{return jetSecondaryVertexFlightDistanceSignificance_;}
