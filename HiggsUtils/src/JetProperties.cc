#include "TopAnalysis/HiggsUtils/interface/JetProperties.h"





JetProperties::JetProperties():
jetChargeGlobalPtWeighted_(-999.),
jetChargeRelativePtWeighted_(-999.),
jetAssociatedPartonPdgId_(-999),
jetAssociatedParton_(math::PtEtaPhiMLorentzVectorD(0.,0.,0.,0.)),
jetSecondaryVertexPtCorrectedMass_(-999.)
{}



JetProperties::JetProperties(const double& jetChargeGlobalPtWeighted, const double& jetChargeRelativePtWeighted,
                             const int& jetAssociatedPartonPdgId, const math::PtEtaPhiMLorentzVectorD& jetAssociatedParton, 
                             const std::vector<math::PtEtaPhiMLorentzVectorD>& jetPfCandidateTrack, const std::vector<int>& jetPfCandidateTrackCharge,
                             const std::vector<int>& jetPfCandidateTrackId, const std::vector<int>& jetPfCandidateTrackRelationToInteractionVertex, 
                             const std::vector<int>& jetSelectedTrackMatchToPfCandidateIndex, const std::vector<math::PtEtaPhiMLorentzVectorD>& jetSelectedTrack, 
                             const std::vector<double>& jetSelectedTrackIPValue, const std::vector<double>& jetSelectedTrackIPSignificance, 
                             const std::vector<int>& jetSelectedTrackCharge, const std::vector<int>& jetSecondaryVertexTrackMatchToSelectedTrackIndex, 
                             const std::vector<int>& jetSecondaryVertexTrackVertexIndex, const std::vector<math::PtEtaPhiMLorentzVectorD>& jetSecondaryVertex, 
                             const std::vector<double>& jetSecondaryVertexFlightDistanceValue, const std::vector<double>& jetSecondaryVertexFlightDistanceSignificance, 
                             const double& jetSecondaryVertexPtCorrectedMass):
jetChargeGlobalPtWeighted_(jetChargeGlobalPtWeighted),
jetChargeRelativePtWeighted_(jetChargeRelativePtWeighted),
jetAssociatedPartonPdgId_(jetAssociatedPartonPdgId),
jetAssociatedParton_(jetAssociatedParton),
jetPfCandidateTrack_(jetPfCandidateTrack),
jetPfCandidateTrackCharge_(jetPfCandidateTrackCharge),
jetPfCandidateTrackId_(jetPfCandidateTrackId),
jetPfCandidateTrackRelationToInteractionVertex_(jetPfCandidateTrackRelationToInteractionVertex),
jetSelectedTrackMatchToPfCandidateIndex_(jetSelectedTrackMatchToPfCandidateIndex),
jetSelectedTrack_(jetSelectedTrack),
jetSelectedTrackIPValue_(jetSelectedTrackIPValue),
jetSelectedTrackIPSignificance_(jetSelectedTrackIPSignificance),
jetSelectedTrackCharge_(jetSelectedTrackCharge),
jetSecondaryVertexTrackMatchToSelectedTrackIndex_(jetSecondaryVertexTrackMatchToSelectedTrackIndex),
jetSecondaryVertexTrackVertexIndex_(jetSecondaryVertexTrackVertexIndex),
jetSecondaryVertex_(jetSecondaryVertex),
jetSecondaryVertexFlightDistanceValue_(jetSecondaryVertexFlightDistanceValue),
jetSecondaryVertexFlightDistanceSignificance_(jetSecondaryVertexFlightDistanceSignificance),
jetSecondaryVertexPtCorrectedMass_(jetSecondaryVertexPtCorrectedMass)
{}



JetProperties::JetProperties(const JetProperties& jetProperties):
jetChargeGlobalPtWeighted_(jetProperties.jetChargeGlobalPtWeighted_),
jetChargeRelativePtWeighted_(jetProperties.jetChargeRelativePtWeighted_),
jetAssociatedPartonPdgId_(jetProperties.jetAssociatedPartonPdgId_),
jetAssociatedParton_(jetProperties.jetAssociatedParton_),
jetPfCandidateTrack_(jetProperties.jetPfCandidateTrack_), 
jetPfCandidateTrackCharge_(jetProperties.jetPfCandidateTrackCharge_),
jetPfCandidateTrackId_(jetProperties.jetPfCandidateTrackId_),
jetPfCandidateTrackRelationToInteractionVertex_(jetProperties.jetPfCandidateTrackRelationToInteractionVertex_),
jetSelectedTrackMatchToPfCandidateIndex_(jetProperties.jetSelectedTrackMatchToPfCandidateIndex_),
jetSelectedTrack_(jetProperties.jetSelectedTrack_),
jetSelectedTrackIPValue_(jetProperties.jetSelectedTrackIPValue_),
jetSelectedTrackIPSignificance_(jetProperties.jetSelectedTrackIPSignificance_),
jetSelectedTrackCharge_(jetProperties.jetSelectedTrackCharge_),
jetSecondaryVertexTrackMatchToSelectedTrackIndex_(jetProperties.jetSecondaryVertexTrackMatchToSelectedTrackIndex_),
jetSecondaryVertexTrackVertexIndex_(jetProperties.jetSecondaryVertexTrackVertexIndex_),
jetSecondaryVertex_(jetProperties.jetSecondaryVertex_),
jetSecondaryVertexFlightDistanceValue_(jetProperties.jetSecondaryVertexFlightDistanceValue_),
jetSecondaryVertexFlightDistanceSignificance_(jetProperties.jetSecondaryVertexFlightDistanceSignificance_),
jetSecondaryVertexPtCorrectedMass_(jetProperties.jetSecondaryVertexPtCorrectedMass_)
{}







const double&
JetProperties::jetChargeGlobalPtWeighted()const{return jetChargeGlobalPtWeighted_;}



const double&
JetProperties::jetChargeRelativePtWeighted()const{return jetChargeRelativePtWeighted_;}



const int&
JetProperties::jetAssociatedPartonPdgId()const{return jetAssociatedPartonPdgId_;}



const math::PtEtaPhiMLorentzVectorD&
JetProperties::jetAssociatedParton()const{return jetAssociatedParton_;}



const std::vector<math::PtEtaPhiMLorentzVectorD>&
JetProperties::jetPfCandidateTrack()const{return jetPfCandidateTrack_;}



const std::vector<int>&
JetProperties::jetPfCandidateTrackCharge()const{return jetPfCandidateTrackCharge_;}



const std::vector<int>&
JetProperties::jetPfCandidateTrackId()const{return jetPfCandidateTrackId_;}



const std::vector<int>&
JetProperties::jetPfCandidateTrackRelationToInteractionVertex()const{return jetPfCandidateTrackRelationToInteractionVertex_;}



const std::vector<int>&
JetProperties::jetSelectedTrackMatchToPfCandidateIndex()const{return jetSelectedTrackMatchToPfCandidateIndex_;}



const std::vector<math::PtEtaPhiMLorentzVectorD>&
JetProperties::jetSelectedTrack()const{return jetSelectedTrack_;}



const std::vector<double>&
JetProperties::jetSelectedTrackIPValue()const{return jetSelectedTrackIPValue_;}



const std::vector<double>&
JetProperties::jetSelectedTrackIPSignificance()const{return jetSelectedTrackIPSignificance_;}



const std::vector<int>&
JetProperties::jetSelectedTrackCharge()const{return jetSelectedTrackCharge_;}



const std::vector<int>&
JetProperties::jetSecondaryVertexTrackMatchToSelectedTrackIndex()const{return jetSecondaryVertexTrackMatchToSelectedTrackIndex_;}



const std::vector<int>&
JetProperties::jetSecondaryVertexTrackVertexIndex()const{return jetSecondaryVertexTrackVertexIndex_;}



const std::vector<math::PtEtaPhiMLorentzVectorD>&
JetProperties::jetSecondaryVertex()const{return jetSecondaryVertex_;}



const std::vector<double>&
JetProperties::jetSecondaryVertexFlightDistanceValue()const{return jetSecondaryVertexFlightDistanceValue_;}



const std::vector<double>&
JetProperties::jetSecondaryVertexFlightDistanceSignificance()const{return jetSecondaryVertexFlightDistanceSignificance_;}


const double&
JetProperties::jetSecondaryVertexPtCorrectedMass()const{return jetSecondaryVertexPtCorrectedMass_;}






