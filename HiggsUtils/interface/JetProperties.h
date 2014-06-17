#ifndef JetProperties_h
#define JetProperties_h



#include <vector>

#include "DataFormats/Math/interface/LorentzVector.h"




class JetProperties{
	
public:
	
    JetProperties();
	
    JetProperties(const double& jetChargeGlobalPtWeighted, const double& jetChargeRelativePtWeighted,
            const int& jetAssociatedPartonPdgId, const math::PtEtaPhiMLorentzVectorD& jetAssociatedParton,
            const std::vector<math::PtEtaPhiMLorentzVectorD>& jetPfCandidateTrack, const std::vector<int>& jetPfCandidateTrackCharge,
            const std::vector<int>& jetPfCandidateTrackId, const std::vector<int>& jetSelectedTrackMatchToPfCandidateIndex,
            const std::vector<math::PtEtaPhiMLorentzVectorD>& jetSelectedTrack, 
            const std::vector<double>& jetSelectedTrackIPValue, const std::vector<double>& jetSelectedTrackIPSignificance,
            const std::vector<int>& jetSelectedTrackCharge, const std::vector<int>& jetSecondaryVertexTrackMatchToSelectedTrackIndex,
            const std::vector<int>& jetSecondaryVertexTrackVertexIndex, const std::vector<math::PtEtaPhiMLorentzVectorD>& jetSecondaryVertex,
            const std::vector<double>& jetSecondaryVertexFlightDistanceValue, const std::vector<double>& jetSecondaryVertexFlightDistanceSignificance,
            const double& jetSecondaryVertexPtCorrectedMass
	);
	
	
	JetProperties(const JetProperties& jetProperties);
    
    double jetChargeGlobalPtWeighted()const;
    double jetChargeRelativePtWeighted()const;
    int jetAssociatedPartonPdgId()const;
    math::PtEtaPhiMLorentzVectorD jetAssociatedParton()const;
    std::vector<math::PtEtaPhiMLorentzVectorD> jetPfCandidateTrack()const;
    std::vector<int> jetPfCandidateTrackCharge()const;
    std::vector<int> jetPfCandidateTrackId()const;
    std::vector<int> jetSelectedTrackMatchToPfCandidateIndex()const;
    std::vector<math::PtEtaPhiMLorentzVectorD> jetSelectedTrack()const;
    std::vector<double> jetSelectedTrackIPValue()const;
    std::vector<double> jetSelectedTrackIPSignificance()const;
    std::vector<int> jetSelectedTrackCharge()const;
    std::vector<int> jetSecondaryVertexTrackMatchToSelectedTrackIndex()const;
    std::vector<int> jetSecondaryVertexTrackVertexIndex()const;
    std::vector<math::PtEtaPhiMLorentzVectorD> jetSecondaryVertex()const;
    std::vector<double> jetSecondaryVertexFlightDistanceValue()const;
    std::vector<double> jetSecondaryVertexFlightDistanceSignificance()const;
    double jetSecondaryVertexPtCorrectedMass()const;
    
private:
	
    double jetChargeGlobalPtWeighted_;
    double jetChargeRelativePtWeighted_;
    int jetAssociatedPartonPdgId_;
    math::PtEtaPhiMLorentzVectorD jetAssociatedParton_;
    std::vector<math::PtEtaPhiMLorentzVectorD> jetPfCandidateTrack_;
    std::vector<int> jetPfCandidateTrackCharge_;
    std::vector<int> jetPfCandidateTrackId_;
    std::vector<int> jetSelectedTrackMatchToPfCandidateIndex_;
    std::vector<math::PtEtaPhiMLorentzVectorD> jetSelectedTrack_;
    std::vector<double> jetSelectedTrackIPValue_;
    std::vector<double> jetSelectedTrackIPSignificance_;
    std::vector<int> jetSelectedTrackCharge_;
    std::vector<int> jetSecondaryVertexTrackMatchToSelectedTrackIndex_;
    std::vector<int> jetSecondaryVertexTrackVertexIndex_;
    std::vector<math::PtEtaPhiMLorentzVectorD> jetSecondaryVertex_;
    std::vector<double> jetSecondaryVertexFlightDistanceValue_;
    std::vector<double> jetSecondaryVertexFlightDistanceSignificance_;
    double jetSecondaryVertexPtCorrectedMass_;
};







#endif // JetProperties_h





