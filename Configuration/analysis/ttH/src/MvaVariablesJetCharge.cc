#include <iostream>
#include <cstdlib>
#include <algorithm>

#include <Math/VectorUtil.h>

#include "MvaVariablesJetCharge.h"
#include "analysisStructs.h"
#include "../../common/include/classes.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/analysisObjectStructs.h"





// ----------------------------------------- Methods for MvaVariablesJetCharge -----------------------------------------------------------



MvaVariablesJetCharge::MvaVariablesJetCharge():
MvaVariablesBase(),
trueBJetId_(MvaVariableInt(name_trueBJetId_)),
thereIsALeadingLepton_(MvaVariableInt(name_thereIsALeadingLepton_)),
thereIsALeadingMuon_(MvaVariableInt(name_thereIsALeadingMuon_)),
thereIsASecondaryVertex_(MvaVariableInt(name_thereIsASecondaryVertex_)),
longChargeJet_(MvaVariableFloat(name_longChargeJet_)),
relChargeJet_(MvaVariableFloat(name_relChargeJet_)),
leadingTrackPtWeightedCharge_(MvaVariableFloat(name_leadingTrackPtWeightedCharge_)),
subleadingTrackPtWeightedCharge_(MvaVariableFloat(name_subleadingTrackPtWeightedCharge_)),
thirdleadingTrackPtWeightedCharge_(MvaVariableFloat(name_thirdleadingTrackPtWeightedCharge_)),
leadingMuonPtWeightedCharge_(MvaVariableFloat(name_leadingMuonPtWeightedCharge_)),
leadingElectronPtWeightedCharge_(MvaVariableFloat(name_leadingElectronPtWeightedCharge_)),
trackNumberWeightedJetPt_(MvaVariableFloat(name_trackNumberWeightedJetPt_)),
chargeWeightedTrackId_(MvaVariableInt(name_chargeWeightedTrackId_)),
svChargeWeightedFlightDistance_(MvaVariableFloat(name_svChargeWeightedFlightDistance_)),
secondaryVertexCharge_(MvaVariableFloat(name_secondaryVertexCharge_)),
ipSignificanceLeadingTrack_(MvaVariableFloat(name_ipSignificanceLeadingTrack_))
{}



MvaVariablesJetCharge::MvaVariablesJetCharge(const int jetIndex, const RecoObjects& recoObjects, const double& eventWeight):
MvaVariablesBase(eventWeight),
trueBJetId_(MvaVariableInt(name_trueBJetId_)),
thereIsALeadingLepton_(MvaVariableInt(name_thereIsALeadingLepton_)),
thereIsALeadingMuon_(MvaVariableInt(name_thereIsALeadingMuon_)),
thereIsASecondaryVertex_(MvaVariableInt(name_thereIsASecondaryVertex_)),
longChargeJet_(MvaVariableFloat(name_longChargeJet_)),
relChargeJet_(MvaVariableFloat(name_relChargeJet_)),
leadingTrackPtWeightedCharge_(MvaVariableFloat(name_leadingTrackPtWeightedCharge_)),
subleadingTrackPtWeightedCharge_(MvaVariableFloat(name_subleadingTrackPtWeightedCharge_)),
thirdleadingTrackPtWeightedCharge_(MvaVariableFloat(name_thirdleadingTrackPtWeightedCharge_)),
leadingMuonPtWeightedCharge_(MvaVariableFloat(name_leadingMuonPtWeightedCharge_)),
leadingElectronPtWeightedCharge_(MvaVariableFloat(name_leadingElectronPtWeightedCharge_)),
trackNumberWeightedJetPt_(MvaVariableFloat(name_trackNumberWeightedJetPt_)),
chargeWeightedTrackId_(MvaVariableInt(name_chargeWeightedTrackId_)),
svChargeWeightedFlightDistance_(MvaVariableFloat(name_svChargeWeightedFlightDistance_)),
secondaryVertexCharge_(MvaVariableFloat(name_secondaryVertexCharge_)),
ipSignificanceLeadingTrack_(MvaVariableFloat(name_ipSignificanceLeadingTrack_))
{
    // Access the reco jets
    const VLV& allJets = *recoObjects.jets_; 
    
    // Specific selected tracks information (tracks with quality requirements applied already at ntuple level) 
    const std::vector<LV>& jetSelectedTrack = *recoObjects.jetSelectedTrack_;
    const std::vector<double>& jetSelectedTrackIPValue = *recoObjects.jetSelectedTrackIPValue_;
    const std::vector<double>& jetSelectedTrackIPSignificance = *recoObjects.jetSelectedTrackIPSignificance_;
    const std::vector<int>& jetSelectedTrackCharge = *recoObjects.jetSelectedTrackCharge_;
    const std::vector<int>& jetSelectedTrackIndex = *recoObjects.jetSelectedTrackIndex_;
    const std::vector<int>& jetSelectedTrackMatchToPfCandidateIndex  = *recoObjects.jetSelectedTrackMatchToPfCandidateIndex_;
    
    // Specific secondary vertex information
    const std::vector<LV>& jetSecondaryVertex = *recoObjects.jetSecondaryVertex_;
    const std::vector<int>& jetSecondaryVertexJetIndex = *recoObjects.jetSecondaryVertexJetIndex_;
    const std::vector<int>& jetSecondaryVertexTrackVertexIndex = *recoObjects.jetSecondaryVertexTrackVertexIndex_;
    const std::vector<int>& jetSecondaryVertexTrackMatchToSelectedTrackIndex = *recoObjects.jetSecondaryVertexTrackMatchToSelectedTrackIndex_;
    const std::vector<double>& jetSecondaryVertexFlightDistanceSignificance = *recoObjects.jetSecondaryVertexFlightDistanceSignificance_;
    
    // Specific track information (from pfcandidates - no quality cuts required at ntuple level)
    const std::vector<int>& jetPfCandidateTrackCharge = *recoObjects.jetPfCandidateTrackCharge_;
    const std::vector<int>& jetPfCandidateTrackId = *recoObjects.jetPfCandidateTrackId_;
    const std::vector<LV>& jetPfCandidateTrack = *recoObjects.jetPfCandidateTrack_;
    const std::vector<int>& jetPfCandidateTrackIndex = *recoObjects.jetPfCandidateTrackIndex_;
    
    bool thereIsASecondaryVertex = false;
    
    LV jet = allJets.at(jetIndex);
    
    double trueBJetPt = jet.pt();
    double jetTrueBPx = jet.px();
    double jetTrueBPy = jet.py();
    double jetTrueBPz = jet.pz();
    
    std::vector<double>sumTrueBMomentum;
    std::vector<double>sumTrueBMomentumQ;
    std::vector<double>sumTrueBMagnitude;
    std::vector<double>sumTrueBMagnitudeQ;
    
    int trueBJetTrackMultiplicity = 0;
    
    for (size_t iSumIni = 0;iSumIni!=10;iSumIni++)
    {
        sumTrueBMomentum.push_back(0);
        sumTrueBMomentumQ.push_back(0);
        sumTrueBMagnitude.push_back(0);
        sumTrueBMagnitudeQ.push_back(0);
    }
    
    double maxPtTrueTrack  = -999.;
    double leadingTrackPt = -999.;
    double leadingTrackCharge = -999.;
    double subleadingTrackPt = -999.;
    double subleadingTrackCharge = -999.;
    double thirdleadingTrackPt = -999.;
    double thirdleadingTrackCharge = -999.;
    
    std::vector<double> trueBJetMuonTracksPt;
    std::vector<LV> trueBJetMuonTracksLV;
    std::vector<double> trueBJetMuonTracksCharge;
    double trueBJetLeadingElectronCharge = 0.;
    
    std::vector<int> trackParticleId;
    std::vector<double> ptRatioValues;
    bool isLeadingMuon = false;
    bool isLeadingElectron = false;
    double ptRatioValuesElectron = -999.;
    
    std::vector<double> impactParameterValue;
    std::vector<double> impactParameterSignificance;
    std::vector<int> impactParameterMatchToPfIndices;
    
    std::vector<double> sumSVPowProduct;
    std::vector<double> sumSVPowProductQ;
    
    for (size_t iSumIni = 0;iSumIni!=10;iSumIni++)
    {
        sumSVPowProduct.push_back(0);
        sumSVPowProductQ.push_back(0);
    }
    
    for (size_t iSelectedTrack=0;iSelectedTrack!=jetSelectedTrack.size();++iSelectedTrack)
    {
        if (jetSelectedTrackIndex.at(iSelectedTrack)!=jetIndex) continue;
        
        // Store for pfCandidates
        if (jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack) != -1) 
        {
            impactParameterValue.push_back(jetSelectedTrackIPValue.at(iSelectedTrack));
            impactParameterSignificance.push_back(jetSelectedTrackIPSignificance.at(iSelectedTrack));
            impactParameterMatchToPfIndices.push_back(jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack));
        }
    }
    
    std::vector<double> secondaryVertexFlightDistanceSignificance;
    unsigned int secondaryVertexMultiplicityPerJet = 0; 
    
    // Secondary vertex multiplicity
    for(size_t iSecondaryVertex=0; iSecondaryVertex<jetSecondaryVertex.size(); ++iSecondaryVertex) 
    {
        if(jetSecondaryVertexJetIndex.at(iSecondaryVertex)!=static_cast<int>(jetIndex)) continue;
        
        for(size_t iSecondaryVertex2=iSecondaryVertex; iSecondaryVertex2<jetSecondaryVertex.size(); ++iSecondaryVertex2)
        {
            if(jetSecondaryVertexJetIndex.at(iSecondaryVertex2)!=jetSecondaryVertexJetIndex.at(iSecondaryVertex)) continue;
            ++secondaryVertexMultiplicityPerJet;
        }
        break; 
    }
    
    std::vector<double> chargeOfSecondaryVerticesForSelectedTracks;
    
    for(size_t jSecondaryVertex=0; jSecondaryVertex<jetSecondaryVertex.size(); ++jSecondaryVertex) 
    {
        // Check that SV belongs to the jet
        if(jetSecondaryVertexJetIndex.at(jSecondaryVertex)!=static_cast<int>(jetIndex)) continue;
        // Ony continue if at least one SV on the jet
        if (secondaryVertexMultiplicityPerJet==0) continue; 
        // Access flight distance information
        secondaryVertexFlightDistanceSignificance.push_back(jetSecondaryVertexFlightDistanceSignificance.at(jSecondaryVertex));
        for (size_t iSelectedTrack=0;iSelectedTrack!=jetSelectedTrack.size();++iSelectedTrack)
        {
            // Check that track belongs to jet
            if (jetSelectedTrackIndex.at(iSelectedTrack)!=jetIndex) continue;
            
            // Check that track belongs to a SV
            std::vector<int>::const_iterator isInVector = std::find(jetSecondaryVertexTrackMatchToSelectedTrackIndex.begin(), jetSecondaryVertexTrackMatchToSelectedTrackIndex.end(),iSelectedTrack);
            if (isInVector ==  jetSecondaryVertexTrackMatchToSelectedTrackIndex.end()) continue;
            
            // Check that track belongs to the SV - if so, sum up for track multiplicity
            if (jetSecondaryVertexTrackVertexIndex.at(isInVector-jetSecondaryVertexTrackMatchToSelectedTrackIndex.begin()) != (int) jSecondaryVertex) continue;
            
            // Calculate secondary vertex charge
            const double svProduct = jetSelectedTrack.at(iSelectedTrack).px()*jetSecondaryVertex.at(jSecondaryVertex).px() +jetSelectedTrack.at(iSelectedTrack).py()*jetSecondaryVertex.at(jSecondaryVertex).py() + jetSelectedTrack.at(iSelectedTrack).pz()*jetSecondaryVertex.at(jSecondaryVertex).pz();
            
            std::vector<double> svProductPow;
            for (double iPow=0.2;iPow<=2.;iPow+=0.2)
            {
                svProductPow.push_back(std::pow(svProduct,iPow));
            }
            
            for (size_t i_sum=0;i_sum!=sumSVPowProduct.size();i_sum++)
            {
                sumSVPowProduct.at(i_sum) += svProductPow.at(i_sum);
                sumSVPowProductQ.at(i_sum) += (svProductPow.at(i_sum))*jetSelectedTrackCharge.at(iSelectedTrack);
            }
            
        }
        const double secondaryVertexChargeTest(sumSVPowProduct.at(4)>0 ? sumSVPowProductQ.at(4)/sumSVPowProduct.at(4) : 0);
        chargeOfSecondaryVerticesForSelectedTracks.push_back(secondaryVertexChargeTest);  
    }
    
    if (secondaryVertexMultiplicityPerJet==1) thereIsASecondaryVertex = true;
    
    std::vector<double> impactParameterValuesForPf;
    double impactParameterSignificanceOfLeadingTrack = -999.;
    
    for (size_t iPfTrack=0;iPfTrack!=jetPfCandidateTrack.size();iPfTrack++)
    {
        double trueBJetPfTrackPt = jetPfCandidateTrack.at(iPfTrack).pt();
        int trueMatched = 0;
        if (jetIndex!=jetPfCandidateTrackIndex.at(iPfTrack)) trueMatched = -1;
        if (trueMatched == -1) continue;
        ++trueBJetTrackMultiplicity;
        
        if(trueBJetPfTrackPt>maxPtTrueTrack) maxPtTrueTrack = trueBJetPfTrackPt;
        
        // Access impact parameter
        double impactParameterValueForPf = 0;
        double impactParameterSignificanceForPf = 0;
        
        std::vector<int>::const_iterator pfIsMatchedToSelTrack = std::find(impactParameterMatchToPfIndices.begin(),impactParameterMatchToPfIndices.end(),iPfTrack);
        if (pfIsMatchedToSelTrack!=impactParameterMatchToPfIndices.end()) 
        {
            impactParameterSignificanceForPf = impactParameterSignificance.at(pfIsMatchedToSelTrack - impactParameterMatchToPfIndices.begin());
            impactParameterValuesForPf.push_back(impactParameterValueForPf);
            if (iPfTrack==0) impactParameterSignificanceOfLeadingTrack = impactParameterSignificanceForPf;
        }
        
        // Check if the track is a lepton or not:
        int particleId = jetPfCandidateTrackId.at(iPfTrack);
        trackParticleId.push_back(particleId);
        
        if (particleId==3) 
        {
            trueBJetMuonTracksPt.push_back(trueBJetPfTrackPt);
            trueBJetMuonTracksCharge.push_back(jetPfCandidateTrackCharge.at(iPfTrack));
        }
        
        double ptRatio = trueBJetPfTrackPt/trueBJetPt;
        ptRatioValues.push_back(ptRatio);
        
        if (ptRatioValues.size()==1)
        {
            leadingTrackPt = jetPfCandidateTrack.at(iPfTrack).pt();
            leadingTrackCharge = jetPfCandidateTrackCharge.at(iPfTrack);
        }
        
        if (ptRatioValues.size()==1&&particleId==2) 
        {
            ptRatioValuesElectron = ptRatio;
            trueBJetLeadingElectronCharge = jetPfCandidateTrackCharge.at(iPfTrack);
            isLeadingElectron = true;
        }
        
        if (ptRatioValues.size()==1&&particleId==3) 
        {
            isLeadingMuon = true;
        }
        
        if (ptRatioValues.size()==2&&particleId==2) 
        {
            subleadingTrackPt = jetPfCandidateTrack.at(iPfTrack).pt();
            subleadingTrackCharge = jetPfCandidateTrackCharge.at(iPfTrack);
        }
        
        if (ptRatioValues.size()==2&&particleId==3) 
        {
            subleadingTrackPt = jetPfCandidateTrack.at(iPfTrack).pt();
            subleadingTrackCharge = jetPfCandidateTrackCharge.at(iPfTrack);
        }
        
        // Calculate the jet c_{rel} 
        const double constituentTrueBPx = jetPfCandidateTrack.at(iPfTrack).px();
        const double constituentTrueBPy = jetPfCandidateTrack.at(iPfTrack).py();
        const double constituentTrueBPz = jetPfCandidateTrack.at(iPfTrack).pz();
        const double trueProduct = constituentTrueBPx*jetTrueBPx + constituentTrueBPy*jetTrueBPy + constituentTrueBPz*jetTrueBPz;
        
        std::vector<double> vectProductMomentumQ;
        double xTrueComponent = (jetTrueBPy*constituentTrueBPz-jetTrueBPz*constituentTrueBPy);
        double yTrueComponent = (jetTrueBPx*constituentTrueBPz-jetTrueBPz*constituentTrueBPx);
        double zTrueComponent = (jetTrueBPx*constituentTrueBPy-jetTrueBPy*constituentTrueBPx); 
        const double trueMagnitude = std::sqrt(xTrueComponent*xTrueComponent+yTrueComponent*yTrueComponent+zTrueComponent*zTrueComponent);
        
        std::vector<double> trueProductPow;
        std::vector<double> trueMagnitudePow;
        
        for (double iPow=0.2;iPow<=2.;iPow+=0.2)
        {
            trueProductPow.push_back(std::pow(trueProduct,iPow));
            trueMagnitudePow.push_back(std::pow(trueMagnitude,iPow));
        }
        
        for (size_t i_sum=0;i_sum!=sumTrueBMomentum.size();i_sum++)
        {
            sumTrueBMomentum.at(i_sum) += trueProductPow.at(i_sum);
            sumTrueBMagnitude.at(i_sum) += trueMagnitudePow.at(i_sum);
            sumTrueBMomentumQ.at(i_sum) += (trueProductPow.at(i_sum))*jetPfCandidateTrackCharge.at(iPfTrack);
            sumTrueBMagnitudeQ.at(i_sum) += (trueMagnitudePow.at(i_sum))*jetPfCandidateTrackCharge.at(iPfTrack);
        }
    }
    
    std::vector<double> trueBJetScalarChargeVector;
    std::vector<double> trueBJetRelChargeVector;
    
    for (size_t iSumHis=0;iSumHis!=sumTrueBMomentum.size();iSumHis++)
    {
        const double trueBJetScalarChargeForX(sumTrueBMomentum.at(iSumHis)>0 ? sumTrueBMomentumQ.at(iSumHis)/sumTrueBMomentum.at(iSumHis) : 0);
        const double trueBJetRelCharge(sumTrueBMagnitude.at(iSumHis)>0 ? sumTrueBMagnitudeQ.at(iSumHis)/sumTrueBMagnitude.at(iSumHis) : 0);
        trueBJetScalarChargeVector.push_back(trueBJetScalarChargeForX);
        trueBJetRelChargeVector.push_back(trueBJetRelCharge);
    }
    
    double leadingTrackPtWeightedCharge = leadingTrackCharge*leadingTrackPt/trueBJetPt;
    double subleadingTrackPtWeightedCharge = subleadingTrackCharge*subleadingTrackPt/trueBJetPt;
    double thirdleadingTrackPtWeightedCharge = thirdleadingTrackCharge*thirdleadingTrackPt/trueBJetPt;
    double leadingMuonPtWeightedCharge = 0.;
    if (isLeadingMuon) leadingMuonPtWeightedCharge = trueBJetMuonTracksCharge.at(0)*trueBJetMuonTracksPt.at(0)/trueBJetPt;
    double leadingElectronPtWeightedCharge = 0.;
    if (isLeadingElectron) leadingElectronPtWeightedCharge = trueBJetLeadingElectronCharge*ptRatioValuesElectron/trueBJetPt;
    double trackNumberWeightedJetPt = leadingTrackCharge*trueBJetPt/trueBJetTrackMultiplicity;
    double chargeWeightedTrackId = leadingTrackCharge;
    if (isLeadingMuon) chargeWeightedTrackId = leadingTrackCharge*3;
    double svChargeWeightedFlightDistance = 0.;
    if (thereIsASecondaryVertex) svChargeWeightedFlightDistance = chargeOfSecondaryVerticesForSelectedTracks.at(0) * jetSecondaryVertexFlightDistanceSignificance.at(0); 
    
    longChargeJet_.value_ = trueBJetScalarChargeVector.at(3);
    relChargeJet_.value_ = trueBJetRelChargeVector.at(3);
    leadingTrackPtWeightedCharge_.value_ = leadingTrackPtWeightedCharge;
    subleadingTrackPtWeightedCharge_.value_ = subleadingTrackPtWeightedCharge;
    thirdleadingTrackPtWeightedCharge_.value_ = thirdleadingTrackPtWeightedCharge;
    leadingMuonPtWeightedCharge_.value_ = leadingMuonPtWeightedCharge;
    leadingElectronPtWeightedCharge_.value_ = leadingElectronPtWeightedCharge;
    trackNumberWeightedJetPt_.value_ = trackNumberWeightedJetPt;
    chargeWeightedTrackId_.value_ = chargeWeightedTrackId;
    svChargeWeightedFlightDistance_.value_ = svChargeWeightedFlightDistance;
    if (impactParameterValuesForPf.size()!=0) ipSignificanceLeadingTrack_.value_ = impactParameterSignificanceOfLeadingTrack;
    else ipSignificanceLeadingTrack_.value_ = 0.;
    if (thereIsASecondaryVertex) secondaryVertexCharge_.value_ = chargeOfSecondaryVerticesForSelectedTracks.at(0);
    else secondaryVertexCharge_.value_ = 0;
    if (isLeadingMuon||isLeadingElectron) thereIsALeadingLepton_.value_ = 1;
    else thereIsALeadingLepton_.value_ = 0;
    if (isLeadingMuon) thereIsALeadingMuon_.value_ = 1;
    else thereIsALeadingMuon_ = 0;
    if (thereIsASecondaryVertex) thereIsASecondaryVertex_.value_ = 1;
    else thereIsASecondaryVertex_.value_ = 0;
}



std::vector<MvaVariablesBase*> MvaVariablesJetCharge::fillVariables(const tth::RecoObjectIndices& recoObjectIndices,
                                                                    const tth::GenObjectIndices&,
                                                                    const RecoObjects& recoObjects,
                                                                    const double& eventWeight)
{
    std::vector<MvaVariablesBase*> result;
    
    const std::vector<int>& jetIndices = recoObjectIndices.jetIndices_;
    
    
    // Loop over all jet pairs
    for(int index : jetIndices){
        // Is it the last entry of an event
        const bool lastInEvent = index == static_cast<int>(jetIndices.size()-1);
        double jetIndex = jetIndices.at(index);
        
        MvaVariablesBase* mvaVariables = new MvaVariablesJetCharge(jetIndex, recoObjects, eventWeight);
        
        result.push_back(mvaVariables);
    }
    
    return result;
}







// ----------------------------------------- Methods for MvaVariablesJetChargePerEvent -----------------------------------------------------------



MvaVariablesJetChargePerEvent::MvaVariablesJetChargePerEvent(const std::vector<MvaVariablesBase*>& v_mvaVariables):
v_mvaVariables_(v_mvaVariables)
{}



MvaVariablesJetChargePerEvent::MvaVariablesJetChargePerEvent(const std::vector<MvaVariablesBase*>& v_mvaVariables,
                                                         const std::map<std::string, std::vector<float> >& m_mvaWeight):
v_mvaVariables_(v_mvaVariables),
m_weight_(m_mvaWeight)
{
    for(const auto& weights : m_weight_){
        if(weights.second.size() != v_mvaVariables_.size()){
            std::cerr<<"ERROR in constructor of MvaVariablesJetChargePerEvent! Vector sizes do not match for weights (variables, weights): "
                     <<v_mvaVariables.size()<<" , "<<weights.second.size()<<"\n...break\n"<<std::endl;
            exit(47);
        }
    }
}



size_t MvaVariablesJetChargePerEvent::maxWeightIndex(const std::string& mvaConfigName)const
{
    const std::vector<float>& v_weight = this->mvaWeights(mvaConfigName);
    return common::extremumIndex(v_weight);
}



float MvaVariablesJetChargePerEvent::maxWeight(const std::string& mvaConfigName)const
{
    const size_t maxWeightIndex = this->maxWeightIndex(mvaConfigName);
    return m_weight_.at(mvaConfigName).at(maxWeightIndex);
}



std::vector<MvaVariablesBase*> MvaVariablesJetChargePerEvent::variables()const
{
    return v_mvaVariables_;
}



std::vector<float> MvaVariablesJetChargePerEvent::mvaWeights(const std::string& mvaConfigName)const
{
    if(m_weight_.find(mvaConfigName) == m_weight_.end()){
        std::cerr<<"ERROR in mvaWeights()! No weights found for mvaConfigName: "
                 <<mvaConfigName<<"\n...break\n"<<std::endl;
        exit(50);
    }
    
    return m_weight_.at(mvaConfigName);
}



std::map<std::string, std::vector<float> > MvaVariablesJetChargePerEvent::mvaWeightsMap()const
{
    return m_weight_;
}







