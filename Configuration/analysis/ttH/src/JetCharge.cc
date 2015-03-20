#include <iostream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <algorithm>
#include <map>
#include <vector>

#include "TString.h"
#include "Math/VectorUtil.h"
#include "TProfile.h"
#include "TH1.h"
#include "TObjString.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "JetCharge.h"
#include "MvaVariablesJetCharge.h"
#include "MvaReaderJetCharge.h"
#include "analysisStructs.h"
#include "JetCategories.h"
#include "higgsUtils.h"
#include "HiggsAnalysis.h"
#include "../../common/include/AnalysisBase.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/classes.h"
#include "../../common/include/storeTemplate.h"
#include "../../common/include/classesFwd.h"



JetCharge::JetCharge(const bool mvaCharge, const bool correction):
mvaReader_(0),
mvaCharge_(mvaCharge),
correction_(correction),
histData_(0),
histMc_(0),
mvaReader2_(0)
{
    std::cout<<"--- Beginning setting up jet charge\n";
    
    if(mvaCharge_){
        const TString weightFilename = "/nfs/dust/cms/user/jgaray/releaseForPseudodata/CMSSW_5_3_18/src/TopAnalysis/Configuration/analysis/ttH/MVA_BDTAdaBoost_1500_beta05.weights.xml";
        
        mvaReader_ = new MvaReaderJetCharge("BDT method");
        mvaReader_->book(weightFilename);
        
        // FIXME: Replace and remove
        //Book the reader
        mvaReader2_ = new TMVA::Reader();
        TMVA::Tools::Instance(); 
        mvaReader2_->AddVariable("longChargeJet", &(jetChargeMvaStruct_.longChargeJet_));
        mvaReader2_->AddVariable("relChargeJet",&(jetChargeMvaStruct_.relChargeJet_));
        mvaReader2_->AddVariable("leadingTrackPtWeightedCharge",&(jetChargeMvaStruct_.leadingTrackPtWeightedCharge_));
        mvaReader2_->AddVariable("subleadingTrackPtWeightedCharge",&(jetChargeMvaStruct_.subleadingTrackPtWeightedCharge_));
        mvaReader2_->AddVariable("thirdleadingTrackPtWeightedCharge",&(jetChargeMvaStruct_.thirdleadingTrackPtWeightedCharge_));
        mvaReader2_->AddVariable("leadingMuonPtWeightedCharge",&(jetChargeMvaStruct_.leadingMuonPtWeightedCharge_));
        mvaReader2_->AddVariable("leadingElectronPtWeightedCharge",&(jetChargeMvaStruct_.leadingElectronPtWeightedCharge_));
        mvaReader2_->AddVariable("chargeWeightedTrackId",&(jetChargeMvaStruct_.chargeWeightedTrackId_));
        mvaReader2_->AddVariable("secondaryVertexCharge",&(jetChargeMvaStruct_.secondaryVertexCharge_));
        mvaReader2_->AddVariable("ipSignificanceLeadingTrack",&(jetChargeMvaStruct_.ipSignificanceLeadingTrack_));
        mvaReader2_->AddSpectator("trueBJetId",&(jetChargeMvaStruct_.trueBJetId_));
        TString case1 = "testingTheReader";
        mvaReader2_->BookMVA(case1, weightFilename);  
    }
    
    if(correction_){
        // FIXME: Read histograms in here, use a clone to modify it, and load it to the memory
        // FIXME: Add protections to check if file exists, and if histos exist
        histData_ = 0;
        histMc_ = 0;
        const double integralData = histData_->Integral(0, histData_->GetNbinsX()+1);
        const double integralMc = histMc_->Integral(0, histMc_->GetNbinsX()+1);
        histData_->Scale(1./integralData);
        histMc_->Scale(1./integralMc);
    }
    
    std::cout<<"=== Finishing setting up jet charge\n\n";
}



double JetCharge::pWeightedCharge(const int jetIndex, const LV& recoJet,
                                  const std::vector<int> pfCandidateTrackIndex, const VLV& pfCandidates,
                                  const std::vector<int> pfCandidateCharge, const std::vector<int>& pfCandidateVertexId,
                                  const double x)const
{
    // FIXME: where is the goodPV selection for tracks??
    
    // Access jet momentum information
    double jetPx = recoJet.px();
    double jetPy = recoJet.py();
    double jetPz = recoJet.pz();
    
    // Define relevant variables for c_{rel} calculation
    double sumMomentum = 0.;
    double sumMomentumQ = 0.;
    
    for (size_t iCandidate=0;iCandidate!=pfCandidates.size();++iCandidate)
    {
        // Check that the pfCandidate corresponds to the jet
        if (jetIndex != pfCandidateTrackIndex.at(iCandidate)) continue;
        
        // Access pfCandidate mometum and charge information
        const double constituentPx = pfCandidates.at(iCandidate).px();
        const double constituentPy = pfCandidates.at(iCandidate).py();
        const double constituentPz = pfCandidates.at(iCandidate).pz();
        const double product = constituentPx*jetPx + constituentPy*jetPy + constituentPz*jetPz;
        
        int charge = pfCandidateCharge.at(iCandidate);
        
        // Sum over all the pfCandidates
        const double productPow = std::pow(product, x);
        sumMomentum += productPow;
        sumMomentumQ += charge*productPow;
    }
    
    // Obtain the jet c_{rel}
    const double ptWeightedJetChargeXValue(sumMomentum>0 ? sumMomentumQ/sumMomentum : 0);
    return ptWeightedJetChargeXValue;
}



double JetCharge::mvaCharge(const int jetIndex, const LV& jet, const RecoObjects& recoObjects)
{
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
    
    double val1 = 0;
    bool thereIsASecondaryVertex = false;
    
    int jetHadronFlavour = -999;
    
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
    
    // Testing mva variable separation
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
    
    
    
    
    
    
    // FIXME: Replace by external variables using function mvaCharges, then remove
    
    // MVA specific variable filling
    if (jetHadronFlavour>0) jetChargeMvaStruct_.trueBJetId_ = -1;
    else if (jetHadronFlavour<0) jetChargeMvaStruct_.trueBJetId_ = 0;
    jetChargeMvaStruct_.relChargeJet_ = trueBJetRelChargeVector.at(3);
    jetChargeMvaStruct_.longChargeJet_ = trueBJetScalarChargeVector.at(3);
    jetChargeMvaStruct_.leadingTrackPtWeightedCharge_ = leadingTrackPtWeightedCharge;
    jetChargeMvaStruct_.subleadingTrackPtWeightedCharge_ = subleadingTrackPtWeightedCharge;
    jetChargeMvaStruct_.thirdleadingTrackPtWeightedCharge_ = thirdleadingTrackPtWeightedCharge;
    jetChargeMvaStruct_.leadingMuonPtWeightedCharge_ = leadingMuonPtWeightedCharge;
    jetChargeMvaStruct_.leadingElectronPtWeightedCharge_ = leadingElectronPtWeightedCharge;
    jetChargeMvaStruct_.trackNumberWeightedJetPt_ = trackNumberWeightedJetPt;
    jetChargeMvaStruct_.chargeWeightedTrackId_ = chargeWeightedTrackId;
    jetChargeMvaStruct_.svChargeWeightedFlightDistance_ = svChargeWeightedFlightDistance;
    if (thereIsASecondaryVertex) jetChargeMvaStruct_.secondaryVertexCharge_ = chargeOfSecondaryVerticesForSelectedTracks.at(0);
    else jetChargeMvaStruct_.secondaryVertexCharge_ = 0.; 
    if (impactParameterValuesForPf.size()!=0) jetChargeMvaStruct_.ipSignificanceLeadingTrack_ = impactParameterSignificanceOfLeadingTrack;
    else (jetChargeMvaStruct_.ipSignificanceLeadingTrack_ = 0.);
    
    // MVA reader variable
    val1  = mvaReader2_->EvaluateMVA ("testingTheReader");
    return val1;
}



void JetCharge::quantileMappingCorrection(double& jetCharge)const
{
    if(!correction_) return;
    
    Int_t bin = histMc_->FindFixBin(jetCharge);
    double integralAtCharge = histMc_->Integral(0, bin);
    // FIXME: put protections against problematic cases, e.g. integral=0 or =1
    if(integralAtCharge>1. || integralAtCharge<0.) std::cout<<"\n\tNormalization didn't work! Please check x!!\n";
    
    Int_t nq = 1;
    Double_t yq [1];
    Double_t xq [1];
    xq[0] = integralAtCharge;
    
    histData_->GetQuantiles(nq, yq, xq);
    
    jetCharge = yq[0];
}



std::vector<float> JetCharge::mvaCharges(const EventMetadata&,
                                         const RecoObjects& recoObjects, const CommonGenObjects&,
                                         const TopGenObjects&, const HiggsGenObjects&,
                                         const KinematicReconstructionSolutions&,
                                         const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                                         const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                         const double& weight)
{
    if(!mvaCharge_){
        std::cerr<<"Error in JetCharge::mvaCharges()! Function called, but MVA charge not set up\n...break\n"<<std::endl;
        exit(102);
    }
    
    std::vector<MvaVariablesBase*> v_mvaVariables = MvaVariablesJetCharge::fillVariables(recoObjectIndices, genObjectIndices, recoObjects, weight);
    return mvaReader_->mvaWeights(v_mvaVariables);
}






















