#include <iostream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <algorithm>
#include <map>
#include <vector>

#include <TString.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <Math/VectorUtil.h>
#include <TProfile.h>
#include <TObjString.h>
#include <TTree.h>
#include <TFile.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "AnalyzerJetCharge.h"
#include "analysisStructs.h"
#include "JetCategories.h"
#include "higgsUtils.h"
#include "HiggsAnalysis.h"
#include "../../common/include/AnalysisBase.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/classes.h"
#include "../../common/include/KinematicReconstructionSolution.h"
#include "../../common/include/storeTemplate.h"

AnalyzerJetCharge::AnalyzerJetCharge(const std::vector<TString>& selectionStepsNoCategories,
                                     const std::vector<TString>& stepsForCategories,
                                     const JetCategories* jetCategories,
                                     const int debug):
AnalyzerBase ("jetCharge_", selectionStepsNoCategories, stepsForCategories, jetCategories),
reader_(0), 
case1_(0),
debug_(0)
{
    debug_ = debug;
    std::cout<<"--- Beginning setting up jetChargeAnalyzer\n";
    
    //reader_ = new TMVA::Reader();
    //TMVA::Tools::Instance(); 
    
    //TString case1_ = "testingTheReader";
    //TString weight1 = "/data/user/jgaray/cmsswFullSetup_14Patch1/CMSSW_5_3_14_patch1/src/TopAnalysis/Configuration/analysis/ttH/factoryOutput/261114/invertSignalAndBkg/weights/MVA_BDTAdaBoost.weights.xml";
    
    
    //reader_->AddVariable("longChargeJet", &(mvaStruct_.longChargeJet_));
    //reader_->AddVariable("relChargeJet",&(mvaStruct_.relChargeJet_));
    //reader_->AddVariable("leadingTrackPtWeightedCharge",&(mvaStruct_.leadingTrackPtWeightedCharge_));
    //reader_->AddVariable("leadingMuonPtWeightedCharge",&(mvaStruct_.leadingMuonPtWeightedCharge_));
    //reader_->AddVariable("trackNumberWeightedJetPt",&(mvaStruct_.trackNumberWeightedJetPt_));
    //reader_->AddVariable("chargeWeightedTrackId",&(mvaStruct_.chargeWeightedTrackId_));
    //reader_->AddVariable("svChargeWeightedFlightDistance",&(mvaStruct_.svChargeWeightedFlightDistance_));
    //reader_->AddVariable("secondaryVertexCharge",&(mvaStruct_.secondaryVertexCharge_));
    //reader_->AddVariable("ipSignificanceLeadingTrack",&(mvaStruct_.ipSignificanceLeadingTrack_));
    
    //reader_->AddSpectator("trueBJetId",&(mvaStruct_.trueBJetId_));
    
    //reader_->BookMVA(case1_, weight1); 
    
    std::cout<<"=== Finishing setting up jetChargeAnalyzer\n\n";
}
                                     


void AnalyzerJetCharge::fillHistos(const EventMetadata& eventMetadata,
                                   const RecoObjects& recoObjects, const CommonGenObjects&,
                                   const TopGenObjects& topGenObjects, const HiggsGenObjects&,
                                   const KinematicReconstructionSolutions&,
                                   const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                                   const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                   const double& weight, const TString&,
                                   std::map<TString, TH1*>& m_histogram)
  
{
    TString name;

    // Stearing options for: true (0), kinReco (1), mlb (2)
    int optionForCalibration = -1; 
    
    // Extracting input data to more comfortable variables
    
    // Gen jets
    const std::vector<std::vector<int> >& genJetBhadronIndices = genObjectIndices.genJetBhadronIndices_; 
    const std::vector<int>& genJetMatchedRecoBjetIndices = genObjectIndices.genJetMatchedRecoBjetIndices_;
    
    // Reco jets
    const VLV& allJets = *recoObjects.jets_; 
    const std::vector<int>& lowerPtCUTJetIdx = recoObjectIndices.jetIndices_;           // Selected jets (point to jets from allJets)
    const std::vector<double>& jetChargeRelativePtWeighted = *recoObjects.jetChargeRelativePtWeighted_;
    //const int& recoBjetFromTopIndex = genObjectIndices.recoBjetFromTopIndex_;
    //const int& recoAntiBjetFromTopIndex = genObjectIndices.recoAntiBjetFromTopIndex_;
    const std::vector<double>& jetBTagCSV = *recoObjects.jetBTagCSV_;
    
    // b-hadron + c-hadron information
    const std::vector<int>& bHadFlavour = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadFlavour_ : std::vector<int>(0);

    // Specific selected tracks information (tracks with quality requirements applied already at ntuple level) 
    const std::vector<LV>& jetSelectedTrack = *recoObjects.jetSelectedTrack_;
    const std::vector<double>& jetSelectedTrackIPValue = *recoObjects.jetSelectedTrackIPValue_;
    const std::vector<double>& jetSelectedTrackIPSignificance = *recoObjects.jetSelectedTrackIPSignificance_;
    const std::vector<int>& jetSelectedTrackCharge = *recoObjects.jetSelectedTrackCharge_;
    const std::vector<int>& jetSelectedTrackIndex = *recoObjects.jetSelectedTrackIndex_;
    const std::vector<int>& jetSelectedTrackMatchToPfCandidateIndex  = *recoObjects.jetSelectedTrackMatchToPfCandidateIndex_;
    
    // Specific track information (from pfcandidates - no quality cuts required at ntuple level)
    const std::vector<int>& jetPfCandidateTrackCharge = *recoObjects.jetPfCandidateTrackCharge_;
    const std::vector<int>& jetPfCandidateTrackId = *recoObjects.jetPfCandidateTrackId_;
    const std::vector<LV>& jetPfCandidateTrack = *recoObjects.jetPfCandidateTrack_;
    const std::vector<int>& jetPfCandidateTrackIndex = *recoObjects.jetPfCandidateTrackIndex_;
    
    // Specific secondary vertex information
    const std::vector<LV>& jetSecondaryVertex = *recoObjects.jetSecondaryVertex_;
    const std::vector<int>& jetSecondaryVertexJetIndex = *recoObjects.jetSecondaryVertexJetIndex_;
    const std::vector<int>& jetSecondaryVertexTrackVertexIndex = *recoObjects.jetSecondaryVertexTrackVertexIndex_;
    const std::vector<int>& jetSecondaryVertexTrackMatchToSelectedTrackIndex = *recoObjects.jetSecondaryVertexTrackMatchToSelectedTrackIndex_;
    const std::vector<double>& jetSecondaryVertexFlightDistanceValue = *recoObjects.jetSecondaryVertexFlightDistanceValue_;
    const std::vector<double>& jetSecondaryVertexFlightDistanceSignificance = *recoObjects.jetSecondaryVertexFlightDistanceSignificance_;
    
    // Vertex information for closesZVertex studies
    const std::vector<int>& jetPfCandidatePrimaryVertexId = *recoObjects.jetPfCandidatePrimaryVertexId_;
    
    // Access the leptons for mlb
    const VLV& isolatedLeptons = *recoObjects.allLeptons_;
    const std::vector<int>& isolatedLeptonsPgdId = *recoObjects.lepPdgId_;
    
    // Event related variables
    const UInt_t& eventNumber = eventMetadata.eventNumber_;

    // =============================================================Now do calculations and filling of histograms===================================================================
  
    
    // Create new indices for the tracks->cut on different p_{T} values
    constexpr double ptTracksCUT = 0.;
    
    // Pt-order the pf tracks (even if they are in principle already ordered)
    std::vector<int> ptOrderedJetTrackIdx = common::initialiseIndices(jetPfCandidateTrack);
    common::orderIndices(ptOrderedJetTrackIdx, jetPfCandidateTrack, common::LVpt);
    
    // Pt-order the selected tracks (they are not p_{T} ordered by default)
    std::vector<int> ptOrderedSelTrackIdx = common::initialiseIndices(jetSelectedTrack);
    common::orderIndices(ptOrderedSelTrackIdx, jetSelectedTrack, common::LVpt);

    std::vector<int> systemLepton0;
    std::vector<int> systemLepton1;
   
    // Mlb related calculations
    if (optionForCalibration==2)
    {
        //if (isolatedLeptons.size() == 2&&bJetsId.size()==2)
        if (isolatedLeptons.size() == 2&&lowerPtCUTJetIdx.size()==2)    
        {
            for (size_t iLepton=0; iLepton!=isolatedLeptons.size(); ++iLepton)
            {
                double massForMlb = 170.;
                //if (iLepton==0) systemLepton0 = leptonToJetMlbCalculator(massForMlb,isolatedLeptons.at(0),isolatedLeptonsPgdId.at(0),allJets,bJetsId);
                //else if (iLepton==1) systemLepton1 = leptonToJetMlbCalculator(massForMlb,isolatedLeptons.at(1),isolatedLeptonsPgdId.at(1),allJets,bJetsId);
                
                if (iLepton==0) systemLepton0 = leptonToJetMlbCalculator(massForMlb,isolatedLeptons.at(0),isolatedLeptonsPgdId.at(0),allJets,lowerPtCUTJetIdx);
                else if (iLepton==1) systemLepton1 = leptonToJetMlbCalculator(massForMlb,isolatedLeptons.at(1),isolatedLeptonsPgdId.at(1),allJets,lowerPtCUTJetIdx);
            }
        }
    }
    
    for(size_t iJet=0;iJet!=lowerPtCUTJetIdx.size();++iJet)
    //for(size_t iJet=0;iJet!=bJetsId.size();++iJet)
    {
        int jetIdx = lowerPtCUTJetIdx.at(iJet);
        //int jetIdx = bJetsId.at(iJet);
        LV jets = allJets.at(jetIdx);
        
        //FIXME weightReweighted is reweighted by the number of tracks
        //double weightReweighted = weight*trackMultiplicityWeight (-0.039234, 1.41540, jetIdx, jetPfCandidateTrackIndex, jetPfCandidatePrimaryVertexId);
        
        if (optionForCalibration==2)
        {
            //if (isolatedLeptons.size() != 2||bJetsId.size()!=2) continue;
            if (isolatedLeptons.size() != 2||lowerPtCUTJetIdx.size()!=2) continue;
            if (systemLepton0.at(0)==0||systemLepton1.at(0)==0) continue;
            if (systemLepton0.at(1) == systemLepton1.at(1)) continue;
        }
        
        // Mva specific tree and booleans to specify event type
        bool fillTree = false;
        bool thereIsASecondaryVertex = false;
        
        //Hyerarchical values for the c_{rel}
        double jetChargeHyerarchicalValue = -999.;
        bool jetChargeHyerarchicalValueWasCalculated = false;
        bool secondaryVertexChargeWasCalculated = false;
        
        int indexOfJetAssociatedToProperIsolatedLepton = -1;
        
        if (optionForCalibration==2)
        {
        // Select only b or bbar
        if (systemLepton0.at(0)<0) indexOfJetAssociatedToProperIsolatedLepton = systemLepton0.at(1);
        if (systemLepton1.at(0)<0) indexOfJetAssociatedToProperIsolatedLepton = systemLepton1.at(1);
        
        //if (systemLepton0.at(0)>0) indexOfJetAssociatedToProperIsolatedLepton = systemLepton0.at(1);
        //if (systemLepton1.at(0)>0) indexOfJetAssociatedToProperIsolatedLepton = systemLepton1.at(1);
        
        if (jetIdx!=indexOfJetAssociatedToProperIsolatedLepton) continue;
        if (indexOfJetAssociatedToProperIsolatedLepton==-1) continue;
        }
        
        //if (optionForCalibration==1)
        //{
            //if (kinRecoAntiBIndex.size()==0) continue;
            //if (jetIdx!=kinRecoAntiBIndex.at(0)) continue;
            
            //if (kinRecoBIndex.size()==0) continue;
            //if (jetIdx!=kinRecoBIndex.at(0)) continue;
        //}
        
        //bool recoBFromTop = false;
        //bool recoAntiBFromTop = false;
        int numHadMatched = -1;
        //bool recoBorAntiBFromAny =  false;
        
        int jetHadronFlavour = -999;
        
        if (optionForCalibration==0)
        {
            //GEN TO RECO matching
            //recoBorAntiBFromAny = std::find(genJetMatchedRecoBjetIndices.begin(),genJetMatchedRecoBjetIndices.end(), jetIdx) != genJetMatchedRecoBjetIndices.end();
            
            //if (recoBjetFromTopIndex==jetIdx) recoBFromTop = true;
            //if (recoAntiBjetFromTopIndex == jetIdx) recoAntiBFromTop = true;
            
            // If the reco jet is not matched to a gen jet, continue
            std::vector <int>::const_iterator recoPositionForGenJet = std::find(genJetMatchedRecoBjetIndices.begin(),genJetMatchedRecoBjetIndices.end(), jetIdx);
            int recoJetPosition = (recoPositionForGenJet == genJetMatchedRecoBjetIndices.end()) ? -1 : recoPositionForGenJet - genJetMatchedRecoBjetIndices.begin();
            if (recoJetPosition==-1) continue;
            
            // Choose only unambiguous cases in which there's only one hadron associated to the genJet
            numHadMatched = genJetBhadronIndices.at(recoJetPosition).size();
            if (numHadMatched>1) continue;
            // Store the hadron associated to the genJet to which the reco jet is associated
           for (size_t iFlavour = 0;iFlavour!=genJetBhadronIndices.at(recoJetPosition).size();++iFlavour)
           {
               int jetHadronFlavourIndex = genJetBhadronIndices.at(recoJetPosition).at(iFlavour);
               jetHadronFlavour = bHadFlavour.at(jetHadronFlavourIndex);
           }
           if (jetHadronFlavour==0) continue;
           //if (jetHadronFlavour<0) continue;
           //if (jetHadronFlavour>0) continue;
        }
        
        double trueBJetScalarCharge = jetChargeRelativePtWeighted.at(jetIdx);
        double trueBJetPt = jets.pt();
        double trueBJetEta = jets.eta();
        double jetTrueBPx = jets.px();
        double jetTrueBPy = jets.py();
        double jetTrueBPz = jets.pz();
        
        // Add some control plots
        m_histogram["h_trueBJetPtInitially"]->Fill(trueBJetPt, weight);
        m_histogram["h_trueBJetCSVvalue"]->Fill(jetBTagCSV.at(jetIdx), weight);
        
        int trueBJetTrackMultiplicity = 0;
        
        std::vector<double>sumTrueBMomentum;
        std::vector<double>sumTrueBMomentumQ;
        std::vector<double>sumTrueBMagnitude;
        std::vector<double>sumTrueBMagnitudeQ;
        
        for (size_t iSumIni = 0;iSumIni!=10;iSumIni++)
        {
            sumTrueBMomentum.push_back(0);
            sumTrueBMomentumQ.push_back(0);
            sumTrueBMagnitude.push_back(0);
            sumTrueBMagnitudeQ.push_back(0);
        }
        
        double maxPtTrueTrack  = -999.;
        double maxTrueMagnitude = -999.;
        double maxTrueProduct = -999.;
        double leadingTrackPt = -999.;
        double leadingTrackCharge = -999.;
        double leadingTrackPx = -999.;
        double leadingTrackPy = -999.;
        double leadingTrackPz = -999.;
        double subleadingTrackPt = -999.;
        double subleadingTrackCharge = -999.;
        double thirdleadingTrackPt = -999.;
        double thirdleadingTrackCharge = -999.;
        
        //lepton-tracks variables
        std::vector<double> trueBJetLeptonTracksPt;
        std::vector<LV> trueBJetLeptonTracksLV;
        std::vector<double> trueBJetLeptonTracksCharge;
        
        std::vector<double> trueBJetMuonTracksPt;
        std::vector<LV> trueBJetMuonTracksLV;
        std::vector<double> trueBJetMuonTracksCharge;
        
        std::vector<double> trueBJetElectronTracksPt;
        std::vector<LV> trueBJetElectronTracksLV;
        std::vector<double> trueBJetElectronTracksCharge;
        
        bool leptonHasNegativeCharge = false;
        bool muonHasNegativeCharge = false;
        bool electronHasNegativeCharge = false;
        std::vector<int> trackParticleId;
        std::vector<double> ptRatioValues;
        
        double ptRatioValuesMuon = -999.;
        double ptRatioValuesElectron = -999.;
        double ptRatioValuesLepton = -999.;
        double ptRatioValuesSubleadingMuon = -999.;
        double ptRatioValuesSubleadingElectron = -999.;
        double ptRatioValuesSubleadingLepton = -999.;
        
        bool isLeadingLepton = false;
        bool isLeadingMuon = false;
        bool isLeadingElectron = false;
        bool isSubleadingLepton = false;
        bool isSubleadingMuon = false;
        bool isSubleadingElectron = false;
        //bool isThirdleadingLepton = false;
        //bool isThirdleadingElectron = false;
        //bool isThirdleadingMuon = false;
        bool isNonLeadingLepton = false;
        
        // Access selected tracks
        double sumSelectedPwoProduct = 0.;
        double sumSelectedPwoProductQ = 0;
        
        std::vector<double> impactParameterValue;
        std::vector<double> impactParameterSignificance;
        std::vector<int> impactParameterMatchToPfIndices;
        double maxImpactParameterValue = -999.;
        double maxImpactParameterValueForPf = -999.;
        double maxImpactParameterSignificance = -999999.;
        double maxImpactParameterSignificanceForPf = -999999.;
        int maxImpactParameterValueForPfIndex = -1;
        
        for (size_t iSelectedTrack=0;iSelectedTrack!=jetSelectedTrack.size();++iSelectedTrack)
        {
            if (jetSelectedTrackIndex.at(iSelectedTrack)!=jetIdx) continue;
            
            // Jet c_{rel} calculation using selectedTracks
            const double product = jetSelectedTrack.at(iSelectedTrack).px()*jetTrueBPx +jetSelectedTrack.at(iSelectedTrack).py()*jetTrueBPy + jetSelectedTrack.at(iSelectedTrack).pz()*jetTrueBPz;
            const double powProduct = std::pow (product,0.8);
            sumSelectedPwoProduct += powProduct;
            sumSelectedPwoProductQ += jetSelectedTrackCharge.at(iSelectedTrack)*powProduct;
            
            // Access Impact parameter
            m_histogram["h_trueBJetTrackIPValueForSelTracks"]->Fill(jetSelectedTrackIPValue.at(iSelectedTrack), weight);
            m_histogram["h_trueBJetTrackIPSignificanceForSelTracks"]->Fill(jetSelectedTrackIPSignificance.at(iSelectedTrack), weight);
            
            if (jetSelectedTrackIPValue.at(iSelectedTrack)>maxImpactParameterValue) maxImpactParameterValue = jetSelectedTrackIPValue.at(iSelectedTrack);
            if (jetSelectedTrackIPSignificance.at(iSelectedTrack)>maxImpactParameterSignificance) maxImpactParameterSignificance = jetSelectedTrackIPSignificance.at(iSelectedTrack);
            
            // Store for pfCandidates
            if (jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack) != -1) 
            {
                impactParameterValue.push_back(jetSelectedTrackIPValue.at(iSelectedTrack));
                impactParameterSignificance.push_back(jetSelectedTrackIPSignificance.at(iSelectedTrack));
                impactParameterMatchToPfIndices.push_back(jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack));
                
                // Max impact parameter for pfCandidates
                if (jetSelectedTrackIPValue.at(iSelectedTrack)>maxImpactParameterValueForPf)
                {
                    maxImpactParameterValueForPf = jetSelectedTrackIPValue.at(iSelectedTrack);
                    //maxImpactParameterSignificanceForPf = jetSelectedTrackIPSignificance.at(iSelectedTrack);
                    maxImpactParameterValueForPfIndex = jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack);
                }
            }
        }
        
        m_histogram["h_trueBJetTrackMaxIPValueForSelTracks"]->Fill(maxImpactParameterValue, weight);
        m_histogram["h_trueBJetTrackMaxIPSignificanceForSelTracks"]->Fill(maxImpactParameterSignificance, weight);
        if (maxImpactParameterValueForPfIndex!=-1) m_histogram["h_trueBJetTrackMaxIPValueForPfCandidates"]->Fill(maxImpactParameterValueForPf, weight);
        if(maxImpactParameterValueForPfIndex!=-1) m_histogram["h_trueBJetTrackMaxIPSignificanceForPfCandidates"]->Fill(maxImpactParameterSignificanceForPf, weight);
        
        const double selectedTrackJetCharge(sumSelectedPwoProduct>0 ? sumSelectedPwoProductQ/sumSelectedPwoProduct : 0);
        m_histogram["h_selectedTrackJetCharge"]->Fill(selectedTrackJetCharge, weight);
        
        // Access secondary vertex information
        unsigned int secondaryVertexMultiplicityPerJet = calculateMultiplicity(jetSecondaryVertexJetIndex, jetIdx);
        
        m_histogram["h_trueBJetTrackSecondaryVertexMultiplicity"]->Fill(secondaryVertexMultiplicityPerJet, weight);
        
        std::vector<double> sumSVPowProduct;
        std::vector<double> sumSVPowProductQ;
        
        for (size_t iSumIni = 0;iSumIni!=10;iSumIni++)
        {
            sumSVPowProduct.push_back(0);
            sumSVPowProductQ.push_back(0);
        }
        
        std::vector<double> secondaryVertexFlightDistanceValue;
        std::vector<double> secondaryVertexFlightDistanceSignificance;
        double minSecondaryVertexFlightDistanceValue = 99999.; 
        
        std::vector<double> chargeOfSecondaryVerticesForSelectedTracks;
        std::vector<double> chargeOfSecondaryVerticesForSelectedTracksPfMatched;
        std::vector<double> chargeOfSecondaryVerticesForSelectedTracksNonPfMatched;
        std::vector<int> multiplicityOfSecondaryVerticesForSelectedTracks;
        std::vector<int> multiplicityOfSecondaryVerticesForSelectedTracksPfMatched;
        std::vector<int> multiplicityOfSecondaryVerticesForSelectedTracksNonPfMatched;
        std::vector<int> nonWeightedSvSelTrackChargeVector;
        
        for(size_t jSecondaryVertex=0; jSecondaryVertex<jetSecondaryVertex.size(); ++jSecondaryVertex) 
        {
            // Check that SV belongs to the jet
            if(jetSecondaryVertexJetIndex.at(jSecondaryVertex)!=static_cast<int>(jetIdx)) continue;
            
            // Ony continue if at least one SV on the jet
            if (secondaryVertexMultiplicityPerJet==0) continue; 
            
            // Access flight distance information
            secondaryVertexFlightDistanceValue.push_back(jetSecondaryVertexFlightDistanceValue.at(jSecondaryVertex));
            m_histogram["h_secondaryVertexFlighDistanceValue"]->Fill(jetSecondaryVertexFlightDistanceValue.at(jSecondaryVertex), weight);
            if (secondaryVertexMultiplicityPerJet == 1) m_histogram["h_secondaryVertexFlighDistanceValueFirstSV"]->Fill(jetSecondaryVertexFlightDistanceValue.at(jSecondaryVertex), weight);
            else if (secondaryVertexMultiplicityPerJet == 2) 
            {
                m_histogram["h_secondaryVertexFlighDistanceValueSecondSV"]->Fill(jetSecondaryVertexFlightDistanceValue.at(jSecondaryVertex), weight);
                m_histogram["h_secondaryVertexFlighDistanceValueFirstSVIfJetTwoSv"]->Fill(jetSecondaryVertexFlightDistanceValue.at(0), weight);
                m_histogram["h_secondaryVertexFlighDistanceValueSecondSVIfJetTwoSv"]->Fill(jetSecondaryVertexFlightDistanceValue.at(1), weight);
                m_histogram["h_secondaryVertexFlighDistanceSignificanceFirstSVIfJetTwoSv"]->Fill(jetSecondaryVertexFlightDistanceSignificance.at(0), weight);
                m_histogram["h_secondaryVertexFlighDistanceSignificanceSecondSVIfJetTwoSv"]->Fill(jetSecondaryVertexFlightDistanceSignificance.at(1), weight);
            }
            else if (secondaryVertexMultiplicityPerJet == 3) 
            {
                m_histogram["h_secondaryVertexFlighDistanceValueThirdSV"]->Fill(jetSecondaryVertexFlightDistanceValue.at(jSecondaryVertex), weight);
                m_histogram["h_secondaryVertexFlighDistanceValueFirstSVIfJetThreeSv"]->Fill(jetSecondaryVertexFlightDistanceValue.at(0), weight);
                m_histogram["h_secondaryVertexFlighDistanceValueSecondSVIfJetThreeSv"]->Fill(jetSecondaryVertexFlightDistanceValue.at(1), weight);
                m_histogram["h_secondaryVertexFlighDistanceValueThirdSVIfJetThreeSv"]->Fill(jetSecondaryVertexFlightDistanceValue.at(2), weight);
                m_histogram["h_secondaryVertexFlighDistanceSignificanceFirstSVIfJetThreeSv"]->Fill(jetSecondaryVertexFlightDistanceSignificance.at(0), weight);
                m_histogram["h_secondaryVertexFlighDistanceSignificanceSecondSVIfJetThreeSv"]->Fill(jetSecondaryVertexFlightDistanceSignificance.at(1), weight);
                m_histogram["h_secondaryVertexFlighDistanceSignificanceThirdSVIfJetThreeSv"]->Fill(jetSecondaryVertexFlightDistanceSignificance.at(2), weight);
            }
            else if (secondaryVertexMultiplicityPerJet == 4)
            {
                m_histogram["h_secondaryVertexFlighDistanceValueFirstSVIfJetFourSv"]->Fill(jetSecondaryVertexFlightDistanceValue.at(0), weight);
                m_histogram["h_secondaryVertexFlighDistanceValueSecondSVIfJetFourSv"]->Fill(jetSecondaryVertexFlightDistanceValue.at(1), weight);
                m_histogram["h_secondaryVertexFlighDistanceValueThirdSVIfJetFourSv"]->Fill(jetSecondaryVertexFlightDistanceValue.at(2), weight);
                m_histogram["h_secondaryVertexFlighDistanceValueFourthSVIfJetFourSv"]->Fill(jetSecondaryVertexFlightDistanceValue.at(3), weight);
                m_histogram["h_secondaryVertexFlighDistanceSignificanceFirstSVIfJetFourSv"]->Fill(jetSecondaryVertexFlightDistanceSignificance.at(0), weight);
                m_histogram["h_secondaryVertexFlighDistanceSignificanceSecondSVIfJetFourSv"]->Fill(jetSecondaryVertexFlightDistanceSignificance.at(1), weight);
                m_histogram["h_secondaryVertexFlighDistanceSignificanceThirdSVIfJetFourSv"]->Fill(jetSecondaryVertexFlightDistanceSignificance.at(2), weight);
                m_histogram["h_secondaryVertexFlighDistanceSignificanceFourthSVIfJetFourSv"]->Fill(jetSecondaryVertexFlightDistanceSignificance.at(3), weight);
            }
            
            
            secondaryVertexFlightDistanceSignificance.push_back(jetSecondaryVertexFlightDistanceSignificance.at(jSecondaryVertex));
            m_histogram["h_secondaryVertexFlighDistanceSignificance"]->Fill(jetSecondaryVertexFlightDistanceSignificance.at(jSecondaryVertex), weight);
            if (secondaryVertexMultiplicityPerJet == 1) m_histogram["h_secondaryVertexFlighDistanceSignificanceFirstSV"]->Fill(jetSecondaryVertexFlightDistanceSignificance.at(jSecondaryVertex), weight);
            else if (secondaryVertexMultiplicityPerJet == 2) m_histogram["h_secondaryVertexFlighDistanceSignificanceSecondSV"]->Fill(jetSecondaryVertexFlightDistanceSignificance.at(jSecondaryVertex), weight);
            else if (secondaryVertexMultiplicityPerJet == 3) m_histogram["h_secondaryVertexFlighDistanceSignificanceThirdSV"]->Fill(jetSecondaryVertexFlightDistanceSignificance.at(jSecondaryVertex), weight);
            
            if (jetSecondaryVertexFlightDistanceValue.at(jSecondaryVertex)<minSecondaryVertexFlightDistanceValue) minSecondaryVertexFlightDistanceValue = jetSecondaryVertexFlightDistanceValue.at(jSecondaryVertex);
            
            double sumSVPowProductPfMatched = 0.;
            double sumSVPowProductQPfMatched = 0;
            
            double sumSVPowProductNonPfMatched = 0.;
            double sumSVPowProductQNonPfMatched = 0;
            
            int secondaryVertexTrackMultiplicity = 0;
            int secondaryVertexTrackMultiplicityPfMatched = 0;
            int secondaryVertexTrackMultiplicityNonPfMatched = 0;
            
            int nonWeightedSvSelTrackCharge = 0;
            
            for (size_t iSelectedTrack=0;iSelectedTrack!=jetSelectedTrack.size();++iSelectedTrack)
            {
                // Check that track belongs to jet
                if (jetSelectedTrackIndex.at(iSelectedTrack)!=jetIdx) continue;
                
                // Check that track belongs to a SV
                std::vector<int>::const_iterator isInVector = std::find(jetSecondaryVertexTrackMatchToSelectedTrackIndex.begin(), jetSecondaryVertexTrackMatchToSelectedTrackIndex.end(),iSelectedTrack);
                if (isInVector ==  jetSecondaryVertexTrackMatchToSelectedTrackIndex.end()) continue;
                
                // Check that track belongs to the SV - if so, sum up for track multiplicity
                if (jetSecondaryVertexTrackVertexIndex.at(isInVector-jetSecondaryVertexTrackMatchToSelectedTrackIndex.begin()) != (int) jSecondaryVertex) continue;
                ++secondaryVertexTrackMultiplicity;
                
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
                
                // Calculate non-weighted secondary vertex charge
                nonWeightedSvSelTrackCharge += jetSelectedTrackCharge.at(iSelectedTrack);
                
                // Secondary vertex charge for selected tracks matched to a pfCandidate - in other words: for pfCandidates
                if (jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack)!=-1)
                {
                    // Calculate secondary vertex charge
                    const double productPfMatched = jetSelectedTrack.at(iSelectedTrack).px()*jetSecondaryVertex.at(jSecondaryVertex).px() +jetSelectedTrack.at(iSelectedTrack).py()*jetSecondaryVertex.at(jSecondaryVertex).py() + jetSelectedTrack.at(iSelectedTrack).pz()*jetSecondaryVertex.at(jSecondaryVertex).pz();
                    const double powProductPfMatched = std::pow (productPfMatched,0.8);
                    sumSVPowProductPfMatched += powProductPfMatched;
                    sumSVPowProductQPfMatched += jetSelectedTrackCharge.at(iSelectedTrack)*powProductPfMatched;
                    
                    // Secondary vertex multiplicity
                    ++secondaryVertexTrackMultiplicityPfMatched;
                    
                    //Secondary vertex tracks Id
                    m_histogram["h_trueBJetTrackSecondaryVertexTrackId"]->Fill(jetPfCandidateTrackId.at(jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack)), weight);
                    m_histogram["h_trueBJetTrackSecondaryVertexTrackPt"]->Fill(jetPfCandidateTrack.at(jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack)).pt(), weight);
                    if (jetPfCandidateTrackId.at(jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack)) == 3)  m_histogram["h_trueBJetTrackSecondaryVertexMuonTrackPt"]->Fill(jetPfCandidateTrack.at(jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack)).pt(), weight);
                    if (jetPfCandidateTrackId.at(jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack)) == 2)  m_histogram["h_trueBJetTrackSecondaryVertexElectronTrackPt"]->Fill(jetPfCandidateTrack.at(jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack)).pt(), weight);
                    if (jetPfCandidateTrackId.at(jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack)) == 3 || jetPfCandidateTrackId.at(jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack)) == 2) m_histogram["h_trueBJetTrackSecondaryVertexLeptonTrackPt"]->Fill(jetPfCandidateTrack.at(jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack)).pt(), weight);
                    if (jetPfCandidateTrackId.at(jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack)) == 1)  m_histogram["h_trueBJetTrackSecondaryVertexNonLeptonTrackPt"]->Fill(jetPfCandidateTrack.at(jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack)).pt(), weight);
                }
                
                // Secondary vertex charge for selected tracks not matched to any pfCandidate - just a test procedure
                else if (jetSelectedTrackMatchToPfCandidateIndex.at(iSelectedTrack)==-1)
                {
                    const double productNonPfMatched = jetSelectedTrack.at(iSelectedTrack).px()*jetSecondaryVertex.at(jSecondaryVertex).px() +jetSelectedTrack.at(iSelectedTrack).py()*jetSecondaryVertex.at(jSecondaryVertex).py() + jetSelectedTrack.at(iSelectedTrack).pz()*jetSecondaryVertex.at(jSecondaryVertex).pz();
                    const double powProductNonPfMatched = std::pow (productNonPfMatched,0.8);
                    sumSVPowProductNonPfMatched += powProductNonPfMatched;
                    sumSVPowProductQNonPfMatched += jetSelectedTrackCharge.at(iSelectedTrack)*powProductNonPfMatched;
                    ++secondaryVertexTrackMultiplicityNonPfMatched;
                }
            }
            const double secondaryVertexChargeTest(sumSVPowProduct.at(4)>0 ? sumSVPowProductQ.at(4)/sumSVPowProduct.at(4) : 0);
            const double secondaryVertexChargePfMatched(sumSVPowProductPfMatched>0 ? sumSVPowProductQPfMatched/sumSVPowProductPfMatched : 0);
            const double secondaryVertexChargeNonPfMatched(sumSVPowProductNonPfMatched>0 ? sumSVPowProductQNonPfMatched/sumSVPowProductNonPfMatched : 0);
            chargeOfSecondaryVerticesForSelectedTracks.push_back(secondaryVertexChargeTest);   
            chargeOfSecondaryVerticesForSelectedTracksPfMatched.push_back(secondaryVertexChargePfMatched);
            chargeOfSecondaryVerticesForSelectedTracksNonPfMatched.push_back(secondaryVertexChargeNonPfMatched);
            multiplicityOfSecondaryVerticesForSelectedTracks.push_back(secondaryVertexTrackMultiplicity);
            multiplicityOfSecondaryVerticesForSelectedTracksPfMatched.push_back(secondaryVertexTrackMultiplicityPfMatched);
            multiplicityOfSecondaryVerticesForSelectedTracksNonPfMatched.push_back(secondaryVertexTrackMultiplicityNonPfMatched);
            nonWeightedSvSelTrackChargeVector.push_back(nonWeightedSvSelTrackCharge);
            
            // Create a vector of histograms: contains histograms for x=0.2 to x=2.0
            std::vector<TString> svChargeHistogramScalarVector;
            
            for (double iHisto=2;iHisto<=20;iHisto+=2)
            {
                std::stringstream ss_histoScalar;
                
                ss_histoScalar<<"h_trueBJetTrackSecondaryVertexSelectedTrackCharge"<<iHisto;
                TString histoSV = ss_histoScalar.str();
                svChargeHistogramScalarVector.push_back(histoSV);
            }
            
            // Fill the histograms with the corresponding value of the charge (we take values from the sumTrueBMagnitude and sumTrueBMomentum vectors)
            std::vector<double> svChargeVector;
            
            for (size_t iSumHis=0;iSumHis!=sumSVPowProduct.size();iSumHis++)
            {
                const double secondaryVertexChargeForX(sumSVPowProduct.at(iSumHis)>0 ? sumSVPowProductQ.at(iSumHis)/sumSVPowProduct.at(iSumHis) : 0);
                TString histoSV = svChargeHistogramScalarVector.at(iSumHis);
                
                m_histogram[histoSV]->Fill(secondaryVertexChargeForX, weight);
                
                svChargeVector.push_back(secondaryVertexChargeForX);
            }
        }
         
        if (secondaryVertexMultiplicityPerJet==1) thereIsASecondaryVertex = true;
        
        // Fill secondary vertex charge and multiplicity information 
        for (size_t iFillMult=0;iFillMult!=chargeOfSecondaryVerticesForSelectedTracks.size();++iFillMult)
        {
            m_histogram["h_trueBJetTrackSecondaryVertexSelectedTrackCharge"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(iFillMult), weight);
            m_histogram["h_trueBJetSecondaryVertexSelectedTrackNonWeightedCharge"]->Fill(nonWeightedSvSelTrackChargeVector.at(iFillMult), weight);
            if (chargeOfSecondaryVerticesForSelectedTracks.at(iFillMult)>0&&nonWeightedSvSelTrackChargeVector.at(iFillMult)>0) m_histogram["h_trueBJetTrackSecondaryVertexSelectedTrackChargeCoincidence"]->Fill(1);
            else if (chargeOfSecondaryVerticesForSelectedTracks.at(iFillMult)<0&&nonWeightedSvSelTrackChargeVector.at(iFillMult)<0) m_histogram["h_trueBJetTrackSecondaryVertexSelectedTrackChargeCoincidence"]->Fill(1);
            else m_histogram["h_trueBJetTrackSecondaryVertexSelectedTrackChargeCoincidence"]->Fill(0);
            m_histogram["h_trueBJetTrackSecondaryVertexSelectedTrackChargePfMatched"]->Fill(chargeOfSecondaryVerticesForSelectedTracksPfMatched.at(iFillMult), weight);
            m_histogram["h_trueBJetTrackSecondaryVertexSelectedTrackChargeNonPfMatched"]->Fill(chargeOfSecondaryVerticesForSelectedTracksNonPfMatched.at(iFillMult), weight);
            m_histogram["h_trueBJetTrackSecondaryVertexSelectedTrackMultiplicity"]->Fill(multiplicityOfSecondaryVerticesForSelectedTracks.at(iFillMult), weight);
            m_histogram["h_trueBJetTrackSecondaryVertexSelectedTrackMultiplicityPfMatched"]->Fill(multiplicityOfSecondaryVerticesForSelectedTracksPfMatched.at(iFillMult), weight);
            m_histogram["h_trueBJetTrackSecondaryVertexSelectedTrackMultiplicityNonPfMatched"]->Fill(multiplicityOfSecondaryVerticesForSelectedTracksNonPfMatched.at(iFillMult), weight);
            if (secondaryVertexMultiplicityPerJet==1) m_histogram["h_trueBJetTrackSecondaryVertexSelectedTrackCharge_IfOneSvInJet"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(iFillMult), weight);
        }
        
        // Sv charge for different sv track multiplicities
        for (size_t svCharge=0;svCharge!=chargeOfSecondaryVerticesForSelectedTracks.size();++svCharge)
        {
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==2) m_histogram["h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity2"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==3) m_histogram["h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity3"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==4) m_histogram["h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity4"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==5) m_histogram["h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity5"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==6) m_histogram["h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity6"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==7) m_histogram["h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity7"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==8) m_histogram["h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity8"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==9) m_histogram["h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity9"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==10) m_histogram["h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity10"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==11) m_histogram["h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity11"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==12) m_histogram["h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity12"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==13) m_histogram["h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity13"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==14) m_histogram["h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity14"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(svCharge), weight);
            
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==2) m_histogram["h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity2"]->Fill(nonWeightedSvSelTrackChargeVector.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==3) m_histogram["h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity3"]->Fill(nonWeightedSvSelTrackChargeVector.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==4) m_histogram["h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity4"]->Fill(nonWeightedSvSelTrackChargeVector.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==5) m_histogram["h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity5"]->Fill(nonWeightedSvSelTrackChargeVector.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==6) m_histogram["h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity6"]->Fill(nonWeightedSvSelTrackChargeVector.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==7) m_histogram["h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity7"]->Fill(nonWeightedSvSelTrackChargeVector.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==8) m_histogram["h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity8"]->Fill(nonWeightedSvSelTrackChargeVector.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==9) m_histogram["h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity9"]->Fill(nonWeightedSvSelTrackChargeVector.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==10) m_histogram["h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity10"]->Fill(nonWeightedSvSelTrackChargeVector.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==11) m_histogram["h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity11"]->Fill(nonWeightedSvSelTrackChargeVector.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==12) m_histogram["h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity12"]->Fill(nonWeightedSvSelTrackChargeVector.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==13) m_histogram["h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity13"]->Fill(nonWeightedSvSelTrackChargeVector.at(svCharge), weight);
            if (multiplicityOfSecondaryVerticesForSelectedTracks.at(svCharge)==14) m_histogram["h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity14"]->Fill(nonWeightedSvSelTrackChargeVector.at(svCharge), weight);
            
        }
        
        // Min. flight distance plots
        if (secondaryVertexMultiplicityPerJet!=0) m_histogram["h_secondaryVertexMinFlighDistanceValue"]->Fill(minSecondaryVertexFlightDistanceValue, weight);
        ((TH2D*)m_histogram["h_secondaryVertexMinFlightDistanceValueVsCSV"])->Fill(minSecondaryVertexFlightDistanceValue,jetBTagCSV.at(jetIdx), weight);
        
        // Secondary vertex information for pfCandidates
        std::vector<int> secondaryVertexSelTrackMatchedToPfCandidateIndex;
        std::vector<int> secondaryVertexSelTrackMatchedToPfCandidateVertexIndex;
        std::vector<int> secondaryVertexSelTrackMatchedToPfCandidateCharge;
        std::vector<LV> secondaryVertexSelTrackMatchedToPfCandidateLV;
        std::vector<int> secondaryVertexSelTrackMatchedToPfCandidateParticleId;
        
        std::vector<double> impactParameterValuesForPf;
        double impactParameterSignificanceOfLeadingTrack = -999.;
        
        for (size_t iPfTrack=0;iPfTrack!=jetPfCandidateTrack.size();++iPfTrack)
        {
            //require only tracks above a certain threshold of pt
            double trueBJetPfTrackPt = jetPfCandidateTrack.at(iPfTrack).pt();
            if (trueBJetPfTrackPt < ptTracksCUT) continue;
            
            //check if the track is matched to a selected jet and in case it is, add one to the multiplicity.
            int trueMatched = 0;
            if (jetIdx!=jetPfCandidateTrackIndex.at(iPfTrack)) trueMatched = -1;
            if (trueMatched == -1) continue;
            
            m_histogram["h_pfCandidateInitialAmount"]->Fill(0., weight);
            if (jetHadronFlavour>0) m_histogram["h_pfCandidateInitialAmountForB"]->Fill(0., weight);
            else if (jetHadronFlavour<0) m_histogram["h_pfCandidateInitialAmountForAntiB"]->Fill(0., weight);
            
            // TEST zvertex studies
            // Check if PfCandidate belongs to primary vertex
            if (jetPfCandidatePrimaryVertexId.at(iPfTrack)==0&&jetPfCandidatePrimaryVertexId.at(iPfTrack)==3) m_histogram["h_pfCandidateNonZeroWeights"]->Fill(0., weight);
            if (jetPfCandidatePrimaryVertexId.at(iPfTrack)==3) m_histogram["h_pfCandidateNonZeroWeightsIndexNonZero"]->Fill(0., weight);
            if (jetPfCandidatePrimaryVertexId.at(iPfTrack)==0) m_histogram["h_pfCandidateNonZeroWeightsIndexZero"]->Fill(0., weight);
            if (jetPfCandidatePrimaryVertexId.at(iPfTrack)==2) m_histogram["h_pfCandidateZeroWeightsIndexNonZero"]->Fill(0., weight);
            if (jetPfCandidatePrimaryVertexId.at(iPfTrack)==1) m_histogram["h_pfCandidateZeroWeightsIndexZero"]->Fill(0., weight);
            if (jetPfCandidatePrimaryVertexId.at(iPfTrack)==-1) std::cout<<"TOO MANY HIGH WEIGHTS!!"<<std::endl;
            
            // Remove tracks not associated to primary vertex 0
            if (jetPfCandidatePrimaryVertexId.at(iPfTrack)==-1|| jetPfCandidatePrimaryVertexId.at(iPfTrack)==2 || jetPfCandidatePrimaryVertexId.at(iPfTrack)==3) continue;
            
            //if (jetHadronFlavour>0) m_histogram["h_test_jetPfCandidateChargeForB"]->Fill(jetPfCandidateTrackCharge.at(iPfTrack), weight);
            //else if (jetHadronFlavour<0) m_histogram["h_test_jetPfCandidateChargeForAntiB"]->Fill(jetPfCandidateTrackCharge.at(iPfTrack), weight);
            
            ++trueBJetTrackMultiplicity;
            
            if(trueBJetPfTrackPt>maxPtTrueTrack) 
            {
                maxPtTrueTrack = trueBJetPfTrackPt;
            }
            
            // Access impact parameter
            double impactParameterValueForPf = 0;
            double impactParameterSignificanceForPf = 0;
            bool ipValueWasDefined = false;
            
            std::vector<int>::const_iterator pfIsMatchedToSelTrack = std::find(impactParameterMatchToPfIndices.begin(),impactParameterMatchToPfIndices.end(),iPfTrack);
            if (pfIsMatchedToSelTrack!=impactParameterMatchToPfIndices.end()) 
            {
                impactParameterValueForPf = impactParameterValue.at(pfIsMatchedToSelTrack - impactParameterMatchToPfIndices.begin());
                impactParameterSignificanceForPf = impactParameterSignificance.at(pfIsMatchedToSelTrack - impactParameterMatchToPfIndices.begin());
                m_histogram["h_trueBJetTrackIPValueForPfCandidates"]->Fill(impactParameterValueForPf, weight);
                m_histogram["h_trueBJetTrackIPSignificanceForPfCandidates"]->Fill(impactParameterSignificanceForPf, weight);
                ipValueWasDefined = true;
                impactParameterValuesForPf.push_back(impactParameterValueForPf);
                if (iPfTrack==0) impactParameterSignificanceOfLeadingTrack = impactParameterSignificanceForPf;
            }
            
            // Check if the track is a lepton or not:
            int particleId = jetPfCandidateTrackId.at(iPfTrack);
            trackParticleId.push_back(particleId);
            
            //if track is a muon (3), fill some muon-specific histograms
            if (particleId==3) 
            {
                trueBJetMuonTracksPt.push_back(trueBJetPfTrackPt);
                trueBJetMuonTracksLV.push_back(jetPfCandidateTrack.at(iPfTrack));
                trueBJetMuonTracksCharge.push_back(jetPfCandidateTrackCharge.at(iPfTrack));
                m_histogram["h_trueBJetMuonTrackPt"]->Fill(trueBJetPfTrackPt, weight);
                m_histogram["h_trueBJetMuonTrackEta"]->Fill(jetPfCandidateTrack.at(iPfTrack).Eta(), weight);
                m_histogram["h_trueBJetMuonTrackCharge"]->Fill(jetPfCandidateTrackCharge.at(iPfTrack), weight);
                ((TH2D*)m_histogram["h_trueBJetMuonChargePt"])->Fill(trueBJetPfTrackPt,jetPfCandidateTrackCharge.at(iPfTrack), weight);
                m_histogram["h_trueBJetToMuonTrackPtRatio"]->Fill(trueBJetPfTrackPt/trueBJetPt, weight);
                if (ipValueWasDefined) m_histogram["h_trueBJetTrackIPValueIfMuon"]->Fill(impactParameterValueForPf, weight);
                if (ipValueWasDefined) m_histogram["h_trueBJetTrackIPSignificanceIfMuon"]->Fill(impactParameterSignificanceForPf, weight);
            }
            
            //if track is an electron (2), fill some electron-specific histograms
            else if (particleId==2)
            {
                trueBJetElectronTracksPt.push_back(trueBJetPfTrackPt);
                trueBJetElectronTracksLV.push_back(jetPfCandidateTrack.at(iPfTrack));
                trueBJetElectronTracksCharge.push_back(jetPfCandidateTrackCharge.at(iPfTrack));
                m_histogram["h_trueBJetElectronTrackPt"]->Fill(trueBJetPfTrackPt, weight);
                m_histogram["h_trueBJetElectronTrackEta"]->Fill(jetPfCandidateTrack.at(iPfTrack).Eta(), weight);
                m_histogram["h_trueBJetElectronTrackCharge"]->Fill(jetPfCandidateTrackCharge.at(iPfTrack), weight);
                ((TH2D*)m_histogram["h_trueBJetElectronChargePt"])->Fill(trueBJetPfTrackPt,jetPfCandidateTrackCharge.at(iPfTrack), weight);
                m_histogram["h_trueBJetToElectronTrackPtRatio"]->Fill(trueBJetPfTrackPt/trueBJetPt, weight);
                if (ipValueWasDefined) m_histogram["h_trueBJetTrackIPValueIfElectron"]->Fill(impactParameterValueForPf, weight);
                if (ipValueWasDefined) m_histogram["h_trueBJetTrackIPSignificanceIfElectron"]->Fill(impactParameterSignificanceForPf, weight);
            }
            
            if (particleId==2 || particleId==3)
            {
                trueBJetLeptonTracksPt.push_back(trueBJetPfTrackPt);
                trueBJetLeptonTracksLV.push_back(jetPfCandidateTrack.at(iPfTrack));
                trueBJetLeptonTracksCharge.push_back(jetPfCandidateTrackCharge.at(iPfTrack));
                m_histogram["h_trueBJetLeptonTrackPt"]->Fill(trueBJetPfTrackPt, weight);
                m_histogram["h_trueBJetLeptonTrackEta"]->Fill(jetPfCandidateTrack.at(iPfTrack).Eta(), weight);
                m_histogram["h_trueBJetLeptonTrackCharge"]->Fill(jetPfCandidateTrackCharge.at(iPfTrack), weight);
                ((TH2D*)m_histogram["h_trueBJetLeptonChargePt"])->Fill(trueBJetPfTrackPt,jetPfCandidateTrackCharge.at(iPfTrack), weight);
                m_histogram["h_trueBJetToLeptonTrackPtRatio"]->Fill(trueBJetPfTrackPt/trueBJetPt, weight);
                if (ipValueWasDefined) m_histogram["h_trueBJetTrackIPValueIfLepton"]->Fill(impactParameterValueForPf, weight);
                if (ipValueWasDefined) m_histogram["h_trueBJetTrackIPSignificanceIfLepton"]->Fill(impactParameterSignificanceForPf, weight);
            }
                
            // Set some boolean to later on identify if the b-quark has the same charge as the lepton with highest p_{T} in the jet.
            if (trueBJetLeptonTracksCharge.size()>0 && trueBJetLeptonTracksCharge.at(0) <0 ) leptonHasNegativeCharge=true;
            if (trueBJetMuonTracksCharge.size()>0 && trueBJetMuonTracksCharge.at(0)<0) muonHasNegativeCharge = true;
            if (trueBJetElectronTracksCharge.size()>0 && trueBJetElectronTracksCharge.at(0)<0) electronHasNegativeCharge = true;
            
            double ptRatio = trueBJetPfTrackPt/trueBJetPt;
            ptRatioValues.push_back(ptRatio);
            
            if (ptRatioValues.size()==1)
            {
                leadingTrackPt = jetPfCandidateTrack.at(iPfTrack).pt();
                leadingTrackCharge = jetPfCandidateTrackCharge.at(iPfTrack);
                leadingTrackPx = jetPfCandidateTrack.at(iPfTrack).px();
                leadingTrackPy = jetPfCandidateTrack.at(iPfTrack).py();
                leadingTrackPz = jetPfCandidateTrack.at(iPfTrack).pz();
            }
            
            if (ptRatioValues.size()==1&&particleId==2) 
            {
                ptRatioValuesElectron = ptRatio;
                ptRatioValuesLepton = ptRatio;
                isLeadingElectron = true;
                isLeadingLepton = true;
                //thereIsALeadingLepton = true;
            }
            
            if (ptRatioValues.size()==1&&particleId==3) 
            {
                ptRatioValuesMuon = ptRatio;
                ptRatioValuesLepton = ptRatio;
                isLeadingMuon = true;
                isLeadingLepton = true;
                //thereIsALeadingMuon = true;
                //thereIsALeadingLepton = true;
            }
            
            if (ptRatioValues.size()==2&&particleId==2) 
            {
                ptRatioValuesSubleadingElectron = ptRatio;
                ptRatioValuesSubleadingLepton = ptRatio;
                subleadingTrackPt = jetPfCandidateTrack.at(iPfTrack).pt();
                subleadingTrackCharge = jetPfCandidateTrackCharge.at(iPfTrack);
                isSubleadingElectron = true;
                isSubleadingLepton = true;
            }
            
            if (ptRatioValues.size()==2&&particleId==3) 
            {
                ptRatioValuesSubleadingMuon = ptRatioValues.at(1);
                ptRatioValuesSubleadingLepton = ptRatioValues.at(1);
                subleadingTrackPt = jetPfCandidateTrack.at(iPfTrack).pt();
                subleadingTrackCharge = jetPfCandidateTrackCharge.at(iPfTrack);
                isSubleadingMuon = true;
                isSubleadingLepton = true;
            }
            
            if (ptRatioValues.size()==3) 
            {
                thirdleadingTrackPt = jetPfCandidateTrack.at(iPfTrack).pt();
                thirdleadingTrackCharge = jetPfCandidateTrackCharge.at(iPfTrack);
                if (particleId==2) isSubleadingElectron = true;
                if (particleId==3) isSubleadingMuon = true;
                if (particleId==2 || particleId==3) isSubleadingLepton = true;
            }
            
            if (ptRatioValues.size()>1&&particleId==3&&isLeadingMuon==false)  isNonLeadingLepton = true;
            
            m_histogram["h_trueBJetToTrackPtRatio"]->Fill(ptRatio, weight);
            m_histogram["h_trueBJetLeptonTracks"]->Fill(particleId, weight);
            m_histogram["h_trueBJetPfTrackPt"]->Fill(trueBJetPfTrackPt, weight);
            
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
            
            m_histogram["h_test_jetPfCandidateCharge"]->Fill(jetPfCandidateTrackCharge.at(iPfTrack), weight); 
            
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
            
            if (trueProduct>maxTrueProduct) maxTrueProduct = trueProduct;
            if (trueMagnitude>maxTrueMagnitude) maxTrueMagnitude = trueMagnitude;
            
            m_histogram["h_trueBJetRelPtTrack"]->Fill(trueMagnitude, weight);
            
            
        } //end of pfCandidates track loop
        
        // Check the particles leading the tracks
        if(trackParticleId.size()>=1) m_histogram["h_trueBJetLeadingTrackParticleId"] -> Fill(trackParticleId.at(0));
        if(trackParticleId.size()>=2) m_histogram["h_trueBJetSubLeadingTrackParticleId"] -> Fill(trackParticleId.at(1));
        if(trackParticleId.size()>=3) m_histogram["h_trueBJetThirdLeadingTrackParticleId"] -> Fill(trackParticleId.at(2));
        if(trackParticleId.size()>=4) m_histogram["h_trueBJetFourthLeadingTrackParticleId"] -> Fill(trackParticleId.at(3));
        
        if (trueBJetLeptonTracksCharge.size()>=1)
        {
            if (trueBJetLeptonTracksCharge.at(0) <0 && trackParticleId.at(0)==2) 
            {
                jetChargeHyerarchicalValue = -1;
                jetChargeHyerarchicalValueWasCalculated = true;
                m_histogram["h_trueBJetLeadingElectronTrackCharge"]->Fill(jetChargeHyerarchicalValue, weight);
                if (optionForCalibration==1 || optionForCalibration==2)
                {
                    if (jetChargeHyerarchicalValueWasCalculated==true) 
                    {
                        m_histogram["h_jetChargeHyerarchicalValue"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueElectrons"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValue_preStep"]->Fill(jetChargeHyerarchicalValue, weight);
                        if(secondaryVertexMultiplicityPerJet==1) ((TH2D*)m_histogram["h_trueBJetLeptonChargeVsSecondaryVertex"])->Fill(-1,nonWeightedSvSelTrackChargeVector.at(0), weight);
                    }
                }
                else if (optionForCalibration == 0)
                {
                    if ( jetHadronFlavour == 6 &&jetChargeHyerarchicalValueWasCalculated==true) 
                    {
                        m_histogram["h_jetChargeHyerarchicalValueBJets"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueBJetsElectrons"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueBJets_preStep"]->Fill(jetChargeHyerarchicalValue, weight);
                        if(secondaryVertexMultiplicityPerJet==1) ((TH2D*)m_histogram["h_trueBJetLeptonChargeVsSecondaryVertex"])->Fill(-1,nonWeightedSvSelTrackChargeVector.at(0), weight);
                    }
                    if ( jetHadronFlavour == -6 &&jetChargeHyerarchicalValueWasCalculated==true) 
                    {
                        m_histogram["h_jetChargeHyerarchicalValueAntiBJets"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueAntiBJetsElectrons"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueAntiBJets_preStep"]->Fill(jetChargeHyerarchicalValue, weight);
                        if(secondaryVertexMultiplicityPerJet==1) ((TH2D*)m_histogram["h_trueBJetLeptonChargeVsSecondaryVertex"])->Fill(-1,nonWeightedSvSelTrackChargeVector.at(0), weight);
                    }
                }
            }
            
            else if (trueBJetLeptonTracksCharge.at(0) >0 && trackParticleId.at(0)==2) 
            {
                jetChargeHyerarchicalValue = 1;
                jetChargeHyerarchicalValueWasCalculated = true;
                m_histogram["h_trueBJetLeadingElectronTrackCharge"]->Fill(jetChargeHyerarchicalValue, weight);
                if (optionForCalibration==1 || optionForCalibration==2)
                {
                    if (jetChargeHyerarchicalValueWasCalculated==true) 
                    {
                        m_histogram["h_jetChargeHyerarchicalValue"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueElectrons"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValue_preStep"]->Fill(jetChargeHyerarchicalValue, weight);
                        if(secondaryVertexMultiplicityPerJet==1) ((TH2D*)m_histogram["h_trueBJetLeptonChargeVsSecondaryVertex"])->Fill(1,nonWeightedSvSelTrackChargeVector.at(0), weight);
                    }
                }
                else if (optionForCalibration == 0)
                {
                    if ( jetHadronFlavour == 6 && jetChargeHyerarchicalValueWasCalculated==true) 
                    {
                        m_histogram["h_jetChargeHyerarchicalValueBJets"]->Fill(jetChargeHyerarchicalValue, weight); 
                        m_histogram["h_jetChargeHyerarchicalValueBJetsElectrons"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueBJets_preStep"]->Fill(jetChargeHyerarchicalValue, weight);
                        if(secondaryVertexMultiplicityPerJet==1) ((TH2D*)m_histogram["h_trueBJetLeptonChargeVsSecondaryVertex"])->Fill(1,nonWeightedSvSelTrackChargeVector.at(0), weight);
                    }
                    if ( jetHadronFlavour == -6 &&jetChargeHyerarchicalValueWasCalculated==true) 
                    {
                        m_histogram["h_jetChargeHyerarchicalValueAntiBJets"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueAntiBJetsElectrons"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueAntiBJets_preStep"]->Fill(jetChargeHyerarchicalValue, weight);
                        if(secondaryVertexMultiplicityPerJet==1) ((TH2D*)m_histogram["h_trueBJetLeptonChargeVsSecondaryVertex"])->Fill(1,nonWeightedSvSelTrackChargeVector.at(0), weight);
                    }
                }
            }
            
            else if (trueBJetLeptonTracksCharge.at(0) <0 && trackParticleId.at(0)==3) 
            {
                jetChargeHyerarchicalValue = - 1;
                jetChargeHyerarchicalValueWasCalculated = true;
                m_histogram["h_trueBJetLeadingMuonTrackCharge"]->Fill(jetChargeHyerarchicalValue, weight);
                if (optionForCalibration==1 || optionForCalibration==2)
                {
                    if (jetChargeHyerarchicalValueWasCalculated==true) 
                    {
                        m_histogram["h_jetChargeHyerarchicalValue"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueMuons"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValue_preStep"]->Fill(jetChargeHyerarchicalValue, weight);
                        if(secondaryVertexMultiplicityPerJet==1) ((TH2D*)m_histogram["h_trueBJetLeptonChargeVsSecondaryVertex"])->Fill(-1,nonWeightedSvSelTrackChargeVector.at(0), weight);
                    }
                }
                else if (optionForCalibration == 0)
                {
                    if ( jetHadronFlavour == 6 &&jetChargeHyerarchicalValueWasCalculated==true) 
                    {
                        m_histogram["h_jetChargeHyerarchicalValueBJets"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueBJetsMuons"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueBJets_preStep"]->Fill(jetChargeHyerarchicalValue, weight);
                        if(secondaryVertexMultiplicityPerJet==1) ((TH2D*)m_histogram["h_trueBJetLeptonChargeVsSecondaryVertex"])->Fill(-1,nonWeightedSvSelTrackChargeVector.at(0), weight);
                    }
                    if ( jetHadronFlavour == -6 &&jetChargeHyerarchicalValueWasCalculated==true) 
                    {
                        m_histogram["h_jetChargeHyerarchicalValueAntiBJets"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueAntiBJetsMuons"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueAntiBJets_preStep"]->Fill(jetChargeHyerarchicalValue, weight);
                        if(secondaryVertexMultiplicityPerJet==1) ((TH2D*)m_histogram["h_trueBJetLeptonChargeVsSecondaryVertex"])->Fill(-1,nonWeightedSvSelTrackChargeVector.at(0), weight);
                    }
                }
            }
            
            else if (trueBJetLeptonTracksCharge.at(0) >0 && trackParticleId.at(0)==3) 
            {
                jetChargeHyerarchicalValue = 1;
                jetChargeHyerarchicalValueWasCalculated = true;
                m_histogram["h_trueBJetLeadingMuonTrackCharge"]->Fill(jetChargeHyerarchicalValue, weight);
                if (optionForCalibration==1 || optionForCalibration==2)
                {
                    if (jetChargeHyerarchicalValueWasCalculated==true) 
                    {
                        m_histogram["h_jetChargeHyerarchicalValue"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueMuons"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValue_preStep"]->Fill(jetChargeHyerarchicalValue, weight);
                        if(secondaryVertexMultiplicityPerJet==1) ((TH2D*)m_histogram["h_trueBJetLeptonChargeVsSecondaryVertex"])->Fill(1,nonWeightedSvSelTrackChargeVector.at(0), weight);
                    }
                }
                else if (optionForCalibration == 0)
                {
                    if ( jetHadronFlavour == 6 && jetChargeHyerarchicalValueWasCalculated==true) 
                    {
                        m_histogram["h_jetChargeHyerarchicalValueBJets"]->Fill(jetChargeHyerarchicalValue, weight); 
                        m_histogram["h_jetChargeHyerarchicalValueBJetsMuons"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueBJets_preStep"]->Fill(jetChargeHyerarchicalValue, weight);
                        if(secondaryVertexMultiplicityPerJet==1) ((TH2D*)m_histogram["h_trueBJetLeptonChargeVsSecondaryVertex"])->Fill(1,nonWeightedSvSelTrackChargeVector.at(0), weight);
                    }
                    if ( jetHadronFlavour == -6 &&jetChargeHyerarchicalValueWasCalculated==true)
                    {
                        m_histogram["h_jetChargeHyerarchicalValueAntiBJets"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueAntiBJetsMuons"]->Fill(jetChargeHyerarchicalValue, weight);
                        m_histogram["h_jetChargeHyerarchicalValueAntiBJets_preStep"]->Fill(jetChargeHyerarchicalValue, weight);
                        if(secondaryVertexMultiplicityPerJet==1) ((TH2D*)m_histogram["h_trueBJetLeptonChargeVsSecondaryVertex"])->Fill(1,nonWeightedSvSelTrackChargeVector.at(0), weight);
                    }
                }
            }
        }
        
        // Location of muons in the track rank
        if (trackParticleId.size()>0) 
        {
            for (size_t iLep = 0;iLep!=trackParticleId.size();iLep++)
            {
                if (trackParticleId.at(iLep)==3) 
                {
                    m_histogram["h_trueBJetMuonRankInTrack"]->Fill(iLep);
                }
            }
        }
        
        if (ptRatioValues.size()!=0) m_histogram["h_trueBJetToLeadingTrackPtRatio"]->Fill(ptRatioValues.at(0), weight);
        if (isLeadingLepton==true) m_histogram["h_trueBJetToLeadingLeptonTrackPtRatio"]->Fill(ptRatioValuesLepton, weight);
        if (isLeadingMuon==true) m_histogram["h_trueBJetToLeadingMuonTrackPtRatio"]->Fill(ptRatioValuesMuon, weight);
        if (isLeadingElectron==true) m_histogram["h_trueBJetToLeadingElectronTrackPtRatio"]->Fill(ptRatioValuesElectron, weight);
        if (isSubleadingLepton==true) m_histogram["h_trueBJetToSubleadingLeptonTrackPtRatio"]->Fill(ptRatioValuesSubleadingLepton, weight);
        if (isSubleadingMuon==true) m_histogram["h_trueBJetToSubleadingMuonTrackPtRatio"]->Fill(ptRatioValuesSubleadingMuon, weight);
        if (isSubleadingElectron==true) m_histogram["h_trueBJetToSubleadingElectronTrackPtRatio"]->Fill(ptRatioValuesSubleadingElectron, weight);
        
        // Lepton histograms
        if (trueBJetLeptonTracksPt.size()>0) 
        {
            m_histogram["h_trueBJetLeptonTrackMultiplicity"] -> Fill(trueBJetLeptonTracksPt.size(), weight);
            ((TH2D*) m_histogram["h_trueBJetLeptonTrackPtMultiplicity"])->Fill(trueBJetPt,trueBJetTrackMultiplicity, weight);
            m_histogram["h_trueBJetTrackMultiplicityIfLepton"]->Fill(trueBJetTrackMultiplicity, weight);
        }
        
        // Hierarchical jet charge studies
        if (secondaryVertexMultiplicityPerJet==1 && multiplicityOfSecondaryVerticesForSelectedTracks.at(0)>=2)  
        {
            secondaryVertexChargeWasCalculated =true;
            
            if (optionForCalibration==0)
            {
                if (jetChargeHyerarchicalValueWasCalculated==false&&secondaryVertexChargeWasCalculated==true&&jetHadronFlavour == 6) m_histogram["h_jetChargeHyerarchicalValueBJets"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(0), weight);
                if (jetChargeHyerarchicalValueWasCalculated==false&&secondaryVertexChargeWasCalculated==true&&jetHadronFlavour == -6) m_histogram["h_jetChargeHyerarchicalValueAntiBJets"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(0), weight);
            }
            
            if (optionForCalibration==1||optionForCalibration==2)
            {
                if (jetChargeHyerarchicalValueWasCalculated==false&&secondaryVertexChargeWasCalculated==true) m_histogram["h_jetChargeHyerarchicalValue"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(0), weight);
            }
        }
        
        // Track general histograms
        m_histogram["h_trueBJetMaxRelPtTrack"]->Fill(maxTrueMagnitude, weight);
        m_histogram["h_trueBJetHighestPtTrack"]->Fill(maxPtTrueTrack, weight);
        
        ((TH2D*)m_histogram["h_trueBJetPtVsMaxPtTrack"])->Fill(maxPtTrueTrack,trueBJetPt, weight);
        ((TH2D*)m_histogram["h_trueBJetPtVsMaxRelPtTrack"])->Fill(maxTrueMagnitude,trueBJetPt, weight);
        ((TH2D*)m_histogram["h_trueBJetPtVsMaxScalarPtTrack"])->Fill(maxTrueProduct,trueBJetPt, weight);
        ((TH2D*)m_histogram["h_trueBJetPtVsNumTracks"])-> Fill(trueBJetTrackMultiplicity,trueBJetPt, weight);
        ((TH2D*)m_histogram["h_trueBJetMaxPtTrackVsNumTracks"])-> Fill(trueBJetTrackMultiplicity,maxPtTrueTrack, weight);
        ((TH2D*)m_histogram["h_trueBJetEtaVsTrackMultiplicity"])-> Fill(trueBJetTrackMultiplicity,trueBJetEta, weight);
        m_histogram["h_trueBJetPt"]->Fill(trueBJetPt, weight);
        
        // Track multiplicity plots: eta and pt slices
        m_histogram["h_trueBJetTrackMultiplicity"]->Fill(trueBJetTrackMultiplicity, weight);
        
        if (std::abs(trueBJetEta)<=0.2) m_histogram["h_trueBJetTrackMultiplicity_eta_0_02"]->Fill(trueBJetTrackMultiplicity, weight);
        else if (std::abs(trueBJetEta)>0.2 && std::abs(trueBJetEta)<=0.8) m_histogram["h_trueBJetTrackMultiplicity_eta_02_08"]->Fill(trueBJetTrackMultiplicity, weight);
        else if (std::abs(trueBJetEta)>0.8 && std::abs(trueBJetEta)<=1.0) m_histogram["h_trueBJetTrackMultiplicity_eta_08_10"]->Fill(trueBJetTrackMultiplicity, weight);
        else if (std::abs(trueBJetEta)>1.0 && std::abs(trueBJetEta)<=1.6) m_histogram["h_trueBJetTrackMultiplicity_eta_10_16"]->Fill(trueBJetTrackMultiplicity, weight);
        else if (std::abs(trueBJetEta)>1.6 && std::abs(trueBJetEta)<=2.0) m_histogram["h_trueBJetTrackMultiplicity_eta_16_20"]->Fill(trueBJetTrackMultiplicity, weight);
        else m_histogram["h_trueBJetTrackMultiplicity_eta_20"]->Fill(trueBJetTrackMultiplicity, weight);
        
        if (trueBJetPt>20.&&trueBJetPt<=30.) m_histogram["h_trueBJetTrackMultiplicity_pt_20_30"]->Fill(trueBJetTrackMultiplicity, weight);
        else if (trueBJetPt>30.&&trueBJetPt<=40.) m_histogram["h_trueBJetTrackMultiplicity_pt_30_40"]->Fill(trueBJetTrackMultiplicity, weight);
        else if (trueBJetPt>40.&&trueBJetPt<=70.) m_histogram["h_trueBJetTrackMultiplicity_pt_40_70"]->Fill(trueBJetTrackMultiplicity, weight);
        else if (trueBJetPt>70.&&trueBJetPt<=100.) m_histogram["h_trueBJetTrackMultiplicity_pt_70_100"]->Fill(trueBJetTrackMultiplicity, weight);
        else if (trueBJetPt>100.&&trueBJetPt<=150.) m_histogram["h_trueBJetTrackMultiplicity_pt_100_150"]->Fill(trueBJetTrackMultiplicity, weight);
        else m_histogram["h_trueBJetTrackMultiplicity_pt_150"]->Fill(trueBJetTrackMultiplicity, weight);
        
       
        // Validation charge histograms
        m_histogram["h_trueBJetScalarChargeValidation"]->Fill(trueBJetScalarCharge, weight);
        
        const double trueBJetScalarCharge10(sumTrueBMomentum.at(4)>0 ? sumTrueBMomentumQ.at(4)/sumTrueBMomentum.at(4) : 0);
        ((TH2D*)m_histogram["h_trueBJetScalarChargeVsMultip"])->Fill(trueBJetTrackMultiplicity,trueBJetScalarCharge10, weight);
        
        //check if lepton track and charge of the jet coincide
        if (trueBJetLeptonTracksPt.size()>0) 
        {
            // For leptons
            if(leptonHasNegativeCharge&&isLeadingLepton&&trueBJetScalarCharge10<0)
            {
                m_histogram["h_trueBJetLeadingLeptonScalarChargeMatch"]->Fill(1);
                m_histogram["h_trueBJetLeadingLeptonRelChargeMatch"]->Fill(1);
            }
            else if(leptonHasNegativeCharge==false&&isLeadingLepton&&trueBJetScalarCharge10>0) 
            {
                m_histogram["h_trueBJetLeadingLeptonScalarChargeMatch"]->Fill(1);
                m_histogram["h_trueBJetLeadingLeptonRelChargeMatch"]->Fill(1);
            }
            else if (leptonHasNegativeCharge&&isLeadingLepton&& trueBJetScalarCharge10>0)
            {
                m_histogram["h_trueBJetLeadingLeptonScalarChargeMatch"]->Fill(0);
                m_histogram["h_trueBJetLeadingLeptonRelChargeMatch"]->Fill(0);
            }
            else if (leptonHasNegativeCharge==false&&isLeadingLepton&&trueBJetScalarCharge10<0) 
            {
                m_histogram["h_trueBJetLeadingLeptonScalarChargeMatch"]->Fill(0);
                m_histogram["h_trueBJetLeadingLeptonRelChargeMatch"]->Fill(0);
            }
            
            // For muons
            if(muonHasNegativeCharge&&isLeadingMuon&&trueBJetScalarCharge10<0) m_histogram["h_trueBJetLeadingMuonScalarChargeMatch"]->Fill(1);
            else if(muonHasNegativeCharge==false&&isLeadingMuon&&trueBJetScalarCharge10>0) m_histogram["h_trueBJetLeadingMuonScalarChargeMatch"]->Fill(1);
            else if (muonHasNegativeCharge&&isLeadingMuon&& trueBJetScalarCharge10>0) m_histogram["h_trueBJetLeadingMuonScalarChargeMatch"]->Fill(0);
            else if (muonHasNegativeCharge==false&&isLeadingMuon&&trueBJetScalarCharge10<0) m_histogram["h_trueBJetLeadingMuonScalarChargeMatch"]->Fill(0);
            
            // For electrons
            if(electronHasNegativeCharge&&isLeadingElectron&&trueBJetScalarCharge10<0) m_histogram["h_trueBJetLeadingElectronScalarChargeMatch"]->Fill(1);
            else if(electronHasNegativeCharge==false&&isLeadingElectron&&trueBJetScalarCharge10>0) m_histogram["h_trueBJetLeadingElectronScalarChargeMatch"]->Fill(1);
            else if (electronHasNegativeCharge&&isLeadingElectron&& trueBJetScalarCharge10>0) m_histogram["h_trueBJetLeadingElectronScalarChargeMatch"]->Fill(0);
            else if (electronHasNegativeCharge==false&&isLeadingElectron&&trueBJetScalarCharge10<0) m_histogram["h_trueBJetLeadingElectronScalarChargeMatch"]->Fill(0);
            
            // For non-leading
            if(leptonHasNegativeCharge&&isNonLeadingLepton&&trueBJetScalarCharge10<0)
            {
                m_histogram["h_trueBJetNonLeadingLeptonScalarChargeMatch"]->Fill(1);
                m_histogram["h_trueBJetNonLeadingLeptonRelChargeMatch"]->Fill(1);
            }
            else if(leptonHasNegativeCharge==false&&isNonLeadingLepton&&trueBJetScalarCharge10>0) 
            {
                m_histogram["h_trueBJetNonLeadingLeptonScalarChargeMatch"]->Fill(1);
                m_histogram["h_trueBJetNonLeadingLeptonRelChargeMatch"]->Fill(1);
            }
            else if (leptonHasNegativeCharge &&isNonLeadingLepton&& trueBJetScalarCharge10>0)
            {
                m_histogram["h_trueBJetNonLeadingLeptonScalarChargeMatch"]->Fill(0);
                m_histogram["h_trueBJetNonLeadingLeptonRelChargeMatch"]->Fill(0);
            }
            else if (leptonHasNegativeCharge==false&&isNonLeadingLepton&&trueBJetScalarCharge10<0) 
            {
                m_histogram["h_trueBJetNonLeadingLeptonScalarChargeMatch"]->Fill(0);
                m_histogram["h_trueBJetNonLeadingLeptonRelChargeMatch"]->Fill(0);
            }
        }
        
        // Create a vector of histograms: contains histograms for x=0.2 to x=2.0
        std::vector<TString> trueBJetHistogramScalarVector;
        std::vector<TString> trueBJetHistogramRelVector;
       
        for (double iHisto=2;iHisto<=20;iHisto+=2)
        {
            std::stringstream ss_histoScalar;
            std::stringstream ss_histoRel;
            
            ss_histoScalar<<"h_trueBJetScalarCharge"<<iHisto;
            ss_histoRel<<"h_trueBJetRelCharge"<<iHisto;
           
            TString histoScalar = ss_histoScalar.str();
            TString histoRel = ss_histoRel.str();
           
            trueBJetHistogramScalarVector.push_back(histoScalar);
            trueBJetHistogramRelVector.push_back(histoRel);
          
        }
        
        // Fill the histograms with the corresponding value of the charge (we take values from the sumTrueBMagnitude and sumTrueBMomentum vectors)
        std::vector<double> trueBJetScalarChargeVector;
        std::vector<double> trueBJetRelChargeVector;
        
        for (size_t iSumHis=0;iSumHis!=sumTrueBMomentum.size();iSumHis++)
        {
            const double trueBJetScalarChargeForX(sumTrueBMomentum.at(iSumHis)>0 ? sumTrueBMomentumQ.at(iSumHis)/sumTrueBMomentum.at(iSumHis) : 0);
            const double trueBJetRelCharge(sumTrueBMagnitude.at(iSumHis)>0 ? sumTrueBMagnitudeQ.at(iSumHis)/sumTrueBMagnitude.at(iSumHis) : 0);
            TString histoScalar = trueBJetHistogramScalarVector.at(iSumHis);
            TString histoRel = trueBJetHistogramRelVector.at(iSumHis);
            
            m_histogram[histoScalar]->Fill(trueBJetScalarChargeForX, weight);
            m_histogram[histoRel]->Fill(trueBJetRelCharge, weight);
            
            trueBJetScalarChargeVector.push_back(trueBJetScalarChargeForX);
            trueBJetRelChargeVector.push_back(trueBJetRelCharge);
        }
        
        // TEST new function added for jet charge calculation, for a given x "squeezing"-parameter
        double jetChargeFromFunction = ptWeightedJetChargeX(jetIdx, jets, 0.8, jetPfCandidateTrackIndex, jetPfCandidateTrack, jetPfCandidateTrackCharge, jetPfCandidatePrimaryVertexId);
        m_histogram["h_jetChargeFromFunction"]->Fill(jetChargeFromFunction, weight);
        
        // Only fill charge if there's ONE secondary vertex in the event
        if(secondaryVertexMultiplicityPerJet==1) m_histogram["h_trueBJetScalarCharge8_IfOneSvInJet"]->Fill(trueBJetScalarChargeVector.at(3), weight);
        if(sumTrueBMomentum.at(3)>0) m_histogram["h_trueBJetScalarCharge8_WithoutZeroValue"]->Fill(trueBJetScalarChargeVector.at(3), weight);
        
        if(secondaryVertexMultiplicityPerJet==1) ((TH2D*)m_histogram["h_trueBJetScalarChargeVsSecondaryVertex"])->Fill(trueBJetScalarChargeVector.at(3),chargeOfSecondaryVerticesForSelectedTracks.at(0), weight);
        
        
        if (optionForCalibration==1 || optionForCalibration==2)
        {
            if (jetChargeHyerarchicalValueWasCalculated==false && secondaryVertexChargeWasCalculated==false) 
            {
                m_histogram["h_jetChargeHyerarchicalValue"]->Fill(trueBJetScalarChargeVector.at(3), weight);
                m_histogram["h_jetChargeHyerarchicalValueElectrons"]->Fill(trueBJetScalarChargeVector.at(3), weight);
                m_histogram["h_jetChargeHyerarchicalValueMuons"]->Fill(trueBJetScalarChargeVector.at(3), weight);
                
                m_histogram["h_jetChargeHyerarchicalValueTruncated"]->Fill(trueBJetScalarChargeVector.at(3), weight);
                m_histogram["h_jetChargeHyerarchicalValueElectronsTruncated"]->Fill(trueBJetScalarChargeVector.at(3), weight);
                m_histogram["h_jetChargeHyerarchicalValueMuonsTruncated"]->Fill(trueBJetScalarChargeVector.at(3), weight);
            }
            
            if (jetChargeHyerarchicalValueWasCalculated==false) m_histogram["h_jetChargeHyerarchicalValue_preStep"]->Fill(trueBJetScalarChargeVector.at(3), weight);
        }
        else if (optionForCalibration==0)
        {
            if (jetChargeHyerarchicalValueWasCalculated==false && secondaryVertexChargeWasCalculated==false &&  jetHadronFlavour == 6) 
            {
                m_histogram["h_jetChargeHyerarchicalValueBJets"]->Fill(trueBJetScalarChargeVector.at(3), weight);
                m_histogram["h_jetChargeHyerarchicalValueBJetsElectrons"]->Fill(trueBJetScalarChargeVector.at(3), weight);
                m_histogram["h_jetChargeHyerarchicalValueBJetsMuons"]->Fill(trueBJetScalarChargeVector.at(3), weight);
                
                m_histogram["h_jetChargeHyerarchicalValueBJetsTruncated"]->Fill(trueBJetScalarChargeVector.at(3), weight);
                m_histogram["h_jetChargeHyerarchicalValueBJetsElectronsTruncated"]->Fill(trueBJetScalarChargeVector.at(3), weight);
                m_histogram["h_jetChargeHyerarchicalValueBJetsMuonsTruncated"]->Fill(trueBJetScalarChargeVector.at(3), weight);
            }
            
            if (jetChargeHyerarchicalValueWasCalculated==false && jetHadronFlavour == 6) m_histogram["h_jetChargeHyerarchicalValueBJets_preStep"]->Fill(trueBJetScalarChargeVector.at(3), weight); 
            
            if (jetChargeHyerarchicalValueWasCalculated==false && secondaryVertexChargeWasCalculated==false &&  jetHadronFlavour == -6) 
            {
                m_histogram["h_jetChargeHyerarchicalValueAntiBJets"]->Fill(trueBJetScalarChargeVector.at(3), weight);
                m_histogram["h_jetChargeHyerarchicalValueAntiBJetsElectrons"]->Fill(trueBJetScalarChargeVector.at(3), weight);
                m_histogram["h_jetChargeHyerarchicalValueAntiBJetsMuons"]->Fill(trueBJetScalarChargeVector.at(3), weight);
                
                m_histogram["h_jetChargeHyerarchicalValueAntiBJetsTruncated"]->Fill(trueBJetScalarChargeVector.at(3), weight);
                m_histogram["h_jetChargeHyerarchicalValueAntiBJetsElectronsTruncated"]->Fill(trueBJetScalarChargeVector.at(3), weight);
                m_histogram["h_jetChargeHyerarchicalValueAntiBJetsMuonsTruncated"]->Fill(trueBJetScalarChargeVector.at(3), weight);
            }
            
            if (jetChargeHyerarchicalValueWasCalculated==false && jetHadronFlavour == -6) m_histogram["h_jetChargeHyerarchicalValueAntiBJets_preStep"]->Fill(trueBJetScalarChargeVector.at(3), weight); 
        }
        
        
        if (optionForCalibration==0)
        {
            if (numHadMatched<1) continue;
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkScalarCharge2"]->Fill(trueBJetScalarChargeVector.at(0), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge2"]->Fill(trueBJetScalarChargeVector.at(0), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkScalarCharge4"]->Fill(trueBJetScalarChargeVector.at(1), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge4"]->Fill(trueBJetScalarChargeVector.at(1), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkScalarCharge6"]->Fill(trueBJetScalarChargeVector.at(2), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge6"]->Fill(trueBJetScalarChargeVector.at(2), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkScalarCharge8"]->Fill(trueBJetScalarChargeVector.at(3), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge8"]->Fill(trueBJetScalarChargeVector.at(3), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkScalarCharge10"]->Fill(trueBJetScalarChargeVector.at(4), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueBJetAntiBQuarkScalarCharge10"]->Fill(trueBJetScalarChargeVector.at(4), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkScalarCharge12"]->Fill(trueBJetScalarChargeVector.at(5), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge12"]->Fill(trueBJetScalarChargeVector.at(5), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkScalarCharge14"]->Fill(trueBJetScalarChargeVector.at(6), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge14"]->Fill(trueBJetScalarChargeVector.at(6), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkScalarCharge16"]->Fill(trueBJetScalarChargeVector.at(7), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge16"]->Fill(trueBJetScalarChargeVector.at(7), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkScalarCharge18"]->Fill(trueBJetScalarChargeVector.at(8), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge18"]->Fill(trueBJetScalarChargeVector.at(8), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkScalarCharge20"]->Fill(trueBJetScalarChargeVector.at(9), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge20"]->Fill(trueBJetScalarChargeVector.at(9), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkRelCharge2"]->Fill(trueBJetRelChargeVector.at(0), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge2"]->Fill(trueBJetRelChargeVector.at(0), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkRelCharge4"]->Fill(trueBJetRelChargeVector.at(1), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge4"]->Fill(trueBJetRelChargeVector.at(1), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkRelCharge6"]->Fill(trueBJetRelChargeVector.at(2), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge6"]->Fill(trueBJetRelChargeVector.at(2), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkRelCharge8"]->Fill(trueBJetRelChargeVector.at(3), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge8"]->Fill(trueBJetRelChargeVector.at(3), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkRelCharge10"]->Fill(trueBJetRelChargeVector.at(4), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueBJetAntiBQuarkRelCharge10"]->Fill(trueBJetRelChargeVector.at(4), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkRelCharge12"]->Fill(trueBJetRelChargeVector.at(5), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge12"]->Fill(trueBJetRelChargeVector.at(5), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkRelCharge14"]->Fill(trueBJetRelChargeVector.at(6), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge14"]->Fill(trueBJetRelChargeVector.at(6), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkRelCharge16"]->Fill(trueBJetRelChargeVector.at(7), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge16"]->Fill(trueBJetRelChargeVector.at(7), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkRelCharge18"]->Fill(trueBJetRelChargeVector.at(8), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge18"]->Fill(trueBJetRelChargeVector.at(8), weight);
            
            if ( jetHadronFlavour>0)  m_histogram["h_trueBJetBQuarkRelCharge20"]->Fill(trueBJetRelChargeVector.at(9), weight);
            else if ( jetHadronFlavour<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge20"]->Fill(trueBJetRelChargeVector.at(9), weight);
            
            // Only coming from top: for comparison with kinReco results
            if ( jetHadronFlavour==6)  m_histogram["h_trueBJetBQuarkFromTopScalarCharge2"]->Fill(trueBJetScalarChargeVector.at(0), weight);
            else if ( jetHadronFlavour==-6)  m_histogram["h_trueAntiBJetBQuarkFromTopScalarCharge2"]->Fill(trueBJetScalarChargeVector.at(0), weight);
            
            if ( jetHadronFlavour==6)  m_histogram["h_trueBJetBQuarkFromTopScalarCharge4"]->Fill(trueBJetScalarChargeVector.at(1), weight);
            else if ( jetHadronFlavour==-6)  m_histogram["h_trueAntiBJetBQuarkFromTopScalarCharge4"]->Fill(trueBJetScalarChargeVector.at(1), weight);
            
            if ( jetHadronFlavour==6)  m_histogram["h_trueBJetBQuarkFromTopScalarCharge6"]->Fill(trueBJetScalarChargeVector.at(2), weight);
            else if ( jetHadronFlavour==-6)  m_histogram["h_trueAntiBJetBQuarkFromTopScalarCharge6"]->Fill(trueBJetScalarChargeVector.at(2), weight);
            
            if ( jetHadronFlavour==6)  m_histogram["h_trueBJetBQuarkFromTopScalarCharge8"]->Fill(trueBJetScalarChargeVector.at(3), weight);
            else if ( jetHadronFlavour==-6)  m_histogram["h_trueAntiBJetBQuarkFromTopScalarCharge8"]->Fill(trueBJetScalarChargeVector.at(3), weight);
            
            if ( jetHadronFlavour==6)  m_histogram["h_trueBJetBQuarkFromTopScalarCharge10"]->Fill(trueBJetScalarChargeVector.at(4), weight);
            else if ( jetHadronFlavour==-6)  m_histogram["h_trueBJetAntiBQuarkFromTopScalarCharge10"]->Fill(trueBJetScalarChargeVector.at(4), weight);
            
            if ( jetHadronFlavour==6)  m_histogram["h_trueBJetBQuarkFromTopScalarCharge12"]->Fill(trueBJetScalarChargeVector.at(5), weight);
            else if ( jetHadronFlavour==-6)  m_histogram["h_trueAntiBJetBQuarkFromTopScalarCharge12"]->Fill(trueBJetScalarChargeVector.at(5), weight);
            
            if ( jetHadronFlavour==6)  m_histogram["h_trueBJetBQuarkFromTopScalarCharge14"]->Fill(trueBJetScalarChargeVector.at(6), weight);
            else if ( jetHadronFlavour==-6)  m_histogram["h_trueAntiBJetBQuarkFromTopScalarCharge14"]->Fill(trueBJetScalarChargeVector.at(6), weight);
            
            if ( jetHadronFlavour==6)  m_histogram["h_trueBJetBQuarkFromTopScalarCharge16"]->Fill(trueBJetScalarChargeVector.at(7), weight);
            else if ( jetHadronFlavour==-6)  m_histogram["h_trueAntiBJetBQuarkFromTopScalarCharge16"]->Fill(trueBJetScalarChargeVector.at(7), weight);
            
            if ( jetHadronFlavour==6)  m_histogram["h_trueBJetBQuarkFromTopScalarCharge18"]->Fill(trueBJetScalarChargeVector.at(8), weight);
            else if ( jetHadronFlavour==-6)  m_histogram["h_trueAntiBJetBQuarkFromTopScalarCharge18"]->Fill(trueBJetScalarChargeVector.at(8), weight);
            
            if ( jetHadronFlavour==6)  m_histogram["h_trueBJetBQuarkFromTopScalarCharge20"]->Fill(trueBJetScalarChargeVector.at(9), weight);
            else if ( jetHadronFlavour==-6)  m_histogram["h_trueAntiBJetBQuarkFromTopScalarCharge20"]->Fill(trueBJetScalarChargeVector.at(9), weight);
        }
        
        // Testing mva variable separation
        double leadingTrackPtWeightedCharge = leadingTrackCharge*leadingTrackPt/trueBJetPt;
        double subleadingTrackPtWeightedCharge = subleadingTrackCharge*subleadingTrackPt/trueBJetPt;
        double thirdleadingTrackPtWeightedCharge = thirdleadingTrackCharge*thirdleadingTrackPt/trueBJetPt;
        double leadingMuonPtWeightedCharge = 0;
        if (isLeadingMuon) leadingMuonPtWeightedCharge = trueBJetMuonTracksCharge.at(0)*trueBJetMuonTracksPt.at(0)/trueBJetPt;
        double leadingElectronPtWeightedCharge = 0;
        if (isLeadingElectron)leadingElectronPtWeightedCharge = trueBJetElectronTracksCharge.at(0)*trueBJetElectronTracksPt.at(0)/trueBJetPt;
        double leadingTrackPtWeightedCharge1 = leadingTrackCharge*(leadingTrackPx*jetTrueBPx + leadingTrackPy*jetTrueBPy + leadingTrackPz*jetTrueBPz);
        double leadingTrackPtWeightedCharge2 = leadingTrackCharge*(std::sqrt(std::pow(jetTrueBPy*leadingTrackPz-jetTrueBPz*leadingTrackPy,2) + std::pow(jetTrueBPx*leadingTrackPz-jetTrueBPz*leadingTrackPx,2) + std::pow(jetTrueBPy*leadingTrackPx-jetTrueBPx*leadingTrackPy,2)));
        double trackNumberWeightedJetPt = leadingTrackCharge*trueBJetPt/trueBJetTrackMultiplicity;
        double chargeWeightedTrackId = leadingTrackCharge;
        if (isLeadingMuon) chargeWeightedTrackId = leadingTrackCharge*3;
        if (isLeadingElectron) chargeWeightedTrackId = leadingTrackCharge*2;
        double svChargeWeightedFlightDistance = 0;
        if (thereIsASecondaryVertex) svChargeWeightedFlightDistance = chargeOfSecondaryVerticesForSelectedTracks.at(0) * jetSecondaryVertexFlightDistanceSignificance.at(0); 
        
        
        m_histogram["h_mva_leadingTrackPtWeightedCharge"]->Fill(leadingTrackPtWeightedCharge, weight);
        m_histogram["h_mva_subleadingTrackPtWeightedCharge"]->Fill(subleadingTrackPtWeightedCharge, weight);
        m_histogram["h_mva_thirdleadingTrackPtWeightedCharge"]->Fill(thirdleadingTrackPtWeightedCharge, weight);
        if (isLeadingMuon) m_histogram["h_mva_leadingMuonPtWeightedCharge"]->Fill(leadingMuonPtWeightedCharge, weight);
        if (isLeadingElectron) m_histogram["h_mva_leadingElectronPtWeightedCharge"]->Fill(leadingElectronPtWeightedCharge, weight);
        m_histogram["h_mva_leadingTrackPtWeightedChargeScalar"]->Fill(leadingTrackPtWeightedCharge1, weight);
        m_histogram["h_mva_leadingTrackPtWeightedChargeVectorial"]->Fill(leadingTrackPtWeightedCharge2, weight);
        m_histogram["h_mva_leadingMuonPtWeightedChargeForAllJets"]->Fill(leadingMuonPtWeightedCharge, weight);
        m_histogram["h_mva_trackNumberWeightedJetPt"]->Fill(trackNumberWeightedJetPt, weight);
        m_histogram["h_mva_chargeWeightedTrackId"]->Fill(chargeWeightedTrackId, weight);
        m_histogram["h_mva_svChargeWeightedFlightDistance"]->Fill(svChargeWeightedFlightDistance, weight);
        if (thereIsASecondaryVertex) m_histogram["h_mva_secondaryVertexWeightedCharge"]->Fill(chargeOfSecondaryVerticesForSelectedTracks.at(0), weight);
        else m_histogram["h_mva_secondaryVertexWeightedCharge"]->Fill(0., weight);
        m_histogram["h_mva_jetCharge_x08"]->Fill(trueBJetScalarChargeVector.at(3), weight);
        fillTree = true;
        
        // MVA specific variable filling
        if (jetHadronFlavour>0) mvaStruct_.trueBJetId_ = -1;
        else if (jetHadronFlavour<0) mvaStruct_.trueBJetId_ = 0;
        mvaStruct_.relChargeJet_ = trueBJetRelChargeVector.at(3);
        mvaStruct_.longChargeJet_ = trueBJetScalarChargeVector.at(3);
        mvaStruct_.leadingTrackPtWeightedCharge_ = leadingTrackPtWeightedCharge;
        mvaStruct_.subleadingTrackPtWeightedCharge_ = subleadingTrackPtWeightedCharge;
        mvaStruct_.thirdleadingTrackPtWeightedCharge_ = thirdleadingTrackPtWeightedCharge;
        mvaStruct_.leadingMuonPtWeightedCharge_ = leadingMuonPtWeightedCharge;
        mvaStruct_.leadingElectronPtWeightedCharge_ = leadingElectronPtWeightedCharge;
        mvaStruct_.trackNumberWeightedJetPt_ = trackNumberWeightedJetPt;
        mvaStruct_.chargeWeightedTrackId_ = chargeWeightedTrackId;
        mvaStruct_.svChargeWeightedFlightDistance_ = svChargeWeightedFlightDistance;
        if (thereIsASecondaryVertex) mvaStruct_.secondaryVertexCharge_ = chargeOfSecondaryVerticesForSelectedTracks.at(0);
        else mvaStruct_.secondaryVertexCharge_ = 0.; 
        if (impactParameterValuesForPf.size()!=0) mvaStruct_.ipSignificanceLeadingTrack_ = impactParameterSignificanceOfLeadingTrack;
        else (mvaStruct_.ipSignificanceLeadingTrack_ = 0.);
        
        //Fill the mva trees
        if (fillTree) 
        {
            if ((eventNumber&1)==0) mvaChargeTestTree_->Fill();
            else mvaChargeTrainTree_->Fill();
        }
        
        //// MVA reader variable
        //double val1 = reader_->EvaluateMVA ("testingTheReader");
        
        //if (recoBjetFromTopIndex==jetIdx) m_histogram["h_testB"]->Fill(val1);
        //if (recoAntiBjetFromTopIndex==jetIdx) m_histogram["h_testAntiB"]->Fill(val1);
        //m_histogram["h_weights"]->Fill(val1);
       
       //if (jetHadronFlavour>0) m_histogram["h_testB_fromAll"]->Fill(val1);
       //if (jetHadronFlavour<0) m_histogram["h_testAntiB_fromAll"]->Fill(val1);
       
       //if (val1<0) m_histogram["h_testB_fromAll_charge"]->Fill(trueBJetScalarChargeVector.at(3));
       //if (val1>0) m_histogram["h_testAntiB_fromAll_charge"]->Fill(trueBJetScalarChargeVector.at(3));
       
       //if (jetHadronFlavour>0) m_histogram["h_trueB_fromAll_charge"]->Fill(trueBJetScalarChargeVector.at(3));
       //if (jetHadronFlavour<0) m_histogram["h_trueAntiB_fromAll_charge"]->Fill(trueBJetScalarChargeVector.at(3));
       
       //if (jetHadronFlavour>0 && val1<0) m_histogram["h_coincidenceTest_flavourB"]->Fill(trueBJetScalarChargeVector.at(3));
       //else if (jetHadronFlavour>0&&val1>0) m_histogram["h_coincidenceTest_flavourB_fail"]->Fill(trueBJetScalarChargeVector.at(3));
       
       //if (jetHadronFlavour<0 && val1>0) m_histogram["h_coincidenceTest_flavourAntiB"]->Fill(trueBJetScalarChargeVector.at(3));
       //else if (jetHadronFlavour<0 && val1<0) m_histogram["h_coincidenceTest_flavourAntiB_fail"]->Fill(trueBJetScalarChargeVector.at(3));
       
       //if (recoBjetFromTopIndex==jetIdx && val1<0) m_histogram["h_coincidenceTest_topB"]->Fill(1);
       //else if (recoBjetFromTopIndex==jetIdx&&val1>0) m_histogram["h_coincidenceTest_topB"]->Fill(0);
       
       //if (recoAntiBjetFromTopIndex==jetIdx && val1>0) m_histogram["h_coincidenceTest_topAntiB"]->Fill(1);
       //else if (recoAntiBjetFromTopIndex==jetIdx && val1<0) m_histogram["h_coincidenceTest_topAntiB"]->Fill(0);
       
       //if (jetHadronFlavour>0 && trueBJetScalarCharge<0) m_histogram["h_coincidenceTest_flavourB_oldDefinition"]->Fill(1);
       //else if (jetHadronFlavour>0 && trueBJetScalarCharge>0) m_histogram["h_coincidenceTest_flavourB_oldDefinition"]->Fill(0);
       
       //if (jetHadronFlavour<0 && trueBJetScalarCharge>0) m_histogram["h_coincidenceTest_flavourAntiB_oldDefinition"]->Fill(1);
       //else if (jetHadronFlavour<0 && trueBJetScalarCharge<0) m_histogram["h_coincidenceTest_flavourAntiB_oldDefinition"]->Fill(0);
       
       //if (recoBjetFromTopIndex==jetIdx&& trueBJetScalarCharge<0) m_histogram["h_coincidenceTest_topB_oldDefinition"]->Fill(1);
       //else if (recoBjetFromTopIndex==jetIdx && trueBJetScalarCharge>0) m_histogram["h_coincidenceTest_topB_oldDefinition"]->Fill(0);
       
       //if (recoAntiBjetFromTopIndex==jetIdx&& trueBJetScalarCharge>0) m_histogram["h_coincidenceTest_topAntiB_oldDefinition"]->Fill(1);
       //else if (recoAntiBjetFromTopIndex==jetIdx && trueBJetScalarCharge<0) m_histogram["h_coincidenceTest_topAntiB_oldDefinition"]->Fill(0);
       
        
    } //end loop over reco jets
    
} //END OF JET CHARGE ANALYZER FUNCTION


//==========================================HISTOGRAM FILLING FUNCTIONS===============================================

void AnalyzerJetCharge::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    bookJetHistos(step, m_histogram);
    bookPfCandidateHistos(step, m_histogram);
    bookSelectedTrackHistos(step, m_histogram);
    bookMvaHistos(step, m_histogram);
    bookOtherHistos(step, m_histogram);
}

void AnalyzerJetCharge::bookMvaHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    mvaChargeTestTree_ = store (new TTree("mvaChargeTestTree_","mvaChargeTestTree_"));
    
    mvaChargeTestTree_->Branch("trueBJetId", &(mvaStruct_.trueBJetId_));
    mvaChargeTestTree_->Branch("longChargeJet", &(mvaStruct_.longChargeJet_));
    mvaChargeTestTree_->Branch("relChargeJet", &(mvaStruct_.relChargeJet_));
    mvaChargeTestTree_->Branch("leadingTrackPtWeightedCharge", &(mvaStruct_.leadingTrackPtWeightedCharge_));
    mvaChargeTestTree_->Branch("subleadingTrackPtWeightedCharge", &(mvaStruct_.subleadingTrackPtWeightedCharge_));
    mvaChargeTestTree_->Branch("thirdleadingTrackPtWeightedCharge", &(mvaStruct_.thirdleadingTrackPtWeightedCharge_));
    mvaChargeTestTree_->Branch("leadingMuonPtWeightedCharge", &(mvaStruct_.leadingMuonPtWeightedCharge_));
    mvaChargeTestTree_->Branch("leadingElectronPtWeightedCharge", &(mvaStruct_.leadingElectronPtWeightedCharge_));
    mvaChargeTestTree_->Branch("trackNumberWeightedJetPt", &(mvaStruct_.trackNumberWeightedJetPt_));
    mvaChargeTestTree_->Branch("chargeWeightedTrackId",&(mvaStruct_.chargeWeightedTrackId_));
    mvaChargeTestTree_->Branch("svChargeWeightedFlightDistance",&(mvaStruct_.svChargeWeightedFlightDistance_));
    mvaChargeTestTree_->Branch("secondaryVertexCharge",&(mvaStruct_.secondaryVertexCharge_));
    mvaChargeTestTree_->Branch("ipSignificanceLeadingTrack",&(mvaStruct_.ipSignificanceLeadingTrack_));
    
    mvaChargeTrainTree_ = store (new TTree("mvaChargeTrainTree_","mvaChargeTrainTree_"));
    
    mvaChargeTrainTree_->Branch("trueBJetId", &(mvaStruct_.trueBJetId_));
    mvaChargeTrainTree_->Branch("longChargeJet", &(mvaStruct_.longChargeJet_));
    mvaChargeTrainTree_->Branch("relChargeJet", &(mvaStruct_.relChargeJet_));
    mvaChargeTrainTree_->Branch("leadingTrackPtWeightedCharge", &(mvaStruct_.leadingTrackPtWeightedCharge_));
    mvaChargeTrainTree_->Branch("subleadingTrackPtWeightedCharge", &(mvaStruct_.subleadingTrackPtWeightedCharge_));
    mvaChargeTrainTree_->Branch("thirdleadingTrackPtWeightedCharge", &(mvaStruct_.thirdleadingTrackPtWeightedCharge_));
    mvaChargeTrainTree_->Branch("leadingMuonPtWeightedCharge", &(mvaStruct_.leadingMuonPtWeightedCharge_));
    mvaChargeTrainTree_->Branch("leadingElectronPtWeightedCharge", &(mvaStruct_.leadingElectronPtWeightedCharge_));
    mvaChargeTrainTree_->Branch("trackNumberWeightedJetPt", &(mvaStruct_.trackNumberWeightedJetPt_));
    mvaChargeTrainTree_->Branch("chargeWeightedTrackId",&(mvaStruct_.chargeWeightedTrackId_));
    mvaChargeTrainTree_->Branch("svChargeWeightedFlightDistance",&(mvaStruct_.svChargeWeightedFlightDistance_));
    mvaChargeTrainTree_->Branch("secondaryVertexCharge",&(mvaStruct_.secondaryVertexCharge_));
    mvaChargeTrainTree_->Branch("ipSignificanceLeadingTrack",&(mvaStruct_.ipSignificanceLeadingTrack_));
    
  
    // MVA variable separation testing
    name = "h_mva_leadingTrackPtWeightedCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of the leading track weighted by the p_{T} ratio with respect to the jet p_{T};c^{track}*p_{T}^{track}*p_{T}^{jet};Jets",80,-1.,1.));
    
    name = "h_mva_leadingTrackPtWeightedChargeScalar";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of the leading track weighted by the scalar product between the jet and the track;c^{track}*p^{track}*p^{jet};Jets",80,-5000.,5000.));
    
    name = "h_mva_leadingTrackPtWeightedChargeVectorial";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of the leading track weighted by the p_{T} ratio with respect to the jet p_{T};c^{track}*p_{T}^{track}*p_{T}^{jet};Jets",80,-1000.,1000.));
    
    name = "h_mva_subleadingTrackPtWeightedCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of the subleading track weighted by the p_{T} ratio with respect to the jet p_{T};c^{track}*p_{T}^{track}/p_{T}^{jet};Jets",80,-1.,1.));
    
    name = "h_mva_thirdleadingTrackPtWeightedCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of the thirdleading track weighted by the p_{T} ratio with respect to the jet p_{T};c^{track}*p_{T}^{track}/p_{T}^{jet};Jets",80,-1.,1.));
    
    name = "h_mva_leadingMuonPtWeightedCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of the leading muon weighted by the p_{T} ratio with respect to the jet p_{T};c^{track}*p_{T}^{track}/p_{T}^{jet} for leading muon;Jets",80,-1.,1.));
    
    name = "h_mva_leadingElectronPtWeightedCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of the leading electron weighted by the p_{T} ratio with respect to the jet p_{T};c^{track}*p_{T}^{track}/p_{T}^{jet} for leading electron;Jets",80,-1.,1.));
    
    name = "h_mva_leadingMuonPtWeightedChargeForAllJets";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of the leading muon weighted by the p_{T} ratio with respect to the jet p_{T} - all jets;c^{track}*p_{T}^{track}/p_{T}^{jet};Jets",80,-1.,1.));
    
    name = "h_mva_trackNumberWeightedJetPt";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Jet p_{T} divided by the track multiplicity;c^{track}*p_{T}^{jet}/(number tracks);Jets",80,-40.,40.));
    
    name = "h_mva_chargeWeightedTrackId";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Id of the leading track multiplied by the charge of it;c^{track}*Id_{track};Jets",8,-4.,4.));
    
    name = "h_mva_secondaryVertexWeightedCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex weighted charge with respect to the jet axis;c_{rel}^{SV};Jets",24,-1.2,1.2));
    
    name = "h_mva_secondaryVertexWeightedCharge_ifSV";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex weighted charge with respect to the jet axis (if there's a SV);c_{rel}^{SV};Jets",24,-1.2,1.2));
    
    name = "h_mva_svChargeWeightedFlightDistance";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex flight distance times secondary vertex weighted charge;c_{rel}^{SV}*(flight distance significance);Jets",90,-50.,50.));
    
    name = "h_mva_jetCharge_x08";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal p_{T} weighted charge with x = 1.0;c_{rel}^{jet};Jets",24,-1.2,1.2));
    
    name = "h_testB";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Weight of b jets from system;weight;Jets",100,-1.,1.));
    
    name = "h_testAntiB";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Weight of anti-b jets from top system;weight;Jets",100,-1.,1.));
    
    name = "h_weights";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Weights coming out from the mva;weight;Jets",100,-1.,1.));
    
    name = "h_testB_fromAll";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Weight of b jets from all;weight;Jets",100,-1.,1.));
    
    name = "h_testAntiB_fromAll";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Weight of anti-b jets from all;weight;Jets",100,-1.,1.));
    
    name = "h_trueB_fromAll_charge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of true b jets from all;c^{jet}_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiB_fromAll_charge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of anti-b jets from all;c^{jet}_{rel};Jets",24,-1.2,1.2));
    
    name = "h_testB_fromAll_charge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of b jets from all if mva weight<0;c^{jet}_{rel};Jets",24,-1.2,1.2));
    
    name = "h_testAntiB_fromAll_charge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of anti-b jets from all if mva weight>0;c^{jet}_{rel};Jets",24,-1.2,1.2));
    
    name = "h_coincidenceTest_flavourB";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge for B hadron flavour when truth and weight output from MVA coincide;c^{jet}_{rel};Jets",24,-1.2,1.2));
    
    name = "h_coincidenceTest_flavourAntiB";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge for anti-B hadron flavour when truth and weight output from MVA coincide;c^{jet}_{rel};Jets",24,-1.2,1.2));
    
    name = "h_coincidenceTest_flavourB_fail";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge for B hadron flavour when truth and weight output from MVA DON'T coincide;c^{jet}_{rel};Jets",24,-1.2,1.2));
    
    name = "h_coincidenceTest_flavourAntiB_fail";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge for anti-B hadron flavour when truth and weight output from MVA DON'T coincide;c^{jet}_{rel};Jets",24,-1.2,1.2));
    
    name = "h_coincidenceTest_topB";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Coincidence between b from top and weight output from MVA;coincidence: 1 agree, 0 fail;Jets",4,0.,3.));
    
    name = "h_coincidenceTest_topAntiB";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Coincidence between anti-b from top and weight output from MVA;coincidence: 1 agree, 0 fail;Jets",4,0.,3.));
    
    name = "h_coincidenceTest_flavourB_oldDefinition";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Coincidence between B hadron flavour and weight output from old definition;coincidence: 1 agree, 0 fail;Jets",4,0.,3.));
    
    name = "h_coincidenceTest_flavourAntiB_oldDefinition";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Coincidence between anti-B hadron flavour and weight output from old definition;coincidence: 1 agree, 0 fail;Jets",4,0.,3.));
    
    name = "h_coincidenceTest_topB_oldDefinition";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Coincidence between b from top and weight output from old definition;coincidence: 1 agree, 0 fail;Jets",4,0.,3.));
    
    name = "h_coincidenceTest_topAntiB_oldDefinition";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Coincidence between anti-b from top and weight output from old definition;coincidence: 1 agree, 0 fail;Jets",4,0.,3.));
}
 
 void AnalyzerJetCharge::bookJetHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
 {
    TString name;
    
    name = "h_trueBJetPt";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet p_{T};Jet p_{T} ;Jets",40,0,200));

    name = "h_trueBJetTrackMultiplicity";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Multiplicity of the tracks in the true b jets;multiplicity;Events",30,0,30));
    
    name = "h_trueBJetTrackMultiplicity_eta_0_02";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Multiplicity of the tracks in the true b jets for eta [0-0.2];multiplicity;Events",30,0,30));
    
    name = "h_trueBJetTrackMultiplicity_eta_02_08";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Multiplicity of the tracks in the true b jets for eta (0.2-0.8];multiplicity;Events",30,0,30));
    
    name = "h_trueBJetTrackMultiplicity_eta_08_10";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Multiplicity of the tracks in the true b jets for eta (0.8-1.0];multiplicity;Events",30,0,30));
    
    name = "h_trueBJetTrackMultiplicity_eta_10_16";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Multiplicity of the tracks in the true b jets for eta (1.0-1.6];multiplicity;Events",30,0,30));
    
    name = "h_trueBJetTrackMultiplicity_eta_16_20";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Multiplicity of the tracks in the true b jets for eta (1.6-2.0];multiplicity;Events",30,0,30));
    
    name = "h_trueBJetTrackMultiplicity_eta_20";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Multiplicity of the tracks in the true b jets for eta >2.0;multiplicity;Events",30,0,30));
    
    name = "h_trueBJetTrackMultiplicity_pt_20_30";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Multiplicity of the tracks in the true b jets for pt [20 - 30];multiplicity;Events",30,0,30));
    
    name = "h_trueBJetTrackMultiplicity_pt_30_40";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Multiplicity of the tracks in the true b jets for pt (30 - 40];multiplicity;Events",30,0,30));
    
    name = "h_trueBJetTrackMultiplicity_pt_40_70";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Multiplicity of the tracks in the true b jets for pt (40 - 70];multiplicity;Events",30,0,30));
    
    name = "h_trueBJetTrackMultiplicity_pt_70_100";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Multiplicity of the tracks in the true b jets for pt (70 - 100];multiplicity;Events",30,0,30));
    
    name = "h_trueBJetTrackMultiplicity_pt_100_150";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Multiplicity of the tracks in the true b jets for pt (100 - 150];multiplicity;Events",30,0,30));
    
    name = "h_trueBJetTrackMultiplicity_pt_150";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Multiplicity of the tracks in the true b jets for pt>150;multiplicity;Events",30,0,30));
    
    name = "h_trueBJetEtaVsTrackMultiplicity";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"True b jet eta to track multiplicity correlation;track multiplicity;jet eta",30,0,30,40,-3.0,3.0));
    
    name = "h_trueBJetScalarChargeVsMultip";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"True b c_{rel} value in function of the jet track multiplicity;track multiplicity;c_{rel}",20,0,40,30, -1.5, 1.5));
    
    name = "h_trueBJetPtInitially";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet p_{T};Jet p_{T} ;Jets",40,0,200));
    
    name = "h_trueBJetCSVvalue";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet CSV;Jet CSV ;Jets",40,-2,2));
    
    name = "h_trueBJetHighestPtTrack";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T} value of the highest p_{T} track of each true b jet;p_{T} track;Jets",40,0,40));
    
    name = "h_trueBJetScalarChargeValidation";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet scalar p_{T}-weighted charge;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge10";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal p_{T} weighted charge with x = 1.0;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge8_IfOneSvInJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet scalar p_{T}-weighted charge if ONE secondary vertex in the event;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge8_WithoutZeroValue";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal p_{T} weighted charge with x = 0.8;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetScalarChargeVsSecondaryVertex";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "Comparison of Inclusive jet charge (x=0.8) vs Secondary Vertex Charge (x=0.8);c_{rel}^{jet};c_{rel}^{SV}",24,-1.2,1.2,24,-1.2,1.2));
    
    name = "h_trueBJetLeptonChargeVsSecondaryVertex";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "Comparison of leading lepton charge vs Secondary Vertex Charge;c^{lepton track};c^{SV}",24,-1.2,1.2,12,-6,6));
    
    name = "h_trueBJetRelCharge10";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet relative p_{T}-weighted charge;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetScalarChargeXWeighted";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal p_{T} weighted charge weighted by 0.7 (as in TOP-11-031);c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge2";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal p_{T} weighted charge with x = 0.2;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge4";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal p_{T} weighted charge with x = 04;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge6";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal p_{T} weighted charge with x = 0.6;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge8";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal p_{T} weighted charge with x = 0.8;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge12";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal p_{T} weighted charge with x = 1.2;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge14";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal p_{T} weighted charge with x = 1.4;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge16";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal p_{T} weighted charge with x = 1.6;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge18";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal p_{T} weighted charge with x = 1.8;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge20";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal p_{T} weighted charge with x = 2.0;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetRelPtTrack";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Relative p_{T} of all tracks in b (and anti-b) jets;Relative p_{T};Jets",40,0.,200.));
    
    name = "h_trueBJetRelCharge2";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse p_{T} weighted charge with x = 0.2;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetRelCharge4";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse p_{T} weighted charge with x = 04;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetRelCharge6";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse p_{T} weighted charge with x = 0.6;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetRelCharge8";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse p_{T} weighted charge with x = 0.8;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetRelCharge12";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse p_{T} weighted charge with x = 1.2;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetRelCharge14";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse p_{T} weighted charge with x = 1.4;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetRelCharge16";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse p_{T} weighted charge with x = 1.6;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetRelCharge18";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse p_{T} weighted charge with x = 1.8;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetRelCharge20";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse p_{T} weighted charge with x = 2.0;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetMaxRelPtTrack";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Relative p_{T} of highest p_{T} track in b (and anti-b) jets;Relative p_{T};Jets",40,0.,200.));
    
    
    name = "h_trueBJetBQuarkScalarCharge2";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal p_{T} weighted charge with x = 0.2;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkScalarCharge4";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal p_{T} weighted charge with x = 04;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkScalarCharge6";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal p_{T} weighted charge with x = 0.6;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkScalarCharge8";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal p_{T} weighted charge with x = 0.8;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkScalarCharge10";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"b quark to longitudinal b c_{rel} correlation ;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkScalarCharge12";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal p_{T} weighted charge with x = 1.2;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkScalarCharge14";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal p_{T} weighted charge with x = 1.4;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkScalarCharge16";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal p_{T} weighted charge with x = 1.6;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkScalarCharge18";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal p_{T} weighted charge with x = 1.8;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkScalarCharge20";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal p_{T} weighted charge with x = 2.0;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkScalarCharge2";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal p_{T} weighted charge with x = 0.2;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkScalarCharge4";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal p_{T} weighted charge with x = 04;;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkScalarCharge6";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal p_{T} weighted charge with x = 0.6;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkScalarCharge8";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal p_{T} weighted charge with x = 0.8;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetAntiBQuarkScalarCharge10";
    m_histogram[name] = store(new TH1D(prefix_+name+step," Anti-b quark to longitudinal b c_{rel} correlation;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkScalarCharge12";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal p_{T} weighted charge with x = 1.2;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkScalarCharge14";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal p_{T} weighted charge with x = 1.4;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkScalarCharge16";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal p_{T} weighted charge with x = 1.6;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkScalarCharge18";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal p_{T} weighted charge with x = 1.8;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkScalarCharge20";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal p_{T} weighted charge with x = 2.0;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkFromTopScalarCharge2";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from b-quark) longitudinal p_{T} weighted charge with x = 0.2;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkFromTopScalarCharge4";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from b-quark) longitudinal p_{T} weighted charge with x = 04;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkFromTopScalarCharge6";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from b-quark) longitudinal p_{T} weighted charge with x = 0.6;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkFromTopScalarCharge8";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from b-quark) longitudinal p_{T} weighted charge with x = 0.8;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkFromTopScalarCharge10";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"b quark to longitudinal b jet from top charge correlation ;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkFromTopScalarCharge12";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from b-quark) longitudinal p_{T} weighted charge with x = 1.2;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkFromTopScalarCharge14";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from b-quark) longitudinal p_{T} weighted charge with x = 1.4;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkFromTopScalarCharge16";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from b-quark) longitudinal p_{T} weighted charge with x = 1.6;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkFromTopScalarCharge18";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from b-quark) longitudinal p_{T} weighted charge with x = 1.8;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkFromTopScalarCharge20";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from b-quark) longitudinal p_{T} weighted charge with x = 2.0;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkFromTopScalarCharge2";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from anti-b quark) longitudinal p_{T} weighted charge with x = 0.2;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkFromTopScalarCharge4";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from anti-b quark) longitudinal p_{T} weighted charge with x = 04;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkFromTopScalarCharge6";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from anti-b quark) longitudinal p_{T} weighted charge with x = 0.6;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkFromTopScalarCharge8";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from anti-b quark) longitudinal p_{T} weighted charge with x = 0.8;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetAntiBQuarkFromTopScalarCharge10";
    m_histogram[name] = store(new TH1D(prefix_+name+step," Anti-b quark to longitudinal b jet from top charge correlation;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkFromTopScalarCharge12";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from anti-b quark) longitudinal p_{T} weighted charge with x = 1.2;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkFromTopScalarCharge14";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from anti-b quark) longitudinal p_{T} weighted charge with x = 1.4;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkFromTopScalarCharge16";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from anti-b quark) longitudinal p_{T} weighted charge with x = 1.6;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkFromTopScalarCharge18";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from anti-b quark) longitudinal p_{T} weighted charge with x = 1.8;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkFromTopScalarCharge20";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet from top (from anti-b quark) longitudinal p_{T} weighted charge with x = 2.0;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkRelCharge2";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse p_{T} weighted charge with x = 0.2;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkRelCharge4";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse p_{T} weighted charge with x = 04;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkRelCharge6";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse p_{T} weighted charge with x = 0.6;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkRelCharge8";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse p_{T} weighted charge with x = 0.8;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkRelCharge10";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"b quark to true b jet transverse charge correlation;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkRelCharge12";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse p_{T} weighted charge with x = 1.2;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkRelCharge14";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse p_{T} weighted charge with x = 1.4;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkRelCharge16";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse p_{T} weighted charge with x = 1.6;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkRelCharge18";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse p_{T} weighted charge with x = 1.8;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetBQuarkRelCharge20";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse p_{T} weighted charge with x = 2.0;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkRelCharge2";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse p_{T} weighted charge with x = 0.2;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkRelCharge4";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse p_{T} weighted charge with x = 04;;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkRelCharge6";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse p_{T} weighted charge with x = 0.6;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkRelCharge8";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse p_{T} weighted charge with x = 0.8;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetAntiBQuarkRelCharge10";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Anti-b quark to true b jet transverse charge correlation;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkRelCharge12";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse p_{T} weighted charge with x = 1.2;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkRelCharge14";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse p_{T} weighted charge with x = 1.4;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkRelCharge16";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse p_{T} weighted charge with x = 1.6;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkRelCharge18";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse p_{T} weighted charge with x = 1.8;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueAntiBJetBQuarkRelCharge20";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse p_{T} weighted charge with x = 2.0;c_{rel};Jets",24,-1.2,1.2));
    
    name = "h_trueBJetPtVsMaxPtTrack";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"true b jet p_{T} to its maximum p_{T} track relation; max track p_{T};jet p_{T} ",40,0.,40., 40, 0.,300));
    
    name = "h_trueBJetPtVsMaxRelPtTrack";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"true b jet p_{T} to its maximum p_{T}-rel track relation; max rel-p_{T} track; jet p_{T}",40,0.,300.,40, 0.,300.));
    
    name = "h_trueBJetPtVsMaxScalarPtTrack";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"true b jet p_{T} to its maximum p_{T}-scalar track scalaration; max scalar-p_{T} track; jet p_{T}",40, 0.,10000.,40,0.,300.));
    
    name = "h_trueBJetPtVsNumTracks";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"true b jet p_{T} to track multiplicity correlation;track multiplicity;jet p_{T}",30,0,30,40,0,300));
    
    name = "h_trueBJetMaxPtTrackVsNumTracks";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"true b jet highest p_{T} track  to track multiplicity correlation;track multiplicity;max track p_{T}",30,0.,30.,40,0.,40.));
    
    //  LEPTON TRACK IDENTIFICATION
    
    name = "h_trueBJetLeptonTracks";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jets containing a lepton among their tracks;PdgId of the track;Tracks",7,0,7));
    
    name = "h_trueBJetLeptonTrackMultiplicity";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Multiplicity of leptons in jets with at least one lepton track;Number of leptons; Tracks",10,0,10));
    
    name = "h_trueBJetLeptonTrackPt";
    m_histogram[name] = store(new TH1D(prefix_+name+step," p_{T}  of the leptons in jets with at least one lepton track;  p_{T}  of the lepton track; Tracks",40,0,40));
    
    name = "h_trueBJetLeptonTrackEta";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Eta of the leptons in jets with at least one lepton track; Eta of the lepton track; Tracks",40,-3,3));
    
    name = "h_trueBJetMuonTrackPt";
    m_histogram[name] = store(new TH1D(prefix_+name+step," p_{T}  of the muons in jets with at least one lepton track;  p_{T}  of the muon track; Tracks",40,0,40));
    
    name = "h_trueBJetMuonTrackEta";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Eta of the muons in jets with at least one lepton track; Eta of the muon track; Tracks",40,-3,3));
    
    name = "h_trueBJetElectronTrackPt";
    m_histogram[name] = store(new TH1D(prefix_+name+step," p_{T}  of the electrons in jets with at least one lepton track;  p_{T}  of the electron track; Tracks",40,0,40));
    
    name = "h_trueBJetElectronTrackEta";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Eta of the electrons in jets with at least one lepton track; Eta of the electron track; Tracks",40,-3,3));
    
    name = "h_trueBJetTrackMultiplicityIfLepton";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Multiplicity of the tracks in events with at least one lepton track;Track multiplicity; Jets",20,0,20));
    
    name = "h_trueBJetLeptonTrackCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of the lepton tracks; Lepton track charge;Tracks",2,-2,2));
    
    name = "h_trueBJetMuonTrackCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of the muon tracks; Muon track charge;Tracks",2,-2,2));
    
    name = "h_trueBJetElectronTrackCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of the electron tracks; Electron track charge;Tracks",2,-2,2));
    
    name = "h_trueBJetLeadingMuonTrackCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of the leading muon tracks; Muon track charge;Tracks",2,-2,2));
    
    name = "h_trueBJetLeadingElectronTrackCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of the leading electron tracks; Electron track charge;Tracks",2,-2,2));
    
    name = "h_trueBJetLeadingLeptonScalarChargeMatch";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Matching between the leading lepton track and jet longitudinal charge; Match of charge for leptons (1=true, 0=false); Tracks",2,0,2));
    
    name = "h_trueBJetLeadingMuonScalarChargeMatch";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Matching between the leading muon track and jet longitudinal charge; Match of charge for muons (1=true, 0=false); Tracks",2,0,2));
    
    name = "h_trueBJetLeadingElectronScalarChargeMatch";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Matching between the leading electron track and jet longitudinal charge; Match of charge for electrons (1=true, 0=false); Tracks",2,0,2));
    
    name = "h_trueBJetLeadingLeptonRelChargeMatch";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Matching between the leading lepton track and jet transverse charge; Match of charge (1=true, 0=false); Tracks",2,0,2));
    
    name = "h_trueBJetNonLeadingLeptonScalarChargeMatch";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Matching between the non-leading lepton track and jet longitudinal charge; Match of charge (1=true, 0=false); Tracks",2,0,2));
    
    name = "h_trueBJetNonLeadingLeptonRelChargeMatch";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Matching between the non-leading lepton track and jet transverse charge; Match of charge (1=true, 0=false); Tracks",2,0,2));
    
    name = "h_trueBJetLeptonTrackPtMultiplicity";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"Relation between track multiplicity and jet track p_{T} for lepton-tracked cases; p_{T}  of the jet; Track multiplicity",40,0.,400.,20,0.,20.));
    
    name = "h_trueBJetLeptonChargePt";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"Relation between p_{T} and charge of the lepton tracks; p_{T}  of the lepton track;Charge of the lepton track",40,0.,400.,2,-2.,2.));
    
    name = "h_trueBJetMuonChargePt";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"Relation between p_{T} and charge of the muon tracks; p_{T}  of the muon track;Charge of the muon track",40,0.,400.,2,-2.,2.));
    
    name = "h_trueBJetElectronChargePt";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"Relation between p_{T} and charge of the electron tracks; p_{T}  of the electron track;Charge of the electron track",40,0.,400.,2,-2.,2.));
    
    name = "h_trueBJetMuonRankInTrack";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Rank of the lepton track among the tracks; Rank of the track; Tracks",30,0,30));
    
    name = "h_trueBJetLeadingTrackParticleId";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"PdgId of the leading track in the jet; PdgId (2=e, 3=mu); Tracks",7,0,7));
    
    name = "h_trueBJetSubLeadingTrackParticleId";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"PdgId of the second-leading track in the jet; PdgId (2=e, 3=mu); Tracks",7,0,7));
    
    name = "h_trueBJetThirdLeadingTrackParticleId";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"PdgId of the third-leading track in the jet; PdgId (2=e, 3=mu); Tracks",7,0,7));
    
    name = "h_trueBJetFourthLeadingTrackParticleId";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"PdgId of the fourth-leading track in the jet; PdgId (2=e, 3=mu); Tracks",7,0,7));
    
    name = "h_trueBJetToTrackPtRatio";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Ratio between p_{T} of jet and p_{T} of track for all tracks;p_{T}Track/p_{T}Jet ;Tracks",40,0.,1.));
    
    name = "h_trueBJetToLeptonTrackPtRatio";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Ratio between p_{T} of jet and p_{T} of track for lepton tracks;p_{T}Track/p_{T}Jet ;Tracks",40,0.,1.));
    
    name = "h_trueBJetToMuonTrackPtRatio";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Ratio between p_{T} of jet and p_{T} of track for muon tracks;p_{T}Track/p_{T}Jet ;Tracks",40,0.,1.));
    
    name = "h_trueBJetToElectronTrackPtRatio";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Ratio between p_{T} of jet and p_{T} of track for lepton tracks;p_{T}Track/p_{T}Jet ;Tracks",40,0.,1.));
    
    name = "h_trueBJetToLeadingTrackPtRatio";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Ratio between p_{T} of jet and p_{T} of track for leading tracks;p_{T}Track/p_{T}Jet ;Tracks",40,0.,1.));
    
    name = "h_trueBJetToLeadingLeptonTrackPtRatio";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Ratio between p_{T} of jet and p_{T} of track for lepton leading tracks;p_{T}Track/p_{T}Jet for leading lepton ;Tracks",40,0.,1.));
    
    name = "h_trueBJetToLeadingMuonTrackPtRatio";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Ratio between p_{T} of jet and p_{T} of track for muon leading tracks;p_{T}Track/p_{T}Jet for leading muon ;Tracks",40,0.,1.));
    
    name = "h_trueBJetToLeadingElectronTrackPtRatio";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Ratio between p_{T} of jet and p_{T} of track for electron leading tracks;p_{T}Track/p_{T}Jet for leading electron;Tracks",40,0.,1.));
    
    name = "h_trueBJetToSubleadingLeptonTrackPtRatio";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Ratio between p_{T} of jet and p_{T} of track for lepton subleading tracks;p_{T}Track/p_{T}Jet for subleading lepton ;Tracks",40,0.,1.));
    
    name = "h_trueBJetToSubleadingMuonTrackPtRatio";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Ratio between p_{T} of jet and p_{T} of track for muon subleading tracks;p_{T}Track/p_{T}Jet for subleading muon ;Tracks",40,0.,1.));
    
    name = "h_trueBJetToSubleadingElectronTrackPtRatio";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Ratio between p_{T} of jet and p_{T} of track for electron subleading tracks;p_{T}Track/p_{T}Jet for subleading electron;Tracks",40,0.,1.));
    
    //TRUE LEPTONS
    
    name = "h_trueBJetGenLeptonPt";
    m_histogram[name] = store(new TH1D(prefix_+name+step," p_{T}  of the gen-leptons associated to a b-hadron; p_{T} ;Leptons",40,0.,400.));
    
    name = "h_trueBJetGenMuonPt";
    m_histogram[name] = store(new TH1D(prefix_+name+step," p_{T}  of the gen-muons associated to a b-hadron; p_{T} ;Leptons",40,0.,400.));
    
    name = "h_trueBJetGenLeptonPdgId";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"PdgId of the gen-leptons associated to a b-hadron;PdgId;Leptons",36,-18,18));
    
    name = "h_trueBJetAgreementMuons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Agreement between reco and gen muons detection;Agreement (ng+nr=0, g+nr=1, ng+r=2, g+r=3);Jets",4,0,4));
    
    name = "h_trueBJetAgreementElectrons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Agreement between reco and gen electron detection;Agreement (ng+nr=0, g+nr=1, ng+r=2, g+r=3);Jets",4,0,4));
    
    name = "h_trueBJetAgreementChargeMuons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Agreement between reco and gen muons charge;Charge (both b, both anti-b, mismatch);Jets",4,0,4));
    
    name = "h_trueBJetAgreementChargeElectrons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Agreement between reco and gen electrons charge;Charge (both b, both anti-b, mismatch);Jets",4,0,4));
    
    name = "h_trueCJetAgreementMuons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Agreement between reco and gen muons detection (for C);Agreement (ng+nr=0, g+nr=1, ng+r=2, g+r=3);Jets",4,0,4));
    
    name = "h_trueCJetAgreementElectrons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Agreement between reco and gen electron detection (for C);Agreement (ng+nr=0, g+nr=1, ng+r=2, g+r=3);Jets",4,0,4));
    
    name = "h_trueCJetAgreementChargeMuons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Agreement between reco and gen muons charge (for C);Charge (both b, both anti-b, mismatch);Jets",4,0,4));
    
    name = "h_trueCJetAgreementChargeElectrons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Agreement between reco and gen electrons charge (for C);Charge (both b, both anti-b, mismatch);Jets",4,0,4));
    
    name = "h_trueBJetGenToRecoPtDifferenceMuons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Difference between the p_{T} of all the gen-muons associated to a b-hadron and all the reco muon; p_{T} difference; Lepton combinations",10,0.,1.));
    
    name = "h_trueBJetGenToRecoPtDifferenceMuonsExtended";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Difference between the p_{T} of all the gen-muons associated to a b-hadron and all the reco muon; p_{T} difference; p_{T} difference; Lepton combinations",30,0.,30.));
    
    name = "h_trueBJetGenToRecoPhiDifferenceMuons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Difference between the Phi of all the gen-muons associated to a b-hadron and all the reco muon;Phi difference; Lepton combinations",40,-6.,2.));
    
    name = "h_trueBJetGenToRecoEtaDifferenceMuons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Difference between the Eta of all the gen-muons associated to a b-hadron and all the reco muon;Eta difference; Lepton combinations",40,0.,2.));
    
    name = "h_trueBJetGenToRecoDeltaRMuons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"DeltaR of all the gen-muons associated to a b-hadron and all the reco muon;DeltaR; Lepton combinations",40,0.,10.));
    
    //leptons
    name = "h_trueBJetBGenLeptonPtIfSameChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step," p_{T}  of the gen leptons coming from a b in a jet with same charge as b;p_{T} of b;Gen leptons",20,0,60));
    
    name = "h_trueBJetBGenLeptonPtIfOppositeChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step," p_{T}  of the gen leptons coming from a b in a jet with opposite charge as b;p_{T} of b;Gen leptons",20,0,60));
    
    name = "h_trueBJetCGenLeptonPtIfSameChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step," p_{T}  of the gen leptons coming from a c in a jet with same charge as c;p_{T} of c;Gen leptons",20,0,60));
    
    name = "h_trueBJetCGenLeptonPtIfOppositeChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step," p_{T}  of the gen leptons coming from a c in a jet with opposite charge as c;p_{T} of c;Gen leptons",20,0,60));
    
    //muons
    name = "h_trueBJetBGenMuonPtIfSameChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step," p_{T}  of the gen muons coming from a b in a jet with same charge as b;p_{T} of b;Gen muons",20,0,60));
    
    name = "h_trueBJetBGenMuonPtIfOppositeChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step," p_{T}  of the gen muons coming from a b in a jet with opposite charge as b;p_{T} of b;Gen muons",20,0,60));
    
    name = "h_trueBJetCGenMuonPtIfSameChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen muons coming from a c in a jet with same charge as c;p_{T} of c;Gen muons",20,0,60));
    
    name = "h_trueBJetCGenMuonPtIfOppositeChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen muons coming from a c in a jet with opposite charge as c;p_{T} of c;Gen muons",20,0,60));
    
    //electrons
    name = "h_trueBJetBGenElectronPtIfSameChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen electrons coming from a b in a jet with same charge as b;p_{T} of b;Gen electrons",20,0,60));
    
    name = "h_trueBJetBGenElectronPtIfOppositeChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen electrons coming from a b in a jet with opposite charge as b;p_{T} of b;Gen electrons",20,0,60));
    
    name = "h_trueBJetCGenElectronPtIfSameChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen electrons coming from a c in a jet with same charge as c;p_{T} of c;Gen electrons",20,0,60));
    
    name = "h_trueBJetCGenElectronPtIfOppositeChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen electrons coming from a c in a jet with opposite charge as c;p_{T} of c;Gen electrons",20,0,60));
    
    //leptons
    name = "h_trueBJetBLeadingGenLeptonPtIfSameChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen leading leptons coming from a b in a jet with same charge as b;p_{T} of b;Gen leptons",20,0,60));
    
    name = "h_trueBJetBLeadingGenLeptonPtIfOppositeChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen leading leptons coming from a b in a jet with opposite charge as b;p_{T} of b;Gen leptons",20,0,60));
    
    name = "h_trueBJetCLeadingGenLeptonPtIfSameChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen leading leptons coming from a c in a jet with same charge as c;p_{T} of c;Gen leptons",20,0,60));
    
    name = "h_trueBJetCLeadingGenLeptonPtIfOppositeChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen leading leptons coming from a c in a jet with opposite charge as c;p_{T} of c;Gen leptons",20,0,60));
    
    //muons
    name = "h_trueBJetBLeadingGenMuonPtIfSameChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen leading muons coming from a b in a jet with same charge as b;p_{T} of b;Gen muons",20,0,60));
    
    name = "h_trueBJetBLeadingGenMuonPtIfOppositeChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen leading muons coming from a b in a jet with opposite charge as b;p_{T} of b;Gen muons",20,0,60));
    
    name = "h_trueBJetCLeadingGenMuonPtIfSameChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen leading muons coming from a c in a jet with same charge as c;p_{T} of c;Gen muons",20,0,60));
    
    name = "h_trueBJetCLeadingGenMuonPtIfOppositeChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen leading muons coming from a c in a jet with opposite charge as c;p_{T} of c;Gen muons",20,0,60));
    
    //electrons
    name = "h_trueBJetBLeadingGenElectronPtIfSameChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen leading electrons coming from a b in a jet with same charge as b;p_{T} of b;Gen electrons",20,0,60));
    
    name = "h_trueBJetBLeadingGenElectronPtIfOppositeChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen leading electrons coming from a b in a jet with opposite charge as b;p_{T} of b;Gen electrons",20,0,60));
    
    name = "h_trueBJetCLeadingGenElectronPtIfSameChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen leading electrons coming from a c in a jet with same charge as c;p_{T} of c;Gen electrons",20,0,60));
    
    name = "h_trueBJetCLeadingGenElectronPtIfOppositeChargeAsJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T}  of the gen leading electrons coming from a c in a jet with opposite charge as c;p_{T} of c;Gen electrons",20,0,60));
    
    //hyerarchical c_{rel} plots
    
    name = "h_jetChargeHyerarchicalValue";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} ;c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueMuons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} (only muons taken as lepton option) ;c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueElectrons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} (only electrons taken as lepton option);c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueTruncated";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} (truncated) ;c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueMuonsTruncated";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} (truncated) (only muons taken as lepton option);c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueElectronsTruncated";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} (truncated) (only electrons taken as lepton option);c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueBJets";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} for b jets ;c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueBJetsMuons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} for b jets (only muons taken as lepton option);c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueBJetsElectrons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} for b jets (only electrons taken as lepton option) ;c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueBJetsTruncated";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} for b jets (truncated) ;c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueBJetsMuonsTruncated";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} for b jets (truncated) (only muons taken as lepton option);c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueBJetsElectronsTruncated";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} for b jets (truncated) (only electrons taken as lepton option);c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueAntiBJets";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} for anti-b jets;c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueAntiBJetsMuons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} for anti-b jets (only muons taken as lepton option);c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueAntiBJetsElectrons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} for anti-b jets (only electrons taken as lepton option);c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueAntiBJetsTruncated";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} for anti-b jets (truncated);c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueAntiBJetsMuonsTruncated";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} for anti-b jets (truncated) (only muons taken as lepton option);c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueAntiBJetsElectronsTruncated";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} for anti-b jets (truncated) (only electrons taken as lepton option);c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValue_preStep";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} - only separating between inclusive and leptonic ;c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueBJets_preStep";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} for b jets - only separating between inclusive and leptonic ;c_{rel}; Jets",24,-1.2,1.2));
    
    name = "h_jetChargeHyerarchicalValueAntiBJets_preStep";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Hyerarchical c_{rel} for anti-b jets - only separating between inclusive and leptonic ;c_{rel}; Jets",24,-1.2,1.2));

     
}


 void AnalyzerJetCharge::bookPfCandidateHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
 {
    TString name;
    
    name = "h_trueBJetPfTrackPt";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet track p_{T};p_{T} track;Tracks",40,0,40));
    
    name = "h_trueBJetTrackIPValueForPfCandidates";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Impact parameter value from PfCandidates for all b jet pfCandidates;Impact parameter value;PfCandidates",40,-0.4,0.4));
    
    name = "h_trueBJetTrackIPSignificanceForPfCandidates";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Impact parameter significance for all b jet pfCandidates;Impact parameter significance;PfCandidates",40,-60.,80.));
    
    name = "h_trueBJetTrackMaxIPValueForPfCandidates";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Max. impact parameter value from PfCandidates for all b jets;Impact parameter value;Jets",40,-0.4,0.4));
    
    name = "h_trueBJetTrackMaxIPSignificanceForPfCandidates";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Max. impact parameter significance for all b jets;Impact parameter significance;Jets",40,-60.,80.));
    
    name = "h_trueBJetTrackIPValueIfLepton";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Impact parameter value for b jet lepton-tracks ;Impact parameter value;Tracks",40,-0.4,0.4));
    
    name = "h_trueBJetTrackIPSignificanceIfLepton";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Impact parameter significance for b jet lepton-tracks ;Impact parameter significance;Tracks",40,-60.,80.));
    
    name = "h_trueBJetTrackIPValueIfMuon";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Impact parameter value for b jet muon-tracks ;Impact parameter value;Tracks",40,-0.4,0.4));
    
    name = "h_trueBJetTrackIPSignificanceIfMuon";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Impact parameter significance for b jet muon-tracks ;Impact parameter significance;Tracks",40,-60.,80.));
    
    name = "h_trueBJetTrackIPValueIfElectron";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Impact parameter value for b jet electron-tracks ;Impact parameter value;Tracks",40,-0.4,0.4));
    
    name = "h_trueBJetTrackIPSignificanceIfElectron";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Impact parameter significance for b jet electron-tracks ;Impact parameter significance;Tracks",40,-60.,80.));
    
    name = "h_trueBJetTrackIPValueVsSignificanceIfLepton";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"Impact parameter value vs significance for b jet muon-tracks;Impact parameter value; Impact parameter significance",40,-0.4,0.4,40,-60.,80.));
    
    name = "h_trueBJetTrackPtVsIPSignificanceIfLepton";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"track p_{T} vs impact parameter significance for all b jet muon-tracks;track p_{T} ; Impact parameter significance",40,0.,10.,40,-60.,80.));
    
    
 }
 
 void AnalyzerJetCharge::bookSelectedTrackHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
 {
    TString name;
    
    name = "h_trueBJetTrackIPValueForSelTracks";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Impact parameter value from selected tracks for all b jet selected tracks;Impact parameter value;Selected Tracks",40,-0.4,0.4));
    
    name = "h_trueBJetTrackIPSignificanceForSelTracks";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Impact parameter significance for all b jet selected tracks;Impact parameter significance;Selected Tracks",40,-60.,80.));
    
    name = "h_trueBJetTrackMaxIPValueForSelTracks";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Max. impact parameter value from selected tracks for all b jets;Impact parameter value;Jets",40,-0.4,0.4));
    
    name = "h_trueBJetTrackMaxIPSignificanceForSelTracks";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Max. impact parameter significance for all b jets;Impact parameter significance;Jets",40,-60.,80.));
    
    name = "h_trueBJetTrackMaxIPValueRank";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Rank of the max. impact parameter value track ;Rank position (p_{T});Jets",30,0,30));
    
    name = "h_trueBJetTrackMaxIPSignificanceRank";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Rank of the max. impact parameter significance track ;Rank position (p_{T});Jets",30,0,30));
    
    name = "h_trueBJetTrackIPValueVsSignificance";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"Impact parameter value vs significance for all b jet-tracks;Impact parameter value; Impact parameter significance",40,-0.4,0.4,40,-60.,80.));
    
    name = "h_trueBJetTrackPtVsIPSignificance";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"track p_{T} vs impact parameter significance for all b jet-tracks;track p_{T} ; Impact parameter significance",40,0.,10.,40,-60.,80.));
    
    name = "h_trueBJetTrackHighestPtVsIPSignificance";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"track with highest p_{T} vs impact parameter significance for all b jet-tracks;track p_{T} ; Impact parameter significance",40,0.,20.,40,-60.,80.));
    
    name = "h_trueBJetTrackPtVsHighestIPSignificance";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"track p_{T} vs highest impact parameter significance for all b jet-tracks;track p_{T} ; Impact parameter significance",40,0.,20.,40,-60.,80.));
    
    name = "h_trueBJetTrackSecondaryVertexMultiplicity";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex multiplicity;sv multiplicity;Jets",8,0.,8.));
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackCharge2";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - x=0.2;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackCharge4";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks- x=0.4;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackCharge6";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks- x=0.6;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackCharge8";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks- x=0.8;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackCharge10";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks- x=1.0;sv charge;Secondary Vertices",24,-1.2,1.2)); 
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackCharge12";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks- x=1.2;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackCharge14";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks- x=1.4;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackCharge16";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks- x=1.6;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackCharge18";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks- x=1.8;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackCharge20";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks- x=2.0;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackChargePfMatched";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - pfMatched;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetTrackSecondaryVertexTrackId";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex pfCandidate Id;sv track Id; PfCandidates",10,0.,10.));
    
    name = "h_trueBJetTrackSecondaryVertexTrackPt";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex pfCandidate p_{T} ;sv track p_{T} ;# pf tracks",40,0.,40.));
    
    name = "h_trueBJetTrackSecondaryVertexLeptonTrackPt";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex pfCandidate p_{T} ;sv lepton track p_{T} ;# pf tracks",40,0.,40.));
    
    name = "h_trueBJetTrackSecondaryVertexMuonTrackPt";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex muon pfCandidate p_{T} ;sv muon track p_{T} ;# pf tracks",40,0.,40.));
    
    name = "h_trueBJetTrackSecondaryVertexElectronTrackPt";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex electron pfCandidate p_{T} ;sv electron track p_{T} ;# pf tracks",40,0.,40.));
    
    name = "h_trueBJetTrackSecondaryVertexNonLeptonTrackPt";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex non-lepton pfCandidate p_{T} ;sv non-lepton track p_{T};# pf tracks",40,0.,40.));
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackChargeNonPfMatched";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - non pfMatched;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackCharge_IfOneSvInJet";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks if ONE sv track in the jet;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackMultiplicity";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex selected track multiplicity;sv track multiplicity;Secondary Vertices",15,0.,15.));
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackMultiplicityPfMatched";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex selected track multiplicity;sv track multiplicity;Secondary Vertices",15,0.,15.));
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackMultiplicityNonPfMatched";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex selected track multiplicity;sv track multiplicity;Secondary Vertices",15,0.,15.));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackNonWeightedCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex non-weighted charge;sv track charge;Secondary Vertices",12,-6,6));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity2";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 2;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity3";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 3;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity4";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 4;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity5";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 5;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity6";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 6;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity7";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 7;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity8";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 8;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity9";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 9;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity10";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 10;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity11";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 11;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity12";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 12;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity13";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 13;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackChargeForMultiplicity14";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 14;sv charge;Secondary Vertices",24,-1.2,1.2));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity2";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 2;sv charge;Secondary Vertices",12,-6,6));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity3";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 3;sv charge;Secondary Vertices",12,-6,6));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity4";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 4;sv charge;Secondary Vertices",12,-6,6));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity5";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 5;sv charge;Secondary Vertices",12,-6,6));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity6";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 6;sv charge;Secondary Vertices",12,-6,6));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity7";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 7;sv charge;Secondary Vertices",12,-6,6));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity8";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 8;sv charge;Secondary Vertices",12,-6,6));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity9";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 9;sv charge;Secondary Vertices",12,-6,6));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity10";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 10;sv charge;Secondary Vertices",12,-6,6));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity11";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 11;sv charge;Secondary Vertices",12,-6,6));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity12";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 12;sv charge;Secondary Vertices",12,-6,6));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity13";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 13;sv charge;Secondary Vertices",12,-6,6));
    
    name = "h_trueBJetSecondaryVertexSelectedTrackNonWeightedChargeForMultiplicity14";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks - multiplicity = 14;sv charge;Secondary Vertices",12,-6,6));
    
    name = "h_trueBJetTrackSecondaryVertexSelectedTrackChargeCoincidence";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex charge for selected tracks coincidence weighted and non-weighted;sv charge;Secondary Vertices",12,-6,6));
    
    // Flight distance ========================
    
    name = "h_secondaryVertexFlighDistanceValue";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex flight distance values;sv flight distance value;Secondary Vertices",45,0.,15.));
    
    name = "h_secondaryVertexMinFlighDistanceValue";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex min. flight distance value;sv flight distance value;Secondary Vertices",45,0.,15.));
    
    name = "h_secondaryVertexFlighDistanceValueFirstSV";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"First secondary vertex flight distance values;sv flight distance value;Secondary Vertices",45,0.,15.));
    
    name = "h_secondaryVertexFlighDistanceValueSecondSV";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Second secondary vertex flight distance values;sv flight distance value;Secondary Vertices",45,0.,15.));
    
    name = "h_secondaryVertexFlighDistanceValueThirdSV";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Third secondary vertex flight distance values;sv flight distance value;Secondary Vertices",45,0.,15.));
    
    name = "h_secondaryVertexFlighDistanceValueFirstSVIfJetTwoSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"First secondary vertex flight distance values if jet has two secondary vertices;sv flight distance value;Secondary Vertices",45,0.,15.));
    
    name = "h_secondaryVertexFlighDistanceValueSecondSVIfJetTwoSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Second secondary vertex flight distance values if jet has two secondary vertices;sv flight distance value;Secondary Vertices",45,0.,15.));
    
    name = "h_secondaryVertexFlighDistanceValueFirstSVIfJetThreeSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"First secondary vertex flight distance values if jet has three secondary vertices;sv flight distance value;Secondary Vertices",45,0.,15.));
    
    name = "h_secondaryVertexFlighDistanceValueSecondSVIfJetThreeSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Second secondary vertex flight distance values if jet has three secondary vertices;sv flight distance value;Secondary Vertices",45,0.,15.));
    
    name = "h_secondaryVertexFlighDistanceValueThirdSVIfJetThreeSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Third secondary vertex flight distance values if jet has three secondary vertices;sv flight distance value;Secondary Vertices",45,0.,15.));
    
    name = "h_secondaryVertexFlighDistanceValueFirstSVIfJetFourSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"First secondary vertex flight distance values if jet has four secondary vertices;sv flight distance value;Secondary Vertices",45,0.,15.));
    
    name = "h_secondaryVertexFlighDistanceValueSecondSVIfJetFourSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Second secondary vertex flight distance values if jet has four secondary vertices;sv flight distance value;Secondary Vertices",45,0.,15.));
    
    name = "h_secondaryVertexFlighDistanceValueThirdSVIfJetFourSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Third secondary vertex flight distance values if jet has four secondary vertices;sv flight distance value;Secondary Vertices",45,0.,15.));
    
    name = "h_secondaryVertexFlighDistanceValueFourthSVIfJetFourSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"First secondary vertex flight distance values if jet has four secondary vertices;sv flight distance value;Secondary Vertices",45,0.,15.));
    
    name = "h_secondaryVertexFlighDistanceSignificance";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Secondary vertex flight distance significances;sv flight distance significance;Secondary Vertices",45,0.,150.));
    
    name = "h_secondaryVertexFlighDistanceSignificanceFirstSV";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"First secondary vertex flight distance significances;sv flight distance significance;Secondary Vertices",45,0.,150.));
    
    name = "h_secondaryVertexFlighDistanceSignificanceSecondSV";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Second secondary vertex flight distance significances;sv flight distance significance;Secondary Vertices",45,0.,150.));
    
    name = "h_secondaryVertexFlighDistanceSignificanceThirdSV";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Third secondary vertex flight distance significances;sv flight distance significance;Secondary Vertices",45,0.,150.));
    
    name = "h_secondaryVertexFlighDistanceSignificanceFirstSVIfJetTwoSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"First secondary vertex flight distance significances if jet has two secondary vertices;sv flight distance significance;Secondary Vertices",45,0.,150.));
    
    name = "h_secondaryVertexFlighDistanceSignificanceSecondSVIfJetTwoSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Second secondary vertex flight distance significances if jet has two secondary vertices;sv flight distance significance;Secondary Vertices",45,0.,150.));
    
    name = "h_secondaryVertexFlighDistanceSignificanceFirstSVIfJetThreeSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"First secondary vertex flight distance significances if jet has three secondary vertices;sv flight distance significance;Secondary Vertices",45,0.,150.));
    
    name = "h_secondaryVertexFlighDistanceSignificanceSecondSVIfJetThreeSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Second secondary vertex flight distance significances if jet has three secondary vertices;sv flight distance significance;Secondary Vertices",45,0.,150.));
    
    name = "h_secondaryVertexFlighDistanceSignificanceThirdSVIfJetThreeSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Third secondary vertex flight distance significances if jet has three secondary vertices;sv flight distance significance;Secondary Vertices",45,0.,150.));
    
    name = "h_secondaryVertexFlighDistanceSignificanceFirstSVIfJetFourSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"First secondary vertex flight distance significances if jet has three secondary vertices;sv flight distance significance;Secondary Vertices",45,0.,150.));
    
    name = "h_secondaryVertexFlighDistanceSignificanceSecondSVIfJetFourSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Second secondary vertex flight distance significances if jet has three secondary vertices;sv flight distance significance;Secondary Vertices",45,0.,150.));
    
    name = "h_secondaryVertexFlighDistanceSignificanceThirdSVIfJetFourSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Third secondary vertex flight distance significances if jet has three secondary vertices;sv flight distance significance;Secondary Vertices",45,0.,150.));
    
    name = "h_secondaryVertexFlighDistanceSignificanceFourthSVIfJetFourSv";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"First secondary vertex flight distance significances if jet has three secondary vertices;sv flight distance significance;Secondary Vertices",45,0.,150.));
    
    name = "h_secondaryVertexMinFlightDistanceValueVsCSV";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "Secondary vertex min. flight distance value vs jet's csv value; sv flight distance value;Jet csv value",40,0.,15.,40,-1.5,1.5));
    
 }
 
 void AnalyzerJetCharge::bookOtherHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
 {
    TString name;
    
    name = "h_trueBJetBQuarkCharge";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "b-quark flavour to reco c_{rel} correlation;c_{rel};b-quark flavour",10,-1.5,1.5,40,-30.,30.));
    
    name = "h_trueBJetTopCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Top b quark to reco c_{rel} correlation;c_{rel};Jets",10,-1.5,1.5));
    
    name = "h_trueBJetAntiTopCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Top anti-b quark to reco c_{rel} correlation;c_{rel};Jets",10,-1.5,1.5));
    
    name = "h_trueBJetHiggsCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Higgs b quark to reco c_{rel} correlation;c_{rel};Jets",10,-1.5,1.5));
    
    name = "h_trueBJetAntiHiggsCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Higgs anti-b quark to reco c_{rel} correlation;c_{rel};Jets",10,-1.5,1.5));
    
    //b-hadron plots
    
    name = "h_trueBJetBPlusCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"B+ hadron to reco c_{rel} correlation;c_{rel};Jets",10,-1.5,1.5));
    
    name = "h_trueBJetBMinusCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"B- hadron to reco c_{rel} correlation;c_{rel};Jets",10,-1.5,1.5));
    
    name = "h_trueBJetBZeroCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"B0 hadron to reco c_{rel} correlation;c_{rel};Jets",10,-1.5,1.5));
    
    name = "h_trueBJetAntiBZeroCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Anti-B0 hadron to reco c_{rel} correlation;c_{rel};Jets",10,-1.5,1.5));
    
    //b-oscillations plots
    
    name = "h_trueBJet_bbarToBPlus";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of bbar going into BPlus",3,0.,2.));
    
    name = "h_trueBJet_bToBPlus";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of b going into BPlus",3,0.,2.));
    
    name = "h_trueBJet_bbarToBZero";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of bbar going into BZero",3,0.,2.));
    
    name = "h_trueBJet_bToBZero";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of b going into BZero",3,0.,2.));
    
    name = "h_trueBJet_bToBMinus";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of b going into BMinus",3,0.,2.));
    
    name = "h_trueBJet_bbarToBMinus";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of bbar going into BMinus",3,0.,2.));
    
    name = "h_trueBJet_bToAntiBZero";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of b going into AntiBZero",3,0.,2.));
    
    name = "h_trueBJet_bbarToAntiBZero";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of bbar going into AntiBZero",3,0.,2.));
    
    name = "h_trueBJet_bbarToBcPlus";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of bbar going into BcPlus",3,0.,2.));
    
    name = "h_trueBJet_bToBcPlus";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of b going into BcPlus",3,0.,2.));
    
    name = "h_trueBJet_bbarToBsZero";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of bbar going into BsZero",3,0.,2.));
    
    name = "h_trueBJet_bToBsZero";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of b going into BsZero",3,0.,2.));
    
    name = "h_trueBJet_bToBcMinus";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of b going into BcMinus",3,0.,2.));
    
    name = "h_trueBJet_bbarToBcMinus";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of bbar going into BcMinus",3,0.,2.));
    
    name = "h_trueBJet_bToAntiBsZero";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of b going into AntiBsZero",3,0.,2.));
    
    name = "h_trueBJet_bbarToAntiBsZero";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of bbar going into AntiBsZero",3,0.,2.));
    
    name = "h_test_jetPfCandidateChargeForB";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of pfCandidates from B;charge;pfCandidates",24,-1.2,1.2));
    
    name = "h_test_jetPfCandidateChargeForAntiB";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of pfCandidates from anti B;charge;pfCandidates",24,-1.2,1.2));
    
    name = "h_test_jetPfCandidateCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of pfCandidates;charge;pfCandidates",24,-1.2,1.2));
    
    name = "h_selectedTrackJetCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Charge of the jet using selTracks;charge;Selected Tracks",24,-1.2,1.2));
    
    // tests for closest z vertex
    name = "h_pfCandidateNonZeroWeights";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Vertex weight if different from zero;vertex weight;pfCandidates",10,0.,1.));
    
    name = "h_pfCandidateNonZeroWeightsIndexZero";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"# of non zero weight vertices = 0;vertex weight;pfCandidates",10,0.,1.));
    
    name = "h_pfCandidateNonZeroWeightsIndexNonZero";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"# of non zero weight vertices != 0;vertex weight;pfCandidates",10,0.,1.));
    
    name = "h_pfCandidateZeroWeightsIndexZero";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"# of zero weight vertices = 0;vertex weight;pfCandidates",10,0.,1.));
    
    name = "h_pfCandidateZeroWeightsIndexNonZero";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"# of zero weight vertices != 0;vertex weight;pfCandidates",10,0.,1.));
    
    name = "h_pfCandidateMaxWeights";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Vertex max. weight;vertex weight;pfCandidates",10,0.,1.));
    
    name = "h_pfCandidateMaxWeightsVsIndex";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"Vertex index to which max. weight corresponds;weight;vertex index",10,0.,1.,22,-2.,20.));
    
    name = "h_pfCandidateInitialAmount";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of PfCandidates before weight check cut ;# pfCandidates; pfCandidate",2,0.,2.));
    
    name = "h_pfCandidateInitialAmountForAntiB";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of PfCandidates from antiB before weight check cut ;# pfCandidates; pfCandidate",2,0.,2.));
    
    name = "h_pfCandidateInitialAmountForB";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Number of PfCandidates from B before weight check cut ;# pfCandidates; pfCandidate",2,0.,2.));
    
    name = "h_pfCandidateMaxWeightIfSeveralFound";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Pf candidate weights if several found ;pfCandidate weight; pfCandidate",10,0.,1.));
    
    name = "h_pfCandidateMaxWeightIfSeveralFoundForB";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Pf candidate weights if several found ;pfCandidate weight; pfCandidate",10,0.,1.));
    
    name = "h_pfCandidateMaxWeightIfSeveralFoundForAntiB";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Pf candidate weights if several found ;pfCandidate weight; pfCandidate",10,0.,1.));
    
    name = "h_pfCandidateMaxWeightIfOnlyOneFound";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Pf candidate weights if only one found;# pfCandidates; pfCandidate",10,0.,1.));
    
    name = "h_pfCandidateMaxWeightIfOnlyOneFoundForB";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Pf candidate weights if only one found;# pfCandidates; pfCandidate",10,0.,1.));
    
    name = "h_pfCandidateMaxWeightIfOnlyOneFoundForAntiB";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Pf candidate weights if only one found;# pfCandidates; pfCandidate",10,0.,1.));
    
    name = "h_pfCandidateMaxWeightIfOnlyNoneFound";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Pf candidate weights if none found;# pfCandidates; pfCandidate",10,0.,1.));
    
    name = "h_pfCandidateMaxWeightIfOnlyNoneFoundForB";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Pf candidate weights if none found;# pfCandidates; pfCandidate",10,0.,1.));
    
    name = "h_pfCandidateMaxWeightIfOnlyNoneFoundForAntiB";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Pf candidate weights if none found;# pfCandidates; pfCandidate",10,0.,1.));
    
    name = "h_pfCandidateMinZ";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Closest vertex z distance;z; pfCandidate",200,0.,20.));
    
    name = "h_pfCandidateZVertexNonZero";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Closest vertex z distance;z; pfCandidate",2,0.,2.));
    
    // TEST jet charge from function
    name = "h_jetChargeFromFunction";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal p_{T} weighted charge with x = 0.8;c_{rel};Jets",24,-1.2,1.2));
    
    
 }

//==========================================OTHER FUNTIONS===============================================

bool AnalyzerJetCharge::checkDecreasingOrderOfAList(const VLV& vector)
{
    bool isNotPtOrdered = false;
    for (size_t i=0;i!=vector.size()-1;++i)
    {
        for(size_t j=i+1;j!=vector.size();++j)
        {
            if (vector.at(j).pt()>vector.at(i).pt()) isNotPtOrdered=true;
            break;
        }
        if (isNotPtOrdered==true) break;
    }
    return isNotPtOrdered;
}

int AnalyzerJetCharge::selTrackMatchToPfIndex (const LV& pfCandidate, const std::vector<LV>& selectedTracks, const int pfCharge, const std::vector<int>& selCharge, const int pfTrackJetIndex, const std::vector<int>& selectedTrackIndex, const std::vector <int>& ptOrderedSelTrackIdx)
{
    int matchedSelectedTrackIdx = -1;
    for (size_t i_jetSelectedTrack=0; i_jetSelectedTrack!=selectedTracks.size();i_jetSelectedTrack++)
    {
        int jetSelectedTrackIdx = ptOrderedSelTrackIdx.at(i_jetSelectedTrack);
        int jetSelectedTrackJetIdx = selectedTrackIndex.at(jetSelectedTrackIdx);
        
        if (pfTrackJetIndex!=jetSelectedTrackJetIdx) continue;
        if(std::abs(pfCandidate.eta()-selectedTracks.at(jetSelectedTrackIdx).eta())>1.e-6 || std::abs(pfCandidate.phi()-selectedTracks.at(jetSelectedTrackIdx).phi())>1.e-6 || std::abs(pfCandidate.pt()-selectedTracks.at(jetSelectedTrackIdx).pt())>1.e-6 || pfCharge!=selCharge.at(jetSelectedTrackIdx)) continue;
        
        matchedSelectedTrackIdx = jetSelectedTrackIdx;
        break;
    }
    return matchedSelectedTrackIdx;
}

std::vector<int> AnalyzerJetCharge::leptonToJetMlbCalculator(const double massThreshold, const LV& lepton, const int leptonCharge, const VLV& recoJets, const std::vector<int>& jetIndices)
{
//     std::cout<<"Entering the function for lepton with index = "<<leptonCharge<<" and an amount of jets = "<<jetIndices.size()<<std::endl;
    std::vector<double> m;
    std::vector<int> mlbOutput;
    
    for(size_t iJet=0; iJet!=jetIndices.size();++iJet)  m.push_back((lepton+recoJets.at(jetIndices.at(iJet))).mass());
    
//     std::cout<<"The systems have the masses "<<m.at(0)<<" and "<<m.at(1)<<std::endl;
    
    // Both masses can't be under threshold mass
    if (m.at(0)<massThreshold && m.at(1)<massThreshold) mlbOutput.push_back(0);
    else if (m.at(0)>massThreshold && m.at(1)>massThreshold) mlbOutput.push_back(0);
    
    else if (m.at(0)<massThreshold&&m.at(1)>massThreshold) 
    {
        // The first entry corresponds to the lepton index
        mlbOutput.push_back(leptonCharge);
        
        // Push back the indices of the jets: 1st the one under the mass threshold, 2nd the one above it
        mlbOutput.push_back(jetIndices.at(0));
        mlbOutput.push_back(jetIndices.at(1));
    }
    
    else if (m.at(0)>massThreshold&&m.at(1)<massThreshold)
    {
        // The first entry corresponds to the lepton index
        mlbOutput.push_back(leptonCharge);
        
        // Push back the indices of the jets: 1st the one under the mass threshold, 2nd the one above it
        mlbOutput.push_back(jetIndices.at(1));
        mlbOutput.push_back(jetIndices.at(0));
    }
    
//     if (mlbOutput.size()==0) std::cout<<"Mlb size is "<<mlbOutput.size()<<" with entries: "<<mlbOutput.at(0)<<std::endl;
//     if (mlbOutput.size()>1) std::cout<<"Mlb size is "<<mlbOutput.size()<<" with entries: "<<mlbOutput.at(0)<<", "<<mlbOutput.at(1)<<", "<<mlbOutput.at(2)<<std::endl;
    return mlbOutput;
}

int AnalyzerJetCharge::genJetIdOfRecoJet ( const LV& recoJet, const VLV& genJets, const float dR_max )
{
    
    float dR_min = 999.9;
    int genJetId = -1;
    // Loop over gen jets to find the closest one to reco jet
    for(size_t iJet=0, nJets=genJets.size(); iJet<nJets; iJet++) {
        float dR = ROOT::Math::VectorUtil::DeltaR(recoJet,genJets.at(iJet));
        if(dR>=dR_min) continue;
        if(dR>=dR_max) continue;
        dR_min = dR;
        genJetId = iJet;
    }
    
    return genJetId;
}


std::vector< int > AnalyzerJetCharge::bHadIdsInGenJet ( const int jetId, const std::vector< int >& hadJetIndices )
{
    
    std::vector<int> hadIndices;
    
    if(jetId<0) return hadIndices;
    
    for(size_t iHad=0, nHads=hadJetIndices.size(); iHad<nHads; iHad++) {
        if(hadJetIndices.at(iHad)!=jetId) continue;
        hadIndices.push_back(iHad);
    }
    
    return hadIndices;
}


std::vector< int > AnalyzerJetCharge::bHadFlavoursInGenJet ( const int jetId, const std::vector< int >& hadJetIndices, const std::vector< int >& hadFlavours, const bool absFlavour )
{
    
    std::vector<int> flavours;
    std::vector<int> hadIndices = bHadIdsInGenJet(jetId, hadJetIndices);
    
    for(size_t iHad=0, nHads=hadIndices.size(); iHad<nHads; iHad++) {
        int hadId = hadIndices.at(iHad);
        int flavour = hadFlavours.at(hadId);
        if(absFlavour) flavour = std::abs(flavour);
        
        putUniquelyInVector(flavours, flavour);
    }
    
    return flavours;
}

bool AnalyzerJetCharge::isInVector(const std::vector<int>& idVector, const int id)
{
    bool isIn = std::find(idVector.begin(), idVector.end(), id) != idVector.end();
    
    return isIn;
}

bool AnalyzerJetCharge::putUniquelyInVector(std::vector<int>& vector, const int id)
{
    if(isInVector(vector, id)) return false;
    
    vector.push_back(id);
    return true;
}

double AnalyzerJetCharge::trackMultiplicityWeight(const double& m, const double& n, int jetIndex, const std::vector<int>& jetPfCandidateTrackIndex, const std::vector<int>& pfCandidateVertexId)
{
    int multiplicity = 0;
    for (size_t i=0;i!=jetPfCandidateTrackIndex.size();++i)
    {
        //check if the track is matched to a selected jet and in case it is, add one to the multiplicity.
        int trueMatched = 0;
        if (jetIndex!=jetPfCandidateTrackIndex.at(i)) trueMatched = -1;
        if (trueMatched == -1) continue;
        // Remove tracks not corresponding to primary vertex
        if (pfCandidateVertexId.at(i) == -1 || pfCandidateVertexId.at(i) == 2 || pfCandidateVertexId.at(i) == 3) continue;
        ++multiplicity;
    }
    
    double y = m*multiplicity+n;
    return y;
}

unsigned int AnalyzerJetCharge::calculateMultiplicity(const std::vector<int>& collection, int jetIndex)
{
    unsigned int multiplicity = 0; 
    
    for(size_t i=0; i<collection.size(); ++i) 
    {
        if(collection.at(i)!=static_cast<int>(jetIndex)) continue;
        
        for(size_t j=i; j<collection.size(); ++j)
        {
            if(collection.at(j)!=collection.at(i)) continue;
            ++multiplicity;
        }
        break; 
    }
    
    return multiplicity;
    
}

double AnalyzerJetCharge::ptWeightedJetChargeX (const int jetId, const LV& recoJet, const double& x, const std::vector<int>& pfCandidateJetIndex, const VLV& pfCandidates, const std::vector<int>& pfCandidateCharge, const std::vector<int>& pfCandidateVertexId)
{
    // Access jet momentum information
    double jetTrueBPx = recoJet.px();
    double jetTrueBPy = recoJet.py();
    double jetTrueBPz = recoJet.pz();
    
    // Define relevant variables for c_{rel} calculation
    double sumMomentum = 0.;
    double sumMomentumQ = 0.;
    
    for (size_t iCandidate=0;iCandidate!=pfCandidates.size();++iCandidate)
    {
        // Check that the pfCandidate corresponds to the jet
        if (jetId!=pfCandidateJetIndex.at(iCandidate)) continue;
        // Remove tracks not corresponding to primary vertex
        if (pfCandidateVertexId.at(iCandidate) == -1 || pfCandidateVertexId.at(iCandidate) == 2 || pfCandidateVertexId.at(iCandidate) == 3) continue;
        
        // Access pfCandidate mometum and charge information
        const double constituentTrueBPx = pfCandidates.at(iCandidate).px();
        const double constituentTrueBPy = pfCandidates.at(iCandidate).py();
        const double constituentTrueBPz = pfCandidates.at(iCandidate).pz();
        const double product = constituentTrueBPx*jetTrueBPx + constituentTrueBPy*jetTrueBPy + constituentTrueBPz*jetTrueBPz;
        
        int charge = pfCandidateCharge.at(iCandidate);
        
        // Sum over all the pfCandidates
        const double productPow = std::pow(product, x);
        sumMomentum += productPow;
        sumMomentumQ += static_cast<double>(charge)*productPow;
    }
    
    // Obtain the jet c_{rel}
    const double ptWeightedJetChargeXValue(sumMomentum>0 ? sumMomentumQ/sumMomentum : 0);
    return ptWeightedJetChargeXValue;
}

double AnalyzerJetCharge::trackMultiplicityWeightPerEvent (const std::vector<int>& jetIndices, const double& m, const double& n, const std::vector<int>& jetPfCandidateTrackIndex, const std::vector<int>& pfCandidateVertexId)
{
    double eventMultiplicity = 1.;
    for (size_t iJet=0; iJet!=jetIndices.size(); ++iJet)
    {
        int jetIndex = jetIndices.at(iJet);
        double jetWeight = trackMultiplicityWeight(m, n, jetIndex, jetPfCandidateTrackIndex, pfCandidateVertexId);
        eventMultiplicity *= jetWeight;
    }
    
    return eventMultiplicity;
}




