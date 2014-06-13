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
#include<TH1I.h>
#include <TH1D.h>
#include <TH2D.h>
#include <Math/VectorUtil.h>
#include <TProfile.h>
#include <TObjString.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>


#include "../include/AnalyzerJetCharge.h"
#include "../include/analysisStructs.h"
#include "../include/JetCategories.h"
#include "../include/higgsUtils.h"
#include "../include/HiggsAnalysis.h"
#include "../../common/include/AnalysisBase.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/classes.h"



AnalyzerJetCharge::AnalyzerJetCharge(const std::vector<TString>& selectionStepsNoCategories,
                                     const std::vector<TString>& stepsForCategories,
                                     const JetCategories* jetCategories,
                                     const int debug):
AnalyzerBase ("jetCharge_", selectionStepsNoCategories, stepsForCategories, jetCategories),debug_(0)
{
    debug_ = debug;
    std::cout<<"--- Beginning setting up jetChargeAnalyzer\n";
    std::cout<<"=== Finishing setting up jetChargeAnalyzer\n\n";
}
                                     


void AnalyzerJetCharge::fillHistos(const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                   const TopGenObjects& topGenObjects, const HiggsGenObjects&,
                                   const KinRecoObjects&,
                                   const tth::RecoObjectIndices&, const tth::GenObjectIndices&,
                                   const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                   const double& weight, const TString&,
                                   std::map<TString, TH1*>& m_histogram)
  
{
    TString name;

    // Extracting input data to more comfortable variables
    
    //   VLV* allGenJets = commonGenObjects.allGenJets_;
    const VLV& genJets = (commonGenObjects.valuesSet_) ? *commonGenObjects.allGenJets_ : VLV(0);
    const VLV& allJets = *recoObjects.jets_; 
    //const std::vector<LV>& bHadPlusMothers = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadPlusMothers_ : std::vector<LV*>(0);
    const std::vector<LV>& bHadPlusMothers = *topGenObjects.genBHadPlusMothers_;
    const std::vector<double>& jetChargeRelativePtWeighted = *recoObjects.jetChargeRelativePtWeighted_;
    
    //b-hadron + c-hadron information
    const std::vector<int>& bHadIndex = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadIndex_ : std::vector<int>(0);
    const std::vector<int>& bHadJetIndex = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadJetIndex_ : std::vector<int>(0);
    const std::vector<int>& bHadFlavour = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadFlavour_ : std::vector<int>(0);
    const std::vector<int>& bHadPlusMothersPdgId = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadPlusMothersPdgId_ : std::vector<int>(0);
    //const std::vector<int>& cHadJetIndex = (topGenObjects.valuesSet_) ? *topGenObjects.genCHadJetIndex_ : std::vector<int>(0);
    //const std::vector<int>& cHadFromBHad = (topGenObjects.valuesSet_) ? *topGenObjects.genCHadFromBHadron_ : std::vector<int>(0);
    
    //genLepton information
    const std::vector<int> genBHadLeptonIndex = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadLeptonIndex_ : std::vector<int>(0);
    const std::vector<int> genBHadLeptonHadronIndex = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadLeptonHadronIndex_ : std::vector<int>(0);
    const std::vector<int> genBHadLeptonViaTau = (topGenObjects.valuesSet_) ? *topGenObjects.genBHadLeptonViaTau_ : std::vector<int>(0);
    const std::vector<int> genCHadLeptonHadronIndex = (topGenObjects.valuesSet_) ? *topGenObjects.genCHadLeptonHadronIndex_ : std::vector<int>(0);
    const std::vector<int> genCHadLeptonViaTau = (topGenObjects.valuesSet_) ? *topGenObjects.genCHadLeptonViaTau_ : std::vector<int>(0);
    
    //specific selected tracks information (tracks with quality requirements applied already at ntuple level) 
    const std::vector<LV>& jetSelectedTrack = *recoObjects.jetSelectedTrack_;
    const std::vector<double>& jetSelectedTrackIPValue = *recoObjects.jetSelectedTrackIPValue_;
    const std::vector<double>& jetSelectedTrackIPSignificance = *recoObjects.jetSelectedTrackIPSignificance_;
    const std::vector<int>& jetSelectedTrackCharge = *recoObjects.jetSelectedTrackCharge_;
    const std::vector<int>& jetSelectedTrackIndex = *recoObjects.jetSelectedTrackIndex_;
    
    //specific track information (from pfcandidates - no quality cuts required at ntuple level)
    const std::vector<int>& jetPfCandidateTrackCharge = *recoObjects.jetPfCandidateTrackCharge_;
    const std::vector<int>& jetPfCandidateTrackId = *recoObjects.jetPfCandidateTrackId_;
    const std::vector<LV>& jetPfCandidateTrack = *recoObjects.jetPfCandidateTrack_;
    const std::vector<int>& jetPfCandidateTrackIndex = *recoObjects.jetPfCandidateTrackIndex_;
    
    
    
    // =============================================================Now do calculations and filling of histograms===================================================================
  
  
    //create a new jet selection-> 20GeV 

    constexpr double lowerJetPtCUT = 20.;
    //NOTE For the moment I'm applying no explicit eta cut because it's applied at ntuple level + no tracking is collected at eta>2.5. But bear it in mind!
    // Get jet indices, apply selection cuts and order them by pt (beginning with the highest value)                                                                                                                
    const VLV& lowerPtCUTjets = *recoObjects.jets_;                                                
    std::vector<int> lowerPtCUTJetIdx = common::initialiseIndices(lowerPtCUTjets);                        
    common::selectIndices(lowerPtCUTJetIdx, lowerPtCUTjets, common::LVpt, lowerJetPtCUT);                
    common::orderIndices(lowerPtCUTJetIdx, lowerPtCUTjets, common::LVpt);                  //pt ordered jets            
    
    //create a new b jet selection-> 20GeV
    constexpr double BtagWP = 0.244;

    const std::vector<double>& jetBTagCSV = *recoObjects.jetBTagCSV_;
    std::vector<int> bjetIndices = lowerPtCUTJetIdx;
    common::selectIndices(bjetIndices, jetBTagCSV, BtagWP);
    common::orderIndices(bjetIndices, jetBTagCSV);

    //create new indices for the tracks->cut on different pt values
    constexpr double ptTracksCUT = 0.;
    
    //pt-order the selected tracks (they are not pt ordered by default)
    std::vector<int> ptOrderedSelTrackIdx = common::initialiseIndices(jetSelectedTrack);
    common::orderIndices(ptOrderedSelTrackIdx, jetSelectedTrack, common::LVpt);
    
    //pt-order the pf tracks (even if they are in principle already ordered)
    std::vector<int> ptOrderedJetTrackIdx = common::initialiseIndices(jetPfCandidateTrack);
    common::orderIndices(ptOrderedJetTrackIdx, jetPfCandidateTrack, common::LVpt);
    
    //define the random generator used later on for splitting test and train trees for the mva.
    TRandom3 *nRandom = new TRandom3(0); 
    UInt_t nTrees = nRandom->Integer(2);
    
    std::vector<int> recoIndex;
    std::vector<int> genIndex;
    std::vector<int> flavour;
    std::vector<int> hadFlav;
    
    std::vector<int> genLeptonFlavour;
    std::vector<double> genLeptonPt;
    std::vector<double> genMuonPt;
    std::vector<LV> genMuonLV;
    
    std::vector <double> charge;
    std::vector<double> chargePerJet;
    
    std::vector<double> trueBJetScalarCharge2V;
    std::vector<double> trueBJetScalarCharge4V;
    std::vector<double> trueBJetScalarCharge6V;
    std::vector<double> trueBJetScalarCharge8V;
    std::vector<double> trueBJetScalarCharge10V;
    std::vector<double> trueBJetScalarCharge12V;
    std::vector<double> trueBJetScalarCharge14V;
    std::vector<double> trueBJetScalarCharge16V;
    std::vector<double> trueBJetScalarCharge18V;
    std::vector<double> trueBJetScalarCharge20V;
    
    std::vector<double> trueBJetRelCharge2V;
    std::vector<double> trueBJetRelCharge4V;
    std::vector<double> trueBJetRelCharge6V;
    std::vector<double> trueBJetRelCharge8V;
    std::vector<double> trueBJetRelCharge10V;
    std::vector<double> trueBJetRelCharge12V;
    std::vector<double> trueBJetRelCharge14V;
    std::vector<double> trueBJetRelCharge16V;
    std::vector<double> trueBJetRelCharge18V;
    std::vector<double> trueBJetRelCharge20V;

    
    //vectors specific for the mva
    std::vector<int> trueBJetMaxChargeTrackV;
    std::vector<double> trueBJetMaxPtTrackV;
    std::vector<int> trueBJetTrackMultiplicityV;
    std::vector<int> trueBJetMaxPtChargeTrackV;
    std::vector<double> trueBJetPtV;
    std::vector<double> trueBJetPtRatioV;
    //std::vector<int> isMuonV;
   
    
    bool fillTree = false;
    

    m_histogram["h_trueBJetLeptonTracks"]->SetMinimum(0);
    m_histogram["h_trueBJetLeptonTrackCharge"]->SetMinimum(0);
    m_histogram["h_trueBJetLeadingLeptonScalarChargeMatch"]->SetMinimum(0);
    m_histogram["h_trueBJetLeadingLeptonRelChargeMatch"]->SetMinimum(0);
    m_histogram["h_trueBJetNonLeadingLeptonScalarChargeMatch"]->SetMinimum(0);
    m_histogram["h_trueBJetNonLeadingLeptonRelChargeMatch"]->SetMinimum(0);
    
    //std::cout<<"======================NEW EVENT========================"<<std::endl;
    
    //RECO JETS (+ GEN JETS) LOOPS ======================================================================================================================
    for(size_t i_jet=0;i_jet!=lowerPtCUTJetIdx.size();i_jet++)
    {
        int jetIdx = lowerPtCUTJetIdx.at(i_jet);
        LV jets = allJets.at(jetIdx);
        if (debug_==2) std::cout<<"jets are = "<<jetIdx<<std::endl;
      
        //std::cout<<"test1"<<std::endl;
        //================================HERE B-QUARK-B-HADRON THINGS================================
        
        int isMatched = -1;
        //bool isGenLepton = false;
        bool isGenMuon = false;
        bool isMuon = false;
        bool isGenPlusMuon = false;
        bool isGenMinMuon = false;
        bool isGenElectron = false;
        bool isElectron = false;
        bool isGenPlusElectron = false;
        bool isGenMinElectron = false;
        
        //matching distance between reco and gen jets
        const float dR_max = 0.5;
        
        //call the matching function between reco and gen jets->returns an index. 
        int genJetIndex = genJetIdOfRecoJet (jets,genJets,dR_max);
        
        //std::cout<<"test1.1"<<std::endl;
        //TEST WHEN COMPARING TO DATA, WE CAN'T REQUIRE GENERATOR LEVEL INFORMATION
        if (genJetIndex==-1) continue;
        
        //std::cout<<"test1.2"<<std::endl;
        //TEST WHEN COMPARING TO DATA, WE CAN'T REQUIRE GENERATOR LEVEL INFORMATION
        m_histogram["h_GenJetPt"]->Fill(genJets.at(genJetIndex).pt(),weight);
        //FIXME should I also require the eta and pt cuts in gen jets??
        
        //find if the jet with the smallest R is associated to a hadron
        std::vector<int> numHadMatched;
        
        for (size_t i_had=0;i_had!=bHadJetIndex.size();++i_had)
        {
            //std::cout<<"test2"<<std::endl;
            //FIXME There are cases in which the same gen jet is matched to two different reco jets... what do I do with this?
            if (genJetIndex!=bHadJetIndex.at(i_had)) continue;
            numHadMatched.push_back(1);
            flavour.push_back(bHadFlavour.at(i_had));
            hadFlav.push_back(bHadPlusMothersPdgId.at(bHadIndex.at(i_had))); //since we have "plus mothers" better use this indexing->see NTupleWriter.cc
            double jetChargeValue=jetChargeRelativePtWeighted.at(jetIdx);
            charge.push_back(jetChargeValue);
            isMatched=1;
            
            //only continue if the jet is a true b-jet
            if (isMatched==-1) continue;
            //std::cout<<"test3"<<std::endl;
            //TEST Here test the genLepton new features. 
            //for each jet check if there's a lepton associated to a hadron associated to this jet. 
            for (size_t i_genLep = 0;i_genLep!=genBHadLeptonIndex.size();i_genLep++)
            {
                //std::cout<<"test4"<<std::endl;
                //if the hadron is not matched to a jet, we skip checking if it has a lepton associated to it 
                if (isMatched==-1) continue;
                //now we check if the index of the hadron is the same of the one of the jet
                if (genBHadLeptonHadronIndex.at(i_genLep) == -1) continue;
                //we only take leptons matched to the hadron
                if(genBHadLeptonHadronIndex.at(i_genLep) != (int)i_had) continue;
                //isGenLepton = true;
                genLeptonFlavour.push_back(bHadPlusMothersPdgId.at(genBHadLeptonIndex.at(i_genLep)));
                genLeptonPt.push_back(bHadPlusMothers.at(genBHadLeptonIndex.at(i_genLep)).pt());
                if (std::abs(bHadPlusMothersPdgId.at(genBHadLeptonIndex.at(i_genLep))) == 13) 
                {
                    genMuonPt.push_back(bHadPlusMothers.at(genBHadLeptonIndex.at(i_genLep)).pt());
                    genMuonLV.push_back(bHadPlusMothers.at(genBHadLeptonIndex.at(i_genLep)));
                }
                if (std::abs(bHadPlusMothersPdgId.at(genBHadLeptonIndex.at(i_genLep))) == 13) isGenMuon = true;
                if (std::abs(bHadPlusMothersPdgId.at(genBHadLeptonIndex.at(i_genLep))) == 11) isGenElectron = true;
                //std::cout<<"Per hadron we have such leptons = "<<genLeptonFlavour.size()<<std::endl;
            }
        }
        
        //choose the lepton with highest pt and take its flavour. Needed to compare charge and pt with the highest pt muon track.
        if (genLeptonPt.size()==1&&genLeptonFlavour.at(0)==-13) isGenPlusMuon = true;
        if (genLeptonPt.size()==1&&genLeptonFlavour.at(0)==13) isGenMinMuon = true;
        
        if (genLeptonPt.size()==1&&genLeptonFlavour.at(0)==-11) isGenPlusElectron = true;
        if (genLeptonPt.size()==1&&genLeptonFlavour.at(0)==11) isGenMinElectron = true;
        
        else if (genLeptonPt.size()>1)
        {
            double maxPtLepton = -999.;
            double maxPtLeptonFlavour = -1;
            for (size_t i_lepPt=0; i_lepPt!=genLeptonPt.size();i_lepPt++)
            {
                if (genLeptonPt.at(i_lepPt)>maxPtLepton) maxPtLepton = genLeptonPt.at(i_lepPt);
                maxPtLeptonFlavour = genLeptonFlavour.at(i_lepPt);
            }
            
            if (maxPtLeptonFlavour == -13) isGenPlusMuon = true;
            else if (maxPtLeptonFlavour == 13) isGenMinMuon = true;
            
            if (maxPtLeptonFlavour == -11) isGenPlusElectron = true;
            else if (maxPtLeptonFlavour == 11) isGenMinElectron = true;
            
        }
        
        //std::cout<<"test5"<<std::endl;
        
        if (debug_==2) std::cout<<"number of hadrons per jet is = "<<numHadMatched.size()<<std::endl;
        
        //if the jet is matched to a hadron, fill the true histograms
        double trueBJetScalarCharge = jetChargeRelativePtWeighted.at(jetIdx);
        double trueBJetPt = jets.pt();
        double etaTrueBJets = jets.eta();
        double jetTrueBPx = jets.px();
        double jetTrueBPy = jets.py();
        double jetTrueBPz = jets.pz();
        
        int trueBJetTrackMultiplicity = 0;
        int trueBJetSelTrackMultiplicity = 0;
        
        std::vector<double>sumTrueBMomentum;
        std::vector<double>sumTrueBMomentumQ;
        std::vector<double>sumTrueBMagnitude;
        std::vector<double>sumTrueBMagnitudeQ;
        
        for (size_t i_sumIni = 0;i_sumIni!=10;i_sumIni++)
        {
            sumTrueBMomentum.push_back(0);
            sumTrueBMomentumQ.push_back(0);
            sumTrueBMagnitude.push_back(0);
            sumTrueBMagnitudeQ.push_back(0);
        }
        
        double maxPtTrueTrack  = -999.;
        int maxChargeTrueTrack = -999;
        double maxPtChargeTrueTrack = -999.;
        double maxTrueMagnitude = -999.;
        double maxTrueProduct = -999.;
        
        //iP-variables
        double trueBJetTrackMaxIPSignificance = -1.; //only take positive values and discard any negative
        double trueBJetTrackMaxIPValue = -1.;
        double trueBJetIPSigForMaxPtTrack = -999.;
        double trueBJetPtForMaxIPSigTrack = -999.;
        std::vector<double> trueBJetTrackIPValueRank;
        std::vector<double> trueBJetTrackIPSignificanceRank;
        
        //lepton-tracks variables
        std::vector<double> trueBJetLeptonTracksPt;
        std::vector<LV> trueBJetLeptonTracksLV;
        std::vector<double> trueBJetLeptonTracksCharge;
        bool isB = false;
        bool isSecB = false;
        std::vector<int> leptonTrackPdgId;
        std::vector<double> ptRatioValues;
        double ptRatioValuesMuon = -999.;
        bool isLeadingMuon = false;
        bool isLeadingElectron = false;
        bool isNonLeadingLepton = false;
        bool isSecLeadingLepton = false;
        
        //std::cout<<"test6"<<std::endl;
        //====================================================================================WE ENTER TRACK LOOP!!=======================================
        for (size_t i_pfTrack=0; i_pfTrack!=jetPfCandidateTrack.size();i_pfTrack++)
        {
            //std::cout<<"test7"<<std::endl;
            //access the pt ordered pftrack collection
            int pfTrackIdx = ptOrderedJetTrackIdx.at(i_pfTrack);
            
            //access the track-LV and put some easy definitions to help.
            LV* pfTracks = new LV();
            *pfTracks = jetPfCandidateTrack.at(pfTrackIdx);
            int pfTrackJetIdx = jetPfCandidateTrackIndex.at(pfTrackIdx);
            double trueBJetPfTrackPt = pfTracks->pt();
            
            //here we require only tracks above a certain threshold of pt
            if (trueBJetPfTrackPt < ptTracksCUT) continue;
            //std::cout<<"test8"<<std::endl;
            
            //here check if the track is matched to a selected jet and in case it is, add one to the multiplicity.
            int trueMatched = 0;
            if (jetIdx!=pfTrackJetIdx) trueMatched = -1;
            if (debug_==2) std::cout<<"jetIdx = "<<jetIdx<<" and "<<"jetPfCandidateTrackIndex at pfTrackIdx = "<<jetPfCandidateTrackIndex.at(pfTrackIdx)<<std::endl;
            if (trueMatched == -1) continue;
            if (debug_==2)  std::cout<<"tracks belong to the jet = "<<jetPfCandidateTrackIndex.at(pfTrackIdx)<<std::endl;  
            trueBJetTrackMultiplicity++;
            //std::cout<<"test9"<<std::endl;
            //double trueBJetPfTrackPx = pfTracks->px();
            //double trueBJetPfTrackPy = pfTracks->py();
            //double trueBJetPfTrackPz = pfTracks->pz();
            
            //check if the pftrack is matched to a sel track. If it is, store it's information. Otherwise, continue to the next track. 
            bool pfMatchedToSelTrack = false;
            double trueBJetTrackIPValue = -999.;
            double trueBJetTrackIPSignificance = -999.;
            std::vector<int> checkHowManySelTracksToOnePfTrack;
            m_histogram["h_checkHowManyTracks"]->Fill(trueBJetPfTrackPt,weight);
            
            for (size_t i_jetSelectedTrack=0; i_jetSelectedTrack!=ptOrderedSelTrackIdx.size();i_jetSelectedTrack++)
            {
                //std::cout<<"test10"<<std::endl;
                int jetSelectedTrackIdx = ptOrderedSelTrackIdx.at(i_jetSelectedTrack);
                
                //access the track-LV and putting some easy definitions to help.
                LV* jetSelectedTracks = new LV();
                *jetSelectedTracks = jetSelectedTrack.at(jetSelectedTrackIdx);
                int jetSelectedTrackJetIdx = jetSelectedTrackIndex.at(jetSelectedTrackIdx);
                
                //double trueBJetSelTrackPx = jetSelectedTracks->px();
                //double trueBJetSelTrackPy = jetSelectedTracks->py();
                //double trueBJetSelTrackPz = jetSelectedTracks->pz();
                
                //match jetSelectedTracks to pfTracks using px, py, pz, charge and jetId
                if (jetSelectedTrackJetIdx!=pfTrackJetIdx) continue;
                //std::cout<<"test11"<<std::endl;
                //if( std::abs(trueBJetPfTrackPx-trueBJetSelTrackPx)>1.e-6 || std::abs(trueBJetPfTrackPy-trueBJetSelTrackPy)>1.e-6 || std::abs(trueBJetPfTrackPz-trueBJetSelTrackPz)>1.e-6 || jetSelectedTrackCharge.at(jetSelectedTrackIdx)!=jetPfCandidateTrackCharge.at(pfTrackIdx)) continue;
                if(std::abs(pfTracks->eta()-jetSelectedTracks->eta())>1.e-6 || std::abs(pfTracks->phi()-jetSelectedTracks->phi())>1.e-6 || std::abs(pfTracks->pt()-jetSelectedTracks->pt())>1.e-6 || jetSelectedTrackCharge.at(jetSelectedTrackIdx)!=jetPfCandidateTrackCharge.at(pfTrackIdx)) continue;
                if (debug_==3) std::cout<<"Pf-track with index "<<pfTrackIdx<<" is matched to sel-track with index "<<jetSelectedTrackIdx<<std::endl;
                //std::cout<<"test12"<<std::endl;
                
                //boolean to determine if the sel-tracks are matched to a pf-track or not
                pfMatchedToSelTrack = true;
                checkHowManySelTracksToOnePfTrack.push_back(1);
               
                //Access the impact parameter
                trueBJetTrackIPValue = jetSelectedTrackIPValue.at(jetSelectedTrackIdx);
                trueBJetTrackIPSignificance = jetSelectedTrackIPSignificance.at(jetSelectedTrackIdx);
                trueBJetTrackIPValueRank.push_back(trueBJetTrackIPValue);
                trueBJetTrackIPSignificanceRank.push_back(trueBJetTrackIPSignificance);
                if (debug_==3) std::cout<<"The impact parameter of this track is ="<<trueBJetTrackIPValue<<std::endl;
                
                if (trueBJetTrackIPValue>trueBJetTrackMaxIPValue)
                {
                    trueBJetTrackMaxIPValue = trueBJetTrackIPValue;
                }
                if (trueBJetTrackIPSignificance>trueBJetTrackMaxIPSignificance)
                {
                    trueBJetTrackMaxIPSignificance = trueBJetTrackIPSignificance;
                    trueBJetPtForMaxIPSigTrack = jetSelectedTracks->pt();
                }
                
               delete  jetSelectedTracks;
                
            }
            
            //check if there's more than one sel-track per pf-track. THIS SHOULD NOT HAPPEN!
            if (checkHowManySelTracksToOnePfTrack.size()>1) std::cout<<"Too many matches!!"<<std::endl; 
            
            //countinue ONLY if the pfTrack was checked to a jetSelectedTrack. Otherwise, no calculation is performed.
            if (pfMatchedToSelTrack==false) continue;
            //std::cout<<"test13"<<std::endl;
            if (debug_==3) std::cout<<"Sanity check: This should only be printed when the pf track is matched to a selected track"<<std::endl;
            
            
            //access the impact parameter
            m_histogram["h_trueBJetTrackIPValue"]->Fill(trueBJetTrackIPValue, weight);
            m_histogram["h_trueBJetTrackIPSignificance"]->Fill(trueBJetTrackIPSignificance, weight);
            ((TH2D*)m_histogram["h_trueBJetTrackIPValueVsSignificance"])->Fill(trueBJetTrackIPValue,trueBJetTrackIPSignificance,weight);
            ((TH2D*)m_histogram["h_trueBJetTrackPtVsIPSignificance"])->Fill(trueBJetPfTrackPt,trueBJetTrackIPSignificance,weight);
            
            //only fill the pt of the tracks that fulfill all the conditions
            if(trueBJetPfTrackPt>maxPtTrueTrack) 
            {
                maxPtTrueTrack = trueBJetPfTrackPt;
                maxChargeTrueTrack = jetPfCandidateTrackCharge.at(pfTrackIdx);
                maxPtChargeTrueTrack = trueBJetPfTrackPt*jetPfCandidateTrackCharge.at(pfTrackIdx);
                trueBJetIPSigForMaxPtTrack = trueBJetTrackIPSignificance;
            }
            //std::cout<<"test14"<<std::endl;
           
            //check if the track is a lepton or not:
            int particleId = jetPfCandidateTrackId.at(pfTrackIdx);
            leptonTrackPdgId.push_back(particleId);
            
            if (debug_==2) std::cout<<"Track pdgId is = "<<particleId<<" for the track "<<pfTrackIdx<<std::endl;
            //std::cout<<"test15"<<std::endl;
            
            //if track is a muon (3), fill some muon-specific histograms
            if (particleId==3) 
            {
                trueBJetLeptonTracksPt.push_back(trueBJetPfTrackPt);
                trueBJetLeptonTracksLV.push_back(*pfTracks);
                trueBJetLeptonTracksCharge.push_back(jetPfCandidateTrackCharge.at(pfTrackIdx));
                m_histogram["h_trueBJetLeptonTrackPt"]->Fill(trueBJetPfTrackPt, weight);
                m_histogram["h_trueBJetLeptonTrackEta"]->Fill(pfTracks->Eta(), weight);
                m_histogram["h_trueBJetLeptonTrackCharge"]->Fill(jetPfCandidateTrackCharge.at(pfTrackIdx), weight);
                ((TH2D*)m_histogram["h_trueBJetLeptonChargePt"])->Fill(trueBJetPfTrackPt,jetPfCandidateTrackCharge.at(pfTrackIdx),weight);
                m_histogram["h_trueBJetToLeptonTrackPtRatio"]->Fill(trueBJetPfTrackPt/trueBJetPt, weight);
                m_histogram["h_trueBJetTrackIPValueIfLepton"]->Fill(trueBJetTrackIPValue, weight);
                m_histogram["h_trueBJetTrackIPSignificanceIfLepton"]->Fill(trueBJetTrackIPSignificance, weight);
                ((TH2D*)m_histogram["h_trueBJetTrackIPValueVsSignificanceIfLepton"])->Fill(trueBJetTrackIPValue,trueBJetTrackIPSignificance, weight);
                ((TH2D*)m_histogram["h_trueBJetTrackPtVsIPSignificanceIfLepton"])->Fill(trueBJetPfTrackPt,trueBJetTrackIPSignificance, weight);
                isMuon = true;
            }
            if (particleId==2) isElectron = true;
            //std::cout<<"test16"<<std::endl;
            
            // We set some boolean to later on identify if the b-quark has the same charge as the lepton with highest pt in the jet.
            if (trueBJetLeptonTracksCharge.size()>0 && trueBJetLeptonTracksCharge.at(0) <0 ) isB=true;
            if (trueBJetLeptonTracksCharge.size()>1 && trueBJetLeptonTracksCharge.at(1) <0 ) isSecB=true;
            //std::cout<<"test17"<<std::endl;
            
            //pt of track to pt of jet relation histograms
            double ptRatio = trueBJetPfTrackPt/trueBJetPt;
            ptRatioValues.push_back(ptRatio);
            //std::cout<<"test18"<<std::endl;
            
            if (ptRatioValues.size()==1&&particleId==2) 
            {
                ptRatioValuesMuon = ptRatio;
                isLeadingElectron = true;
            }
            
            if (ptRatioValues.size()==1&&particleId==3) 
            {
                ptRatioValuesMuon = ptRatio;
                isLeadingMuon = true;
            }
            else if (ptRatioValues.size()>1&&particleId==3&&isLeadingMuon==false)  isNonLeadingLepton = true;
            else if (ptRatioValues.size()>1&&particleId==3&&isLeadingMuon)  isSecLeadingLepton = true;
               
            //std::cout<<"test19"<<std::endl;
            
            m_histogram["h_trueBJetToTrackPtRatio"]->Fill(ptRatio, weight);
            
            m_histogram["h_trueBJetLeptonTracks"]->Fill(particleId);
            
            m_histogram["h_trueBJetPfTrackPt"]->Fill(trueBJetPfTrackPt, weight);
            
            //calculate the jet charge in the same way it is done at ntuple level to cross-check the information 
            //obtain from the trees is the same obtained directly from the RECO objects
            const double constituentTrueBPx = pfTracks->px();
            const double constituentTrueBPy = pfTracks->py();
            const double constituentTrueBPz = pfTracks->pz();
            const double trueProduct = constituentTrueBPx*jetTrueBPx + constituentTrueBPy*jetTrueBPy + constituentTrueBPz*jetTrueBPz;
            //std::cout<<"test20"<<std::endl;
            
            std::vector<double> vectProductMomentumQ;
            double xTrueComponent = (jetTrueBPy*constituentTrueBPz-jetTrueBPz*constituentTrueBPy);
            double yTrueComponent = (jetTrueBPx*constituentTrueBPz-jetTrueBPz*constituentTrueBPx);
            double zTrueComponent = (jetTrueBPx*constituentTrueBPy-jetTrueBPy*constituentTrueBPx); 
            const double trueMagnitude = std::sqrt(xTrueComponent*xTrueComponent+yTrueComponent*yTrueComponent+zTrueComponent*zTrueComponent);
            //std::cout<<"test21"<<std::endl;
            
            std::vector<double> trueProductPow;
            std::vector<double> trueMagnitudePow;
            
            for (double i_pow=0.2;i_pow<=2.;i_pow+=0.2)
            {
                trueProductPow.push_back(std::pow(trueProduct,i_pow));
                trueMagnitudePow.push_back(std::pow(trueMagnitude,i_pow));
            }
            //std::cout<<"test22"<<std::endl;
            
            for (size_t i_sum=0;i_sum!=sumTrueBMomentum.size();i_sum++)
            {
                sumTrueBMomentum.at(i_sum) += trueProductPow.at(i_sum);
                sumTrueBMagnitude.at(i_sum) += trueMagnitudePow.at(i_sum);
                sumTrueBMomentumQ.at(i_sum) += (trueProductPow.at(i_sum))*jetPfCandidateTrackCharge.at(pfTrackIdx);
                sumTrueBMagnitudeQ.at(i_sum) += (trueMagnitudePow.at(i_sum))*jetPfCandidateTrackCharge.at(pfTrackIdx);
            }
            
            //std::cout<<"test23"<<std::endl;
            if (trueProduct>maxTrueProduct) maxTrueProduct = trueProduct;
            if (trueMagnitude>maxTrueMagnitude) maxTrueMagnitude = trueMagnitude;
            
            m_histogram["h_trueBJetRelPtTrack"]->Fill(trueMagnitude, weight);
            
            if (debug_==3) std::cout<<"------------------------"<<std::endl;
            
            delete pfTracks;
            
        } //end of sel track loop
        
        //check sel-pf matching -> If greater than one, then some selected track was not matched
        m_histogram["h_trueBJetSelTrackMultiplicity"]->Fill(trueBJetSelTrackMultiplicity, weight);
        ((TH2D*)m_histogram["h_pfSizeVsSelSizeNoMatching"])->Fill(trueBJetSelTrackMultiplicity,trueBJetTrackMultiplicity,weight);
        //((TH2D*) m_histogram["h_trueBJetLeptonTrackPtMultiplicity"])->Fill(trueBJetPt,trueBJetTrackMultiplicity,weight);
         //std::cout<<"test24"<<std::endl;
        
        //Which particle leads our tracks?
        if(leptonTrackPdgId.size()>=1) m_histogram["h_trueBJetLeptonLeadingTrackPdgId"] -> Fill(leptonTrackPdgId.at(0));
        if(leptonTrackPdgId.size()>=2) m_histogram["h_trueBJetLeptonSubLeadingTrackPdgId"] -> Fill(leptonTrackPdgId.at(1));
        if(leptonTrackPdgId.size()>=3) m_histogram["h_trueBJetLeptonThirdLeadingTrackPdgId"] -> Fill(leptonTrackPdgId.at(2));
        if(leptonTrackPdgId.size()>=4) m_histogram["h_trueBJetLeptonFourthLeadingTrackPdgId"] -> Fill(leptonTrackPdgId.at(3));
        
        //std::cout<<"test25"<<std::endl;
        //Where are our muons in the track list? 
        if (leptonTrackPdgId.size()>0) 
        {
            for (size_t i_lep = 0;i_lep!=leptonTrackPdgId.size();i_lep++)
            {
                if (leptonTrackPdgId.at(i_lep)==3) 
                {
                    m_histogram["h_trueBJetLeptonRankInTrack"]->Fill(i_lep);
                }
            }
        }
        
        //std::cout<<"test26"<<std::endl;
        if (ptRatioValues.size()!=0) m_histogram["h_trueBJetToLeadingTrackPtRatio"]->Fill(ptRatioValues.at(0),weight);
        //std::cout<<"test27"<<std::endl;
        if (isLeadingMuon==true) m_histogram["h_trueBJetToLeadingLeptonTrackPtRatio"]->Fill(ptRatioValuesMuon, weight);
        //std::cout<<"test28"<<std::endl;
        
        //lepton histograms
        if (trueBJetLeptonTracksPt.size()>0) 
        {
            m_histogram["h_trueBJetLeptonTrackMultiplicity"] -> Fill(trueBJetLeptonTracksPt.size(), weight);
            ((TH2D*) m_histogram["h_trueBJetLeptonTrackPtMultiplicity"])->Fill(trueBJetPt,trueBJetTrackMultiplicity,weight);
            m_histogram["h_trueBJetTrackMultiplicityIfLepton"]->Fill(trueBJetTrackMultiplicity, weight);
        }
        
        //genLepton to recoLepton comparison plots
        if (isGenMuon==false && isMuon==false) m_histogram["h_trueBJetAgreementMuons"]->Fill(0);
        else if (isGenMuon==true && isMuon==false) m_histogram["h_trueBJetAgreementMuons"]->Fill(1);
        else if (isGenMuon==false && isMuon==true) m_histogram["h_trueBJetAgreementMuons"]->Fill(2);
        else if (isGenMuon==true && isMuon==true) m_histogram["h_trueBJetAgreementMuons"]->Fill(3);
        
        if (isGenElectron==false && isElectron==false) m_histogram["h_trueBJetAgreementElectrons"]->Fill(0);
        else if (isGenElectron==true && isElectron==false) m_histogram["h_trueBJetAgreementElectrons"]->Fill(1);
        else if (isGenElectron==false && isElectron==true) m_histogram["h_trueBJetAgreementElectrons"]->Fill(2);
        else if (isGenElectron==true && isElectron==true) m_histogram["h_trueBJetAgreementElectrons"]->Fill(3);
        
        if (isLeadingMuon&&isB&&isMuon&&isGenMinMuon) m_histogram["h_trueBJetAgreementChargeMuons"]->Fill(0);
        if (isLeadingMuon&&isB==false&&isMuon&&isGenPlusMuon) m_histogram["h_trueBJetAgreementChargeMuons"]->Fill(1);
        if (isLeadingMuon&&isB&&isMuon&&isGenPlusMuon) m_histogram["h_trueBJetAgreementChargeMuons"]->Fill(2);
        if (isLeadingMuon&&isB==false&&isMuon&&isGenMinMuon) m_histogram["h_trueBJetAgreementChargeMuons"]->Fill(3);
        
        if (isLeadingElectron&&isB&&isElectron&&isGenMinElectron) m_histogram["h_trueBJetAgreementChargeElectrons"]->Fill(0);
        if (isLeadingElectron&&isB==false&&isElectron&&isGenPlusElectron) m_histogram["h_trueBJetAgreementChargeElectrons"]->Fill(1);
        if (isLeadingElectron&&isB&&isElectron&&isGenPlusElectron) m_histogram["h_trueBJetAgreementChargeElectrons"]->Fill(2);
        if (isLeadingElectron&&isB==false&&isElectron&&isGenMinElectron) m_histogram["h_trueBJetAgreementChargeElectrons"]->Fill(3);
        
        //find mean distances on pt between all reco and gen leptons
        if (genMuonPt.size()>0 && trueBJetLeptonTracksPt.size()>0)
        {
            
            for (size_t i_recoLep = 0; i_recoLep!=trueBJetLeptonTracksPt.size(); i_recoLep++)
            {
                std::cout<<"The reco lep has a pt = "<<trueBJetLeptonTracksPt.at(i_recoLep)<<std::endl;
                std::cout<<"The reco lep has a phi = "<<trueBJetLeptonTracksLV.at(i_recoLep).phi()<<std::endl;
                std::cout<<"The reco lep has an eta = "<<trueBJetLeptonTracksLV.at(i_recoLep).eta()<<std::endl;
                for (size_t i_genLep = 0; i_genLep!=genMuonPt.size(); i_genLep++)
                {
                    std::cout<<"The gen lep has a pt = "<<genMuonPt.at(i_genLep)<<std::endl;
                    std::cout<<"The gen lep has a phi = "<<genMuonLV.at(i_genLep).phi()<<std::endl;
                    std::cout<<"The gen lep has an eta = "<<genMuonLV.at(i_genLep).eta()<<std::endl;
                    //pt difference between gen and reco muons
                    m_histogram["h_trueBJetGenToRecoPtDifferenceMuons"]->Fill(std::abs(trueBJetLeptonTracksPt.at(i_recoLep)-genMuonPt.at(i_genLep)));
                    m_histogram["h_trueBJetGenToRecoPtDifferenceMuonsExtended"]->Fill(std::abs(trueBJetLeptonTracksPt.at(i_recoLep)-genMuonPt.at(i_genLep)));
                    std::cout<<"the deltaPt is = "<<std::abs(trueBJetLeptonTracksPt.at(i_recoLep)-genMuonPt.at(i_genLep))<<std::endl;
                    
                    //phi difference between reco and gen muons
                    m_histogram["h_trueBJetGenToRecoPhiDifferenceMuons"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(trueBJetLeptonTracksLV.at(i_recoLep), genMuonLV.at(i_genLep)));
                    std::cout<<"the deltaPhi is = "<<ROOT::Math::VectorUtil::DeltaPhi(trueBJetLeptonTracksLV.at(i_recoLep), genMuonLV.at(i_genLep))<<std::endl;;
                    
                    //eta difference between reco and gen muons
                    m_histogram["h_trueBJetGenToRecoEtaDifferenceMuons"]->Fill(std::abs(trueBJetLeptonTracksLV.at(i_recoLep).eta()-genMuonLV.at(i_genLep).eta()));
                    std::cout<<"the deltaEta is = "<<std::abs(trueBJetLeptonTracksLV.at(i_recoLep).eta()-genMuonLV.at(i_genLep).eta())<<std::endl;
                    
                    //deltaR between reco and gen muons
                    m_histogram["h_trueBJetGenToRecoDeltaRMuons"]->Fill(ROOT::Math::VectorUtil::DeltaR(trueBJetLeptonTracksLV.at(i_recoLep), genMuonLV.at(i_genLep)));
                    std::cout<<"the deltaR is = "<<ROOT::Math::VectorUtil::DeltaR(trueBJetLeptonTracksLV.at(i_recoLep), genMuonLV.at(i_genLep))<<std::endl;
                }
            }
            std::cout<<"==========================================================="<<std::endl;
        }
        
        
        //IP histograms (only filled if there was a matching to a selected track!
        if (trueBJetTrackMaxIPValue>=0.) m_histogram["h_trueBJetTrackMaxIPValue"]->Fill(trueBJetTrackMaxIPValue, weight);
        if (trueBJetTrackMaxIPSignificance>=0.) m_histogram["h_trueBJetTrackMaxIPSignificance"]->Fill(trueBJetTrackMaxIPSignificance, weight);
        if (trueBJetTrackMaxIPValue>=0.&&trueBJetTrackMaxIPSignificance>=0.) ((TH2D*) m_histogram["h_trueBJetTrackPtVsHighestIPSignificance"])->Fill(trueBJetPtForMaxIPSigTrack,trueBJetTrackMaxIPSignificance,weight);
        if (trueBJetTrackMaxIPSignificance>=0.) ((TH2D*) m_histogram["h_trueBJetTrackHighestPtVsIPSignificance"])->Fill(maxPtTrueTrack, trueBJetIPSigForMaxPtTrack, weight);
        
        for (size_t i_val=0;i_val!=trueBJetTrackIPValueRank.size();i_val++)
        {
            if (trueBJetTrackIPValueRank.at(i_val) == trueBJetTrackMaxIPValue) m_histogram["h_trueBJetTrackMaxIPValueRank"]->Fill(i_val, weight);
            if (trueBJetTrackIPSignificanceRank.at(i_val) == trueBJetTrackMaxIPSignificance) m_histogram["h_trueBJetTrackMaxIPSignificanceRank"]->Fill(i_val,weight);
        }
        
        //track general histograms
        m_histogram["h_trueBJetMaxRelPtTrack"]->Fill(maxTrueMagnitude, weight);
        m_histogram["h_trueBJetHighestPtTrack"]->Fill(maxPtTrueTrack, weight);
        
        ((TH2D*)m_histogram["h_trueBJetPtVsMaxPtTrack"])->Fill(maxPtTrueTrack,trueBJetPt,weight);
        ((TH2D*)m_histogram["h_trueBJetPtVsMaxRelPtTrack"])->Fill(maxTrueMagnitude,trueBJetPt,weight);
        ((TH2D*)m_histogram["h_trueBJetPtVsMaxScalarPtTrack"])->Fill(maxTrueProduct,trueBJetPt,weight);
        ((TH2D*)m_histogram["h_trueBJetPtVsNumTracks"])-> Fill(trueBJetTrackMultiplicity,trueBJetPt,weight);
        ((TH2D*)m_histogram["h_trueBJetMaxPtTrackVsNumTracks"])-> Fill(trueBJetTrackMultiplicity,maxPtTrueTrack,weight);
        ((TH2D*)m_histogram["h_trueBJetEtaVsTrackMultiplicity"])-> Fill(trueBJetTrackMultiplicity,etaTrueBJets,weight);
        m_histogram["h_trueBJetPt"]->Fill(trueBJetPt, weight);
        
        if (debug_==2) std::cout<<"The trueBJetTrackMultiplicity is = "<<trueBJetTrackMultiplicity<<std::endl;
        m_histogram["h_trueBJetTrackMultiplicity"]->Fill(trueBJetTrackMultiplicity,weight);
       
        //charge histograms
        m_histogram["h_trueBJetScalarChargeValidation"]->Fill(trueBJetScalarCharge, weight);
         //std::cout<<"test29"<<std::endl;
        
        const double trueBJetScalarCharge10(sumTrueBMomentum.at(4)>0 ? sumTrueBMomentumQ.at(4)/sumTrueBMomentum.at(4) : 0);
        ((TH2D*)m_histogram["h_trueBJetScalarChargeVsMultip"])->Fill(trueBJetTrackMultiplicity,trueBJetScalarCharge10,weight);
        //std::cout<<"test30"<<std::endl;
        
        
        //check if lepton track and charge of the jet coincide
        if (trueBJetLeptonTracksPt.size()>0) 
            //std::cout<<"test31"<<std::endl;
        
        {
            if(isB&&isLeadingMuon&&trueBJetScalarCharge10<0)
            {
                m_histogram["h_trueBJetLeadingLeptonScalarChargeMatch"]->Fill(1);
                m_histogram["h_trueBJetLeadingLeptonRelChargeMatch"]->Fill(1);
            }
            else if(isB==false&&isLeadingMuon&&trueBJetScalarCharge10>0) 
            {
                m_histogram["h_trueBJetLeadingLeptonScalarChargeMatch"]->Fill(1);
                m_histogram["h_trueBJetLeadingLeptonRelChargeMatch"]->Fill(1);
                
            }
            else if (isB &&isLeadingMuon&& trueBJetScalarCharge10>0)
            {
                m_histogram["h_trueBJetLeadingLeptonScalarChargeMatch"]->Fill(0);
                m_histogram["h_trueBJetLeadingLeptonRelChargeMatch"]->Fill(0);
            }
            else if (isB==false&&isLeadingMuon&&trueBJetScalarCharge10<0) 
            {
                m_histogram["h_trueBJetLeadingLeptonScalarChargeMatch"]->Fill(0);
                m_histogram["h_trueBJetLeadingLeptonRelChargeMatch"]->Fill(0);
            }
            else if(isB&&isNonLeadingLepton&&trueBJetScalarCharge10<0)
            {
                m_histogram["h_trueBJetNonLeadingLeptonScalarChargeMatch"]->Fill(1);
                m_histogram["h_trueBJetNonLeadingLeptonRelChargeMatch"]->Fill(1);
            }
            else if(isB==false&&isNonLeadingLepton&&trueBJetScalarCharge10>0) 
            {
                m_histogram["h_trueBJetNonLeadingLeptonScalarChargeMatch"]->Fill(1);
                m_histogram["h_trueBJetNonLeadingLeptonRelChargeMatch"]->Fill(1);
            }
            else if (isB &&isNonLeadingLepton&& trueBJetScalarCharge10>0)
            {
                m_histogram["h_trueBJetNonLeadingLeptonScalarChargeMatch"]->Fill(0);
                m_histogram["h_trueBJetNonLeadingLeptonRelChargeMatch"]->Fill(0);
            }
            else if (isB==false&&isNonLeadingLepton&&trueBJetScalarCharge10<0) 
            {
                m_histogram["h_trueBJetNonLeadingLeptonScalarChargeMatch"]->Fill(0);
                m_histogram["h_trueBJetNonLeadingLeptonRelChargeMatch"]->Fill(0);
            }
            else if(isSecB&&isSecLeadingLepton&&trueBJetScalarCharge10<0)
            {
                m_histogram["h_trueBJetSecLeadingLeptonScalarChargeMatch"]->Fill(1);
                m_histogram["h_trueBJetSecLeadingLeptonRelChargeMatch"]->Fill(1);
            }
            else if(isSecB==false&&isSecLeadingLepton&&trueBJetScalarCharge10>0) 
            {
                m_histogram["h_trueBJetSecLeadingLeptonScalarChargeMatch"]->Fill(1);
                m_histogram["h_trueBJetSecLeadingLeptonRelChargeMatch"]->Fill(1);
            }
            else if (isSecB &&isSecLeadingLepton&& trueBJetScalarCharge10>0)
            {
                m_histogram["h_trueBJetSecLeadingLeptonScalarChargeMatch"]->Fill(0);
                m_histogram["h_trueBJetSecLeadingLeptonRelChargeMatch"]->Fill(0);
            }
            else if (isSecB==false&&isSecLeadingLepton&&trueBJetScalarCharge10<0) 
            {
                m_histogram["h_trueBJetSecLeadingLeptonScalarChargeMatch"]->Fill(0);
                m_histogram["h_trueBJetSecLeadingLeptonRelChargeMatch"]->Fill(0);
            }
        }
        //std::cout<<"test32"<<std::endl;
        
        
        //create a vector of histograms: contains histograms for x=0.2 to x=2.0
        std::vector<TString> trueBJetHistogramScalarVector;
        std::vector<TString> trueBJetHistogramRelVector;
       
        for (double i_histo=2;i_histo<=20;i_histo+=2)
        {
            std::stringstream ss_histoScalar;
            std::stringstream ss_histoRel;
            
            ss_histoScalar<<"h_trueBJetScalarCharge"<<i_histo;
            ss_histoRel<<"h_trueBJetRelCharge"<<i_histo;
           
            TString histoScalar = ss_histoScalar.str();
            TString histoRel = ss_histoRel.str();
           
            trueBJetHistogramScalarVector.push_back(histoScalar);
            trueBJetHistogramRelVector.push_back(histoRel);
          
        }
        //std::cout<<"test33"<<std::endl;
        
        //fill the histograms with the corresponding value of the charge (we take values from the sumTrueBMagnitude and sumTrueBMomentum vectors)
        std::vector<double> trueBJetScalarChargeVector;
        std::vector<double> trueBJetRelChargeVector;
        
        for (size_t i_sumHis=0;i_sumHis!=sumTrueBMomentum.size();i_sumHis++)
        {
            const double trueBJetScalarCharge(sumTrueBMomentum.at(i_sumHis)>0 ? sumTrueBMomentumQ.at(i_sumHis)/sumTrueBMomentum.at(i_sumHis) : 0);
            const double trueBJetRelCharge(sumTrueBMagnitude.at(i_sumHis)>0 ? sumTrueBMagnitudeQ.at(i_sumHis)/sumTrueBMagnitude.at(i_sumHis) : 0);
            
            TString histoScalar = trueBJetHistogramScalarVector.at(i_sumHis);
            TString histoRel = trueBJetHistogramRelVector.at(i_sumHis);
            
            m_histogram[histoScalar]->Fill(trueBJetScalarCharge);
            m_histogram[histoRel]->Fill(trueBJetRelCharge);
            
            trueBJetScalarChargeVector.push_back(trueBJetScalarCharge);
            trueBJetRelChargeVector.push_back(trueBJetRelCharge);
        }
        //std::cout<<"test34"<<std::endl;
        
        
        //fill a vector with the charge value for the different x taking into account that a jet might have more than one hadron associated
        for (size_t i_hadCorr = 0; i_hadCorr!=numHadMatched.size();i_hadCorr++)
        {
            trueBJetScalarCharge2V.push_back(trueBJetScalarChargeVector.at(0));
            trueBJetScalarCharge4V.push_back(trueBJetScalarChargeVector.at(1));
            trueBJetScalarCharge6V.push_back(trueBJetScalarChargeVector.at(2));
            trueBJetScalarCharge8V.push_back(trueBJetScalarChargeVector.at(3));
            trueBJetScalarCharge10V.push_back(trueBJetScalarChargeVector.at(4));
            trueBJetScalarCharge12V.push_back(trueBJetScalarChargeVector.at(5));
            trueBJetScalarCharge14V.push_back(trueBJetScalarChargeVector.at(6));
            trueBJetScalarCharge16V.push_back(trueBJetScalarChargeVector.at(7));
            trueBJetScalarCharge18V.push_back(trueBJetScalarChargeVector.at(8));
            trueBJetScalarCharge20V.push_back(trueBJetScalarChargeVector.at(9));
            
            trueBJetRelCharge2V.push_back(trueBJetRelChargeVector.at(0));
            trueBJetRelCharge4V.push_back(trueBJetRelChargeVector.at(1));
            trueBJetRelCharge6V.push_back(trueBJetRelChargeVector.at(2));
            trueBJetRelCharge8V.push_back(trueBJetRelChargeVector.at(3));
            trueBJetRelCharge10V.push_back(trueBJetRelChargeVector.at(4));
            trueBJetRelCharge12V.push_back(trueBJetRelChargeVector.at(5));
            trueBJetRelCharge14V.push_back(trueBJetRelChargeVector.at(6));
            trueBJetRelCharge16V.push_back(trueBJetRelChargeVector.at(7));
            trueBJetRelCharge18V.push_back(trueBJetRelChargeVector.at(8));
            trueBJetRelCharge20V.push_back(trueBJetRelChargeVector.at(9));
            //std::cout<<"test35"<<std::endl;
            
            //for mva exclusive:
            //if (isMuon == false)  
            //{
                trueBJetMaxPtTrackV.push_back(maxPtTrueTrack);
                trueBJetMaxChargeTrackV.push_back(maxChargeTrueTrack);
                trueBJetTrackMultiplicityV.push_back(trueBJetTrackMultiplicity);
                trueBJetMaxPtChargeTrackV.push_back(maxPtChargeTrueTrack);
                trueBJetPtV.push_back(trueBJetPt);
                //FIXME This needs to be understood. But for the moment should work
                if (ptRatioValues.size()>0) trueBJetPtRatioV.push_back(ptRatioValues.at(0));
                else trueBJetPtRatioV.push_back(-999.);
                //if (isMuon) isMuonV.push_back(1);
                //else if (isMuon==false) isMuonV.push_back(0);
            //}   
                //std::cout<<"test36"<<std::endl;
                
            // else if (isMuon == true)
            //   {
            //        ////std::cout<<"test25"<<std::endl;
            //       trueBJetMaxPtTrackV.push_back(-999.);
            //       trueBJetMaxChargeTrackV.push_back(-999.);
            //       trueBJetTrackMultiplicityV.push_back(-999.);
            //       trueBJetMaxPtChargeTrackV.push_back(-999.);
            //       trueBJetPtV.push_back(-999.);
            //       if (ptRatioValues.size()>0) trueBJetPtRatioV.push_back(-999.);
            //       else trueBJetPtRatioV.push_back(-999.);
            //    }
        }
        if (debug_==3) std::cout<<"====================================================="<<std::endl;

    } //end loop over reco jets
    
   
    //if (charge.size() != trueBJetRelChargeV.size()) std::cout<<"Size of the new charge = "<<trueBJetRelChargeV.size()<<" and charge size = "<<charge.size()<<std::endl;
  
    for (size_t i_fill=0;i_fill!=charge.size();++i_fill) //loop for charge.size()=flavour.size()->one entry per hadron/bquark, not per jet
    {
        //std::cout<<"test37"<<std::endl;
        ((TH2D*)m_histogram["h_trueBJetBQuarkCharge"])->Fill(charge.at(i_fill),flavour.at(i_fill), weight);
      
        //FIXME this would only work for Pythia. Herwig contains flavour oscillations!! Check from which hadron it's coming to take into account hadron oscillations
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkScalarCharge2"]->Fill(trueBJetScalarCharge2V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge2"]->Fill(trueBJetScalarCharge2V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkScalarCharge4"]->Fill(trueBJetScalarCharge4V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge4"]->Fill(trueBJetScalarCharge4V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkScalarCharge6"]->Fill(trueBJetScalarCharge6V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge6"]->Fill(trueBJetScalarCharge6V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkScalarCharge8"]->Fill(trueBJetScalarCharge8V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge8"]->Fill(trueBJetScalarCharge8V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkScalarCharge10"]->Fill(trueBJetScalarCharge10V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueBJetAntiBQuarkScalarCharge10"]->Fill(trueBJetScalarCharge10V.at(i_fill), weight);
      
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkScalarCharge12"]->Fill(trueBJetScalarCharge12V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge12"]->Fill(trueBJetScalarCharge12V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkScalarCharge14"]->Fill(trueBJetScalarCharge14V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge14"]->Fill(trueBJetScalarCharge14V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkScalarCharge16"]->Fill(trueBJetScalarCharge16V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge16"]->Fill(trueBJetScalarCharge16V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkScalarCharge18"]->Fill(trueBJetScalarCharge18V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge18"]->Fill(trueBJetScalarCharge18V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkScalarCharge20"]->Fill(trueBJetScalarCharge20V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkScalarCharge20"]->Fill(trueBJetScalarCharge20V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkRelCharge2"]->Fill(trueBJetRelCharge2V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge2"]->Fill(trueBJetRelCharge2V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkRelCharge4"]->Fill(trueBJetRelCharge4V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge4"]->Fill(trueBJetRelCharge4V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkRelCharge6"]->Fill(trueBJetRelCharge6V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge6"]->Fill(trueBJetRelCharge6V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkRelCharge8"]->Fill(trueBJetRelCharge8V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge8"]->Fill(trueBJetRelCharge8V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkRelCharge10"]->Fill(trueBJetRelCharge10V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueBJetAntiBQuarkRelCharge10"]->Fill(trueBJetRelCharge10V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkRelCharge12"]->Fill(trueBJetRelCharge12V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge12"]->Fill(trueBJetRelCharge12V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkRelCharge14"]->Fill(trueBJetRelCharge14V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge14"]->Fill(trueBJetRelCharge14V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkRelCharge16"]->Fill(trueBJetRelCharge16V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge16"]->Fill(trueBJetRelCharge16V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkRelCharge18"]->Fill(trueBJetRelCharge18V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge18"]->Fill(trueBJetRelCharge18V.at(i_fill),weight);
        
        if (flavour.at(i_fill)>0)  m_histogram["h_trueBJetBQuarkRelCharge20"]->Fill(trueBJetRelCharge20V.at(i_fill),weight);
        else if (flavour.at(i_fill)<0)  m_histogram["h_trueAntiBJetBQuarkRelCharge20"]->Fill(trueBJetRelCharge20V.at(i_fill),weight);
        
        //mva variables filling
        if (trueBJetPtV.at(i_fill)!=-999.)
        {
            //std::cout<<"test38"<<std::endl;
            if (flavour.at(i_fill)<0) mvaStruct.trueBJetId_ .push_back(-1);
            else if (flavour.at(i_fill)>0) mvaStruct.trueBJetId_.push_back(0);
            
            mvaStruct.relChargeJet_.push_back(trueBJetRelCharge8V.at(i_fill));
            mvaStruct.longChargeJet_.push_back(trueBJetScalarCharge8V.at(i_fill));
            mvaStruct.leadingTrackCharge_.push_back(trueBJetMaxChargeTrackV.at(i_fill));
            mvaStruct.leadingTrackPtCharge_.push_back(trueBJetMaxPtChargeTrackV.at(i_fill));
            mvaStruct.numTracks_.push_back(trueBJetTrackMultiplicityV.at(i_fill));
            mvaStruct.leadingTrackPt_.push_back(trueBJetMaxPtTrackV.at(i_fill));
            mvaStruct.trueBJetPt_.push_back(trueBJetPtV.at(i_fill));
            mvaStruct.ptRatioTrackJet_.push_back(trueBJetPtRatioV.at(i_fill));
            //mvaStruct.isMuonEvent_.push_back(isMuonV.at(i_fill));
            fillTree = true;
        }
        
        //std::cout<<"test39"<<std::endl;
        
        if (flavour.at(i_fill) == 6)  m_histogram["h_trueBJetTopCharge"]->Fill(trueBJetScalarCharge10V.at(i_fill),weight);
        else if (flavour.at(i_fill) == -6) m_histogram["h_trueBJetAntiTopCharge"]->Fill(trueBJetScalarCharge10V.at(i_fill),weight);
        else if (flavour.at(i_fill) == 25) m_histogram["h_trueBJetHiggsCharge"]->Fill(trueBJetScalarCharge10V.at(i_fill),weight);
        else if (flavour.at(i_fill) == -25) m_histogram["h_trueBJetAntiHiggsCharge"]->Fill(trueBJetScalarCharge10V.at(i_fill),weight);
        
        if (hadFlav.at(i_fill) == 521) m_histogram["h_trueBJetBPlusCharge"]->Fill(trueBJetScalarCharge10V.at(i_fill),weight);
        else if (hadFlav.at(i_fill) == -521) m_histogram["h_trueBJetBMinusCharge"]->Fill(trueBJetScalarCharge10V.at(i_fill),weight);
        else if (hadFlav.at(i_fill) == 511) m_histogram["h_trueBJetBZeroCharge"]->Fill(trueBJetScalarCharge10V.at(i_fill),weight);
        else if (hadFlav.at(i_fill) == -511) m_histogram["h_trueBJetAntiBZeroCharge"]->Fill(trueBJetScalarCharge10V.at(i_fill),weight);
        else m_histogram["checkBHadOthers"]->Fill(trueBJetScalarCharge10V.at(i_fill),weight);
        
        
        //hadron oscillations
        if (flavour.at(i_fill)<0 && hadFlav.at(i_fill) == 521) m_histogram["h_trueBJet_bbarToBPlus"]->Fill(1.,weight);
        if (flavour.at(i_fill)>0 && hadFlav.at(i_fill) == 521) m_histogram["h_trueBJet_bToBPlus"]->Fill(1.,weight);
        if (flavour.at(i_fill)<0 && hadFlav.at(i_fill) == 511) m_histogram["h_trueBJet_bbarToBZero"]->Fill(1.,weight);
        if (flavour.at(i_fill)>0 && hadFlav.at(i_fill) == 511) m_histogram["h_trueBJet_bToBZero"]->Fill(1.,weight);
        
        if (flavour.at(i_fill)>0 && hadFlav.at(i_fill) == -521) m_histogram["h_trueBJet_bToBMinus"]->Fill(1.,weight);
        if (flavour.at(i_fill)<0 && hadFlav.at(i_fill) == -521) m_histogram["h_trueBJet_bbarToBMinus"]->Fill(1.,weight);
        if (flavour.at(i_fill)>0 && hadFlav.at(i_fill) == -511) m_histogram["h_trueBJet_bToAntiBZero"]->Fill(1.,weight);
        if (flavour.at(i_fill)<0 && hadFlav.at(i_fill) == -511) m_histogram["h_trueBJet_bbarToAntiBZero"]->Fill(1.,weight);
        
        if (flavour.at(i_fill)<0 && hadFlav.at(i_fill) == 541) m_histogram["h_trueBJet_bbarToBcPlus"]->Fill(1.,weight);
        if (flavour.at(i_fill)>0 && hadFlav.at(i_fill) == 541) m_histogram["h_trueBJet_bToBcPlus"]->Fill(1.,weight);
        if (flavour.at(i_fill)<0 && hadFlav.at(i_fill) == 531) m_histogram["h_trueBJet_bbarToBsZero"]->Fill(1.,weight);
        if (flavour.at(i_fill)>0 && hadFlav.at(i_fill) == 531) m_histogram["h_trueBJet_bToBsZero"]->Fill(1.,weight);
        
        if (flavour.at(i_fill)>0 && hadFlav.at(i_fill) == -541) m_histogram["h_trueBJet_bToBcMinus"]->Fill(1.,weight);
        if (flavour.at(i_fill)<0 && hadFlav.at(i_fill) == -541) m_histogram["h_trueBJet_bbarToBcMinus"]->Fill(1.,weight);
        if (flavour.at(i_fill)>0 && hadFlav.at(i_fill) == -531) m_histogram["h_trueBJet_bToAntiBsZero"]->Fill(1.,weight);
        if (flavour.at(i_fill)<0 && hadFlav.at(i_fill) == -531) m_histogram["h_trueBJet_bbarToAntiBsZero"]->Fill(1.,weight);
    }
    //Fill the mva trees
    
    //std::cout<<"test40"<<std::endl;
    if (fillTree == true) 
    {
        if (nTrees == 0 ) mvaChargeTestTree->Fill();
        else mvaChargeTrainTree->Fill();
    }
    mvaStruct.relChargeJet_.clear();
    mvaStruct.longChargeJet_.clear();
    mvaStruct.leadingTrackCharge_.clear();
    mvaStruct.leadingTrackPtCharge_.clear();
    mvaStruct.leadingTrackPt_.clear();
    mvaStruct.numTracks_.clear();
    mvaStruct.trueBJetId_.clear();
    mvaStruct.trueBJetPt_.clear();
    mvaStruct.ptRatioTrackJet_.clear();
    mvaStruct.isMuonEvent_.clear();
        
} //END OF JET CHARGE ANALYZER FUNCTION

void AnalyzerJetCharge::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    mvaChargeTestTree = store (new TTree("mvaChargeTestTree","mvaChargeTestTree"));
    
    mvaChargeTestTree->Branch("trueBJetId", &(mvaStruct.trueBJetId_));
    mvaChargeTestTree->Branch("longChargeJet", &(mvaStruct.longChargeJet_));
    mvaChargeTestTree->Branch("leadingTrackCharge", &(mvaStruct.leadingTrackCharge_));
    mvaChargeTestTree->Branch("leadingTrackPt", &(mvaStruct.leadingTrackPt_));
    mvaChargeTestTree->Branch("leadingTrackPtCharge", &(mvaStruct.leadingTrackPtCharge_));
    mvaChargeTestTree->Branch("numTracks", &(mvaStruct.numTracks_));
    mvaChargeTestTree->Branch("relChargeJet", &(mvaStruct.relChargeJet_));
    mvaChargeTestTree->Branch("BJetPt",&(mvaStruct.trueBJetPt_));
    mvaChargeTestTree->Branch("ptRatioTrackJet",&(mvaStruct.ptRatioTrackJet_));
    mvaChargeTestTree->Branch("isMuonEvent",&(mvaStruct.isMuonEvent_));
    
    mvaChargeTrainTree = store (new TTree("mvaChargeTrainTree","mvaChargeTrainTree"));
    
    mvaChargeTrainTree->Branch("trueBJetId", &(mvaStruct.trueBJetId_));
    mvaChargeTrainTree->Branch("longChargeJet", &(mvaStruct.longChargeJet_));
    mvaChargeTrainTree->Branch("leadingTrackCharge", &(mvaStruct.leadingTrackCharge_));
    mvaChargeTrainTree->Branch("leadingTrackPt", &(mvaStruct.leadingTrackPt_));
    mvaChargeTrainTree->Branch("leadingTrackPtCharge", &(mvaStruct.leadingTrackPtCharge_));
    mvaChargeTrainTree->Branch("numTracks", &(mvaStruct.numTracks_));
    mvaChargeTrainTree->Branch("relChargeJet", &(mvaStruct.relChargeJet_));
    mvaChargeTrainTree->Branch("BJetPt",&(mvaStruct.trueBJetPt_));
    mvaChargeTrainTree->Branch("ptRatioTrackJet",&(mvaStruct.ptRatioTrackJet_));
    mvaChargeTrainTree->Branch("isMuonEvent",&(mvaStruct.isMuonEvent_));
    
    //==================TRUE LEVEL INFORMATION===========================================================================================
    
    //***ALL THE SELECTED JETS*******************************************************************************************************
    //b-quark plots
    
    name = "h_trueBJetBQuarkCharge";
    m_histogram[name] = store(new TH2D(prefix_+name+step, "b-quark flavour to reco jet charge correlation;jet charge;b-quark flavour",10,-1.5,1.5,40,-30.,30.));
    
    name = "h_trueBJetTopCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Top b quark to reco jet charge correlation;jet charge;Number of jets",10,-1.5,1.5));
    
    name = "h_trueBJetAntiTopCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Top anti-b quark to reco jet charge correlation;jet charge;Number of jets",10,-1.5,1.5));
    
    name = "h_trueBJetHiggsCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Higgs b quark to reco jet charge correlation;jet charge;Number of jets",10,-1.5,1.5));
    
    name = "h_trueBJetAntiHiggsCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Higgs anti-b quark to reco jet charge correlation;jet charge;Number of jets",10,-1.5,1.5));
    
    //b-hadron plots
    
    name = "h_trueBJetBPlusCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"B+ hadron to reco jet charge correlation;jet charge;Number of jets",10,-1.5,1.5));
    
    name = "h_trueBJetBMinusCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"B- hadron to reco jet charge correlation;jet charge;Number of jets",10,-1.5,1.5));
    
    name = "h_trueBJetBZeroCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"B0 hadron to reco jet charge correlation;jet charge;Number of jets",10,-1.5,1.5));
    
    name = "h_trueBJetAntiBZeroCharge";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Anti-B0 hadron to reco jet charge correlation;jet charge;Number of jets",10,-1.5,1.5));
    
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
    
    //===============RECO LEVEL INFORMATION================================================================================================
    
    //***ONLY TRUE-SELECTED B-JETS*******************************************************************************************************
    
    name = "h_trueBJetTrackMultiplicity";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Multiplicity of the tracks in the true b jets;multiplicity;Events",30,0,30));
    
    name = "h_trueBJetEtaVsTrackMultiplicity";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"True b jet eta to track multiplicity correlation;track multiplicity;jet eta",30,0,30,40,-3.0,3.0));
    
    name = "h_trueBJetScalarChargeVsMultip";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"True b jet charge value in function of the jet track multiplicity;track multiplicity;jet charge",20,0,40,30, -1.5, 1.5));
    
    //some control plots
    name = "h_trueBJetPt";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"True b jet p_{T};jet p_{T};Number of jets",40,0,200));
    
    name = "h_trueBJetPfTrackPt";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"True b jet track p_{T};p_{T} track;Number of tracks",40,0,40));
    
    name = "h_trueBJetHighestPtTrack";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T} value of the highest p_{T} track of each true b jet;p_{T} track;Number of jets",40,0,40));
    
    //charge validation plots
    
    name = "h_trueBJetScalarChargeValidation";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"True b jet scalar p_{T}-weighted charge;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge10";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"True b-Jet longitudinal charge calculation;jet charge;Number of jets",40,-1.2,1.2));
    
    //NEW CHARGE DEFINITIONS===============================================================================================
    
    name = "h_trueBJetRelCharge10";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet relative p_{T}-weighted charge;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetScalarChargeXWeighted";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal pT weighted charge weighted by 0.7 (as in TOP-11-031);Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge2";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal pT weighted charge with x = 0.2;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge4";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal pT weighted charge with x = 04;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge6";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal pT weighted charge with x = 0.6;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge8";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal pT weighted charge with x = 0.8;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge12";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal pT weighted charge with x = 1.2;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge14";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal pT weighted charge with x = 1.4;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge16";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal pT weighted charge with x = 1.6;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge18";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal pT weighted charge with x = 1.8;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetScalarCharge20";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet longitudinal pT weighted charge with x = 2.0;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetRelPtTrack";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Relative p_{T} of all tracks in b (and anti-b) jets;Relative p_{T};Number of jets",40,0.,200.));
    
    name = "h_trueBJetRelCharge2";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse pT weighted charge with x = 0.2;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetRelCharge4";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse pT weighted charge with x = 04;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetRelCharge6";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse pT weighted charge with x = 0.6;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetRelCharge8";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse pT weighted charge with x = 0.8;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetRelCharge12";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse pT weighted charge with x = 1.2;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetRelCharge14";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse pT weighted charge with x = 1.4;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetRelCharge16";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse pT weighted charge with x = 1.6;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetRelCharge18";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse pT weighted charge with x = 1.8;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetRelCharge20";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet transverse pT weighted charge with x = 2.0;Jet charge;Number of jets",40,-1.2,1.2));
    
    name = "h_trueBJetMaxRelPtTrack";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Relative p_{T} of highest p_{T} track in b (and anti-b) jets;Relative p_{T};Number of jets",40,0.,200.));
    
    
    //EXTRA PLOTS FOR MVA USE================================================================================================================
    
    name = "h_trueBJetBQuarkScalarCharge2";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal pT weighted charge with x = 0.2;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkScalarCharge4";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal pT weighted charge with x = 04;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkScalarCharge6";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal pT weighted charge with x = 0.6;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkScalarCharge8";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal pT weighted charge with x = 0.8;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkScalarCharge10";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"b quark to longitudinal b jet charge correlation ;c_{rel};Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkScalarCharge12";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal pT weighted charge with x = 1.2;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkScalarCharge14";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal pT weighted charge with x = 1.4;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkScalarCharge16";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal pT weighted charge with x = 1.6;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkScalarCharge18";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal pT weighted charge with x = 1.8;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkScalarCharge20";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) longitudinal pT weighted charge with x = 2.0;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkScalarCharge2";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal pT weighted charge with x = 0.2;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkScalarCharge4";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal pT weighted charge with x = 04;;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkScalarCharge6";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal pT weighted charge with x = 0.6;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkScalarCharge8";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal pT weighted charge with x = 0.8;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetAntiBQuarkScalarCharge10";
    m_histogram[name] = store(new TH1D(prefix_+name+step," Anti-b quark to longitudinal b jet charge correlation;c_{rel};Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkScalarCharge12";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal pT weighted charge with x = 1.2;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkScalarCharge14";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal pT weighted charge with x = 1.4;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkScalarCharge16";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal pT weighted charge with x = 1.6;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkScalarCharge18";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal pT weighted charge with x = 1.8;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkScalarCharge20";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) longitudinal pT weighted charge with x = 2.0;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkRelCharge2";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse pT weighted charge with x = 0.2;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkRelCharge4";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse pT weighted charge with x = 04;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkRelCharge6";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse pT weighted charge with x = 0.6;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkRelCharge8";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse pT weighted charge with x = 0.8;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkRelCharge10";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"b quark to true b jet transverse charge correlation;jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkRelCharge12";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse pT weighted charge with x = 1.2;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkRelCharge14";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse pT weighted charge with x = 1.4;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkRelCharge16";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse pT weighted charge with x = 1.6;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkRelCharge18";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse pT weighted charge with x = 1.8;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetBQuarkRelCharge20";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from b-quark) transverse pT weighted charge with x = 2.0;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkRelCharge2";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse pT weighted charge with x = 0.2;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkRelCharge4";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse pT weighted charge with x = 04;;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkRelCharge6";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse pT weighted charge with x = 0.6;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkRelCharge8";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse pT weighted charge with x = 0.8;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueBJetAntiBQuarkRelCharge10";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Anti-b quark to true b jet transverse charge correlation;jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkRelCharge12";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse pT weighted charge with x = 1.2;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkRelCharge14";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse pT weighted charge with x = 1.4;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkRelCharge16";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse pT weighted charge with x = 1.6;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkRelCharge18";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse pT weighted charge with x = 1.8;Jet charge;Number of jets",40,-1.5,1.5));
    
    name = "h_trueAntiBJetBQuarkRelCharge20";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"True b jet (from anti-b quark) transverse pT weighted charge with x = 2.0;Jet charge;Number of jets",40,-1.5,1.5));
    
    // JET-TRACK RELATION PLOTS===============================================================================================================================
    
    name = "h_trueBJetPtVsMaxPtTrack";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"true b jet p_{T} to its maximum pT track relation; max track p_{T};jet p_{T} ",40,0.,40., 40, 0.,300));
    
    name = "h_trueBJetPtVsMaxRelPtTrack";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"true b jet p_{T} to its maximum pT-rel track relation; max rel-p_{T} track; jet p_{T}",40,0.,300.,40, 0.,300.));
    
    name = "h_trueBJetPtVsMaxScalarPtTrack";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"true b jet p_{T} to its maximum pT-scalar track scalaration; max scalar-p_{T} track; jet p_{T}",40, 0.,10000.,40,0.,300.));
    
    name = "h_trueBJetPtVsNumTracks";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"true b jet p_{T} to track multiplicity correlation;track multiplicity;jet p_{T}",30,0,30,40,0,300));
    
    name = "h_trueBJetMaxPtTrackVsNumTracks";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"true b jet highest p_{T} track  to track multiplicity correlation;track multiplicity;max track p_{T}",30,0.,30.,40,0.,40.));
    
    //  LEPTON TRACK IDENTIFICATION
    
    name = "h_trueBJetLeptonTracks";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"True b jets containing a lepton among their tracks;PdgId of the track;Number of tracks",7,0,7));
    
    name = "h_trueBJetLeptonTrackMultiplicity";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Multiplicity of leptons in jets with at least one lepton track;Number of leptons; Number of tracks",10,0,10));
    
    name = "h_trueBJetLeptonTrackPt";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Pt of the leptons in jets with at least one lepton track; Pt of the lepton track; Number of tracks",40,0,40));
    
    name = "h_trueBJetLeptonTrackEta";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Eta of the leptons in jets with at least one lepton track; Eta of the lepton track; Number of tracks",40,-3,3));
    
    name = "h_trueBJetTrackMultiplicityIfLepton";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Multiplicity of the tracks in events with at least one lepton track;Track multiplicity; Number of jets",20,0,20));
    
    name = "h_trueBJetLeptonTrackCharge";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Charge of the lepton tracks; Lepton track charge;Number of tracks",2,-2,2));
    
    name = "h_trueBJetLeadingLeptonScalarChargeMatch";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Matching between the leading lepton track and jet longitudinal charge; Match of charge (1=true, 0=false); Number of tracks",2,0,2));
    
    name = "h_trueBJetLeadingLeptonRelChargeMatch";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Matching between the leading lepton track and jet transverse charge; Match of charge (1=true, 0=false); Number of tracks",2,0,2));
    
    name = "h_trueBJetNonLeadingLeptonScalarChargeMatch";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Matching between the non-leading lepton track and jet longitudinal charge; Match of charge (1=true, 0=false); Number of tracks",2,0,2));
    
    name = "h_trueBJetNonLeadingLeptonRelChargeMatch";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Matching between the non-leading lepton track and jet transverse charge; Match of charge (1=true, 0=false); Number of tracks",2,0,2));
    
    name = "h_trueBJetSecLeadingLeptonScalarChargeMatch";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Matching between the second-leading lepton track and jet longitudinal charge; Match of charge (1=true, 0=false); Number of tracks",2,0,2));
    
    name = "h_trueBJetSecLeadingLeptonRelChargeMatch";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Matching between the second-leading lepton track and jet transverse charge; Match of charge (1=true, 0=false); Number of tracks",2,0,2));
    
    name = "h_trueBJetLeptonTrackPtMultiplicity";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"Relation between track multiplicity and jet track pt for lepton-tracked cases;Pt of the jet; Track multiplicity",40,0.,400.,20,0.,20.));
    
    name = "h_trueBJetLeptonChargePt";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"Relation between pt and charge of the lepton tracks;Pt of the lepton track;Charge of the lepton track",40,0.,400.,2,-2.,2.));
    
    name = "h_trueBJetLeptonRankInTrack";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Rank of the lepton track among the tracks; Rank of the track; Number of tracks",30,0,30));
    
    name = "h_trueBJetLeptonLeadingTrackPdgId";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"PdgId of the leading track in the jet; PdgId (2=e, 3=mu); Number of tracks",7,0,7));
    
    name = "h_trueBJetLeptonSubLeadingTrackPdgId";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"PdgId of the second-leading track in the jet; PdgId (2=e, 3=mu); Number of tracks",7,0,7));
    
    name = "h_trueBJetLeptonThirdLeadingTrackPdgId";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"PdgId of the third-leading track in the jet; PdgId (2=e, 3=mu); Number of tracks",7,0,7));
    
    name = "h_trueBJetLeptonFourthLeadingTrackPdgId";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"PdgId of the fourth-leading track in the jet; PdgId (2=e, 3=mu); Number of tracks",7,0,7));
    
    name = "h_trueBJetToTrackPtRatio";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Ratio between pt of jet and pt of track for all tracks;ptTrack/ptJet ;Number of tracks",40,0.,1.));
    
    name = "h_trueBJetToLeptonTrackPtRatio";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Ratio between pt of jet and pt of track for muon tracks;ptTrack/ptJet ;Number of tracks",40,0.,1.));
    
    name = "h_trueBJetToLeadingTrackPtRatio";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Ratio between pt of jet and pt of track for leading tracks;ptTrack/ptJet ;Number of tracks",40,0.,1.));
    
    name = "h_trueBJetToLeadingLeptonTrackPtRatio";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Ratio between pt of jet and pt of track for muon leading tracks;ptTrack/ptJet ;Number of tracks",40,0.,1.));
    
    //TRUE LEPTONS
    
    name = "h_trueBJetGenLeptonPt";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Pt of the gen-leptons associated to a b-hadron;Pt;# of leptons",40,0.,400.));
    
    name = "h_trueBJetGenMuonPt";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Pt of the gen-muons associated to a b-hadron;Pt;# of leptons",40,0.,400.));
    
    name = "h_trueBJetGenLeptonPdgId";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"PdgId of the gen-leptons associated to a b-hadron;PdgId;# of leptons",36,-18,18));
    
    name = "h_trueBJetAgreementMuons";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Agreement between reco and gen muons detection;Agreement (ng+nr=0, g+nr=1, ng+r=2, g+r=3);#jets",4,0,4));
    
    name = "h_trueBJetAgreementElectrons";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Agreement between reco and gen electron detection;Agreement (ng+nr=0, g+nr=1, ng+r=2, g+r=3);#jets",4,0,4));
    
    name = "h_trueBJetAgreementChargeMuons";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Agreement between reco and gen muons charge;Charge (both b, both anti-b, mismatch;#jets",4,0,4));
    
    name = "h_trueBJetAgreementChargeElectrons";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Agreement between reco and gen electrons charge;Charge (both b, both anti-b, mismatch;#jets",4,0,4));
    
    name = "h_trueBJetGenToRecoPtDifferenceMuons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Difference between the Pt of all the gen-muons associated to a b-hadron and all the reco muon;Pt difference;# of lepton combinations",10,0.,1.));
    
    name = "h_trueBJetGenToRecoPtDifferenceMuonsExtended";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Difference between the Pt of all the gen-muons associated to a b-hadron and all the reco muon;Pt difference;Pt difference;# of lepton combinations",30,0.,30.));
    
    name = "h_trueBJetGenToRecoPhiDifferenceMuons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Difference between the Phi of all the gen-muons associated to a b-hadron and all the reco muon;Phi difference;# of lepton combinations",40,-6.,2.));
    
    name = "h_trueBJetGenToRecoEtaDifferenceMuons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Difference between the Eta of all the gen-muons associated to a b-hadron and all the reco muon;Eta difference;# of lepton combinations",40,0.,2.));
    
    name = "h_trueBJetGenToRecoDeltaRMuons";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"DeltaR of all the gen-muons associated to a b-hadron and all the reco muon;DeltaR;# of lepton combinations",40,0.,10.));
    
    //Impact parameter====================================================
    
    name = "h_trueBJetTrackIPValue";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Impact parameter value for all b jet-tracks;Impact parameter value;#tracks",40,-0.4,0.4));
    
    name = "h_trueBJetTrackIPSignificance";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Impact parameter significance for all b jet-tracks;Impact parameter significance;#tracks",40,-60.,80.));
    
    name = "h_trueBJetTrackMaxIPValue";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Max. impact parameter value for all b jets;Impact parameter value;#Jets",40,-0.4,0.4));
    
    name = "h_trueBJetTrackMaxIPSignificance";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Max. impact parameter significance for all b jets;Impact parameter significance;#Jets",40,-60.,80.));
    
    name = "h_trueBJetTrackMaxIPValueRank";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Rank of the max. impact parameter value track ;Rank position (pT);#Jets",30,0,30));
    
    name = "h_trueBJetTrackMaxIPSignificanceRank";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Rank of the max. impact parameter significance track ;Rank position (pT);#Jets",30,0,30));
    
    name = "h_trueBJetTrackIPValueVsSignificance";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"Impact parameter value vs significance for all b jet-tracks;Impact parameter value; Impact parameter significance",40,-0.4,0.4,40,-60.,80.));
    
    name = "h_trueBJetTrackPtVsIPSignificance";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"track pt vs impact parameter significance for all b jet-tracks;track pt; Impact parameter significance",40,0.,10.,40,-60.,80.));
    
    name = "h_trueBJetTrackHighestPtVsIPSignificance";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"track with highest pt vs impact parameter significance for all b jet-tracks;track pt; Impact parameter significance",40,0.,20.,40,-60.,80.));
    
    name = "h_trueBJetTrackPtVsHighestIPSignificance";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"track pt vs highest impact parameter significance for all b jet-tracks;track pt; Impact parameter significance",40,0.,20.,40,-60.,80.));
    
    //If leptons-----------------------------------//
    
    name = "h_trueBJetTrackIPValueIfLepton";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Impact parameter value for b jet muon-tracks ;Impact parameter value;#tracks",40,-0.4,0.4));
    
    name = "h_trueBJetTrackIPSignificanceIfLepton";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Impact parameter significance for b jet muon-tracks ;Impact parameter significance;#tracks",40,-60.,80.));
    
    name = "h_trueBJetTrackIPValueVsSignificanceIfLepton";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"Impact parameter value vs significance for b jet muon-tracks;Impact parameter value; Impact parameter significance",40,-0.4,0.4,40,-60.,80.));
    
    name = "h_trueBJetTrackPtVsIPSignificanceIfLepton";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"track pt vs impact parameter significance for all b jet muon-tracks;track pt; Impact parameter significance",40,0.,10.,40,-60.,80.));
    
    
    //Test PLOTS==================================================================================================================
    
    name = "checkBHadOthers";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"charge of b jets coming from other B (no Top or Higgs);charge;#jets",40,-1.2,1.2));
    
    name = "h_GenJetPt";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"p_{T} of GenJets;p_{T};#jets genjets",40,0,200));
    
    name = "h_checkPfSelMatching";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Number of selected tracks - number of matched selected tracks;# selTracks-# matched selTracks;#jets",20,0,20));
   
    name = "h_trueBJetSelTrackMultiplicity";
    m_histogram[name] = store(new TH1I(prefix_+name+step,"Multiplicity of the selected tracks;# selTrack multiplicity;#jets",30,0,200));
    
    name = "h_pfSizeVsSelSizeNoMatching";
    m_histogram[name] = store(new TH2D(prefix_+name+step,"Size of the pfCandidates vs selectedTracks;pfCandidates;selectedTracks",150,0.,150.,150,0.,150.));
    
    name = "h_checkHowManyTracks";
    m_histogram[name] = store(new TH1D(prefix_+name+step,"Size of the pfCandidates;pfCandidates;#tracks",40,0.,100.));
    
    }

//==========================================HERE MY OTHER FUNTIONS===============================================

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


