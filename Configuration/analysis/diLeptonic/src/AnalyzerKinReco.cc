#include <map>
#include <vector>
#include <iostream>

#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <Math/VectorUtil.h>

#include "AnalyzerKinReco.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"
#include "../../common/include/KinematicReconstructionSolution.h"
#include "../../common/include/KinematicReconstruction_LSroutines.h"

#include "../../common/include/KinematicReconstruction_MeanSol.h"




constexpr double TopMASS = 172.5;





AnalyzerKinReco::AnalyzerKinReco(const std::vector<TString>& selectionStepsNoCategories):
AnalyzerBaseClass("KinReco_", selectionStepsNoCategories)
{
    std::cout<<"--- Beginning setting up basic histograms\n";
    std::cout<<"=== Finishing setting up basic histograms\n\n";
}



void AnalyzerKinReco::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
        name = "nRecoEvt_vs_JetMult";
        m_histogram[name]    = store(new TH1D ( prefix_+name+step, "N reco events vs Jet Multiplicity", 10, -0.5, 9.5 ));
        
        name = "nRecoEvt_vs_LepEta";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "lead-LeptonEta - Kin. Reco. behaviour;#eta;eff.", 10, -2.5, 2.5 ));
        name = "nRecoEvt_vs_LepEta2";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "nLead-LeptonEta - Kin. Reco. behaviour;#eta;eff.", 10, -2.5, 2.5 ));
        
        name = "nRecoEvt_vs_LeppT";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "lead-LeptonpT - Kin. Reco. behaviour;p_{T}[GeV];eff.", 20, 0, 400 ));
        name = "nRecoEvt_vs_LeppT2";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "nLead-LeptonpT - Kin. Reco. behaviour;p_{T}[GeV];eff.", 20, 0, 400 ));
        name = "nRecoEvt_vs_LeppT_10Bins";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "lead-LeptonpT - Kin. Reco. behaviour;p_{T}[GeV];eff.", 10, 0, 400 ));
        
        
        name = "nRecoEvt_vs_JetEta";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "JetEta - Kin. Reco. behaviour;#eta;eff.", 10, -2.5, 2.5 ));
        name = "nRecoEvt_vs_JetEta2";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "JetEta - Kin. Reco. behaviour;#eta;eff.", 10, -2.5, 2.5 ));
        
        name = "nRecoEvt_vs_JetpT";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "JetpT - Kin. Reco. behaviour;p_{T}[GeV];eff.", 20, 0, 400 ));
        name = "nRecoEvt_vs_JetpT2";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "JetpT - Kin. Reco. behaviour;p_{T}[GeV];eff.", 20, 0, 400 ));
        
        name = "nRecoEvt_vs_LeptonEta";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "LeptonEta - Kin. Reco. behaviour;#eta;eff.", 10, -2.5, 2.5 ));
        name = "nRecoEvt_vs_AntiLeptonEta";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "AntiLeptonEta - Kin. Reco. behaviour;#eta;eff.", 10, -2.5, 2.5 ));
        
        name = "nRecoEvt_vs_LeptonpT";
        m_histogram[name]   = store(new TH1D ( prefix_+name+step, "LeptonpT - Kin. Reco. behaviour;p_{T}[GeV];eff.", 20, 0, 400 ));
        name = "nRecoEvt_vs_AntiLeptonpT";
        m_histogram[name]   = store(new TH1D ( prefix_+name+step, "AntiLeptonpT - Kin. Reco. behaviour;p_{T}[GeV];eff.", 20, 0, 400 ));
        
        name = "nRecoEvt_vs_MET";
        m_histogram[name]   = store(new TH1D ( prefix_+name+step, "MET - Kin. Reco. behaviour;MET;eff.", 10, 0, 400 ));
        
        name = "nRecoEvt_vs_top_pt";
        m_histogram[name]   = store(new TH1D ( prefix_+name+step, "p_{T}^{t} - Kin. Reco. efficiency ;p_{T}^{t};eff.", 10, 0, 1000 ));
        
        name = "nRecoEvt_Eff";
        m_histogram[name]   = store(new TH1D ( prefix_+name+step, "", 1, 0, 2 ));
        

        name = "signalTopEvents_vs_JetMult";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "N top events vs Jet Multiplicity", 10, -0.5, 9.5 ));
        name = "MatchedJets_vs_JetMult";
        m_histogram[name]  = store(new TH1D ( prefix_+name+step, "N events with 2 jets from t#bar{t}", 10, -0.5, 9.5 ));
    
        name = "nSolTtJets_vs_JetMult";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "N sol with 2 jets from t#bar{t}", 10, -0.5, 9.5 ));
        name = "nSolCorrSignJets_vs_JetMult";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "N sol with correct sign jets", 10, -0.5, 9.5 ));
        
        
        //Histograms for kinReco input, step7 needed
        name = "d_angle_jet";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, ";#alpha, rad;Entries", 200, 0, 0.5 ));
        name = "fE_jet";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "Jet Energy Correction Factor;#frac{E^{true}_{jet}}{E^{reco}_{jet}};Entries", 100, 0, 4 ));
        name = "d_angle_lep";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, ";#alpha, rad;Entries", 200, 0, 0.2 ));
        name = "fE_lep";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "Lepton Energy Correction Factor;#frac{E^{true}_{lep}}{E^{reco}_{lep}};Entries", 200, 0.5, 2.5 ));
        
        //Histograms for kinReco input, step0 needed
        name = "mbl_true";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "mbl_true", 100, 0, 180 ));
        name = "mbl_true_wrong";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "mbl_true_wrong", 100, 0, 180 ));
        
        //KinReco input W mass
        name="W_mass";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "W mass", 250, 30, 130 ));
        
        //boosted top analysis (scatter plots)
        name = "mtt_true_vs_reco_pTt300";
        m_histogram[name] = this->store(new TH2D (prefix_+name+step,"mtt true vs reco (pTt>300); Mtt(reco) [GeV];Mtt(true), [GeV]",2000,0,2000,2000,0,2000));
        name = "pTt_true_vs_reco";
        m_histogram[name] = this->store(new TH2D (prefix_+name+step,"pTt true vs reco ; pTt(reco), [GeV];pTt(true), [GeV]",1500,0,1500,1500,0,1500));
	name = "pTt_true_vs_reco_pTt300";
        m_histogram[name] = this->store(new TH2D (prefix_+name+step,"pTt true vs reco (pTt>300); pTt(reco), [GeV];pTt(true), [GeV]",1500,0,1500,1500,0,1500));
        
        //2d analysis (scatter plots)
        name = "mtt_true_vs_reco";
        m_histogram[name] = this->store(new TH2D (prefix_+name+step,"M^{t#bar{t}}_{ True}  vs  M^{t#bar{t}}_{ Reco}; M^{t#bar{t}}_{ Reco} , GeV;M^{t#bar{t}}_{ True} , GeV",2000,0,2000,2000,0,2000));
        name = "ytt_true_vs_reco";
        m_histogram[name] = this->store(new TH2D (prefix_+name+step,"y^{tt} true vs reco; y^{tt}(reco);y^{tt}(true)",2000,-5,5,2000,-5,5));
        name = "pTtt_true_vs_reco";
        m_histogram[name] = this->store(new TH2D (prefix_+name+step,"pT^{tt} true vs reco ; pT^{tt}(reco), [GeV];pT^{tt}(true), [GeV]",1500,0,1500,1500,0,1500));
        name = "yt_true_vs_reco";
        m_histogram[name] = this->store(new TH2D (prefix_+name+step,"y^{t} true vs reco; y^{t}(reco);y^{t}(true)",2000,-5,5,2000,-5,5));
        
        
        //Number of solution per gen ttbar pair.
        name = "Nsolutions";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, ";N(solutions);Fraction of events", 17, -0.5, 16.5 ));
        name = "ptGen_dN";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "ptGen - ptSol for generator ttbar input;ptGen - ptSol;Nevents", 2000, -1000, 1000 ));
        name = "ptGen_dR";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "ptGen - ptSol for generator ttbar input;ptGen - ptSol;Nevents", 2000, -1000, 1000 ));
        name = "ptGen_dRN";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "ptGen - ptSol for generator ttbar input;ptGen- ptSol;Nevents", 2000, -1000, 1000 ));
        
        name = "mttPosition";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "mtt, correct solution position", 5, -0.5, 4.5 ));
        
        name = "mttPositionMttLow";// < 500
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "mtt, correct solution position, after gen level cut", 5, -0.5, 4.5 ));
        name = "mttPositionMttMiddle";//500 800
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "mtt, correct solution position, after gen level cut", 5, -0.5, 4.5 ));
        name = "mttPositionMttHigh";// > 800
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "mtt, correct solution position, after gen level cut", 5, -0.5, 4.5 ));
        
        name = "mtt";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "mtt", 2000, 0, 2000 ));
        
        //RMS vs topPt, TopRapidity, tt-mass
        name = "RMSvsGenToppT";
        m_histogram[name] = store(new TH2D (prefix_+name+step , "RMS vs Gen", 250, 0, 1000, 250, -500, 500 ));
        name = "RMSvsGenTopRapidity";
        m_histogram[name] = store(new TH2D ( prefix_+name+step, "RMS vs Gen", 150, -3, 3, 150, -3, 3 ));
        name = "RMSvsGenTTBarMass";
        m_histogram[name] = store(new TH2D ( prefix_+name+step, "RMS vs Gen", 2000, 0, 2000, 4000, -2000, 2000 ));
        
        name = "RMSvsGenToppT_flat";
        m_histogram[name] = store(new TH2D (prefix_+name+step , "RMS vs Gen", 250, 0, 1000, 250, -500, 500 ));
        name = "RMSvsGenTopRapidity_flat";
        m_histogram[name] = store(new TH2D ( prefix_+name+step, "RMS vs Gen", 150, -3, 3, 150, -3, 3 ));
        name = "RMSvsGenTTBarMass_flat";
        m_histogram[name] = store(new TH2D ( prefix_+name+step, "RMS vs Gen", 2000, 0, 2000, 4000, -2000, 2000 ));
        
        name = "RMSvsGenToppT_posWeight";
        m_histogram[name] = store(new TH2D (prefix_+name+step , "RMS vs Gen", 250, 0, 1000, 250, -500, 500 ));
        name = "RMSvsGenTopRapidity_posWeight";
        m_histogram[name] = store(new TH2D ( prefix_+name+step, "RMS vs Gen", 150, -3, 3, 150, -3, 3 ));
        name = "RMSvsGenTTBarMass_posWeight";
        m_histogram[name] = store(new TH2D ( prefix_+name+step, "RMS vs Gen", 2000, 0, 2000, 4000, -2000, 2000 ));
        
        name = "ABSvsGenToppT";
        m_histogram[name] = store(new TProfile (prefix_+name+step, "ABS vs Gen", 200, 0, 1000));
        name = "ABSvsGenTopRapidity";
        m_histogram[name] = store(new TProfile (prefix_+name+step, "ABS vs Gen", 200, -5, 5));
        name = "ABSvsGenTTBarMass";
        m_histogram[name] = store(new TProfile (prefix_+name+step, "ABS vs Gen", 200, 0, 2000));
        
        name = "ABSvsGenToppT_flat";
        m_histogram[name] = store(new TProfile (prefix_+name+step, "ABS vs Gen", 200, 0, 1000));
        name = "ABSvsGenTopRapidity_flat";
        m_histogram[name] = store(new TProfile (prefix_+name+step, "ABS vs Gen", 200, -5, 5));
        name = "ABSvsGenTTBarMass_flat";
        m_histogram[name] = store(new TProfile (prefix_+name+step, "ABS vs Gen", 200, 0, 2000));
        
        name = "ABSvsGenToppT_posWeight";
        m_histogram[name] = store(new TProfile (prefix_+name+step, "ABS vs Gen", 200, 0, 1000));
        name = "ABSvsGenTopRapidity_posWeight";
        m_histogram[name] = store(new TProfile (prefix_+name+step, "ABS vs Gen", 200, -5, 5));
        name = "ABSvsGenTTBarMass_posWeight";
        m_histogram[name] = store(new TProfile (prefix_+name+step, "ABS vs Gen", 200, 0, 2000));
        
        name = "recoABSvsGenToppT";
        m_histogram[name] = store(new TProfile (prefix_+name+step, "reco ABS vs Gen", 200, 0.1, 1000.1));
        name = "recoABSvsGenTopRapidity";
        m_histogram[name] = store(new TProfile (prefix_+name+step, "reco ABS vs Gen", 200, -5, 5));
        name = "recoABSvsGenTTBarMass";
        m_histogram[name] = store(new TProfile (prefix_+name+step, "reco ABS vs Gen", 200, 0, 2000));
        
}



void AnalyzerKinReco::fillHistos(const EventMetadata& eventMetadata,
                                 const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                      const TopGenObjects& topGenObjects,
                                      const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                                      const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices&,
                                      const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                                      const double& weight, const TString&,
                                      std::map< TString, TH1* >& m_histogram)
{
    
    //Basic kinReco histograms
//     if(topGenObjects.valuesSet_ ){
//         
//         LV gen_top = (*topGenObjects.GenTop_);
//         LV gen_topbar = (*topGenObjects.GenAntiTop_);
//         LV gen_ttbar = gen_top + gen_topbar;
//         LV gen_b = (*topGenObjects.GenB_);
//         LV gen_bbar = (*topGenObjects.GenAntiB_);
//         LV gen_Wp = (*topGenObjects.GenWPlus_);
//         LV gen_Wm = (*topGenObjects.GenWMinus_);
//         LV gen_al = (*topGenObjects.GenAntiLepton_);
//         LV gen_l = (*topGenObjects.GenLepton_);
//         LV gen_neutrino = (*topGenObjects.GenNeutrino_);
//         LV gen_antineutrino = (*topGenObjects.GenAntiNeutrino_);        
//         
//         LV tempTop;
//         LV tempAntiTop;
//         LV tempNeutrino;
//         LV tempAntiNeutrino;
//         
//         KinematicReconstruction_LSroutines topPair(gen_top.M(), gen_topbar.M(),gen_b.M(),gen_bbar.M(), gen_Wp.M(), gen_Wm.M(), gen_al.M(), gen_l.M());
//         topPair.setConstraints(gen_al, gen_l, gen_b, gen_bbar, (gen_neutrino+gen_antineutrino).Px(), (gen_neutrino+gen_antineutrino).Py());
//         topPair.setTrueInfo(gen_top,gen_topbar,gen_neutrino,gen_antineutrino);
//         m_histogram["Nsolutions"]->Fill(topPair.getNsol());
//         
//         std::vector<double> mttPosWeight = {0.610,0.336,0.034,0.019 };
//         
//         
//         if(topPair.getNsol()>0 && topPair.getNsol()<5){
//             
//             if(gen_b.Pt()>30 && gen_bbar.Pt()>30 && gen_al.Pt()>20 && gen_l.Pt()>20 && fabs(gen_b.Eta())<2.4&& fabs(gen_bbar.Eta())<2.4&& fabs(gen_al.Eta())<2.4&& fabs(gen_l.Eta())<2.4){
//                 KinematicReconstruction_MeanSol meanSolPosWeight(TopMASS);
//                 KinematicReconstruction_MeanSol meanSolFlat(TopMASS);
//                 
//                                 for(int i=0;i<topPair.getNsol();++i)
//                                 {
//                                     meanSolFlat.add(topPair.getTtSol()->at(i).top,topPair.getTtSol()->at(i).topbar,topPair.getTtSol()->at(i).neutrino,topPair.getTtSol()->at(i).neutrinobar,1);
//                                     meanSolPosWeight.add(topPair.getTtSol()->at(i).top,topPair.getTtSol()->at(i).topbar,topPair.getTtSol()->at(i).neutrino,topPair.getTtSol()->at(i).neutrinobar,mttPosWeight.at(i));
//                                 }
//                 
//                 m_histogram["RMSvsGenToppT"]->Fill(gen_top.Pt(),topPair.getTtSol()->at(0).top.Pt()-gen_top.Pt());
//                 m_histogram["ABSvsGenToppT"]->Fill(gen_top.Pt(),std::fabs(topPair.getTtSol()->at(0).top.Pt()-gen_top.Pt())/gen_top.Pt());
//                 
//                 m_histogram["RMSvsGenTopRapidity"]->Fill(gen_top.Rapidity(),topPair.getTtSol()->at(0).top.Rapidity()-gen_top.Rapidity());
//                 m_histogram["ABSvsGenTopRapidity"]->Fill(gen_top.Rapidity(),std::fabs(topPair.getTtSol()->at(0).top.Rapidity()-gen_top.Rapidity())/std::fabs(gen_top.Rapidity()));
//                 
//                 m_histogram["RMSvsGenTTBarMass"]->Fill(gen_ttbar.M(),topPair.getTtSol()->at(0).mtt-gen_ttbar.M());
//                 m_histogram["ABSvsGenTTBarMass"]->Fill(gen_ttbar.M(),std::fabs(topPair.getTtSol()->at(0).mtt-gen_ttbar.M())/gen_ttbar.M());
//         
//                 meanSolFlat.meanSolution(tempTop, tempAntiTop, tempNeutrino, tempAntiNeutrino);
//                 m_histogram["RMSvsGenToppT_flat"]->Fill(gen_top.Pt(),tempTop.Pt()-gen_top.Pt());
//                 m_histogram["ABSvsGenToppT_flat"]->Fill(gen_top.Pt(),std::fabs(tempTop.Pt()-gen_top.Pt())/gen_top.Pt());
//                 
//                 m_histogram["RMSvsGenTopRapidity_flat"]->Fill(gen_top.Rapidity(),tempTop.Rapidity()-gen_top.Rapidity());
//                 m_histogram["ABSvsGenTopRapidity_flat"]->Fill(gen_top.Rapidity(),std::fabs(tempTop.Rapidity()-gen_top.Rapidity())/std::fabs(gen_top.Rapidity()));
//                 
//                 m_histogram["RMSvsGenTTBarMass_flat"]->Fill(gen_ttbar.M(),(tempTop+tempAntiTop).M()-gen_ttbar.M());
//                 m_histogram["ABSvsGenTTBarMass_flat"]->Fill(gen_ttbar.M(),std::fabs((tempTop+tempAntiTop).M()-gen_ttbar.M())/gen_ttbar.M());
//                 
//                 meanSolPosWeight.meanSolution(tempTop, tempAntiTop, tempNeutrino, tempAntiNeutrino);
//                 m_histogram["RMSvsGenToppT_posWeight"]->Fill(gen_top.Pt(),tempTop.Pt()-gen_top.Pt());
//                 m_histogram["ABSvsGenToppT_posWeight"]->Fill(gen_top.Pt(),std::fabs(tempTop.Pt()-gen_top.Pt())/gen_top.Pt());
//                 
//                 m_histogram["RMSvsGenTopRapidity_posWeight"]->Fill(gen_top.Rapidity(),tempTop.Rapidity()-gen_top.Rapidity());
//                 m_histogram["ABSvsGenTopRapidity_posWeight"]->Fill(gen_top.Rapidity(),std::fabs(tempTop.Rapidity()-gen_top.Rapidity())/std::fabs(gen_top.Rapidity()));
//                 
//                 m_histogram["RMSvsGenTTBarMass_posWeight"]->Fill(gen_ttbar.M(),(tempTop+tempAntiTop).M()-gen_ttbar.M());
//                 m_histogram["ABSvsGenTTBarMass_posWeight"]->Fill(gen_ttbar.M(),std::fabs((tempTop+tempAntiTop).M()-gen_ttbar.M())/gen_ttbar.M());
//             }
//                 
//                 
//                 
//                 topPair.sortBy("dN");
//                     m_histogram["ptGen_dN"]->Fill((gen_top.Pt() - ((topPair.getTtSol()->at(0)).top.Pt()))*1000);
//                         
//                     int posNumer=1;
//                     std::vector<double> mttVal;
//                     for(auto sol : (*topPair.getTtSol()))
//                         mttVal.push_back(sol.mtt);
//                     for(int i=1;i<(int)mttVal.size() ; ++i)
//                         if(mttVal.at(0)>mttVal.at(i))posNumer++;
//                     
//                     
//                     if(gen_b.Pt()>30 && gen_bbar.Pt()>30 && gen_al.Pt()>20 && gen_l.Pt()>20 && fabs(gen_b.Eta())<2.4&& fabs(gen_bbar.Eta())<2.4&& fabs(gen_al.Eta())<2.4&& fabs(gen_l.Eta())<2.4){
//                         
//                             m_histogram["mtt"]->Fill((gen_top+gen_topbar).M());
//                             m_histogram["mttPosition"]->Fill(posNumer);
//                             if(gen_ttbar.M()<500)m_histogram["mttPositionMttLow"]->Fill(posNumer);
//                             if(gen_ttbar.M()>500&&gen_ttbar.M()<800)m_histogram["mttPositionMttMiddle"]->Fill(posNumer);
//                             if(gen_ttbar.M()>800)m_histogram["mttPositionMttHigh"]->Fill(posNumer);
//                     }
// //                         
//                 topPair.sortBy("dR");
//                     m_histogram["ptGen_dR"]->Fill((gen_top.Pt() - ((topPair.getTtSol()->at(0)).top.Pt()))*1000);
//                 topPair.sortBy("dRN");
//                     m_histogram["ptGen_dRN"]->Fill((gen_top.Pt() - ((topPair.getTtSol()->at(0)).top.Pt()))*1000);
//             }
//             
//             
//         
//     }
    
    //boosted top analysis
    if(topGenObjects.valuesSet_ && kinematicReconstructionSolutions.numberOfSolutions()){
        TLorentzVector gentop = common::LVtoTLV((*topGenObjects.GenTop_));
        TLorentzVector gentopbar = common::LVtoTLV((*topGenObjects.GenAntiTop_));
        TLorentzVector recotop = common::LVtoTLV(kinematicReconstructionSolutions.solution().top());
        TLorentzVector recotopbar = common::LVtoTLV(kinematicReconstructionSolutions.solution().antiTop());
        
        TLorentzVector genttbar = gentop + gentopbar;
        TLorentzVector recottbar = recotop + recotopbar;
        
        if(std::fabs(recotop.Pt()-gentop.Pt())/gentop.Pt()<1000000) 
            m_histogram["recoABSvsGenToppT"]->Fill(gentop.Pt(),std::fabs(recotop.Pt()-gentop.Pt())/gentop.Pt());
        if(std::fabs(recotop.Rapidity()-gentop.Rapidity())/std::fabs(gentop.Rapidity())<1000000)
            m_histogram["recoABSvsGenTopRapidity"]->Fill(gentop.Rapidity(),std::fabs(recotop.Rapidity()-gentop.Rapidity())/std::fabs(gentop.Rapidity()));
        if(std::fabs(recottbar.M()-genttbar.M())/genttbar.M()<1000000)
            m_histogram["recoABSvsGenTTBarMass"]->Fill(genttbar.M(),std::fabs(recottbar.M()-genttbar.M())/genttbar.M());

        //if((gentop.Pt()>300 && recotop.Pt()>300) || (gentopbar.Pt()>300 && recotopbar.Pt()>300)){
	if( recotop.Pt()>300 || recotopbar.Pt()>300 ){
            ((TH2D*)m_histogram["mtt_true_vs_reco_pTt300"])->Fill(recottbar.M(),genttbar.M(),weight);
	    ((TH2D*)m_histogram["pTt_true_vs_reco_pTt300"])->Fill(recotop.Pt(),gentop.Pt(),weight);
        }
        ((TH2D*)m_histogram["pTt_true_vs_reco"])->Fill(recotop.Pt(),gentop.Pt());
        ((TH2D*)m_histogram["yt_true_vs_reco"])->Fill(recotop.Rapidity(),gentop.Rapidity());
        ((TH2D*)m_histogram["pTtt_true_vs_reco"])->Fill(recottbar.Pt(),genttbar.Pt());
        ((TH2D*)m_histogram["mtt_true_vs_reco"])->Fill(recottbar.M(),genttbar.M());
        ((TH2D*)m_histogram["ytt_true_vs_reco"])->Fill(recottbar.Rapidity(),genttbar.Rapidity());
    }
    // ...
    
    if(recoObjects.valuesSet_){
        double weight7 = (recoLevelWeights.weightNoPileup_)*(genLevelWeights.weightPileup_);
        weight7 *= recoLevelWeights.weightBtagSF_;
    
        int Njets = recoObjectIndices.jetIndices_.size();
    
        m_histogram["nRecoEvt_vs_JetMult"]->Fill(Njets,weight7);
        
        m_histogram["nRecoEvt_vs_LepEta"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.leadingLeptonIndex_).Eta(),weight7);
        m_histogram["nRecoEvt_vs_LepEta2"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.nLeadingLeptonIndex_).Eta(),weight7);
        m_histogram["nRecoEvt_vs_LeppT"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.leadingLeptonIndex_).Pt(),weight7);
        m_histogram["nRecoEvt_vs_LeppT_10Bins"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.leadingLeptonIndex_).Pt(),weight7);
        m_histogram["nRecoEvt_vs_LeppT2"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.nLeadingLeptonIndex_).Pt(),weight7);
        
        m_histogram["nRecoEvt_vs_LeptonEta"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.leptonIndex_).Eta(),weight7);
        m_histogram["nRecoEvt_vs_AntiLeptonEta"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.antiLeptonIndex_).Eta(),weight7);
        m_histogram["nRecoEvt_vs_LeptonpT"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.leptonIndex_).Pt(),weight7);
        m_histogram["nRecoEvt_vs_AntiLeptonpT"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.antiLeptonIndex_).Pt(),weight7);
        
        m_histogram["nRecoEvt_vs_MET"]->Fill( (*recoObjects.met_).Pt(),weight7);
        m_histogram["nRecoEvt_vs_JetEta"]->Fill((*recoObjects.jets_).at(0).Eta(),weight7);
        m_histogram["nRecoEvt_vs_JetEta2"]->Fill((*recoObjects.jets_).at(1).Eta(),weight7);
        m_histogram["nRecoEvt_vs_JetpT"]->Fill((*recoObjects.jets_).at(0).Pt(),weight7);
        m_histogram["nRecoEvt_vs_JetpT2"]->Fill((*recoObjects.jets_).at(1).Pt(),weight7);
        
        m_histogram["nRecoEvt_Eff"]->Fill(1,weight7);
        
        if(topGenObjects.valuesSet_)m_histogram["nRecoEvt_vs_top_pt"]->Fill((*topGenObjects.GenTop_).Pt(),weight7);
        
    }
        
        
    if(topGenObjects.valuesSet_ && kinematicReconstructionSolutions.numberOfSolutions()){
        
        
        int Njets = recoObjectIndices.jetIndices_.size();
        m_histogram["signalTopEvents_vs_JetMult"]->Fill(Njets,weight);
        
     
        int bjet_index=-1;
        int bbarjet_index=-1;
        int b_matched_jetIndex=-1,bbar_matched_jetIndex=-1;
        
            for(int i=0;i<(int)((*topGenObjects.genBHadFlavour_).size());i++)
            {
                if(((*topGenObjects.genBHadFlavour_).at(i))==6)bjet_index =  (*topGenObjects.genBHadJetIndex_).at(((*topGenObjects.genBHadIndex_).at(i)));
                if(((*topGenObjects.genBHadFlavour_).at(i))==-6)bbarjet_index =  (*topGenObjects.genBHadJetIndex_).at(((*topGenObjects.genBHadIndex_).at(i)));
            }
        if(bjet_index>=0&&bbarjet_index>=0&&bjet_index!=bbarjet_index){
            
            TLorentzVector truejetB = common::LVtoTLV((*topGenObjects.allGenJets_).at(bjet_index));
            TLorentzVector truejetBbar = common::LVtoTLV((*topGenObjects.allGenJets_).at(bbarjet_index));
            TLorentzVector bjetHyp=common::LVtoTLV(kinematicReconstructionSolutions.solution().bjet());
            TLorentzVector bBarjetHyp=common::LVtoTLV(kinematicReconstructionSolutions.solution().antiBjet());
            
            double dR_b_b=truejetB.DeltaR(bjetHyp);
            double dR_bar_bar=truejetBbar.DeltaR(bBarjetHyp);
            double dR_b_bar=truejetB.DeltaR(bBarjetHyp);
            double dR_bar_b=truejetBbar.DeltaR(bjetHyp);
            
                if((dR_b_b<0.3&&dR_bar_bar<0.3)||(dR_b_bar<0.3&&dR_bar_b<0.3))
                {
                        m_histogram["nSolTtJets_vs_JetMult"]->Fill(Njets,weight);
                        if((dR_b_b<0.3&&dR_bar_bar<0.3))
                        {
                            m_histogram["nSolCorrSignJets_vs_JetMult"]->Fill(Njets,weight);
                        }
                }
            
         for(const int index : (recoObjectIndices.jetIndices_))
         { 
            double dr=100;
            TLorentzVector recojet = common::LVtoTLV((*recoObjects.jets_).at(index));
            dr=recojet.DeltaR(common::LVtoTLV((*topGenObjects.allGenJets_).at(bjet_index)));
            if(dr<0.3)b_matched_jetIndex=index;
            dr=recojet.DeltaR(common::LVtoTLV((*topGenObjects.allGenJets_).at(bbarjet_index)));
            if(dr<0.3)bbar_matched_jetIndex=index;
         }
         if(b_matched_jetIndex>=0 && bbar_matched_jetIndex>=0 && b_matched_jetIndex!=bbar_matched_jetIndex)
         {
            m_histogram["MatchedJets_vs_JetMult"]->Fill(Njets,weight);
         }
        
        }
    }
    
    //fill histograms which are using in KinematicReconstruction.cc, loadData() function
    if(topGenObjects.valuesSet_ ){
        
        TLorentzVector b_true,bbar_true,l_true,al_true;
        b_true.SetXYZM((*topGenObjects.GenB_).Px(),(*topGenObjects.GenB_).Py(),(*topGenObjects.GenB_).Pz(),(*topGenObjects.GenB_).M());
        bbar_true.SetXYZM((*topGenObjects.GenAntiB_).Px(),(*topGenObjects.GenAntiB_).Py(),(*topGenObjects.GenAntiB_).Pz(),(*topGenObjects.GenAntiB_).M());
        l_true.SetXYZM((*topGenObjects.GenLepton_).Px(),(*topGenObjects.GenLepton_).Py(),(*topGenObjects.GenLepton_).Pz(),(*topGenObjects.GenLepton_).M());
        al_true.SetXYZM((*topGenObjects.GenAntiLepton_).Px(),(*topGenObjects.GenAntiLepton_).Py(),(*topGenObjects.GenAntiLepton_).Pz(),(*topGenObjects.GenAntiLepton_).M());
        
        //mbl_true, mbl_true_wrong
        m_histogram["mbl_true"]->Fill((b_true+al_true).M(),(genLevelWeights.trueLevelWeight_));
        m_histogram["mbl_true"]->Fill((bbar_true+l_true).M(),(genLevelWeights.trueLevelWeight_));
        m_histogram["mbl_true_wrong"]->Fill((b_true+l_true).M(),(genLevelWeights.trueLevelWeight_));
        m_histogram["mbl_true_wrong"]->Fill((bbar_true+al_true).M(),(genLevelWeights.trueLevelWeight_));
        
        //w mass
        //m_histogram["W_mass"]->Fill((*topGenObjects.GenWPlus_).M(),(genLevelWeights.trueLevelWeight_));
        //m_histogram["W_mass"]->Fill((*topGenObjects.GenWMinus_).M(),(genLevelWeights.trueLevelWeight_));
        
        
    if(recoObjects.valuesSet_){
        //jet
        for(const int reco_j : (recoObjectIndices.bjetIndices_))
        { 
            TLorentzVector recojet;
            recojet.SetPxPyPzE((*recoObjects.jets_).at(reco_j).Px(),(*recoObjects.jets_).at(reco_j).Py(),(*recoObjects.jets_).at(reco_j).Pz(),(*recoObjects.jets_).at(reco_j).E());
            TLorentzVector asedjet;
            asedjet.SetPxPyPzE((*commonGenObjects.associatedGenJet_).at(reco_j).Px(),(*commonGenObjects.associatedGenJet_).at(reco_j).Py(),(*commonGenObjects.associatedGenJet_).at(reco_j).Pz(),(*commonGenObjects.associatedGenJet_).at(reco_j).E());
            if(asedjet.Px()==0&&asedjet.Py()==0&&asedjet.Pz()==0&&asedjet.E()==0)continue;
            
            for(int gen_j=0;gen_j<(int)(*topGenObjects.allGenJets_).size();gen_j++)
            {
                TLorentzVector genjet;
                genjet.SetPxPyPzE((*topGenObjects.allGenJets_).at(gen_j).Px(),(*topGenObjects.allGenJets_).at(gen_j).Py(),(*topGenObjects.allGenJets_).at(gen_j).Pz(),(*topGenObjects.allGenJets_).at(gen_j).E());
                  if(genjet.Px()==asedjet.Px()&&genjet.Py()==asedjet.Py()&&genjet.Pz()==asedjet.Pz()&&genjet.E()==asedjet.E())
                  {
                        for (int IDjb=0; IDjb<(int)(*topGenObjects.genBHadJetIndex_).size(); IDjb++)
                        {
                            if(gen_j==(*topGenObjects.genBHadJetIndex_).at(IDjb))
                            {
                                if((*topGenObjects.genBHadFlavour_).at(IDjb)==6)
                                {
                                        m_histogram["fE_jet"]->Fill(b_true.E()/recojet.E());
                                        m_histogram["d_angle_jet"]->Fill(b_true.Angle(recojet.Vect()));
                                }   
                                if((*topGenObjects.genBHadFlavour_).at(IDjb)==-6)
                                {
                                        m_histogram["fE_jet"]->Fill(bbar_true.E()/recojet.E());
                                        m_histogram["d_angle_jet"]->Fill(bbar_true.Angle(recojet.Vect()));
                                } 
                            }
                        }
                  }
            }
            
        }
        //...
        //lep
        
                const int leptonIndex = recoObjectIndices.leptonIndices_.size()>0 ? recoObjectIndices.leptonIndices_.at(0) : -1;
                const int antiLeptonIndex = recoObjectIndices.antiLeptonIndices_.size()>0 ? recoObjectIndices.antiLeptonIndices_.at(0) : -1;
                const bool hasLeptonPair = (leptonIndex!=-1 && antiLeptonIndex!=-1);
                
                // Leading lepton and antilepton
                if(hasLeptonPair){

                           TLorentzVector l_temp,al_temp;
                           l_temp.SetPxPyPzE((*recoObjects.allLeptons_).at(leptonIndex).Px(),(*recoObjects.allLeptons_).at(leptonIndex).Py(),(*recoObjects.allLeptons_).at(leptonIndex).Pz(),(*recoObjects.allLeptons_).at(leptonIndex).E());
                           al_temp.SetPxPyPzE((*recoObjects.allLeptons_).at(antiLeptonIndex).Px(),(*recoObjects.allLeptons_).at(antiLeptonIndex).Py(),(*recoObjects.allLeptons_).at(antiLeptonIndex).Pz(),(*recoObjects.allLeptons_).at(antiLeptonIndex).E());
                           if(l_true.Mag()>0&&al_true.Mag()>0){

                           double dR_l=l_true.DeltaR(l_temp);
                           double dR_al=al_true.DeltaR(al_temp);
                           
                               if(dR_l<0.3&&dR_al<0.3)
                               {
                                   m_histogram["fE_lep"]->Fill(l_true.E()/l_temp.E());
                                   m_histogram["d_angle_lep"]->Fill(l_true.Angle(l_temp.Vect()));
                                   
                                   m_histogram["fE_lep"]->Fill(al_true.E()/al_temp.E());
                                   m_histogram["d_angle_lep"]->Fill(al_true.Angle(al_temp.Vect()));
                               }
                           }
                }
        //...
        }
        
    }

    
    
    
    
    
    
}








