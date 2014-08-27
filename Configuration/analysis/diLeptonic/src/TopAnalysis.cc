#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <cmath>
#include <limits>
#include <iomanip>

#include <TString.h>
#include <TMath.h>
#include <TSystem.h>
#include <Math/VectorUtil.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TIterator.h>
#include <TObject.h>

#include "TopAnalysis.h"
#include "HistoListReader.h"
#include "../../common/include/utils.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"
#include "../../common/include/ScaleFactors.h"
#include "../../common/include/analysisObjectStructs.h"
#include "analysisStructs.h"
#include "AnalyzerBaseClass.h" //FIXME: rename to AnalyzerBase.
#include "TreeHandlerBase.h"





///top production xsec in pb
constexpr double TOPXSEC = 234.;

/// Luminosity in 1/fb
constexpr double LUMI = 19.712;

///do we want to run the sync excercise?
constexpr bool RUNSYNC = false;

/// Lepton eta selection (absolute value)
constexpr double LeptonEtaCUT = 2.4;

/// Lepton pt selection in GeV
constexpr double LeptonPtCut = 20.;

/// Cut value for jet eta
constexpr double JetEtaCUT = 2.4;

/// Cut value for jet pt
constexpr double JetPtCUT = 30.;

/// Cut value for jet pt
constexpr double JetPtCUT2 = 30.;

/// B-tag algorithm and working point
constexpr Btag::Algorithm BtagALGO = Btag::csv;
constexpr Btag::WorkingPoint BtagWP = Btag::L;

/// PF MET or MVA MET
constexpr bool MvaMET = false;

/// MET selection for same-flavour channels (ee, mumu)
constexpr double MetCUT = 40.;

/// Dxy(vertex) cut for electrons
constexpr double DVertex = 0.04;

/// rho_s (mass measurement)
constexpr double rho0 = 340.0;

/// DeltaR(l,jet) cut for tt+jets
constexpr bool isttdiffXS = true;







TopAnalysis::TopAnalysis():
kinRecoOnTheFly_(false),
doClosureTest_(false),
pdf_no_(-1),
closureFunction_(nullptr),
binnedControlPlots_(0)
{
    if(MvaMET) this->mvaMet();
}




void TopAnalysis::Begin(TTree*)
{
    // Defaults from AnalysisBase
    AnalysisBase::Begin(0);
    
    // Set b-tagging working point
    this->setBtagAlgorithmAndWorkingPoint(BtagALGO, BtagWP);
    
    // Set up selection steps of tree handlers
    for(TreeHandlerBase* treeHandler : v_treeHandler_){
        if(treeHandler) treeHandler->book();
    }
    
}



void TopAnalysis::Terminate()
{
    // Produce b-tag efficiencies if required for given correction mode
    this->produceBtagEfficiencies();
    
    // Do everything needed for TTree
    for(TreeHandlerBase* treeHandler : v_treeHandler_){
        if(treeHandler){
            // Produce and write tree
            treeHandler->writeTrees(this->outputFilename(), this->channel(), this->systematic());
            //treeHandler->writeTrees(fOutput);

            // Cleanup
            treeHandler->clear();
        }
    }
    
    
    // Calculate an overall weight due to the shape reweighting, and apply it
    const double globalNormalisationFactor = this->overallGlobalNormalisationFactor();
    TIterator* it = fOutput->MakeIterator();
    while (TObject* obj = it->Next()) {
        TH1 *hist = dynamic_cast<TH1*>(obj);
        if (hist) hist->Scale(globalNormalisationFactor); 
    }
    
    // Defaults from AnalysisBase
    AnalysisBase::Terminate();
}



void TopAnalysis::SlaveBegin(TTree*)
{
    // Defaults from AnalysisBase
    AnalysisBase::SlaveBegin(0);

    //h_jetMultiAll = store(new TH1D ( "HypjetMultiAll", "Jet Multiplicity (AllJets)", 10, -0.5, 9.5 ));
    h_jetMultiXSec = store(new TH1D ( "HypjetMultiXSec", "Jet Multiplicity (for cross-section)", 10, -0.5, 9.5 ));
    h_jetMulti = store(new TH1D ( "HypjetMulti", "Jet Multiplicity", 10, -0.5, 9.5 ));
    h_jetMulti_diLep = store(new TH1D ( "HypjetMulti_diLep", "Jet Multiplicity (after dilepton)", 10, -0.5, 9.5 ));
    h_jetMulti_noBTag = store(new TH1D ( "HypjetMulti_noBTag", "Jet Multiplicity (after dilepton before btag)", 10, -0.5, 9.5 ));
    h_jetMultiNoPU = store(new TH1D ( "HypjetMultiNoPU", "Jet Multiplicity (No Pileup or lumi weight)", 10, -0.5, 9.5 ));
//     h_jetMultiVisTop = store(new TH1D ( "HypjetMultiVisTop", "Jet Multiplicity for Visible Top (No Pileup or lumi Weight)", 10, -0.5, 9.5 ));
    h_BjetMulti = store(new TH1D ( "HypBjetMulti", "B-Jet Multiplicity", 10, -0.5, 9.5 ));
    h_BjetMulti_noBTag = store(new TH1D ( "HypBjetMulti_noBTag", "B-Jet Multiplicity before B-Tag requirement", 10, -0.5, 9.5 ));

    h_HypTTBarRapidity = store(new TH1D ( "HypTTBarRapidity", "Rapidity of TTbar System (HYP)", 100, -5, 5 ));
    h_HypTTBarpT = store(new TH1D ( "HypTTBarpT", "pT of TTbar System (HYP)", 500, 0, 500 ));
    h_HypTTBarMass = store(new TH1D ( "HypTTBarMass", "Mass of TTbar System (HYP)", 2000, 0, 2000 ));
    h_HypLLBarMass = store(new TH1D ( "HypLLBarMass", "Mass of LLbar System (HYP)", 500, 0, 1000 ));
    h_HypLLBarpT = store(new TH1D ( "HypLLBarpT", "pT of LLbar System (HYP)", 200, 0, 1000 ));

    h_VisGenTTBarMass = store(new TH1D ( "VisGenTTBarMass", "Mass of TTbar System(VisGEN)", 2000, 0, 2000 ));
    h_VisGenTTBarRapidity = store(new TH1D ( "VisGenTTBarRapidity", "Rapidity of TTbar System(VisGEN)", 100, -5, 5 ));
    h_VisGenTTBarpT = store(new TH1D ( "VisGenTTBarpT", "pT of TTbar System(VisGEN)", 1200, 0, 1200 ));
    h_VisGenTopRapidity = store(new TH1D ( "VisGenTopRapidity", "Rapidity of Top(VisGEN)", 100, -5, 5 ));
    h_VisGenAntiTopRapidity = store(new TH1D ( "VisGenAntiTopRapidity", "Rapidity of AntiTop(VisGEN)", 100, -5, 5 ));

    h_VisGenLLBarpT = store(new TH1D ( "VisGenLLBarpT", "pT of LLbar System(VisGEN)", 200, 0, 1000 ));
    h_VisGenLLBarMass = store(new TH1D ( "VisGenLLBarMass", "Mass of LLbar System(VisGEN)", 500, 0, 1000 ));

    h_RecoTTBarMass = store(new TH1D ( "RecoTTBarMass","Mass of TTbar System (HYP)",2000,0,2000 ));
    h_RecoTTBarRapidity = store(new TH1D ( "RecoTTBarRapidity","Rapidity of TTbar System (HYP)",100,-5,5 ));
    h_RecoTTBarpT = store(new TH1D ( "RecoTTBarpT","pT of TTbar System (HYP)",1200,0,1200 ));
    h_RecoToppT = store(new TH1D ( "RecoToppT","pT of Top (HYP)",1000,0,1000 ));
    h_RecoAntiToppT = store(new TH1D ( "RecoAntiToppT","pT of AntiTop (HYP)",1000,0,1000 ));
    h_RecoTopRapidity = store(new TH1D ( "RecoTopRapidity","Rapidity of Top (HYP)",100,-5,5 ));
    h_RecoAntiTopRapidity = store(new TH1D ( "RecoAntiTopRapidity","Rapidity of AntiTop (HYP)",100,-5,5 ));

    h_RecoBJetpT = store(new TH1D ( "RecoBJetpT","pT of BJet (HYP)",80,0,400 ));
    h_RecoAntiBJetpT = store(new TH1D ( "RecoAntiBJetpT","pT of AntiBJet (HYP)",80,0,400 ));
    h_RecoBJetRapidity = store(new TH1D ( "RecoBJetRapidity","Rapidity of BJet (HYP)",100,-5,5 ));
    h_RecoAntiBJetRapidity = store(new TH1D ( "RecoAntiBJetRapidity","Rapidity of AntiBJet (HYP)",100,-5,5 ));
    h_RecoBJetEta = store(new TH1D ( "RecoBJetEta","#eta of BJet (HYP)",100,-5,5 ));
    h_RecoAntiBJetEta = store(new TH1D ( "RecoAntiBJetEta","#eta of AntiBJet (HYP)",100,-5,5 ));

    h_RecoLLBarMass = store(new TH1D ( "RecoLLBarMass","Mass of LLbar System (HYP)",500,0,1000 ));
    h_RecoLLBarpT = store(new TH1D ( "RecoLLBarpT","pT of LLbar System (HYP)",200,0,1000 ));

    h_VisGenAll = store(new TH1D ( "VisGenAll", "All Visible Generated particles (IM)", 40, 0, 400 ));
    h_VisGenAll_noweight = store(new TH1D ( "VisGenAll_noweight", "All Visible Generated particles (IM)", 40, 0, 400 ));
    h_GenAll = store(new TH1D ( "GenAll", "AllGenerated particles (IM)", 40, 0, 400 ));         h_GenAll->Sumw2();
    h_GenAll_noweight = store(new TH1D ( "GenAll_noweight", "AllGenerated particles (IM)", 40, 0, 400 ));         h_GenAll_noweight->Sumw2();
    h_GenAll_RecoCuts = store(new TH1D ( "GenAll_RecoCuts", "AllGenerated particles (IM)", 40, 0, 400 ));         h_GenAll_RecoCuts->Sumw2();
    h_GenAll_RecoCuts_noweight= store(new TH1D ( "GenAll_RecoCuts_noweight", "AllGenerated particles (IM)", 40, 0, 400 ));         h_GenAll_RecoCuts_noweight->Sumw2();
    
    h_diLepMassFull = store(new TH1D ( "DIMFull", "DiLepton Mass (Full Range)", 100, 0, 300 ));
    h_diLepMassFull_fullSel = store(new TH1D ( "DIMFull_fullSel", "DiLepton Mass (Full Range)", 100, 0, 300 ));
    
    h_vertMulti = store(new TH1D ( "vertMulti", "Primary Vertex Multiplicity", 30, 0, 30 ));
    h_vertMulti_noPU = store(new TH1D ( "vertMulti_noPU", "Primary Vertex Multiplicity (no Pileup)", 30, 0, 30 ));
    h_MET = store(new TH1D ( "MET", "Missing Transverse Energy", 80, 0, 400 ));
    h_jetpT = store(new TH1D ( "jetpT", "jet pT", 80, 0, 400 ));
    h_jetHT = store(new TH1D ( "jetHT", "jet HT", 80, 0, 1000 ));

    h_MuonpT = store(new TH1D ( "MuonpT", "Muon pT", 80, 0, 400 ));
    h_MuonEta = store(new TH1D ( "MuonEta", "Muon Eta", 100, -5, 5 ));
    h_ElectronpT = store(new TH1D ( "ElectronpT", "Electron pT", 80, 0, 400 ));
    h_ElectronEta = store(new TH1D ( "ElectronEta", "Electron Eta", 100, -5, 5 ));

    h_leptonPtBeforeKinReco = store(new TH1D ( "LeptonpTbkr", "Lepton pT (before kin reco)", 80, 0, 400 ));
    h_leptonPtAfterKinReco = store(new TH1D ( "LeptonpTakr", "Lepton pT (after kin reco)", 80, 0, 400 ));
    h_leptonEtaBeforeKinReco = store(new TH1D ( "LeptonEtabkr", "Lepton #eta (before kin reco)", 80, -2.5, 2.5 ));
    h_leptonEtaAfterKinReco = store(new TH1D ( "LeptonEtaakr", "Lepton #eta (after kin reco)", 80, -2.5, 2.5 ));
    h_METBeforeKinReco = store(new TH1D ( "METbkr", "Missing Transverse Energy (before kin reco)", 80, 0, 400 ));
    h_METAfterKinReco = store(new TH1D ( "METakr", "Missing Transverse Energy (after kin reco)", 80, 0, 400 ));
    h_bjetetaBeforeKinReco = store(new TH1D ( "BjetEtabkr", "b-jet eta (before kin reco)", 80, -2.5, 2.5 ));
    h_bjetetaAfterKinReco = store(new TH1D ( "BjetEtaakr", "b-jet eta (after kin reco)", 80, -2.5, 2.5 ));
    
    
    h_LeptonpT_postMETcut = store(new TH1D ( "LeptonpT_postMETcut", "Lepton pT (post MET cut)", 80, 0, 400 ));
    h_LeptonEta_postMETcut = store(new TH1D ( "LeptonEta_postMETcut", "Lepton Eta (post MET cut)", 100, -5, 5 ));
    h_AntiLeptonpT_postMETcut = store(new TH1D ( "AntiLeptonpT_postMETcut", "AntiLepton pT (post MET cut)", 80, 0, 400 ));
    h_AntiLeptonEta_postMETcut = store(new TH1D ( "AntiLeptonEta_postMETcut", "AntiLepton Eta (post MET cut)", 100, -5, 5 ));
    h_MuonpT_postMETcut = store(new TH1D ( "MuonpT_postMETcut", "Muon pT (post MET cut)", 80, 0, 400 ));
    h_MuonEta_postMETcut = store(new TH1D ( "MuonEta_postMETcut", "Muon Eta (post MET cut)", 100, -5, 5 ));
    h_ElectronpT_postMETcut = store(new TH1D ( "ElectronpT_postMETcut", "Electron pT (post MET cut)", 80, 0, 400 ));
    h_ElectronEta_postMETcut = store(new TH1D ( "ElectronEta_postMETcut", "Electron Eta (post MET cut)", 100, -5, 5 ));
    
    h_LeptonpT = store(new TH1D ( "LeptonpT", "Lepton pT", 80, 0, 400 ));
    h_LeptonEta = store(new TH1D ( "LeptonEta", "Lepton Eta", 100, -5, 5 ));
    h_LeptonpT_diLep = store(new TH1D ( "LeptonpT_diLep", "Lepton pT (after dilepton cut)", 80, 0, 400 ));
    h_LeptonEta_diLep = store(new TH1D ( "LeptonEta_diLep", "Lepton Eta (after dilepton cut)", 100, -5, 5 ));

    h_AntiLeptonpT = store(new TH1D ( "AntiLeptonpT", "AntiLepton pT", 80, 0, 400 ));
    h_AntiLeptonEta = store(new TH1D ( "AntiLeptonEta", "AntiLepton Eta", 100, -5, 5 ));
    h_AntiLeptonpT_diLep = store(new TH1D ( "AntiLeptonpT_diLep", "Lepton pT (after dilepton cut)", 80, 0, 400 ));
    h_AntiLeptonEta_diLep = store(new TH1D ( "AntiLeptonEta_diLep", "Lepton Eta (after dilepton cut)", 100, -5, 5 ));

    h_RecoToppTTTRestFrame = store(new TH1D ( "RecoToppTTTRestFrame", "Top pT in TTBar Rest Frame",1000,0,1000 ));
    h_HypToppTTTRestFrame = store(new TH1D ( "HypToppTTTRestFrame", "Top pT in TTBar Rest Frame (HYP)",1000,0,1000 ));
    h_VisGenToppTTTRestFrame = store(new TH1D ( "VisGenToppTTTRestFrame", "Top pT in TTBar Rest Frame (VisGEN)",1000,0,1000 ));
    h_GenRecoToppTTTRestFrame = store(new TH2D ( "GenRecoToppTTTRestFrame", "Gen/Reco (Top pT in TTBar Rest Frame)",1000,0,1000,1000,0,1000));

    h_RecoAntiToppTTTRestFrame = store(new TH1D ( "RecoAntiToppTTTRestFrame", "AntiTop pT in TTBar Rest Frame",1000,0,1000 ));
    h_HypAntiToppTTTRestFrame = store(new TH1D ( "HypAntiToppTTTRestFrame", "AntiTop pT in TTBar Rest Frame (HYP)",1000,0,1000 ));
    h_VisGenAntiToppTTTRestFrame = store(new TH1D ( "VisGenAntiToppTTTRestFrame", "AntiTop pT in TTBar Rest Frame (VisGEN)",1000,0,1000 ));
    h_GenRecoAntiToppTTTRestFrame = store(new TH2D ( "GenRecoAntiToppTTTRestFrame", "Gen/Reco (AntiTop pT in TTBar Rest Frame)",1000,0,1000,1000,0,1000));

    h_HypToppT = store(new TH1D ( "HypToppT", "Top pT",1000,0,1000 ));
    h_HypTopEta = store(new TH1D ( "HypTopEta", "Top #eta", 100, -5, 5 ));
    h_HypTopMass = store(new TH1D ( "HypTopMass", "Top Mass", 80, 0, 400 ));
    h_HypTopRapidity = store(new TH1D ( "HypTopRapidity", "Top Rapidity", 100, -5, 5 ));
    
    h_HypAntiToppT = store(new TH1D ( "HypAntiToppT", "AntiTop pT",1000,0,1000 ));
    h_HypAntiTopEta = store(new TH1D ( "HypAntiTopEta", "AntiTop #eta", 100, -5, 5 ));
    h_HypAntiTopMass = store(new TH1D ( "HypAntiTopMass", "AntiTop Mass", 80, 0, 400 ));
    h_HypAntiTopRapidity = store(new TH1D ( "HypAntiTopRapidity", "Top Rapidity", 100, -5, 5 ));

    h_GenRecoLeptonpT = store(new TH2D ( "GenRecoLeptonpT", "Gen/Reco Matching", 240,0,1200, 240,0,1200 ));
    h_GenRecoAntiLeptonpT = store(new TH2D ( "GenRecoAntiLeptonpT", "Gen/Reco Matching", 240,0,1200, 240,0,1200 ));
    h_HypLeptonpT = store(new TH1D ( "HypLeptonpT", "Lepton Hypothesis pT", 240,0,1200 ));
    h_HypAntiLeptonpT = store(new TH1D ( "HypAntiLeptonpT", "AntiLepton Hypothesis pT", 240,0,1200 ));
    h_VisGenLeptonpT = store(new TH1D ( "VisGenLeptonpT", "Lepton VisGenothesis pT", 240,0,1200 ));
    h_VisGenAntiLeptonpT = store(new TH1D ( "VisGenAntiLeptonpT", "AntiLepton VisGenothesis pT", 240,0,1200 ));
    h_RecoLeptonpT = store(new TH1D ( "RecoLeptonpT","pT of Lepton (HYP)",240,0,1200 ));
    h_RecoAntiLeptonpT = store(new TH1D ( "RecoAntiLeptonpT","pT of AntiLepton (HYP)",240,0,1200 ));

    h_GenRecoLeptonEta = store(new TH2D ( "GenRecoLeptonEta", "Gen/Reco Matching", 200,-5,5, 200,-5,5 ));
    h_GenRecoAntiLeptonEta = store(new TH2D ( "GenRecoAntiLeptonEta", "Gen/Reco Matching", 200,-5,5, 200,-5,5 ));
    h_HypLeptonEta = store(new TH1D ( "HypLeptonEta", "Lepton Eta", 200,-5,5 ));
    h_HypAntiLeptonEta = store(new TH1D ( "HypAntiLeptonEta", "AntiLepton Hypothesis Eta", 200,-5,5 ));
    h_VisGenLeptonEta = store(new TH1D ( "VisGenLeptonEta", "Lepton Eta", 200,-5,5 ));
    h_VisGenAntiLeptonEta = store(new TH1D ( "VisGenAntiLeptonEta", "AntiLepton VisGenothesis Eta", 200,-5,5 ));
    h_RecoLeptonEta = store(new TH1D ( "RecoLeptonEta","Eta of Lepton (HYP)",200,-5,5 ));
    h_RecoAntiLeptonEta = store(new TH1D ( "RecoAntiLeptonEta","Eta of AntiLepton (HYP)",200,-5,5 ));

    h_HypBJetpT = store(new TH1D ( "HypBJetpT", "B Hypothesis pT", 80, 0, 400 ));
    h_HypBJetEta = store(new TH1D ( "HypBJetEta", "B Hypothesis Eta", 100, -5, 5 ));
    h_HypBJetRapidity = store(new TH1D ( "HypBJetRapidity", "B Hypothesis Eta", 100, -5, 5 ));

    h_HypAntiBJetpT = store(new TH1D ( "HypAntiBJetpT", "AntiB Hypothesis pT", 80, 0, 400 ));
    h_HypAntiBJetEta = store(new TH1D ( "HypAntiBJetEta", "AntiB Hypothesis Eta", 100, -5, 5 ));
    h_HypAntiBJetRapidity = store(new TH1D ( "HypAntiBJetRapidity", "AntiB Hypothesis Eta", 100, -5, 5 ));

    h_VisGenToppT = store(new TH1D ( "VisGenToppT", "Top pT (VisGen)", 1000,0,1000 ));
    h_VisGenTopEta = store(new TH1D ( "VisGenTopEta", "Top Eta (VisGen)", 100, -5, 5 ));

    h_VisGenAntiToppT = store(new TH1D ( "VisGenAntiToppT", "AntiTop pT (VisGen)", 1000,0,1000 ));
    h_VisGenAntiTopEta = store(new TH1D ( "VisGenAntiTopEta", "AntiTop pT (VisGen)", 100, -5, 5 ));

    h_VisGenBJetpT = store(new TH1D ( "VisGenBJetpT", "B VisGenothesis pT", 80, 0, 400 ));
    h_VisGenBJetEta = store(new TH1D ( "VisGenBJetEta", "B VisGenothesis Eta", 100, -5, 5 ));
    h_VisGenBJetRapidity = store(new TH1D ( "VisGenBJetRapidity", "B VisGenothesis Rapidity", 100, -5, 5 ));

    h_VisGenAntiBJetpT = store(new TH1D ( "VisGenAntiBJetpT", "AntiB VisGenothesis pT", 80, 0, 400 ));
    h_VisGenAntiBJetEta = store(new TH1D ( "VisGenAntiBJetEta", "AntiB VisGenothesis Eta", 100, -5, 5 ));
    h_VisGenAntiBJetRapidity = store(new TH1D ( "VisGenAntiBJetRapidity", "AntiB VisGenothesis Rapidity", 100, -5, 5 ));

    h_GenRecoBJetpT = store(new TH2D ( "GenRecoBJetpT", "Gen/Reco Matching", 80, 0, 400, 80, 0, 400 ));
    h_GenRecoBJetEta = store(new TH2D ( "GenRecoBJetEta", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 ));
    h_GenRecoBJetRapidity = store(new TH2D ( "GenRecoBJetRapidity", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 ));

    h_GenRecoAntiBJetpT = store(new TH2D ( "GenRecoAntiBJetpT", "Gen/Reco Matching", 80, 0, 400, 80, 0, 400 ));
    h_GenRecoAntiBJetEta = store(new TH2D ( "GenRecoAntiBJetEta", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 ));
    h_GenRecoAntiBJetRapidity = store(new TH2D ( "GenRecoAntiBJetRapidity", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 ));

    h_GenRecoTopRapidity = store(new TH2D ( "GenRecoTopRapidity", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 ));
    h_GenRecoAntiTopRapidity = store(new TH2D ( "GenRecoAntiTopRapidity", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 ));
    h_GenRecoToppT = store(new TH2D ( "GenRecoToppT", "Gen/Reco Matching", 1000,0,1000, 1000,0,1000 ));
    h_GenRecoAntiToppT = store(new TH2D ( "GenRecoAntiToppT", "Gen/Reco Matching", 1000,0,1000, 1000,0,1000 ));

    h_GenRecoMet = store(new TH2D("GenRecoMet", "Missing ET in the event", 500, 0, 500, 500, 0, 500));
    h_VisGenMet = store(new TH1D("VisGenMet", "MET (VisGEN)", 500, 0, 500));
    h_RecoMet = store(new TH1D("RecoMet", "Reconstructed MET", 500, 0, 500));
    h_HypMet = store(new TH1D("HypMet","MET ", 500, 0, 500));

    h_GenRecoHT = store(new TH2D("GenRecoHT", "HT in the event", 800, 0, 800, 800, 0, 800));
    h_VisGenHT = store(new TH1D("VisGenHT", "HT (VisGEN)", 800, 0, 800));
    h_RecoHT = store(new TH1D("RecoHT", "Reconstructed HT", 800, 0, 800));
    h_HypHT = store(new TH1D("HypHT", "HT", 800, 0, 800));

    h_GenRecoNeutrinopT = store(new TH2D("GenRecoNeutrinopT", "Gen/Reco nu pt", 80, 0, 400, 80, 0, 400));
    h_VisGenNeutrinopT = store(new TH1D("VisGenNeutrinopT", "Nu pT (VisGEN)", 80, 0, 400));
    h_RecoNeutrinopT = store(new TH1D("RecoNeutrinopT", "reco nu pT", 80, 0, 400));
    h_HypNeutrinopT = store(new TH1D("HypNeutrinopT", "hyp nu pT", 80, 0, 400));

    h_GenRecoAntiNeutrinopT = store(new TH2D("GenRecoAntiNeutrinopT", "Gen/Reco nubar pt", 80, 0, 400, 80, 0, 400));
    h_VisGenAntiNeutrinopT = store(new TH1D("VisGenAntiNeutrinopT", "Nubar pT (VisGEN)", 80, 0, 400));
    h_RecoAntiNeutrinopT = store(new TH1D("RecoAntiNeutrinopT", "reco nubar pT", 80, 0, 400));
    h_HypAntiNeutrinopT = store(new TH1D("HypAntiNeutrinopT", "hyp nubar pT", 80, 0, 400));

    h_GenRecoTTBarDeltaPhi = store( new TH2D("GenRecoTTBarDeltaPhi", "Gen/Reco #Delta#Phi (ttbar)", 200, 0, 3.5, 200, 0, 3.5));
    h_RecoTTBarDeltaPhi = store( new TH1D("RecoTTBarDeltaPhi", "#Delta#Phi ofTTBar (RECO)", 200, 0, 3.5));
    h_HypTTBarDeltaPhi = store( new TH1D("HypTTBarDeltaPhi", "#Delta#Phi ofTTBar (HYP)", 200, 0, 3.5));
    h_VisGenTTBarDeltaPhi = store( new TH1D("VisGenTTBarDeltaPhi", "#Delta#Phi ofTTBar (VisGen)", 200, 0, 3.5));

    h_GenRecoTTBarDeltaRapidity = store( new TH2D("GenRecoTTBarDeltaRapidity", "Gen/Reco |y^{t}| - |y^{#bar{t}}|", 500, -5, 5, 500, -5, 5));
    h_RecoTTBarDeltaRapidity = store( new TH1D("RecoTTBarDeltaRapidity", "|y^{t}| - |y^{#bar{t}}| (RECO)", 500, -5, 5));
    h_HypTTBarDeltaRapidity = store( new TH1D("HypTTBarDeltaRapidity", "|y^{t}| - |y^{#bar{t}}| (HYP)", 500, -5, 5));
    h_VisGenTTBarDeltaRapidity = store( new TH1D("VisGenTTBarDeltaRapidity", "|y^{t}| - |y^{#bar{t}}| (VisGen)", 500, -5, 5));
    
    h_GenRecoBBBarpT = store( new TH2D("GenRecoBBBarpT", "Gen/Reco p_{T} (bbbar)", 400, 0, 800, 400, 0, 800));
    h_RecoBBBarpT = store( new TH1D("RecoBBBarpT", "p_{T} of bbbar (RECO)", 400, 0, 800));
    h_HypBBBarpT = store( new TH1D("HypBBBarpT", "p_{T} of bbbarpT (HYP)", 400, 0, 800));
    h_VisGenBBBarpT = store( new TH1D("VisGenBBBarpT", "p_{T} of bbbarpT (VisGen)", 400, 0, 800));

    h_GenRecoBBBarMass = store( new TH2D("GenRecoBBBarMass", "Gen/Reco Mass (bbbar)", 400, 0, 800, 400, 0, 800));
    h_RecoBBBarMass = store( new TH1D("RecoBBBarMass", "Mass of bbbar (RECO)", 400, 0, 800));
    h_HypBBBarMass = store( new TH1D("HypBBBarMass", "Mass of bbbarMass (HYP)", 400, 0, 800));
    h_VisGenBBBarMass = store( new TH1D("VisGenBBBarMass", "Mass of bbbarMass (VisGen)", 400, 0, 800));

    h_GenRecoTTBarRapidity = store(new TH2D ( "GenRecoTTBarRapidity", "Rapidity of TTbar System (HYP)", 100, -5, 5, 100, -5, 5 ));
    h_GenRecoTTBarpT = store(new TH2D ( "GenRecoTTBarpT", "pT of TTbar System (HYP)", 500, 0, 500, 500, 0, 500 ));
    h_GenRecoTTBarMass = store(new TH2D ( "GenRecoTTBarMass", "Mass of TTbar System (HYP)", 2000, 0, 2000, 2000, 0, 2000 ));
    h_GenRecoLLBarMass = store(new TH2D ( "GenRecoLLBarMass", "Mass of LLbar System (HYP)", 500, 0, 1000, 500, 0, 1000 ));
    h_GenRecoLLBarpT = store(new TH2D ( "GenRecoLLBarpT", "pT of LLbar System (HYP)", 200, 0, 1000, 200, 0, 1000 ));

    h_NJetMatching = store(new TH1D ( "NJetMatching", "NJet Gen/Reco Matching", 5, 0, 5 ));

    h_GenRecoLLBarDPhi = store(new TH2D ( "GenRecoLLBarDPhi", "Gen/Reco Matching", 320, 0., 3.2, 320, 0., 3.2 ));
    h_GenRecoLeptonantiBjetMass = store(new TH2D ( "GenRecoLeptonBjetMass", "Gen/Reco Matching", 500, 0, 1000, 500, 0, 1000 ));
    h_GenRecoAntiLeptonBjetMass = store(new TH2D ( "GenRecoAntiLeptonBjetMass", "Gen/Reco Matching", 500, 0, 1000, 500, 0, 1000 ));
    h_GenRecoJetMult = store(new TH2D ( "GenRecoJetMult", "Gen/Reco Matching", 26, -0.5, 25.5, 26, -0.5, 25.5 ));

    h_HypLLBarDPhi = store(new TH1D ( "HypLLBarDPhi", "#Delta#phi(Lep, AntiLep) (HYP)",320, 0., 3.2 ));
    h_HypLeptonantiBjetMass = store(new TH1D ( "HypLeptonBjetMass", "Mass(Lep, AntiBJet) (HYP)", 500, 0, 1000 ));
    h_HypAntiLeptonBjetMass = store(new TH1D ( "HypAntiLeptonBjetMass", "Mass(AntiLep, BJet) (HYP)", 500, 0, 1000 ));
    h_HypJetMult = store(new TH1D ( "HypJetMult", "Jet Multiplicity (HYP)", 26, -0.5, 25.5 ));

    h_VisGenLLBarDPhi = store(new TH1D ( "VisGenLLBarDPhi", "#Delta #Phi (Lep, AntiLep) (VisGEN)", 320, 0., 3.2 ));
    h_VisGenLeptonantiBjetMass = store(new TH1D ( "VisGenLeptonBjetMass", "M(Lep, AntiBJet) (VisGEN)", 500, 0, 1000 ));
    h_VisGenAntiLeptonBjetMass = store(new TH1D ( "VisGenAntiLeptonBjetMass", "M(AntiLep, BJet) (VisGEN)", 500, 0, 1000 ));
    h_VisGenJetMult = store(new TH1D ( "VisGenJetMult", "Jet Multiplicty (VisGEN)", 26, -0.5, 25.5 ));

    h_RecoLLBarDPhi = store(new TH1D ( "RecoLLBarDPhi", "#Delta #Phi (Lep, AntiLep) (Reco)", 320, 0., 3.2 ));
    h_RecoLeptonantiBjetMass = store(new TH1D ( "RecoLeptonBjetMass", "M(Lep, AntiBJet) (Reco)", 500, 0, 1000 ));
    h_RecoAntiLeptonBjetMass = store(new TH1D ( "RecoAntiLeptonBjetMass", "M(AntiLep, BJet) (Reco)", 500, 0, 1000 ));
    h_RecoJetMult = store(new TH1D ( "RecoJetMult", "Jet Multiplicty (Reco)", 26, -0.5, 25.5 ));

    h_HypToppTLead = store(new TH1D ( "HypToppTLead","Leading pT Top pT",1000,0,1000 ));
    h_RecoToppTLead = store(new TH1D ( "RecoToppTLead","Leading pT Top pT",1000,0,1000 ));
    h_VisGenToppTLead = store(new TH1D ( "VisGenToppTLead","Leading pT Top pT",1000,0,1000 ));
    h_GenRecoToppTLead = store(new TH2D ( "GenRecoToppTLead", "Gen/Reco Matching", 1000,0,1000,1000,0,1000 ));

    h_HypToppTNLead = store(new TH1D ( "HypToppTNLead","NLeading pT Top pT",1000,0,1000 ));
    h_RecoToppTNLead = store(new TH1D ( "RecoToppTNLead","NLeading pT Top pT",1000,0,1000 ));
    h_VisGenToppTNLead = store(new TH1D ( "VisGenToppTNLead","NLeading pT Top pT",1000,0,1000 ));
    h_GenRecoToppTNLead = store(new TH2D ( "GenRecoToppTNLead", "Gen/Reco Matching", 1000,0,1000,1000,0,1000 ));

    h_HypTopRapidityLead = store(new TH1D ( "HypTopRapidityLead","Leading pT Top Rapidity",100,-5,5 ));
    h_RecoTopRapidityLead = store(new TH1D ( "RecoTopRapidityLead","Leading pT Top Rapidity",100,-5,5 ));
    h_VisGenTopRapidityLead = store(new TH1D ( "VisGenTopRapidityLead","Leading pT Top Rapidity",100,-5,5 ));
    h_GenRecoTopRapidityLead = store(new TH2D ( "GenRecoTopRapidityLead", "Gen/Reco Matching", 100,-5,5,100,-5,5 ));

    h_HypTopRapidityNLead = store(new TH1D ( "HypTopRapidityNLead","NLeading pT Top Rapidity",100,-5,5 ));
    h_RecoTopRapidityNLead = store(new TH1D ( "RecoTopRapidityNLead","NLeading pT Top Rapidity",100,-5,5 ));
    h_VisGenTopRapidityNLead = store(new TH1D ( "VisGenTopRapidityNLead","NLeading pT Top Rapidity",100,-5,5 ));
    h_GenRecoTopRapidityNLead = store(new TH2D ( "GenRecoTopRapidityNLead", "Gen/Reco Matching", 100,-5,5,100,-5,5 ));

    h_HypTopMassLead = store(new TH1D ( "HypTopMassLead","Leading pT Top Mass",80,0,400 ));
    h_RecoTopMassLead = store(new TH1D ( "RecoTopMassLead","Leading pT Top Mass",80,0,400 ));
    h_VisGenTopMassLead = store(new TH1D ( "VisGenTopMassLead","Leading pT Top Mass",80,0,400 ));
    h_GenRecoTopMassLead = store(new TH2D ( "GenRecoTopMassLead", "Gen/Reco Matching", 80,0,400,80,0,400 ));

    h_HypTopMassNLead = store(new TH1D ( "HypTopMassNLead","NLeading pT Top Mass",80,0,400 ));
    h_RecoTopMassNLead = store(new TH1D ( "RecoTopMassNLead","NLeading pT Top Mass",80,0,400 ));
    h_VisGenTopMassNLead = store(new TH1D ( "VisGenTopMassNLead","NLeading pT Top Mass",80,0,400 ));
    h_GenRecoTopMassNLead = store(new TH2D ( "GenRecoTopMassNLead", "Gen/Reco Matching", 80,0,400,80,0,400 ));


    h_HypLeptonpTLead = store(new TH1D ( "HypLeptonpTLead","Leading pT Lepton pT",400,0,400 ));
    h_RecoLeptonpTLead = store(new TH1D ( "RecoLeptonpTLead","Leading pT Lepton pT",400,0,400 ));
    h_VisGenLeptonpTLead = store(new TH1D ( "VisGenLeptonpTLead","Leading pT Lepton pT",400,0,400 ));
    h_GenRecoLeptonpTLead = store(new TH2D ( "GenRecoLeptonpTLead", "Gen/Reco Matching", 400,0,400,400,0,400 ));

    h_HypLeptonpTNLead = store(new TH1D ( "HypLeptonpTNLead","NLeading pT Lepton pT",400,0,400 ));
    h_RecoLeptonpTNLead = store(new TH1D ( "RecoLeptonpTNLead","NLeading pT Lepton pT",400,0,400 ));
    h_VisGenLeptonpTNLead = store(new TH1D ( "VisGenLeptonpTNLead","NLeading pT Lepton pT",400,0,400 ));
    h_GenRecoLeptonpTNLead = store(new TH2D ( "GenRecoLeptonpTNLead", "Gen/Reco Matching", 400,0,400,400,0,400 ));

    h_HypLeptonEtaLead = store(new TH1D ( "HypLeptonEtaLead","Leading pT Lepton Eta",100,-5,5 ));
    h_RecoLeptonEtaLead = store(new TH1D ( "RecoLeptonEtaLead","Leading pT Lepton Eta",100,-5,5 ));
    h_VisGenLeptonEtaLead = store(new TH1D ( "VisGenLeptonEtaLead","Leading pT Lepton Eta",100,-5,5 ));
    h_GenRecoLeptonEtaLead = store(new TH2D ( "GenRecoLeptonEtaLead", "Gen/Reco Matching", 100,-5,5,100,-5,5 ));

    h_HypLeptonEtaNLead = store(new TH1D ( "HypLeptonEtaNLead","NLeading pT Lepton Eta",100,-5,5 ));
    h_RecoLeptonEtaNLead = store(new TH1D ( "RecoLeptonEtaNLead","NLeading pT Lepton Eta",100,-5,5 ));
    h_VisGenLeptonEtaNLead = store(new TH1D ( "VisGenLeptonEtaNLead","NLeading pT Lepton Eta",100,-5,5 ));
    h_GenRecoLeptonEtaNLead = store(new TH2D ( "GenRecoLeptonEtaNLead", "Gen/Reco Matching", 100,-5,5,100,-5,5 ));

    h_HypBJetpTLead = store(new TH1D ( "HypBJetpTLead","Leading pT BJet pT",400,0,400 ));
    h_RecoBJetpTLead = store(new TH1D ( "RecoBJetpTLead","Leading pT BJet pT",400,0,400 ));
    h_VisGenBJetpTLead = store(new TH1D ( "VisGenBJetpTLead","Leading pT BJet pT",400,0,400 ));
    h_GenRecoBJetpTLead = store(new TH2D ( "GenRecoBJetpTLead", "Gen/Reco Matching", 400,0,400,400,0,400 ));

    h_HypBJetpTNLead = store(new TH1D ( "HypBJetpTNLead","NLeading pT BJet pT",400,0,400 ));
    h_RecoBJetpTNLead = store(new TH1D ( "RecoBJetpTNLead","NLeading pT BJet pT",400,0,400 ));
    h_VisGenBJetpTNLead = store(new TH1D ( "VisGenBJetpTNLead","NLeading pT BJet pT",400,0,400 ));
    h_GenRecoBJetpTNLead = store(new TH2D ( "GenRecoBJetpTNLead", "Gen/Reco Matching", 400,0,400,400,0,400 ));

    h_HypBJetEtaLead = store(new TH1D ( "HypBJetEtaLead","Leading pT BJet Eta",100,-5,5 ));
    h_RecoBJetEtaLead = store(new TH1D ( "RecoBJetEtaLead","Leading pT BJet Eta",100,-5,5 ));
    h_VisGenBJetEtaLead = store(new TH1D ( "VisGenBJetEtaLead","Leading pT BJet Eta",100,-5,5 ));
    h_GenRecoBJetEtaLead = store(new TH2D ( "GenRecoBJetEtaLead", "Gen/Reco Matching", 100,-5,5,100,-5,5 ));

    h_HypBJetEtaNLead = store(new TH1D ( "HypBJetEtaNLead","NLeading pT BJet Eta",100,-5,5 ));
    h_RecoBJetEtaNLead = store(new TH1D ( "RecoBJetEtaNLead","NLeading pT BJet Eta",100,-5,5 ));
    h_VisGenBJetEtaNLead = store(new TH1D ( "VisGenBJetEtaNLead","NLeading pT BJet Eta",100,-5,5 ));
    h_GenRecoBJetEtaNLead = store(new TH2D ( "GenRecoBJetEtaNLead", "Gen/Reco Matching", 100,-5,5,100,-5,5 ));
    h_RecoExtraJetpT   = store(new TH1D("RecoExtraJetpT","pT of additional Jet (HYP)",80,0,400));
    h_RecoExtraJetEta  = store(new TH1D("RecoExtraJetEta","#eta of additional Jet (HYP)",100,-5,5));
    h_RecoExtraJetpT2  = store(new TH1D("RecoExtraJetpT2","pT of additional Jet (HYP)",80,0,400));
    h_RecoExtraJetEta2 = store(new TH1D("RecoExtraJetEta2","#eta of additional Jet (HYP)",100,-5,5));
    h_RecoExtraJetpT3  = store(new TH1D("RecoExtraJetpT3","pT of additional Jet (HYP)",80,0,400));
    h_RecoExtraJetEta3 = store(new TH1D("RecoExtraJetEta3","#eta of additional Jet (HYP)",100,-5,5));
    h_RecoExtraJetpT4  = store(new TH1D("RecoExtraJetpT4","pT of additional Jet (HYP)",80,0,400));
    h_RecoExtraJetEta4 = store(new TH1D("RecoExtraJetEta4","#eta of additional Jet (HYP)",100,-5,5));

    h_HypExtraJetpT   = store(new TH1D("HypExtraJetpT","pT of additional Jet",80,0,400));
    h_HypExtraJetEta  = store(new TH1D("HypExtraJetEta","#eta of additional Jet",100,-5,5));
    h_HypExtraJetpT2  = store(new TH1D("HypExtraJetpT2","pT of additional Jet",80,0,400));
    h_HypExtraJetEta2 = store(new TH1D("HypExtraJetEta2","#eta of additional Jet",100,-5,5));
    h_HypExtraJetpT3  = store(new TH1D("HypExtraJetpT3","pT of additional Jet",80,0,400));
    h_HypExtraJetEta3 = store(new TH1D("HypExtraJetEta3","#eta of additional Jet",100,-5,5));
    h_HypExtraJetpT4  = store(new TH1D("HypExtraJetpT4","pT of additional Jet",80,0,400));
    h_HypExtraJetEta4 = store(new TH1D("HypExtraJetEta4","#eta of additional Jet",100,-5,5));

    h_VisGenExtraJetpT   = store(new TH1D("VisGenExtraJetpT","pT of gen additional Jet",80,0,400));
    h_VisGenExtraJetEta  = store(new TH1D("VisGenExtraJetEta","#eta of gen additional Jet",100,-5,5));
    h_VisGenExtraJetpT2  = store(new TH1D("VisGenExtraJetpT2","pT of gen additional Jet",80,0,400));
    h_VisGenExtraJetEta2 = store(new TH1D("VisGenExtraJetEta2","#eta of gen additional Jet",100,-5,5));
    h_VisGenExtraJetpT3  = store(new TH1D("VisGenExtraJetpT3","pT of gen additional Jet",80,0,400));
    h_VisGenExtraJetEta3 = store(new TH1D("VisGenExtraJetEta3","#eta of gen additional Jet",100,-5,5));
    h_VisGenExtraJetpT4  = store(new TH1D("VisGenExtraJetpT4","pT of gen additional Jet",80,0,400));
    h_VisGenExtraJetEta4 = store(new TH1D("VisGenExtraJetEta4","#eta of gen additional Jet",100,-5,5));

    h_GenRecoExtraJetpT   = store(new TH2D("GenRecoExtraJetpT","Gen/Reco pT of additional Jet",80,0,400,80,0,400));
    h_GenRecoExtraJetEta  = store(new TH2D("GenRecoExtraJetEta","Gen/Reco #eta of additional Jet",100,-5,5,100,-5,5));
    h_GenRecoExtraJetpT2  = store(new TH2D("GenRecoExtraJetpT2","Gen/Reco pT of additional Jet",80,0,400,80,0,400));
    h_GenRecoExtraJetEta2 = store(new TH2D("GenRecoExtraJetEta2","Gen/Reco #eta of additional Jet",100,-5,5,100,-5,5));
    h_GenRecoExtraJetpT3  = store(new TH2D("GenRecoExtraJetpT3","Gen/Reco pT of additional Jet",80,0,400,80,0,400));
    h_GenRecoExtraJetEta3 = store(new TH2D("GenRecoExtraJetEta3","Gen/Reco #eta of additional Jet",100,-5,5,100,-5,5));
    h_GenRecoExtraJetpT4  = store(new TH2D("GenRecoExtraJetpT4","Gen/Reco pT of additional Jet",80,0,400,80,0,400));
    h_GenRecoExtraJetEta4 = store(new TH2D("GenRecoExtraJetEta4","Gen/Reco #eta of additional Jet",100,-5,5,100,-5,5));

    h_GenRecoJetMultpt30 = store(new TH2D("GenRecoJetMultpt30", "Gen/Reco Matching",10,-0.5,9.5,10,-0.5,9.5));
    h_RecoJetMultpt30    = store(new TH1D("RecoJetMultpt30", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_HypJetMultpt30     = store(new TH1D("HypJetMultpt30", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_VisGenJetMultpt30  = store(new TH1D("VisGenJetMultpt30", "Jet Multiplicty (VisGEN)",10,-0.5,9.5));

    h_GenRecoJetMultpt40 = store(new TH2D("GenRecoJetMultpt40", "Gen/Reco Matching",10,-0.5,9.5,10,-0.5,9.5));
    h_RecoJetMultpt40    = store(new TH1D("RecoJetMultpt40", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_HypJetMultpt40     = store(new TH1D("HypJetMultpt40", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_VisGenJetMultpt40  = store(new TH1D("VisGenJetMultpt40", "Jet Multiplicty (VisGEN)",10,-0.5,9.5));

    h_GenRecoJetMultpt60 = store(new TH2D("GenRecoJetMultpt60", "Gen/Reco Matching",10,-0.5,9.5,10,-0.5,9.5));
    h_RecoJetMultpt60    = store(new TH1D("RecoJetMultpt60", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_HypJetMultpt60     = store(new TH1D("HypJetMultpt60", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_VisGenJetMultpt60  = store(new TH1D("VisGenJetMultpt60", "Jet Multiplicty (VisGEN)",10,-0.5,9.5));

    h_GenRecoJetMultpt100 = store(new TH2D("GenRecoJetMultpt100", "Gen/Reco Matching",10,-0.5,9.5,10,-0.5,9.5));
    h_RecoJetMultpt100    = store(new TH1D("RecoJetMultpt100", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_HypJetMultpt100     = store(new TH1D("HypJetMultpt100", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_VisGenJetMultpt100  = store(new TH1D("VisGenJetMultpt100", "Jet Multiplicty (VisGEN)",10,-0.5,9.5));

    h_GenRecoDeltaRExtraJet12 = store(new TH2D("GenRecoDeltaRExtraJet12","Gen/Reco #Delta R of additional Jets",100,0.,10,100,0,10));
    h_RecoDeltaRExtraJet12 = store(new TH1D("RecoDeltaRExtraJet12","Reco #Delta R of additional Jets",100,0,10));
    h_HypDeltaRExtraJet12 = store(new TH1D("HypDeltaRExtraJet12","Hyp #DeltaR of additional Jets",100,0,10));
    h_VisGenDeltaRExtraJet12 = store(new TH1D("VisGenDeltaRExtraJet12","#Delta R of additional Jets",100,0,10));

    h_GenRecoMassExtraJet12 = store(new TH2D("GenRecoMassExtraJet12","Gen/Reco #Delta R of additional Jets",80,0,400,80,0,400));
    h_RecoMassExtraJet12 = store(new TH1D("RecoMassExtraJet12","Reco #Delta R of additional Jets",80,0,400));
    h_HypMassExtraJet12 = store(new TH1D("HypMassExtraJet12","Hyp #DeltaR of additional Jets",80,0,400));
    h_VisGenMassExtraJet12 = store(new TH1D("VisGenMassExtraJet12","#Delta R of additional Jets",80,0,400));


    h_RecoDeltaPhiExtraJet12 = store(new TH1D("RecoDeltaPhiExtraJet12","Reco #Delta#Phi of additional Jets",100,-TMath::Pi(), TMath::Pi()));
    h_HypDeltaPhiExtraJet12 = store(new TH1D("HypDeltaPhiExtraJet12","Hyp #Delta#Phi of additional Jets",100,-TMath::Pi(), TMath::Pi()));
    h_GenRecoDeltaPhiExtraJet12 = store(new TH2D("GenRecoDeltaPhiExtraJet12","Gen/Reco #Delta #phi of additional Jets",100,-TMath::Pi(), TMath::Pi(),100,0,10));
    h_VisGenDeltaPhiExtraJet12 = store(new TH1D("VisGenDeltaPhiExtraJet12","#Delta #phi of additional Jets",100,-TMath::Pi(), TMath::Pi()));

    h_RecoPhiExtraJet12 = store(new TH1D("RecoPhiExtraJet12","Reco #sum#Phi of additional Jets",100,-TMath::Pi(), TMath::Pi()));
    h_HypPhiExtraJet12 = store(new TH1D("HypPhiExtraJet12","Hyp #sum#Phi of additional Jets",100,-TMath::Pi(), TMath::Pi()));
    h_GenRecoPhiExtraJet12 = store(new TH2D("GenRecoPhiExtraJet12","Gen/Reco #Sum #phi of additional Jets",100,-TMath::Pi(), TMath::Pi(),100,-TMath::Pi(), TMath::Pi()));
    h_VisGenPhiExtraJet12 = store(new TH1D("VisGenPhiExtraJet12","#Sum #phi of additional Jets",100, -TMath::Pi(), TMath::Pi()));

    h_VisGenTTBar1stJetMass = store(new TH1D("VisGenTTBar1stJetMass","TTBar1stJetMass (VisGEN)",100,0,1));
    h_GenRecoTTBar1stJetMass = store(new TH2D("GenRecoTTBar1stJetMass","TTBar1stJetMass Gen/Reco",100,0,1,100,0,1));
    h_RecoTTBar1stJetMass = store(new TH1D("RecoTTBar1stJetMass", "TTBar1stJetMass (Reco)",100,0,1));
    h_HypTTBar1stJetMass = store(new TH1D("HypTTBar1stJetMass","TTBar1stJetMass (HYP)",100,0,1));

    h_VisGenTTBar0Mass = store(new TH1D("VisGenTTBar0Mass","TTBar0Mass (VisGEN)",200,0,2));
    h_GenRecoTTBar0Mass = store(new TH2D("GenRecoTTBar0Mass","TTBar0Mass Gen/Reco",200,0,2,200,0,2));
    h_RecoTTBar0Mass = store(new TH1D("RecoTTBar0Mass", "TTBar0Mass (Reco)",200,0,2));
    h_HypTTBar0Mass = store(new TH1D("HypTTBar0Mass","TTBar0Mass (HYP)",200,0,2));

    h_HypTopPartonFraction = store(new TH1D("HypTopPartonFraction","Parton Momentum Fraction (HYP)",100,0,1));
    h_VisGenTopPartonFraction = store(new TH1D("VisGenTopPartonFraction","Parton Momentum Fraction (VisGEN)",100,0,1));
    h_RecoTopPartonFraction = store(new TH1D("RecoTopPartonFraction","Parton Momentum Fraction (reco)",100,0,1));
    h_GenRecoTopPartonFraction = store(new TH2D("GenRecoTopPartonFraction","Parton Momentum Fraction (Gen/Reco)",100,0,1, 100, 0,1));

    h_HypAntiTopPartonFraction = store(new TH1D("HypAntiTopPartonFraction","AntiParton Momentum Fraction (HYP)",100,0,1));
    h_VisGenAntiTopPartonFraction = store(new TH1D("VisGenAntiTopPartonFraction","AntiParton Momentum Fraction (VisGEN)",100,0,1));
    h_RecoAntiTopPartonFraction = store(new TH1D("RecoAntiTopPartonFraction","AntiParton Momentum Fraction (reco)",100,0,1));
    h_GenRecoAntiTopPartonFraction = store(new TH2D("GenRecoAntiTopPartonFraction","AntiParton Momentum Fraction (Gen/Reco)",100,0,1, 100, 0,1));

//    double xbin[20]={30.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,320.,360.,400.};
    double xbin[21]={20.,30.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,320.,360.,400.};
    int nbin = 20;
    h_RecoJetMultQ0   = store(new TH1D("RecoJetMultQ0","Gap fraction Q0",nbin,xbin));
    h_HypJetMultQ0  = store(new TH1D("HypJetMultQ0",         "Gap fraction Q0 (HYP)", nbin,xbin));
    h_VisGenJetMultQ0 = store(new TH1D("VisGenJetMultQ0",         "Gap fraction Q0 (VisGEN)", nbin,xbin));
    h_GenRecoJetMultQ0 = store(new TH2D("GenRecoJetMultQ0","Gap fraction Q0", nbin,xbin,nbin,xbin));

    h_RecoJetExtra2Q0   = store(new TH1D("RecoJet2Q0","Gap fraction 2nd jet",nbin,xbin));
    h_HypJetExtra2Q0  = store(new TH1D("HypJet2Q0",         "Gap fraction 2nd jet (HYP)", nbin,xbin));
    h_VisGenJetExtra2Q0 = store(new TH1D("VisGenJet2Q0",         "Gap fraction 2nd jet (VisGEN)", nbin,xbin));
    h_GenRecoJetExtra2Q0 = store(new TH2D("GenRecoJet2Q0","Gap fraction 2nd jet", nbin,xbin,nbin,xbin));

    h_RecoJetMultQsum   = store(new TH1D("RecoJetMultQsum","Gap fraction Qsum",nbin,xbin));
    h_HypJetMultQsum  = store(new TH1D("HypJetMultQsum",         "Gap fraction Qsum (HYP)", nbin,xbin));
    h_VisGenJetMultQsum = store(new TH1D("VisGenJetMultQsum",         "Gap fraction Qsum (VisGEN)", nbin,xbin));
    h_GenRecoJetMultQsum = store(new TH2D("GenRecoJetMultQsum","Gap fraction Qsum", nbin,xbin,nbin,xbin));

    h_RecoJetMultTotal  = store(new TH1D("RecoJetMultTotal",         "Jet Multiplicity (HYP)", nbin,xbin));
    h_HypJetMultTotal  = store(new TH1D("HypJetMultTotal",         "Jet Multiplicity (HYP)", nbin,xbin));
    h_VisGenJetMultTotal = store(new TH1D("VisGenJetMultTotal",         "Gap fraction Qsum (VisGEN)", nbin,xbin));
    h_GenRecoJetMultTotal = store(new TH2D("GenRecoJetMultTotal","Gap fraction Qsum", nbin,xbin,nbin,xbin)); 

    h_ClosureTotalWeight = store(new TH1D("ClosureTotalWeight", "Total Weights from closure test",1,0,2));
    h_PDFTotalWeight = store(new TH1D("PDFTotalWeight", "PDF Weights",1,0,2));
    
    // Book histograms for b-tagging efficiencies if required for given correction mode
    this->bookBtagEfficiencyHistos();
    
    h_PUSF = store(new TH1D("PUSF", "PU SF per event", 200, 0.5, 1.5));
    h_TrigSF = store(new TH1D("TrigSF", "Trigger SF per event", 200, 0.5, 1.5));
    h_LepSF = store(new TH1D("LepSF", "Lep. Id and Isol. SF per event", 200, 0.75, 1.25));
    h_BTagSF = store(new TH1D("BTagSF", "BTagging SF per event", 200 , 0.95, 1.15 ));
    h_BTagSF->Sumw2();
    h_KinRecoSF = store(new TH1D("KinRecoSF", "Kinematic Reco. SF per event", 200, 0.5, 1.5));
    h_EventWeight = store(new TH1D("EventWeight", "Event SF", 600, 0, 3));
    
    
    // Map for binned control plots
    binnedControlPlots_ = new std::map<std::string, std::pair<TH1*, std::vector<std::map<std::string, TH1*> > > >;
    
    CreateBinnedControlPlots(h_HypToppT, h_LeptonpT);
    CreateBinnedControlPlots(h_HypToppT, h_LeptonEta);
    CreateBinnedControlPlots(h_HypToppT, h_MET);
    CreateBinnedControlPlots(h_HypToppT, h_diLepMassFull);
    
    CreateBinnedControlPlots(h_HypTopRapidity, h_LeptonpT);
    CreateBinnedControlPlots(h_HypTopRapidity, h_LeptonEta);
    CreateBinnedControlPlots(h_HypTopRapidity, h_MET);
    CreateBinnedControlPlots(h_HypTopRapidity, h_diLepMassFull);
    
    // Ievgen
        
        h_RMSvsGenToppT = store(new TH2D ( "RMSvsGenToppT", "RMS vs Gen", 500, 0, 500, 1000, -500, 500 ));
        h_RMSvsGenTopRapidity = store(new TH2D ( "RMSvsGenTopRapidity", "RMS vs Gen", 400, -5, 5, 400, -5, 5 ));
        h_RMSvsGenTTBarMass = store(new TH2D ( "RMSvsGenTTBarMass", "RMS vs Gen", 2000, 0, 2000, 4000, -2000, 2000 ));
        
        //2d cs
        
        h_HypTTBarRapidityvsTTBarpT = store(new TH2D ("HypTTBarRapidityvsTTBarpT","TTBarRapidity vs TTBarpT;p_{T}^{t#bar{t}} [GeV];y(t#bar{t})",400,0,400,100,-2.4,2.4));
        h_VisGenTTBarRapidityvsTTBarpT = store(new TH2D ("VisGenTTBarRapidityvsTTBarpT","TTBarRapidity vs TTBarpT;p_{T}^{t#bar{t}} [GeV];y(t#bar{t})",400,0,400,100,-2.4,2.4));
    // ...
    
    // Book histograms of all analyzers
    this->bookAll();
}



void TopAnalysis::SlaveTerminate()
{
    
    this->clearAll();
    
    for (auto it = binnedControlPlots_->begin(); it != binnedControlPlots_->end(); ++it) {
        delete (*it).second.first;
    }
    delete binnedControlPlots_;
    
    // Defaults from AnalysisBase
    AnalysisBase::SlaveTerminate();
}



Bool_t TopAnalysis::Process ( Long64_t entry )
{
    // Defaults from AnalysisBase
    if(!AnalysisBase::Process(entry)) return kFALSE;
    
    
    // Use utilities without namespaces
    using ROOT::Math::VectorUtil::DeltaPhi;
    using ROOT::Math::VectorUtil::DeltaR;
    using namespace common;
    
    
    // Entry for object structs are not yet read, so reset
    this->resetObjectStructEntry();


    // Define the selection steps as strings
    std::string selectionStep("");
    
    //===CUT===
    // select events on generator level and access true level weights

    // Separate DY dilepton decays in lepton flavours
    if(this->failsDrellYanGeneratorSelection(entry)) return kTRUE;

    // Separate dileptonic ttbar decays via tau
    if(this->failsTopGeneratorSelection(entry)) return kTRUE;
    
    // Count events for closure test here, where no more taus are available
    if (doClosureTest_) {
        static int closureTestEventCounter = 0;
        if (++closureTestEventCounter > closureMaxEvents_) return kTRUE;
    }
    
    // Correct for the MadGraph branching fraction being 1/9 for dileptons (PDG average is .108)
    const double weightMadgraphCorrection = this->madgraphWDecayCorrection(entry);

    // Weight due to PDF variation systematics
    const double pdfWeight = this->weightPdf(entry, pdf_no_);
    h_PDFTotalWeight->Fill(1, pdfWeight);
    
    // Get weight due to pileup reweighting
    const double weightPU = this->weightPileup(entry);
    
    // Get weight due to generator weights
    const double weightGenerator = this->weightGenerator(entry);

    // Get weight due to top-pt reweighting
    const double weightTopPt = this->weightTopPtReweighting(entry);
    
    // Get true level weights
    const double trueLevelWeightNoPileupNoClosure = weightGenerator*weightMadgraphCorrection*pdfWeight*weightTopPt;
    const double trueLevelWeightNoPileup = doClosureTest_ ? this->calculateClosureTestWeight(entry) : trueLevelWeightNoPileupNoClosure;
    const double trueLevelWeight = trueLevelWeightNoPileup*weightPU;
    
    const ttbar::GenLevelWeights genLevelWeights(weightMadgraphCorrection, weightPU, weightGenerator,
                                                 trueLevelWeightNoPileup, trueLevelWeight);
    
    const RecoObjects recoObjectsDummy;
    const CommonGenObjects commonGenObjectsDummy;

    // Access MC general generator info
    const CommonGenObjects& commonGenObjects = this->getCommonGenObjects(entry); 
    
    // Access Top signal generator info
    const TopGenObjects& topGenObjects = this->getTopGenObjects(entry);

    // Get indices of B and anti-B hadrons steming from ttbar system
    int BHadronIndex=-1;
    int AntiBHadronIndex=-1;
    this->bHadronIndices(BHadronIndex, AntiBHadronIndex, topGenObjects);
    //std::cout<<"\nINDICES: "<<BHadronIndex<<" , "<<AntiBHadronIndex<<"\n";

    // Access ttbar dilepton generator event
    LV LeadGenTop, NLeadGenTop;
    LV LeadGenLepton, NLeadGenLepton;
    LV LeadGenBJet, NLeadGenBJet;
    double genHT = -1;
    this->generatorTopEvent(LeadGenTop, NLeadGenTop,
                            LeadGenLepton, NLeadGenLepton,
                            LeadGenBJet, NLeadGenBJet,
                            genHT,
                            BHadronIndex, AntiBHadronIndex,
                            trueLevelWeightNoPileup, trueLevelWeight, topGenObjects);

    double jetHTGen = 0.;
    int GenJets_cut = -1000, GenJets_cut40 = -1000, GenJets_cut60 = -1000, GenJets_cut100 = -1000, jetnum = -1; 
    double extragenjet[4]={0.,0.,0.,0.};
    this->generatorTTbarjetsEvent(jetHTGen,
                                  BHadronIndex, AntiBHadronIndex,
                                  trueLevelWeight,
                                  GenJets_cut, GenJets_cut40, GenJets_cut60, GenJets_cut100, jetnum,
                                  extragenjet,topGenObjects);
    
    std::vector<int> genVisJetIndices;
    this->generatorVisJets(topGenObjects,genVisJetIndices);
    
    selectionStep = "0";
    
    const KinRecoObjects kinRecoObjectsDummy;
    const ttbar::GenObjectIndices genObjectIndices(-1, -1, -1, -1, -1, -1, -1, -1,genVisJetIndices);
    const ttbar::RecoObjectIndices recoObjectIndicesDummy({0},{0},{0},0,0,0,0,0,0,{0},{0});
    ttbar::RecoLevelWeights recoLevelWeightsDummy(0,0,0,0,0);
    
    this->fillAll(selectionStep,
                  recoObjectsDummy, commonGenObjects,
                  topGenObjects,
                  kinRecoObjectsDummy,
                  genObjectIndices, recoObjectIndicesDummy,
                  genLevelWeights, recoLevelWeightsDummy,
                  1.);
    
    //===CUT===
    selectionStep = "1";
    // check if event was triggered
    if(this->failsDileptonTrigger(entry)) return kTRUE;

    
    
    // === FULL RECO OBJECT SELECTION === (can thus be used at each selection step)
    
    // Access objects info
    const RecoObjects& recoObjects = this->getRecoObjects(entry);
    
    // Get allLepton indices, apply selection cuts and order them by pt (beginning with the highest value)
    const VLV& allLeptons = *recoObjects.allLeptons_;
    const std::vector<int>& lepPdgId = *recoObjects.lepPdgId_;
    std::vector<int> allLeptonIndices = initialiseIndices(allLeptons);
    selectIndices(allLeptonIndices, allLeptons, LVeta, LeptonEtaCUT, false);
    selectIndices(allLeptonIndices, allLeptons, LVeta, -LeptonEtaCUT);
    selectIndices(allLeptonIndices, allLeptons, LVpt, LeptonPtCut);
    // for electrons apply D0 cut
    selectIndices(allLeptonIndices, lepPdgId, *recoObjects.lepDxyVertex0_, DVertex, false);
    selectIndices(allLeptonIndices, lepPdgId, *recoObjects.lepDxyVertex0_, -DVertex, true);
    orderIndices(allLeptonIndices, allLeptons, LVpt);

    // Get indices of leptons and antiLeptons separated by charge, and get the leading ones if they exist
    std::vector<int> leptonIndices = allLeptonIndices;
    std::vector<int> antiLeptonIndices = allLeptonIndices;
    selectIndices(leptonIndices, lepPdgId, 0);
    selectIndices(antiLeptonIndices, lepPdgId, 0, false);
    const int numberOfLeptons = leptonIndices.size();
    const int numberOfAntiLeptons = antiLeptonIndices.size();
    const int leptonIndex = numberOfLeptons>0 ? leptonIndices.at(0) : -1;
    const int antiLeptonIndex = numberOfAntiLeptons>0 ? antiLeptonIndices.at(0) : -1;
    
    // In case of an existing opposite-charge dilepton system, get indices for leading and next-to-leading lepton
    int leadingLeptonIndex(-1);
    int nLeadingLeptonIndex(-1);
    if(numberOfLeptons>0 && numberOfAntiLeptons>0){
        leadingLeptonIndex = leptonIndex;
        nLeadingLeptonIndex = antiLeptonIndex;
        orderIndices(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, LVpt);
    }
    const bool hasLeptonPair = this->hasLeptonPair(leadingLeptonIndex, nLeadingLeptonIndex, lepPdgId);

    // Get two indices of the two leptons in the right order for trigger scale factor, if existing
    int leptonXIndex(leadingLeptonIndex);
    int leptonYIndex(nLeadingLeptonIndex);
    if(hasLeptonPair){
        //in ee and mumu channel leptonX must be the highest pt lepton, i.e. this is already correct
        // in emu channel leptonX must be electron
        if(std::abs(lepPdgId.at(leptonXIndex)) != std::abs(lepPdgId.at(leptonYIndex))){
            orderIndices(leptonYIndex, leptonXIndex, lepPdgId, true);
        }
    }
    
    // Get dilepton system, if existing
    const LV dummyLV(0.,0.,0.,0.);
    const LV dilepton(hasLeptonPair ? allLeptons.at(leadingLeptonIndex)+allLeptons.at(nLeadingLeptonIndex) : dummyLV);
    
    // Get jet indices, apply selection cuts and order them by pt (beginning with the highest value)
    const VLV& jets = *recoObjects.jets_;
    std::vector<int> jetIndices = initialiseIndices(jets);
    selectIndices(jetIndices, jets, LVeta, JetEtaCUT, false);
    selectIndices(jetIndices, jets, LVeta, -JetEtaCUT);
    selectIndices(jetIndices, jets, LVpt, JetPtCUT);
    orderIndices(jetIndices, jets, LVpt);
    const int numberOfJets = jetIndices.size();
    const bool has2Jets = numberOfJets > 1;

    // Get b-jet indices, apply selection cuts
    // and apply b-tag efficiency MC correction using random number based tag flipping (if requested correction mode is applied)
    // and order b-jets by btag discriminator (beginning with the highest value)
    const std::vector<double>& jetBTagCSV = *recoObjects.jetBTagCSV_;
    const std::vector<int>& jetPartonFlavour = *commonGenObjects.jetPartonFlavour_;
    std::vector<int> bjetIndices = jetIndices;
    selectIndices(bjetIndices, jetBTagCSV, this->btagCutValue());
    this->retagJets(bjetIndices, jetIndices, jets, jetPartonFlavour, jetBTagCSV);
    orderIndices(bjetIndices, jetBTagCSV);
    const int numberOfBjets = bjetIndices.size();
    const bool hasBtag = numberOfBjets > 0;
    
    // Get MET, and in case of MVA MET apply recoil correction for Drell-Yan sample
    this->correctMvaMet(dilepton, numberOfJets, entry);
    const LV& met = *recoObjects.met_;
    const bool hasMetOrEmu = this->channel()==Channel::emu || met.Pt()>MetCUT;
    
    const ttbar::RecoObjectIndices recoObjectIndices(allLeptonIndices,
                                                     leptonIndices, antiLeptonIndices,
                                                     leptonIndex, antiLeptonIndex,
                                                     leadingLeptonIndex, nLeadingLeptonIndex,
                                                     leptonXIndex, leptonYIndex,
                                                     jetIndices, bjetIndices);

    // Determine all reco level weights
    const double weightLeptonSF = this->weightLeptonSF(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId);
    const double weightTriggerSF = this->weightTriggerSF(leptonXIndex, leptonYIndex, allLeptons);
    const double weightNoPileup = trueLevelWeightNoPileup*weightTriggerSF*weightLeptonSF;
    const double weightBtagSF = this->weightBtagSF(jetIndices, jets, jetPartonFlavour, jetBTagCSV);
    const double weightKinReco = this->weightKinReco();
    
    // The weight to be used for filling the histograms
    double weight = weightNoPileup*weightPU;
    
    ttbar::RecoLevelWeights recoLevelWeights(weightLeptonSF, weightTriggerSF, weightBtagSF,
                                             weightNoPileup, weight);
    
    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjects,
                  kinRecoObjectsDummy,
                  genObjectIndices, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  1.);
    
    h_PUSF->Fill(weightPU, 1);
    
    //===CUT===
    selectionStep = "2";
    // we need an OS lepton pair
    if (! hasLeptonPair) return kTRUE;
    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjects,
                  kinRecoObjectsDummy,
                  genObjectIndices, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  1.);
    
    //===CUT===
    selectionStep = "3";
    // with at least 20 GeV invariant mass
    if (dilepton.M() < 20.) return kTRUE;
    
    // weight even without PU reweighting
    h_vertMulti_noPU->Fill(recoObjects.vertMulti_, weightNoPileup);
    
    h_vertMulti->Fill(recoObjects.vertMulti_, weight);
    
    h_TrigSF->Fill(weightTriggerSF, 1.);
    h_LepSF->Fill(weightLeptonSF, 1.);
    
    h_jetMulti_diLep->Fill(numberOfJets, weight);
    h_diLepMassFull->Fill(dilepton.M(), weight);
    
    if(topGenObjects.valuesSet_ && hasLeptonPair && has2Jets)
    {// Set of histograms needed to estimate the efficiency and acceptance requested by the TopXSection conveners
        h_GenAll_RecoCuts_noweight->Fill((*topGenObjects.GenTop_).M(), trueLevelWeightNoPileup);
        h_GenAll_RecoCuts->Fill((*topGenObjects.GenTop_).M(), weight);
    }

    // Access kinematic reconstruction info
    //const KinRecoObjects& kinRecoObjects = this->getKinRecoObjects(entry);
    const KinRecoObjects& kinRecoObjects = !this->makeBtagEfficiencies() ? this->getKinRecoObjectsOnTheFly(leptonIndex, antiLeptonIndex, jetIndices, bjetIndices, allLeptons, jets, jetBTagCSV, met) : kinRecoObjectsDummy;
    bool hasSolution = kinRecoObjects.valuesSet_;
    
    
    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjects,
                  kinRecoObjects,
                  genObjectIndices, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);
    
    
    
    //****************************************
    // Handle inverted Z cut
    const bool isZregion = dilepton.M() > 76. && dilepton.M() < 106.;
    if ( isZregion ) {
        double fullWeights = weight;
        selectionStep = "4zWindow";
        this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjects,
                  kinRecoObjects,
                  genObjectIndices, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  fullWeights);
        
        if ( has2Jets ) {
            selectionStep = "5zWindow";
            this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjects,
                  kinRecoObjects,
                  genObjectIndices, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  fullWeights);
            
            if ( hasMetOrEmu ) {
                selectionStep = "6zWindow";
                this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjects,
                  kinRecoObjects,
                  genObjectIndices, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  fullWeights);
                
                if ( hasBtag ) {
                    selectionStep = "7zWindow";
                    fullWeights *= weightBtagSF;
                    this->fillAll(selectionStep,
                        recoObjects, commonGenObjects,
                        topGenObjects,
                        kinRecoObjects,
                        genObjectIndices, recoObjectIndices,
                        genLevelWeights, recoLevelWeights,
                        fullWeights);
                    
                    if ( hasSolution ) {
                        fullWeights *= weightKinReco;
                        selectionStep = "8zWindow";
                        this->fillAll(selectionStep,
                            recoObjects, commonGenObjects,
                            topGenObjects,
                            kinRecoObjects,
                            genObjectIndices, recoObjectIndices,
                            genLevelWeights, recoLevelWeights,
                            fullWeights);
                    }
                }
            }
        }
    }
    
    //=== CUT ===
    selectionStep = "4";
    //Exclude the Z window
    if (this->channel() != Channel::emu && isZregion) return kTRUE;
    
    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjects,
                  kinRecoObjects,
                  genObjectIndices, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);
    
    h_LeptonpT_diLep->Fill((*recoObjects.allLeptons_).at(leptonIndex).Pt(), weight);
    h_AntiLeptonpT_diLep->Fill((*recoObjects.allLeptons_).at(antiLeptonIndex).Pt(), weight);
    h_LeptonEta_diLep->Fill((*recoObjects.allLeptons_).at(leptonIndex).Eta(), weight);
    h_AntiLeptonEta_diLep->Fill((*recoObjects.allLeptons_).at(antiLeptonIndex).Eta(), weight);
    
    h_MET->Fill(met.Pt(), weight);
    //loop over both leptons
    for (const int index : {leadingLeptonIndex, nLeadingLeptonIndex}) {
        if ( std::abs((*recoObjects.lepPdgId_).at(index)) == 11 ) {
            h_ElectronpT->Fill((*recoObjects.allLeptons_).at(index).Pt(), weight);
            h_ElectronEta->Fill((*recoObjects.allLeptons_).at(index).Eta(), weight);
        }
        else if ( std::abs((*recoObjects.lepPdgId_).at(index)) == 13 ) {
            h_MuonpT->Fill((*recoObjects.allLeptons_).at(index).Pt(), weight);
            h_MuonEta->Fill((*recoObjects.allLeptons_).at(index).Eta(), weight);
        }
    }
    
    //=== CUT ===
    selectionStep = "5";
    //Require at least two jets > 30 GeV (check for > 30 needed because we might have 20 GeV jets in our NTuple)
    if(!has2Jets) return kTRUE;
    
    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjects,
                  kinRecoObjects,
                  genObjectIndices, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);
    
    
    //=== CUT ===
    selectionStep = "6";
    //Require MET > 30 GeV in non-emu channels
    if (!hasMetOrEmu) return kTRUE;
    
    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjects,
                  kinRecoObjects,
                  genObjectIndices, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);
 
    h_LeptonpT_postMETcut->Fill((*recoObjects.allLeptons_).at(leptonIndex).Pt(), weight);
    h_AntiLeptonpT_postMETcut->Fill((*recoObjects.allLeptons_).at(antiLeptonIndex).Pt(), weight);
    h_LeptonEta_postMETcut->Fill((*recoObjects.allLeptons_).at(leptonIndex).Eta(), weight);
    h_AntiLeptonEta_postMETcut->Fill((*recoObjects.allLeptons_).at(antiLeptonIndex).Eta(), weight);

    //loop over both leptons
    for(const int index : {leadingLeptonIndex, nLeadingLeptonIndex}){
        if ( std::abs((*recoObjects.lepPdgId_).at(index)) == 11 ) {
            h_ElectronpT_postMETcut->Fill((*recoObjects.allLeptons_).at(index).Pt(), weight);
            h_ElectronEta_postMETcut->Fill((*recoObjects.allLeptons_).at(index).Eta(), weight);
        }
        else if ( std::abs((*recoObjects.lepPdgId_).at(index)) == 13 ) {
            h_MuonpT_postMETcut->Fill((*recoObjects.allLeptons_).at(index).Pt(), weight);
            h_MuonEta_postMETcut->Fill((*recoObjects.allLeptons_).at(index).Eta(), weight);
        }
    }
    
    h_jetMulti_noBTag->Fill(numberOfJets, weight);
    h_BjetMulti_noBTag->Fill(numberOfBjets, weight);

    double jetHT = getJetHT(jetIndices, (*recoObjects.jets_));
    h_jetHT->Fill(jetHT, weight);

    for ( size_t i = 0; i < 2; ++i ) {
        const int index = jetIndices.at(i);
        h_jetpT->Fill((*recoObjects.jets_).at(index).Pt(), weight);
    }
    
    // Fill b-tagging efficiencies if required for given correction mode, and in case do not process further steps
    this->fillBtagEfficiencyHistos(jetIndices, jetBTagCSV, jets, jetPartonFlavour, weight);
    if(this->makeBtagEfficiencies()) return kTRUE;

    //=== CUT ===
    selectionStep = "7";
    //Require at least one b tagged jet
    if (!hasBtag) return kTRUE;

    weight *= weightBtagSF;
    
    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjects,
                  kinRecoObjects,
                  genObjectIndices, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);
    
    h_BTagSF->Fill(weightBtagSF);
    
    if (RUNSYNC) {
        static int fullSelectionCounter = 0;
        if (fullSelectionCounter == 0)
            std::cout << "Selected#\tRun\tEvent\tlep+\tlep-\tMll\tNJets\tjet0\tjet1\tNTags\tGenJet1\tGenJet2\tMet\tGenMet\tt/tbar_decay\n"
            << std::setprecision(2) << std::fixed;
            std::cout << "Event#" << ++fullSelectionCounter << ":\t" << recoObjects.runNumber_ << "\t" << recoObjects.eventNumber_ << "\t" << (*recoObjects.allLeptons_).at(antiLeptonIndex) << "\t" << (*recoObjects.allLeptons_).at(leptonIndex) << "\t"
            << dilepton.M() << "\t" << numberOfJets << "\t"
            << (*recoObjects.jets_).at(jetIndices.at(0)) << "\t" << (*recoObjects.jets_).at(jetIndices.at(1)) << "\t" << numberOfBjets << "\t"
            << (*commonGenObjects.associatedGenJet_).at(jetIndices.at(0)) << "\t" << (*commonGenObjects.associatedGenJet_).at(jetIndices.at(1)) << "\t"
            << met.Pt() << "\t" << (*topGenObjects.GenMet_).Pt() << "\t"
            << topDecayModeString()
            << "\n";
    }
    
    h_BjetMulti->Fill(numberOfBjets, weight);
    h_jetMulti->Fill(numberOfJets, weight);
    
    h_leptonPtBeforeKinReco->Fill((*recoObjects.allLeptons_).at(leptonIndex).Pt(), weight);
    h_leptonPtBeforeKinReco->Fill((*recoObjects.allLeptons_).at(antiLeptonIndex).Pt(), weight);
    h_leptonEtaBeforeKinReco->Fill((*recoObjects.allLeptons_).at(leptonIndex).Eta(), weight);
    h_leptonEtaBeforeKinReco->Fill((*recoObjects.allLeptons_).at(antiLeptonIndex).Eta(), weight);
    h_METBeforeKinReco->Fill(met.Pt(), weight);
    for (const int index : bjetIndices)
        h_bjetetaBeforeKinReco->Fill((*recoObjects.jets_).at(index).Eta(), weight);
    
    //...
    //=== CUT ===
    selectionStep = "8";
    //Require at least one solution for the kinematic event reconstruction
    if (!hasSolution) return kTRUE;
    
    weight *= weightKinReco;
    
    this->fillAll(selectionStep,
                  recoObjects, commonGenObjects,
                  topGenObjects,
                  kinRecoObjects,
                  genObjectIndices, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight);
    
    if(topGenObjects.valuesSet_){//Ievgen
        h_RMSvsGenToppT->Fill((*topGenObjects.GenTop_).Pt(),(*topGenObjects.GenTop_).Pt()-(*kinRecoObjects.HypTop_).at(0).Pt());
        h_RMSvsGenToppT->Fill((*topGenObjects.GenAntiTop_).Pt(),(*topGenObjects.GenAntiTop_).Pt()-(*kinRecoObjects.HypAntiTop_).at(0).Pt());
   
        h_RMSvsGenTopRapidity->Fill((*topGenObjects.GenTop_).Rapidity(),(*topGenObjects.GenTop_).Rapidity()-((*kinRecoObjects.HypTop_).at(0)).Rapidity());
        h_RMSvsGenTopRapidity->Fill((*topGenObjects.GenAntiTop_).Rapidity(),(*topGenObjects.GenAntiTop_).Rapidity()-((*kinRecoObjects.HypAntiTop_).at(0)).Rapidity());
    
        h_RMSvsGenTTBarMass->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).M(),((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).M()-((*kinRecoObjects.HypTop_).at(0)+(*kinRecoObjects.HypAntiTop_).at(0)).M());
    }
    
    
    h_leptonPtAfterKinReco->Fill((*recoObjects.allLeptons_).at(leptonIndex).Pt(), weight);
    h_leptonPtAfterKinReco->Fill((*recoObjects.allLeptons_).at(antiLeptonIndex).Pt(), weight);
    h_leptonEtaAfterKinReco->Fill((*recoObjects.allLeptons_).at(leptonIndex).Eta(), weight);
    h_leptonEtaAfterKinReco->Fill((*recoObjects.allLeptons_).at(antiLeptonIndex).Eta(), weight);
    h_METAfterKinReco->Fill(met.Pt(), weight);
    for (const int index : bjetIndices)
        h_bjetetaAfterKinReco->Fill((*recoObjects.jets_).at(index).Eta(), weight);
    
    h_KinRecoSF->Fill(weightKinReco, 1);
    h_EventWeight->Fill(weight, 1);
    
    
    h_jetMultiXSec->Fill(numberOfJets, weight);
    h_jetMultiNoPU->Fill(numberOfJets, weight / weightPU );
    h_diLepMassFull_fullSel->Fill(dilepton.M(), weight);
        
    //create helper variables
    size_t solutionIndex = 0; //always zero!
    
    // Find 1st (and 2nd) leading pT particles: Top, Lepton, BJetIndex
    LV LeadHypTop, NLeadHypTop;
    LV LeadHypLepton, NLeadHypLepton;
    LV LeadHypBJet, NLeadHypBJet;
    orderLV(LeadHypTop, NLeadHypTop, (*kinRecoObjects.HypTop_).at(solutionIndex), (*kinRecoObjects.HypAntiTop_).at(solutionIndex), LVpt);
    orderLV(LeadHypLepton, NLeadHypLepton, (*kinRecoObjects.HypLepton_).at(solutionIndex), (*kinRecoObjects.HypAntiLepton_).at(solutionIndex), LVpt);
    orderLV(LeadHypBJet, NLeadHypBJet, (*kinRecoObjects.HypBJet_).at(solutionIndex), (*kinRecoObjects.HypAntiBJet_).at(solutionIndex), LVpt);

    //create ll, bb and tt system 
    LV hypllbar((*kinRecoObjects.HypLepton_).at(solutionIndex) + (*kinRecoObjects.HypAntiLepton_).at(solutionIndex));
    LV hypbbbar((*kinRecoObjects.HypBJet_).at(solutionIndex) + (*kinRecoObjects.HypAntiBJet_).at(solutionIndex));
    LV hypttbar((*kinRecoObjects.HypTop_).at(solutionIndex)+(*kinRecoObjects.HypAntiTop_).at(solutionIndex));

    // create top/antitop quark in the ttbar rest frame
    ROOT::Math::Boost CoMBoostHypTtbar (hypttbar.BoostToCM());
    LV top = (*kinRecoObjects.HypTop_).at(solutionIndex);
    LV antitop = (*kinRecoObjects.HypAntiTop_).at(solutionIndex);
    top = CoMBoostHypTtbar(top);
    antitop = CoMBoostHypTtbar(antitop);

    //First fill the reco histograms (which have no scaling factors applied)
    const double recoWeight = trueLevelWeight;
    h_RecoTTBarDeltaRapidity->Fill(std::abs((*kinRecoObjects.HypTop_).at(solutionIndex).Rapidity()) - std::abs((*kinRecoObjects.HypAntiTop_).at(solutionIndex).Rapidity()), recoWeight);
    h_RecoTTBarDeltaPhi->Fill(std::abs(DeltaPhi((*kinRecoObjects.HypTop_).at(solutionIndex), (*kinRecoObjects.HypAntiTop_).at(solutionIndex))), recoWeight);
    h_RecoTTBarMass->Fill(hypttbar.M(), recoWeight);
    h_RecoTTBarRapidity->Fill(hypttbar.Rapidity(), recoWeight);
    h_RecoTTBarpT->Fill(hypttbar.Pt(), recoWeight);
    h_RecoToppT->Fill((*kinRecoObjects.HypTop_).at(solutionIndex).Pt(), recoWeight);
    h_RecoAntiToppT->Fill((*kinRecoObjects.HypAntiTop_).at(solutionIndex).Pt(), recoWeight);
    h_RecoTopRapidity->Fill((*kinRecoObjects.HypTop_).at(solutionIndex).Rapidity(), recoWeight);
    h_RecoAntiTopRapidity->Fill((*kinRecoObjects.HypAntiTop_).at(solutionIndex).Rapidity(), recoWeight);

    h_RecoToppTTTRestFrame->Fill(top.Pt(), recoWeight);
    h_RecoAntiToppTTTRestFrame->Fill(antitop.Pt(), recoWeight);
    
    h_RecoLLBarMass->Fill(hypllbar.M(), recoWeight);
    h_RecoLLBarpT->Fill(hypllbar.Pt(), recoWeight);
    h_RecoLeptonpT->Fill((*kinRecoObjects.HypLepton_).at(solutionIndex).Pt(), recoWeight);
    h_RecoAntiLeptonpT->Fill((*kinRecoObjects.HypAntiLepton_).at(solutionIndex).Pt(), recoWeight);
    h_RecoLeptonEta->Fill((*kinRecoObjects.HypLepton_).at(solutionIndex).Eta(), recoWeight);
    h_RecoAntiLeptonEta->Fill((*kinRecoObjects.HypAntiLepton_).at(solutionIndex).Eta(), recoWeight);
    
    h_RecoMet->Fill(met.Pt(), recoWeight);
    h_RecoHT->Fill(jetHT, recoWeight);
    
    h_RecoNeutrinopT->Fill((*kinRecoObjects.HypNeutrino_).at(solutionIndex).Pt(), recoWeight);
    h_RecoAntiNeutrinopT->Fill((*kinRecoObjects.HypAntiNeutrino_).at(solutionIndex).Pt(), recoWeight);

    h_RecoBBBarMass->Fill(hypbbbar.M(), recoWeight);
    h_RecoBBBarpT->Fill(hypbbbar.Pt(), recoWeight);

    h_RecoBJetpT->Fill((*kinRecoObjects.HypBJet_).at(solutionIndex).Pt(), recoWeight);
    h_RecoAntiBJetpT->Fill((*kinRecoObjects.HypAntiBJet_).at(solutionIndex).Pt(), recoWeight);
    h_RecoBJetRapidity->Fill((*kinRecoObjects.HypBJet_).at(solutionIndex).Rapidity(), recoWeight);
    h_RecoAntiBJetRapidity->Fill((*kinRecoObjects.HypAntiBJet_).at(solutionIndex).Rapidity(), recoWeight);
    h_RecoBJetEta->Fill((*kinRecoObjects.HypBJet_).at(solutionIndex).Eta(), recoWeight);
    h_RecoAntiBJetEta->Fill((*kinRecoObjects.HypAntiBJet_).at(solutionIndex).Eta(), recoWeight);

    h_RecoLLBarDPhi->Fill(std::fabs ( DeltaPhi ( (*kinRecoObjects.HypLepton_).at(solutionIndex), (*kinRecoObjects.HypAntiLepton_).at(solutionIndex) ) ), recoWeight);
    h_RecoLeptonantiBjetMass->Fill(( (*kinRecoObjects.HypLepton_).at(solutionIndex)+(*kinRecoObjects.HypAntiBJet_).at(solutionIndex) ).M(), recoWeight);
    h_RecoAntiLeptonBjetMass->Fill(( (*kinRecoObjects.HypAntiLepton_).at(solutionIndex)+(*kinRecoObjects.HypBJet_).at(solutionIndex) ).M(), recoWeight);

    h_RecoToppTLead->Fill(LeadHypTop.Pt(), recoWeight);
    h_RecoToppTNLead->Fill(NLeadHypTop.Pt(), recoWeight);
    h_RecoTopRapidityLead->Fill(LeadHypTop.Rapidity(), recoWeight);
    h_RecoTopRapidityNLead->Fill(NLeadHypTop.Rapidity(), recoWeight);
    h_RecoTopMassLead->Fill(LeadHypTop.M(), recoWeight);
    h_RecoTopMassNLead->Fill(NLeadHypTop.M(), recoWeight);
    
    h_RecoLeptonpTLead->Fill(LeadHypLepton.Pt(), recoWeight);
    h_RecoLeptonpTNLead->Fill(NLeadHypLepton.Pt(), recoWeight);
    h_RecoLeptonEtaLead->Fill(LeadHypLepton.Eta(), recoWeight);
    h_RecoLeptonEtaNLead->Fill(NLeadHypLepton.Eta(), recoWeight);
    
    h_RecoBJetpTLead->Fill(LeadHypBJet.Pt(), recoWeight);
    h_RecoBJetpTNLead->Fill(NLeadHypBJet.Pt(), recoWeight);
    h_RecoBJetEtaLead->Fill(LeadHypBJet.Eta(), recoWeight);
    h_RecoBJetEtaNLead->Fill(NLeadHypBJet.Eta(), recoWeight);
    
    //now go to the plots 
    h_HypTTBarMass->Fill(hypttbar.M(), weight);
    h_HypTTBarRapidity->Fill(hypttbar.Rapidity(), weight);
    h_HypTTBarpT->Fill(hypttbar.Pt(), weight);
    h_HypTTBarDeltaPhi->Fill(std::abs(DeltaPhi((*kinRecoObjects.HypTop_).at(solutionIndex), (*kinRecoObjects.HypAntiTop_).at(solutionIndex))), weight);
    h_HypTTBarDeltaRapidity->Fill(std::abs((*kinRecoObjects.HypTop_).at(solutionIndex).Rapidity()) - std::abs((*kinRecoObjects.HypAntiTop_).at(solutionIndex).Rapidity()), weight);

    h_HypBBBarMass->Fill(hypbbbar.M(), weight);
    h_HypBBBarpT->Fill(hypbbbar.Pt(), weight);

    h_HypLLBarMass->Fill(hypllbar.M(), weight);
    h_HypLLBarpT->Fill(hypllbar.Pt(), weight);

    h_HypToppTTTRestFrame->Fill(top.Pt(), weight);
    h_HypAntiToppTTTRestFrame->Fill(antitop.Pt(), weight);

    h_HypMet->Fill(met.Pt(), weight);
    h_HypHT->Fill(jetHT, weight);

    h_HypTopMass->Fill((*kinRecoObjects.HypTop_).at(solutionIndex).M(), weight);
    h_HypAntiTopMass->Fill((*kinRecoObjects.HypAntiTop_).at(solutionIndex).M(), weight);
    h_HypToppT->Fill((*kinRecoObjects.HypTop_).at(solutionIndex).Pt(), weight);
    h_HypAntiToppT->Fill((*kinRecoObjects.HypAntiTop_).at(solutionIndex).Pt(), weight);
        
    h_HypLeptonpT->Fill((*kinRecoObjects.HypLepton_).at(solutionIndex).Pt(), weight);
    h_HypAntiLeptonpT->Fill((*kinRecoObjects.HypAntiLepton_).at(solutionIndex).Pt(), weight);

    h_HypBJetpT->Fill((*kinRecoObjects.HypBJet_).at(solutionIndex).Pt(), weight);
    h_HypAntiBJetpT->Fill((*kinRecoObjects.HypAntiBJet_).at(solutionIndex).Pt(), weight);
    h_HypBJetRapidity->Fill((*kinRecoObjects.HypBJet_).at(solutionIndex).Rapidity(), weight);
    h_HypAntiBJetRapidity->Fill((*kinRecoObjects.HypAntiBJet_).at(solutionIndex).Rapidity(), weight);

    h_HypTopRapidity->Fill((*kinRecoObjects.HypTop_).at(solutionIndex).Rapidity(), weight);
    h_HypAntiTopRapidity->Fill((*kinRecoObjects.HypAntiTop_).at(solutionIndex).Rapidity(), weight);
    
    h_HypNeutrinopT->Fill((*kinRecoObjects.HypNeutrino_).at(solutionIndex).Pt(), weight);
    h_HypAntiNeutrinopT->Fill((*kinRecoObjects.HypAntiNeutrino_).at(solutionIndex).Pt(), weight);
    
    h_HypTopEta->Fill((*kinRecoObjects.HypTop_).at(solutionIndex).Eta(), weight);
    h_HypAntiTopEta->Fill((*kinRecoObjects.HypAntiTop_).at(solutionIndex).Eta(), weight);
    h_HypBJetEta->Fill((*kinRecoObjects.HypBJet_).at(solutionIndex).Eta(), weight);
    h_HypAntiBJetEta->Fill((*kinRecoObjects.HypAntiBJet_).at(solutionIndex).Eta(), weight);
    h_HypLeptonEta->Fill((*kinRecoObjects.HypLepton_).at(solutionIndex).Eta(), weight);

    h_HypAntiLeptonEta->Fill((*kinRecoObjects.HypAntiLepton_).at(solutionIndex).Eta(), weight);

    h_HypLLBarDPhi->Fill(std::fabs ( DeltaPhi ( (*kinRecoObjects.HypLepton_).at(solutionIndex), (*kinRecoObjects.HypAntiLepton_).at(solutionIndex) ) ), weight);
    h_HypLeptonantiBjetMass->Fill(( (*kinRecoObjects.HypLepton_).at(solutionIndex) + (*kinRecoObjects.HypAntiBJet_).at(solutionIndex) ).M(), weight);
    h_HypAntiLeptonBjetMass->Fill(( (*kinRecoObjects.HypAntiLepton_).at(solutionIndex) + (*kinRecoObjects.HypBJet_).at(solutionIndex) ).M(), weight);

    h_HypToppTLead->Fill(LeadHypTop.Pt(), weight);
    h_HypToppTNLead->Fill(NLeadHypTop.Pt(), weight);
    h_HypTopRapidityLead->Fill(LeadHypTop.Rapidity(), weight);
    h_HypTopRapidityNLead->Fill(NLeadHypTop.Rapidity(), weight);
    h_HypTopMassLead->Fill(LeadHypTop.M(), weight);
    h_HypTopMassNLead->Fill(NLeadHypTop.M(), weight);

    h_HypLeptonpTLead->Fill(LeadHypLepton.Pt(), weight);
    h_HypLeptonpTNLead->Fill(NLeadHypLepton.Pt(), weight);
    h_HypLeptonEtaLead->Fill(LeadHypLepton.Eta(), weight);
    h_HypLeptonEtaNLead->Fill(NLeadHypLepton.Eta(), weight);

    h_HypBJetpTLead->Fill(LeadHypBJet.Pt(), weight);
    h_HypBJetpTNLead->Fill(NLeadHypBJet.Pt(), weight);
    h_HypBJetEtaLead->Fill(LeadHypBJet.Eta(), weight);
    h_HypBJetEtaNLead->Fill(NLeadHypBJet.Eta(), weight);

    //Ievgen 
    
    h_HypTTBarRapidityvsTTBarpT->Fill(hypttbar.Pt(),hypttbar.Rapidity(),weight);
    
    // ...
    // Parton momentum fraction
    double RecoPartonMomFraction = ((*kinRecoObjects.HypTop_).at(solutionIndex).energy() - (*kinRecoObjects.HypTop_).at(solutionIndex).Pz() + (*kinRecoObjects.HypAntiTop_).at(solutionIndex).energy() - (*kinRecoObjects.HypAntiTop_).at(solutionIndex).Pz()) / (2 * 4000);
    double RecoAntipartonMomFraction = ((*kinRecoObjects.HypTop_).at(solutionIndex).energy() + (*kinRecoObjects.HypTop_).at(solutionIndex).Pz() + (*kinRecoObjects.HypAntiTop_).at(solutionIndex).energy() + (*kinRecoObjects.HypAntiTop_).at(solutionIndex).Pz()) / (2 * 4000);
    h_HypTopPartonFraction->Fill(RecoPartonMomFraction, weight);
    h_HypAntiTopPartonFraction->Fill(RecoAntipartonMomFraction, weight);
    h_RecoTopPartonFraction->Fill(RecoPartonMomFraction, recoWeight);
    h_RecoAntiTopPartonFraction->Fill(RecoAntipartonMomFraction, recoWeight);

    //New plots from Carmen: Begin
    h_RecoTTBar0Mass->Fill(rho0/(hypttbar).M(), recoWeight);
    h_HypTTBar0Mass->Fill(rho0/(hypttbar).M(), weight);

    int extrarecojet[4] = {0,0,0,0};
    int jetnumReco = -1;
    double jetHTreco = 0;
    int RecoJets = 0, RecoJets_cut40 = 0, RecoJets_cut60 = 0, RecoJets_cut100 = 0;
//    double cbin[19]={35.,45.,55.,65.,75.,85.,95.,110.,130.,150.,170.,190.,210.,230.,250.,270.,300.,340.,380.};
    double cbin[20]={25.,35.,45.,55.,65.,75.,85.,95.,110.,130.,150.,170.,190.,210.,230.,250.,270.,300.,340.,380.};
    for(int k=0; k<(int)(*recoObjects.jets_).size(); k++){
            //if(std::fabs(jets->at(k).Eta())>1.5 || std::fabs(jets->at(k).Eta())<0.8) continue;
        if(std::fabs((*recoObjects.jets_).at(k).Eta()) > 2.4) continue;//carmen eta cuts
//        if((*recoObjects.jets_).at(k).Pt()< JetPtCUT) {continue;}
    // changed to fill Mult***
        if ((*recoObjects.jets_).at(k).Pt()> JetPtCUT2)
        {
            RecoJets++;
            if(std::fabs((*kinRecoObjects.HypAntiBJet_).at(solutionIndex).Pt() - (*recoObjects.jets_).at(k).Pt())>0.1 && std::fabs((*kinRecoObjects.HypBJet_).at(solutionIndex).Pt() - (*recoObjects.jets_).at(k).Pt())>0.1 && jetnumReco<3) {

            jetHTreco+=(*recoObjects.jets_).at(k).Pt();
            jetnumReco++;
            extrarecojet[jetnumReco]= k;
            }
        }
        if((*recoObjects.jets_).at(k).Pt()>40.) RecoJets_cut40++;
        if((*recoObjects.jets_).at(k).Pt()>60.) RecoJets_cut60++;
        if((*recoObjects.jets_).at(k).Pt()>100.) RecoJets_cut100++;
    }

    h_RecoJetMultpt30->Fill(RecoJets,recoWeight);
    h_RecoJetMultpt40->Fill(RecoJets_cut40,recoWeight);
    h_RecoJetMultpt60->Fill(RecoJets_cut60,recoWeight);
    h_RecoJetMultpt100->Fill(RecoJets_cut100,recoWeight);
    h_HypJetMultpt30->Fill(RecoJets,weight);
    h_HypJetMultpt40->Fill(RecoJets_cut40,weight);
    h_HypJetMultpt60->Fill(RecoJets_cut60,weight);
    h_HypJetMultpt100->Fill(RecoJets_cut100,weight);
    
    int first= -1, second=-1, third=-1, fourth=-1;
    double ptjet = 0;
    for(int ord = 0; ord <= jetnumReco; ord++)
    {
        if((*recoObjects.jets_).at(extrarecojet[ord]).Pt()> ptjet) {
            first = ord; ptjet=(*recoObjects.jets_).at(extrarecojet[ord]).Pt();
        }
    }
    ptjet = 0;
    for(int ord = 0; ord <= jetnumReco && jetnumReco>0; ord++)
    {
        if((*recoObjects.jets_).at(extrarecojet[ord]).Pt()> ptjet && (*recoObjects.jets_).at(extrarecojet[ord]).Pt()< (*recoObjects.jets_).at(extrarecojet[first]).Pt()) {
            second = ord; ptjet=(*recoObjects.jets_).at(extrarecojet[ord]).Pt();
        }
    }
    ptjet = 0;
    for(int ord = 0; ord <= jetnumReco && jetnumReco>1; ord++)
    {
        if((*recoObjects.jets_).at(extrarecojet[ord]).Pt()> ptjet && (*recoObjects.jets_).at(extrarecojet[ord]).Pt()< (*recoObjects.jets_).at(extrarecojet[second]).Pt() ) {
            third = ord; ptjet=(*recoObjects.jets_).at(extrarecojet[ord]).Pt();
        }
    }
    ptjet = 0;
    for(int ord = 0; ord <= jetnumReco && jetnumReco>2; ord++)
    {
        if((*recoObjects.jets_).at(extrarecojet[ord]).Pt()> ptjet && (*recoObjects.jets_).at(extrarecojet[ord]).Pt()< (*recoObjects.jets_).at(extrarecojet[third]).Pt() ) {
            fourth = ord; ptjet=(*recoObjects.jets_).at(extrarecojet[ord]).Pt();
        }
    }
    if(jetnumReco>2)
    {
        h_RecoExtraJetpT4->Fill((*recoObjects.jets_).at(extrarecojet[fourth]).Pt(),recoWeight);
        h_RecoExtraJetEta4->Fill((*recoObjects.jets_).at(extrarecojet[fourth]).Eta(),recoWeight);
        h_HypExtraJetpT4->Fill((*recoObjects.jets_).at(extrarecojet[fourth]).Pt(),weight);
        h_HypExtraJetEta4->Fill((*recoObjects.jets_).at(extrarecojet[fourth]).Eta(),weight);
    }
    if (jetnumReco>1)
    {
        h_RecoExtraJetpT3->Fill((*recoObjects.jets_).at(extrarecojet[third]).Pt(),recoWeight);
        h_RecoExtraJetEta3->Fill((*recoObjects.jets_).at(extrarecojet[third]).Eta(),recoWeight); 
        h_HypExtraJetpT3->Fill((*recoObjects.jets_).at(extrarecojet[third]).Pt(),weight);
        h_HypExtraJetEta3->Fill((*recoObjects.jets_).at(extrarecojet[third]).Eta(),weight);
    }
    if (jetnumReco>0)
    {
        h_RecoExtraJetpT2->Fill((*recoObjects.jets_).at(extrarecojet[second]).Pt(),recoWeight);
        h_RecoExtraJetEta2->Fill((*recoObjects.jets_).at(extrarecojet[second]).Eta(),recoWeight);
        h_HypExtraJetpT2->Fill((*recoObjects.jets_).at(extrarecojet[second]).Pt(),weight);
        h_HypExtraJetEta2->Fill((*recoObjects.jets_).at(extrarecojet[second]).Eta(),weight); 
        if(second >-1 && first > -1 )
        {
            h_RecoDeltaRExtraJet12->Fill(std::fabs(DeltaR((*recoObjects.jets_).at(extrarecojet[second]),(*recoObjects.jets_).at(extrarecojet[first]))),recoWeight);
            h_HypDeltaRExtraJet12->Fill(std::fabs(DeltaR((*recoObjects.jets_).at(extrarecojet[second]),(*recoObjects.jets_).at(extrarecojet[first]))),weight);
            h_RecoMassExtraJet12->Fill(((*recoObjects.jets_).at(extrarecojet[second])+(*recoObjects.jets_).at(extrarecojet[first])).M(),weight);
            h_HypMassExtraJet12->Fill(((*recoObjects.jets_).at(extrarecojet[second])+(*recoObjects.jets_).at(extrarecojet[first])).M(),weight);
         if((*recoObjects.jets_).at(extrarecojet[first]).Pt()< 60. && (*recoObjects.jets_).at(extrarecojet[second]).Pt()< 60. && TMath::Abs(DeltaR((*recoObjects.jets_).at(extrarecojet[second]),(*recoObjects.jets_).at(extrarecojet[first]))) > 0.6 && (*recoObjects.jets_).at(extrarecojet[first]).Eta() >0.0 && (*recoObjects.jets_).at(extrarecojet[second]).Eta() < 0.0 ) {
         h_RecoDeltaPhiExtraJet12->Fill(DeltaPhi((*recoObjects.jets_).at(extrarecojet[first]),(*recoObjects.jets_).at(extrarecojet[second])),recoWeight);
         h_HypDeltaPhiExtraJet12->Fill(DeltaPhi((*recoObjects.jets_).at(extrarecojet[first]),(*recoObjects.jets_).at(extrarecojet[second])),weight);
         h_RecoPhiExtraJet12->Fill((*recoObjects.jets_).at(extrarecojet[second]).phi()+(*recoObjects.jets_).at(extrarecojet[first]).phi(),recoWeight);
         h_HypPhiExtraJet12->Fill((*recoObjects.jets_).at(extrarecojet[second]).phi()+(*recoObjects.jets_).at(extrarecojet[first]).phi(),weight);
         }
        }
    }
    if (jetnumReco > -1)
    {
        h_RecoExtraJetpT->Fill((*recoObjects.jets_).at(extrarecojet[first]).Pt(),recoWeight);
        h_RecoExtraJetEta->Fill((*recoObjects.jets_).at(extrarecojet[first]).Eta(),recoWeight);
        h_HypExtraJetpT->Fill((*recoObjects.jets_).at(extrarecojet[first]).Pt(),weight);
        h_HypExtraJetEta->Fill((*recoObjects.jets_).at(extrarecojet[first]).Eta(),weight);

        h_RecoTTBar1stJetMass->Fill(rho0/(hypttbar+(*recoObjects.jets_).at(extrarecojet[first])).M(), recoWeight);
        h_HypTTBar1stJetMass->Fill(rho0/(hypttbar+(*recoObjects.jets_).at(extrarecojet[first])).M(),weight);

    }
    for(int q0 = 0; q0<20;q0++){
        h_RecoJetMultTotal->Fill(cbin[q0],recoWeight);
        h_HypJetMultTotal->Fill(cbin[q0],weight);
        if((first >-1 && (*recoObjects.jets_).at(extrarecojet[first]).Pt()<= cbin[q0] )|| jetnumReco <0) {h_RecoJetMultQ0->Fill(cbin[q0],recoWeight);h_HypJetMultQ0->Fill(cbin[q0],weight);}
        if((second >-1 && (*recoObjects.jets_).at(extrarecojet[second]).Pt()<= cbin[q0]) || jetnumReco <1 ) {h_RecoJetExtra2Q0->Fill(cbin[q0],recoWeight);h_HypJetExtra2Q0->Fill(cbin[q0],weight);}
        if(jetHTreco<=cbin[q0]) {h_RecoJetMultQsum->Fill(cbin[q0],recoWeight); h_HypJetMultQsum->Fill(cbin[q0],weight);}
    }
//     //New plots from Carmen: End

    //make sure you have called CreateBinnedControlPlots in the SlaveBegin first
    for (const auto& i : { (*kinRecoObjects.HypTop_).at(solutionIndex), (*kinRecoObjects.HypAntiTop_).at(solutionIndex) } ) {
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_LeptonpT, (*recoObjects.allLeptons_).at(leptonIndex).Pt(), weight);
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_LeptonpT, (*recoObjects.allLeptons_).at(antiLeptonIndex).Pt(), weight);
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_diLepMassFull, dilepton.M(), weight);
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_LeptonEta, (*recoObjects.allLeptons_).at(leptonIndex).Eta(), weight);
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_LeptonEta, (*recoObjects.allLeptons_).at(antiLeptonIndex).Eta(), weight);
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_MET, met.Pt(), weight);
        
        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_LeptonpT, (*recoObjects.allLeptons_).at(leptonIndex).Pt(), weight);
        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_LeptonpT, (*recoObjects.allLeptons_).at(antiLeptonIndex).Pt(), weight);
        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_diLepMassFull, dilepton.M(), weight);
        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_LeptonEta, (*recoObjects.allLeptons_).at(leptonIndex).Eta(), weight);
        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_LeptonEta, (*recoObjects.allLeptons_).at(antiLeptonIndex).Eta(), weight);
        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_MET, met.Pt(), weight);
    }

    //=== CUT ===
    //Following histograms only filled for the signal sample
    if (!topGenObjects.valuesSet_) return kTRUE;

    // top quark properties
    h_GenRecoToppT->Fill((*kinRecoObjects.HypTop_).at(solutionIndex).Pt(), (*topGenObjects.GenTop_).Pt(), weight );
    h_GenRecoAntiToppT->Fill((*kinRecoObjects.HypAntiTop_).at(solutionIndex).Pt(), (*topGenObjects.GenAntiTop_).Pt(), weight );
    h_GenRecoToppTLead->Fill(LeadHypTop.Pt(), LeadGenTop.Pt(), weight);
    h_GenRecoTopRapidity->Fill((*kinRecoObjects.HypTop_).at(solutionIndex).Rapidity(), (*topGenObjects.GenTop_).Rapidity(), weight );
    h_GenRecoAntiTopRapidity->Fill((*kinRecoObjects.HypAntiTop_).at(solutionIndex).Rapidity(), (*topGenObjects.GenAntiTop_).Rapidity(), weight );
    h_GenRecoToppTNLead->Fill(NLeadHypTop.Pt(), NLeadGenTop.Pt(), weight);
    h_GenRecoTopRapidityLead->Fill(LeadHypTop.Rapidity(), LeadGenTop.Rapidity(), weight);
    h_GenRecoTopRapidityNLead->Fill(NLeadHypTop.Rapidity(), NLeadGenTop.Rapidity(), weight);
    h_GenRecoTopMassLead->Fill(LeadHypTop.M(), LeadGenTop.M(), weight);
    h_GenRecoTopMassNLead->Fill(NLeadHypTop.M(), NLeadGenTop.M(), weight);

    // ttbar properties
    LV genttbar((*topGenObjects.GenTop_) + (*topGenObjects.GenAntiTop_));
    h_GenRecoTTBarMass->Fill(hypttbar.M(), genttbar.M(), weight );
    h_GenRecoTTBarpT->Fill(hypttbar.Pt(), genttbar.Pt(), weight );
    h_GenRecoTTBarRapidity->Fill(hypttbar.Rapidity(), genttbar.Rapidity(), weight );
    h_GenRecoTTBarDeltaPhi->Fill(std::abs(DeltaPhi((*kinRecoObjects.HypTop_).at(solutionIndex), (*kinRecoObjects.HypAntiTop_).at(solutionIndex))),
                                 std::abs(DeltaPhi((*topGenObjects.GenTop_), (*topGenObjects.GenAntiTop_))), weight);
    h_GenRecoTTBarDeltaRapidity->Fill(std::abs((*kinRecoObjects.HypTop_).at(solutionIndex).Rapidity()) - std::abs((*kinRecoObjects.HypAntiTop_).at(solutionIndex).Rapidity()),
                                     std::abs((*topGenObjects.GenTop_).Rapidity()) - std::abs((*topGenObjects.GenAntiTop_).Rapidity()),
                                      weight);

    // create top/antitop quark in the ttbar rest frame
    ROOT::Math::Boost CoMBoostGenTtbar (genttbar.BoostToCM());
    LV gentop ((*topGenObjects.GenTop_));
    LV genantitop ((*topGenObjects.GenAntiTop_));
    gentop = CoMBoostGenTtbar(gentop);
    genantitop = CoMBoostGenTtbar(genantitop);
    h_GenRecoToppTTTRestFrame->Fill(top.Pt(), gentop.Pt(), weight);
    h_GenRecoAntiToppTTTRestFrame->Fill(antitop.Pt(), genantitop.Pt(), weight);

    // Parton momentum fraction as defined by Olaf
    double GenPartonMomFraction = ((*topGenObjects.GenTop_).energy() - (*topGenObjects.GenTop_).Pz() + (*topGenObjects.GenAntiTop_).energy() - (*topGenObjects.GenAntiTop_).Pz()) / (2 * 4000);
    double GenAntipartonMomFraction = ((*topGenObjects.GenTop_).energy() + (*topGenObjects.GenTop_).Pz() + (*topGenObjects.GenAntiTop_).energy() + (*topGenObjects.GenAntiTop_).Pz()) / (2 * 4000);
    h_GenRecoTopPartonFraction->Fill(RecoPartonMomFraction,GenPartonMomFraction, weight);
    h_GenRecoAntiTopPartonFraction->Fill(RecoAntipartonMomFraction, GenAntipartonMomFraction, weight);


    // fill object in visible phase-space
    double noDeltaRleptjet = -0.1;
    if(!isttdiffXS) noDeltaRleptjet = 0.4; 

    if( (*topGenObjects.GenLepton_).Pt()> LeptonPtCut && std::fabs((*topGenObjects.GenLepton_).Eta()) <LeptonEtaCUT &&
       (*topGenObjects.GenAntiLepton_).Pt()> LeptonPtCut && std::fabs((*topGenObjects.GenAntiLepton_).Eta()) <LeptonEtaCUT &&
       BHadronIndex >= 0 && AntiBHadronIndex >= 0 && (*topGenObjects.allGenJets_).at(BHadronIndex).pt() > JetPtCUT && std::fabs ( (*topGenObjects.allGenJets_).at(BHadronIndex).eta() ) < JetEtaCUT && (*topGenObjects.allGenJets_).at(AntiBHadronIndex).pt() > JetPtCUT && std::fabs ( (*topGenObjects.allGenJets_).at(AntiBHadronIndex).Eta() ) < JetEtaCUT && std::fabs(DeltaR((*topGenObjects.GenLepton_), (*topGenObjects.allGenJets_).at(BHadronIndex))) > noDeltaRleptjet  && std::fabs(DeltaR((*topGenObjects.GenLepton_), (*topGenObjects.allGenJets_).at(AntiBHadronIndex))) > noDeltaRleptjet  && std::fabs(DeltaR((*topGenObjects.GenAntiLepton_), (*topGenObjects.allGenJets_).at(BHadronIndex))) > noDeltaRleptjet  && std::fabs(DeltaR((*topGenObjects.GenAntiLepton_), (*topGenObjects.allGenJets_).at(AntiBHadronIndex))) > noDeltaRleptjet) 
            {
            // extra objects: met, ht, ...
            h_GenRecoMet->Fill(met.Pt(), (*topGenObjects.GenMet_).Pt(), weight);
            h_GenRecoHT->Fill(jetHT, genHT, weight);
            h_GenRecoNeutrinopT->Fill((*kinRecoObjects.HypNeutrino_).at(solutionIndex).Pt(), (*topGenObjects.GenNeutrino_).Pt(), weight);
            h_GenRecoAntiNeutrinopT->Fill((*kinRecoObjects.HypAntiNeutrino_).at(solutionIndex).Pt(), (*topGenObjects.GenAntiNeutrino_).Pt(), weight);

            // lepton distributions
            h_GenRecoLeptonpT->Fill((*kinRecoObjects.HypLepton_).at(solutionIndex).Pt(), (*topGenObjects.GenLepton_).Pt(), weight );
            h_GenRecoAntiLeptonpT->Fill((*kinRecoObjects.HypAntiLepton_).at(solutionIndex).Pt(), (*topGenObjects.GenAntiLepton_).Pt(), weight );
            h_GenRecoLeptonEta->Fill((*kinRecoObjects.HypLepton_).at(solutionIndex).Eta(), (*topGenObjects.GenLepton_).Eta(), weight );
            h_GenRecoAntiLeptonEta->Fill((*kinRecoObjects.HypAntiLepton_).at(solutionIndex).Eta(), (*topGenObjects.GenAntiLepton_).Eta(), weight );
            h_GenRecoLeptonpTLead->Fill(LeadHypLepton.Pt(), LeadGenLepton.Pt(), weight);
            h_GenRecoLeptonpTNLead->Fill(NLeadHypLepton.Pt(), NLeadGenLepton.Pt(), weight);
            h_GenRecoLeptonEtaLead->Fill(LeadHypLepton.Eta(), LeadGenLepton.Eta(), weight);
            h_GenRecoLeptonEtaNLead->Fill(NLeadHypLepton.Eta(), NLeadGenLepton.Eta(), weight);

            // lepton-pair distributions
            LV genllbar((*topGenObjects.GenLepton_) + (*topGenObjects.GenAntiLepton_));
            h_GenRecoLLBarMass->Fill(hypllbar.M(), genllbar.M(), weight );
            h_GenRecoLLBarpT->Fill(hypllbar.Pt(), genllbar.Pt(), weight );
            h_GenRecoLLBarDPhi->Fill(
                std::fabs( DeltaPhi( (*kinRecoObjects.HypLepton_).at(solutionIndex), (*kinRecoObjects.HypAntiLepton_).at(solutionIndex) ) ),
                std::fabs( DeltaPhi( (*topGenObjects.GenLepton_), (*topGenObjects.GenAntiLepton_) ) ),
                weight );

            // letpon-b-jet mass
            h_GenRecoLeptonantiBjetMass->Fill(( (*kinRecoObjects.HypLepton_).at(solutionIndex) + (*kinRecoObjects.HypAntiBJet_).at(solutionIndex) ).M(), ( (*topGenObjects.GenLepton_)+(*topGenObjects.allGenJets_).at(AntiBHadronIndex) ).M(), weight );
            h_GenRecoAntiLeptonBjetMass->Fill(( (*kinRecoObjects.HypAntiLepton_).at(solutionIndex)+(*kinRecoObjects.HypBJet_).at(solutionIndex) ).M(), ( (*topGenObjects.GenAntiLepton_)+(*topGenObjects.allGenJets_).at(BHadronIndex) ).M(), weight );

            // b-jet-pair distributions
            LV genbbbar ((*topGenObjects.allGenJets_).at(BHadronIndex) + (*topGenObjects.allGenJets_).at(AntiBHadronIndex));
            h_GenRecoBBBarpT->Fill(hypbbbar.Pt(), genbbbar.Pt(), weight);
            h_GenRecoBBBarMass->Fill(hypbbbar.M(), genbbbar.M(), weight);

            // bjet distributions
            h_GenRecoBJetpT->Fill((*kinRecoObjects.HypBJet_).at(solutionIndex).Pt(), (*topGenObjects.allGenJets_).at(BHadronIndex).Pt(), weight );
            h_GenRecoAntiBJetpT->Fill((*kinRecoObjects.HypAntiBJet_).at(solutionIndex).Pt(), (*topGenObjects.allGenJets_).at(AntiBHadronIndex).Pt(), weight );
            h_GenRecoBJetEta->Fill((*kinRecoObjects.HypBJet_).at(solutionIndex).Eta(), (*topGenObjects.allGenJets_).at(BHadronIndex).Eta(), weight );
            h_GenRecoAntiBJetEta->Fill((*kinRecoObjects.HypAntiBJet_).at(solutionIndex).Eta(), (*topGenObjects.allGenJets_).at(AntiBHadronIndex).Eta(), weight );
            h_GenRecoBJetRapidity->Fill((*kinRecoObjects.HypBJet_).at(solutionIndex).Rapidity(), (*topGenObjects.allGenJets_).at(BHadronIndex).Rapidity(), weight );
            h_GenRecoAntiBJetRapidity->Fill((*kinRecoObjects.HypAntiBJet_).at(solutionIndex).Rapidity(), (*topGenObjects.allGenJets_).at(AntiBHadronIndex).Rapidity(), weight );
            h_GenRecoBJetpTLead->Fill(LeadHypBJet.Pt(), LeadGenBJet.Pt(), weight);
            h_GenRecoBJetpTNLead->Fill(NLeadHypBJet.Pt(), NLeadGenBJet.Pt(), weight);
            h_GenRecoBJetEtaLead->Fill(LeadHypBJet.Eta(), LeadGenBJet.Eta(), weight);
            h_GenRecoBJetEtaNLead->Fill(NLeadHypBJet.Eta(), NLeadGenBJet.Eta(), weight);

            //rho/masstt 
            h_GenRecoTTBar0Mass->Fill(rho0/hypttbar.M(),rho0/genttbar.M(),weight);
            // jet multiplicities
            h_GenRecoJetMult->Fill(numberOfJets, (*topGenObjects.allGenJets_).size(), weight );
            h_GenRecoJetMultpt30->Fill(RecoJets,GenJets_cut,weight);
            h_GenRecoJetMultpt40->Fill(RecoJets_cut40,GenJets_cut40,weight);
            h_GenRecoJetMultpt60->Fill(RecoJets_cut60,GenJets_cut60,weight);
            h_GenRecoJetMultpt100->Fill(RecoJets_cut100,GenJets_cut100,weight);
            double ptaddgen = -1000.;
            double etaaddgen = -1000.;
            double DeltaRgen = -1000.;
            double massgen = -1000.;
            if(jetnumReco>2){
                if(jetnum >2) {ptaddgen = (*topGenObjects.allGenJets_).at(extragenjet[3]).Pt(); etaaddgen = (*topGenObjects.allGenJets_).at(extragenjet[3]).Eta();}
                h_GenRecoExtraJetpT4->Fill((*recoObjects.jets_).at(extrarecojet[fourth]).Pt(),ptaddgen,weight);
                h_GenRecoExtraJetEta4->Fill((*recoObjects.jets_).at(extrarecojet[fourth]).Eta(),etaaddgen,weight);
            }
            if (jetnumReco>1){
                if(jetnum >1) {ptaddgen = (*topGenObjects.allGenJets_).at(extragenjet[2]).Pt(); etaaddgen = (*topGenObjects.allGenJets_).at(extragenjet[2]).Eta();}
                h_GenRecoExtraJetpT3->Fill((*recoObjects.jets_).at(extrarecojet[third]).Pt(),ptaddgen,weight);
                h_GenRecoExtraJetEta3->Fill((*recoObjects.jets_).at(extrarecojet[third]).Eta(),etaaddgen,weight);
            }
            if (jetnumReco>0){
                if(jetnum >0) {
                ptaddgen = (*topGenObjects.allGenJets_).at(extragenjet[1]).Pt(); etaaddgen = (*topGenObjects.allGenJets_).at(extragenjet[1]).Eta();
                DeltaRgen = std::fabs(DeltaR((*topGenObjects.allGenJets_).at(extragenjet[1]),(*topGenObjects.allGenJets_).at(extragenjet[0])));               
                massgen = ((*topGenObjects.allGenJets_).at(extragenjet[1])+(*topGenObjects.allGenJets_).at(extragenjet[0])).M();
                }
                h_GenRecoExtraJetpT2->Fill((*recoObjects.jets_).at(extrarecojet[second]).Pt(),ptaddgen,weight);
                h_GenRecoExtraJetEta2->Fill((*recoObjects.jets_).at(extrarecojet[second]).Eta(),etaaddgen,weight);
                h_GenRecoDeltaRExtraJet12->Fill(std::fabs(DeltaR((*recoObjects.jets_).at(extrarecojet[second]),(*recoObjects.jets_).at(extrarecojet[first]))),DeltaRgen,weight);
                h_GenRecoMassExtraJet12->Fill(((*recoObjects.jets_).at(extrarecojet[second])+(*recoObjects.jets_).at(extrarecojet[first])).M(),massgen, weight);
                if((*recoObjects.jets_).at(extrarecojet[first]).Pt()< 60. && (*recoObjects.jets_).at(extrarecojet[first]).Eta() >0.0 &&
                (*recoObjects.jets_).at(extrarecojet[second]).Pt()< 60. && (*recoObjects.jets_).at(extrarecojet[second]).Eta() < 0.0 &&
                TMath::Abs(DeltaR((*recoObjects.jets_).at(extrarecojet[second]),(*recoObjects.jets_).at(extrarecojet[first]))) > 0.6 )
                {
                    h_GenRecoDeltaPhiExtraJet12->Fill(DeltaPhi((*recoObjects.jets_).at(extrarecojet[first]),(*recoObjects.jets_).at(extrarecojet[second])),DeltaPhi((*topGenObjects.allGenJets_).at(extragenjet[1]),(*topGenObjects.allGenJets_).at(extragenjet[0])),weight);
                    h_GenRecoPhiExtraJet12->Fill((*recoObjects.jets_).at(extrarecojet[first]).Phi()+(*recoObjects.jets_).at(extrarecojet[second]).Phi(),(*topGenObjects.allGenJets_).at(extragenjet[0]).Phi()+(*topGenObjects.allGenJets_).at(extragenjet[1]).Phi(),weight);
                }
            }
            if (jetnumReco >-1){
//               if(jetnum > -1 && std::fabs(DeltaR((*topGenObjects.allGenJets_).at(extrarecojet[0]),(*topGenObjects.allGenJets_).at(extragenjet[0])))<0.5 ) {
               if(jetnum > -1) {
                    ptaddgen = (*topGenObjects.allGenJets_).at(extragenjet[0]).Pt(); etaaddgen = (*topGenObjects.allGenJets_).at(extragenjet[0]).Eta();
                    h_GenRecoTTBar1stJetMass->Fill(rho0/(hypttbar+(*recoObjects.jets_).at(extrarecojet[first])).M(),rho0/(genttbar+(*topGenObjects.allGenJets_).at(extragenjet[0])).M(),weight);
               } else h_GenRecoTTBar1stJetMass->Fill(rho0/(hypttbar+(*recoObjects.jets_).at(extrarecojet[first])).M(),-1000.,weight);

                h_GenRecoExtraJetpT->Fill((*recoObjects.jets_).at(extrarecojet[first]).Pt(),ptaddgen,weight); 
                h_GenRecoExtraJetEta->Fill((*recoObjects.jets_).at(extrarecojet[first]).Eta(),etaaddgen,weight);
//                h_GenRecoTTBar1stJetMass->Fill(rho0/(hypttbar+(*recoObjects.jets_).at(extrarecojet[first])).M(),rho0/(genttbar+(*topGenObjects.allGenJets_).at(extragenjet[0])).M(),weight);
            }
    }
    else{// fill underflow/overflow for reco objects not in vis. phase space
        // extra objects: met, ht, ...
        h_GenRecoMet->Fill(met.Pt(), -1000, weight);
        h_GenRecoHT->Fill(jetHT, -1000, weight);
        h_GenRecoNeutrinopT->Fill((*kinRecoObjects.HypNeutrino_).at(solutionIndex).Pt(), -1000, weight);
        h_GenRecoAntiNeutrinopT->Fill((*kinRecoObjects.HypAntiNeutrino_).at(solutionIndex).Pt(), -1000, weight);

        // lepton distributions
        h_GenRecoLeptonpT->Fill((*kinRecoObjects.HypLepton_).at(solutionIndex).Pt(), -1000, weight );
        h_GenRecoAntiLeptonpT->Fill((*kinRecoObjects.HypAntiLepton_).at(solutionIndex).Pt(), -1000, weight );
        h_GenRecoLeptonEta->Fill((*kinRecoObjects.HypLepton_).at(solutionIndex).Eta(), -1000, weight );
        h_GenRecoAntiLeptonEta->Fill((*kinRecoObjects.HypAntiLepton_).at(solutionIndex).Eta(), -1000, weight );
        h_GenRecoLeptonpTLead->Fill(LeadHypLepton.Pt(), -1000, weight);
        h_GenRecoLeptonpTNLead->Fill(NLeadHypLepton.Pt(), -1000, weight);
        h_GenRecoLeptonEtaLead->Fill(LeadHypLepton.Eta(), -1000, weight);
        h_GenRecoLeptonEtaNLead->Fill(NLeadHypLepton.Eta(), -1000, weight);

        // lepton-pair distributions
        h_GenRecoLLBarMass->Fill(hypllbar.M(), -1000, weight );
        h_GenRecoLLBarpT->Fill(hypllbar.Pt(), -1000, weight );
        h_GenRecoLLBarDPhi->Fill(
            std::fabs( DeltaPhi( (*kinRecoObjects.HypLepton_).at(solutionIndex), (*kinRecoObjects.HypAntiLepton_).at(solutionIndex) ) ), -1000, weight );

        // letpon-b-jet mass
        h_GenRecoLeptonantiBjetMass->Fill(( (*kinRecoObjects.HypLepton_).at(solutionIndex) + (*kinRecoObjects.HypAntiBJet_).at(solutionIndex) ).M(), -1000, weight );
        h_GenRecoAntiLeptonBjetMass->Fill(( (*kinRecoObjects.HypAntiLepton_).at(solutionIndex)+(*kinRecoObjects.HypBJet_).at(solutionIndex) ).M(), -1000, weight );

        // b-jet-pair distributions
        h_GenRecoBBBarpT->Fill(hypbbbar.Pt(), -1000, weight);
        h_GenRecoBBBarMass->Fill(hypbbbar.M(), -1000, weight);

        // bjet distributions
        h_GenRecoBJetpT->Fill((*kinRecoObjects.HypBJet_).at(solutionIndex).Pt(), -1000, weight );
        h_GenRecoAntiBJetpT->Fill((*kinRecoObjects.HypAntiBJet_).at(solutionIndex).Pt(), -1000, weight );
        h_GenRecoBJetEta->Fill((*kinRecoObjects.HypBJet_).at(solutionIndex).Eta(), -1000, weight );
        h_GenRecoAntiBJetEta->Fill((*kinRecoObjects.HypAntiBJet_).at(solutionIndex).Eta(), -1000, weight );
        h_GenRecoBJetRapidity->Fill((*kinRecoObjects.HypBJet_).at(solutionIndex).Rapidity(), -1000, weight );
        h_GenRecoAntiBJetRapidity->Fill((*kinRecoObjects.HypAntiBJet_).at(solutionIndex).Rapidity(), -1000, weight );
        h_GenRecoBJetpTLead->Fill(LeadHypBJet.Pt(), -1000, weight);
        h_GenRecoBJetpTNLead->Fill(NLeadHypBJet.Pt(), -1000, weight);
        h_GenRecoBJetEtaLead->Fill(LeadHypBJet.Eta(), -1000, weight);
        h_GenRecoBJetEtaNLead->Fill(NLeadHypBJet.Eta(), -1000, weight);

        //rho/masstt 
        h_GenRecoTTBar0Mass->Fill(rho0/hypttbar.M(), -1000, weight);

        // jet multiplicities
        h_GenRecoJetMult->Fill(numberOfJets, -1000, weight );
        h_GenRecoJetMultpt30->Fill(RecoJets,-1000,weight);
        h_GenRecoJetMultpt40->Fill(RecoJets_cut40,-1000,weight);
        h_GenRecoJetMultpt60->Fill(RecoJets_cut60,-1000,weight);
        h_GenRecoJetMultpt100->Fill(RecoJets_cut100,-1000,weight);
        if(jetnumReco>2){
            h_GenRecoExtraJetpT4->Fill((*recoObjects.jets_).at(extrarecojet[fourth]).Pt(),-1000.,weight);
            h_GenRecoExtraJetEta4->Fill((*recoObjects.jets_).at(extrarecojet[fourth]).Eta(),-1000.,weight);
        }
        if(jetnumReco>1){
            h_GenRecoExtraJetpT3->Fill((*recoObjects.jets_).at(extrarecojet[third]).Pt(),-1000.,weight);
            h_GenRecoExtraJetEta3->Fill((*recoObjects.jets_).at(extrarecojet[third]).Eta(),-1000.,weight);
        }
        if (jetnumReco>0){
            h_GenRecoExtraJetpT2->Fill((*recoObjects.jets_).at(extrarecojet[second]).Pt(),-1000.,weight);
            h_GenRecoExtraJetEta2->Fill((*recoObjects.jets_).at(extrarecojet[second]).Eta(),-1000.,weight);
            h_GenRecoDeltaRExtraJet12->Fill(std::fabs(DeltaR((*recoObjects.jets_).at(extrarecojet[second]),(*recoObjects.jets_).at(extrarecojet[first]))),-1000.,weight);
            h_GenRecoMassExtraJet12->Fill(((*recoObjects.jets_).at(extrarecojet[second])+(*recoObjects.jets_).at(extrarecojet[first])).M(),-1000.,weight);
            h_GenRecoDeltaPhiExtraJet12->Fill(DeltaPhi((*recoObjects.jets_).at(extrarecojet[first]),(*recoObjects.jets_).at(extrarecojet[second])),-1000,weight);
            h_GenRecoPhiExtraJet12->Fill((*recoObjects.jets_).at(extrarecojet[first]).Phi()+(*recoObjects.jets_).at(extrarecojet[second]).Phi(),-1000,weight);
        }
        if (jetnumReco >-1 ){
            h_GenRecoExtraJetpT->Fill((*recoObjects.jets_).at(extrarecojet[first]).Pt(),-1000.,weight);
            h_GenRecoExtraJetEta->Fill((*recoObjects.jets_).at(extrarecojet[first]).Eta(),-1000.,weight);
            h_GenRecoTTBar1stJetMass->Fill(rho0/(hypttbar+(*recoObjects.jets_).at(extrarecojet[first])).M(),-1000.,weight);
        }
    }

    // ...

    return kTRUE;
}



void TopAnalysis::SetPDF(int pdf_no)
{
    this->pdf_no_ = pdf_no;
}



void TopAnalysis::SetClosureTest(TString closure, double slope)
{
    if (closure == "") {
        doClosureTest_ = false;
    } else {
        doClosureTest_ = true;
        if (closure == "pttop") {
            closureFunction_ = [&,slope](Long64_t entry) -> double {
                const TopGenObjects& topGenObjects = this->getTopGenObjects(entry);
                return std::max((1.+(topGenObjects.GenTop_->Pt()-100.)*slope)
                               *(1.+(topGenObjects.GenAntiTop_->Pt()-100.)*slope) , 0.1);
            };
        } else if (closure == "ytop") {
            closureFunction_ = [&,slope](Long64_t entry) -> double {
                const TopGenObjects& topGenObjects = this->getTopGenObjects(entry);
                return std::max((1.+(std::fabs(topGenObjects.GenTop_->Rapidity())-1.)*slope)
                               *(1.+(std::fabs(topGenObjects.GenAntiTop_->Rapidity()-1.))*slope) , 0.1);
            };
        } else if (closure == "nominal") {
            closureFunction_ = [](Long64_t) -> double {return 1.;};
        } else {
            std::cerr << "invalid closure test function\n";
            exit(1);
        }
        
        TString outputFilename = this->outputFilename();
        if (closure != "nominal") {
            outputFilename.ReplaceAll(".root", TString::Format("_fakerun_%s%.3f.root", closure.Data(), slope));
        } else {
            outputFilename.ReplaceAll(".root", TString::Format("_fakerun_%s.root", closure.Data()));
        }
        this->SetOutputfilename(outputFilename);
        std::cout << "<<< Closure test. Writing to: " << outputFilename << "\n";
        
        //BRANCHING FRACTION
        double br = 0.;
        const int channelPdgIdProduct = this->channelPdgIdProduct();
        if (channelPdgIdProduct == -11*11 || channelPdgIdProduct == -13*13) br = 0.01166;
        else if (channelPdgIdProduct == -11*13) br = 0.02332;
        else {
            std::cerr << "closure test channel invalid\n"; exit(1);
        }
        closureMaxEvents_ = TOPXSEC * 1000 * LUMI * br;
        TString samplename = this->samplename();
        samplename.Append("_fakedata");
        this->SetSamplename(samplename);
    }
}



double TopAnalysis::calculateClosureTestWeight(const Long64_t& entry)
{
    if(!doClosureTest_ || !this->isTopSignal()) return 1.;
    const double weight = closureFunction_(entry);
    h_ClosureTotalWeight->Fill(1, weight);
    return weight;
}



double TopAnalysis::globalNormalisationFactorClosureTest()
{
    double globalNormalisationFactor(1);
    if (doClosureTest_) {
        TH1 *total = dynamic_cast<TH1*>(fOutput->FindObject("ClosureTotalWeight"));
        if (!total) {
            std::cerr << "ClosureTotalWeight histogram is missing!\n"; exit(1);
        }
        globalNormalisationFactor *= total->GetEntries() / total->GetBinContent(1);
        std::cout << "Global normalisation factor of closure test: " << globalNormalisationFactor << "\n";
    }
    return globalNormalisationFactor;
}



double TopAnalysis::globalNormalisationFactorPDF()
{
    double globalNormalisationFactor(1);
    if (pdf_no_ >= 0) {
        TH1 *total = dynamic_cast<TH1*>(fOutput->FindObject("PDFTotalWeight"));
        if (!total) {
            std::cerr << "PDFTotalWeight histogram is missing!\n"; exit(1);
        }
        globalNormalisationFactor *= total->GetEntries() / total->GetBinContent(1);
        std::cout << "PDF Weight Normalisation = " << globalNormalisationFactor << "\n";
    }
    return globalNormalisationFactor;
}



double TopAnalysis::overallGlobalNormalisationFactor()
{
    double globalNormalisationFactor(1);
    globalNormalisationFactor *= globalNormalisationFactorClosureTest();
    globalNormalisationFactor *= globalNormalisationFactorPDF();
    return globalNormalisationFactor;
}



void TopAnalysis::bHadronIndices(int& bHadronIndex, int& antiBHadronIndex, const TopGenObjects& topGenObjects)
{
    int& BHadronIndex(bHadronIndex);
    int& AntiBHadronIndex(antiBHadronIndex);
    BHadronIndex = -1;
    AntiBHadronIndex = -1;
    
    if(!topGenObjects.valuesSet_) return;
    
    
    std::vector<size_t> idx_leadbHadJet;
    std::vector<size_t> idx_nleadbHadJet;
    //To avoid recopying may code lines, we select HERE the BHadron JET Indices to cut on.

    //time to choose which genJet we actually want

    idx_leadbHadJet.insert ( idx_leadbHadJet.begin(), 4, -1 );
    idx_nleadbHadJet.insert ( idx_nleadbHadJet.begin(), 4, -1 );
    /*
      idx_bHadJet will have 4 jet indices
      [0] is the highest pT jet with a B-Hadron
      [1] is the highest pT jet with a B-Hadron also matched to a top quark
      [2] highest pT jet of those matched closest (in DeltaR) to the B-Hadron
      [3] highest pT jet of those matched closest (in DeltaR) to the B-Hadron also matched to a top quark
    */

    bool LeadBHadhighpTjet = false;
    bool LeadBHadhighpTjetfromtop = false;
    bool NLeadBHadhighpTjet = false;
    bool NLeadBHadhighpTjetfromtop = false;

    size_t hadron_index = std::numeric_limits<size_t>::max();
    size_t antihadron_index = std::numeric_limits<size_t>::max();
    size_t hadrontop_index = std::numeric_limits<size_t>::max();
    size_t antihadrontop_index = std::numeric_limits<size_t>::max();
    
    //Case 1: highest pT genJet matched to a BHadron
    //need to remove jets from the genJetCollection which are below the JETPTCUT
    //while (allGenJets->size() > 0 && allGenJets->back().Pt() < JETPTCUT) allGenJets->pop_back();
    //while (jet->size() > 0 && jet->back().Pt() < JETPTCUT) jet->pop_back();
    
    for ( size_t genJet = 0;
          genJet < (*topGenObjects.allGenJets_).size() && (*topGenObjects.allGenJets_).at(genJet).pt() >= JetPtCUT && std::fabs((*topGenObjects.allGenJets_).at(genJet).eta()) < JetEtaCUT; 
          ++genJet )
    {
        
        for ( size_t bHadron=0; bHadron < (*topGenObjects.BHadrons_).size(); bHadron++ ) {
            if ( (*topGenObjects.BHadronVsJet_)[genJet*(*topGenObjects.BHadrons_).size()+bHadron]==1 &&
                 (!LeadBHadhighpTjet || !LeadBHadhighpTjetfromtop || !NLeadBHadhighpTjet || !NLeadBHadhighpTjetfromtop) )
            {
                if ( LeadBHadhighpTjet==false ) {
                    idx_leadbHadJet[0] = genJet;
                    LeadBHadhighpTjet = true;
                    hadron_index = bHadron;
                    if ( ( *topGenObjects.BHadronFromTopB_ ) [bHadron] == true ) {
                        idx_leadbHadJet[1] = genJet;
                        LeadBHadhighpTjetfromtop = true;
                        hadrontop_index = bHadron;
                    }
                } else if ( LeadBHadhighpTjetfromtop == false ) {
                    if ( ( *topGenObjects.BHadronFromTopB_ ) [bHadron] == true ) {
                        idx_leadbHadJet[1] = genJet;
                        LeadBHadhighpTjetfromtop = true;
                        hadrontop_index = bHadron;
                    }
                } else if ( NLeadBHadhighpTjet==false && bHadron!=hadron_index && idx_leadbHadJet[0] != genJet ) {
                    idx_nleadbHadJet[0] = genJet;
                    NLeadBHadhighpTjet = true;
                    if ( ( *topGenObjects.BHadronFromTopB_ ) [bHadron] == true && bHadron!=hadrontop_index && idx_leadbHadJet[1] != genJet ) {
                        idx_nleadbHadJet[1] = genJet;
                        NLeadBHadhighpTjetfromtop = true;
                    }
                } else if ( NLeadBHadhighpTjetfromtop == false && bHadron!=hadrontop_index && idx_leadbHadJet[1] != genJet ) {
                    if ( ( *topGenObjects.BHadronFromTopB_ ) [bHadron] == true ) {
                        idx_nleadbHadJet[1] = genJet;
                        LeadBHadhighpTjetfromtop = true;
                    }
                }//series of if statements to find highest pT jet
            }
        }
        for ( size_t antibHadron=0; antibHadron < (*topGenObjects.AntiBHadrons_).size(); antibHadron++ ) {
            if ( (*topGenObjects.AntiBHadronVsJet_)[genJet*(*topGenObjects.AntiBHadrons_).size()+antibHadron]==1 &&
                 (LeadBHadhighpTjet ==false || LeadBHadhighpTjetfromtop == false || NLeadBHadhighpTjet ==false || NLeadBHadhighpTjetfromtop == false ) &&
                 idx_leadbHadJet[0] != genJet ) 
            {
                if ( LeadBHadhighpTjet==false ) {
                    idx_leadbHadJet[0] = genJet;
                    LeadBHadhighpTjet = true;
                    antihadron_index = antibHadron;
                    if ( ( *topGenObjects.AntiBHadronFromTopB_ ) [antibHadron] == true ) {
                        idx_leadbHadJet[1] = genJet;
                        LeadBHadhighpTjetfromtop = true;
                        antihadrontop_index = antibHadron;
                    }
                } else if ( LeadBHadhighpTjetfromtop == false ) {
                    if ( ( *topGenObjects.AntiBHadronFromTopB_ ) [antibHadron] == true ) {
                        idx_leadbHadJet[1] = genJet;
                        LeadBHadhighpTjetfromtop = true;
                        antihadrontop_index = antibHadron;
                    }
                } else if ( NLeadBHadhighpTjet==false && antibHadron!=antihadron_index && idx_leadbHadJet[0] != genJet ) {
                    idx_nleadbHadJet[0] = genJet;
                    NLeadBHadhighpTjet = true;
                    if ( ( *topGenObjects.AntiBHadronFromTopB_ ) [antibHadron] == true && antibHadron!=antihadrontop_index && idx_leadbHadJet[1] != genJet ) {
                        idx_nleadbHadJet[1] = genJet;
                        NLeadBHadhighpTjetfromtop = true;
                    }
                } else if ( NLeadBHadhighpTjetfromtop == false && antibHadron!=antihadrontop_index && idx_leadbHadJet[1] != genJet ) {
                    if ( ( *topGenObjects.AntiBHadronFromTopB_ ) [antibHadron] == true ) {
                        idx_nleadbHadJet[1] = genJet;
                        LeadBHadhighpTjetfromtop = true;
                    }
                }
            }
        }
    }

    //Case 2: highest pT genJets matched closest to a BHadron
    //BHadJetIndex: vector containing the GetJet indices matched, in DeltaR, to a BHadron. Starting from the highest pT jet.
    if ( (*topGenObjects.BHadJetIndex_).size() != 0 ) idx_leadbHadJet[2] = ( *topGenObjects.BHadJetIndex_ ) [0];
    for ( size_t i=0; i < (*topGenObjects.BHadJetIndex_).size(); ++i ) {
        //Only search for those jets matched in DeltaR with a BHadron
        for ( size_t j=0; j<(*topGenObjects.BHadrons_).size() ; ++j ) {
            if ( ( *topGenObjects.BHadronVsJet_ ) [ ( ( *topGenObjects.BHadJetIndex_ ) [i] ) * (*topGenObjects.BHadrons_).size()+j] == 1 && ( *topGenObjects.BHadronFromTopB_ ) [j] == true ) {
                idx_leadbHadJet[3] = ( *topGenObjects.BHadJetIndex_ ) [i];
            }
        }
    }

    //AntiBHadJetIndex: vector containing the GetJet indices matched, in DeltaR, to a AntiBHadron. Starting from the highest pT jet.
    if ( (*topGenObjects.AntiBHadJetIndex_).size() != 0 ) idx_nleadbHadJet[2] = ( *topGenObjects.AntiBHadJetIndex_ ) [0];
    for ( size_t i=0; i < (*topGenObjects.AntiBHadJetIndex_).size(); ++i ) {
        //Only search for those jets matched in DeltaR with a AntiBHadron
        for ( size_t j=0; j < (*topGenObjects.AntiBHadrons_).size() ; ++j ) {
            //if ((*AntiBHadronVsJet)[i*AntiBHadrons_+j] == 1 && (*AntiBHadronFromTopB)[j] == true) {idx_antibHadJet[3] = (*AntiBHadJetIndex)[i];}
            if ( ( *topGenObjects.AntiBHadronVsJet_ ) [ ( ( *topGenObjects.AntiBHadJetIndex_ ) [i] ) * (*topGenObjects.AntiBHadrons_).size()+j] == 1 && ( *topGenObjects.AntiBHadronFromTopB_ ) [j] == true ) {
                idx_nleadbHadJet[3] = ( *topGenObjects.AntiBHadJetIndex_ ) [i];
            }
        }
    }
    
    
    //To avoid recopying many code lines, we select HERE the BHadron JET Indices to cut on.
    //int BHadronIndex;
    //int AntiBHadronIndex;
    //Case 1A: highest pT genJet matched to a BHadron
    BHadronIndex = idx_leadbHadJet[0];
    AntiBHadronIndex = idx_nleadbHadJet[0];
    //
    //   //Case 1B: highest pT genJet matched to a BHadron from Top
    //BHadronIndex = idx_bHadJet[1];
    //AntiBHadronIndex = idx_antibHadJet[1];
    //
    //   //Case 2A: highest pT genJets matched closest to a BHadron
    //BHadronIndex = idx_bHadJet[2];
    //AntiBHadronIndex = idx_antibHadJet[2];
    //
    //   //Case 2B: highest pT genJets matched closest to a BHadron from Top
    //    BHadronIndex = idx_leadbHadJet[3];
    //AntiBHadronIndex = idx_nleadbHadJet[3];
}

void TopAnalysis::generatorTopEvent(LV& leadGenTop, LV& nLeadGenTop,
                                    LV& leadGenLepton, LV& nLeadGenLepton,
                                    LV& leadGenBJet, LV& nLeadGenBJet,
                                    double& genHT,
                                    const int bHadronIndex, const int antiBHadronIndex,
                                    const double trueLevelWeightNoPileup, const double trueLevelWeight, const TopGenObjects& topGenObjects)
{
    // Use utilities without namespaces
    using namespace common;
    using ROOT::Math::VectorUtil::DeltaPhi;
    
    LV& LeadGenTop(leadGenTop);
    LV& NLeadGenTop(nLeadGenTop);
    LV& LeadGenLepton(leadGenLepton);
    LV& NLeadGenLepton(nLeadGenLepton);
    LV& LeadGenBJet(leadGenBJet);
    LV& NLeadGenBJet(nLeadGenBJet);
    
    genHT = -1.;
    
    if(!topGenObjects.valuesSet_) return;
    
    const int BHadronIndex(bHadronIndex);
    const int AntiBHadronIndex(antiBHadronIndex);
    
    h_GenAll->Fill((*topGenObjects.GenTop_).M(), trueLevelWeight);
    h_GenAll_noweight->Fill((*topGenObjects.GenTop_).M(), trueLevelWeightNoPileup);
    
    //Begin: Select & Fill histograms with Leading pT and 2nd Leading pT: Lepton and BJet
    orderLV(LeadGenLepton, NLeadGenLepton, (*topGenObjects.GenLepton_), (*topGenObjects.GenAntiLepton_), LVpt);
    
    if (BHadronIndex != -1 && AntiBHadronIndex != -1) {
        orderLV(LeadGenBJet, NLeadGenBJet, (*topGenObjects.allGenJets_).at(BHadronIndex), (*topGenObjects.allGenJets_).at(AntiBHadronIndex), LVpt);
    }

    if ( (*topGenObjects.GenLepton_).pt() > LeptonPtCut && (*topGenObjects.GenAntiLepton_).pt() > LeptonPtCut &&
        std::fabs( (*topGenObjects.GenLepton_).eta() ) < LeptonEtaCUT && std::fabs ( (*topGenObjects.GenAntiLepton_).eta() ) < LeptonEtaCUT) {

        if ( BHadronIndex != -1 && AntiBHadronIndex != -1 ) {
            if ( (*topGenObjects.allGenJets_).at(BHadronIndex).pt() > JetPtCUT && std::fabs ( (*topGenObjects.allGenJets_).at(BHadronIndex).eta() ) < JetEtaCUT &&
                 (*topGenObjects.allGenJets_).at(AntiBHadronIndex).pt() > JetPtCUT && std::fabs ( (*topGenObjects.allGenJets_).at(AntiBHadronIndex).Eta() ) < JetEtaCUT )
            {

                h_VisGenAll->Fill((*topGenObjects.GenTop_).M(), trueLevelWeight);
                h_VisGenAll_noweight->Fill((*topGenObjects.GenTop_).M(), trueLevelWeightNoPileup);

                h_VisGenLLBarpT->Fill(( (*topGenObjects.GenLepton_) + (*topGenObjects.GenAntiLepton_) ).Pt(), trueLevelWeight );
                h_VisGenLLBarMass->Fill(( (*topGenObjects.GenLepton_) + (*topGenObjects.GenAntiLepton_) ).M(), trueLevelWeight );

                h_VisGenLeptonpT->Fill((*topGenObjects.GenLepton_).pt(), trueLevelWeight );
                h_VisGenAntiLeptonpT->Fill((*topGenObjects.GenAntiLepton_).Pt(), trueLevelWeight );

                h_VisGenLeptonEta->Fill((*topGenObjects.GenLepton_).Eta(), trueLevelWeight );
                h_VisGenAntiLeptonEta->Fill((*topGenObjects.GenAntiLepton_).Eta(), trueLevelWeight );

                h_VisGenBJetEta->Fill((*topGenObjects.allGenJets_).at(BHadronIndex).Eta(), trueLevelWeight );
                h_VisGenAntiBJetEta->Fill((*topGenObjects.allGenJets_).at(AntiBHadronIndex).Eta(), trueLevelWeight );
                h_VisGenBJetRapidity->Fill((*topGenObjects.allGenJets_).at(BHadronIndex).Rapidity(), trueLevelWeight );
                h_VisGenAntiBJetRapidity->Fill((*topGenObjects.allGenJets_).at(AntiBHadronIndex).Rapidity(), trueLevelWeight );
                h_VisGenBJetpT->Fill((*topGenObjects.allGenJets_).at(BHadronIndex).Pt(), trueLevelWeight );
                h_VisGenAntiBJetpT->Fill((*topGenObjects.allGenJets_).at(AntiBHadronIndex).Pt(), trueLevelWeight );
                h_VisGenMet->Fill((*topGenObjects.GenMet_).Pt(), trueLevelWeight);

                //for HT, count only >= JetPtCUT
                std::vector<int> genJetIndices = initialiseIndices(*topGenObjects.allGenJets_);
                selectIndices(genJetIndices, *topGenObjects.allGenJets_, LVpt, JetPtCUT);
                genHT = getJetHT(genJetIndices, *topGenObjects.allGenJets_);
                h_VisGenHT->Fill(genHT, trueLevelWeight);

                h_VisGenLLBarDPhi->Fill(std::fabs( DeltaPhi((*topGenObjects.GenLepton_), (*topGenObjects.GenAntiLepton_))), trueLevelWeight );
                h_VisGenLeptonantiBjetMass->Fill(( (*topGenObjects.GenLepton_) + (*topGenObjects.allGenJets_).at(AntiBHadronIndex) ).M(), trueLevelWeight );
                h_VisGenAntiLeptonBjetMass->Fill(( (*topGenObjects.GenAntiLepton_) + (*topGenObjects.allGenJets_).at(BHadronIndex) ).M(), trueLevelWeight );
                h_VisGenJetMult->Fill((*topGenObjects.allGenJets_).size(), trueLevelWeight );

                h_VisGenBBBarpT->Fill(((*topGenObjects.allGenJets_).at(BHadronIndex) + (*topGenObjects.allGenJets_).at(AntiBHadronIndex)).Pt(), trueLevelWeight );
                h_VisGenBBBarMass->Fill(((*topGenObjects.allGenJets_).at(BHadronIndex) + (*topGenObjects.allGenJets_).at(AntiBHadronIndex)).M(), trueLevelWeight );

                //Begin: Select & Fill histograms with Leading pT and 2nd Leading pT: Lepton and BJet
                h_VisGenLeptonpTLead->Fill(LeadGenLepton.Pt(), trueLevelWeight);
                h_VisGenLeptonpTNLead->Fill(NLeadGenLepton.Pt(), trueLevelWeight);
                h_VisGenLeptonEtaLead->Fill(LeadGenLepton.Eta(), trueLevelWeight);
                h_VisGenLeptonEtaNLead->Fill(NLeadGenLepton.Eta(), trueLevelWeight);

                h_VisGenBJetpTLead->Fill(LeadGenBJet.Pt(), trueLevelWeight);
                h_VisGenBJetpTNLead->Fill(NLeadGenBJet.Pt(), trueLevelWeight);
                h_VisGenBJetEtaLead->Fill(LeadGenBJet.Eta(), trueLevelWeight);
                h_VisGenBJetEtaNLead->Fill(NLeadGenBJet.Eta(), trueLevelWeight);
                //End: Select & Fill histograms with Leading pT and 2nd Leading pT: Lepton and BJet
            }
        }
    }
    
    LV genttbar((*topGenObjects.GenTop_) + (*topGenObjects.GenAntiTop_));
    h_VisGenTTBarMass->Fill(genttbar.M(), trueLevelWeight );
    h_VisGenTTBarRapidity->Fill(genttbar.Rapidity(), trueLevelWeight );
    h_VisGenTTBarpT->Fill(genttbar.Pt(), trueLevelWeight );
    
    h_VisGenTTBarRapidityvsTTBarpT->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Pt(),((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Rapidity(),trueLevelWeight);
    
    h_VisGenToppT->Fill((*topGenObjects.GenTop_).Pt(), trueLevelWeight );
    h_VisGenAntiToppT->Fill((*topGenObjects.GenAntiTop_).Pt(), trueLevelWeight );
    h_VisGenTopRapidity->Fill((*topGenObjects.GenTop_).Rapidity(), trueLevelWeight );
    h_VisGenAntiTopRapidity->Fill((*topGenObjects.GenAntiTop_).Rapidity(), trueLevelWeight );
    h_VisGenTopEta->Fill((*topGenObjects.GenTop_).Eta(), trueLevelWeight );
    h_VisGenAntiTopEta->Fill((*topGenObjects.GenAntiTop_).Eta(), trueLevelWeight );
    h_VisGenTTBarDeltaPhi->Fill(std::abs(DeltaPhi((*topGenObjects.GenTop_), (*topGenObjects.GenAntiTop_))), trueLevelWeight);
    h_VisGenTTBarDeltaRapidity->Fill(std::abs((*topGenObjects.GenTop_).Rapidity()) - std::abs((*topGenObjects.GenAntiTop_).Rapidity()), trueLevelWeight);

    h_VisGenNeutrinopT->Fill((*topGenObjects.GenNeutrino_).Pt(), trueLevelWeight);
    h_VisGenAntiNeutrinopT->Fill((*topGenObjects.GenAntiNeutrino_).Pt(), trueLevelWeight);

    /// Parton momentum fraction
    double partonMomFraction = ((*topGenObjects.GenTop_).energy() - (*topGenObjects.GenTop_).Pz() + (*topGenObjects.GenAntiTop_).energy() - (*topGenObjects.GenAntiTop_).Pz()) / (2 * 4000);
    double antipartonMomFraction = ((*topGenObjects.GenTop_).energy() + (*topGenObjects.GenTop_).Pz() + (*topGenObjects.GenAntiTop_).energy() + (*topGenObjects.GenAntiTop_).Pz()) / (2 * 4000);
    h_VisGenTopPartonFraction->Fill(partonMomFraction, trueLevelWeight);
    h_VisGenAntiTopPartonFraction->Fill(antipartonMomFraction, trueLevelWeight);

    //Begin: Fill histograms with Leading pT and 2nd Leading pT: Top
    orderLV(LeadGenTop, NLeadGenTop, (*topGenObjects.GenTop_), (*topGenObjects.GenAntiTop_), LVpt);
    h_VisGenToppTLead->Fill(LeadGenTop.Pt(), trueLevelWeight);
    h_VisGenToppTNLead->Fill(NLeadGenTop.Pt(), trueLevelWeight);
    h_VisGenTopRapidityLead->Fill(LeadGenTop.Rapidity(), trueLevelWeight);
    h_VisGenTopRapidityNLead->Fill(NLeadGenTop.Rapidity(), trueLevelWeight);
    h_VisGenTopMassLead->Fill(LeadGenTop.M(), trueLevelWeight);
    h_VisGenTopMassNLead->Fill(NLeadGenTop.M(), trueLevelWeight);
    //End: Fill histograms with Leading pT and 2nd Leading pT: Top
    
    // create top/antitop quark in the ttbar rest frame
    ROOT::Math::Boost CoMBoostGenTtbar (genttbar.BoostToCM());
    LV gentop ((*topGenObjects.GenTop_));
    LV genantitop ((*topGenObjects.GenAntiTop_));
    gentop = CoMBoostGenTtbar(gentop);
    genantitop = CoMBoostGenTtbar(genantitop);

    h_VisGenToppTTTRestFrame->Fill(gentop.Pt(), trueLevelWeight);
    h_VisGenAntiToppTTTRestFrame->Fill(genantitop.Pt(), trueLevelWeight);

}

void TopAnalysis::generatorTTbarjetsEvent(double& jetHTGen,
                                          const int bHadronIndex, const int antiBHadronIndex,
                                          const double trueLevelWeight,
                                          int& GenJets_cut, int& GenJets_cut40, int& GenJets_cut60, int& GenJets_cut100, int& jetnum,
                                          double extragenjet[4],
                                          const TopGenObjects& topGenObjects)
{
    // Use utilities without namespaces
    using namespace common;
    using ROOT::Math::VectorUtil::DeltaPhi;
    using ROOT::Math::VectorUtil::DeltaR;

//    double cbin[19]={35.,45.,55.,65.,75.,85.,95.,110.,130.,150.,170.,190.,210.,230.,250.,270.,300.,340.,380.};
    double cbin[20]={25.,35.,45.,55.,65.,75.,85.,95.,110.,130.,150.,170.,190.,210.,230.,250.,270.,300.,340.,380.};
    jetnum = -1;
    jetHTGen = 0.;
    
    if(!topGenObjects.valuesSet_) return;


    const int BHadronIndex(bHadronIndex);
    const int AntiBHadronIndex(antiBHadronIndex);
    
    if ( (*topGenObjects.GenLepton_).pt() > LeptonPtCut && (*topGenObjects.GenAntiLepton_).pt() > LeptonPtCut &&
        std::fabs( (*topGenObjects.GenLepton_).eta() ) < LeptonEtaCUT && std::fabs ( (*topGenObjects.GenAntiLepton_).eta() ) < LeptonEtaCUT ) {
        if ( BHadronIndex != -1 && AntiBHadronIndex != -1 ) {
            if ( (*topGenObjects.allGenJets_).at(BHadronIndex).Pt() > JetPtCUT && std::fabs ( (*topGenObjects.allGenJets_).at(BHadronIndex).Eta() ) < JetEtaCUT &&
                (*topGenObjects.allGenJets_).at(AntiBHadronIndex).Pt() > JetPtCUT && std::fabs ( (*topGenObjects.allGenJets_).at(AntiBHadronIndex).Eta() ) < JetEtaCUT && std::fabs(DeltaR((*topGenObjects.GenLepton_), (*topGenObjects.allGenJets_).at(BHadronIndex))) > 0.4  && std::fabs(DeltaR((*topGenObjects.GenLepton_), (*topGenObjects.allGenJets_).at(AntiBHadronIndex))) > 0.4  && std::fabs(DeltaR((*topGenObjects.GenAntiLepton_), (*topGenObjects.allGenJets_).at(BHadronIndex))) > 0.4  && std::fabs(DeltaR((*topGenObjects.GenAntiLepton_), (*topGenObjects.allGenJets_).at(AntiBHadronIndex))) > 0.4)
            {    
                GenJets_cut = GenJets_cut40 = GenJets_cut60 = GenJets_cut100 = 0;
                //New plots from Carmen: Begin
                for(int genJet=0; genJet<(int)(*topGenObjects.allGenJets_).size(); genJet++)
                {
                    if(std::fabs((*topGenObjects.allGenJets_).at(genJet).Eta() ) > JetEtaCUT || (*topGenObjects.allGenJets_).at(genJet).Pt() <= JetPtCUT2 || (std::fabs(DeltaR((*topGenObjects.GenLepton_), (*topGenObjects.allGenJets_).at(genJet))) < 0.4 && (*topGenObjects.allGenJets_).at(BHadronIndex) != (*topGenObjects.allGenJets_).at(genJet)  && (*topGenObjects.allGenJets_).at(AntiBHadronIndex) != (*topGenObjects.allGenJets_).at(genJet)) || (std::fabs(DeltaR((*topGenObjects.GenAntiLepton_), (*topGenObjects.allGenJets_).at(genJet))) < 0.4 && (*topGenObjects.allGenJets_).at(BHadronIndex) != (*topGenObjects.allGenJets_).at(genJet)  && (*topGenObjects.allGenJets_).at(AntiBHadronIndex) != (*topGenObjects.allGenJets_).at(genJet) ) )
                    {
                        continue;
                    }
//                    if((*topGenObjects.allGenJets_).at(genJet).Pt()> JetPtCUT2 && std::fabs((*topGenObjects.allGenJets_).at(genJet).Eta())<2.4) {//CARMEN ETA CUT
                        GenJets_cut++;
                        if((*topGenObjects.allGenJets_).at(BHadronIndex) != (*topGenObjects.allGenJets_).at(genJet) && (*topGenObjects.allGenJets_).at(AntiBHadronIndex) != (*topGenObjects.allGenJets_).at(genJet))
                        {
                            jetHTGen+=(*topGenObjects.allGenJets_).at(genJet).Pt();
                            if(jetnum < 3)
                            {
                                jetnum++;
                                extragenjet[jetnum] = genJet;
                            }
                        }
                        if((*topGenObjects.allGenJets_).at(genJet).Pt()> 40.0) GenJets_cut40++;
                        if((*topGenObjects.allGenJets_).at(genJet).Pt()> 60.0) GenJets_cut60++;
                        if((*topGenObjects.allGenJets_).at(genJet).Pt()> 100.0) GenJets_cut100++;
//                    }
                }
                h_VisGenJetMultpt30->Fill(GenJets_cut,trueLevelWeight);
                h_VisGenJetMultpt40->Fill(GenJets_cut40,trueLevelWeight);
                h_VisGenJetMultpt60->Fill(GenJets_cut60,trueLevelWeight);
                h_VisGenJetMultpt100->Fill(GenJets_cut100,trueLevelWeight);

                LV genttbar((*topGenObjects.GenTop_) + (*topGenObjects.GenAntiTop_));
                h_VisGenTTBar0Mass->Fill(rho0/genttbar.M(),trueLevelWeight);

                if(jetnum>2){
                    h_VisGenExtraJetpT4->Fill((*topGenObjects.allGenJets_).at(extragenjet[3]).Pt(),trueLevelWeight);
                    h_VisGenExtraJetEta4->Fill((*topGenObjects.allGenJets_).at(extragenjet[3]).Eta(),trueLevelWeight);
                } if(jetnum>1){
                    h_VisGenExtraJetpT3->Fill((*topGenObjects.allGenJets_).at(extragenjet[2]).Pt(),trueLevelWeight);
                    h_VisGenExtraJetEta3->Fill((*topGenObjects.allGenJets_).at(extragenjet[2]).Eta(),trueLevelWeight);
                } if(jetnum>0){
                    h_VisGenExtraJetpT2->Fill((*topGenObjects.allGenJets_).at(extragenjet[1]).Pt(),trueLevelWeight);
                    h_VisGenExtraJetEta2->Fill((*topGenObjects.allGenJets_).at(extragenjet[1]).Eta(),trueLevelWeight);
                    h_VisGenMassExtraJet12->Fill(((*topGenObjects.allGenJets_).at(extragenjet[0])+(*topGenObjects.allGenJets_).at(extragenjet[1])).M(),trueLevelWeight);
                    h_VisGenDeltaRExtraJet12->Fill(std::fabs(DeltaR((*topGenObjects.allGenJets_).at(extragenjet[0]),(*topGenObjects.allGenJets_).at(extragenjet[1]))),trueLevelWeight);
                } if(jetnum > -1 ){
                    h_VisGenExtraJetpT->Fill((*topGenObjects.allGenJets_).at(extragenjet[0]).Pt(),trueLevelWeight);
                    h_VisGenExtraJetEta->Fill((*topGenObjects.allGenJets_).at(extragenjet[0]).Eta(),trueLevelWeight);
                    LV genttbar((*topGenObjects.GenTop_) + (*topGenObjects.GenAntiTop_));
                    h_VisGenTTBar1stJetMass->Fill(rho0/(genttbar+(*topGenObjects.allGenJets_).at(extragenjet[0])).M(),trueLevelWeight);  
                }

                for(int q0 = 0; q0<20; q0++){
                    h_VisGenJetMultTotal->Fill(cbin[q0],trueLevelWeight);
                    if((*topGenObjects.allGenJets_).at(extragenjet[0]).Pt()<=cbin[q0] || jetnum<0 ) h_VisGenJetMultQ0->Fill(cbin[q0],trueLevelWeight);
                    if((*topGenObjects.allGenJets_).at(extragenjet[1]).Pt()<=cbin[q0] || jetnum<1 ) h_VisGenJetExtra2Q0->Fill(cbin[q0],trueLevelWeight);
                    if(jetHTGen<=cbin[q0]) h_VisGenJetMultQsum->Fill(cbin[q0],trueLevelWeight);
                }
            }
        }
    }
}


void TopAnalysis::generatorVisJets(const TopGenObjects& topGenObjects,std::vector<int>& genVisJetIndices)
{
    if(topGenObjects.valuesSet_)
    {
        for(int genJet=0; genJet<(int)(*topGenObjects.allGenJets_).size(); genJet++)
        {
            if(std::fabs((*topGenObjects.allGenJets_).at(genJet).Eta() ) < JetEtaCUT && (*topGenObjects.allGenJets_).at(genJet).Pt() <= JetPtCUT  )
            {
                genVisJetIndices.push_back(genJet);
            }
        }
    }
}


void TopAnalysis::CreateBinnedControlPlots(TH1* h_differential, TH1* h_control, const bool fromHistoList)
{
    auto &pair = (*binnedControlPlots_)[h_differential->GetName()];
    if(fromHistoList){
        HistoListReader histoList("HistoList");
        if(histoList.IsZombie()) { std::cout << "Need a HistoList to create binned control plots!\n"; exit(273); }
        pair.first = histoList.getPlotProperties(h_differential->GetName()).getClonedHistogram();
    }
    else{
        bool old = TH1::AddDirectoryStatus();
        TH1::AddDirectory(false);
        TH1* clone = static_cast<TH1*>(h_differential->Clone());
        TH1::AddDirectory(old);
        pair.first = clone;
    }
    std::string name = "bcp_";
    name.append(h_differential->GetName()).append("_bin_");
    //create maps if we are called for the first time with a certain h_differential
    if (pair.second.size() == 0) {
        for (int i = 0; i <= pair.first->GetNbinsX() + 1; ++i)
            pair.second.push_back(std::map<std::string, TH1*>());
    }
    //now really create the histograms
    for (int i = 0; i <= pair.first->GetNbinsX() + 1; ++i) {
        std::string binning = 
            i == 0 ? "underflow" :
            i == pair.first->GetNbinsX() + 1 ? "overflow" :
            common::d2s(pair.first->GetBinLowEdge(i)) + " to " + common::d2s(pair.first->GetBinLowEdge(i+1));
        binning = std::string(" (") + h_differential->GetName() + " " + binning + ")";
        std::string n = name + std::to_string(i) + "_" + h_control->GetName();
        pair.second[i][h_control->GetName()] = store(
            new TH1D(n.c_str(), (std::string(h_control->GetName())+binning).c_str(), 
                     h_control->GetNbinsX(), h_control->GetBinLowEdge(1), 
                     h_control->GetBinLowEdge(h_control->GetNbinsX()+1)));
    }
}



void TopAnalysis::FillBinnedControlPlot(TH1* h_differential, double binvalue, 
                                        TH1* h_control, double value, double weight)
{
    auto pair = (*binnedControlPlots_)[h_differential->GetName()];
    auto bin = pair.first->FindBin(binvalue);
    auto m = pair.second.at(bin);
    TH1* h = m[h_control->GetName()];
    if (!h) { 
        std::cerr << "Error: please CreateBinnedControlPlots for " << h_differential->GetName() 
            << " and " << h_control->GetName() << std::endl;
        exit(911);
    }
    h->Fill(value, weight);
}

void TopAnalysis::SetAllAnalyzers(std::vector<AnalyzerBaseClass*> v_analyzer)
{
    v_analyzer_ = v_analyzer;
}

void TopAnalysis::SetAllTreeHandlers(std::vector<TreeHandlerBase*> v_treeHandler)
{
    v_treeHandler_ = v_treeHandler;
}

void TopAnalysis::fillAll(const std::string& selectionStep,
                            const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                            const TopGenObjects& topGenObjects,
                            const KinRecoObjects& kinRecoObjects,
                            const ttbar::GenObjectIndices& genObjectIndices, const ttbar::RecoObjectIndices& recoObjectIndices,
                            const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                            const double& defaultWeight)const
{
    // In case b-tag efficiencies are produced, analysis output is not
    if(this->makeBtagEfficiencies()) return;
    
    for(AnalyzerBaseClass* analyzer : v_analyzer_){
        if(analyzer) analyzer->fill(recoObjects, commonGenObjects,
                                                   topGenObjects,
                                                   kinRecoObjects,
                                                   recoObjectIndices, genObjectIndices,
                                                   genLevelWeights, recoLevelWeights,
                                                   defaultWeight, selectionStep);
    }
    
    for(TreeHandlerBase* treeHandler : v_treeHandler_){
        if(treeHandler) treeHandler->fill(recoObjects, commonGenObjects,
                                                topGenObjects,
                                                kinRecoObjects,
                                                recoObjectIndices, genObjectIndices,
                                                genLevelWeights, recoLevelWeights,
                                                defaultWeight, selectionStep);
    }
    
}

void TopAnalysis::bookAll()
{
    for(AnalyzerBaseClass* analyzer : v_analyzer_){
        if(analyzer) analyzer->book(fOutput);
    }
}

void TopAnalysis::clearAll()
{
    for(AnalyzerBaseClass* analyzer : v_analyzer_){
        if(analyzer) analyzer->clear();
    }
}






