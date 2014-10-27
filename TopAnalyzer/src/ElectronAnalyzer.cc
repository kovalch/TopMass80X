#include "TMath.h"
#include "TopAnalysis/TopAnalyzer/interface/ElectronAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TopAnalysis/TopAnalyzer/interface/DileptonEventWeight.h"

ElectronAnalyzer::ElectronAnalyzer(const edm::ParameterSet& cfg):
  electrons_  ( cfg.getParameter<edm::InputTag>     ( "electrons"  )),
  jets_       ( cfg.getParameter<edm::InputTag>     ( "jets"       )),
  puWeight_   ( cfg.getParameter<edm::InputTag>     ( "weightPU"   )),
  lepSfWeight_( cfg.getParameter<edm::InputTag>     ( "weightLepSF")),    
  verbosity_  ( cfg.getParameter<bool>              ( "verbosity"  )),
  fromTo_     ( cfg.getParameter<std::vector<int> > ( "from_to"    ))
{

}

ElectronAnalyzer::~ElectronAnalyzer()
{
}

void
ElectronAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{
  double weight = getDileptonEventWeight(evt, puWeight_, lepSfWeight_);
  if(verbosity_){
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "      Run: " << evt.id().run()              << std::endl;
    std::cout << "LumiBlock: " << evt.id().luminosityBlock()  << std::endl;
    std::cout << "    Event: " << evt.id().event()            << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
  }

  edm::Handle<PatElectronCollection> electrons;
  evt.getByLabel(electrons_, electrons);

  edm::Handle<PatJetCollection> jets;
  evt.getByLabel(jets_, jets);

  multi_->Fill(electrons->size(), weight);

  int i=0;
  for(PatElectronCollection::const_iterator electron = electrons->begin(); electron!= electrons->end(); ++electron) {

    // only take the electrons from the selected range
    if(i<fromTo_[0]) continue;
    if(i>fromTo_[1]) break;

    ++i;
    pt_      ->Fill(electron->pt(), weight);
    eta_     ->Fill(electron->eta(), weight);
    phi_->Fill(electron->phi(), weight); /////

    supercluster_energy_->Fill(electron->superCluster()->energy() , weight);
    // hcaloverecal_->Fill(electron->hcalOverEcal() , weight);
    hcaloverecal_->Fill(electron->hadronicOverEm(), weight);
    fbrem_->Fill(electron->fbrem(), weight);
    // eoverp_->Fill(electron->eSuperClusterOverp() , weight);


    //if(electron->isEE()){
    if( fabs(electron->superCluster()->eta()) > 1.479 ){
      iso_combEE_->Fill(
                        ((electron->dr03TkSumPt()+electron->dr03EcalRecHitSumEt()+electron->dr03HcalTowerSumEt())/std::max<double>(20.,electron->pt())),
                        weight);
     iso_combEE2_->Fill(
                        ((electron->dr03TkSumPt()+electron->dr03EcalRecHitSumEt()+electron->dr03HcalTowerSumEt())/std::max<double>(20.,electron->pt())),
                        weight);
    }

    //if(electron->isEB()){
    if( fabs(electron->superCluster()->eta()) < 1.479 ){
      iso_combEB_->Fill(
                        ((electron->dr03TkSumPt()+std::max(0.,electron->dr03EcalRecHitSumEt()-1.)+electron->dr03HcalTowerSumEt())/std::max<double>(20.,electron->pt())),
                        weight);
      iso_combEB2_->Fill(
                        ((electron->dr03TkSumPt()+std::max(0.,electron->dr03EcalRecHitSumEt()-1.)+electron->dr03HcalTowerSumEt())/std::max<double>(20.,electron->pt())),
                        weight);
    }
    pfIso_->Fill((electron->chargedHadronIso()+electron->neutralHadronIso()+electron->photonIso()) / electron->pt());

    dB_->Fill(electron->dB(), weight);

    nlost_->Fill(electron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(), weight);

    if (jets.failedToGet()) {
        jet_dist_->Fill(jet_dist_->GetBinLowEdge(jet_dist_->GetNbinsX()), weight);
    } else {
        double minDr = 999.;
        for(PatJetCollection::const_iterator jet = jets->begin(); jet!= jets->end(); ++jet) {
            double dr = reco::deltaR<const pat::Electron, const pat::Jet>(*electron, *jet);
            if(dr<minDr) minDr = dr;
        }
        jet_dist_->Fill(minDr, weight);
    }

    d0_   ->Fill(electron->gsfTrack()->d0(), weight); /////

    deltaphi_->Fill(electron->deltaPhiSuperClusterTrackAtVtx(), weight);
    deltaeta_->Fill(electron->deltaEtaSuperClusterTrackAtVtx(), weight);
    sigmaietaieta_->Fill(electron->sigmaIetaIeta(), weight);

    double mapID= electron->electronID("simpleEleId90cIso");
    electronmap_->Fill(mapID, weight);

    convDcot_->Fill(electron->convDcot(), weight);
    convDist_->Fill(electron->convDist(), weight);

    if(verbosity_){
      std::cout << "-------------------------------------------" << std::endl;
      std::cout << "ELECTRON " << i-1 << ":" << std::endl;
      std::cout << "        pt = " << electron->pt() << std::endl;
      std::cout << "       eta = " << electron->eta() << std::endl;
      std::cout << "       phi = " << electron->phi() << std::endl;
      std::cout << "        dB = " << electron->dB() << std::endl;
    }}

  if(verbosity_){
    std::cout << "-------------------------------------------" << std::endl;
  }
}


void
ElectronAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;
  if( !fs ){
    throw edm::Exception( edm::errors::Configuration,
                          "TFile Service is not registered in cfg file" );
  }

  multi_= fs->make<TH1D>( "multi" , "N_{el} per Event", 10, -0.5, 9.5);
  multi_->GetXaxis()->SetTitle("N_{el}");
  multi_->GetYaxis()->SetTitle("N_{evts}");

  pt_= fs->make<TH1D>( "pt" , "Electron Transverse Momentum", 200, 0., 400.);
  pt_->GetXaxis()->SetTitle("p_{t,el} [GeV]");
  pt_->GetYaxis()->SetTitle("Electrons");

  eta_= fs->make<TH1D>( "eta", "Electron Pseudorapidity", 80,-4., 4.);
  eta_->GetXaxis()->SetTitle("#eta_{el}");
  eta_->GetYaxis()->SetTitle("Electrons");

  phi_ = fs->make<TH1D>("phi", "Electron Phi", 48, -TMath::Pi(), TMath::Pi());
  phi_->GetXaxis()->SetTitle("#phi_{e}");
  phi_->GetYaxis()->SetTitle("N per 7.5 degree");

  supercluster_energy_= fs->make<TH1D>("supercluster_energy", "Electron supercluster Energy", 300, 0., 300.);
  supercluster_energy_->GetXaxis()->SetTitle("E supercluster [GeV]");
  supercluster_energy_->GetYaxis()->SetTitle("Electrons");

  hcaloverecal_= fs->make<TH1D>("hcaloverecal", "HCAL over ECAL", 100, 0., 0.2);
  hcaloverecal_->GetXaxis()->SetTitle("HCAL over ECAL [GeV]");
  hcaloverecal_->GetYaxis()->SetTitle("Electrons");

  fbrem_= fs->make<TH1D>("fbrem", "fbrem", 100, 0., 1.);
  fbrem_->GetXaxis()->SetTitle("fbrem");
  fbrem_->GetYaxis()->SetTitle("Electrons");

  eoverp_= fs->make<TH1D>("eoverp", "Energy over P", 500, 0., 10.);
  eoverp_->GetXaxis()->SetTitle("E/P");
  eoverp_->GetYaxis()->SetTitle("Electrons");

  //iso_combEE_= fs->make<TH1D>("iso_combEE", "Electron Combined Isolation in Endcaps", 200, 0.,  1.);
  iso_combEE_= fs->make<TH1D>("iso_combEE", "Electron Combined Isolation in Endcaps", 300, 0.,  15.);
  iso_combEE_->GetXaxis()->SetTitle("I_{combEE}");
  iso_combEE_->GetYaxis()->SetTitle("Electrons");

  iso_combEE2_= fs->make<TH1D>("iso_combEE2", "Electron Combined Isolation in Endcaps", 200, 0.,  2.);
  iso_combEE2_->GetXaxis()->SetTitle("I_{combEE}");
  iso_combEE2_->GetYaxis()->SetTitle("Electrons");

  //iso_combEB_= fs->make<TH1D>("iso_combEB", "Electron Combined Isolation in Barrel", 200, 0.,  1.);
  iso_combEB_= fs->make<TH1D>("iso_combEB", "Electron Combined Isolation in Barrel", 300, 0.,  15.);
  iso_combEB_->GetXaxis()->SetTitle("I_{combEB}");
  iso_combEB_->GetYaxis()->SetTitle("Electrons");

  iso_combEB2_= fs->make<TH1D>("iso_combEB2", "Electron Combined Isolation in Barrel", 200, 0.,  2.);
  iso_combEB2_->GetXaxis()->SetTitle("I_{combEB}");
  iso_combEB2_->GetYaxis()->SetTitle("Electrons");

  pfIso_= fs->make<TH1D>("pf_iso", "Particle Flow Isolation", 200, 0.,  0.5);
  pfIso_->GetXaxis()->SetTitle("I_{PF}");
  pfIso_->GetYaxis()->SetTitle("N");

  d0_= fs->make<TH1D>( "d0", "D0 of Inner Track", 100, 0., 10.);
  d0_->GetXaxis()->SetTitle("d0 [cm]");
  d0_->GetYaxis()->SetTitle("N");

  //dB_= fs->make<TH1D>( "dB", "DB of Electron Track", 100, 0., 0.5);
  dB_= fs->make<TH1D>( "dB", "DB of Electron Track", 100, 0., 10.);
  dB_->GetXaxis()->SetTitle("dB [cm]");
  dB_->GetYaxis()->SetTitle("Electrons");

  jet_dist_= fs->make<TH1D>("jet_dist", "Distance to closest Jet", 100, 0.,  10.);
  jet_dist_->GetXaxis()->SetTitle("#Delta r");
  jet_dist_->GetYaxis()->SetTitle("N");

  deltaphi_= fs->make<TH1D>("deltaphi", "Deltaphi supercluster track", 400, -0.5,  0.5);
  deltaphi_->GetXaxis()->SetTitle("#Delta #phi");
  deltaphi_->GetYaxis()->SetTitle("Electrons");

  deltaeta_= fs->make<TH1D>("deltaeta", "Deltaeta supercluster track", 400, -0.02,  0.02);
  deltaeta_->GetXaxis()->SetTitle("#Delta #eta");
  deltaeta_->GetYaxis()->SetTitle("Electrons");

  sigmaietaieta_= fs->make<TH1D>("sigmaietaieta", "Sigmaietaieta", 200, 0.,  0.1);
  sigmaietaieta_->GetXaxis()->SetTitle("#sigmaI#etaI#eta");
  sigmaietaieta_->GetYaxis()->SetTitle("Electrons");

  electronmap_= fs->make<TH1D>("electronmap", "Electron Map", 10, 0.,  10);
  electronmap_->GetXaxis()->SetTitle("MapID");
  electronmap_->GetYaxis()->SetTitle("Electrons");

  nlost_= fs->make<TH1D>("nlost", "Conversions Nlost hits", 10, 0.,  10);
  nlost_->GetXaxis()->SetTitle("Nlost");
  nlost_->GetYaxis()->SetTitle("Electrons");

  convDcot_ = fs->make<TH1D>( "convDcot","convDcot", 200, 0., 1.);
  convDcot_->GetXaxis()->SetTitle("convDcot");
  convDcot_->GetYaxis()->SetTitle("Electrons");

  convDist_ = fs->make<TH1D>( "convDist","convDist", 200, 0., 1.);
  convDist_->GetXaxis()->SetTitle("convDist");
  convDist_->GetYaxis()->SetTitle("Electrons");

}


void
ElectronAnalyzer::endJob()
{
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( ElectronAnalyzer );
