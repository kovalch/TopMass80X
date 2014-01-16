#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TopAnalysis/TopAnalyzer/plugins/DiLeptonAnalyzer.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "TopAnalysis/TopAnalyzer/interface/DileptonEventWeight.h"

using namespace std;
using namespace edm;

DiLeptonAnalyzer::DiLeptonAnalyzer(const edm::ParameterSet& cfg):

  fileOutput_      (cfg.getParameter<bool>("fileOutput")),
  outputFile_      (cfg.getUntrackedParameter<string>("outputFile")),

  vertex_          (cfg.getParameter<edm::InputTag>("vertexCollection")),
  vertex_X_cut_    (cfg.getParameter<double>("vertex_X_cut")),
  vertex_Y_cut_    (cfg.getParameter<double>("vertex_Y_cut")),
  vertex_Z_cut_    (cfg.getParameter<double>("vertex_Z_cut")),

  loose_muons_     (cfg.getParameter<edm::InputTag>("loose_muons")),
  muons_           (cfg.getParameter<edm::InputTag>("muons")),
  elecs_           (cfg.getParameter<edm::InputTag>("elecs")),
  jets_            (cfg.getParameter<edm::InputTag>("jets")),
  puWeight_        (cfg.getParameter<edm::InputTag>("weightPU")),
  lepSfWeight_     (cfg.getParameter<edm::InputTag>("weightLepSF")),   

  muon_iso_cut_    (cfg.getParameter<double>("muon_iso_cut")),
  elec_iso_cut_    (cfg.getParameter<double>("elec_iso_cut")),

  MassWindow_up_   (cfg.getParameter<double>("MassWindow_up")),
  MassWindow_down_ (cfg.getParameter<double>("MassWindow_down"))
{
  if( fileOutput_ ) {
    const char *fileName = outputFile_.c_str();
    outfile.open (fileName);
  }
}


DiLeptonAnalyzer::~DiLeptonAnalyzer() {

}


void DiLeptonAnalyzer::beginJob() {

  edm::Service<TFileService> fs;

  if( !fs ){
    throw edm::Exception( edm::errors::Configuration,
                          "TFile Service is not registered in cfg file" );
  }

  Events_RC = fs->make<TH1F>("Events_RC", "Isolated dilepton events (RC)", 5, 0., 5.);
  Events_WC = fs->make<TH1F>("Events_WC", "Isolated dilepton events (WC)", 5, 0., 5.);

  TAxis* axis = Events_RC->GetXaxis();
  axis->SetBinLabel( 2, "#mu #mu");
  axis->SetBinLabel( 3, "#mu e");
  axis->SetBinLabel( 4, "e e");

  Events_iso = fs->make<TH2F>("Events_iso", "Isolated muons AND electrons", 5,  0., 5., 5, 0., 5.);

  for( int i = 0; i < 5; i++ ) {

    stringstream ss;

    ss << "N = " << i;
    const string &label = ss.str();

    Events_iso->GetXaxis()->SetBinLabel( i+1,label.c_str() );
    Events_iso->GetYaxis()->SetBinLabel( i+1,label.c_str() );

  }

  Ntracks       = fs->make<TH1F>("Ntracks",       "Number of tracks",              50, -0.5, 49.5);
  VxVy_muons    = fs->make<TH2F>("VxVy_muons",    "Vertex x-y-positon (global)",   40,  -1.,   1., 40 , -1., 1.);
  Vz_muons      = fs->make<TH1F>("Vz_muons",      "Vertex z-positon (global)",     40, -20.,  20.);

  Nmuons        = fs->make<TH1F>("Nmuons",        "Number of muons",               10, -0.5,  9.5);
  Nmuons_iso    = fs->make<TH1F>("Nmuons_iso",    "Number of isolated muons",      10, -0.5,  9.5);
  Nmuons_charge = fs->make<TH1F>("Nmuons_charge", "Number of muons * moun charge", 19, -10.,  10.);

  Nelecs        = fs->make<TH1F>("Nelecs",        "Number of electrons",           10, -0.5,  9.5);
  Nelecs_iso    = fs->make<TH1F>("Nelecs_iso",    "Number of isolated electrons",  10, -0.5,  9.5);
  Nelecs_charge = fs->make<TH1F>("Nelecs_charge", "Number of elecs * elec charge", 19, -10.,  10.);

  Muon_CombRelIso03 = fs->make<TH1F>("Muon_CombRelIso03", "Muon CombRelIso03",     50,  0.,  1.);
  Elec_CombRelIso03 = fs->make<TH1F>("Elec_CombRelIso03", "Electron CombRelIso03", 50,  0.,  1.);

  DiLeptonMass = fs->make<TH1F>("DiLeptonMass", "M^{l^{+}l^{-}} #left[#frac{GeV}{c^{2}}#right]", 250, 0., 500.);
  DiLeptonPt   = fs->make<TH1F>("DiLeptonPt",   "p_{T}^{l^{+}l^{-}} #left[#frac{GeV}{c}#right]", 250, 0., 500.);  
  LeptonMinPt = fs->make<TH1F>("LeptonMinPt", "min p_{T}^{l} #left[#frac{GeV}{c}#right]", 250, 0., 500.);
  LeptonMaxPt = fs->make<TH1F>("LeptonMaxPt", "max p_{T}^{l} #left[#frac{GeV}{c}#right]", 250, 0., 500.);    
  LeptonMinEta = fs->make<TH1F>("LeptonMinEta", "min #eta^{l}", 100, -5., 5.);
  LeptonMaxEta = fs->make<TH1F>("LeptonMaxEta", "max #eta^{l}", 100, -5., 5.);    

  const int nbins = 100;

  double logmin = 1.;
  double logmax = log(7000)/log(10);  // 100 bins from 1 GeV to 7 TeV

  float bins[nbins+1];

  for (int i = 0; i <= nbins; i++) {
    double log = logmin + (logmax-logmin)*i/nbins;
    bins[i] = std::pow(10.0, log);
  }

  DimassRC_MM = fs->make<TH1F>("DimassRC_MM", "Dilepton mass (mu mu) RC", 150, 0., 300.);
  DimassRC_EE = fs->make<TH1F>("DimassRC_EE", "Dilepton mass (e e)   RC", 150, 0., 300.);
  DimassRC_ME = fs->make<TH1F>("DimassRC_ME", "Dilepton mass (mu e)  RC", 150, 0., 300.);
  DimassWC_MM = fs->make<TH1F>("DimassWC_MM", "Dilepton mass (mu mu) WC", 150, 0., 300.);
  DimassWC_EE = fs->make<TH1F>("DimassWC_EE", "Dilepton mass (e e)   WC", 150, 0., 300.);
  DimassWC_ME = fs->make<TH1F>("DimassWC_ME", "Dilepton mass (mu e)  WC", 150, 0., 300.);

  DimassRC_LOGX_MM = fs->make<TH1F>("DimassRC_LOGX_MM", "Dilepton mass (mu mu) RC", nbins, &bins[0]);
  DimassRC_LOGX_EE = fs->make<TH1F>("DimassRC_LOGX_EE", "Dilepton mass (e e)   RC", nbins, &bins[0]);
  DimassRC_LOGX_ME = fs->make<TH1F>("DimassRC_LOGX_ME", "Dilepton mass (mu e)  RC", nbins, &bins[0]);
  DimassWC_LOGX_MM = fs->make<TH1F>("DimassWC_LOGX_MM", "Dilepton mass (mu mu) WC", nbins, &bins[0]);
  DimassWC_LOGX_EE = fs->make<TH1F>("DimassWC_LOGX_EE", "Dilepton mass (e e)   WC", nbins, &bins[0]);
  DimassWC_LOGX_ME = fs->make<TH1F>("DimassWC_LOGX_ME", "Dilepton mass (mu e)  WC", nbins, &bins[0]);

  D_R_muon_jet  = fs->make<TH1F>("D_R_muon_jet",  "(#Delta R)_{min} muon-jets",  100, -1., 9.);
  D_R_elec_muon = fs->make<TH1F>("D_R_elec_muon", "(#Delta R)_{min} elec-muons", 100, -1., 9.);
  D_R_elec_jet  = fs->make<TH1F>("D_R_elec_jet",  "(#Delta R)_{min} elec-jets",  100, -1., 9.);

  D_eta_muonsRC = fs->make<TH1F>("D_eta_muonsRC", "#Delta eta_muons RC", 50, -5., 5.);
  D_eta_muonsWC = fs->make<TH1F>("D_eta_muonsWC", "#Delta eta_muons WC", 50, -5., 5.);
  D_phi_muonsRC = fs->make<TH1F>("D_phi_muonsRC", "#Delta phi_muons RC", 50, -5., 5.);
  D_phi_muonsWC = fs->make<TH1F>("D_phi_muonsWC", "#Delta phi_muons WC", 50, -5., 5.);

  D_eta_elecsRC = fs->make<TH1F>("D_eta_elecsRC", "#Delta eta_elecs RC", 50, -5., 5.);
  D_eta_elecsWC = fs->make<TH1F>("D_eta_elecsWC", "#Delta eta_elecs WC", 50, -5., 5.);
  D_phi_elecsRC = fs->make<TH1F>("D_phi_elecsRC", "#Delta phi_elecs RC", 50, -5., 5.);
  D_phi_elecsWC = fs->make<TH1F>("D_phi_elecsWC", "#Delta phi_elecs WC", 50, -5., 5.);

  D_eta_leptsRC = fs->make<TH1F>("D_eta_leptsRC", "#Delta eta_lepts RC", 50, -5., 5.);
  D_eta_leptsWC = fs->make<TH1F>("D_eta_leptsWC", "#Delta eta_lepts WC", 50, -5., 5.);
  D_phi_leptsRC = fs->make<TH1F>("D_phi_leptsRC", "#Delta phi_lepts RC", 50, -5., 5.);
  D_phi_leptsWC = fs->make<TH1F>("D_phi_leptsWC", "#Delta phi_lepts WC", 50, -5., 5.);

  corr_pt = fs->make<TH2F>("corr_pt", "Lepton+pt vs Lepton-pt", 250, 0, 500, 250, 0, 500);
}


void DiLeptonAnalyzer::beginRun(const edm::Run& r, const edm::EventSetup& context) {

}

void DiLeptonAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& context) {

  double weight = getDileptonEventWeight(evt, puWeight_, lepSfWeight_);

  // ------------------------
  //  Global Event Variables
  // ------------------------

  int N_run   = (evt.id()).run();
  int N_event = (evt.id()).event();
  int N_lumi  =  evt.luminosityBlock();

  int N_iso_mu    = 0;
  int N_iso_el    = 0;
  int N_iso_lep   = 0;
  int N_jets_Pt30 = 0;

  reco::MuonCollection        Isolated_muons;
  reco::GsfElectronCollection Isolated_elecs;

  TLorentzVector diLepLVector;
  double DilepMass = 0.;
  double DilepPt   = 0.;

  double vertex_X  = 100.;
  double vertex_Y  = 100.;
  double vertex_Z  = 100.;

  double min_dR_muon_jets  = 999.;
  double min_dR_elec_jets  = 999.;
  double min_dR_elec_muons = 999.;

  edm::Handle< vector<pat::Jet> > jets;
  evt.getByLabel(jets_, jets);

  vector<pat::Jet>::const_iterator jet;

  if( !jets.failedToGet() ) {

    for(jet = jets->begin(); jet!= jets->end(); ++jet) {

      if( (jet->pt() > 30.) && (abs(jet->eta()) < 2.5) )  ++N_jets_Pt30;

    }

  }

  // ------------------------
  //  Analyze Primary Vertex
  // ------------------------

  edm::Handle<reco::VertexCollection> vertexs;
  evt.getByLabel(vertex_, vertexs);

  if( vertexs.failedToGet() ) {

    cout << endl << "----------------------------" << endl;
    cout << "--- NO PRIMARY VERTEX !! ---" << endl;
    cout << "----------------------------" << endl << endl;

  }

  if( !vertexs.failedToGet() ) {

    reco::Vertex primaryVertex = vertexs->front();

    int numberTracks = primaryVertex.tracksSize();
    //    double ndof      = primaryVertex.ndof();
    bool fake        = primaryVertex.isFake();

    Ntracks->Fill(numberTracks, weight);

    if( !fake && numberTracks > 3 ) {

      vertex_X = primaryVertex.x();
      vertex_Y = primaryVertex.y();
      vertex_Z = primaryVertex.z();

    }
  }

  // ------------------------
  //  Analyze Muon Isolation
  // ------------------------

  edm::Handle< vector<pat::Muon> > muons;
  evt.getByLabel(muons_, muons);

  vector<pat::Muon>::const_iterator muon;

  if( muons.failedToGet() ) {

    cout << endl << "-----------------------" << endl;
    cout << "--- NO PAT MUONS !! ---" << endl;
    cout << "-----------------------" << endl << endl;

  }

  if( !muons.failedToGet() ) {

    Nmuons->Fill( muons->size(), weight);

    for(muon = muons->begin(); muon!= muons->end(); ++muon) {

      float  N_muons = muons->size();
      float  Q_muon  = muon->charge();

      Nmuons_charge->Fill(N_muons*Q_muon, weight);

      double track_X = 100.;
      double track_Y = 100.;
      double track_Z = 100.;

      if( muon->isGlobalMuon() && muon->isTrackerMuon() ) {

	reco::TrackRef track = muon->globalTrack();

	track_X = track->vx();
	track_Y = track->vy();
	track_Z = track->vz();

	VxVy_muons->Fill(track_X, track_Y, weight);
	Vz_muons->Fill(track_Z, weight);

      }

      // Vertex and kinematic cuts

      if(       track_X > vertex_X_cut_ )  continue;
      if(       track_Y > vertex_Y_cut_ )  continue;
      if(       track_Z > vertex_Z_cut_ )  continue;

      // Find closest jet to muon

      if( !jets.failedToGet() ) {

	for(jet = jets->begin(); jet!= jets->end(); ++jet) {

	  double deltaR_muon_jets = deltaR( muon->eta(), muon->phi(), jet->eta(), jet->phi() );

	  if( deltaR_muon_jets < min_dR_muon_jets )  min_dR_muon_jets = deltaR_muon_jets;

	  //	  if( jet->pt() > 30. )  ++N_jets_Pt30;

	}

	D_R_muon_jet->Fill(min_dR_muon_jets, weight);

      }

      // Muon isolation

      reco::MuonIsolation muIso03 = muon->isolationR03();

      double muonCombRelIso = 999.;

      if ( muon->pt() != 0. )
	muonCombRelIso = ( muIso03.emEt + muIso03.hadEt + muIso03.hoEt + muIso03.sumPt ) / muon->pt();

      Muon_CombRelIso03->Fill( muonCombRelIso, weight );

      if( muonCombRelIso < muon_iso_cut_ ) {

	++N_iso_mu;
	Isolated_muons.push_back(*muon);

      }

    }

    Nmuons_iso->Fill(N_iso_mu, weight);

  }

  // ----------------------------
  //  Analyze Electron Isolation
  // ----------------------------

  edm::Handle< vector<pat::Electron> > elecs;
  evt.getByLabel(elecs_, elecs);

  vector<pat::Electron>::const_iterator elec;

  edm::Handle< vector<pat::Muon> > loose_muons;
  evt.getByLabel(loose_muons_, loose_muons);

  vector<pat::Muon>::const_iterator loose_muon;

  if( elecs.failedToGet() ) {

    cout << endl << "---------------------------" << endl;
    cout << "--- NO PAT ELECTRONS !! ---" << endl;
    cout << "---------------------------" << endl << endl;

  }

  if( !elecs.failedToGet() ) {

    Nelecs->Fill( elecs->size(), weight );

    for(elec = elecs->begin(); elec!= elecs->end(); ++elec) {

      float N_elecs = elecs->size();
      float Q_elec  = elec->charge();

      Nelecs_charge->Fill(N_elecs*Q_elec, weight);

      double track_X = 100.;
      double track_Y = 100.;
      double track_Z = 100.;

      reco::GsfTrackRef track = elec->gsfTrack();

      track_X = track->vx();
      track_Y = track->vy();
      track_Z = track->vz();

      // Vertex and kinematic cuts

      if( track_X > vertex_X_cut_ )  continue;
      if( track_Y > vertex_Y_cut_ )  continue;
      if( track_Z > vertex_Z_cut_ )  continue;

      // Find closest jet to elec

      if( !jets.failedToGet() ) {

	for(jet = jets->begin(); jet!= jets->end(); ++jet) {

	  double deltaR_elec_jets = deltaR( elec->eta(), elec->phi(), jet->eta(), jet->phi() );

	  if( deltaR_elec_jets < min_dR_elec_jets )  min_dR_elec_jets = deltaR_elec_jets;

	}

	D_R_elec_jet->Fill(min_dR_elec_jets, weight);

      }

      // Find closest loose muon to elec

      if( !loose_muons.failedToGet() ) {

	for(loose_muon = loose_muons->begin(); loose_muon!= loose_muons->end(); ++loose_muon) {

	  double deltaR_elec_muons = deltaR( elec->eta(), elec->phi(), loose_muon->eta(), loose_muon->phi() );

	  if( deltaR_elec_muons < min_dR_elec_muons )  min_dR_elec_muons = deltaR_elec_muons;

	}

	D_R_elec_muon->Fill(min_dR_elec_muons, weight);

      }

      // Electron isolation

      //reco::GsfElectron::IsolationVariables elecIso = elec->dr03IsolationVariables();

      double elecCombRelIso = 999.;

      if ( elec->et() > 20. && elec->superCluster()->eta() < 1.479 )
	elecCombRelIso = ( elec->dr03TkSumPt() + max(0.,elec->dr03EcalRecHitSumEt()-1.) + elec->dr03HcalTowerSumEt() ) / elec->et();

      if ( elec->et() > 20. && elec->superCluster()->eta() > 1.479 )
	elecCombRelIso = ( elec->dr03TkSumPt() + elec->dr03EcalRecHitSumEt() + elec->dr03HcalTowerSumEt() ) / elec->et();


      Elec_CombRelIso03->Fill( elecCombRelIso, weight );

      if( elecCombRelIso < elec_iso_cut_ ) {

	++N_iso_el;
	Isolated_elecs.push_back(*elec);

      }

    }

    Nelecs_iso->Fill(N_iso_el, weight);

  }

  Events_iso->Fill(N_iso_mu, N_iso_el, weight);

  N_iso_lep = N_iso_el + N_iso_mu;

  if (Isolated_elecs.size() + Isolated_muons.size() ==2) {
    const reco::RecoCandidate *lep1, *lep2;
    if (Isolated_elecs.size() == 0) {
        lep1 = &(Isolated_muons.at(0));
        lep2 = &(Isolated_muons.at(1));
    } else if (Isolated_elecs.size() == 1) {
        lep1 = &(Isolated_elecs.at(0));
        lep2 = &(Isolated_muons.at(0));
    } else {
        lep1 = &(Isolated_elecs.at(0));
        lep2 = &(Isolated_elecs.at(1));
    }
    if (lep1->charge() > 0) 
        corr_pt->Fill(lep1->pt(), lep2->pt(), weight);
    else
        corr_pt->Fill(lep2->pt(), lep1->pt(), weight);
  }
  
  
  // --------------------
  //  TWO Isolated MUONS
  // --------------------

  if( Isolated_muons.size() > 1 ) {

    // Vertex cut

    if( vertex_X < vertex_X_cut_ && vertex_Y < vertex_Y_cut_ && vertex_Z < vertex_Z_cut_ ) {

      for( int i = 0; i < (static_cast<int>(Isolated_muons.size()) - 1); ++i ) {

	for( int j = i+1; j < static_cast<int>(Isolated_muons.size()); ++j ) {

	  reco::MuonCollection::const_reference mu1 = Isolated_muons.at(i);
	  reco::MuonCollection::const_reference mu2 = Isolated_muons.at(j);

          LeptonMinPt->Fill(min(mu1.pt(),mu2.pt()));
	  LeptonMaxPt->Fill(max(mu1.pt(),mu2.pt()));
	  
          LeptonMinEta->Fill(fabs(mu1.eta())<=fabs(mu2.eta()) ? mu1.eta() : mu2.eta());
	  LeptonMaxEta->Fill(fabs(mu1.eta())> fabs(mu2.eta()) ? mu1.eta() : mu2.eta());	  
	  
          diLepLVector = TLorentzVector(mu1.px()+mu2.px(), mu1.py()+mu2.py(), 
			                mu1.pz()+mu2.pz(), mu1.energy()+mu2.energy());

	  DilepMass = diLepLVector.M();
	  DilepPt   = diLepLVector.Pt();
	  
	  DiLeptonMass->Fill(DilepMass, weight);
	  DiLeptonPt  ->Fill(DilepPt,   weight);
	  
	  if( DilepMass < 1. ) {

	    if( i > 0 ) {

	      Isolated_muons.erase(Isolated_muons.begin()+i);
	      --i;

	    }

	    continue;

	  }

	  // Opposite muon charges -> Right Charge (RC)

 	  if( mu1.charge()*mu2.charge() < 0. ) {

	    DimassRC_MM->Fill(      DilepMass, weight );
	    DimassRC_LOGX_MM->Fill( DilepMass, weight );

	    if( DilepMass > 12. ) {

	      Events_RC->Fill(1., weight);

	      if( fileOutput_ ) {

 		outfile <<  N_run                       << "\t";
		outfile <<  N_lumi                      << "\t";
		outfile <<  N_event                     << "\t";
 		outfile << "N_Lep_iso: " << N_iso_lep   << "\t";
 		outfile << "RC mu mu"                   << "\t";
		outfile << "DiLepMass: " << DilepMass   << "\t";
 		outfile << "N_jets: "    << N_jets_Pt30 << "\n";

// 		outfile << "--------------------"      << "\n";
// 		outfile << "      Run : " << N_run     << "\n";
// 		outfile << "    Event : " << N_event   << "\n";
// 		outfile << "LumiBlock : " << N_lumi    << "\n";
// 		outfile << "N_Lep_iso : " << N_iso_lep << "\n";
// 		outfile << "     Type : mu mu"         << "\n";
// 		outfile << "--------------------"      << "\n";
// 		outfile << "DilepMass : "       << DilepMass << "\n";
// 		outfile << "Mu"<<i<<" Pt    : " << mu1.pt()  << "\n";
// 		outfile << "Mu"<<i<<" Eta   : " << mu1.eta() << "\n";
// 		outfile << "Mu"<<i<<" Phi   : " << mu1.phi() << "\n";
// 		outfile << "Mu"<<j<<" Pt    : " << mu2.pt()  << "\n";
// 		outfile << "Mu"<<j<<" Eta   : " << mu2.eta() << "\n";
// 		outfile << "Mu"<<j<<" Phi   : " << mu2.phi() << "\n";
// 		outfile << "--------------------"            << "\n";

	      }

	    }

	    if( DilepMass > MassWindow_down_ && DilepMass < MassWindow_up_ ) {

	      D_eta_muonsRC->Fill(mu1.eta()-mu2.eta(), weight);
	      D_phi_muonsRC->Fill(deltaPhi(mu1.phi(),mu2.phi()), weight);

	    }

	  }

	  // Same muon charges -> Wrong Charge (WC)

	  if( mu1.charge()*mu2.charge() > 0. ) {

	    DimassWC_MM->Fill(      DilepMass, weight );
	    DimassWC_LOGX_MM->Fill( DilepMass, weight );

	    if( DilepMass > 12. ) {

	      Events_WC->Fill(1., weight);

	      if( fileOutput_ ) {

 		outfile <<  N_run                       << "\t";
		outfile <<  N_lumi                      << "\t";
		outfile <<  N_event                     << "\t";
 		outfile << "N_Lep_iso: " << N_iso_lep   << "\t";
 		outfile << "WC mu mu"                   << "\t";
		outfile << "DiLepMass: " << DilepMass   << "\t";
 		outfile << "N_jets: "    << N_jets_Pt30 << "\n";

	      }

	    }

	    if( DilepMass > MassWindow_down_ && DilepMass < MassWindow_up_ ) {

	      D_eta_muonsWC->Fill(mu1.eta()-mu2.eta(), weight);
	      D_phi_muonsWC->Fill(deltaPhi(mu1.phi(),mu2.phi()), weight);

	    }

	  }

	}

      }

    }

  }

  // -----------------------------
  //  TWO Isolated LEPTONS (mu/e)
  // -----------------------------

  if( Isolated_muons.size() > 0 && Isolated_elecs.size() > 0 ) {

    // Vertex cut

    if( vertex_X < vertex_X_cut_ && vertex_Y < vertex_Y_cut_ && vertex_Z < vertex_Z_cut_ ) {

      for( int i = 0; i < static_cast<int>(Isolated_muons.size()); ++i ) {

	for( int j = 0; j < static_cast<int>(Isolated_elecs.size()); ++j ) {

	  reco::MuonCollection::const_reference        mu1 = Isolated_muons.at(i);
	  reco::GsfElectronCollection::const_reference el1 = Isolated_elecs.at(j);

          LeptonMinPt->Fill(min(mu1.pt(),el1.pt()));
	  LeptonMaxPt->Fill(max(mu1.pt(),el1.pt()));
	  
          LeptonMinEta->Fill(fabs(mu1.eta())<=fabs(el1.eta()) ? mu1.eta() : el1.eta());
	  LeptonMaxEta->Fill(fabs(mu1.eta())> fabs(el1.eta()) ? mu1.eta() : el1.eta());	

          diLepLVector = TLorentzVector(mu1.px()+el1.px(), mu1.py()+el1.py(), 
			                mu1.pz()+el1.pz(), mu1.energy()+el1.energy());

	  DilepMass = diLepLVector.M();
	  DilepPt   = diLepLVector.Pt();
	  
	  DiLeptonMass->Fill(DilepMass, weight);
	  DiLeptonPt  ->Fill(DilepPt,   weight);

	  if( DilepMass < 1. )  continue;

	  // Opposite lepton charges -> Right Charge (RC)

	  if( mu1.charge()*el1.charge() < 0. ) {

	    DimassRC_ME->Fill(      DilepMass, weight );
	    DimassRC_LOGX_ME->Fill( DilepMass, weight );

	    if( DilepMass > 12. ) {

	      Events_RC->Fill(2., weight);

	      if( fileOutput_ ) {

//  		cout <<  N_run                       << "\t";
// 		cout <<  N_lumi                      << "\t";
// 		cout <<  N_event                     << "\t";
//  		cout << "N_Lep_iso: " << N_iso_lep   << "\t";
//  		cout << "RC mu el"                   << "\t";
// 		cout << "DiLepMass: " << DilepMass   << "\t";
//  		cout << "N_jets: "    << N_jets_Pt30 << endl;

 		outfile <<  N_run                       << "\t";
		outfile <<  N_lumi                      << "\t";
		outfile <<  N_event                     << "\t";
 		outfile << "N_Lep_iso: " << N_iso_lep   << "\t";
 		outfile << "RC mu el"                   << "\t";
		outfile << "DiLepMass: " << DilepMass   << "\t";
 		outfile << "N_jets: "    << N_jets_Pt30 << "\n";

// 		outfile << "--------------------"      << "\n";
// 		outfile << "      Run : " << N_run     << "\n";
// 		outfile << "    Event : " << N_event   << "\n";
// 		outfile << "LumiBlock : " << N_lumi    << "\n";
// 		outfile << "N_Lep_iso : " << N_iso_lep << "\n";
// 		outfile << "     Type : mu el"         << "\n";
// 		outfile << "--------------------"      << "\n";
// 		outfile << "DilepMass : "       << DilepMass << "\n";
// 		outfile << "Mu"<<i<<" Pt    : " << mu1.pt()  << "\n";
// 		outfile << "Mu"<<i<<" Eta   : " << mu1.eta() << "\n";
// 		outfile << "Mu"<<i<<" Phi   : " << mu1.phi() << "\n";
// 		outfile << "El"<<j<<" Pt    : " << el1.pt()  << "\n";
// 		outfile << "El"<<j<<" Eta   : " << el1.eta() << "\n";
// 		outfile << "El"<<j<<" Phi   : " << el1.phi() << "\n";
// 		outfile << "--------------------"            << "\n";

	      }

	    }

	    if( DilepMass > MassWindow_down_ && DilepMass < MassWindow_up_ ) {

	      D_eta_leptsRC->Fill(mu1.eta()-el1.eta(), weight);
	      D_phi_leptsRC->Fill(deltaPhi(mu1.phi(),el1.phi()), weight);

	    }

	  }

	  // Same muon charges -> Wrong Charge (WC)

	  if( mu1.charge()*el1.charge() > 0. ) {

	    DimassWC_ME->Fill(      DilepMass, weight );
	    DimassWC_LOGX_ME->Fill( DilepMass, weight );

	    if( DilepMass > 12. ) {

	      Events_WC->Fill(2., weight);

	      if( fileOutput_ ) {

 		outfile <<  N_run                       << "\t";
		outfile <<  N_lumi                      << "\t";
		outfile <<  N_event                     << "\t";
 		outfile << "N_Lep_iso: " << N_iso_lep   << "\t";
 		outfile << "WC mu el"                   << "\t";
		outfile << "DiLepMass: " << DilepMass   << "\t";
 		outfile << "N_jets: "    << N_jets_Pt30 << "\n";

	      }

	    }

	    if( DilepMass > MassWindow_down_ && DilepMass < MassWindow_up_ ) {

	      D_eta_leptsWC->Fill(mu1.eta()-el1.eta(), weight);
	      D_phi_leptsWC->Fill(deltaPhi(mu1.phi(),el1.phi()), weight);

	    }

	  }

	}

      }

    }

  }

  // ------------------------
  //  TWO Isolated ELECTRONS
  // ------------------------

  if( Isolated_elecs.size() > 1 ) {

    // Vertex cut

    if( vertex_X < vertex_X_cut_ && vertex_Y < vertex_Y_cut_ && vertex_Z < vertex_Z_cut_ ) {

      for( int i = 0; i < (static_cast<int>(Isolated_elecs.size()) - 1); ++i ) {

	for( int j = i+1; j < static_cast<int>(Isolated_elecs.size()); ++j ) {

	  reco::GsfElectronCollection::const_reference el1 = Isolated_elecs.at(i);
	  reco::GsfElectronCollection::const_reference el2 = Isolated_elecs.at(j);

          LeptonMinPt->Fill(min(el1.pt(),el2.pt()));
	  LeptonMaxPt->Fill(max(el1.pt(),el2.pt()));
	  
          LeptonMinEta->Fill(fabs(el1.eta())<=fabs(el2.eta()) ? el1.eta() : el2.eta());
	  LeptonMaxEta->Fill(fabs(el1.eta())> fabs(el2.eta()) ? el1.eta() : el2.eta());	

          diLepLVector = TLorentzVector(el1.px()+el2.px(), el1.py()+el2.py(), 
			                el1.pz()+el2.pz(), el1.energy()+el2.energy());

	  DilepMass = diLepLVector.M();
	  DilepPt   = diLepLVector.Pt();
	  
	  DiLeptonMass->Fill(DilepMass, weight);
	  DiLeptonPt  ->Fill(DilepPt,   weight);

	  if( DilepMass < 1. )  continue;

	  // Opposite lepton charges -> Right Charge (RC)

	  if( el1.charge()*el2.charge() < 0. ) {

	    DimassRC_EE->Fill(      DilepMass, weight );
	    DimassRC_LOGX_EE->Fill( DilepMass, weight );

	    if( DilepMass > 12. ) {

	      Events_RC->Fill(3., weight);

	      if( fileOutput_ ) {

 		outfile <<  N_run                       << "\t";
		outfile <<  N_lumi                      << "\t";
		outfile <<  N_event                     << "\t";
 		outfile << "N_Lep_iso: " << N_iso_lep   << "\t";
 		outfile << "RC el el"                   << "\t";
		outfile << "DiLepMass: " << DilepMass   << "\t";
 		outfile << "N_jets: "    << N_jets_Pt30 << "\n";

// 		outfile << "--------------------"      << "\n";
// 		outfile << "      Run : " << N_run     << "\n";
// 		outfile << "    Event : " << N_event   << "\n";
// 		outfile << "LumiBlock : " << N_lumi    << "\n";
// 		outfile << "N_Lep_iso : " << N_iso_lep << "\n";
// 		outfile << "     Type : el el"         << "\n";
// 		outfile << "--------------------"      << "\n";
// 		outfile << "DilepMass : "       << DilepMass << "\n";
// 		outfile << "El"<<i<<" Pt    : " << el1.pt()  << "\n";
// 		outfile << "El"<<i<<" Eta   : " << el1.eta() << "\n";
// 		outfile << "El"<<i<<" Phi   : " << el1.phi() << "\n";
// 		outfile << "El"<<j<<" Pt    : " << el2.pt()  << "\n";
// 		outfile << "El"<<j<<" Eta   : " << el2.eta() << "\n";
// 		outfile << "El"<<j<<" Phi   : " << el2.phi() << "\n";
// 		outfile << "--------------------"            << "\n";

	      }

	    }

	    if( DilepMass > MassWindow_down_ && DilepMass < MassWindow_up_ ) {

	      D_eta_elecsRC->Fill(el1.eta()-el2.eta(), weight);
	      D_phi_elecsRC->Fill(deltaPhi(el1.phi(),el2.phi()), weight);

	    }

	  }

	  // Same muon charges -> Wrong Charge (WC)

	  if( el1.charge()*el2.charge() > 0. ) {

	    DimassWC_EE->Fill(      DilepMass, weight );
	    DimassWC_LOGX_EE->Fill( DilepMass, weight );

	    if( DilepMass > 12. ) {

	      Events_WC->Fill(3., weight);

	      if( fileOutput_ ) {

 		outfile <<  N_run                       << "\t";
		outfile <<  N_lumi                      << "\t";
		outfile <<  N_event                     << "\t";
 		outfile << "N_Lep_iso: " << N_iso_lep   << "\t";
 		outfile << "WC el el"                   << "\t";
		outfile << "DiLepMass: " << DilepMass   << "\t";
 		outfile << "N_jets: "    << N_jets_Pt30 << "\n";

	      }

	    }

	    if( DilepMass > MassWindow_down_ && DilepMass < MassWindow_up_ ) {

	      D_eta_elecsWC->Fill(el1.eta()-el2.eta(), weight);
	      D_phi_elecsWC->Fill(deltaPhi(el1.phi(),el2.phi()), weight);

	    }

	  }

	}

      }

    }

  }

}


void DiLeptonAnalyzer::endRun(const edm::Run& r, const edm::EventSetup& context) {

}

void DiLeptonAnalyzer::endJob() {

  if(fileOutput_)  outfile.close();

}

DEFINE_FWK_MODULE( DiLeptonAnalyzer );
