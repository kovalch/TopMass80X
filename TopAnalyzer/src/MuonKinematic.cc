#include "TopAnalysis/TopAnalyzer/interface/MuonKinematic.h"


/// constructor for FWLite analyzer
MuonKinematic::MuonKinematic(double dRMax):
  jets_(), towers_(), tracks_(), dRMax_(dRMax)
{ 
}

/// constructor for full FW analyzer
MuonKinematic::MuonKinematic(const edm::ParameterSet& cfg):
  jets_  ( cfg.getParameter<edm::InputTag>( "jets" ) ),
  towers_( cfg.getParameter<edm::InputTag>( "towers" ) ),
  tracks_( cfg.getParameter<edm::InputTag>( "tracks" ) ),  
  dRMax_ ( cfg.getParameter<double>( "dRMax" ) )  
{
}

/// fill interface for full FW analyzer
void
MuonKinematic::fill(const edm::Event& evt, const std::vector<pat::Muon>& muons, const double& weight=1.)
{
  if(muons.begin()==muons.end()) return;
  
  // jet isolation
  edm::Handle<std::vector<pat::Jet> > jets; 
  evt.getByLabel(jets_, jets);

  // track isolation
  edm::Handle<reco::TrackCollection> tracks;
  evt.getByLabel(tracks_, tracks);
  
  // calo isolation
  edm::Handle<reco::CandidateCollection> towers;
  evt.getByLabel(towers_, towers);
  
  fill(*jets, *tracks, *towers, muons, weight);
}

/// fill interface for FWLite analyzer
void
MuonKinematic::fill(const std::vector<pat::Jet>& jets, 
		    const reco::TrackCollection& tracks, 
		    const reco::CandidateCollection& towers,
		    const std::vector<pat::Muon>& muons, const double& weight=1.)
{
  std::vector<pat::Muon>::const_iterator muon=muons.begin();
  if(muon!=muons.end()){
    //---------------------------------------------
    // basic kinematics
    //---------------------------------------------
    en_ ->Fill( muon->energy(), weight );
    pt_ ->Fill( muon->et(),     weight );    
    eta_->Fill( muon->eta(),    weight );
    phi_->Fill( muon->phi(),    weight );

    //---------------------------------------------
    // jet isolation
    //---------------------------------------------
    double minDR   =-1.;
    double minDR_5 =-1.;
    double minDR_10=-1.;
    double minDR_15=-1.;
    double minDR_20=-1.;
    for(std::vector<pat::Jet>::const_iterator jet = jets.begin(); 
	jet!=jets.end(); ++jet) {
      double dR=deltaR(muon->eta(), muon->phi(), jet->eta(), jet->phi());
      if(  minDR   <0 || dR<minDR   )                   minDR   =dR;
      if( (minDR_5 <0 || dR<minDR_5 ) && jet->pt()>5  ) minDR_5 =dR;
      if( (minDR_10<0 || dR<minDR_10) && jet->pt()>10 ) minDR_10=dR;
      if( (minDR_15<0 || dR<minDR_15) && jet->pt()>15 ) minDR_15=dR;
      if( (minDR_20<0 || dR<minDR_20) && jet->pt()>20 ) minDR_20=dR;
    }
    if( minDR   >=0 ) isoJet_  ->Fill( minDR   , weight );
    if( minDR_5 >=0 ) isoJet5_ ->Fill( minDR_5 , weight );
    if( minDR_10>=0 ) isoJet10_->Fill( minDR_10, weight );
    if( minDR_15>=0 ) isoJet15_->Fill( minDR_15, weight );
    if( minDR_20>=0 ) isoJet20_->Fill( minDR_20, weight );

    //---------------------------------------------
    // calo isolation
    //---------------------------------------------
    
    // count number of towers in cone
    int nTowers = 0;
    for(reco::CandidateCollection::const_iterator tower = towers.begin(); 
	tower != towers.end(); ++tower) {
      double dR=deltaR( muon->eta(), muon->phi(), tower->eta(), tower->phi() );
      dRCalN_ ->Fill(dR);
      dRCalPt_->Fill(dR, tower->et());
      if(dR<dRMax_)++nTowers;
    }
    isoCalN_->Fill(nTowers, weight);
    
    //---------------------------------------------
    // track isolation
    //---------------------------------------------

    // search for track closest to the muon
    double dRmin = -1.;
    for(reco::TrackCollection::const_iterator track = tracks.begin(); 
	track != tracks.end(); ++track) {
      double dR=deltaR( muon->eta(), muon->phi(), track->eta(), track->phi() );
      if(dRmin<0 || dR<dRmin) dRmin=dR;
    }
    closestCtf_->Fill(TMath::Log10(dRmin), weight);
    
    // count number of tracks in cone
    int nTracks = 0;
    for(reco::TrackCollection::const_iterator track = tracks.begin(); 
	track != tracks.end(); ++track) {
      double dR=deltaR( muon->eta(), muon->phi(), track->eta(), track->phi() );
      if( dR==dRmin) 
	continue;
      dRTrkN_->Fill(dR);
      if(track->pt()<100.)
	dRTrkPt_->Fill(dR, track->pt());    
      if(dR<dRMax_) ++nTracks;
    }
    isoTrkN_->Fill(nTracks, weight);

    //test
    double offset=-0.1;
    for(int i=0; i<42; ++i){
      offset+=i*(2.1/42);
      std::cout << "isoDeposit: " << offset << " : " << muon->trackerIsoDeposit()-> depositWithin(offset);
    }
    
    //---------------------------------------------
    // fill std isolation plots
    //---------------------------------------------
    isoTrk_ ->Fill( muon->trackIso(), weight );
    isoCal_ ->Fill( muon->caloIso (), weight );
    isoEcal_->Fill( muon->ecalIso (), weight );
    isoHcal_->Fill( muon->hcalIso (), weight );
    
    ptVsTrkIso_ ->Fill( muon->pt(), muon->trackIso() );
    ptVsCalIso_ ->Fill( muon->pt(), muon->caloIso () );
    ptVsEcalIso_->Fill( muon->pt(), muon->ecalIso () );
    ptVsHcalIso_->Fill( muon->pt(), muon->hcalIso () );
  }
}

/// book for FWLite
void 
MuonKinematic::book()
{
  NameScheme kin("kin");
  en_ = new TH1F(kin.name( "en"  ), kin.name("en"  ), 60,   0., 300.);
  pt_ = new TH1F(kin.name( "pt"  ), kin.name("pt"  ), 30,   0., 150.);
  eta_= new TH1F(kin.name( "eta" ), kin.name("eta" ), 35, -3.5,  3.5);
  phi_= new TH1F(kin.name( "phi" ), kin.name("phi" ), 35, -3.5,  3.5);

  NameScheme iso("iso");
  isoJet_    = new TH1F(iso.name( "isoJet"  ), iso.name("isoJet"  ), 80,   0.,  4.);
  isoJet5_   = new TH1F(iso.name( "isoJet5" ), iso.name("isoJet5" ), 80,   0.,  4.);
  isoJet10_  = new TH1F(iso.name( "isoJet10"), iso.name("isoJet10"), 80,   0.,  4.);
  isoJet15_  = new TH1F(iso.name( "isoJet15"), iso.name("isoJet15"), 80,   0.,  4.);
  isoJet20_  = new TH1F(iso.name( "isoJet20"), iso.name("isoJet20"), 80,   0.,  4.);
  isoTrk_    = new TH1F(iso.name( "isoTrk"  ), iso.name("isoTrk"  ), 60,  -1.,  5.);
  isoCal_    = new TH1F(iso.name( "isoCal"  ), iso.name("isoCal"  ), 40, -10., 30.);
  isoEcal_   = new TH1F(iso.name( "isoEcal" ), iso.name("isoEcal" ), 40, -10., 30.);
  isoHcal_   = new TH1F(iso.name( "isoHcal" ), iso.name("isoHcal" ), 40, -10., 30.);
  dRTrkPt_   = new TH1F(iso.name( "dRTrkPt" ), iso.name("dRTrkPt" ), 42, -0.1,  2.);
  dRTrkN_    = new TH1F(iso.name( "dRTrkN"  ), iso.name("dRTrkN"  ), 42, -0.1,  2.);
  dRCalPt_   = new TH1F(iso.name( "dRCalPt" ), iso.name("dRCalPt" ), 42, -0.1,  2.);
  dRCalN_    = new TH1F(iso.name( "dRCalN"  ), iso.name("dRCalN"  ), 42, -0.1,  2.);
  isoTrkN_   = new TH1F(iso.name( "isoTrkN" ), iso.name("isoTrkN" ), 21,  -1., 20.);
  isoCalN_   = new TH1F(iso.name( "isoCalN" ), iso.name("isoCaloN"), 31,  -1., 30.);
  closestCtf_= new TH1F(iso.name( "closestCtf" ), iso.name("closestCtf" ), 25, -6., -1.);

  ptVsTrkIso_ = new TH2F(iso.name( "ptVsTrkIso" ), iso.name( "ptVsTrkIso" ), 100, 0., 100., 50,   0., 25.);
  ptVsCalIso_ = new TH2F(iso.name( "ptVsCalIso" ), iso.name( "ptVsCalIso" ), 100, 0., 100., 50, -10., 25.);
  ptVsEcalIso_= new TH2F(iso.name( "ptVsEcalIso"), iso.name( "ptVsEcalIso"), 100, 0., 100., 50, -10., 25.);
  ptVsHcalIso_= new TH2F(iso.name( "ptVsHcalIso"), iso.name( "ptVsHcalIso"), 100, 0., 100., 50, -10., 25.);
}

/// book for full FW
void 
MuonKinematic::book(edm::Service<TFileService>& fs)
{
  NameScheme kin("kin");
  en_ = fs->make<TH1F>(kin.name( "en"  ), kin.name("en"  ), 60,   0., 300.);
  pt_ = fs->make<TH1F>(kin.name( "pt"  ), kin.name("pt"  ), 30,   0., 150.);
  eta_= fs->make<TH1F>(kin.name( "eta" ), kin.name("eta" ), 35, -3.5,  3.5);
  phi_= fs->make<TH1F>(kin.name( "phi" ), kin.name("phi" ), 35, -3.5,  3.5);

  NameScheme iso("iso");
  isoJet_  = fs->make<TH1F>(iso.name( "isoJet"  ), iso.name("isoJet"  ), 80,   0.,  4.);
  isoJet5_ = fs->make<TH1F>(iso.name( "isoJet5" ), iso.name("isoJet5" ), 80,   0.,  4.);
  isoJet10_= fs->make<TH1F>(iso.name( "isoJet10"), iso.name("isoJet10"), 80,   0.,  4.);
  isoJet15_= fs->make<TH1F>(iso.name( "isoJet15"), iso.name("isoJet15"), 80,   0.,  4.);
  isoJet20_= fs->make<TH1F>(iso.name( "isoJet20"), iso.name("isoJet20"), 80,   0.,  4.);
  isoTrk_  = fs->make<TH1F>(iso.name( "isoTrk"  ), iso.name("isoTrk"  ), 60,  -1.,  5.);
  isoCal_  = fs->make<TH1F>(iso.name( "isoCal"  ), iso.name("isoCal"  ), 40, -10., 30.);
  isoEcal_ = fs->make<TH1F>(iso.name( "isoEcal" ), iso.name("isoEcal" ), 40, -10., 30.);
  isoHcal_ = fs->make<TH1F>(iso.name( "isoHcal" ), iso.name("isoHcal" ), 40, -10., 30.);
  dRTrkPt_ = fs->make<TH1F>(iso.name( "dRTrkPt" ), iso.name("dRTrkPt" ), 42, -0.1,  2.);
  dRTrkN_  = fs->make<TH1F>(iso.name( "dRTrkN"  ), iso.name("dRTrkN"  ), 42, -0.1,  2.);
  dRCalPt_ = fs->make<TH1F>(iso.name( "dRCalPt" ), iso.name("dRCalPt" ), 42, -0.1,  2.);
  dRCalN_  = fs->make<TH1F>(iso.name( "dRCalN"  ), iso.name("dRCalN"  ), 42, -0.1,  2.);
  isoTrkN_ = fs->make<TH1F>(iso.name( "isoTrkN" ), iso.name("isoTrkN" ), 21,  -1., 20.);
  isoCalN_ = fs->make<TH1F>(iso.name( "isoCalN" ), iso.name("isoCaloN"), 31,  -1., 30.);
  closestCtf_= fs->make<TH1F>(iso.name( "closestCtf" ), iso.name("closestCtf" ), 25, -6., -1.);

  ptVsTrkIso_ = fs->make<TH2F>(iso.name( "ptVsTrkIso" ), iso.name( "ptVsTrkIso" ), 100, 0., 100., 50,   0., 25.);
  ptVsCalIso_ = fs->make<TH2F>(iso.name( "ptVsCalIso" ), iso.name( "ptVsCalIso" ), 100, 0., 100., 50, -10., 25.);
  ptVsEcalIso_= fs->make<TH2F>(iso.name( "ptVsEcalIso"), iso.name( "ptVsEcalIso"), 100, 0., 100., 50, -10., 25.);
  ptVsHcalIso_= fs->make<TH2F>(iso.name( "ptVsHcalIso"), iso.name( "ptVsHcalIso"), 100, 0., 100., 50, -10., 25.);
}

/// book for full FW with output stream
void 
MuonKinematic::book(edm::Service<TFileService>& fs, ofstream& file)
{
  NameScheme kin("kin");
  en_ = fs->make<TH1F>(kin.name( file, "en"  ), kin.name("en"  ), 60,   0., 300.);
  pt_ = fs->make<TH1F>(kin.name( file, "pt"  ), kin.name("pt"  ), 30,   0., 150.);
  eta_= fs->make<TH1F>(kin.name( file, "eta" ), kin.name("eta" ), 35, -3.5,  3.5);
  phi_= fs->make<TH1F>(kin.name( file, "phi" ), kin.name("phi" ), 35, -3.5,  3.5);

  NameScheme iso("iso");
  isoJet_  = fs->make<TH1F>(iso.name( file, "isoJet"  ), iso.name("isoJet"  ), 80,   0.,  4.);
  isoJet5_ = fs->make<TH1F>(iso.name( file, "isoJet5" ), iso.name("isoJet5" ), 80,   0.,  4.);
  isoJet10_= fs->make<TH1F>(iso.name( file, "isoJet10"), iso.name("isoJet10"), 80,   0.,  4.);
  isoJet15_= fs->make<TH1F>(iso.name( file, "isoJet15"), iso.name("isoJet15"), 80,   0.,  4.);
  isoJet20_= fs->make<TH1F>(iso.name( file, "isoJet20"), iso.name("isoJet20"), 80,   0.,  4.);
  isoTrk_  = fs->make<TH1F>(iso.name( file, "isoTrk"  ), iso.name("isoTrk"  ), 60,  -1.,  5.);
  isoCal_  = fs->make<TH1F>(iso.name( file, "isoCal"  ), iso.name("isoCal"  ), 40, -10., 30.);
  isoEcal_ = fs->make<TH1F>(iso.name( file, "isoEcal" ), iso.name("isoEcal" ), 40, -10., 30.);
  isoHcal_ = fs->make<TH1F>(iso.name( file, "isoHcal" ), iso.name("isoHcal" ), 40, -10., 30.);
  dRTrkPt_ = fs->make<TH1F>(iso.name( file, "dRTrkPt" ), iso.name("dRTrkPt" ), 42, -0.1,  2.);
  dRTrkN_  = fs->make<TH1F>(iso.name( file, "dRTrkN"  ), iso.name("dRTrkN"  ), 42, -0.1,  2.);
  dRCalPt_ = fs->make<TH1F>(iso.name( file, "dRCalPt" ), iso.name("dRCalPt" ), 42, -0.1,  2.);
  dRCalN_  = fs->make<TH1F>(iso.name( file, "dRCalN"  ), iso.name("dRCalN"  ), 42, -0.1,  2.);
  isoTrkN_ = fs->make<TH1F>(iso.name( file, "isoTrkN" ), iso.name("isoTrkN" ), 21,  -1., 20.);
  isoCalN_ = fs->make<TH1F>(iso.name( file, "isoCalN" ), iso.name("isoCaloN"), 31,  -1., 30.);
  closestCtf_= fs->make<TH1F>(iso.name( file, "closestCtf" ), iso.name("closestCtf" ), 25, -6., -1.);

  ptVsTrkIso_ = fs->make<TH2F>(iso.name( "ptVsTrkIso" ), iso.name( "ptVsTrkIso" ), 100, 0., 100., 50,   0., 25.);
  ptVsCalIso_ = fs->make<TH2F>(iso.name( "ptVsCalIso" ), iso.name( "ptVsCalIso" ), 100, 0., 100., 50, -10., 25.);
  ptVsEcalIso_= fs->make<TH2F>(iso.name( "ptVsEcalIso"), iso.name( "ptVsEcalIso"), 100, 0., 100., 50, -10., 25.);
  ptVsHcalIso_= fs->make<TH2F>(iso.name( "ptVsHcalIso"), iso.name( "ptVsHcalIso"), 100, 0., 100., 50, -10., 25.);
}

/// write to file and free allocated space for FWLite
void 
MuonKinematic::write(const char* filename, const char* directory)
{
  /// save histograms to file
  TFile outFile( filename, "recreate" );
  outFile.mkdir( directory );
  outFile.cd( directory );

  /// basic kinematic
  en_ ->Write( );
  pt_ ->Write( );
  eta_->Write( );
  phi_->Write( );

  /// isolation
  isoJet_  ->Write( );
  isoJet5_ ->Write( );
  isoJet10_->Write( );
  isoJet15_->Write( );
  isoJet20_->Write( );
  isoTrk_  ->Write( );
  isoCal_  ->Write( );
  isoEcal_ ->Write( );
  isoHcal_ ->Write( );
  dRTrkPt_ ->Write( );
  dRTrkN_  ->Write( );
  dRCalPt_ ->Write( );
  dRCalN_  ->Write( );
  closestCtf_ ->Write( );

  /// correlations
  ptVsTrkIso_ ->Write( );
  ptVsCalIso_ ->Write( );
  ptVsEcalIso_->Write( );
  ptVsHcalIso_->Write( );

  outFile.Close();

  // free allocated space
  delete en_;
  delete pt_;
  delete eta_;
  delete phi_;

  /// isolation
  delete isoJet_;
  delete isoJet5_;
  delete isoJet10_;
  delete isoJet15_;
  delete isoJet20_;
  delete isoTrk_;
  delete isoCal_;
  delete isoEcal_;
  delete isoHcal_;
  delete dRTrkPt_;
  delete dRTrkN_;
  delete dRCalPt_;
  delete dRCalN_;
  delete closestCtf_;

  /// correlations
  delete ptVsTrkIso_;
  delete ptVsCalIso_;
  delete ptVsEcalIso_;
  delete ptVsHcalIso_;
}
