
#include "TopAnalysis/TopAnalyzer/interface/MixedObjectsAnalyzer.h"

#include "stdio.h"
#include <TString.h>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

/// default constructor
MixedObjectsAnalyzer::MixedObjectsAnalyzer(const edm::ParameterSet& cfg) :
  JetSrc_            (cfg.getParameter<edm::InputTag>("JetSrc"      )),
  METSrc_            (cfg.getParameter<edm::InputTag>("METSrc"      )),
  MuonSrc_           (cfg.getParameter<edm::InputTag>("MuonSrc"     )),
  ElectronSrc_       (cfg.getParameter<edm::InputTag>("ElectronSrc" )),
  GenJetSrc_         (cfg.getParameter<edm::InputTag>("GenJetSrc"   )),
  addGenJetSrc_      (cfg.getParameter<edm::InputTag>("addGenJetSrc")),
  GenMETSrc_         (cfg.getParameter<edm::InputTag>("GenMETSrc"   )),
  GenLepSrc_         (cfg.getParameter<edm::InputTag>("GenLepSrc"   )),
  ingenPS_           (cfg.getParameter<edm::InputTag>("ingenPS"     )),
  weight_            (cfg.getParameter<edm::InputTag>("weight"      )),
  VertexSrc_         (cfg.getParameter<edm::InputTag>("VertexSrc"   )),
  semiLepEvt_        (cfg.getParameter<edm::InputTag>("semiLepEvent")),
  hypoKey_           (cfg.getParameter<std::string>  ("hypoKey"     )),
  btagAlgo_          (cfg.getParameter<std::string>  ("btagAlgo"    )),
  btagDiscr_         (cfg.getParameter<double>       ("btagDiscr"   )),
  addJetPt_          (cfg.getParameter<double>       ("addJetPt"    ))
{
}

/// default destructor
MixedObjectsAnalyzer::~MixedObjectsAnalyzer()
{
}

/// initiate histograms
void
MixedObjectsAnalyzer::beginJob()
{
  // define TFileService which takes care of writing contents to root file
  edm::Service<TFileService> fs;
  if( !fs ) throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

  // define the tree and make it known to the TFileService
  tree = fs->make<TTree>("tree","tree",0);

  // event identifiers
  runNumber = 0;
  tree->Branch("runNumber", &runNumber, "runNumber/i");
  luminosityBlockNumber = 0;
  tree->Branch("luminosityBlockNumber", &luminosityBlockNumber, "luminosityBlockNumber/i");
  eventNumber= 0;
  tree->Branch("eventNumber", &eventNumber, "eventNumber/i");
  // PDF related info
  Q   =  0.;
  id1 = -42;
  id2 = -42;
  x1  =  0.;
  x2  =  0.;
  tree->Branch("Q"  , &Q  , (std::string("Q"  ) + "/F").c_str());
  tree->Branch("id1", &id1, (std::string("id1") + "/I").c_str());
  tree->Branch("id2", &id2, (std::string("id2") + "/I").c_str());
  tree->Branch("x1" , &x1 , (std::string("x1" ) + "/D").c_str());
  tree->Branch("x2" , &x2 , (std::string("x2" ) + "/D").c_str());

  // event weight
  tree->Branch("weight", &weight, "weight/F");

  // parton level quark vs reco jet assignment
  tree->Branch("qAssignment" , &valueAssignment , "qAssignment/F" );

  // different invariant masses
  MuNu4J=-1000;
  ElNu4J=-1000;
  mJJ=-1000; 
  mWJJ=-1000;
  mWFitJJ=-1000;
  mHbb=-1000;         
  tree->Branch("MWJJ"      , &mWJJ      , "MWJJ/F"      );
  tree->Branch("MWFitJJ"   , &mWFitJJ   , "MWFitJJ/F"   );
  tree->Branch("MuNu4J"    , &MuNu4J    , "MuNu4J/F"    );
  tree->Branch("ElNu4J"    , &ElNu4J    , "MuNu4J/F"    );                
  tree->Branch("mHbb"      , &mHbb      , "MHbb/F"      ); 

  // multiplicities
  Nbjets=-1000;
  Njets=-1000;
  NjetsTrue=-1000;
  tree->Branch("Nbjets"    , &Nbjets    , "Nbjets/F"    );
  tree->Branch("Njets"     , &Njets     , "Njets/F"     );
  tree->Branch("NjetsTrue" , &NjetsTrue , "NjetsTrue/F" );
 
  // leading jet which is not associated with the ttbar system
  tree->Branch("leadNonttjetPt" , &leadNonttjetPt , "leadNonttjetPt/F" );   
  tree->Branch("leadNonttjetY"  , &leadNonttjetY  , "leadNonttjetY/F"  );   
  tree->Branch("leadNonttjetEta", &leadNonttjetEta, "leadNonttjetEta/F"); 
  tree->Branch("leadNonttjetPtTrue" , &leadNonttjetPtTrue , "leadNonttjetPtTrue/F" );   
  tree->Branch("leadNonttjetYTrue"  , &leadNonttjetYTrue  , "leadNonttjetYTrue/F"  );   
  tree->Branch("leadNonttjetEtaTrue", &leadNonttjetEtaTrue, "leadNonttjetEtaTrue/F");  
  // m(ttbar+jet)
  tree->Branch("ttbarJetMass"    , &ttbarJetMass    , "ttbarJetMass/F"    ); 
  tree->Branch("ttbarJetMassTrue", &ttbarJetMassTrue, "ttbarJetMassTrue/F"); 
  tree->Branch("rhos"            , &rhos            , "rhos/F"            ); 
  tree->Branch("rhosTrue"        , &rhosTrue        , "rhosTrue/F"        ); 

  // object kinematics to check kinfit shift
  tree->Branch("bqhadPtPre"   ,&bqhadPtPre , "bqhadPtPre/F" ); 
  tree->Branch("bqhadEtaPre"  ,&bqhadEtaPre, "bqhadEtaPre/F"); 
  tree->Branch("bqhadPhiPre"  ,&bqhadPhiPre, "bqhadEtaPre/F"); 
  tree->Branch("bqlepPtPre"   ,&bqlepPtPre , "bqlepPtPre/F" ); 
  tree->Branch("bqlepEtaPre"  ,&bqlepEtaPre, "bqlepEtaPre/F"); 
  tree->Branch("bqlepPhiPre"  ,&bqlepPhiPre, "bqlepPhiPre/F"); 
  tree->Branch("lqPtPre"      ,&lqPtPre    , "lqPtPre/F"    ); 
  tree->Branch("lqEtaPre"     ,&lqEtaPre   , "lqEtaPre/F"   ); 
  tree->Branch("lqPhiPre"     ,&lqPhiPre   , "lqPhiPre/F"   ); 
  tree->Branch("lqbarPtPre"   ,&lqbarPtPre , "lqbarPtPre/F" ); 
  tree->Branch("lqbarEtaPre"  ,&lqbarEtaPre, "lqbarEtaPre/F"); 
  tree->Branch("lqbarPhiPre"  ,&lqbarPhiPre, "lqbarPhiPre/F"); 
  tree->Branch("bqhadPtFit"   ,&bqhadPtFit , "bqhadPtFit/F" ); 
  tree->Branch("bqhadEtaFit"  ,&bqhadEtaFit, "bqhadEtaFit/F"); 
  tree->Branch("bqhadPhiFit"  ,&bqhadPhiFit, "bqhadPhiFit/F"); 
  tree->Branch("bqlepPtFit"   ,&bqlepPtFit , "bqlepPtFit/F" ); 
  tree->Branch("bqlepEtaFit"  ,&bqlepEtaFit, "bqlepEtaFit/F"); 
  tree->Branch("bqlepPhiFit"  ,&bqlepPhiFit, "bqlepPhiFit/F"); 
  tree->Branch("lqPtFit"      ,&lqPtFit    , "lqPtFit/F"    ); 
  tree->Branch("lqEtaFit"     ,&lqEtaFit   , "lqEtaFit/F"   ); 
  tree->Branch("lqPhiFit"     ,&lqPhiFit   , "lqPhiFit/F"   ); 
  tree->Branch("lqbarPtFit"   ,&lqbarPtFit , "lqbarPtFit/F" ); 
  tree->Branch("lqbarEtaFit"  ,&lqbarEtaFit, "lqbarEtaFit/F"); 
  tree->Branch("lqbarPhiFit"  ,&lqbarPhiFit, "lqbarPhiFit/F"); 
  tree->Branch("nuPtPre"      ,&nuPtPre    , "nuPtPre/F"    ); 
  tree->Branch("nuEtaPre"     ,&nuEtaPre   , "nuEtaPre/F"   ); 
  tree->Branch("nuPhiPre"     ,&nuPhiPre   , "nuEtaPre/F"   ); 
  tree->Branch("nuPtTrue"     ,&nuPtTrue   , "nuPtTrue/F"   ); 
  tree->Branch("nuEtaTrue"    ,&nuEtaTrue  , "nuEtaTrue/F"  ); 
  tree->Branch("nuPhiTrue"    ,&nuPhiTrue  , "nuEtaTrue/F"  ); 
  tree->Branch("nuPtFit"      ,&nuPtFit    , "nuPtFit/F"    ); 
  tree->Branch("nuEtaFit"     ,&nuEtaFit   , "nuEtaFit/F"   ); 
  tree->Branch("nuPhiFit"     ,&nuPhiFit   , "nuEtaFit/F"   ); 
  tree->Branch("lepPtPre"     ,&lepPtPre   , "lepPtPre/F"   ); 
  tree->Branch("lepEtaPre"    ,&lepEtaPre  , "lepEtaPre/F"  ); 
  tree->Branch("lepPhiPre"    ,&lepPhiPre  , "lepEtaPre/F"  );
  tree->Branch("lepPtTrue"    ,&lepPtTrue  , "lepPtTrue/F"  ); 
  tree->Branch("lepEtaTrue"   ,&lepEtaTrue , "lepEtaTrue/F" ); 
  tree->Branch("lepPhiTrue"   ,&lepPhiTrue , "lepEtaTrue/F" );
  tree->Branch("lepPtFit"     ,&lepPtFit   , "lepPtFit/F"   ); 
  tree->Branch("lepEtaFit"    ,&lepEtaFit  , "lepEtaFit/F"  ); 
  tree->Branch("lepPhiFit"    ,&lepPhiFit  , "lepEtaFit/F"  );
  tree->Branch("sumEtPre"     ,&sumEtPre   , "sumEtPre/F"   ); 
  tree->Branch("sumEtTrue"    ,&sumEtTrue  , "sumEtTrue/F"  ); 
  // top and ttbar quantities
  tree->Branch("topPtLepFit"  ,&topPtLepFit , "topPtLepFit/F" );
  tree->Branch("topPtHadFit"  ,&topPtHadFit , "topPtHadFit/F" );
  tree->Branch("topYLepFit"   ,&topYLepFit  , "topYLepFit/F"  );
  tree->Branch("topYHadFit"   ,&topYHadFit  , "topYHadFit/F"  );
  tree->Branch("ttbarMassFit" ,&ttbarMassFit, "ttbarMassFit/F");
  tree->Branch("ttbarYFit"    ,&ttbarYFit   , "ttbarYFit/F"   );
  tree->Branch("ttbarPtFit"   ,&ttbarPtFit  , "ttbarPtFit/F"  );
  tree->Branch("topPtLepTrue"  ,&topPtLepTrue , "topPtLepTrue/F" );
  tree->Branch("topPtHadTrue"  ,&topPtHadTrue , "topPtHadTrue/F" );
  tree->Branch("topYLepTrue"   ,&topYLepTrue  , "topYLepTrue/F"  );
  tree->Branch("topYHadTrue"   ,&topYHadTrue  , "topYHadTrue/F"  );
  tree->Branch("ttbarMassTrue" ,&ttbarMassTrue, "ttbarMassTrue/F");
  tree->Branch("ttbarYTrue"    ,&ttbarYTrue   , "ttbarYTrue/F"   );
  tree->Branch("ttbarPtTrue"   ,&ttbarPtTrue  , "ttbarPtTrue/F"  );

  // chi2
  tree->Branch("chi2", &chi2, "chi2/F"  );
  // nPV
  nPV=-1000;
  tree->Branch("nPV", &nPV, "nPV/F");
  // within visible phase space?
  tree->Branch("inVisPS"       , &inVisPS       , "inVisPS/O"      );

  // histos
  // 1D
  hists_["MuNu4J" ] = fs->make<TH1F>( "MuNu4J"  , "MuNu4J"  , 2500,  0. , 2500 );
  hists_["ElNu4J" ] = fs->make<TH1F>( "ElNu4J"  , "ElNu4J"  , 2500,  0. , 2500 );
  hists_["MJJ"    ] = fs->make<TH1F>( "MJJ"     , "MJJ"     , 2500,  0. , 2500 );        
  hists_["MWJJ"   ] = fs->make<TH1F>( "MWJJ"    , "MWJJ"    , 160 ,  0. , 160  );        
  hists_["MWFitJJ"] = fs->make<TH1F>( "MWFitJJ" , "MWFitJJ" , 160 ,  0. , 160  );       
  hists_["mHbb"   ] = fs->make<TH1F>( "mHbb"    , "mHbb"    , 500 ,  0. , 500  );      
  hists_["mbb"    ] = fs->make<TH1F>( "mbb"     , "mbb"     , 500 ,  0. , 500  );  
  hists_["Nbjets" ] = fs->make<TH1F>( "Nbjets"  , "Nbjets"  , 10  ,  -0.5, 9.5 ); 
  hists_["Njets"  ] = fs->make<TH1F>( "Njets"   , "Njets"   , 15  ,  -0.5, 14.5);  
  hists_["Ngenjets"]= fs->make<TH1F>( "Ngenjets", "Ngenjets", 15  ,  -0.5, 14.5);  
  hists_["leadNonttjetPt" ] = fs->make<TH1F>( "leadNonttjetPt" , "leadNonttjetPt" , 1500, 0 , 1500.);
  hists_["leadNonttjetY"  ] = fs->make<TH1F>( "leadNonttjetY"  , "leadNonttjetY"  , 100 ,-5., 5    );
  hists_["leadNonttjetEta"] = fs->make<TH1F>( "leadNonttjetEta", "leadNonttjetEta", 100 ,-5., 5    );
  hists_["leadNonttjetPtTrue" ] = fs->make<TH1F>( "leadNonttjetPtTrue" , "leadNonttjetPtTrue" , 1500, 0 , 1500.);
  hists_["leadNonttjetYTrue"  ] = fs->make<TH1F>( "leadNonttjetYTrue"  , "leadNonttjetYTrue"  , 100 ,-5., 5    );
  hists_["leadNonttjetEtaTrue"] = fs->make<TH1F>( "leadNonttjetEtaTrue", "leadNonttjetEtaTrue", 100 ,-5., 5    );
  hists_["shiftBqPt"    ] = fs->make<TH1F>( "shiftBqPt" , "shiftBqPt"   , 2000,  -100. , 100 );
  hists_["shiftLqPt"    ] = fs->make<TH1F>( "shiftLqPt" , "shiftLqPt"   , 2000,  -100. , 100 );
  hists_["shiftBqEta"   ] = fs->make<TH1F>( "shiftBqEta", "shiftBqEta"  , 3000,  -0.3  , 0.3 );
  hists_["shiftLqEta"   ] = fs->make<TH1F>( "shiftLqEta", "shiftLqEta"  , 3000,  -0.3  , 0.3 );
  hists_["shiftBqPhi"   ] = fs->make<TH1F>( "shiftBqPhi", "shiftBqPhi"  , 2000,  -0.2  , 0.2 );
  hists_["shiftLqPhi"   ] = fs->make<TH1F>( "shiftLqPhi", "shiftLqPhi"  , 2000,  -0.2  , 0.2 );
  hists_["shiftLepPt"   ] = fs->make<TH1F>( "shiftLepPt" , "shiftLepPt" , 2000,  -2.   , 2.  );
  hists_["shiftLepEta"  ] = fs->make<TH1F>( "shiftLepEta", "shiftLepEta", 1000,  -0.000005, 0.000005);
  hists_["shiftLepPhi"  ] = fs->make<TH1F>( "shiftLepPhi", "shiftLepPhi", 1000,  -0.000002, 0.000002);
  hists_["shiftNuPt"    ] = fs->make<TH1F>( "shiftNuPt" , "shiftNuPt"   , 2000,  -200. , 200 );
  hists_["shiftNuEta"   ] = fs->make<TH1F>( "shiftNuEta", "shiftNuEta"  , 1100,  -5.5  , 5.5 );
  hists_["shiftNuPhi"   ] = fs->make<TH1F>( "shiftNuPhi", "shiftNuPhi"  , 1500,  -1.5  , 1.5 );
  hists_["nPV"          ] = fs->make<TH1F>("nPV"          , "nPV"          , 50, 0., 50.);
  hists_["nPVunweighted"] = fs->make<TH1F>("nPVunweighted", "nPVunweighted", 50, 0., 50.);
  hists_["nPV"          ]->GetXaxis()->SetTitle(VertexSrc_.label().c_str());
  hists_["nPVunweighted"]->GetXaxis()->SetTitle(VertexSrc_.label().c_str());
  hists_["rhos"         ] = fs->make<TH1F>("rhos"   , "rhos"   , 110,  0. , 1.1);
  hists_["rhosGen"      ] = fs->make<TH1F>("rhosGen", "rhosGen", 110,  0. , 1.1);

  // 2D
  hists2D_["Njets_"  ] = fs->make<TH2F>( "Njets_" , "Njets_" , 15, -0.5, 14.5, 15, -0.5, 14.5);  
  hists2D_["rhos_"   ] = fs->make<TH2F>( "rhos_"  , "rhos_"  , 110,  0. , 1.1, 110,  0. , 1.1);
}

/// fill variables
void
MixedObjectsAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
  // activate debug printout
  bool debug=false;

  // Event identifier quantities
  edm::EventAuxiliary aux = event.eventAuxiliary();
  runNumber             = aux.run();
  luminosityBlockNumber = aux.luminosityBlock();
  eventNumber           = aux.event();

  if(debug) std::cout << "run:lumi:evt=" << runNumber << ":" << luminosityBlockNumber << ":" << eventNumber << std::endl;

  // PDF related quantities
  edm::Handle<GenEventInfoProduct> genInfo; 
  event.getByLabel("generator", genInfo);  
  Q   =  0.; id1 = -42; id2 = -42; x1  =  0.; x2  =  0.;
  if(genInfo.isValid()){
    Q     = genInfo->pdf()->scalePDF;
    id1   = genInfo->pdf()->id.first;
    id2   = genInfo->pdf()->id.second;
    x1    = genInfo->pdf()->x.first;
    x2    = genInfo->pdf()->x.second;
  }

  // other collections
  edm::Handle<edm::View< pat::MET > > MET_h;
  event.getByLabel(METSrc_, MET_h);
  
  edm::Handle<edm::View< pat::Muon > > muons_h;
  event.getByLabel(MuonSrc_, muons_h);
  
  edm::Handle<edm::View< pat::Electron > > electrons_h;
  event.getByLabel(ElectronSrc_, electrons_h);
  
  edm::Handle<edm::View< pat::Jet > > jets_h;
  event.getByLabel(JetSrc_, jets_h);
  
  edm::Handle<edm::View< reco::GenJet > > genJets_h;
  event.getByLabel(GenJetSrc_, genJets_h);

  edm::Handle<edm::View< reco::GenJet > > addGenJets_h;
  event.getByLabel(addGenJetSrc_, addGenJets_h);

  edm::Handle<edm::View< reco::GenMET > > genMET_h;
  event.getByLabel(GenMETSrc_, genMET_h);
  
  edm::Handle<edm::View< reco::GenParticle > > genLep_h;
  event.getByLabel(GenLepSrc_, genLep_h);
  
  edm::Handle<std::vector<reco::Vertex> > vertecies_h;
  event.getByLabel(VertexSrc_, vertecies_h);
  
  edm::Handle<double> weight_h;
  event.getByLabel(weight_, weight_h);

  edm::Handle<TtSemiLeptonicEvent> semiLepEvt_h;
  event.getByLabel(semiLepEvt_, semiLepEvt_h);

  edm::Handle<bool > visPS;
  if(genLep_h.isValid()) event.getByLabel(ingenPS_, visPS);
  inVisPS=false;
  if(genLep_h.isValid()&&!visPS.failedToGet()&&visPS.isValid()&&*visPS) inVisPS=true;
  
  // event weight
  weight = weight_h.isValid() ? *weight_h : 1.;

  // information concerning PU/NPV   
  nPV = vertecies_h.isValid() ? vertecies_h->size() : -1000;
  hists_["nPV"          ]->Fill(nPV, weight);
  hists_["nPVunweighted"]->Fill(nPV, 1.    );
  
  // collect information from KinFitHypothesis
  int lepBIndex=-1000, hadBIndex=-1000, lightQIndex=-1000, lightQBarIndex=-1000;
  bqhadPtFit=bqhadEtaFit=bqhadPhiFit=bqlepPtFit=bqlepEtaFit=bqlepPhiFit=lqPtFit=lqEtaFit=lqPhiFit=lqbarPtFit=lqbarEtaFit=lqbarPhiFit=chi2=-1000;
  lepPtFit=lepEtaFit=lepPhiFit=lepPtTrue=lepEtaTrue=lepPhiTrue=nuPtFit=nuEtaFit=nuPhiFit=nuPtTrue=nuEtaTrue=nuPhiTrue=-1000;
  if( semiLepEvt_h.isValid()&&semiLepEvt_h->isHypoValid(hypoKey_) ){
    // jet indices
    lepBIndex      = semiLepEvt_h->jetLeptonCombination(hypoKey_)[TtSemiLepEvtPartons::LepB     ];
    hadBIndex      = semiLepEvt_h->jetLeptonCombination(hypoKey_)[TtSemiLepEvtPartons::HadB     ];
    lightQIndex    = semiLepEvt_h->jetLeptonCombination(hypoKey_)[TtSemiLepEvtPartons::LightQ   ];
    lightQBarIndex = semiLepEvt_h->jetLeptonCombination(hypoKey_)[TtSemiLepEvtPartons::LightQBar];
    // jet kinematics after fit
    bqhadPtFit  =semiLepEvt_h->hadronicDecayB(hypoKey_)->pt();
    bqhadEtaFit =semiLepEvt_h->hadronicDecayB(hypoKey_)->eta();
    bqhadPhiFit =semiLepEvt_h->hadronicDecayB(hypoKey_)->phi();
    bqlepPtFit  =semiLepEvt_h->leptonicDecayB(hypoKey_)->pt();
    bqlepEtaFit =semiLepEvt_h->leptonicDecayB(hypoKey_)->eta();
    bqlepPhiFit =semiLepEvt_h->leptonicDecayB(hypoKey_)->phi();
    lqPtFit     =semiLepEvt_h->hadronicDecayQuark(hypoKey_)->pt();
    lqEtaFit    =semiLepEvt_h->hadronicDecayQuark(hypoKey_)->eta();
    lqPhiFit    =semiLepEvt_h->hadronicDecayQuark(hypoKey_)->phi();
    lqbarPtFit  =semiLepEvt_h->hadronicDecayQuarkBar(hypoKey_)->pt();
    lqbarEtaFit =semiLepEvt_h->hadronicDecayQuarkBar(hypoKey_)->eta();
    lqbarPhiFit =semiLepEvt_h->hadronicDecayQuarkBar(hypoKey_)->phi();
    lepPtFit  =semiLepEvt_h->singleLepton(hypoKey_)->pt();
    lepEtaFit =semiLepEvt_h->singleLepton(hypoKey_)->eta();
    lepPhiFit =semiLepEvt_h->singleLepton(hypoKey_)->phi();
    nuPtFit   =semiLepEvt_h->singleNeutrino(hypoKey_)->pt();
    nuEtaFit  =semiLepEvt_h->singleNeutrino(hypoKey_)->eta();
    nuPhiFit  =semiLepEvt_h->singleNeutrino(hypoKey_)->phi();
    chi2=semiLepEvt_h->fitChi2();
  }
  if( semiLepEvt_h.isValid()&&semiLepEvt_h->singleLepton()&&semiLepEvt_h->singleNeutrino() ){
    lepPtTrue =semiLepEvt_h->singleLepton()->pt();
    lepEtaTrue=semiLepEvt_h->singleLepton()->eta();
    lepPhiTrue=semiLepEvt_h->singleLepton()->phi();
    nuPtTrue  =semiLepEvt_h->singleNeutrino()->pt();
    nuEtaTrue =semiLepEvt_h->singleNeutrino()->eta();
    nuPhiTrue =semiLepEvt_h->singleNeutrino()->phi();
  }
  // parton level quark vs reco jet assignment
  valueAssignment=-1000.;
  if(genLep_h.isValid()&&!visPS.failedToGet()&&visPS.isValid()&&semiLepEvt_h.isValid()) valueAssignment=checkPartonAssignment(semiLepEvt_h, 5); // FIXME: N(jets) in KinFit hardcoded for the moment
  if(debug) std::cout << "reco jet assignment (kinFit) wrt reco jet-quark matching: " << valueAssignment << std::endl;
  // invariant masses (measure for mttbar)
  MuNu4J=-1000;
  ElNu4J=-1000;
  if(jets_h.isValid()&&MET_h.isValid()&&muons_h.isValid()&&electrons_h.isValid()){
    if(jets_h->size()>3&&MET_h->size()>0){
      if(muons_h->size()>0){
	MuNu4J=(muons_h->at(0).p4()+MET_h->at(0).p4()+jets_h->at(0).p4()+jets_h->at(1).p4()+jets_h->at(2).p4()+jets_h->at(3).p4()).mass();
	hists_["MuNu4J"]->Fill(MuNu4J, weight);
	if(debug) std::cout << "mu: " << MuNu4J << std::endl;
      }
      if(electrons_h->size()>0){ 
	ElNu4J=(electrons_h->at(0).p4()+MET_h->at(0).p4()+jets_h->at(0).p4()+jets_h->at(1).p4()+jets_h->at(2).p4()+jets_h->at(3).p4()).mass();
	hists_["ElNu4J"]->Fill(ElNu4J, weight);
	if(debug) std::cout << "el: " << ElNu4J << std::endl;
      }
    }
  }
  mJJ=-1000;
  mWJJ=-1000;
  mWFitJJ=-1000;
  mHbb=-1000;
  BindexA=-1000;
  BindexB=-1000;
  BindexC=-1000;
  BindexD=-1000;
  Nbjets=-1000;
  Njets=-1000;
  leadNonttjetIX=-1000;
  leadNonttjetPt=leadNonttjetY=leadNonttjetEta=-1000;
  bqhadPtPre=bqhadEtaPre=bqhadPhiPre=bqlepPtPre=bqlepEtaPre=bqlepPhiPre=lqPtPre=lqEtaPre=lqPhiPre=lqbarPtPre=lqbarEtaPre=lqbarPhiPre=-1000;
  // loop jets
  if(jets_h.isValid()){
    Nbjets=0;
    Njets=0;
    for(edm::View< pat::Jet>::const_iterator jetsA=jets_h->begin(); jetsA!=jets_h->end(); ++jetsA){
      int tempindex=jetsA-jets_h->begin();
      // check if kinfit information is available
      if(semiLepEvt_h.isValid()&&semiLepEvt_h->isHypoValid(hypoKey_)){
	// jet kinematics before the fit
	if(tempindex==hadBIndex){
	  bqhadPtPre =jetsA->pt();
	  bqhadEtaPre=jetsA->eta();
	  bqhadPhiPre=jetsA->phi();
	}
	else if(tempindex==lepBIndex){
	  bqlepPtPre =jetsA->pt();
	  bqlepEtaPre=jetsA->eta();
	  bqlepPhiPre=jetsA->phi();
	}
	else if(tempindex==lightQIndex){
	  lqPtPre =jetsA->pt();
	  lqEtaPre=jetsA->eta();
	  lqPhiPre=jetsA->phi();
	}
	else if(tempindex==lightQBarIndex){
	  lqbarPtPre =jetsA->pt();
	  lqbarEtaPre=jetsA->eta();
	  lqbarPhiPre=jetsA->phi();
	}
	// search for leading non ttbar jet within the 5 leading jets
	if(tempindex<=5&&tempindex!=lepBIndex&&tempindex!=hadBIndex&&tempindex!=lightQIndex&&tempindex!=lightQBarIndex){
	  leadNonttjetIX=tempindex;
	  leadNonttjetPt=jetsA->pt();
	  leadNonttjetY =jetsA->rapidity();
	  leadNonttjetEta =jetsA->eta();
	}
      }
      // count jet multiplicity
      ++Njets;
      // count b jet multiplicity (CSVM algo)
      if(jetsA->bDiscriminator(btagAlgo_) >btagDiscr_){
	++Nbjets;
	// store b-jet index 
	if     (BindexA<=-1) BindexA=tempindex;
	else if(BindexB<=-1) BindexB=tempindex;
	else if(BindexC<=-1) BindexC=tempindex;
	else if(BindexD<=-1) BindexD=tempindex;
      }
      // 2nd loop for jet-jet masses
      for(edm::View< pat::Jet>::const_iterator jetsB=jets_h->begin(); jetsB!=jets_h->end(); ++jetsB){      
	// take care of double counting
	if(jetsB>jetsA){
	  // a) invariant jet-jet mass for all combinations
	  mJJ=(jetsA->p4()+jetsB->p4()).mass();
	  // b) for combination closest to w mass
	  if(mWJJ<=-1) mWJJ=mJJ;
	  else if(std::abs(mJJ-80.4)<std::abs(mWJJ-80.4)) mWJJ=mJJ;
	  // c) for invariant bjet-bjet masses
	  if(jetsA->bDiscriminator(btagAlgo_) >btagDiscr_&&jetsB->bDiscriminator(btagAlgo_) >btagDiscr_){
	    hists_["mbb"]->Fill((jetsA->p4()+jetsB->p4()).mass(), weight);
	  }
	  // d) for all invariant masses hists_["MJJ"]->Fill(mJJ, weight);
	}
      }
    }
  }
  //std::cout << "loping genjets" << std::endl;

  // calculate N(genjets)
  NjetsTrue=-1000;
  if(genJets_h.isValid()&&inVisPS){
    NjetsTrue=0;
    // loop gen Jets
    for(edm::View< reco::GenJet >::const_iterator genJets=genJets_h->begin(); genJets!=genJets_h->end(); ++genJets){
      if(genJets->pt()>30&&genJets->eta()<2.4){
	++NjetsTrue;
      }
    }
  }

  // fill multiplicity plots
  if(Nbjets!=-1000   ) hists_["Nbjets"]->Fill(Nbjets, weight);
  hists_["Nbjets"]->GetXaxis()->SetTitle("N("+TString(btagAlgo_)+" b-jets)");
  if(Njets!=-1000    ) hists_["Njets" ]->Fill(Njets , weight);
  if(debug) std::cout << "fill gen jet plot" << std::endl;
  if(NjetsTrue!=-1000) hists_["Ngenjets"]->Fill(NjetsTrue, weight); 
  if(debug) std::cout << "fill gen-reco jet plot" << std::endl;
  if(Njets!=-1000    ) hists2D_["Njets_"]->Fill(NjetsTrue, Njets, weight); 
  if(debug) std::cout << "done" << std::endl;

  // lead non ttbar jet
  if(debug) std::cout << "fill leading non ttbar jet plot" << std::endl;
  if(leadNonttjetIX!=-1000){
    hists_["leadNonttjetPt" ]->Fill(leadNonttjetPt , weight);
    hists_["leadNonttjetY"  ]->Fill(leadNonttjetY  , weight);
    hists_["leadNonttjetEta"]->Fill(leadNonttjetEta, weight);
  }
  if(debug) std::cout << "done" << std::endl;

  // jet-jet invariant mass closest to w mass
  if(debug) std::cout << "fill m(WJJ)" << std::endl;
  if(mWJJ>0   ) hists_["MWJJ"   ]->Fill(mWJJ, weight);   
  if(debug) std::cout << "done" << std::endl;

  // jet-jet invariant mass from jets associated to W decay by kinfit
  if(debug) std::cout << "fill m(WJJ) from fit" << std::endl;
  if(lightQIndex!=-1000&&lightQBarIndex!=-1000) mWFitJJ=(jets_h->at(lightQIndex).p4()+jets_h->at(lightQBarIndex).p4()).mass();
  if(mWFitJJ>0) hists_["MWFitJJ"]->Fill(mWFitJJ, weight);  
  if(debug) std::cout << "done" << std::endl;

  // bjet-bjet invariant mass for bjet pair comming NOT from ttbar decay (using kinfit information) in 4 b events
  if(debug) std::cout << "fill m(bb) not from ttbar" << std::endl;
  if(Nbjets>=4&&lepBIndex!=-1000&&hadBIndex!=-1000){
    // search non ttbar bjets
    int HbbindexA=-1000, HbbindexB=-1000;
    if(debug){
      std::cout << "lepBIndex: " << lepBIndex << std::endl;
      std::cout << "hadBIndex: " << hadBIndex << std::endl;
      std::cout << "BindexA: " << BindexA << std::endl;
      std::cout << "BindexB: " << BindexB << std::endl;
      std::cout << "BindexC: " << BindexC << std::endl;
      std::cout << "BindexD: " << BindexD << std::endl;
    }
    if(BindexA!=lepBIndex&&BindexA!=hadBIndex){
      if     (HbbindexA<=-1000) HbbindexA=BindexA;
      else if(HbbindexB<=-1000) HbbindexB=BindexA;
    }
    if(BindexB!=lepBIndex&&BindexB!=hadBIndex){
      if     (HbbindexA<=-1000) HbbindexA=BindexB;
      else if(HbbindexB<=-1000) HbbindexB=BindexB;
    }
    if(BindexC!=lepBIndex&&BindexC!=hadBIndex){
      if     (HbbindexA<=-1000) HbbindexA=BindexC;
      else if(HbbindexB<=-1000) HbbindexB=BindexC;
    }
    if(BindexD!=lepBIndex&&BindexD!=hadBIndex){
      if(     HbbindexA<=-1000) HbbindexA=BindexD;
      else if(HbbindexB<=-1000) HbbindexB=BindexD;
    }
    if(debug){
      std::cout << "HbbindexA: " << HbbindexA << std::endl;
      std::cout << "HbbindexB: " << HbbindexB << std::endl;
    }
    if(HbbindexA!=-1000&&HbbindexB!=-1000) mHbb=(jets_h->at(HbbindexA).p4()+jets_h->at(HbbindexB).p4()).mass();
    if(mHbb!=-1000) hists_["mHbb"]->Fill(mHbb, weight); 
  }
  if(debug) std::cout << "done" << std::endl;
    
  // fill kinfit shift histos
  if(debug) std::cout << "fill kinfit shift histos" << std::endl;
  if(bqhadPtFit!=-1000&&bqhadPtPre!=-1000) hists_["shiftBqPt" ]->Fill( bqhadPtFit-bqhadPtPre, weight);
  if(bqlepPtFit!=-1000&&bqlepPtPre!=-1000) hists_["shiftBqPt" ]->Fill( bqlepPtFit-bqlepPtPre, weight);
  if(lqPtFit!=-1000&&lqPtPre!=-1000&&lqbarPtFit!=-1000&&lqbarPtPre!=-1000){
    // consider ambiguity for light jets
    if(std::abs(lqPtFit-lqPtPre)+std::abs(lqbarPtFit-lqbarPtPre)<std::abs(lqPtFit-lqbarPtPre)+std::abs(lqbarPtFit-lqPtPre)){
      hists_["shiftLqPt"    ]->Fill( lqPtFit-lqPtPre, weight);
      hists_["shiftLqPt"    ]->Fill( lqbarPtFit-lqbarPtPre, weight);
    }
    else{
      hists_["shiftLqPt"    ]->Fill( lqPtFit-lqbarPtPre, weight);
      hists_["shiftLqPt"    ]->Fill( lqbarPtFit-lqPtPre, weight);
    }
  }
  if(bqhadEtaFit!=-1000&&bqhadEtaPre!=-1000) hists_["shiftBqEta" ]->Fill( bqhadEtaFit-bqhadEtaPre, weight);
  if(bqlepEtaFit!=-1000&&bqlepEtaPre!=-1000) hists_["shiftBqEta" ]->Fill( bqlepEtaFit-bqlepEtaPre, weight);
  if(lqEtaFit!=-1000&&lqEtaPre!=-1000&&lqbarEtaFit!=-1000&&lqbarEtaPre!=-1000){
    // consider ambiguity for light jets
    if(std::abs(lqEtaFit-lqEtaPre)+std::abs(lqbarEtaFit-lqbarEtaPre)<std::abs(lqEtaFit-lqbarEtaPre)+std::abs(lqbarEtaFit-lqEtaPre)){
      hists_["shiftLqEta"    ]->Fill( lqEtaFit-lqEtaPre, weight);
      hists_["shiftLqEta"    ]->Fill( lqbarEtaFit-lqbarEtaPre, weight);
    }
    else{
      hists_["shiftLqEta"    ]->Fill( lqEtaFit-lqbarEtaPre, weight);
      hists_["shiftLqEta"    ]->Fill( lqbarEtaFit-lqEtaPre, weight);
    }
  }
  if(bqhadPhiFit!=-1000&&bqhadPhiPre!=-1000) hists_["shiftBqPhi" ]->Fill( modTwoPi(bqhadPhiFit-bqhadPhiPre), weight);
  if(bqlepPhiFit!=-1000&&bqlepPhiPre!=-1000) hists_["shiftBqPhi" ]->Fill( modTwoPi(bqlepPhiFit-bqlepPhiPre), weight);
  if(lqPhiFit!=-1000&&lqPhiPre!=-1000&&lqbarPhiFit!=-1000&&lqbarPhiPre!=-1000){
    // consider ambiguity for light jets
    if(std::abs(modTwoPi(lqPhiFit-lqPhiPre))+std::abs(modTwoPi(lqbarPhiFit-lqbarPhiPre))<std::abs(modTwoPi(lqPhiFit-lqbarPhiPre))+std::abs(modTwoPi(lqbarPhiFit-lqPhiPre))){
      hists_["shiftLqPhi"    ]->Fill( modTwoPi(lqPhiFit-lqPhiPre), weight);
      hists_["shiftLqPhi"    ]->Fill( modTwoPi(lqbarPhiFit-lqbarPhiPre), weight);
    }
    else{
      hists_["shiftLqPhi"    ]->Fill( modTwoPi(lqPhiFit-lqbarPhiPre), weight);
      hists_["shiftLqPhi"    ]->Fill( modTwoPi(lqbarPhiFit-lqPhiPre), weight);
    }
  }
  lepPtPre=lepEtaPre=lepPhiPre=nuPtPre=nuEtaPre=nuPhiPre=sumEtPre=sumEtTrue=-1000;
  if(muons_h.isValid()&&(muons_h->end()-muons_h->begin()>0)){
    lepPtPre  =muons_h->at(0).pt();
    lepEtaPre =muons_h->at(0).eta();
    lepPhiPre =muons_h->at(0).phi();
  }
  else if(electrons_h.isValid()&&(electrons_h->end()-electrons_h->begin()>0)){
    lepPtPre  =electrons_h->at(0).pt();
    lepEtaPre =electrons_h->at(0).eta();
    lepPhiPre =electrons_h->at(0).phi();
  }
  if(MET_h.isValid()&&(MET_h->end()-MET_h->begin()>0)){
    nuPtPre   =MET_h->at(0).pt();
    nuEtaPre  =MET_h->at(0).eta();
    nuPhiPre  =MET_h->at(0).phi();
    sumEtPre  =MET_h->at(0).sumEt();
    if(MET_h->begin()->genMET()&&inVisPS){
      sumEtTrue=MET_h->begin()->genMET()->sumEt();
    }
  }
  if(lepPtFit !=-1000&&lepPtPre !=-1000) hists_["shiftLepPt" ]->Fill( lepPtFit -lepPtPre , weight);  
  if(nuPtFit  !=-1000&&nuPtPre  !=-1000) hists_["shiftNuPt"  ]->Fill( nuPtFit  -nuPtPre  , weight);
  if(lepEtaFit!=-1000&&lepEtaPre!=-1000) hists_["shiftLepEta"]->Fill( lepEtaFit-lepEtaPre, weight);  
  if(nuEtaFit !=-1000&&nuEtaPre !=-1000) hists_["shiftNuEta" ]->Fill( nuEtaFit -nuEtaPre , weight);
  if(lepPhiFit!=-1000&&lepPhiPre!=-1000) hists_["shiftLepPhi"]->Fill( modTwoPi(lepPhiFit-lepPhiPre), weight);  
  if(nuPhiFit !=-1000&&nuPhiPre !=-1000) hists_["shiftNuPhi" ]->Fill( modTwoPi( nuPhiFit -nuPhiPre), weight);
  if(debug) std::cout << "done" << std::endl;

  // top and ttbar distributions
  // a) reco
  if(debug) std::cout << "fill reco top histos" << std::endl;
  topPtLepFit=topPtHadFit=topYLepFit=topYHadFit=ttbarMassFit=ttbarYFit=ttbarPtFit=-1000;
  if(semiLepEvt_h.isValid()&&semiLepEvt_h->isHypoValid(hypoKey_) ){
    topPtLepFit=semiLepEvt_h->leptonicDecayTop(hypoKey_)->pt();
    topPtHadFit=semiLepEvt_h->hadronicDecayTop(hypoKey_)->pt();
    topYLepFit =semiLepEvt_h->leptonicDecayTop(hypoKey_)->rapidity();
    topYHadFit =semiLepEvt_h->hadronicDecayTop(hypoKey_)->rapidity();
    ttbarMassFit=(semiLepEvt_h->hadronicDecayTop(hypoKey_)->p4()+semiLepEvt_h->leptonicDecayTop(hypoKey_)->p4()).mass();
    ttbarYFit   =(semiLepEvt_h->hadronicDecayTop(hypoKey_)->p4()+semiLepEvt_h->leptonicDecayTop(hypoKey_)->p4()).Rapidity();
    ttbarPtFit  =(semiLepEvt_h->hadronicDecayTop(hypoKey_)->p4()+semiLepEvt_h->leptonicDecayTop(hypoKey_)->p4()).pt();
  }
  if(debug) std::cout << "done" << std::endl;

  // b) gen
  if(debug) std::cout << "fill gen top histos" << std::endl;
  topPtLepTrue=topPtHadTrue=topYLepTrue=topYHadTrue=ttbarMassTrue=ttbarYTrue=ttbarPtTrue=-1000;
  if(semiLepEvt_h.isValid()&&semiLepEvt_h->leptonicDecayTop()&&semiLepEvt_h->hadronicDecayTop()){  
    topPtLepTrue=semiLepEvt_h->leptonicDecayTop()->pt();
    topPtHadTrue=semiLepEvt_h->hadronicDecayTop()->pt();
    topYLepTrue =semiLepEvt_h->leptonicDecayTop()->rapidity();
    topYHadTrue =semiLepEvt_h->hadronicDecayTop()->rapidity();
    ttbarMassTrue=(semiLepEvt_h->hadronicDecayTop()->p4()+semiLepEvt_h->leptonicDecayTop()->p4()).mass();
    ttbarYTrue   =(semiLepEvt_h->hadronicDecayTop()->p4()+semiLepEvt_h->leptonicDecayTop()->p4()).Rapidity();
    ttbarPtTrue  =(semiLepEvt_h->hadronicDecayTop()->p4()+semiLepEvt_h->leptonicDecayTop()->p4()).pt();
  }
  if(debug) std::cout << "done" << std::endl;

  // tt+jet distribution
  // scale value chosen
  double mtop=170.0; //172.5
  // initialize tree entries
  ttbarJetMass=ttbarJetMassTrue=-1000;
  rhos=rhosTrue=-1000;
  // calculate quantity
  if(debug) std::cout << "calculate m(ttbar+1jet)" << std::endl;
  if(semiLepEvt_h.isValid()){
    if(debug) std::cout << "valid top hypothesis found" << std::endl;
    // gen level
    if(debug) std::cout << "a) gen level" << std::endl;
    if(semiLepEvt_h->hadronicDecayTop()&&semiLepEvt_h->leptonicDecayTop()&&inVisPS){
      if(addGenJets_h.isValid()&&(addGenJets_h->size()>0)&&inVisPS){
	ttbarJetMassTrue=(semiLepEvt_h->hadronicDecayTop()->p4()+semiLepEvt_h->leptonicDecayTop()->p4()+(addGenJets_h->at(0)).p4()).mass();
	rhosTrue=2*mtop/ttbarJetMassTrue;
	if(debug) std::cout << "ttbarJetMassTrue=" << ttbarJetMassTrue << std::endl;
      }
      else if(debug){
	std::cout << "addGenJets_h.isValid(): " << addGenJets_h.isValid() << std::endl;
	std::cout << "addGenJets_h->size(): " << addGenJets_h->size() << std::endl;
	if(addGenJets_h->size()>0) std::cout << "addGenJets_h->at(0).pt(): " << addGenJets_h->at(0).pt() << std::endl;
      }
    }
    // reco level
    if(debug) std::cout << "b) reco level" << std::endl;
    if(semiLepEvt_h->hadronicDecayTop(hypoKey_)&&semiLepEvt_h->leptonicDecayTop(hypoKey_)){
      if( (semiLepEvt_h->isHypoValid(hypoKey_))&&(leadNonttjetIX!=-1000) ){
	ttbarJetMass=(semiLepEvt_h->hadronicDecayTop(hypoKey_)->p4()+semiLepEvt_h->leptonicDecayTop(hypoKey_)->p4()+(jets_h->at(leadNonttjetIX)).p4()).mass();
	rhos=2*mtop/ttbarJetMass;
	if(debug) std::cout << "ttbarJetMass=" << ttbarJetMass << std::endl;
      }
      else if(debug){
	std::cout << "semiLepEvt_h->isHypoValid(hypoKey_): " << semiLepEvt_h->isHypoValid(hypoKey_) << std::endl;
	std::cout << "leadNonttjetIX: "   << leadNonttjetIX   << std::endl;
	std::cout << "leadNonttjetPt: " << leadNonttjetPt << std::endl;
      }
    }
  }

  // get leading jet kinematics
  if(debug) std::cout << "get gen kinematics of additional jet" << std::endl;
  leadNonttjetPtTrue=leadNonttjetYTrue=leadNonttjetEtaTrue=-1000;
  if(addGenJets_h.isValid()&&(addGenJets_h->size()>0)&&inVisPS){
    leadNonttjetPtTrue=addGenJets_h->at(0).pt();
    leadNonttjetEtaTrue=addGenJets_h->at(0).eta();
    leadNonttjetYTrue=addGenJets_h->at(0).rapidity();
  }
  if(leadNonttjetPt    <=addJetPt_)rhos    =-1000;
  if(leadNonttjetPtTrue<=addJetPt_)rhosTrue=-1000;

  // fill rhos plots
  if(debug) std::cout << "fill ttbar+1jet plots" << std::endl;
  if(rhos    !=-1000&&leadNonttjetPt    >addJetPt_) hists_  ["rhos"   ]->Fill( rhos    , weight);
  if(rhosTrue!=-1000&&leadNonttjetPtTrue>addJetPt_) hists_  ["rhosGen"]->Fill( rhosTrue, weight);
  if(rhos    !=-1000&&leadNonttjetPt    >addJetPt_) hists2D_["rhos_"  ]->Fill( rhosTrue, rhos, weight );
  if(debug) std::cout << "done" << std::endl;

  // lead non ttbar jet GEN
  if(debug) std::cout << "fill leading additional jet plots" << std::endl;
  if(leadNonttjetPtTrue>=0){
    hists_["leadNonttjetPtTrue" ]->Fill(leadNonttjetPtTrue , weight);
    hists_["leadNonttjetYTrue"  ]->Fill(leadNonttjetYTrue  , weight);
    hists_["leadNonttjetEtaTrue"]->Fill(leadNonttjetEtaTrue, weight);
  }
  if(debug) std::cout << "done" << std::endl;

  if(debug) std::cout << "fill tree" << std::endl;
  tree->Fill();
  if(debug) std::cout << "done" << std::endl;

}

void
MixedObjectsAnalyzer::endJob()
{

}

/// helper funtion to take care that delta phi is < pi
double
MixedObjectsAnalyzer::modTwoPi(double DeltaPhi){
  const double pi=TMath::Pi();
  if(DeltaPhi>pi){
    double output=std::abs(2*pi-std::abs(DeltaPhi));
    if(DeltaPhi>0) output*=-1.;
  }
  return DeltaPhi;
}


bool
MixedObjectsAnalyzer::findAncestor(const reco::Candidate* cand, TString decaychain){
  //static int calls = 0;
  //int calls = decaychain.CountChar('<');
  if(cand->numberOfMothers()==0){
    decaychain.Remove(decaychain.Sizeof()-3, 2);
    //resetCount=true;
    //std::cout << decaychain << std::endl;
    // check for light quark candiates hadrons
    TString decaychain2=decaychain;
    int index=decaychain.Index("6");
    decaychain2.Remove(index+1,decaychain.Sizeof()-index-1);
    if((decaychain.Contains("1<-24<-6")&&!decaychain2.Contains("21"))||
       (decaychain.Contains("2<-24<-6")&&!decaychain2.Contains("21"))||
       (decaychain.Contains("3<-24<-6")&&!decaychain2.Contains("21"))||
       (decaychain.Contains("4<-24<-6")&&!decaychain2.Contains("21"))||
       (decaychain.Contains("5<-6")&&!decaychain2.Contains("21")
	&&(decaychain.Contains("511")||
	   decaychain.Contains("521")||
	   decaychain.Contains("513")||
	   decaychain.Contains("523")||
	   decaychain.Contains("515")||
	   decaychain.Contains("525")||
	   decaychain.Contains("531")||
	   decaychain.Contains("533")||
	   decaychain.Contains("535")||
	   decaychain.Contains("541")||
	   decaychain.Contains("543")||
	   decaychain.Contains("545")||
	   decaychain.Contains("551")||
	   decaychain.Contains("553")||
	   decaychain.Contains("555")||
	   decaychain.Contains("557")||
	   decaychain.Contains("5122")||
	   decaychain.Contains("5112")||
	   decaychain.Contains("5212")||
	   decaychain.Contains("5222")||
	   decaychain.Contains("5114")||
	   decaychain.Contains("5214")||
	   decaychain.Contains("5132")||
	   decaychain.Contains("5232")||
	   decaychain.Contains("5312")||
	   decaychain.Contains("5322")||
	   decaychain.Contains("5314")||
	   decaychain.Contains("5324")||
	   decaychain.Contains("5332")||
	   decaychain.Contains("5334")||
	   decaychain.Contains("5142")||
	   decaychain.Contains("5242")||
	   decaychain.Contains("5412")||
	   decaychain.Contains("5422")||
	   decaychain.Contains("5414")||
	   decaychain.Contains("5424")||
	   decaychain.Contains("5342")||
	   decaychain.Contains("5432")||
	   decaychain.Contains("5434")||
	   decaychain.Contains("5442")||
	   decaychain.Contains("5444")||
	   decaychain.Contains("5512")||
	   decaychain.Contains("5522")||
	   decaychain.Contains("5514")||
	   decaychain.Contains("5524")||
	   decaychain.Contains("5532")||
	   decaychain.Contains("5534")||
	   decaychain.Contains("5542")||
	   decaychain.Contains("5544")||
	   decaychain.Contains("5554")
	   ))
       ){
      //std::cout << "possible ttbar jet!" << std::endl;
      std::cout << decaychain << std::endl;
      return true;
    }
  }
  
  // pdgid original particle
  //std::cout.width(5*calls);
  //if(calls==0) std::cout << " " << cand->pdgId() << "(" << cand->numberOfMothers() << " ancestors)" << std::endl;
  // check mothers
  if(cand->numberOfMothers()>0){
    
    // pdgid mothers
    //for(unsigned int numMother=0; numMother<cand->numberOfMothers(); ++numMother){
      //if(numMother==0) { std::cout.width(5*calls); std::cout << " "; std::cout << " <- " << cand->mother(numMother)->pdgId(); }
      //else             { std::cout << ", " << cand->mother(numMother)->pdgId(); }
      //if(numMother==cand->numberOfMothers()-1) std::cout << std::endl; 
      
    //}
    // bubble up mothers
    for(unsigned int numMother=0; numMother<cand->numberOfMothers(); ++numMother){
      TString newchain=decaychain;
      newchain+=std::abs(cand->mother(numMother)->pdgId());
      newchain+="<-";
      findAncestor(cand->mother(numMother), newchain);
    }
  }
  return false;
}

double  
MixedObjectsAnalyzer::checkPartonAssignment(const edm::Handle<TtSemiLeptonicEvent> semiLepEvt, int maxNJets){
  double assignment=-1;
  // ---
  //    check quark assignment
  // ---
  // no parton jet parton match exists
  if(semiLepEvt->isHypoValid("kKinFit")&&(!semiLepEvt->isHypoValid("kGenMatch"))){
    assignment= 9.;
  }
  // if jet parton match exists:
  if(semiLepEvt->isHypoValid("kKinFit")&& semiLepEvt->isHypoValid("kGenMatch")){
    // indices for all quarks from Kinfit Hypothesis and genmatch
    int lepBIndex         = semiLepEvt->jetLeptonCombination("kKinFit"  )[TtSemiLepEvtPartons::LepB     ];
    int hadBIndex         = semiLepEvt->jetLeptonCombination("kKinFit"  )[TtSemiLepEvtPartons::HadB     ];
    int lightQIndex       = semiLepEvt->jetLeptonCombination("kKinFit"  )[TtSemiLepEvtPartons::LightQ   ];
    int lightQBarIndex    = semiLepEvt->jetLeptonCombination("kKinFit"  )[TtSemiLepEvtPartons::LightQBar];
    int lepBIndexGen      = semiLepEvt->jetLeptonCombination("kGenMatch")[TtSemiLepEvtPartons::LepB     ];
    int hadBIndexGen      = semiLepEvt->jetLeptonCombination("kGenMatch")[TtSemiLepEvtPartons::HadB     ];
    int lightQIndexGen    = semiLepEvt->jetLeptonCombination("kGenMatch")[TtSemiLepEvtPartons::LightQ   ];
    int lightQBarIndexGen = semiLepEvt->jetLeptonCombination("kGenMatch")[TtSemiLepEvtPartons::LightQBar];
    // calculate permutation
    // 0: nothing wrong
    if((lepBIndex==lepBIndexGen)&&(hadBIndex==hadBIndexGen)&&
       (((lightQIndex==lightQIndexGen   )&&(lightQBarIndex==lightQBarIndexGen))||
	((lightQIndex==lightQBarIndexGen)&&(lightQBarIndex==lightQIndexGen   )))) assignment=0;
    else{
      // 1: b quarks switched
      if((lepBIndex==hadBIndexGen)&&(hadBIndex==lepBIndexGen)) assignment=1;
      // 2: leptonic b and light quark switched
      if(((lepBIndex==lightQIndexGen)||(lepBIndex==lightQBarIndexGen))&&
	 (((lightQIndex==lepBIndexGen)||(lightQBarIndex==lepBIndexGen)))) assignment=2;
      // 3: hadronic b and light quark switched/
      if(((hadBIndex==lightQIndexGen)||(hadBIndex==lightQBarIndexGen))&&
	 (((lightQIndex==hadBIndexGen)||(lightQBarIndex==hadBIndexGen)))) assignment=3;
      // 4: light quark->leptonic b & leptonic b->hadronic b & hadronic b-> light quark
      if(((lepBIndex==lightQIndexGen)||(lepBIndex==lightQBarIndexGen))&&(hadBIndex==lepBIndexGen)&&
	 ((lightQIndex==hadBIndexGen)||(lightQBarIndex==hadBIndexGen))) assignment=4;
      // 5: light quark->hadronic b & hadronic b->leptonic b & leptonic b-> light quark
      if(((hadBIndex==lightQIndexGen)||(hadBIndex==lightQBarIndexGen))&&(lepBIndex==hadBIndexGen)&&
	 ((lightQIndex==lepBIndexGen)||(lightQBarIndex==lepBIndexGen))) assignment=5;
      // 6: hadronic/leptonic b-> light quarks &  light quarks->hadronic/leptonic b
      if(((hadBIndex     ==lightQIndexGen)||(hadBIndex     ==lightQBarIndexGen))&&
	 ((lepBIndex     ==lightQIndexGen)||(lepBIndex     ==lightQBarIndexGen))&&
	 ((lightQIndex   ==lepBIndexGen  )||(lightQIndex   ==hadBIndexGen     ))&&
	 ((lightQBarIndex==lepBIndexGen  )||(lightQBarIndex==hadBIndexGen     ))) assignment=6;
      // make sure that no relevant jet is missing
      std::vector<int> genJets_, recoJets_;
      // list of genJets
      genJets_ .push_back(lepBIndexGen     );
      genJets_ .push_back(hadBIndexGen     );
      genJets_ .push_back(lightQIndexGen   );
      genJets_ .push_back(lightQBarIndexGen);
      std::sort( genJets_.begin(), genJets_.end());
      // list of recoJets
      recoJets_.push_back(lepBIndex);
      recoJets_.push_back(hadBIndex);
      recoJets_.push_back(lightQIndex);
      recoJets_.push_back(lightQBarIndex);
      std::sort( recoJets_.begin(), recoJets_.end());
      // compare recoJets and genJets
      for(unsigned int i=0; i<recoJets_.size(); ++i){
	if( recoJets_[i]!=genJets_[i] ){ 
	  if(maxNJets<4){
	    std::cout << "ERROR: number of conidered jets can not be smaller than 4" << std::endl;
	    exit(1);
	  }

	  // 7: jet is missing
	  if( genJets_.back()>maxNJets-1 ) assignment=7;
	  // 8: wrong jet chosen (only valid if kinFitTtSemiLepEventHypothesis.maxNJets>4)
	  // e.g. took the wrong 4 out of 5 jets 
	  else{
	    if(assignment<0) assignment=8;
	  }
	  break;
	}
      }
    }
  }
  return assignment;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( MixedObjectsAnalyzer );
