#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "AnalysisDataFormats/TopObjects/interface/TtSemiLeptonicEvent.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiLepEvtPartons.h"

#include "TopAnalysis/TopAnalyzer/plugins/SemiLepLeptonAnalyzer.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/View.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// default constructor
SemiLepLeptonAnalyzer::SemiLepLeptonAnalyzer(const edm::ParameterSet& cfg):
  semiLepEvt_(cfg.getParameter<edm::InputTag>("semiLepEvent")),
  hypoKey_   (cfg.getParameter<std::string>  ("hypoKey"     )),
  recLeptons_(cfg.getParameter<edm::InputTag>("recLeptons"  )),
  genLeptons_(cfg.getParameter<edm::InputTag>("genLeptons"  )),
  verbose    (cfg.getParameter<int>          ("output"      )),
  weight_    (cfg.getParameter<edm::InputTag>("weight"      )),
  genPlots_  (cfg.getParameter<bool>         ("genPlots"    )),
  ingenPS_   (cfg.getParameter<edm::InputTag>("ingenPS"     )),
  recPlots_  (cfg.getParameter<bool>         ("recPlots"    )),
  preLep     (cfg.getParameter<bool>         ("useRecLeptonKinematicsBeforeFit")),
  useTree_   (cfg.getParameter<bool>         ("useTree"     )),
  valueLepPtRec (-1000),
  valueLepPtGen (-1000),
  valueLepEtaRec(-1000),
  valueLepEtaGen(-1000),
  valueLepYRec  (-1000),
  valueLepYGen  (-1000)
{
}

void
SemiLepLeptonAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup)
{

  // ---
  //     Z: Event information
  // ---
  if(useTree_){
    // event identification
    edm::EventAuxiliary aux = event.eventAuxiliary();
    runNumber             = aux.run();
    luminosityBlockNumber = aux.luminosityBlock();
    eventNumber           = aux.event();
    // PDF related info
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
  }

  // ---
  //     A: configure input
  // ---
  if(verbose>0) std::cout << "prepare inputs" << std::endl;
  if(verbose>1) std::cout << "a) event weight" << std::endl;
   // a) prepare the event weight
  weight = 1;
  // get weight from the CMSSW event
  edm::Handle<double> wgt;
  event.getByLabel(weight_, wgt);
  // ignore non existing weights
  if(wgt.isValid()) weight=*wgt;
  else{
    if(verbose>0) std::cout << "eventWeight not found" << std::endl;
    edm::LogInfo("weightNotFound") << "eventWeight not found";
  }
  if(verbose>1) std::cout << "weight=" << weight << std::endl;
  // b) get a handle for the TtSemiLepEvent 
  if(verbose>1) std::cout << "b) TtSemiLepEvent" << std::endl;
  edm::Handle<TtSemiLeptonicEvent> semiLepEvt;
  if(recPlots_) event.getByLabel(semiLepEvt_, semiLepEvt);
  // c) get a handle for genLepton collection
  if(verbose>1) std::cout << "c) genLeptonCollection" << std::endl;
  // get lepton collection
  edm::Handle<std::vector<reco::GenParticle> > genLeptons;
  if(genPlots_) event.getByLabel(genLeptons_, genLeptons);
  // d) get visible phase space information
  edm::Handle<bool > visPS;
  if(genPlots_) event.getByLabel(ingenPS_, visPS);
  inVisPS=!genPlots_||!*visPS ? false: true;

  // ---
  //     B: rec level plots
  // ---
  if(verbose>0) std::cout << std::endl << "create rec level plots" << std::endl;
  edm::Handle<edm::View<reco::Candidate> > recLeptons;
  if(recPlots_) event.getByLabel(recLeptons_, recLeptons);
  // get rec lepton (from kinfit hypothesis or before kinFit)
  const reco::Candidate* zero=0;
  if(!recPlots_) preLep=false;
  const reco::Candidate* recLepton = !recPlots_ ? zero : (preLep ? &(recLeptons->front()) : getLepton(semiLepEvt, hypoKey_));
  // fill rec histograms
  valueLepPtRec =-1000;
  valueLepEtaRec=-1000;
  valueLepYRec  =-1000; 
  if(recLepton){
    valueLepPtRec =recLepton->pt();
    valueLepEtaRec=recLepton->eta();
    valueLepYRec  =recLepton->rapidity(); 
    if(verbose>1) std::cout << "do filling" << std::endl;
    lepPtRec ->Fill( recLepton->pt()         , weight);
    lepEtaRec->Fill( recLepton->eta()        , weight);
    lepYRec  ->Fill( recLepton->rapidity()   , weight);
  }
  else if(verbose>1) std::cout << "no filling done" << std::endl;
 
  // ---
  //     C: gen level plots
  // ---
  if(verbose>0) std::cout << std::endl << "create gen level plots" << std::endl;
  // get generator level  jet collection using indices from B Hadron identification procedure

  //edm::Handle<std::vector<reco::GenParticle> > genLeptons;
  //if(genPlots_) event.getByLabel(genLeptons_, genLeptons);
  const reco::GenParticle* zeroL=0;
  const reco::GenParticle* genLepton = genPlots_&&inVisPS ? &genLeptons->at(0) : zeroL;

  // fill gen histograms
  valueLepPtGen =-1000;
  valueLepEtaGen=-1000;
  valueLepYGen  =-1000;
  if(genLepton){
    valueLepPtGen =genLepton->pt();
    valueLepEtaGen=genLepton->eta();
    valueLepYGen  =genLepton->rapidity();
    if(verbose>1) std::cout << "do filling" << std::endl;
    lepPtGen ->Fill( genLepton->pt()         , weight);
    lepEtaGen->Fill( genLepton->eta()        , weight);
    lepYGen  ->Fill( genLepton->rapidity()   , weight);
  }
  else if(verbose>1) std::cout << "no filling done" << std::endl;

  // ---
  //     D: gen rec correlation plots
  // ---  
  if(verbose>0) std::cout << std::endl << "correlation plots" << std::endl;
  if(recLepton){
    if(verbose>1) std::cout << "do filling" << std::endl;
    // fill correlation histograms
    lepPt_ ->Fill( valueLepPtGen , recLepton->pt()      , weight);
    lepEta_->Fill( valueLepEtaGen, recLepton->eta()     , weight);
    lepY_  ->Fill( valueLepYGen  , recLepton->rapidity(), weight);
  }
  else if(verbose>1) std::cout << "no filling done" << std::endl;
 

  if(useTree_) tree->Fill();
}

void 
SemiLepLeptonAnalyzer::beginJob()
{
  // check file service 
  edm::Service<TFileService> fs;
  if( !fs ) throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

  // print information
  if(verbose>0){ 
    std::cout << "running SemiLepLeptonAnalyzer" << std::endl;
    std::cout << "produce gen plots? " << genPlots_ << std::endl;
    std::cout << "produce rec plots? " << recPlots_ << std::endl;
  }
  
  //---
  //   book histograms
  //---
  // pt
  if(recPlots_) lepPtRec  = fs->make<TH1F>("lepPtRec" , "p_{t}^{lep} (rec) [GeV]", 1200, 0., 1200.);
  if(genPlots_) lepPtGen  = fs->make<TH1F>("lepPtGen" , "p_{t}^{lep} (gen) [GeV]", 1200, 0., 1200.);
  // eta
  if(recPlots_) lepEtaRec = fs->make<TH1F>("lepEtaRec", "#eta^{lep} (rec)"       , 100, -5., 5.);
  if(genPlots_) lepEtaGen = fs->make<TH1F>("lepEtaGen", "#eta^{lep} (gen)"       , 100, -5., 5.);
  // rapidity
  if(recPlots_) lepYRec   = fs->make<TH1F>("lepYRec"  , "y^{lep} (rec)"          , 100, -5., 5.);
  if(genPlots_) lepYGen   = fs->make<TH1F>("lepYGen"  , "y^{lep} (gen)"          , 100, -5., 5.);
  // correlation
  if(recPlots_&&genPlots_){
    lepPt_  = fs->make<TH2F>("lepPt_" , "p_{t}^{lep} (gen vs rec) [GeV]", 1200,  0., 1200., 1200,  0., 1200.);
    lepEta_ = fs->make<TH2F>("lepEta_", "#eta^{lep} (gen vs rec)"       ,  100, -5.,    5.,  100, -5.,    5.);
    lepY_   = fs->make<TH2F>("lepY_"  , "y^{lep} (gen vs rec)"          ,  100, -5.,    5.,  100, -5.,    5.);
  }
  //---
  //   create branches
  //---
  if(useTree_){
    // define the tree and make it known to the TFileService
    tree = fs->make<TTree>("tree","tree",0);
    // weight
    tree->Branch("weight" , &weight , "weight/F" );
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
    // variables
    tree->Branch("lepPtRec" , &valueLepPtRec , "lepPtRec/F" );
    tree->Branch("lepPtGen" , &valueLepPtGen , "lepPtGen/F" );
    tree->Branch("lepEtaRec", &valueLepEtaRec, "lepEtaRec/F");
    tree->Branch("lepEtaGen", &valueLepEtaGen, "lepEtaGen/F");
    tree->Branch("lepYRec"  , &valueLepYRec  , "lepYRec/F"  );
    tree->Branch("lepYGen"  , &valueLepYGen  , "lepYGen/F"  );
    // within visible phase space?
    tree->Branch("inVisPS"       , &inVisPS       , "inVisPS/O"      );
  }
}

void
SemiLepLeptonAnalyzer::endJob() 
{
}

const reco::Candidate*
SemiLepLeptonAnalyzer::getLepton(const edm::Handle<TtSemiLeptonicEvent> semiLepEvt, const std::string hypoKey_){
  const reco::Candidate* zero=0;
  if(verbose>1) std::cout << "get lepton candidate from kinfit" << std::endl;
  // check if hypothesis is available and valid in this event
  if(!semiLepEvt->isHypoValid(hypoKey_)){
    if(verbose>0) std::cout << "Hypothesis not valid for this event" << std::endl;
  }
  else{
    // check if decay is semileptonic
    if( !semiLepEvt->genEvent().isAvailable()                  || 
	(semiLepEvt->genEvent()->isSemiLeptonic(WDecay::kMuon) || 
	 semiLepEvt->genEvent()->isSemiLeptonic(WDecay::kElec) ||
	 semiLepEvt->genEvent()->isSemiLeptonic(WDecay::kTau )   ) ){
      if(verbose>1) std::cout << "semilep. decay or no ttgenevent available" << std::endl;
      // get reconstructed lepton from the hypothesis 
      const reco::Candidate* lepton  = semiLepEvt->singleLepton(hypoKey_);
      // return lepton
      if(verbose>0) std::cout << "rec lepton candidate found!" << std::endl;
      return lepton;
    }
  }
  // otherwise: return NULL
  if(verbose>0) std::cout << "no rec lepton found" << std::endl;
  return zero;
}
