#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "AnalysisDataFormats/TopObjects/interface/TtSemiLeptonicEvent.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiLepEvtPartons.h"

#include "TopAnalysis/TopAnalyzer/plugins/SemiLepBjetAnalyzer.h"
#include "DataFormats/Math/interface/deltaR.h"

// default constructor
SemiLepBjetAnalyzer::SemiLepBjetAnalyzer(const edm::ParameterSet& cfg):
  semiLepEvt_(cfg.getParameter<edm::InputTag>("semiLepEvent")),
  hypoKey_   (cfg.getParameter<std::string>  ("hypoKey"     )),
  genJets_   (cfg.getParameter<edm::InputTag>("genJets"     )),
  genLeptons_(cfg.getParameter<edm::InputTag>("genLeptons"  )),

  bJetCollection_(cfg.getParameter<bool>     ("bJetCollection")),
  recoJets_  (cfg.getParameter<edm::InputTag>("recoJets"    )),
  verbose    (cfg.getParameter<int>          ("output"      )),
  weight_    (cfg.getParameter<edm::InputTag>("weight"      )),
  genPlots_  (cfg.getParameter<bool>         ("genPlots"    )),
  ingenPS_   (cfg.getParameter<edm::InputTag>("ingenPS"     )),
  recPlots_  (cfg.getParameter<bool>         ("recPlots"    )),
  preBjets   (cfg.getParameter<bool>         ("useRecBjetsKinematicsBeforeFit")),
  bHadJetIdx_    (cfg.getParameter<edm::InputTag> ("BHadJetIndex"    )),
  antibHadJetIdx_(cfg.getParameter<edm::InputTag> ("AntiBHadJetIndex")),
  useClosestDrBs_(cfg.getParameter<bool>   ("useClosestDrBs")),
  useTree_   (cfg.getParameter<bool>("useTree")),
  valueBqPtRec(-1000),
  valueBqPtGen(-1000),
  valueLeadBqPtRec(-1000),
  valueLeadBqPtGen(-1000),
  valueSubLeadBqPtRec(-1000),
  valueSubLeadBqPtGen(-1000),
  valueBqEtaRec(-1000),
  valueBqEtaGen(-1000),
  valueBqYRec(-1000),
  valueBqYGen(-1000),
  valueBbarqPtRec(-1000),
  valueBbarqPtGen(-1000),
  valueBbarqEtaRec(-1000),
  valueBbarqEtaGen(-1000),
  valueBbarqYRec(-1000),
  valueBbarqYGen(-1000),
  valueBbbarPtRec(-1000),
  valueBbbarPtGen(-1000),
  valueBbbarEtaRec(-1000),
  valueBbbarEtaGen(-1000),
  valueBbbarYRec(-1000),
  valueBbbarYGen(-1000),
  valueBbbarMassRec(-1000),
  valueBbbarMassGen(-1000),
  valueLbMassRec(-1000),
  valueLbMassGen(-1000),
  //valuexBHadRec(-1000),
  //valuexBHadGen(-1000),
  //valuexBLepRec(-1000),
  //valuexBLepGen(-1000),
  bbSwapBetter(false),
  valueAssignment(-1000)
{
}

void
SemiLepBjetAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup)
{

  // ---
  //     Z: Event information
  // ---
  if(useTree_){
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

  if(verbose>0){
    std::cout << "processing SemiLepBjetAnalyzer" << std::endl;
    std::cout << "- prepare inputs" << std::endl;
  }
  if(verbose>1){
    std::cout << "  genPlots_: " << genPlots_ << std::endl;
    std::cout << "  recPlots_: " << recPlots_ << std::endl;
    std::cout << "  preBjets : " << preBjets  << std::endl;
  }
  if(verbose>1) std::cout << "a) event weight: ";
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
  if(verbose>1) std::cout << weight << std::endl;
  // b) get a handle for the TtSemiLepEvent 
  if(verbose>1) std::cout << "b) TtSemiLepEvent" << std::endl;
  edm::Handle<TtSemiLeptonicEvent> semiLepEvt;
  if(recPlots_) event.getByLabel(semiLepEvt_, semiLepEvt);
  // c1) get a handle for genjet collection
  if(verbose>1) std::cout << "c1) genJetCollection" << std::endl;
  // get jet collection
  edm::Handle<reco::GenJetCollection> genJets;
  if(genPlots_) event.getByLabel(genJets_, genJets);
  // c2) get a handle for reco jet collection
  if(verbose>1) std::cout << "c2) recoJetCollection" << std::endl;
  if(!recPlots_) preBjets=false;
  edm::Handle<pat::JetCollection> recoJets;
  if(recPlots_) event.getByLabel(recoJets_, recoJets);
  // d) position of the jet within genjet collection identified as b-hadron jet
  if(verbose>1) std::cout << "d) b-jet indices produced by TopAnalysis/TopUtils/GenLevelBJetProducer" << std::endl;
  int bIX   =-1;
  int bbarIX=-1;
  // e) get lepton collection
  edm::Handle<std::vector<reco::GenParticle> > genLeptons;
  if(genPlots_) event.getByLabel(genLeptons_, genLeptons);
  if(verbose>1) std::cout << "e) gen lepton collection" << std::endl;
  // f) visible hadron level gen phase space indicator
  edm::Handle<bool > visPS;
  if(genPlots_) event.getByLabel(ingenPS_, visPS);
  inVisPS=false;
  if(!visPS.failedToGet()&&visPS.isValid()&&*visPS) inVisPS=true;
  if(verbose>1) std::cout << "f) visible phase space indicator: " << inVisPS << std::endl;
  if(genPlots_){
    // take leading two jets assuming that collection genJets consists of b-jets only in the order 1st entry: b-jet, 2nd entry: anti-bjet
    if(bJetCollection_==true){
      bIX   =0;
      bbarIX=1;
    }
    else{
      // get indices produced by TopAnalysis/TopUtils/GenLevelBJetProducer.cc
      std::vector<int> VBHadJetIdx, VAntiBHadJetIdx;
      edm::Handle<std::vector<int> > BHadJetIndex, AntiBHadJetIndex;
      event.getByLabel(bHadJetIdx_    , BHadJetIndex    );
      event.getByLabel(antibHadJetIdx_, AntiBHadJetIndex);
      // save them in vector
      for (int i=0; i<(int) BHadJetIndex->size()    ; ++i){  VBHadJetIdx    .push_back(BHadJetIndex->at(i)    );};
      for (int i=0; i<(int) AntiBHadJetIndex->size(); ++i){  VAntiBHadJetIdx.push_back(AntiBHadJetIndex->at(i));};
      // use leading vector entries finally 
      if(VBHadJetIdx    .size()>0) bIX   =VBHadJetIdx    [0];
      if(VAntiBHadJetIdx.size()>0) bbarIX=VAntiBHadJetIdx[0];
      // output for user information
      if(verbose>1){
	std::cout << "N(b jets): "    << VBHadJetIdx    .size() << std::endl;
	for(unsigned int idx=0; idx<VBHadJetIdx    .size(); ++idx){
	  if(idx==0) std::cout << "("; 
	  std::cout << VBHadJetIdx[idx];
	  if(idx<VBHadJetIdx.size()-1) std::cout << ", ";
	  else std::cout << ")" << std::endl;
	}
	std::cout << "N(bbar jets): " << VAntiBHadJetIdx.size() << std::endl;
	for(unsigned int idx=0; idx<VAntiBHadJetIdx.size(); ++idx){
	  if(idx==0) std::cout << "("; 
	  std::cout << VAntiBHadJetIdx[idx];
	  if(idx<VAntiBHadJetIdx.size()-1) std::cout << ", ";
	  else std::cout << ")" << std::endl;
	}
      }
    }
  }
  // ---
  //     B: rec level plots
  // ---
  if(verbose>0) std::cout << std::endl << "create rec level plots" << std::endl;
  // get rec b and bbar (from kinfit hypothesis)
  const reco::Candidate* zero=0;
  std::pair<const reco::Candidate*,const reco::Candidate*> recBJets = recPlots_ ? getbJets(semiLepEvt, hypoKey_) : std::make_pair(zero,zero);
  // get rec b and bbar (assignment from kinfit hypothesis but no shifted kinematics)
  // (not used for xB)
  if(preBjets) recBJets = recPlots_ ? getPreBJets(semiLepEvt, hypoKey_, recoJets) : std::make_pair(zero,zero);
  const reco::Candidate* b   = recBJets.first ;
  const reco::Candidate* bbar= recBJets.second;

  //double recxBhad = -1; 
  //double recxBlep = -1;

  valueBqPtRec        =-1000;
  valueLeadBqPtRec    =-1000;
  valueSubLeadBqPtRec =-1000;
  valueBqEtaRec       =-1000;
  valueBqYRec         =-1000;
  valueBbarqPtRec     =-1000;
  valueBbarqEtaRec    =-1000;
  valueBbarqYRec      =-1000;
  valueBbbarPtRec     =-1000;
  valueBbbarEtaRec    =-1000;
  valueBbbarYRec      =-1000;
  valueBbbarMassRec   =-1000;
  valueLbMassRec      =-1000;
  //valuexBHadRec    =recxBhad;
  //valuexBLepRec    =recxBlep;
  
  // fill rec histograms
  if(b&&bbar){
    // reconstructed top and W
    
    //reco::Particle::LorentzVector hadWRec=semiLepEvt->hadronicDecayW(hypoKey_)->p4();
    //reco::Particle::LorentzVector lepWRec=semiLepEvt->leptonicDecayW(hypoKey_)->p4();
    //reco::Particle::LorentzVector hadTopRec=semiLepEvt->hadronicDecayTop(hypoKey_)->p4();
    //reco::Particle::LorentzVector lepTopRec=semiLepEvt->leptonicDecayTop(hypoKey_)->p4();
    
    // pick b(lepton/hadron) wrt lepton charge
    double lepChargeRec= semiLepEvt->singleLepton(hypoKey_)->charge();
    // q(lep) < 0 -> q(leptonic bjet)>0 -> bjet(leptonic)=bbar
    reco::Particle::LorentzVector lepBRec = lepChargeRec < 0 ? bbar->p4() : b   ->p4();
    //reco::Particle::LorentzVector hadBRec = lepChargeRec < 0 ? b   ->p4() : bbar->p4();

    // fill tree variables   
    valueBqPtRec        =b->pt();
    valueLeadBqPtRec    =b->pt()>bbar->pt() ? b->pt() : bbar->pt();
    valueSubLeadBqPtRec =b->pt()<bbar->pt() ? b->pt() : bbar->pt();
    valueBqEtaRec       =b->eta();
    valueBqYRec         =b->rapidity(); 
    valueBbarqPtRec     =bbar->pt();
    valueBbarqEtaRec    =bbar->eta();
    valueBbarqYRec      =bbar->rapidity();
    valueBbbarPtRec     =(b->p4()+bbar->p4()).pt();
    valueBbbarEtaRec    =(b->p4()+bbar->p4()).eta();
    valueBbbarYRec      =(b->p4()+bbar->p4()).Rapidity();
    valueBbbarMassRec   =(b->p4()+bbar->p4()).mass();
    //valuexBHadRec    =recxBhad;
    //valuexBLepRec    =recxBlep;

    // calculate more complex variables
    valueLbMassRec=(semiLepEvt->singleLepton(hypoKey_)->p4()+lepBRec).mass();
    //recxBhad=(1/(1-(hadWRec.mass2())/(hadTopRec.mass2()))) * 2*(hadBRec.Dot(hadTopRec)/(hadTopRec.mass2()));
    //recxBlep=(1/(1-(lepWRec.mass2())/(lepTopRec.mass2()))) * 2*(lepBRec.Dot(lepTopRec)/(lepTopRec.mass2()));

    // debug output
    if(verbose>1) std::cout << "do filling" << std::endl;
    // fill plots
    bqPtRec        ->Fill( valueBqPtRec       , weight );
    bqPtRec        ->Fill( valueBbarqPtRec    , weight );
    bqPtLeadRec    ->Fill( valueLeadBqPtRec   , weight);
    bqPtSubLeadRec ->Fill( valueSubLeadBqPtRec, weight);
    bqEtaRec       ->Fill( valueBqEtaRec      , weight);
    bqEtaRec       ->Fill( valueBbarqEtaRec   , weight);
    bqYRec         ->Fill( valueBqYRec        , weight);
    bqYRec         ->Fill( valueBbarqYRec     , weight);
    bbbarPtRec     ->Fill( valueBbbarPtRec    , weight);
    bbbarEtaRec    ->Fill( valueBbbarEtaRec   , weight);
    bbbarYRec      ->Fill( valueBbbarYRec     , weight);
    bbbarMassRec   ->Fill( valueBbbarMassRec  , weight);
    lbMassRec      ->Fill( valueLbMassRec     , weight);
    //xBLepRec->Fill(recxBlep, weight);
    //xBHadRec->Fill(recxBhad, weight);
    //xBRec   ->Fill(recxBlep, weight);
    //xBRec   ->Fill(recxBhad, weight);
  }
  else if(verbose>1) std::cout << "no filling done" << std::endl;
 
  // ---
  //     C: gen level plots
  // ---
  if(verbose>0) std::cout << std::endl << "create gen level plots" << std::endl;
  // ensure there are indices for both
  if(bIX<-1||bbarIX<-1){
    if(genPlots_&&verbose>0) std::cout << "WARNING: gen bjet indices incomplete!" << std::endl;
    genPlots_=false;
  }
  // get generator level b-Jets from gen jet collection using indices from B Hadron identification procedure
  const reco::GenJet*      genb      = genPlots_&&inVisPS ? getJetFromCollection(*genJets,bIX   ) : 0;
  const reco::GenJet*      genbbar   = genPlots_&&inVisPS ? getJetFromCollection(*genJets,bbarIX) : 0;
  const reco::GenParticle* genLepton = genPlots_&&inVisPS ? &genLeptons->at(0) : 0;

  //double genxBhad=-1;
  //double genxBlep=-1;
  // fill gen histograms
  valueBqPtGen       =-1000.;
  valueLeadBqPtGen   =-1000.;
  valueSubLeadBqPtGen=-1000.;
  valueBqEtaGen      =-1000.;
  valueBqYGen        =-1000.;
  valueBbarqPtGen    =-1000.;
  valueBbarqEtaGen   =-1000.;
  valueBbarqYGen     =-1000.;
  valueBbbarPtGen    =-1000.;
  valueBbbarEtaGen   =-1000.;
  valueBbbarYGen     =-1000.;
  valueBbbarMassGen  =-1000.;
  valueLbMassGen     =-1000.;
  if(genb&&genbbar){
    // top and W from parton truth
    
    //reco::Particle::LorentzVector hadWGen=semiLepEvt->hadronicDecayW()->p4();
    //reco::Particle::LorentzVector lepWGen=semiLepEvt->leptonicDecayW()->p4();
    //reco::Particle::LorentzVector hadTopGen=semiLepEvt->hadronicDecayTop()->p4();
    //reco::Particle::LorentzVector lepTopGen=semiLepEvt->leptonicDecayTop()->p4();
    
    // pick b(lepton/hadron) wrt lepton charge
    double lepChargeGen=genLepton->charge();
    // q(lep) < 0 -> bjet(leptonic)>0 -> bjet(leptonic)=bbar
    reco::Particle::LorentzVector lepBGen = lepChargeGen < 0 ? genbbar->p4() : genb   ->p4();
    //reco::Particle::LorentzVector hadBGen = lepChargeGen < 0 ? genb   ->p4() : genbbar->p4();

    // set tree gen variables
    valueBqPtGen =genb->pt();
    valueLeadBqPtGen    = genb->pt()>genbbar->pt() ? genb->pt() : genbbar->pt();
    valueSubLeadBqPtGen = genb->pt()<genbbar->pt() ? genb->pt() : genbbar->pt();
    valueBqEtaGen   =genb->eta();
    valueBqYGen     =genb->rapidity();
    valueBbarqPtGen =genbbar->pt();
    valueBbarqEtaGen=genbbar->eta();
    valueBbarqYGen  =genbbar->rapidity();
    valueBbbarPtGen  =(genb->p4()+genbbar->p4()).pt();
    valueBbbarEtaGen =(genb->p4()+genbbar->p4()).eta();
    valueBbbarYGen   =(genb->p4()+genbbar->p4()).Rapidity();
    valueBbbarMassGen=(genb->p4()+genbbar->p4()).mass();
    //valuexBHadGen=genxBhad;
    //valuexBLepGen=genxBlep;

    // calculate more complex variables
    valueLbMassGen = (genLepton->p4()+lepBGen).mass();
    //genxBhad=(1/(1-(hadWGen.mass2())/(hadTopGen.mass2()))) * 2*(hadBGen.Dot(hadTopGen)/(hadTopGen.mass2()));
    //genxBlep=(1/(1-(lepWGen.mass2())/(lepTopGen.mass2()))) * 2*(lepBGen.Dot(lepTopGen)/(lepTopGen.mass2()));

    // debug output
    if(!recPlots_&&valueBqPtGen>0&&valueBbarqPtGen>0&&verbose>1){ 
      std::cout << "gen truth filling" << std::endl;
      std::cout << "b(pt, eta)=("    << valueBqPtGen    << ", " << valueBqEtaGen    << ")" << std::endl;
      std::cout << "bbar(pt, eta)=(" << valueBbarqPtGen << ", " << valueBbarqEtaGen << ")" << std::endl;
    }
    if(verbose>1) std::cout << "do filling" << std::endl;
    // fill plots
    bqPtGen        ->Fill( valueBqPtGen       , weight );
    bqPtGen        ->Fill( valueBbarqPtGen    , weight );
    bqPtLeadGen    ->Fill( valueLeadBqPtGen   , weight );
    bqPtSubLeadGen ->Fill( valueSubLeadBqPtGen, weight );
    bqEtaGen       ->Fill( valueBqEtaGen      , weight );
    bqEtaGen       ->Fill( valueBbarqEtaGen   , weight );
    bqYGen         ->Fill( valueBqYGen        , weight );
    bqYGen         ->Fill( valueBbarqYGen     , weight );
    bbbarPtGen     ->Fill( valueBbbarPtGen    , weight );
    bbbarEtaGen    ->Fill( valueBbarqEtaGen   , weight );
    bbbarYGen      ->Fill( valueBbbarYGen     , weight );
    bbbarMassGen   ->Fill( valueBbbarMassGen  , weight );
    lbMassGen      ->Fill( valueLbMassGen     , weight );
    //xBLepGen->Fill(genxBlep, weight);
    //xBHadGen->Fill(genxBhad, weight);
    //xBGen   ->Fill(genxBlep, weight);
    //xBGen   ->Fill(genxBhad, weight);
  }
  else if(verbose>1) std::cout << "no filling done" << std::endl;

  // ---
  //     D: gen rec correlation plots
  // ---  
  if(verbose>0) std::cout << std::endl << "correlation plots" << std::endl;
  // check ONLY if RECO b-jets exist in the event
  bbSwapBetter = false;
  if(b&&bbar){
    // parton level quark vs reco jet assignment
    valueAssignment=checkPartonAssignment(semiLepEvt, 5); // FIXME: N(jets) in KinFit hardcoded for the moment
    if(verbose>1){
      std::cout << "mlb(gen correct)=" << valueLbMassGen  << std::endl;
      std::cout << "mlb(rec KinFit )=" << valueLbMassRec  << std::endl;
      std::cout << "reco jet assignment (kinFit) wrt reco jet-quark matching: " << valueAssignment << std::endl;
    }
    if(verbose>1) std::cout << "do filling" << std::endl;
    // fill lead b-jet and bbbar correlation histograms 
    // (those have only one entry per event and are independend of the rec-gen b-jet association)
    if(verbose>1){
      std::cout << "valueLeadBqPtGen: " << valueLeadBqPtGen << std::endl; 
      std::cout << "valueLeadBqPtRec: " << valueLeadBqPtRec << std::endl;
      std::cout << "valueSubLeadBqPtGen: " << valueSubLeadBqPtGen << std::endl; 
      std::cout << "valueSubLeadBqPtRec: " << valueSubLeadBqPtRec << std::endl;
      std::cout << "valueBbbarPtGen: " << valueBbbarPtGen << std::endl; 
      std::cout << "valueBbarqEtaGen: " << valueBbbarPtRec << std::endl;
      std::cout << "valueBbarqEtaGen: " << valueBbarqEtaGen << std::endl; 
      std::cout << "valueBbbarEtaRec: " << valueBbbarEtaRec << std::endl;
      std::cout << "valueBbbarYGen: " << valueBbbarYGen << std::endl; 
      std::cout << "valueBbbarYRec: " << valueBbbarYRec << std::endl;
      std::cout << "valueBbbarMassGen: " << valueBbbarMassGen   << std::endl; 
      std::cout << "valueBbbarMassRec: " << valueBbbarMassRec   << std::endl;
    }
    bqPtLead_   ->Fill(valueLeadBqPtGen   , valueLeadBqPtRec   , weight);
    bqPtSubLead_->Fill(valueSubLeadBqPtGen, valueSubLeadBqPtRec, weight);
    bbbarPt_    ->Fill(valueBbbarPtGen    , valueBbbarPtRec    , weight);
    bbbarEta_   ->Fill(valueBbarqEtaGen   , valueBbbarEtaRec   , weight);
    bbbarY_     ->Fill(valueBbbarYGen     , valueBbbarYRec     , weight);
    bbbarMass_  ->Fill(valueBbbarMassGen  , valueBbbarMassRec  , weight);
    if(valueLbMassRec>=0) lbMass_->Fill( valueLbMassGen, valueLbMassRec, weight );
    //if(genxBlep>=0&&recxBlep>=0) xB_->Fill( genxBlep, recxBlep, weight );
    //if(genxBhad>=0&&recxBhad>=0) xB_->Fill( genxBhad, recxBhad, weight );
    
    // for plots with combined b-jet and bbar-jet measurement:
    // take into account that both values end up in the same histogram
    // -> search for best gen-reco association based on minimal dR 

    // start with gen: association from B/AntiB, reco: association from KinFit
    double bpt    =valueBqPtRec;         
    double bbarpt =valueBbarqPtRec; 
    double beta   =valueBqEtaRec;        
    double bbareta=valueBbarqEtaRec; 
    double bY     =valueBqYRec;   
    double bbarY  =valueBbarqYRec;
    bool swapped = false;
    // check if dR is smaller for swopping gen-rec assignment
    if(genb&&genbbar){
      double valueBqPhiRec   =b   ->phi();
      double valueBbarqPhiRec=bbar->phi();
      double valueBqPhiGen   =genb   ->phi();
      double valueBbarqPhiGen=genbbar->phi();
      double dROri =deltaR( valueBqEtaRec   , valueBqPhiRec   , valueBqEtaGen   , valueBqPhiGen   )+
	            deltaR( valueBbarqEtaRec, valueBbarqPhiRec, valueBbarqEtaGen, valueBbarqPhiGen);
      double dRSwap=deltaR( valueBbarqEtaRec, valueBbarqPhiRec, valueBqEtaGen   , valueBqPhiGen   )+ 
	deltaR( valueBqEtaRec   , valueBqPhiRec   , valueBbarqEtaGen, valueBbarqPhiGen);
      if(dROri > dRSwap){
	bbSwapBetter = true;
	// use the better dR combination only if parameter useClosestDrBs_ is true
	if(useClosestDrBs_){
	  swapped = true;
	  bY     =valueBbarqYRec; 
	  bbarY  =valueBqYRec;
	  bpt    =valueBbarqPtRec; 
	  bbarpt =valueBqPtRec;
	  beta   =valueBbarqEtaRec; 
	  bbareta=valueBqEtaRec;
	}
	if(verbose>1){
	  std::cout << "dR unswapped = " << dROri  << std::endl;
	  std::cout << "dR swapped   = " << dRSwap << std::endl;
	}
      }
    }
    // fill combined b-jet and bbar-jet correlation histograms 
    // (those have both b-jets as entry per event and depend on rec-gen b-jet association)
    bqPt_ ->Fill( valueBqPtGen    , bpt    , weight);
    bqPt_ ->Fill( valueBbarqPtGen , bbarpt , weight);
    bqEta_->Fill( valueBqEtaGen   , beta   , weight);
    bqEta_->Fill( valueBbarqEtaGen, bbareta, weight);
    bqY_  ->Fill( valueBqYGen     , bY     , weight);
    bqY_  ->Fill( valueBbarqYGen  , bbarY  , weight);
    // output
    if(verbose>1){
      std::cout << "swapped: " << swapped << "; bbSwapBetter: " << bbSwapBetter << std::endl;
      std::cout << "pt: genb   = " << valueBqPtGen    << "; recb   = " << bpt    << "-- eta: genb   = " << valueBqEtaGen    << "; recb   = " << beta << std::endl;
      std::cout << "pt: genbbar= " << valueBbarqPtGen << "; recbbar= " << bbarpt << "-- eta: genbbar= " << valueBbarqEtaGen << "; recbbar= " << bbareta << std::endl;
    }
  }
  else if(verbose>1) std::cout << "no filling done" << std::endl;

  if(useTree_) tree->Fill();
}

void 
SemiLepBjetAnalyzer::beginJob()
{
  // check file service 
  edm::Service<TFileService> fs;
  if( !fs ) throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

  // print information
  if(verbose>0){ 
    std::cout << "running SemiLepBjetAnalyzer" << std::endl;
    std::cout << "produce gen plots? " << genPlots_ << std::endl;
    std::cout << "produce rec plots? " << recPlots_ << std::endl;
  }
  
  //---
  //   book histograms
  //---
  // 1D rec und gen
  // A) b-quark quantities
  // pt
  if(recPlots_) bqPtRec  = fs->make<TH1F>("bqPtRec" , "p_{t}^{b and #bar{b}} (rec) [GeV]", 1200, 0., 1200.);
  if(genPlots_) bqPtGen  = fs->make<TH1F>("bqPtGen" , "p_{t}^{b and #bar{b}} (gen) [GeV]", 1200, 0., 1200.);
  if(recPlots_) bqPtLeadRec  = fs->make<TH1F>("bqPtLeadRec" , "p_{t}^{lead b} (rec) [GeV]", 1200, 0., 1200.);
  if(genPlots_) bqPtLeadGen  = fs->make<TH1F>("bqPtLeadGen" , "p_{t}^{lead b} (gen) [GeV]", 1200, 0., 1200.);
  if(recPlots_) bqPtSubLeadRec  = fs->make<TH1F>("bqPtSubLeadRec" , "p_{t}^{sublead b} (rec) [GeV]", 1200, 0., 1200.);
  if(genPlots_) bqPtSubLeadGen  = fs->make<TH1F>("bqPtSubLeadGen" , "p_{t}^{sublead b} (gen) [GeV]", 1200, 0., 1200.);
  // eta
  if(recPlots_) bqEtaRec = fs->make<TH1F>("bqEtaRec", "#eta^{b and #bar{b}} (rec)"       , 100, -5., 5.);
  if(genPlots_) bqEtaGen = fs->make<TH1F>("bqEtaGen", "#eta^{b and #bar{b}} (gen)"       , 100, -5., 5.);
  // rapidity
  if(recPlots_) bqYRec   = fs->make<TH1F>("bqYRec"  , "y^{b and #bar{b}} (rec)"          , 100, -5., 5.);
  if(genPlots_) bqYGen   = fs->make<TH1F>("bqYGen"  , "y^{b and #bar{b}} (gen)"          , 100, -5., 5.);
  // B) bbbar-system quantities
  // pt
  if(recPlots_) bbbarPtRec  = fs->make<TH1F>("bbbarPtRec"   , "p_{t}^{b#bar{b}} (rec) [GeV]", 1200,  0., 1200.);
  if(genPlots_) bbbarPtGen  = fs->make<TH1F>("bbbarPtGen"   , "p_{t}^{b#bar{b}} (gen) [GeV]", 1200,  0., 1200.);
  // eta
  if(recPlots_) bbbarEtaRec = fs->make<TH1F>("bbbarEtaRec"  , "#eta^{b#bar{b}} (rec)"       ,  100, -5.,    5.);
  if(genPlots_) bbbarEtaGen = fs->make<TH1F>("bbbarEtaGen"  , "#eta^{b#bar{b}} (gen)"       ,  100, -5.,    5.);
  // rapidity
  if(recPlots_) bbbarYRec   = fs->make<TH1F>("bbbarYRec"    , "y^{b#bar{b}} (rec)"          ,  100, -5.,    5.);
  if(genPlots_) bbbarYGen   = fs->make<TH1F>("bbbarYGen"    , "y^{b#bar{b}} (gen)"          ,  100, -5.,    5.);
  // mass
  if(recPlots_) bbbarMassRec= fs->make<TH1F>("bbbarMassRec" , "m^{b#bar{b}} (rec) [GeV]"    , 1200,  0., 1200.);
  if(genPlots_) bbbarMassGen= fs->make<TH1F>("bbbarMassGen" , "m^{b#bar{b}} (gen) [GeV]"    , 1200,  0., 1200.);
  // C) lb-system quantities
  if(recPlots_) lbMassRec= fs->make<TH1F>("lbMassRec" , "m^{lb} (rec) [GeV]"    , 500,  0., 500.);
  if(genPlots_) lbMassGen= fs->make<TH1F>("lbMassGen" , "m^{lb} (gen) [GeV]"    , 500,  0., 500.);
  // D) hadronization variable
  //if(recPlots_) xBRec= fs->make<TH1F>("xBRec" , "x_B", 200, -0.5, 1.5 );
  //if(genPlots_) xBGen= fs->make<TH1F>("xBGen" , "x_B", 200, -0.5, 1.5 );
  //if(recPlots_) xBHadRec= fs->make<TH1F>("xBHadRec" , "x_B (hadronic branch)", 200, -0.5, 1.5 );
  //if(genPlots_) xBHadGen= fs->make<TH1F>("xBHadGen" , "x_B (hadronic branch)", 200, -0.5, 1.5 );
  //if(recPlots_) xBLepRec= fs->make<TH1F>("xBLepRec" , "x_B (leptonic branch)", 200, -0.5, 1.5 );
  //if(genPlots_) xBLepGen= fs->make<TH1F>("xBLepGen" , "x_B (leptonic branch)", 200, -0.5, 1.5 );

  // 2D correlation
  if(recPlots_&&genPlots_){
    // A) b-jet quantities
    bqPt_       = fs->make<TH2F>("bqPt_" , "p_{t}^{b and #bar{b}} (gen vs rec) [GeV]", 1200, 0., 1200., 1200, 0., 1200.);
    bqPtLead_   = fs->make<TH2F>("bqPtLead_" , "p_{t}^{lead b} (gen vs rec) [GeV]"   , 1200, 0., 1200., 1200, 0., 1200.);
    bqPtSubLead_= fs->make<TH2F>("bqPtSubLead_" , "p_{t}^{lead b} (gen vs rec) [GeV]", 1200, 0., 1200., 1200, 0., 1200.);
    bqEta_ = fs->make<TH2F>("bqEta_", "#eta^{b and #bar{b}} (gen vs rec)"       ,  100, -5.,    5.,  100, -5.,    5.);
    bqY_   = fs->make<TH2F>("bqY_"  , "y^{b and #bar{b}} (gen vs rec)"          ,  100, -5.,    5.,  100, -5.,    5.);
    // B) bbbar-system quantities
    bbbarPt_  = fs->make<TH2F>("bbbarPt_"   , "p_{t}^{b#bar{b}} (gen vs rec) [GeV]", 1200,  0., 1200., 1200,  0., 1200.);
    bbbarEta_ = fs->make<TH2F>("bbbarEta_"  , "#eta^{b#bar{b}} (gen vs rec)"       ,  100, -5.,    5.,  100, -5.,    5.);
    bbbarY_   = fs->make<TH2F>("bbbarY_"    , "y^{b#bar{b}} (gen vs rec)"          ,  100, -5.,    5.,  100, -5.,    5.);
    bbbarMass_= fs->make<TH2F>("bbbarMass_" , "m^{b#bar{b}} (gen vs rec) [GeV]"    , 1200,  0., 1200., 1200,  0., 1200.);
    // C) lb-system quantities
    lbMass_   = fs->make<TH2F>("lbMass_"    , "m^{lb} (gen vs rec) [GeV]"          ,  500,  0.,  500.,  500,  0.,  500.);
    // D) hadronization variable
    //xB_= fs->make<TH2F>("xB_" , "xB_", 200, -0.5, 1.5, 200, -0.5, 1.5 );
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
    // variables
    tree->Branch("bqPtRec" , &valueBqPtRec , "bqPtRec/F" );
    tree->Branch("bqPtGen" , &valueBqPtGen , "bqPtGen/F" );
    tree->Branch("bqEtaRec", &valueBqEtaRec, "bqEtaRec/F");
    tree->Branch("bqEtaGen", &valueBqEtaGen, "bqEtaGen/F");
    tree->Branch("bqYRec"  , &valueBqYRec  , "bqYRec/F"  );
    tree->Branch("bqYGen"  , &valueBqYGen  , "bqYGen/F"  );
    tree->Branch("bbarqPtRec" , &valueBbarqPtRec , "bbarqPtRec/F" );
    tree->Branch("bbarqPtGen" , &valueBbarqPtGen , "bbarqPtGen/F" );
    tree->Branch("bbarqEtaRec", &valueBbarqEtaRec, "bbarqEtaRec/F");
    tree->Branch("bbarqEtaGen", &valueBbarqEtaGen, "bbarqEtaGen/F");
    tree->Branch("bbarqYRec"  , &valueBbarqYRec  , "bbarqYRec/F"  );
    tree->Branch("bbarqYGen"  , &valueBbarqYGen  , "bbarqYGen/F"  );
    tree->Branch("bbbarMassRec", &valueBbbarMassRec, "bbbarMassRec/F");
    tree->Branch("bbbarMassGen", &valueBbbarMassGen, "bbbarMassGen/F");
    tree->Branch("lbMassRec"   , &valueLbMassRec   , "lbMassRec/F" );
    tree->Branch("lbMassGen"   , &valueLbMassGen   , "lbMassGen/F" );
    tree->Branch("bbbarPtRec"  , &valueBbbarPtRec  , "bbbarPtRec/F");
    tree->Branch("bbbarPtGen"  , &valueBbbarPtGen  , "bbbarPtGen/F");
    tree->Branch("bbbarYRec"   , &valueBbbarYRec   , "bbbarYRec/F" );
    tree->Branch("bbbarYGen"   , &valueBbbarYGen   , "bbbarYGen/F" );
    tree->Branch("bbbarEtaRec" , &valueBbbarEtaRec , "bbbarEtaRec/F");
    tree->Branch("bbbarEtaGen" , &valueBbbarEtaGen , "bbbarEtaGen/F");
    tree->Branch("bqPtLeadRec" , &valueLeadBqPtRec , "bqPtLeadRec/F");
    tree->Branch("bqPtLeadGen" , &valueLeadBqPtGen , "bqPtLeadGen/F");
    tree->Branch("bqPtSubLeadRec", &valueSubLeadBqPtRec, "bqPtSubLeadRec/F");
    tree->Branch("bqPtSubLeadGen", &valueSubLeadBqPtGen, "bqPtSubLeadGen/F");
    // boolean = true if swapping gives better results
    tree->Branch("bbSwapBetter"  , &bbSwapBetter  , "bbSwapBetter/O" );
    tree->Branch("inVisPS"       , &inVisPS       , "inVisPS/O"      );
    // parton level quark vs reco jet assignment
    if(recPlots_&&genPlots_) tree->Branch("qAssignment" , &valueAssignment , "qAssignment/F" );
  }
}

void
SemiLepBjetAnalyzer::endJob() 
{
}

std::pair<const reco::Candidate*,const reco::Candidate*> 
SemiLepBjetAnalyzer::getbJets(const edm::Handle<TtSemiLeptonicEvent> semiLepEvt, const std::string hypoKey_){
  const reco::Candidate* zero=0;
  if(verbose>1) std::cout << "get b jet candidates from kinfit" << std::endl;
  // check if hypothesis is available and valid in this event
  if(!semiLepEvt->isHypoValid(hypoKey_)){
    if(verbose>0) std::cout << "Hypothesis not valid for this event" << std::endl;
  }
  else{
    // check if decay is semileptonic
    if(  semiLepEvt->genEvent().isAvailable()                    && 
	(semiLepEvt->genEvent()->isSemiLeptonic(WDecay::kMuon) || 
	 semiLepEvt->genEvent()->isSemiLeptonic(WDecay::kElec) ||
	 semiLepEvt->genEvent()->isSemiLeptonic(WDecay::kTau )   ) ){
      if(verbose>1) std::cout << "semilep. decay or no ttgenevent available" << std::endl;
      // get reconstructed b quarks (jets) from the hypothesis 
      const reco::Candidate* bhad  = semiLepEvt->hadronicDecayB(hypoKey_);
      const reco::Candidate* blep  = semiLepEvt->leptonicDecayB(hypoKey_);
      if(bhad&&blep){
	// destinguish b and bbar using the lepton charge
	// (tbar->W-&bbar, W- -> l- nu, therefore: l- -> bbar; l+ -> b)
	bool lminus = ((reco::LeafCandidate*)(semiLepEvt->singleLepton(hypoKey_)))->charge()<0 ? true : false;
	const reco::Candidate* b    = lminus ? bhad : blep;
	const reco::Candidate* bbar = lminus ? blep : bhad;   
	if(verbose>2) std::cout << "q(lepRec)<0? " << lminus << std::endl;
	// return b jets
	if(verbose>0) std::cout << "rec bjet candidates found!" << std::endl;
	return std::make_pair(b,bbar);
      }
    }
  }
  // otherwise: return NULL
  if(verbose>0) std::cout << "no rec bjets found" << std::endl;
  return std::make_pair(zero,zero);
}

std::pair<const reco::Candidate*,const reco::Candidate*> 
SemiLepBjetAnalyzer::getPreBJets(const edm::Handle<TtSemiLeptonicEvent> semiLepEvt, const std::string hypoKey_, const edm::Handle<std::vector<pat::Jet> > jets){
  const reco::Candidate* zero=0;
  if(verbose>1) std::cout << "get b jet candidate kinematics before KinFit" << std::endl;
  // check if kinfit information is available
  if(semiLepEvt.isValid()&&semiLepEvt->isHypoValid(hypoKey_)){
    if(verbose>0) std::cout << "Hypothesis is valid" << std::endl;
    // get indices
    int hadBIndex=semiLepEvt->jetLeptonCombination(hypoKey_)[TtSemiLepEvtPartons::HadB];
    int lepBIndex=semiLepEvt->jetLeptonCombination(hypoKey_)[TtSemiLepEvtPartons::LepB];
    if(verbose>1) std::cout << "recoLepBIX: " << lepBIndex << ", recoHadBIX: " << hadBIndex << std::endl;
    // get jets before the fit
    const reco::Candidate* bhad = &(jets->at(hadBIndex));
    const reco::Candidate* blep = &(jets->at(lepBIndex));
    // get charge of lepton 
    bool lminus = ((reco::LeafCandidate*)(semiLepEvt->singleLepton(hypoKey_)))->charge()<0 ? true : false;
    // return candidates
    // pick b/bbar wrt lepton charge
    // q(lep) < 0 -> bjet(leptonic)>0 -> bjet(leptonic)=bbar
    if(bhad&&blep){
      if(verbose>0) std::cout << "rec bjet candidates before KinFit found!" << std::endl;
      return std::make_pair((lminus ? bhad : blep),(lminus ? blep : bhad));
    }
  }
  return std::make_pair(zero,zero);
}

const reco::GenJet* 
SemiLepBjetAnalyzer::getJetFromCollection(const reco::GenJetCollection& genJets, int IX)
{
  const reco::GenJet* zero=0;
  int index=0;
  // loop collection
  for (size_t iJet = 0; iJet < genJets.size(); ++iJet) {
    const reco::GenJet* thisJet = &(genJets[iJet]);
    // return jet with index "IX" 
    if(index==IX){
      if(verbose>1) std::cout << "jet with index " << IX << " found!" << std::endl;
      return thisJet;
    }
    ++ index;
  }
  // if not found: return NULL
  if(verbose>0){
    std::cout << "WARNING: index " << IX << " is larger than the size";
    std::cout << " of the jet collection (" << genJets.size() << ")" << std::endl;
  }
  return zero;
}

double  
SemiLepBjetAnalyzer::checkPartonAssignment(const edm::Handle<TtSemiLeptonicEvent> semiLepEvt, int maxNJets){
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
