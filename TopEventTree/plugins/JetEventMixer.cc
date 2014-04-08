#include "FWCore/Framework/interface/InputSourceDescription.h"
#include "FWCore/Framework/interface/EventPrincipal.h"
#include "FWCore/Sources/interface/VectorInputSourceFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/TypeID.h"
//#include "FWCore/Utilities/interface/EDMException.h"

#include <algorithm>
#include <csignal>
#include "boost/bind.hpp"

#include "TRandom3.h"

#include "TopMass/TopEventTree/plugins/JetEventMixer.h"

JetEventMixer::JetEventMixer(const edm::ParameterSet& cfg) :
eventSrc_(0),
nMix_   (cfg.getParameter<int>("nMix")),
nMixMin_(cfg.getParameter<int>("nMixMin")),
speedUp_(cfg.getParameter<int>("speedUp")),
comboIndex_(0),
combos_(0),
validCombos_(0),
oriPatJets_(0),
oriPatJetsCalo_(0),
//oriGenJets_(0),
oriPFCandidates_(0),
oriPUInfos_(0)
{
  if(nMix_ > 10){
    throw edm::Exception( edm::errors::Configuration, "Cannot run JetEventMixer with with nMix > 10" );
  }

  eventSrc_ = edm::VectorInputSourceFactory::get()->makeVectorInputSource(cfg.getParameterSet("input"), edm::InputSourceDescription());
  // DOES NOT WORK, skip() IS NOT IMPLEMENTED IN edm::VectorInputSource
  //if(cfg.getParameterSet("input").exists("skipEvents")){
  //  unsigned int offset = cfg.getParameterSet("input").getUntrackedParameter<unsigned int>("skipEvents");
  //  std::cout << offset << std::endl;
  //  eventSrc_->skipEvents(offset);
  //}

  produces<std::vector<pat::Jet> >();
  produces<std::vector<pat::Jet> >("calo");
  //produces<std::vector<reco::GenJet> > ("genJets");
  //produces<std::vector<CaloTower>  > ("caloTowers");
  produces<std::vector<reco::PFCandidate> > ("pfCandidates");
  //produces<edm::OwnVector<reco::BaseTagInfo> > ("tagInfos");
  produces<std::vector<PileupSummaryInfo> >("addPileupInfo");
}

JetEventMixer::~JetEventMixer()
{
}

void
JetEventMixer::produce(edm::Event& evt, const edm::EventSetup& setup)
{
  // fetch the events once, before starting to iterate the combinations
  if(!oriPatJets_.size())
    getEvents(evt);

  putOneEvent(evt);
  ++comboIndex_;

  // clean up everything before putting new stuff into the event
  if(comboIndex_ == combos_.size()){
    cleanUp();
    getEvents(evt);
  }
}

void
JetEventMixer::getEvents(edm::Event& evt)
{
  while(oriPatJets_.size() < nMix_){
  //for(unsigned int iEvt = 0; iEvt < nMix_; ++iEvt){
    //std::vector< std::string > wantedBranches;
    //eventSrc_->dropUnwantedBranches(wantedBranches);

    // exit program gracefully once less then nMix events are available
    if(!eventSrc_->loopSequential(1, boost::bind(&JetEventMixer::processOneEvent, this, _1, boost::ref(evt)))){
      cleanUp();
      //std::stringstream errorMessage;
      //errorMessage << "Less than nMix (" << nMix_ << ") events left in source to be processed." << "\n" << "Terminating program!";
      //throw edm::Exception( edm::errors::UnimplementedFeature, errorMessage.str() );
      //throw NotEnoughEventsLeftException("NotEnoughEventsLeftException", errorMessage.str());
      edm::LogWarning("NotEnoughEventsLeftException") << __FILE__ << ":" << "\n"
          << __FUNCTION__ << ":" << __LINE__ << "\n"
          << "Less than nMix (" << nMix_ << ") events left in source to be processed." << "\n"
          << "Terminating program!";
      kill(getpid(),SIGUSR2);
      //kill(getpid(),SIGINT);
      return;
    }
  }
  fillCombos();
}

void
JetEventMixer::processOneEvent(edm::EventPrincipal& eventPrincipal, edm::Event& evt)
{
  //std::cout << "Run: "     << eventPrincipal.aux().run()             << " -> " << evt.eventAuxiliary().run()             ;
  //std::cout << ", Lumi: "  << eventPrincipal.aux().luminosityBlock() << " -> " << evt.eventAuxiliary().luminosityBlock() ;
  //std::cout << ", Event: " << eventPrincipal.aux().event()           << " -> " << evt.eventAuxiliary().event()           << std::endl;

  //for(edm::Principal::const_iterator it = eventPrincipal.begin(), itEnd = eventPrincipal.end(); it != itEnd; ++it){
  //  std::cout << (*it)->moduleLabel() << " : " << (*it)->productInstanceName() << " : " << (*it)->processName() << " : " << (*it)->productType().className() << std::endl;
  //}

  size_t cachedOffset2 = 0;
  int fillCount2 = 0;
  edm::BasicHandle hPatJetsCalo = eventPrincipal.getByLabel(edm::TypeID(typeid(std::vector<pat::Jet>)), "selectedPatJetsAK5Calo", "", "FullHadTreeWriter", cachedOffset2, fillCount2);
  edm::Wrapper<std::vector<pat::Jet> > const* wPatJetsCalo = static_cast<edm::Wrapper<std::vector<pat::Jet> > const*>(hPatJetsCalo.wrapper());
  std::vector<pat::Jet> pPatJetsCalo = *wPatJetsCalo->product();
  if(pPatJetsCalo.size() > nMix_) pPatJetsCalo.resize(nMix_);
  if(pPatJetsCalo.size() < 4 || pPatJetsCalo[3].pt() < 60) return;
  oriPatJetsCalo_.push_back(pPatJetsCalo);
  //delete wPatJetsCalo;

  size_t cachedOffset = 0;
  int fillCount = 0;
  edm::BasicHandle hPatJets = eventPrincipal.getByLabel(edm::TypeID(typeid(std::vector<pat::Jet>)), "tightLeadingJets", "", "FullHadTreeWriter", cachedOffset, fillCount);
  //std::cout << "cachedOffset: " << cachedOffset << ", fillCount: " << fillCount << std::endl;
  //std::cout << "failedToGet: " << hPatJets.failedToGet() << std::endl;
  edm::Wrapper<std::vector<pat::Jet> > const* wPatJets = static_cast<edm::Wrapper<std::vector<pat::Jet> > const*>(hPatJets.wrapper());
  //std::cout << "wPatJets: " << wPatJets << std::endl;
  std::vector<pat::Jet> pPatJets = *wPatJets->product();
  if(pPatJets.size() > nMix_) pPatJets.resize(nMix_);
  //std::cout << "pPatJets: " << pPatJets << std::endl;
  //std::cout << "pPatJets: " << pPatJets->at(0).pt() << std::endl;
  oriPatJets_.push_back(pPatJets);
  //oriPatJets_.push_back(*(static_cast<edm::Wrapper<std::vector<pat::Jet> > const*>(hPatJets.wrapper())->product()));
  //delete wPatJets;

  ////edm::BasicHandle hGenJets = eventPrincipal.getByLabel(edm::TypeID(typeid(std::vector<reco::GenJet>)), "tightLeadingJets", "genJets", "FullHadTreeWriter", cachedOffset, fillCount);
  //edm::BasicHandle hGenJets = eventPrincipal.getByType(edm::TypeID(typeid(std::vector<reco::GenJet>)));
  //edm::Wrapper<std::vector<reco::GenJet> > const* wGenJets = static_cast<edm::Wrapper<std::vector<reco::GenJet> > const*>(hGenJets.wrapper());
  //std::vector<reco::GenJet> pGenJets = *wGenJets->product();
  //oriGenJets_.push_back(pGenJets);
  //delete wGenJets;

  //edm::BasicHandle hPFCandidates = eventPrincipal.getByLabel(edm::TypeID(typeid(std::vector<reco::PFCandidate>)), "tightLeadingJets", "pfCandidates", "FullHadTreeWriter", cachedOffset, fillCount);
  edm::BasicHandle hPFCandidates = eventPrincipal.getByType(edm::TypeID(typeid(std::vector<reco::PFCandidate>)));
  edm::Wrapper<std::vector<reco::PFCandidate> > const* wPFCandidates = static_cast<edm::Wrapper<std::vector<reco::PFCandidate> > const*>(hPFCandidates.wrapper());
  std::vector<reco::PFCandidate> pPFCandidates = *wPFCandidates->product();
  oriPFCandidates_.push_back(pPFCandidates);
  //delete wPFCandidates;

  //edm::BasicHandle hPileupSummaryInfos = eventPrincipal.getByLabel(edm::TypeID(typeid(std::vector<reco::PileupSummaryInfo>)), "tightLeadingJets", "pfCandidates", "FullHadTreeWriter", cachedOffset, fillCount);
  edm::BasicHandle hPileupSummaryInfos = eventPrincipal.getByType(edm::TypeID(typeid(std::vector<PileupSummaryInfo>)));
  //if(hPileupSummaryInfos.isValid()){
  if(false){
    edm::Wrapper<std::vector<PileupSummaryInfo> > const* wPileupSummaryInfos = static_cast<edm::Wrapper<std::vector<PileupSummaryInfo> > const*>(hPileupSummaryInfos.wrapper());
    std::vector<PileupSummaryInfo> pPileupSummaryInfos = *wPileupSummaryInfos->product();
    oriPUInfos_.push_back(pPileupSummaryInfos);
    //delete wPileupSummaryInfos;
  }
  eventPrincipal.clearEventPrincipal();
}

int factorial(int n)
{
  if(n > 12 || n < 0) return -1;

  int fact=1;
  for (int i=1; i<=n; ++i)
    fact*=i;
  return fact;
}

void
JetEventMixer::fillCombos()
{
  comboIndex_ = 0;
  validCombos_.clear();
  //static unsigned int maxCombosToSave = factorial(nMixMin_) / speedUp_;
  static TRandom3 rand(0);

  std::vector<char> nJets(0);
  for(std::vector<std::vector<pat::Jet> >::const_iterator jets = oriPatJets_.begin(), endJets = oriPatJets_.end(); jets != endJets; ++jets)
    nJets.push_back(jets->size());

  //std::cout << "maxCombosToSave: " << maxCombosToSave << ", rand: " << rand.Uniform() << ", nJets.size(): " << nJets.size() << ": ";
  int maxNJets = 0;
  for(std::vector<char>::const_iterator n = nJets.begin(); n < nJets.end(); ++n){
    //std::cout << int(*n) << ", ";
    if(*n > maxNJets) maxNJets = *n;
  }
  maxNJets = std::min(maxNJets, int(nMix_));
  //std::cout << maxNJets << std::endl;
  unsigned int maxCombosToSave = factorial(maxNJets) / speedUp_;

  std::vector<char> comb(0);
  for(unsigned int i = 0; i < nMix_; ++i)
    comb.push_back(i);
  comb.erase(comb.begin()+maxNJets,comb.end());

  do{
    for(char i = 0; i < maxNJets; ++i){
      if     (comb[i] >= nJets[comb[i]] && i+1 <  int(nMixMin_)) continue;
      else if(comb[i] >= nJets[comb[i]] && i+1 >= int(nMixMin_)){
        std::vector<char> tmpComb(comb.begin(), comb.begin()+i+1);
        validCombos_.push_back(tmpComb);
        break;
      }
      else if(maxNJets-1 == i){
        validCombos_.push_back(comb);
      }
    }
  } while(std::next_permutation(comb.begin(), comb.end()));

  //std::cout << "before: " << validCombos_.size() << std::flush;
  do{
    int size = validCombos_.size();
    int idx  = rand.Integer(size);
    combos_.push_back(validCombos_[idx]);
    std::swap(validCombos_[idx], validCombos_[size-1]);
    validCombos_.pop_back();
  } while(combos_.size() < maxCombosToSave);
  //std::cout << ", after: " << validCombos_.size() << ", combos_.size(): " << combos_.size() << std::endl;
}

void
JetEventMixer::putOneEvent(edm::Event& evt)
{
  std::vector<char> combo = combos_[comboIndex_];
  std::auto_ptr< std::vector<pat::Jet> > patJets ( new std::vector<pat::Jet>() );
  std::auto_ptr< std::vector<pat::Jet> > patJetsCalo ( new std::vector<pat::Jet>() );

  //std::auto_ptr<reco::GenJetCollection > genJetsOut ( new reco::GenJetCollection() );
  //std::auto_ptr<std::vector<CaloTower>  >  caloTowersOut( new std::vector<CaloTower> () );
  std::auto_ptr<reco::PFCandidateCollection > pfCandidatesOut( new reco::PFCandidateCollection() );
  //std::auto_ptr<edm::OwnVector<reco::BaseTagInfo> > tagInfosOut ( new edm::OwnVector<reco::BaseTagInfo>() );

  //edm::RefProd<reco::GenJetCollection > h_genJetsOut = evt.getRefBeforePut<reco::GenJetCollection >( "genJets" );
  //edm::RefProd<std::vector<CaloTower>  >  h_caloTowersOut = evt.getRefBeforePut<std::vector<CaloTower>  > ( "caloTowers" );
  edm::RefProd<reco::PFCandidateCollection > h_pfCandidatesOut = evt.getRefBeforePut<reco::PFCandidateCollection > ( "pfCandidates" );
  //edm::RefProd<edm::OwnVector<reco::BaseTagInfo> > h_tagInfosOut = evt.getRefBeforePut<edm::OwnVector<reco::BaseTagInfo> > ( "tagInfos" );

  std::auto_ptr< std::vector<PileupSummaryInfo> > puInfos ( new std::vector<PileupSummaryInfo>() );

  //int nPU0 = 0, nPU1 = 0, nPU2 = 0;
  //double nPUTrue0 = 0, nPUTrue1 = 0, nPUTrue2 = 0;
  for(unsigned int i = 0; i < nMix_; ++i){

    // update the forward pointers to the PFCandidates of the jet
    unsigned int jetIdx = 0;
    unsigned int pfCandidateIndex = 0;
    for(std::vector<reco::PFCandidate>::const_iterator iPart = oriPFCandidates_[i].begin(), partEnd = oriPFCandidates_[i].end(); iPart != partEnd; ++iPart, ++pfCandidateIndex){
      //std::cout << pfCandidateIndex << " " << iPart - oriPFCandidates_[i].begin() << " " << oriPatJets_[i][jetIdx].pfCandidatesFwdPtr().size() << " " << oriPFCandidates_[i].size() << std::endl;
      //std::vector<pat::Jet>& myVPJs = oriPatJets_[i];
      //pat::Jet& myPJ = myVPJs[jetIdx];
      //const std::vector<reco::PFCandidateFwdPtr> myVFP = myPJ.pfCandidatesFwdPtr();
      //unsigned int mySize = myVFP.size();
      //if(pfCandidateIndex == mySize){
      if(pfCandidateIndex == oriPatJets_[i][jetIdx].pfCandidatesFwdPtr().size()){
        ++jetIdx;
        pfCandidateIndex = 0;
        //std::cout << "NEXT JET!!!" << std::endl;
      }
      if(jetIdx >= nMix_) continue;
      if(jetIdx >= oriPatJets_[i].size()){
        std::cout << "Something went wrong, jetIdx (" << jetIdx << ") >= nJets (" << oriPatJets_[i].size() << ")!" << std::endl;
        break;
      }
      //if( combo[i] == jetIdx ){
      pfCandidatesOut->push_back( *iPart );
      edm::Ref<reco::PFCandidateCollection> pfCollectionRef( h_pfCandidatesOut, pfCandidatesOut->size() - 1);
      edm::Ptr<reco::PFCandidate> pfForwardRef( h_pfCandidatesOut.id(), pfCollectionRef.key(),  h_pfCandidatesOut.productGetter() );
      oriPatJets_[i][jetIdx].updateFwdPFCandidateFwdPtr(pfCandidateIndex, pfForwardRef);
      std::vector<edm::FwdPtr<reco::PFCandidate> > emptyDummy(0);
      oriPatJets_[i][jetIdx].setPFCandidates(emptyDummy);
      //}
    }

    //// update the forward references to the GenJet of the jet
    //jetIdx = 0;
    //unsigned int genJetIndex = 0;
    //for(std::vector<reco::GenJet>::const_iterator iGenJet = oriGenJets_[i].begin(), genJetsEnd = oriGenJets_[i].end(); iGenJet != genJetsEnd; ++iGenJet, ++genJetIndex, ++jetIdx){
    //  if(oriPatJets_[i][jetIdx].genJetFwdRef().isNonnull()){
    //    ++jetIdx;
    //  }
    //  if( i == 0 /*jetIdx*/ ){
    //    genJetsOut->push_back( *iGenJet );
    //    edm::Ref<std::vector<reco::GenJet> > genRef( h_genJetsOut, genJetsOut->size() - 1 );
    //    ////edm::Ref<std::vector<reco::GenJet> > genRefZero;
    //    ////edm::FwdRef<std::vector<reco::GenJet> > genForwardRef( genRef, genRefZero );
    //    //oriPatJets_[i][jetIdx].updateFwdGenJetFwdRef(genRef);
    //    ////oriPatJets_[i][jetIdx].setGenJetRef(genForwardRef);
    //  }
    //}
    for(std::vector<pat::Jet>::iterator iJet = oriPatJets_[i].begin(), patJetsEnd = oriPatJets_[i].end(); iJet != patJetsEnd; ++iJet){
      edm::FwdRef<std::vector<reco::GenJet> > zeroGenJetForwardRef;
      iJet->setGenJetRef(zeroGenJetForwardRef);
    }

    // create new PileupSummaryInfo for the combined event
    //if(oriPUInfos_.size()){
    //  for(std::vector<PileupSummaryInfo>::const_iterator iterPU = oriPUInfos_[i].begin(), iterEnd = oriPUInfos_[i].end(); iterPU != iterEnd; ++iterPU){
    //    // vector size is 3
    //    // -1: previous BX, 0: current BX,  1: next BX
    //    if     (iterPU->getBunchCrossing()==-1) { nPU0 += iterPU->getPU_NumInteractions(); nPUTrue0 += iterPU->getTrueNumInteractions(); }
    //    else if(iterPU->getBunchCrossing()== 0) { nPU1 += iterPU->getPU_NumInteractions(); nPUTrue1 += iterPU->getTrueNumInteractions(); }
    //    else if(iterPU->getBunchCrossing()== 1) { nPU2 += iterPU->getPU_NumInteractions(); nPUTrue2 += iterPU->getTrueNumInteractions(); }
    //  }
    //}
  }
  //if(oriPUInfos_.size()){
  //  nPU0 /= nMix_;
  //  nPU1 /= nMix_;
  //  nPU2 /= nMix_;
  //  nPUTrue0 /= nMix_;
  //  nPUTrue1 /= nMix_;
  //  nPUTrue2 /= nMix_;
  //  std::vector< float > zpositions, sumpT_lowpT, sumpT_highpT;
  //  std::vector< int > ntrks_lowpT, ntrks_highpT;
  //  puInfos->push_back(PileupSummaryInfo(nPU0, zpositions, sumpT_lowpT, sumpT_highpT, ntrks_lowpT, ntrks_highpT, -1, nPUTrue0));
  //  puInfos->push_back(PileupSummaryInfo(nPU1, zpositions, sumpT_lowpT, sumpT_highpT, ntrks_lowpT, ntrks_highpT,  0, nPUTrue1));
  //  puInfos->push_back(PileupSummaryInfo(nPU2, zpositions, sumpT_lowpT, sumpT_highpT, ntrks_lowpT, ntrks_highpT, +1, nPUTrue2));
  //}

  if(combo.size() > 0 && oriPatJets_[combo[0]].size() > 0) { patJets->push_back(oriPatJets_[combo[0]][0]); if(oriPatJetsCalo_[combo[0]].size() > 0) patJetsCalo->push_back(oriPatJetsCalo_[combo[0]][0]); }
  if(combo.size() > 1 && oriPatJets_[combo[1]].size() > 1) { patJets->push_back(oriPatJets_[combo[1]][1]); if(oriPatJetsCalo_[combo[1]].size() > 1) patJetsCalo->push_back(oriPatJetsCalo_[combo[1]][1]); }
  if(combo.size() > 2 && oriPatJets_[combo[2]].size() > 2) { patJets->push_back(oriPatJets_[combo[2]][2]); if(oriPatJetsCalo_[combo[2]].size() > 2) patJetsCalo->push_back(oriPatJetsCalo_[combo[2]][2]); }
  if(combo.size() > 3 && oriPatJets_[combo[3]].size() > 3) { patJets->push_back(oriPatJets_[combo[3]][3]); if(oriPatJetsCalo_[combo[3]].size() > 3) patJetsCalo->push_back(oriPatJetsCalo_[combo[3]][3]); }
  if(combo.size() > 4 && oriPatJets_[combo[4]].size() > 4) { patJets->push_back(oriPatJets_[combo[4]][4]); if(oriPatJetsCalo_[combo[4]].size() > 4) patJetsCalo->push_back(oriPatJetsCalo_[combo[4]][4]); }
  if(combo.size() > 5 && oriPatJets_[combo[5]].size() > 5) { patJets->push_back(oriPatJets_[combo[5]][5]); if(oriPatJetsCalo_[combo[5]].size() > 5) patJetsCalo->push_back(oriPatJetsCalo_[combo[5]][5]); }
  if(combo.size() > 6 && oriPatJets_[combo[6]].size() > 6) { patJets->push_back(oriPatJets_[combo[6]][6]); if(oriPatJetsCalo_[combo[6]].size() > 6) patJetsCalo->push_back(oriPatJetsCalo_[combo[6]][6]); }
  if(combo.size() > 7 && oriPatJets_[combo[7]].size() > 7) { patJets->push_back(oriPatJets_[combo[7]][7]); if(oriPatJetsCalo_[combo[7]].size() > 7) patJetsCalo->push_back(oriPatJetsCalo_[combo[7]][7]); }
  if(combo.size() > 8 && oriPatJets_[combo[8]].size() > 8) { patJets->push_back(oriPatJets_[combo[8]][8]); if(oriPatJetsCalo_[combo[8]].size() > 8) patJetsCalo->push_back(oriPatJetsCalo_[combo[8]][8]); }
  if(combo.size() > 9 && oriPatJets_[combo[9]].size() > 9) { patJets->push_back(oriPatJets_[combo[9]][9]); if(oriPatJetsCalo_[combo[9]].size() > 9) patJetsCalo->push_back(oriPatJetsCalo_[combo[9]][9]); }
  evt.put(patJets, "");
  evt.put(patJetsCalo, "calo");
  //evt.put(genJetsOut, "genJets");
  //evt.put(caloTowersOut, "caloTowers");
  evt.put(pfCandidatesOut, "pfCandidates");
  //evt.put(tagInfosOut, "tagInfos");
  if(oriPUInfos_.size()) evt.put(puInfos, "addPileupInfo");
}

void
JetEventMixer::cleanUp()
{
  //std::cout << "CLEANING UP !!!" << std::endl;
  //if(oriPatJets_.size()){
    oriPatJets_     .clear();
    oriPatJetsCalo_ .clear();
    //oriGenJets_     .clear();
    oriPFCandidates_.clear();
    oriPUInfos_     .clear();
  //}
  combos_.clear();
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetEventMixer);
