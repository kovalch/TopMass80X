#include "FWCore/Framework/interface/InputSourceDescription.h"
#include "FWCore/Framework/interface/EventPrincipal.h"
#include "FWCore/Sources/interface/VectorInputSourceFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/TypeID.h"

#include <algorithm>
#include "boost/bind.hpp"

#include "TopMass/TopEventTree/plugins/JetEventMixer.h"

JetEventMixer::JetEventMixer(const edm::ParameterSet& cfg) :
eventSrc_(0),
nMix_(cfg.getParameter<int>("nMix")),
combo_(0),
oriPatJets_(0),
//oriGenJets_(0),
oriPFCandidates_(0),
oriPUInfos_(0)
{
  if(nMix_ > 12){
    throw edm::Exception( edm::errors::Configuration, "Cannot run JetEventMixer with with nMix > 12" );
  }

  eventSrc_ = edm::VectorInputSourceFactory::get()->makeVectorInputSource(cfg.getParameterSet("input"), edm::InputSourceDescription()).release();
  // DOES NOT WORK, skip() IS NOT IMPLEMENTED IN edm::VectorInputSource
  //if(cfg.getParameterSet("input").exists("skipEvents")){
  //  unsigned int offset = cfg.getParameterSet("input").getUntrackedParameter<unsigned int>("skipEvents");
  //  std::cout << offset << std::endl;
  //  eventSrc_->skipEvents(offset);
  //}

  produces<std::vector<pat::Jet> >();
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

  unsigned int cSize = combo_.size();
  std::cout << cSize << ", Mix: ";
  for(unsigned int i = 0; i < cSize; ++i)
    std::cout << combo_[i] << " ";
  std::cout << std::endl;

  putOneEvent(evt);

  bool nextCombo = std::next_permutation(combo_.begin(), combo_.end());

  // clean up everything before putting new stuff into the event
  if(!nextCombo){
    cleanUp();
    getEvents(evt);
  }
}

void
JetEventMixer::getEvents(edm::Event& evt)
{
  combo_.clear();
  for(unsigned int i = 0; i < nMix_; ++i)
    combo_.push_back(i);
  for(unsigned int iEvt = 0; iEvt < nMix_; ++iEvt){
    //std::vector< std::string > wantedBranches;
    //eventSrc_->dropUnwantedBranches(wantedBranches);
    //int i =
    eventSrc_->loopSequential(1, boost::bind(&JetEventMixer::processOneEvent, this, _1, boost::ref(evt)));
    //std::cout << iEvt << " " << i << std::endl;
  }
  for(unsigned int iEvt = 0; iEvt < nMix_; ++iEvt){
    if(oriPatJets_[iEvt].size() < iEvt+1){
      combo_.erase(combo_.begin()+iEvt, combo_.end());
      break;
    }
  }
}

void
JetEventMixer::processOneEvent(edm::EventPrincipal const& eventPrincipal, edm::Event& evt)
{
  std::cout << "Run: "     << eventPrincipal.aux().run()             << " -> " << evt.eventAuxiliary().run()             ;
  std::cout << ", Lumi: "  << eventPrincipal.aux().luminosityBlock() << " -> " << evt.eventAuxiliary().luminosityBlock() ;
  std::cout << ", Event: " << eventPrincipal.aux().event()           << " -> " << evt.eventAuxiliary().event()           << std::endl;

  //for(edm::Principal::const_iterator it = eventPrincipal.begin(), itEnd = eventPrincipal.end(); it != itEnd; ++it){
  //  std::cout << (*it)->moduleLabel() << " : " << (*it)->productInstanceName() << " : " << (*it)->processName() << " : " << (*it)->productType().className() << std::endl;
  //}

  size_t cachedOffset;
  int fillCount;
  edm::BasicHandle hPatJets = eventPrincipal.getByLabel(edm::TypeID(typeid(std::vector<pat::Jet>)), "tightLeadingJets", "", "FullHadTreeWriter", cachedOffset, fillCount);
  //std::cout << "cachedOffset: " << cachedOffset << ", fillCount: " << fillCount << std::endl;
  //std::cout << "failedToGet: " << hPatJets.failedToGet() << std::endl;
  edm::Wrapper<std::vector<pat::Jet> > const* wPatJets = static_cast<edm::Wrapper<std::vector<pat::Jet> > const*>(hPatJets.wrapper());
  //std::cout << "wPatJets: " << wPatJets << std::endl;
  std::vector<pat::Jet> pPatJets = *wPatJets->product();
  //std::cout << "pPatJets: " << pPatJets << std::endl;
  //std::cout << "pPatJets: " << pPatJets->at(0).pt() << std::endl;
  oriPatJets_.push_back(pPatJets);
  //oriPatJets_.push_back(*(static_cast<edm::Wrapper<std::vector<pat::Jet> > const*>(hPatJets.wrapper())->product()));

  ////edm::BasicHandle hGenJets = eventPrincipal.getByLabel(edm::TypeID(typeid(std::vector<reco::GenJet>)), "tightLeadingJets", "genJets", "FullHadTreeWriter", cachedOffset, fillCount);
  //edm::BasicHandle hGenJets = eventPrincipal.getByType(edm::TypeID(typeid(std::vector<reco::GenJet>)));
  //edm::Wrapper<std::vector<reco::GenJet> > const* wGenJets = static_cast<edm::Wrapper<std::vector<reco::GenJet> > const*>(hGenJets.wrapper());
  //std::vector<reco::GenJet> pGenJets = *wGenJets->product();
  //oriGenJets_.push_back(pGenJets);

  //edm::BasicHandle hPFCandidates = eventPrincipal.getByLabel(edm::TypeID(typeid(std::vector<reco::PFCandidate>)), "tightLeadingJets", "pfCandidates", "FullHadTreeWriter", cachedOffset, fillCount);
  edm::BasicHandle hPFCandidates = eventPrincipal.getByType(edm::TypeID(typeid(std::vector<reco::PFCandidate>)));
  edm::Wrapper<std::vector<reco::PFCandidate> > const* wPFCandidates = static_cast<edm::Wrapper<std::vector<reco::PFCandidate> > const*>(hPFCandidates.wrapper());
  std::vector<reco::PFCandidate> pPFCandidates = *wPFCandidates->product();
  oriPFCandidates_.push_back(pPFCandidates);

  //edm::BasicHandle hPileupSummaryInfos = eventPrincipal.getByLabel(edm::TypeID(typeid(std::vector<reco::PileupSummaryInfo>)), "tightLeadingJets", "pfCandidates", "FullHadTreeWriter", cachedOffset, fillCount);
  edm::BasicHandle hPileupSummaryInfos = eventPrincipal.getByType(edm::TypeID(typeid(std::vector<PileupSummaryInfo>)));
  edm::Wrapper<std::vector<PileupSummaryInfo> > const* wPileupSummaryInfos = static_cast<edm::Wrapper<std::vector<PileupSummaryInfo> > const*>(hPileupSummaryInfos.wrapper());
  std::vector<PileupSummaryInfo> pPileupSummaryInfos = *wPileupSummaryInfos->product();
  oriPUInfos_.push_back(pPileupSummaryInfos);
}

void
JetEventMixer::putOneEvent(edm::Event& evt)
{
  std::auto_ptr< std::vector<pat::Jet> > patJets ( new std::vector<pat::Jet>() );

  //std::auto_ptr<reco::GenJetCollection > genJetsOut ( new reco::GenJetCollection() );
  //std::auto_ptr<std::vector<CaloTower>  >  caloTowersOut( new std::vector<CaloTower> () );
  std::auto_ptr<reco::PFCandidateCollection > pfCandidatesOut( new reco::PFCandidateCollection() );
  //std::auto_ptr<edm::OwnVector<reco::BaseTagInfo> > tagInfosOut ( new edm::OwnVector<reco::BaseTagInfo>() );

  //edm::RefProd<reco::GenJetCollection > h_genJetsOut = evt.getRefBeforePut<reco::GenJetCollection >( "genJets" );
  //edm::RefProd<std::vector<CaloTower>  >  h_caloTowersOut = evt.getRefBeforePut<std::vector<CaloTower>  > ( "caloTowers" );
  edm::RefProd<reco::PFCandidateCollection > h_pfCandidatesOut = evt.getRefBeforePut<reco::PFCandidateCollection > ( "pfCandidates" );
  //edm::RefProd<edm::OwnVector<reco::BaseTagInfo> > h_tagInfosOut = evt.getRefBeforePut<edm::OwnVector<reco::BaseTagInfo> > ( "tagInfos" );

  std::auto_ptr< std::vector<PileupSummaryInfo> > puInfos ( new std::vector<PileupSummaryInfo>() );

  int nPU0 = 0, nPU1 = 0, nPU2 = 0;
  double nPUTrue0 = 0, nPUTrue1 = 0, nPUTrue2 = 0;
  for(unsigned int i = 0; i < nMix_; ++i){

    // update the forward pointers to the PFCandidates of the jet
    unsigned int jetIdx = 0;
    unsigned int pfCandidateIndex = 0;
    for(std::vector<reco::PFCandidate>::const_iterator iPart = oriPFCandidates_[i].begin(), partEnd = oriPFCandidates_[i].end(); iPart != partEnd; ++iPart, ++pfCandidateIndex){
      //std::cout << pfCandidateIndex << " " << iPart - oriPFCandidates_[i].begin() << " " << oriPatJets_[i][jetIdx].pfCandidatesFwdPtr().size() << " " << oriPFCandidates_[i].size() << std::endl;
      if(pfCandidateIndex == oriPatJets_[i][jetIdx].pfCandidatesFwdPtr().size()){
        ++jetIdx;
        pfCandidateIndex = 0;
        //std::cout << "NEXT JET!!!" << std::endl;
      }
      if(jetIdx >= oriPatJets_[i].size()){
        std::cout << "Something went wrong, jetIdx (" << jetIdx << ") >= nJets (" << oriPatJets_[i].size() << ")!" << std::endl;
        break;
      }
      //if( combo_[i] == jetIdx ){
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
    for(std::vector<PileupSummaryInfo>::const_iterator iterPU = oriPUInfos_[i].begin(), iterEnd = oriPUInfos_[i].end(); iterPU != iterEnd; ++iterPU){
      // vector size is 3
      // -1: previous BX, 0: current BX,  1: next BX
      if     (iterPU->getBunchCrossing()==-1) { nPU0 += iterPU->getPU_NumInteractions(); nPUTrue0 += iterPU->getTrueNumInteractions(); }
      else if(iterPU->getBunchCrossing()== 0) { nPU1 += iterPU->getPU_NumInteractions(); nPUTrue1 += iterPU->getTrueNumInteractions(); }
      else if(iterPU->getBunchCrossing()== 1) { nPU2 += iterPU->getPU_NumInteractions(); nPUTrue2 += iterPU->getTrueNumInteractions(); }
    }
  }
  nPU0 /= nMix_;
  nPU1 /= nMix_;
  nPU2 /= nMix_;
  nPUTrue0 /= nMix_;
  nPUTrue1 /= nMix_;
  nPUTrue2 /= nMix_;
  std::vector< float > zpositions, sumpT_lowpT, sumpT_highpT;
  std::vector< int > ntrks_lowpT, ntrks_highpT;
  puInfos->push_back(PileupSummaryInfo(nPU0, zpositions, sumpT_lowpT, sumpT_highpT, ntrks_lowpT, ntrks_highpT, -1, nPUTrue0));
  puInfos->push_back(PileupSummaryInfo(nPU1, zpositions, sumpT_lowpT, sumpT_highpT, ntrks_lowpT, ntrks_highpT,  0, nPUTrue1));
  puInfos->push_back(PileupSummaryInfo(nPU2, zpositions, sumpT_lowpT, sumpT_highpT, ntrks_lowpT, ntrks_highpT, +1, nPUTrue2));

  if(combo_.size() > 0 && oriPatJets_[combo_[0]].size() > 0) patJets->push_back(oriPatJets_[combo_[0]][0]);
  if(combo_.size() > 1 && oriPatJets_[combo_[1]].size() > 1) patJets->push_back(oriPatJets_[combo_[1]][1]);
  if(combo_.size() > 2 && oriPatJets_[combo_[2]].size() > 2) patJets->push_back(oriPatJets_[combo_[2]][2]);
  if(combo_.size() > 3 && oriPatJets_[combo_[3]].size() > 3) patJets->push_back(oriPatJets_[combo_[3]][3]);
  if(combo_.size() > 4 && oriPatJets_[combo_[4]].size() > 4) patJets->push_back(oriPatJets_[combo_[4]][4]);
  if(combo_.size() > 5 && oriPatJets_[combo_[5]].size() > 5) patJets->push_back(oriPatJets_[combo_[5]][5]);
  if(combo_.size() > 6 && oriPatJets_[combo_[6]].size() > 6) patJets->push_back(oriPatJets_[combo_[6]][6]);
  if(combo_.size() > 7 && oriPatJets_[combo_[7]].size() > 7) patJets->push_back(oriPatJets_[combo_[7]][7]);
  if(combo_.size() > 8 && oriPatJets_[combo_[8]].size() > 8) patJets->push_back(oriPatJets_[combo_[8]][8]);
  if(combo_.size() > 9 && oriPatJets_[combo_[9]].size() > 9) patJets->push_back(oriPatJets_[combo_[9]][9]);
  evt.put(patJets, "");
  //evt.put(genJetsOut, "genJets");
  //evt.put(caloTowersOut, "caloTowers");
  evt.put(pfCandidatesOut, "pfCandidates");
  //evt.put(tagInfosOut, "tagInfos");
  evt.put(puInfos, "addPileupInfo");
}

void
JetEventMixer::cleanUp()
{
  //std::cout << "CLEANING UP !!!" << std::endl;
  if(oriPatJets_.size()){
    oriPatJets_     .clear();
    //oriGenJets_     .clear();
    oriPFCandidates_.clear();
    oriPUInfos_     .clear();
  }
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetEventMixer);
