/*
 * JetEventAnalyzer.cc
 *
 *  Created on: Feb 6, 2013
 *      Author: eschliec
 */

//#include <memory>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"


#include "TopMass/TopEventTree/interface/TreeRegistryService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputerWrapper.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputerRecord.h"
#include "RecoBTag/SecondaryVertex/interface/CombinedSVComputer.h"

#include "TopMass/TopEventTree/plugins/JetEventAnalyzer.h"

JetEventAnalyzer::JetEventAnalyzer(const edm::ParameterSet& cfg):
jets_           (cfg.getParameter<edm::InputTag>("jets")),
alternativeJets_(cfg.getParameter<edm::InputTag>("alternativeJets")),
met_            (cfg.getParameter<edm::InputTag>("met")),
gluonTagName_   (cfg.getParameter<edm::InputTag>("gluonTagSrc").encode()),
kJetMAX_(cfg.getParameter<int>("maxNJets")),
jet(0),
checkedIsPFJet(false), checkedJERSF(false), checkedJESSF(false), checkedTotalSF(false), checkedQGTag(false), checkedBReg(false),
       isPFJet(false),     hasJERSF(false),     hasJESSF(false),     hasTotalSF(false),     hasQGTag(false),     hasBReg(false),
alternativeJetsAvailable(alternativeJets_.label().size())
{
}

void
JetEventAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{
  //////////////////////////////////////////////////////////////////////////////
  // INIT JetEvent
  ////////////////////////////////////////////////////////////////////////////

  jet->init();

  //////////////////////////////////////////////////////////////////////////
  // JETS
  ////////////////////////////////////////////////////////////////////////

  edm::Handle<std::vector<pat::Jet> > jets;
  evt.getByLabel(jets_, jets);

  //////////////////////////////////////////////////////////////////////// 
  // instantiate a tagging variable computer for unification of some calculations like vertex mass corrections
  ////////////////////////////////////////////////////////////////////////
  //edm::ESHandle<JetTagComputer> computerHandle;;
  //setup.get<JetTagComputerRecord>().get( "combinedSecondaryVertex", computerHandle );
  //const GenericMVAJetTagComputer *computer = dynamic_cast<const GenericMVAJetTagComputer*>( computerHandle.product() );
  //if (!computer){
  // edm::LogError("DataLost")<<"computer missing !!!"<<std::endl;
  //  exit(1);
  //}
  //computer->passEventSetup(setup);
  //use TagInfos in future!!!!

  unsigned short jetIndex = 0;
  for(std::vector< pat::Jet >::const_iterator ijet = jets->begin(); ijet != jets->end(); ++ijet, ++jetIndex) {
    // write only kJetMAX_ jets into the event
    if(jetIndex == kJetMAX_) break;

    jet->jet.push_back(TLorentzVector(ijet->px(), ijet->py(), ijet->pz(), ijet->energy()));
    if(ijet->genJet()) jet->genJet.push_back(TLorentzVector(ijet->genJet()->px(), ijet->genJet()->py(), ijet->genJet()->pz(), ijet->genJet()->energy()));
    else               jet->genJet.push_back(TLorentzVector(0,0,0,-1));
    if(ijet->genParton()) jet->genParton.push_back(TLorentzVector(ijet->genParton()->px(), ijet->genParton()->py(), ijet->genParton()->pz(), ijet->genParton()->energy()));
    else                  jet->genParton.push_back(TLorentzVector(0,0,0,-1));

    // check only once per module run if the needed collections are available
    if(!checkedIsPFJet) { checkedIsPFJet = true; isPFJet = ijet->isPFJet(); }
    if(!checkedQGTag  ) { checkedQGTag   = true; if(ijet->hasUserFloat(gluonTagName_)) hasQGTag   = true; }
    if(!checkedJERSF  ) { checkedJERSF   = true; if(ijet->hasUserFloat("jerSF"      )) hasJERSF   = true; }
    if(!checkedJESSF  ) { checkedJESSF   = true; if(ijet->hasUserFloat("jesSF"      )) hasJESSF   = true; }
    if(!checkedTotalSF) { checkedTotalSF = true; if(ijet->hasUserFloat("totalSF"    )) hasTotalSF = true; }
    if(!checkedBReg   ) { checkedBReg    = true; if(ijet->hasUserFloat("BRegResult" )) hasBReg    = true; }

    if(isPFJet){
      jet->fChargedHadron .push_back(ijet->chargedHadronEnergyFraction());
      jet->fNeutralHadron .push_back(ijet->neutralHadronEnergyFraction());
      jet->fElectron      .push_back(ijet->electronEnergyFraction());
      jet->fPhoton        .push_back(ijet->photonEnergyFraction());
      jet->fMuon          .push_back(ijet->muonEnergyFraction());
      jet->nConstituents  .push_back(ijet->nConstituents());
      jet->nChargedHadrons.push_back(ijet->chargedHadronMultiplicity());
      jet->nNeutralHadrons.push_back(ijet->neutralHadronMultiplicity());
      jet->nElectrons     .push_back(ijet->electronMultiplicity());
      jet->nPhotons       .push_back(ijet->photonMultiplicity());
      jet->nMuons         .push_back(ijet->muonMultiplicity());
      // cross-check (would fail in HF as of now)
      // why does the code have to crash, when at least
      // one constituent of a jet is in the HF ???
      // as this may happen for jets at eta = 2.4 ...
      //assert(jet->nConstituents.back()== jet->nChargedHadrons.back() + jet->nNeutralHadrons.back() + jet->nElectrons.back() + jet->nPhotons.back() + jet->nMuons.back());

    }
    jet->charge .push_back(ijet->jetCharge());
    jet->flavour.push_back(ijet->partonFlavour());
    jet->bTagCSV.push_back(ijet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

    const reco::SecondaryVertexTagInfo *svTagInfo = ijet->tagInfoSecondaryVertex();

    if(svTagInfo && svTagInfo->nVertices()>0){
	  jet->nSV.push_back(svTagInfo->nVertices());
      jet->SVChi2.push_back(svTagInfo->secondaryVertex(0).chi2());
      jet->SV3DLength.push_back      ( svTagInfo->flightDistance(0).value());
      jet->SV3DLengthError.push_back ( svTagInfo->flightDistance(0).error());
      //use https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideCMSDataAnalysisSchool2015BTaggingExercise#Accessing_additional_b_tag_infor
      //std::vector<const reco::BaseTagInfo*>  baseTagInfos;
      //baseTagInfos.push_back( ijet->tagInfoTrackIP ("impactParameter"));
      //baseTagInfos.push_back( ijet->tagInfoSecondaryVertex("secondaryVertex"));
      //JetTagComputer::TagInfoHelper helper(baseTagInfos);
      //reco::TaggingVariableList vars = computer->taggingVariables(helper);
      TLorentzVector svmom;   
      //if(vars.checkTag(reco::btau::vertexMass)) {
	//const reco::Vertex &vertex = svTagInfo.secondaryVertex(0);
	//svmom.SetPtEtaPhiM(vertex.p4().pt(), vertex.p4().eta(), vertex.p4().phi(), vars.get(reco::btau::vertexMass));
      //}
      jet->SVMomentum.push_back(svmom);
    } else {
      jet->nSV.push_back(0);
      jet->SVChi2.push_back(-1);
      jet->SV3DLength.push_back      (0);
      jet->SV3DLengthError.push_back (0 );
      jet->SVMomentum.push_back(TLorentzVector(0,0,0,0));
    }

    const reco::SecondaryVertexTagInfo &svTagInfo = *ijet->tagInfoSecondaryVertex();
    jet->nSV.push_back(svTagInfo.nVertices());
    if(svTagInfo.nVertices()>0){
      jet->SVChi2.push_back(svTagInfo.secondaryVertex(0).chi2());
      jet->SV3DLength.push_back      ( svTagInfo.flightDistance(0).value());
      jet->SV3DLengthError.push_back ( svTagInfo.flightDistance(0).error());
      std::vector<const reco::BaseTagInfo*>  baseTagInfos;
      baseTagInfos.push_back( ijet->tagInfoTrackIP ("impactParameter"));
      baseTagInfos.push_back( ijet->tagInfoSecondaryVertex("secondaryVertex"));
      JetTagComputer::TagInfoHelper helper(baseTagInfos);
      reco::TaggingVariableList vars = computer->taggingVariables(helper);
      TLorentzVector svmom;   
      if(vars.checkTag(reco::btau::vertexMass)) {
	const reco::Vertex &vertex = svTagInfo.secondaryVertex(0);
	svmom.SetPtEtaPhiM(vertex.p4().pt(), vertex.p4().eta(), vertex.p4().phi(), vars.get(reco::btau::vertexMass));
      }
      jet->SVMomentum.push_back(svmom);
    } else {
      jet->SVChi2.push_back(-1);
      jet->SV3DLength.push_back      (0);
      jet->SV3DLengthError.push_back (0 );
      jet->SVMomentum.push_back(TLorentzVector(0,0,0,0));
    }

    if(hasQGTag  ) jet->gluonTag.push_back(ijet->userFloat(gluonTagName_));
    if(hasJERSF  ) jet->jerSF   .push_back(ijet->userFloat("jerSF"      ));
    if(hasJESSF  ) jet->jesSF   .push_back(ijet->userFloat("jesSF"      ));
    if(hasTotalSF) jet->totalSF .push_back(ijet->userFloat("totalSF"    ));
    if(hasBReg   ) jet->breg    .push_back(ijet->userFloat("BRegResult" ));

    std::pair<TVector2, TVector2> pulls = getPullVector( ijet );
    jet->pull       .push_back(pulls.first );
    jet->pullCharged.push_back(pulls.second);

    std::pair<TVector2, TVector2> genPulls = getGenPullVector( ijet );
    jet->genPull       .push_back(genPulls.first );
    jet->genPullCharged.push_back(genPulls.second);
  }

  if(alternativeJetsAvailable){
    edm::Handle<std::vector<pat::Jet> > alternativeJets;
    evt.getByLabel(alternativeJets_, alternativeJets);

    unsigned short alternativeJetIndex = 0;
    for(const auto &altJet : *alternativeJets) {
      jet->alternativeJet.push_back(TLorentzVector(altJet.px(), altJet.py(), altJet.pz(), altJet.energy()));
      if(alternativeJetIndex == 5) break;
      ++alternativeJetIndex;
    }
  }
  edm::Handle<std::vector<pat::MET> > met;
  evt.getByLabel(met_, met);
  if(met.isValid()){
    jet->met = TLorentzVector(met->at(0).px(), met->at(0).py(), met->at(0).pz(), met->at(0).energy());
    jet->sumEt = met->at(0).sumEt();
  }
  else{
    jet->met = TLorentzVector(0,0,0,0);
    jet->sumEt = -1.;
  }

  trs->Fill();
}

void
JetEventAnalyzer::beginJob()
{
  if( !trs ) throw edm::Exception( edm::errors::Configuration, "TreeRegistryService is not registered in cfg file!" );

  jet = new JetEvent();
  trs->Branch("jet.", jet);
}

void
JetEventAnalyzer::endJob()
{
}

JetEventAnalyzer::~JetEventAnalyzer()
{
}

std::pair<TVector2, TVector2>
JetEventAnalyzer::getPullVector( std::vector<pat::Jet>::const_iterator patJet )
{
  TVector2 null(0,0);

  if (patJet->isPFJet() == false) {
    return std::make_pair(null,null);
  }

  //re-reconstruct the jet direction with the charged tracks
  TLorentzVector chargedJet(0,0,0,0);
  TLorentzVector constituent(0,0,0,0);
  unsigned int nCharged = 0;

  for(size_t idx = 0, length = patJet->numberOfDaughters() ; idx < length; ++idx){
    const reco::Candidate *cand = patJet->daughter(idx);
    if( cand->charge() != 0 ){
      constituent.SetPtEtaPhiE( cand->pt(), cand->eta(), cand->phi(), cand->energy() );
      chargedJet += constituent;
      ++nCharged;
    }
  }

  double jetPt        = patJet   ->pt(), jetPhi        = patJet   ->phi(), jetRapidity        = patJet   ->rapidity();
  double jetPtCharged = chargedJet.Pt(), jetPhiCharged = chargedJet.Phi(), jetRapidityCharged = chargedJet.Rapidity();
  TVector2 r(0,0);
  TVector2 pullAll(0,0);
  TVector2 pullCharged(0,0);

  for(size_t idx = 0, length = patJet->numberOfDaughters() ; idx < length; ++idx){
    const reco::Candidate *cand = patJet->daughter(idx);
    double constituentPt       = cand->pt();
    double constituentPhi      = cand->phi();
    double constituentRapidity = cand->rapidity();
    r.Set( constituentRapidity - jetRapidity, TVector2::Phi_mpi_pi( constituentPhi - jetPhi ) );
    pullAll += ( constituentPt / jetPt ) * r.Mod() * r;
    //calculate TVector using only charged tracks
    if( cand->charge() != 0  )
      r.Set( constituentRapidity - jetRapidityCharged, TVector2::Phi_mpi_pi( constituentPhi - jetPhiCharged ) );
    pullCharged += ( constituentPt / jetPtCharged ) * r.Mod() * r;
  }

  // if there are less than two charged tracks do not calculate the pull (there is not enough info), return null vector
  //TODO really needed???????
  if( nCharged < 2 )
    pullCharged = null;

  return std::make_pair(pullAll, pullCharged);
}

// TODO Sorry for duplicating this function. Maybe there is a nicer option?
std::pair<TVector2, TVector2>
JetEventAnalyzer::getGenPullVector( std::vector<pat::Jet>::const_iterator patJet )
{
  TVector2 null(0,0);

  if ((patJet->genJet()) == false) {
    return std::make_pair(null,null);
  }

  //re-reconstruct the jet direction with the charged tracks
  TLorentzVector chargedJet(0,0,0,0);
  TLorentzVector constituent(0,0,0,0);
  unsigned int nCharged = 0;

  for(size_t idx = 0, length = patJet->genJet()->numberOfDaughters() ; idx < length; ++idx){
    const reco::Candidate *cand = patJet->genJet()->daughter(idx);
    if( cand->charge() != 0 ){
      constituent.SetPtEtaPhiE( cand->pt(), cand->eta(), cand->phi(), cand->energy() );
      chargedJet += constituent;
      ++nCharged;
    }
  }

  double jetPt        = patJet->genJet()->pt(), jetPhi = patJet->genJet()->phi(), jetRapidity = patJet->genJet()->rapidity();
  double jetPtCharged = chargedJet.Pt(), jetPhiCharged = chargedJet.Phi(), jetRapidityCharged = chargedJet.Rapidity();
  TVector2 r(0,0);
  TVector2 pullAll(0,0);
  TVector2 pullCharged(0,0);

  for(size_t idx = 0, length = patJet->genJet()->numberOfDaughters() ; idx < length; ++idx){
    const reco::Candidate *cand = patJet->genJet()->daughter(idx);
    double constituentPt       = cand->pt();
    double constituentPhi      = cand->phi();
    double constituentRapidity = cand->rapidity();
    r.Set( constituentRapidity - jetRapidity, TVector2::Phi_mpi_pi( constituentPhi - jetPhi ) );
    pullAll += ( constituentPt / jetPt ) * r.Mod() * r;
    //calculate TVector using only charged tracks
    if( cand->charge() != 0  )
      r.Set( constituentRapidity - jetRapidityCharged, TVector2::Phi_mpi_pi( constituentPhi - jetPhiCharged ) );
    pullCharged += ( constituentPt / jetPtCharged ) * r.Mod() * r;
  }

  // if there are less than two charged tracks do not calculate the pull (there is not enough info), return null vector
  //TODO really needed???????
  if( nCharged < 2 )
    pullCharged = null;

  return std::make_pair(pullAll, pullCharged);
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetEventAnalyzer);
