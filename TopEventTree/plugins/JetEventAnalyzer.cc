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

#include "TopMass/TopEventTree/interface/TreeRegistryService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "TopMass/TopEventTree/plugins/JetEventAnalyzer.h"

JetEventAnalyzer::JetEventAnalyzer(const edm::ParameterSet& cfg):
jets_        (cfg.getParameter<edm::InputTag>("jets")),
//allJets_     (cfg.getParameter<edm::InputTag>("allJets")),
//noPtEtaJets_ (cfg.getParameter<edm::InputTag>("noPtEtaJets")),
gluonTagSrc_ (cfg.getParameter<edm::InputTag>("gluonTagSrc")),
kJetMAX_(cfg.getParameter<int>("maxNJets")),
checkedJERSF(false), checkedJESSF(false), checkedTotalSF(false), checkedQGTag(false),
    hasJERSF(false),     hasJESSF(false),     hasTotalSF(false),     hasQGTag(false)
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

  edm::Handle<edm::ValueMap<float> >  QGTagsHandleLikelihood;
  evt.getByLabel(gluonTagSrc_, QGTagsHandleLikelihood);

  unsigned short jetIndex = 0;
  for(std::vector< pat::Jet >::const_iterator ijet = jets->begin(); ijet != jets->end(); ++ijet, ++jetIndex) {
    // write only kJetMAX_ jets into the event
    if(jetIndex == kJetMAX_) break;

    jet->jet.push_back(TLorentzVector(ijet->px(), ijet->py(), ijet->pz(), ijet->energy()));

    jet->charge .push_back(ijet->jetCharge());
    jet->flavour.push_back(ijet->partonFlavour());
    jet->bTagCSV.push_back(ijet->bDiscriminator("combinedSecondaryVertexBJetTags"));

    // check only once per module run if the needed collections are available
    if(!checkedQGTag  ) { checkedQGTag   = true; if(QGTagsHandleLikelihood.isValid()) hasQGTag   = true; }
    if(!checkedJERSF  ) { checkedJERSF   = true; if(ijet->hasUserFloat("jerSF"     )) hasJERSF   = true; }
    if(!checkedJESSF  ) { checkedJESSF   = true; if(ijet->hasUserFloat("jesSF"     )) hasJESSF   = true; }
    if(!checkedTotalSF) { checkedTotalSF = true; if(ijet->hasUserFloat("totalSF"   )) hasTotalSF = true; }

    if(hasQGTag){
      edm::RefToBase<pat::Jet> jetRef(edm::Ref<std::vector<pat::Jet> >(jets, jetIndex));
      jet->gluonTag.push_back((*QGTagsHandleLikelihood)[jetRef]);
    }
    if(hasJERSF  ) jet->jerSF  .push_back(ijet->userFloat("jerSF"  ));
    if(hasJESSF  ) jet->jesSF  .push_back(ijet->userFloat("jesSF"  ));
    if(hasTotalSF) jet->totalSF.push_back(ijet->userFloat("totalSF"));
  }

  trs->Fill();
}

void
JetEventAnalyzer::beginJob()
{
  if( !trs ) throw edm::Exception( edm::errors::Configuration, "TreeRegistryService is not registered in cfg file!" );

  jet = new JetEvent();
  trs->Branch("jet", jet);
}

void
JetEventAnalyzer::endJob()
{
}

JetEventAnalyzer::~JetEventAnalyzer()
{
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetEventAnalyzer);
