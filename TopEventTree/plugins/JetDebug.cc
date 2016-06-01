
//
// Original Author:  Christoph Garbers
//         Created: 17 Oct 2016
//

#include <memory>
#include <iostream>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "AnalysisDataFormats/TopObjects/interface/TtSemiLeptonicEvent.h"

//
// class declaration
//

class JetDebug : public edm::EDProducer {
 public:
  explicit JetDebug(const edm::ParameterSet&);
  ~JetDebug();

 private:
  edm::EDGetTokenT<std::vector<pat::Jet>> jets_;
  edm::EDGetTokenT<std::vector<TtSemiLeptonicEvent>> evtSol_;
  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
};

JetDebug::JetDebug(const edm::ParameterSet& iConfig) {
  jets_ = consumes<std::vector<pat::Jet>>(
      iConfig.getParameter<edm::InputTag>("jets"));
  evtSol_ = mayConsume<std::vector<TtSemiLeptonicEvent>>(
      iConfig.getParameter<edm::InputTag>("evtSolLabel"));
  produces<std::vector<double>>();
}

JetDebug::~JetDebug() {}

//
// member functions
//

// ------------ method called to produce the data  ------------
void JetDebug::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  Handle<std::vector<pat::Jet>> jets;
  iEvent.getByToken(jets_, jets);

  Handle<std::vector<TtSemiLeptonicEvent>> evtSols;
  iEvent.getByToken(evtSol_, evtSols);

  std::auto_ptr<std::vector<double>> fDoutput(new std::vector<double>);
  unsigned int kJetMAX_ = 4;
  unsigned short jetIndex = 0;
  for (std::vector<pat::Jet>::const_iterator ijet = jets->begin();
       ijet != jets->end(); ++ijet, ++jetIndex) {
    if (jetIndex == kJetMAX_) break;
    if (ijet == jets->begin()) {
      std::cout << "Flavour outOf FlavourInfo "
                << ijet->jetFlavourInfo().getPartonFlavour() << std::endl;
      /*std::cout<<"JetEnergy ueber p4 "<<ijet->p4().E()<<std::endl;
      std::vector<std::string> sets = ijet->availableJECSets();
      for(unsigned int i = 0; i<sets.size();i++){
              std::cout<<sets.at(i)<<std::endl;
              std::vector<std::string> labels =ijet->availableJECLevels(i);
              for(unsigned int j = 0; j<labels.size();j++){
                      std::cout<<"	"<<labels.at(j)<<std::endl;
                      std::cout<<"		"<<ijet->jecFactor(j,
      ijet->currentJECFlavor(), i)<<std::endl;
                      std::cout<<"		correctedJet "<<ijet->correctedJet(j,
      ijet->currentJECFlavor(), i)<<std::endl;
                      //std::cout<<"		correctedJet b
      "<<ijet->correctedJet(labels.at(j), "BOTTOM", sets.at(i))<<std::endl;
                      //std::cout<<"		correctedJet light
      "<<ijet->correctedJet(labels.at(j), "UDS", sets.at(i))<<std::endl;
              }
      }
      */
      std::cout << std::endl;
    }
    fDoutput->push_back(ijet->energy());
  }

  /*std::vector<TtSemiLeptonicEvent>::const_iterator iEvtSol = evtSols->begin();
  std::cout<<"TtSemiEvtSol: "<<std::endl;
  std::cout<<iEvtSol->hadronicDecayB("ttSemiLepHypHitFit") <<std::endl;
  std::cout<<iEvtSol->hadronicDecayB("ttSemiLepHypHitFit")->energy()
  <<std::endl;*/

  iEvent.put(fDoutput);
}

// ------------ method called once each job just before starting event loop
// ------------
void JetDebug::beginJob() {}

// ------------ method called once each job just after ending the event loop
// ------------
void JetDebug::endJob() {}

// define this as a plug-in
DEFINE_FWK_MODULE(JetDebug);
