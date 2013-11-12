#include <algorithm>

#include "TopMass/TopEventTree/plugins/BRegProducer.h"



#include "CommonTools/Utils/interface/PtComparator.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


BRegProducer::BRegProducer(const edm::ParameterSet& cfg):
  inputJets_           (cfg.getParameter<edm::InputTag>("jets"           )),
  kJetMAX_(cfg.getParameter<int>("maxNJets"))
{
	if(inputJets_.instance()==""){
		outputJets_ = inputJets_.label();// + "DummyBReg";
	}
	else outputJets_ = inputJets_.instance();// + "DummyBReg";
	std::cout << "inputJets_.label() " << inputJets_.label() << " inputJets_.instance(): " << inputJets_.instance() << " outputJets; " << outputJets_.c_str() << std::endl;

	bRegAnalyzer = new BRegJetEventAnalyzer(cfg);
	std::cout << "defined a BRegJetEventAnalyzer object (with all necessary configs)" << std::endl;


	produces<std::vector<pat::Jet> >(outputJets_);

}

void
BRegProducer::beginJob()
{ 
	bRegAnalyzer->iniBRegEvent();
}

///function to sort any auto_ptr to a vector of anything which has a pt function
template<typename T>
void sortByPt(std::auto_ptr<std::vector<T> > &collection) {
  std::sort(collection->begin(), collection->end(), GreaterByPt<T>());
}

void
BRegProducer::produce(edm::Event& event, const edm::EventSetup& setup)
{


  // access jets
	edm::Handle<edm::View<pat::Jet> > jets;
//  edm::Handle<std::vector<pat::Jet> > jets;
  event.getByLabel(inputJets_, jets);

//  std::cout << "inputJets_.label() " << inputJets_.label() << " inputJets_.instance() " << inputJets_.instance() << " inputJets_.process() " << inputJets_.process() << std::endl;

  std::auto_ptr<std::vector<pat::Jet> > pJets(new std::vector<pat::Jet>);

    std::vector<double> tempBRegGBRTrainResultCrossCheck = 	bRegAnalyzer->fillBRegJetAndReturnGBRResults(event,  setup);
//  std::vector<double> tempBRegGBRTrainResultCrossCheck = 	bRegAnalyzer->fillBRegJetAndReturnTMVAResults(event,  setup);

unsigned int jet_i =0;
for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){

//  for(std::vector<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
    pat::Jet scaledJet = *jet;




//std::cout << "jet_i" << jet_i << "kJetMAX_" << kJetMAX_ << "tempBRegGBRTrainResultCrossCheck.size() " << tempBRegGBRTrainResultCrossCheck.size() << std::endl;

    if(jet_i>=(unsigned int) kJetMAX_ || jet_i>=tempBRegGBRTrainResultCrossCheck.size()){
    	if((unsigned int) kJetMAX_!=tempBRegGBRTrainResultCrossCheck.size()){
    	edm::LogWarning msg("Max. no of jets");
	    msg << "The maximum number of jets has been reached. This means, that there are no more BReg-Results available to add to the jets. The maximum number of jets does not match the size of the results vector. Check inputs!\n";
    	}
//	    throw cms::Exception("Configuration Error");
    	continue;
    }
//	std::cout <<   tempBRegGBRTrainResultCrossCheck.at(jet_i) << std::endl;
	scaledJet.addUserFloat("BRegResult" , tempBRegGBRTrainResultCrossCheck.at(jet_i));
//	std::cout << "trying to add ref to jet before" << std::endl;
	edm::RefToBase<pat::Jet> refToJetWithValueMaps =  jets->refAt(jet_i);
//	std::cout << "adduserdatatrying to add ref to jet before" << std::endl;
	scaledJet.addUserData< edm::RefToBase<pat::Jet> >( "refToJetWithValueMaps", refToJetWithValueMaps );
//	std::cout << "did it adduserdatatrying to add ref to jet before" << std::endl;

    pJets->push_back( scaledJet );
    jet_i++;
  }
//
//  //p4 changes might have changed the pt order, so need to sort the new collections
//  sortByPt(pJets);

//std::cout << "------------------------------ pJets.size() "<< pJets->size() << "--------------------------------------" << std::endl;
//std::cout << "------------------------------trying to put new jets " << outputJets_ << " into event----------------------------" << std::endl;
event.put(pJets, outputJets_);

}


