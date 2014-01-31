#include "FWCore/Utilities/interface/EDMException.h"
#include "Boolcombiner.h"


Boolcombiner::Boolcombiner(const edm::ParameterSet& cfg):
debug(cfg.getParameter< bool >("debug") ),
mustbetrue_ ( cfg.getParameter<std::vector<edm::InputTag> >("mustBeTrue") ),
mustbefalse_( cfg.getParameter<std::vector<edm::InputTag> >("mustBeFalse") )
{
	produces<bool>();
	if(debug){
		std::cout << "Boolcombiner: " << std::endl;
		std::cout << "will require the following to be true: " <<std::endl;
		for(size_t i=0;i<mustbetrue_.size();i++){
			std::cout << mustbetrue_.at(i) << ' ';
		}
		std::cout << std::endl;
		std::cout << "will require the following to be false: " <<std::endl;
		for(size_t i=0;i<mustbefalse_.size();i++){
			std::cout << mustbefalse_.at(i) << ' ' ;
		}
		std::cout << std::endl;
	}

}

void
Boolcombiner::produce(edm::Event& evt, const edm::EventSetup& setup)
{
	//---------------------------------------------
	// get event weight
	//---------------------------------------------
	//edm::Handle<double> weight;
	//evt.getByLabel(weight_, weight);

	std::auto_ptr<bool> combined(new bool(true));
	edm::Handle<bool> filterbool;





	for(size_t i=0;i<mustbefalse_.size();i++){
		evt.getByLabel(mustbefalse_.at(i), filterbool);
		if(debug) std::cout << *filterbool << ' ';
		if(*filterbool){
			*combined=false;
			if(!debug) break;
		}
	}
	if(debug) std::cout << std::endl;
	for(size_t i=0;i<mustbetrue_.size();i++){
		evt.getByLabel(mustbetrue_.at(i), filterbool);
		if(debug) std::cout << *filterbool << ' ';
		if(! *filterbool){
			*combined=false;
			if(!debug) break;
		}
	}

	if(debug){
		std::cout << std::endl;
		std::cout << (*combined);
		std::cout << std::endl;}
	evt.put(combined);

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( Boolcombiner );
