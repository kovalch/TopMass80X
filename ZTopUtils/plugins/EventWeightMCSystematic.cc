#include "FWCore/Utilities/interface/EDMException.h"
#include "TopAnalysis/ZTopUtils/plugins/EventWeightMCSystematic.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "FWCore/Framework/interface/MakerMacros.h"

EventWeightMCSystematic::EventWeightMCSystematic(const edm::ParameterSet& cfg):
    parameterSet_(cfg) 
{  
    produces<double>();
}

EventWeightMCSystematic::~EventWeightMCSystematic(){}

void EventWeightMCSystematic::produce(edm::Event& evt, const edm::EventSetup& setup)
{
    // get final event weight
    std::string weightId(parameterSet_.getParameter<std::string>("weightId"));
    edm::InputTag genEventInfoTag("generator");
    edm::Handle<GenEventInfoProduct> evt_info;
    evt.getByLabel(genEventInfoTag, evt_info);

    // get systematic variation at LHE level  
    edm::InputTag lheEventInfoTag("externalLHEProducer");
    edm::Handle<LHEEventProduct> lhe_info;
    evt.getByLabel(lheEventInfoTag, lhe_info);

    double weight(evt_info->weight());
    bool foundid(false);	
    for (size_t iwgt=0; iwgt<lhe_info->weights().size(); ++iwgt) {
	const LHEEventProduct::WGT& wgt = lhe_info->weights().at(iwgt);
	if(wgt.id==weightId){
	    weight=evt_info->weight()*wgt.wgt/lhe_info->originalXWGTUP(); 
	    foundid=true;
	    break;
	}	
    }
    if (!foundid)
	edm::LogWarning ("EventWeightMCSystematic") << "##### Specified weight ID not available! Take default MC weight!#####" ;
    
    std::auto_ptr<double> eventWeight(new double);
    *eventWeight = weight ;
    evt.put(eventWeight);
}

// Define this as a plug-in                                                                             
DEFINE_FWK_MODULE(EventWeightMCSystematic);              
