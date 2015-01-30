#include "FWCore/Utilities/interface/EDMException.h"
#include "TopAnalysis/ZTopUtils/plugins/EventWeightMCSystematic.h"
//#include "CommonTools/UtilAlgos/interface/TFileService.h"
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
    edm::InputTag lheEventInfoTag("source");
    edm::Handle<LHEEventProduct> lhe_info;
    evt.getByLabel(lheEventInfoTag, lhe_info);
    weight_=evt_info->weight();
    // weightid=0 gives default MC weight, same in case id is not available for this MC  
    if(weightId!="0"){
	for (size_t iwgt=0; iwgt<lhe_info->weights().size(); ++iwgt) {
	    const LHEEventProduct::WGT& wgt = lhe_info->weights().at(iwgt);
	    std::cout << "agrohsje weight id is " << wgt.id <<" with value  "<< wgt.wgt << std::endl;    
	    if(wgt.id==weightId)
		weight_=evt_info->weight()*wgt.wgt/lhe_info->originalXWGTUP(); 
	}
    }
    std::auto_ptr<double> eventWeight(new double);
    *eventWeight = weight_ ;
    evt.put(eventWeight);
}

// Define this as a plug-in                                                                             
DEFINE_FWK_MODULE(EventWeightMCSystematic);              
