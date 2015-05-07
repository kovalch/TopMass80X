#include "FWCore/Utilities/interface/EDMException.h"
#include "TopAnalysis/ZTopUtils/plugins/EventWeightMCSystematic.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "FWCore/Framework/interface/MakerMacros.h"

/////////////////////////////////////////////////////////////////////////////////

EventWeightMCSystematic::EventWeightMCSystematic(const edm::ParameterSet& cfg):
    parameterSet_(cfg) 
{  
    produces<double>();
}

/////////////////////////////////////////////////////////////////////////////////

EventWeightMCSystematic::~EventWeightMCSystematic(){}

/////////////////////////////////////////////////////////////////////////////////

void EventWeightMCSystematic::produce(edm::Event& evt, const edm::EventSetup& setup)
{
    // get final event weight
    std::string weightId(parameterSet_.getParameter<std::string>("weightId"));
    edm::InputTag genEventInfoTag("generator");
    edm::Handle<GenEventInfoProduct> evt_info;
    try {evt.getByLabel(genEventInfoTag, evt_info);}
    catch (...) {;}
    // get systematic variation at LHE level  
    edm::InputTag lheEventInfoTag("externalLHEProducer");
    edm::Handle<LHEEventProduct> lhe_info;
    try {evt.getByLabel(lheEventInfoTag, lhe_info);} 
    catch (...) {;}
    double weight(1.0);        
    // calculate new weight 
    if ( evt_info.isValid() && lhe_info.isValid() ) { 
	weight = evt_info->weight();
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
	    edm::LogWarning ("EventWeightMCSystematic") << "Specified weight ID not available! Take default MC weight!" ;
    }else{
	edm::LogWarning ("EventWeightMCSystematic") << "Can't get weight information! Take 1 as weight!" ;
    }
    std::auto_ptr<double> eventWeight(new double);
    *eventWeight = weight ;
    evt.put(eventWeight);
}

/////////////////////////////////////////////////////////////////////////////////

void EventWeightMCSystematic::endRun(edm::Run const& run, const edm::EventSetup& setup)
{
    bool printLHE(parameterSet_.getParameter<bool>("printLHE"));
    if(!printLHE) return;
    edm::Handle<LHERunInfoProduct> lheRunInfoTag;
    typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
    try {run.getByLabel("externalLHEProducer", lheRunInfoTag);}
    catch (...) {;}
    if(lheRunInfoTag.isValid()){
	LHERunInfoProduct lheRunInfo = *(lheRunInfoTag.product());
	for (headers_const_iterator iter=lheRunInfo.headers_begin(); iter!=lheRunInfo.headers_end(); iter++){
	    std::cout << iter->tag() << std::endl;
	    std::vector<std::string> lines = iter->lines();
	    for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
		std::cout << lines.at(iLine);
	    }
	}
    }else{
	edm::LogWarning ("EventWeightMCSystematic") << "Can't get LHE header!" ;
    }
}

// Define this as a plug-in                                                                             
DEFINE_FWK_MODULE(EventWeightMCSystematic);              
