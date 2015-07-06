#include "FWCore/Utilities/interface/EDMException.h"
#include "TopAnalysis/ZTopUtils/plugins/EventWeightMCSystematic.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "FWCore/Framework/interface/MakerMacros.h"

/////////////////////////////////////////////////////////////////////////////////

EventWeightMCSystematic::EventWeightMCSystematic(const edm::ParameterSet& cfg):
genEventInfoTag_(cfg.getParameter<edm::InputTag>("genEventInfoTag")), //generator
lheEventInfoTag_(cfg.getParameter<edm::InputTag>("lheEventInfoTag")), //externalLHEProducer
weightID_(cfg.getParameter<std::string>("weightID")), 
printLHE_(cfg.getParameter<bool>("printLHE"))
{  

    produces<double>();
}

/////////////////////////////////////////////////////////////////////////////////

EventWeightMCSystematic::~EventWeightMCSystematic(){}

/////////////////////////////////////////////////////////////////////////////////

void EventWeightMCSystematic::produce(edm::Event& evt, const edm::EventSetup& setup)
{
    if(evt.isRealData()) return;
    // get nominal event weight
    edm::Handle<GenEventInfoProduct> evt_info;
    try {evt.getByLabel(genEventInfoTag_, evt_info);}
    catch (...) {;}
    // get systematic variation at LHE level  
    edm::Handle<LHEEventProduct> lhe_info;
    try {evt.getByLabel(lheEventInfoTag_, lhe_info);} 
    catch (...) {;}
    double weight(1.0);
    bool allinfo(false);
    bool foundid(false);	
    if(weightID_=="0000"){
	if(evt_info.isValid()){
	    weight=evt_info->weight(); 
	    allinfo=true;
	    //std::cout<<" fill in nominal weight " << weight << std::endl; 
	}
    }else{ 
	if(lhe_info.isValid()){ 
	    allinfo=true;
	    for(size_t iwgt=0; iwgt<lhe_info->weights().size(); ++iwgt){
		const LHEEventProduct::WGT& wgt = lhe_info->weights().at(iwgt);
		if(wgt.id==weightID_){
		    weight=wgt.wgt/lhe_info->originalXWGTUP(); 
		    //std::cout<<" fill in weight with ID=="<< weightID_ <<" : " << weight << std::endl; 
		    foundid=true;
		    break;
		}	
	    }
	}
    }
    
    if(!allinfo)
	edm::LogWarning("EventWeightMCSystematic")<<"Can't get weight information! Take "<<weight<<" as weight!";
    
    if(weightID_!="0000"&&!foundid)
	edm::LogWarning("EventWeightMCSystematic")<<"Specified weight ID not available! Take "<<weight<<" as weight!";
    
    std::auto_ptr<double> eventWeight(new double);
    *eventWeight = weight ;
    evt.put(eventWeight);
}

/////////////////////////////////////////////////////////////////////////////////
#ifndef CMSSW_LEQ_5
void EventWeightMCSystematic::endRun(edm::Run const& run, const edm::EventSetup& setup)
{
    if(!printLHE_) return;
    edm::Handle<LHERunInfoProduct> lheRunInfo;
    typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
    try {run.getByLabel(lheEventInfoTag_, lheRunInfo);}
    catch (...) {;}
    if(lheRunInfo.isValid()){
	LHERunInfoProduct lheRunInfoProd = *(lheRunInfo.product());
	for (headers_const_iterator iter=lheRunInfoProd.headers_begin(); iter!=lheRunInfoProd.headers_end(); iter++){
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
#endif
// Define this as a plug-in                                                                             
DEFINE_FWK_MODULE(EventWeightMCSystematic);              
