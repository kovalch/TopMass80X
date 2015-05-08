#include "FWCore/Utilities/interface/EDMException.h"
#include "TopAnalysis/ZTopUtils/plugins/EventWeightPDF.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <stdlib.h>

/////////////////////////////////////////////////////////////////////////////////

EventWeightPDF::EventWeightPDF(const edm::ParameterSet& cfg):
genEventInfoTag_(cfg.getParameter<edm::InputTag>("genEventInfoTag")), //generator
lheEventInfoTag_(cfg.getParameter<edm::InputTag>("lheEventInfoTag")), //externalLHEProducer
pdfSetNames_(cfg.getParameter<std::vector<std::string> >("PDFSetNames")), 
beginWeightID_(cfg.getParameter<std::vector<std::string> >("beginWeightID")), 
printLHE_(cfg.getParameter<bool>("printLHE"))
{  
    if(pdfSetNames_.size()!=beginWeightID_.size()){
	throw cms::Exception("EventWeightPDF")<<"PDF and WeightID vectors different in size!";
    }
    for(size_t k=0; k<pdfSetNames_.size(); ++k){     
	if(atoi(beginWeightID_.at(k).c_str())%1000!=1) 
	    throw cms::Exception("EventWeightPDF")<<"The array should start at x001, no?";
    }
    for(unsigned int k=0; k<pdfSetNames_.size(); ++k){
	size_t dot = pdfSetNames_.at(k).find_first_of('.');
	size_t underscore = pdfSetNames_.at(k).find_first_of('_');
	if(underscore<dot){
	    pdfShortNames_.push_back(pdfSetNames_.at(k).substr(0,underscore));
	}else{
	    pdfShortNames_.push_back(pdfSetNames_.at(k).substr(0,dot));
	}
	produces<std::vector<double> >(pdfShortNames_.at(k).data());
    }
}

/////////////////////////////////////////////////////////////////////////////////

EventWeightPDF::~EventWeightPDF(){}

/////////////////////////////////////////////////////////////////////////////////

void EventWeightPDF::produce(edm::Event& evt, const edm::EventSetup& setup)
{

    if (evt.isRealData()) return;
    // get final event weight
    edm::Handle<GenEventInfoProduct> evt_info;
    try {evt.getByLabel(genEventInfoTag_, evt_info);}
    catch (...) {;}
    // get systematic variation at LHE level  
    edm::Handle<LHEEventProduct> lhe_info;
    try {evt.getByLabel(lheEventInfoTag_, lhe_info);} 
    catch (...) {;}
    // put PDF weights in the event
    for (size_t k=0; k<pdfSetNames_.size(); ++k) {    
	std::auto_ptr<std::vector<double> > weights (new std::vector<double>);
	int begin(-1);
	int end(-1);
	if (evt_info.isValid()&&lhe_info.isValid()){ 
	    for(size_t iwgt=0; iwgt<lhe_info->weights().size(); ++iwgt){
		if(lhe_info->weights().at(iwgt).id==beginWeightID_.at(k)) 
		    begin=iwgt;
		if(atoi(lhe_info->weights().at(iwgt).id.c_str())-atoi(beginWeightID_.at(k).c_str())<999) 
		    end=iwgt;
	    }
	}
	if(begin>-1 && end>-1 && end>begin){
	    unsigned nPDFvars=end-begin+1;
	    weights->reserve(nPDFvars);
	    for (unsigned i=0; i<nPDFvars; ++i) {	    
		//std::cout<<"Fill in weight at point "<< begin+i << " with id " << lhe_info->weights().at(begin+i).id << std::endl;
		const LHEEventProduct::WGT& wgt = lhe_info->weights().at(begin+i);
		weights->push_back(evt_info->weight()*wgt.wgt/lhe_info->originalXWGTUP()); 	
	    }
	}
	else{
	    double defweight=(evt_info.isValid()) ? evt_info->weight() : 1.;	  
	    edm::LogWarning("EventWeightPDF")<<"Specified weight ID for PDF("<<k<<")not available or screwed up! Use "<<defweight<<" as weight!";
	    weights->reserve(1);
	    weights->push_back(defweight);
	}
	evt.put(weights,pdfShortNames_[k]);
    }
}

/////////////////////////////////////////////////////////////////////////////////

void EventWeightPDF::endRun(edm::Run const& run, const edm::EventSetup& setup)
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
	edm::LogWarning ("EventWeightPDF") << "Can't get LHE header!" ;
    }
}


// Define this as a plug-in                                                                             
DEFINE_FWK_MODULE(EventWeightPDF);              
