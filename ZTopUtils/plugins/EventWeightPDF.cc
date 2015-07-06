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
lheEventInfoTag_(cfg.getParameter<edm::InputTag>("lheEventInfoTag")), //externalLHEProducer
PDFName_(cfg.getParameter<std::vector<std::string> >("PDFName")), 
nominalWeightID_(cfg.getParameter<std::vector<std::string> >("nominalWeightID")), 
beginWeightID_(cfg.getParameter<std::vector<std::string> >("beginWeightID")), 
endWeightID_(cfg.getParameter<std::vector<std::string> >("endWeightID")), 
printLHE_(cfg.getParameter<bool>("printLHE"))
{  
    if(PDFName_.size()!=nominalWeightID_.size() ||
       PDFName_.size()!=beginWeightID_.size() ||
       PDFName_.size()!=endWeightID_.size()){
	throw cms::Exception("EventWeightPDF")<<"PDFName, nominalWeightID, beginWeightID, endWeightID vectors different in size!";
    }

    for(unsigned int k=0; k<PDFName_.size(); ++k){
	size_t dot = PDFName_.at(k).find_first_of('.');
	size_t underscore = PDFName_.at(k).find_first_of('_');
	if(underscore<dot){
	    shortPDFName_.push_back(PDFName_.at(k).substr(0,underscore));
	}else{
	    shortPDFName_.push_back(PDFName_.at(k).substr(0,dot));
	}
	produces<std::vector<double> >(shortPDFName_.at(k).data());
    }
}

/////////////////////////////////////////////////////////////////////////////////

EventWeightPDF::~EventWeightPDF(){}

/////////////////////////////////////////////////////////////////////////////////

void EventWeightPDF::produce(edm::Event& evt, const edm::EventSetup& setup)
{

    if (evt.isRealData()) return;
    // get systematic variation at LHE level  
    edm::Handle<LHEEventProduct> lhe_info;
    try {evt.getByLabel(lheEventInfoTag_, lhe_info);} 
    catch (...) {;}
    // put PDF weights in the event
    for (size_t k=0; k<PDFName_.size(); ++k) {    
	std::auto_ptr<std::vector<double> > weights (new std::vector<double>);
	int nominal(-1);
	int begin(-1);
	int end(-1);
	if(lhe_info.isValid()){ 
	    for(size_t iwgt=0; iwgt<lhe_info->weights().size(); ++iwgt){
		if(lhe_info->weights().at(iwgt).id==nominalWeightID_.at(k)) 
		    nominal=iwgt;
		if(lhe_info->weights().at(iwgt).id==beginWeightID_.at(k))
		    begin=iwgt;
		if(lhe_info->weights().at(iwgt).id==endWeightID_.at(k))
		    end=iwgt;
		if(nominal>-1 && begin>-1 && end>-1) break;
	    }
	}
	//atoi(lhe_info->weights().at(iwgt).id.c_str())-atoi(beginWeightID_.at(k).c_str())<999) 
	if(nominal>-1 && begin>-1 && end>-1 && end>begin){
	    unsigned nPDFvars=end-begin+1;
	    weights->reserve(nPDFvars+1);
	    weights->push_back(lhe_info->weights().at(nominal).wgt/lhe_info->originalXWGTUP()); 	
	    //std::cout<<"Fill in weight at point "<< nominal << " with id " << lhe_info->weights().at(nominal).id << " "<< 
	    //lhe_info->weights().at(nominal).wgt/lhe_info->originalXWGTUP()<< std::endl;
	    for (unsigned i=0; i<nPDFvars; ++i) {	    
		const LHEEventProduct::WGT& wgt = lhe_info->weights().at(begin+i);
		weights->push_back(wgt.wgt/lhe_info->originalXWGTUP()); 	
		//std::cout<<"Fill in weight at point "<< begin+i << " with id " << lhe_info->weights().at(begin+i).id << " "<< 
		//wgt.wgt/lhe_info->originalXWGTUP()<<std::endl;
		
	    }
	}
	else{
	    edm::LogWarning("EventWeightPDF")<<"Specified weight ID for PDF("<<k<<")not available or screwed up! Use 1.0 as weight!";
	    weights->reserve(1);
	    weights->push_back(1.0);
	}
	evt.put(weights,shortPDFName_[k]);
    }
}

/////////////////////////////////////////////////////////////////////////////////
#ifndef CMSSW_LEQ_5
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
#endif

// Define this as a plug-in                                                                             
DEFINE_FWK_MODULE(EventWeightPDF);              
