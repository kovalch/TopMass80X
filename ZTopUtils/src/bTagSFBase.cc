/*
 * bTagSFBase.cc
 *
 *  Created on: Apr 7, 2015
 *      Author: kiesej
 */

#include "../interface/bTagSFBase.h"
#include <cmath> //std::abs
#include "TRandom3.h"


namespace ztop{




bool bTagSFBase::debug=false;




bTagSFBase::bTagSFBase():bTagEfficiency(), calib_(0),fullcalibration_(0),bccalibration_(0), udsgcalibration_(0), syst_(nominal),wpval_(0){

}
bTagSFBase::~bTagSFBase(){
	clear();
}

void bTagSFBase::loadSF(const std::string& filename, BTagEntry::OperatingPoint op,
		const std::string& tagger, std::string measurementType,std::string systsourceup, std::string  systsourcedown){

	if(syst_!=nominal && (systsourceup.length()<1 || systsourcedown.length()< 1)){
		throw std::logic_error("bTagSFBase::loadSF: systematic variation required without loading systematic variations for scale factors");
	}

	loadSFToCalib(&fullcalibration_,filename,op,tagger,measurementType,systsourceup,systsourcedown);

	//fast check
	try{
		jetSF(BTagEntry::FLAV_B,40, 0,0.5);
	}
	catch(std::exception& e){
		throw std::runtime_error("bTagSFBase::loadSF: payload misses b-SF");
	}
	try{
		jetSF(BTagEntry::FLAV_C,40, 0,0.5);
	}
	catch(...){
		throw std::runtime_error("bTagSFBase::loadSF: payload misses c-SF");
	}
	try{
		jetSF(BTagEntry::FLAV_UDSG,40, 0,0.5);
	}
	catch(...){
		throw std::runtime_error("bTagSFBase::loadSF: payload misses light-SF");
	}

	if(debug)
		std::cout << "bTagSFBase::loadSF: all flavours available" << std::endl;

}

void bTagSFBase::loadBCSF(const std::string& filename, BTagEntry::OperatingPoint op,
		const std::string& tagger, std::string measurementType,std::string systsourceup, std::string  systsourcedown){

	if(syst_!=nominal && (systsourceup.length()<1 || systsourcedown.length()< 1)){
		throw std::logic_error("bTagSFBase::loadSF: systematic variation required without loading systematic variations for scale factors");
	}

	loadSFToCalib(&bccalibration_,filename,op,tagger,measurementType,systsourceup,systsourcedown);

	//fast check
	try{
		jetSF(BTagEntry::FLAV_B,40, 0,0.5);
	}
	catch(...){
		throw std::runtime_error("bTagSFBase::loadBCSF: payload misses b-SF");
	}
	try{
		jetSF(BTagEntry::FLAV_C,40, 0,0.5);
	}
	catch(...){
		throw std::runtime_error("bTagSFBase::loadBCSF: payload misses c-SF");
	}

}

void bTagSFBase::loadUDSGSF(const std::string& filename, BTagEntry::OperatingPoint op,
		const std::string& tagger, std::string measurementType,std::string systsourceup, std::string  systsourcedown){

	if(syst_!=nominal && (systsourceup.length()<1 || systsourcedown.length()< 1)){
		throw std::logic_error("bTagSFBase::loadSF: systematic variation required without loading systematic variations for scale factors");
	}

	loadSFToCalib(&udsgcalibration_,filename,op,tagger,measurementType,systsourceup,systsourcedown);

	//fast check
	try{
		jetSF(BTagEntry::FLAV_UDSG,40, 0,0.5);
	}
	catch(...){
		throw std::runtime_error("bTagSFBase::loadSF: payload misses light-SF");
	}

}






void bTagSFBase::setSystematics(systematics sys){
	bool sysfullavail=fullcalibration_ && (fullcalibration_->up() && fullcalibration_->down());
	bool sysbcaval=bccalibration_ && (bccalibration_->up() && bccalibration_->down());
	bool sysudsgaval=udsgcalibration_ && (udsgcalibration_->up() && udsgcalibration_->down());
	bool sysavail=sysfullavail || (sysbcaval && sysudsgaval) ||
			(sysfullavail && sysbcaval)||
			(sysfullavail && sysudsgaval);

	if(sys!=nominal && !sysavail){
		throw std::logic_error("bTagSFBase::setSystematics: systematic variation required without loading systematic variations for scale factors");
	}
	syst_=sys;
}

float bTagSFBase::getWPDiscrValue(const std::string& tagger, BTagEntry::OperatingPoint op)const{
	if(tagger == "csv"){
		if(op == BTagEntry::OP_LOOSE)
			return 0.244;
		else if(op == BTagEntry::OP_MEDIUM)
			return 0.679;
		else if(op == BTagEntry::OP_TIGHT)
			return 0.898;
	}
	//pre 50 ns  csvv2 taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagging for TTH synchronisation
	else if(tagger == "csvv2"){
		if(op == BTagEntry::OP_LOOSE)
			return 0.423;
		else if(op == BTagEntry::OP_MEDIUM)
			return 0.814;
		else if(op == BTagEntry::OP_TIGHT)
			return 0.941;
	}
        //50 ns csvv2 btagging working points
        else if(tagger == "csvv2_50ns"){
                if(op == BTagEntry::OP_LOOSE)
                        return 0.605;
                else if(op == BTagEntry::OP_MEDIUM)
                        return 0.890;
                else if(op == BTagEntry::OP_TIGHT)
                        return 0.970;
        }

	else{
		std::string errstr="bTagSFBase::getWPDiscrValue: tagger string unknown\navailable:\n";
		errstr+="csv\ncsvv1\n";
		throw std::runtime_error(errstr);
	}
	return 0;//shape reweighting mode
}
std::string bTagSFBase::getWorkingPointString()const{
	return wpstring_;
}



bool bTagSFBase::jetIsTagged(const float & pt, const float& eta, const int & genPartonFlavor,
		const float & tagValue, const unsigned int & seed) const{


	BTagEntry::JetFlavor jettype=jetFlavor(genPartonFlavor);

	const bool isBTagged = tagValue > wpval_;
	if (genPartonFlavor == 0)
		return isBTagged;

	if(makeeffs_) return isBTagged;
	const double Btag_eff = jetEff(pt, fabs(eta), jettype);
	const double Btag_SF = jetSF(jettype,pt, eta,tagValue);
	const float coin = TRandom3(seed).Uniform(1.0);

	if ( Btag_SF > 1. ) {  // use this if SF>1
		if ( !isBTagged ) {
			// fraction of jets that need to be upgraded
			const float mistagPercent = (1.0 - Btag_SF) / (1.0 - (1.0/Btag_eff) );

			//upgrade to tagged
			if( coin < mistagPercent ) return true;
		}
	}
	else if ( Btag_SF < 1. ) {  // use this if SF<1
		// downgrade tagged to untagged
		if ( isBTagged && coin > Btag_SF ) return false;
	}
	else {  // no change if exactly SF==1
		return isBTagged;
	}
	return isBTagged;

}


float bTagSFBase::getJetDiscrShapeWeight(const float & pt, const float& eta,
		const int & genPartonFlavor, const float& jetdiscr)const{
	if(debug) std::cout << "bTagSFBase::getJetDiscrShapeWeight " << std::endl;
	if(makeeffs_) return 1.;
	if (genPartonFlavor == 0)
		return 1;

	return jetSF(jetFlavor(genPartonFlavor),pt,eta,jetdiscr);

}




void bTagSFBase::clear(){

	if(calib_) delete calib_;
	calib_=0;
	if(fullcalibration_) delete fullcalibration_;
	fullcalibration_=0;
	if(bccalibration_) delete bccalibration_;
	bccalibration_=0;
	if(udsgcalibration_) delete udsgcalibration_;
	udsgcalibration_=0;

	wpval_=0;
	wpstring_="";
	syst_=nominal;
	filename_="";
}





void bTagSFBase::loadSFToCalib(calibrations** cal, const std::string& filename, BTagEntry::OperatingPoint op,
		const std::string& tagger, const std::string measurementType="",
		const std::string systsourceup="", const std::string  systsourcedown=""){

	if(*cal) delete *cal;
	*cal=0; //except safe
	if(fullcalibration_ ){///set up once
		//check for consist
		if(filename_ != filename){
			throw std::logic_error("bTagSFBase::loadSFToCalib: Do not use calibrations from different files");
		}
	}
	else{
		filename_= filename;
		try{
			calib_ = new BTagCalibration(tagger,filename);
		}
		catch(...){
			calib_=0;
			throw std::runtime_error("bTagSFBase::loadSFToCalib: problem loading payload");
		}
	}
	try{
		*cal = new calibrations(calib_,op,measurementType,systsourceup,systsourcedown);
	}
	catch(...){
		*cal=0;
		throw std::runtime_error("bTagSFBase::loadSFToCalib: problem getting SFs");
	}
	std::string tmpwpstr;
	if(op == BTagEntry::OP_LOOSE)
		tmpwpstr="loose";
	else if(op == BTagEntry::OP_MEDIUM)
		tmpwpstr="medium";
	else if(op == BTagEntry::OP_TIGHT)
		tmpwpstr="tight";
	else
		tmpwpstr="shape";

	if(wpstring_.length() > 0 &&tmpwpstr != wpstring_){
		throw std::logic_error("bTagSFBase::loadSFToCalib: working points not identical for different flavours!");
	}else{
		wpstring_=tmpwpstr;
	}

	float wpvaltmp=getWPDiscrValue(tagger,op);

	if(wpval_ && wpval_!=wpvaltmp){
		throw std::logic_error("bTagSFBase::loadSFToCalib: taggers not identical for different flavours!");
	}else{
		wpval_=wpvaltmp;
	}

	if(debug)
		std::cout << "bTagSFBase::loadSFToCalib: set up SF for "<< tagger <<", "<< measurementType << std::endl;

}






bTagSFBase::calibrations::calibrations()
:up_(0),down_(0),central_(0)
{}
bTagSFBase::calibrations::calibrations(BTagCalibration *c, BTagEntry::OperatingPoint op,
		std::string measurementType,std::string systsourceup, std::string  systsourcedown)
:up_(0),down_(0),central_(0)
{
	central_=new BTagCalibrationReader(c,op,measurementType,"central");
	if(systsourceup.length()>0)
		up_=new BTagCalibrationReader(c,op,measurementType,systsourceup);
	if(systsourcedown.length()>0)
		down_=new BTagCalibrationReader(c,op,measurementType,systsourcedown);
}
bTagSFBase::calibrations::~calibrations(){
	if(central_) delete central_;
	if(up_) delete up_;
	if(down_) delete down_;
}



float bTagSFBase::jetSF(const BTagEntry::JetFlavor& jetflav, const float& pt, const float& eta, const float& discr)const{

	const BTagCalibrationReader * usecalib=0;
	const calibrations * flavcalib=fullcalibration_;
	float abseta=std::abs(eta);
	if(jetflav == BTagEntry::FLAV_B || jetflav == BTagEntry::FLAV_C){
		if(bccalibration_)  flavcalib=bccalibration_;
		if(!flavcalib)
			throw std::logic_error("bTagSFBase::jetSF: calibration not loaded");

		if (syst_ == heavyup ||
				(syst_ == heavyuppt && pt <= getMedian(med_bpt)) ||
				(syst_ == heavydownpt && pt >= getMedian(med_bpt)) ||
				(syst_ == heavyupeta &&  abseta <= getMedian(med_beta)) ||
				(syst_ == heavydowneta && abseta >= getMedian(med_beta)) ){
			usecalib=flavcalib->up();
		}
		else if (syst_ == heavydown ||
				(syst_ == heavyuppt && pt > getMedian(med_bpt)) ||
				(syst_ == heavydownpt && pt < getMedian(med_bpt)) ||
				(syst_ == heavyupeta &&  abseta > getMedian(med_beta)) ||
				(syst_ == heavydowneta && abseta < getMedian(med_beta)) ){
			usecalib=flavcalib->down();
		}
		else if (syst_ == allup)
			usecalib=flavcalib->up();
		else if (syst_ == alldown)
			usecalib=flavcalib->down();
		else
			usecalib=flavcalib->central();
	}

	else if(jetflav == BTagEntry::FLAV_UDSG){

		if(udsgcalibration_)  flavcalib=udsgcalibration_;
		if(!flavcalib)
			throw std::logic_error("bTagSFBase::jetSF: calibration not loaded");

		// For the time being use the median points of b-falvour jets!!!
		if(syst_ == lightup ||
				(syst_ == lightuppt && pt >= getMedian(med_bpt)) ||
				(syst_ == lightdownpt && pt <= getMedian(med_bpt)) ||
				(syst_ == lightupeta && abseta >= getMedian(med_beta)) ||
				(syst_ == lightdowneta && abseta <= getMedian(med_beta))){
			usecalib=flavcalib->up();
		}
		else if(syst_ == lightdown ||
				(syst_ == lightuppt && pt < getMedian(med_bpt)) ||
				(syst_ == lightdownpt && pt > getMedian(med_bpt)) ||
				(syst_ == lightupeta && abseta < getMedian(med_beta)) ||
				(syst_ == lightdowneta && abseta > getMedian(med_beta))) {
			usecalib=flavcalib->down();
		}
		else if (syst_ == allup)
			usecalib=flavcalib->up();
		else if (syst_ == alldown)
			usecalib=flavcalib->down();
		else
			usecalib=flavcalib->central();
	}

	//everything went fine
	if(usecalib){
		return usecalib->eval(jetflav,eta,pt,discr);
	}
	throw std::logic_error("bTagSFBase::jetSF: systematic variation that was requested was never loaded");

}


} //namespace
