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




bTagSFBase::bTagSFBase():bTagEfficiency(), calib_(0),currentcalibration_(0),syst_(nominal),wpval_(0){

}
bTagSFBase::~bTagSFBase(){
	if(calib_) delete calib_;
	if(currentcalibration_) delete currentcalibration_;
}

void bTagSFBase::loadSF(const std::string& filename, BTagEntry::OperatingPoint op,
		const std::string& tagger, std::string measurementType,std::string systsourceup, std::string  systsourcedown){

	if(syst_!=nominal && (systsourceup.length()<1 || systsourcedown.length()< 1)){
		throw std::logic_error("bTagSFBase::loadSF: systematic variation required without loading systematic variations for scale factors");
	}

	if(calib_) delete calib_;
	if(currentcalibration_) delete currentcalibration_;
	try{
		calib_ = new BTagCalibration(tagger,filename);
		currentcalibration_ = new calibrations(calib_,op,measurementType,systsourceup,systsourcedown);
	}
	catch(...){
		throw std::runtime_error("bTagSFBase::loadSF: problem loading payload");
	}


	if(op == BTagEntry::OP_LOOSE)
		wpstring_="loose";
	else if(op == BTagEntry::OP_MEDIUM)
		wpstring_="medium";
	else if(op == BTagEntry::OP_TIGHT)
		wpstring_="tight";
	else
		wpstring_="shape";

	wpval_=getWPDiscrValue(tagger,op);

	if(debug)
		std::cout << "bTagSFBase::loadSF: set up SF for "<< tagger <<", "<< measurementType << std::endl;


	//fast check
	try{
		jetSF(BTagEntry::FLAV_B,40, 0,0.5);
	}
	catch(...){
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





void bTagSFBase::setSystematics(systematics sys){
	if(sys!=nominal && (! currentcalibration_->up() || ! currentcalibration_->down())){
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
	else if(tagger == "csvv1"){
		if(op == BTagEntry::OP_LOOSE)
			return -1;
		else if(op == BTagEntry::OP_MEDIUM)
			return -1;
		else if(op == BTagEntry::OP_TIGHT)
			return -1;
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
	float abseta=std::abs(eta);
	if(jetflav == BTagEntry::FLAV_B || jetflav == BTagEntry::FLAV_C){
		if (syst_ == heavyup ||
				(syst_ == heavyuppt && pt < getMedian(med_bpt)) || (syst_ == heavydownpt && pt > getMedian(med_bpt)) ||
				(syst_ == heavyupeta &&  abseta < getMedian(med_beta)) || (syst_ == heavydowneta && abseta > getMedian(med_beta)) ){
			usecalib=currentcalibration_->up();
		}
		else if (syst_ == heavydown ||
				(syst_ == heavyuppt && pt > getMedian(med_bpt)) || (syst_ == heavydownpt && pt < getMedian(med_bpt)) ||
				(syst_ == heavyupeta &&  abseta > getMedian(med_beta)) || (syst_ == heavydowneta && abseta < getMedian(med_beta)) ){
			usecalib=currentcalibration_->down();
		}
		else
			usecalib=currentcalibration_->central();
	}

	else if(jetflav == BTagEntry::FLAV_UDSG){
		// For the time being use the median points of b-falvour jets!!!
		if((syst_ == lightuppt && pt > getMedian(med_bpt)) ||
				(syst_ == lightdownpt && pt < getMedian(med_bpt)) ||
				(syst_ == lightupeta && abseta > getMedian(med_beta)) ||
				(syst_ == lightdowneta && abseta < getMedian(med_beta))){
			usecalib=currentcalibration_->up();
		}
		else if((syst_ == lightuppt && pt < getMedian(med_bpt)) ||
				(syst_ == lightdownpt && pt > getMedian(med_bpt)) ||
				(syst_ == lightupeta && abseta < getMedian(med_beta)) ||
				(syst_ == lightdowneta && abseta > getMedian(med_beta))) {
			usecalib=currentcalibration_->down();
		}
		else
			usecalib=currentcalibration_->central();
	}
	if(syst_ == allup)
		usecalib=currentcalibration_->up();
	else if (syst_ == alldown)
		usecalib=currentcalibration_->down();

	//everything went fine
	if(usecalib){
		return usecalib->eval(jetflav,eta,pt,discr);
	}
	throw std::logic_error("bTagSFBase::jetSF: systematic variation that was requested was never loaded");

}


} //namespace
