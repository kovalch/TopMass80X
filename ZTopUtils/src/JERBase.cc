#include "../interface/JERBase.h"
#include <stdexcept>
namespace ztop {

void JERBase::setSystematics(std::string type) {
	resranges_.clear();

	// these are the ones from 2011 onwards
	/*
    resranges_ << 0 << 0.5 << 1.1 << 1.7 << 2.3 << 5;

    if (type == "def") { //standard
        resfactors_.clear();
        resfactors_ << 1.052 << 1.057 << 1.096 << 1.134 << 1.288;
        std::cout << "JER set to default" << std::endl;
    } else if (type == "down") {
        resfactors_.clear();
        resfactors_ << 1.115 << 1.114 << 1.162 << 1.228 << 1.488;
        std::cout << "JER set to syst down" << std::endl;
    } else if (type == "up") {
        resfactors_.clear();
        resfactors_ << 0.990 << 1.001 << 1.032 << 1.042 << 1.089;
        std::cout << "JER set to syst up" << std::endl;
    }
	 */

	//new ones from August 2014 (for 2012 data)
	resranges_     << 0 << 0.5    << 1.1   << 1.7    << 2.3   << 2.8    << 3.2  << 5.0;

	if (type == "def") { //standard
		resfactors_.clear();
		resfactors_ << 1.079 << 1.099 << 1.121 << 1.208 << 1.254 << 1.395 << 1.056;
		std::cout << "JER set to default" << std::endl;
	} else if (type == "down") {
		resfactors_.clear();
		resfactors_ << 1.053 << 1.071 << 1.092 << 1.162 << 1.192 << 1.332 << 0.865;
		std::cout << "JER set to syst down" << std::endl;
	} else if (type == "up") {
		resfactors_.clear();
		resfactors_ << 1.105 << 1.127 << 1.150 << 1.254 << 1.316 << 1.458 << 1.247 ;
		std::cout << "JER set to syst up" << std::endl;
	}

	// when setting ranges and factors here, we assume that the sizes
	// are coded in the right format, check is only done if sizes have set up
	// externally
	checksizesonfirstcall_=false;

}

void JERBase::correctP4(float & recopt, float& recoabs_eta, float & recophi, float & recom, //full lorentzvector
		const float & genpt) {
	if (genpt < 1)
		return;

	if(checksizesonfirstcall_){
		if(resranges_.size() != resfactors_.size()+1){
			throw std::out_of_range("JERBase: Resolution ranges and factors don't match in size. The range vector defines\
					 bin boundaries, hence must have one more entry than the factors vector!");
		}

		checksizesonfirstcall_=false;
	}

	std::vector<float>::const_iterator it=std::lower_bound(resranges_.begin(),
			resranges_.end(), recoabs_eta);
	size_t etabin=0;
	if(recoabs_eta==*it)
		etabin= it-resranges_.begin();
	else
		etabin= it-resranges_.begin()-1;
	if(etabin==resranges_.size()-1) //last entry
		etabin--;

	double deltapt = (1. - resfactors_.at(etabin))* (recopt - genpt);
	double scale = std::max(0., recopt + deltapt) / recopt;
	recopt*=scale;
	recom*=scale;
}
} //namespace
