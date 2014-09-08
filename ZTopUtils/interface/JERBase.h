#ifndef JERBASE_H
#define JERBASE_H


#include <vector>
#include "../interface/miscUtils.h"
#include <algorithm>
#include <iostream>

// following recommendation https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution;   https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefSyst

namespace ztop {

// typedef ROOT::Math::PolarLorentzVector<ROOT::Math::PxPyPzE4D<double> > PolarLorentzVector;
/**
 *
 WHATEVER you add as functions, please don't use exit() in case an error occurs.
 replace it with either:
 - throw an exception (throw std::logic_error("sometext") or std::runtime_error("");)
 - return something (-1 or another int for dubugging)

 */

class JERBase {
public:
	JERBase():checksizesonfirstcall_(true) {
		setSystematics("def");
	}
	~JERBase() {
	}
	/**
	 * Provide bin boundaries. That means, this vector's size should
	 * exceed the factors vector size by exactly one entry
	 * These will overwrite anything
	 * automatically set by setSystematics()
	 */
	void setResolutionEtaRanges(std::vector<float> ranges) {
		resranges_ = ranges;
		checksizesonfirstcall_=true;
	}
	/**
	 * Provide resolution factors. These will overwrite anything
	 * automatically set by setSystematics()
	 */
	void setResolutionFactors(std::vector<float> factors) {
		resfactors_ = factors;
		checksizesonfirstcall_=true;
	}

	void setSystematics(std::string type);

	void correctP4(float & recopt, float& recoabs_eta, float & recophi, float & recom, //full lorentzvector
			const float & genpt) ;

protected:

	std::vector<float> resranges_;
	std::vector<float> resfactors_;

	bool checksizesonfirstcall_;

};
}
#endif
