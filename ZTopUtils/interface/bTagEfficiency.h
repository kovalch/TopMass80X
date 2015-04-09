/*
 * bTagEfficiency.h
 *
 *  Created on: Apr 7, 2015
 *      Author: kiesej
 */

#ifndef BTAGEFFICIENCY_H_
#define BTAGEFFICIENCY_H_
#include "TH2D.h"
#include <string>
#include <vector>
#include "../ext/BTagCalibrationStandalone.h"
#include "TString.h"
#include "TFile.h"
#include <stdexcept>

namespace ztop{

/**
 * Class that handles the production of MC efficiencies for b-tagging.
 * It also incorporates input and output to TFiles.
 * The methods to apply and calculate scale factors are handled by a different class:
 * bTagSFBase
 *
 */
class bTagEfficiency {

public:
	bTagEfficiency();
	~bTagEfficiency();


    void setMakeEff(bool makee) {makeeffs_ = makee; }
    const bool& getMakeEff()const { return makeeffs_;}

	void fillEff(const float& pt, const float&abs_eta,
			const int& genpartflav, const bool &tagged,
			const float &puweight);




	//for IO: Fully root based!
	//produces efficiency histos and writes them to file
	void writeToTFile(TFile *f);
	void readFromTFile(TFile * f);

	static bool debug;

protected:

	enum btag_medians {
		med_bpt=0,
		med_beta=1,
		med_cpt=2,
		med_ceta=3,
		med_lpt=4,
		med_leta=5,
		med_length_medians=6};

	float getMedian(btag_medians med)const;
	float jetEff(const float &pt, const float& abs_eta,const BTagEntry::JetFlavor& jetflav) const;


private:
	static std::vector<bTagEfficiency*> all_;

	void initHistos();
	//enumeration scheme is crucial
	void makeEff();

	enum histoNames{
		hist_bjets_h=0,
		hist_btagged_h=1,
		hist_cjets_h=2,
		hist_ctagged_h=3,
		hist_ljets_h=4,
		hist_ltagged_h=5,
		hist_length_histoNames=6};

	enum effHistoNames{
		effhist_b=0,
		effhist_c=1,
		effhist_udsg=2,
		effhist_length_effHistoNames=3};


	static std::vector<std::string> getJetHistoOrderedNames();
	static std::vector<std::string> getEffHistoOrderedNames();


	histoNames assoFlavorToHist(const BTagEntry::JetFlavor& jetflav,const bool & tagged=false)const;
	effHistoNames assoFlavorToEff(const BTagEntry::JetFlavor& jetflav)const;

	float produceMedian(TH1 *)const ;
	template<class T>
	T  tryToGet(TFile * f, TString name)const;


	std::vector<TH2D> histos_;      //! bjets, btagged, cjets, ctagged, ljets, ltagged
	std::vector<TH2D> effhistos_;   //! beff, ceff, leff
	std::vector<float>  medians_;

	bool init_;
protected:
	bool makeeffs_;

	BTagEntry::JetFlavor jetFlavor(const int & genPartonFlavor)const;

};


inline
bTagEfficiency::histoNames bTagEfficiency::assoFlavorToHist(const BTagEntry::JetFlavor& jetflav,const bool & tagged)const{
	if(!tagged){
		if(jetflav ==  BTagEntry::FLAV_B) return hist_bjets_h;
		if(jetflav ==  BTagEntry::FLAV_C) return hist_cjets_h;
		if(jetflav ==  BTagEntry::FLAV_UDSG) return hist_ljets_h;
	}
	else{
		if(jetflav ==  BTagEntry::FLAV_B) return hist_btagged_h;
		if(jetflav ==  BTagEntry::FLAV_C) return hist_ctagged_h;
		if(jetflav ==  BTagEntry::FLAV_UDSG) return hist_ltagged_h;
	}
	return hist_length_histoNames;
}
inline
bTagEfficiency::effHistoNames bTagEfficiency::assoFlavorToEff(const BTagEntry::JetFlavor& jetflav)const{
	if(jetflav ==  BTagEntry::FLAV_B) return effhist_b;
	if(jetflav ==  BTagEntry::FLAV_C) return effhist_c;
	if(jetflav ==  BTagEntry::FLAV_UDSG) return effhist_udsg;
	return effhist_length_effHistoNames;
}

template<class T>
inline
T bTagEfficiency::tryToGet(TFile * f, TString name)const{
	T  * h= (T  *)f->Get(name);
	if(!h || h->IsZombie()){
		throw std::runtime_error("bTagEfficiency::tryToGet: Histogram not found " );
	}
	T out=*h;
	delete h;
	return out;
}

}


#endif /* BTAGEFFICIENCY_H_ */
