/*
 * bTagSFBase.h
 *
 *  Created on: Apr 7, 2015
 *      Author: kiesej
 */

#ifndef BTAGSFBASE_H_
#define BTAGSFBASE_H_

#include "../ext/BTagCalibrationStandalone.h"
#include <fstream>
#include "bTagEfficiency.h"

namespace ztop{


/**
 * Warning: if jets are out of range (20 < pt < ~500) (|eta|< 2.4), the SF will be 0!
 *
 * Class to apply b-tag SF
 * The most important function to set up the scale factors is
 * loadSF(..)
 * It loads and configures the official tools BTagCalibration and BTagCalibrationReader
 * For information on the options, please refer to:
 *  https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagCalibration
 *
 * Also shape reweighting is supported, as well as random-tagging. Event-wise
 * scale factors are not supported anymore, but could be implemented on demand
 *
 * A small test program can be found in ../test/testBTag.cc
 * Whenever you test against the old implementation, please keep in mind
 * that the old implementation used b-tag SF without ttbar (now called mujets),
 * but used ttbar data for the light-jet SF. With this new tool this mixing is
 * not possible (right now).
 *
 * The IO for efficiency histograms is inherited from the bTagEfficiency class.
 * This has also changed wrt to the old bTagBase where this was up to the user.
 *
 * The input files can be found at ../data (right now only 8 TeV)
 *
 * Instances of this class can not be copied or overwritten.
 *
 * Basic workflow for random tagging (different from bTagBase)
 *
 *   bTagSFBase bsf;
 *   bsf.loadSF(...,...)
 *   bsf.setMakeEff(true);
 *
 *   eventloop for efficiencies{
 *      for each jet{
 *        bsf.fillEff(...)
 *      }
 *    }
 *
 *    bsf.writeToTFile(...)
 *
 *  ////////
 *    bsf.readFromTFile(...)
 *    eventloop for analysis{
 *      for each jet{
 *        bool istagged=bsf.jetIsTagged(...)
 *      }
 *    }
 *
 * Basic workflow for disc shape reweighting (similar to before)
 *
 *    bTagSFBase bsf;
 *    bsf.loadSF(...,...)
 *    eventloop {
 *      for each b-jet{
 *        float jetSF=getJetDiscrShapeWeight(....)
 *      }
 *    }
 *
 * The following functionality has been successfully tested:
 *  -> Random-tagging consistency with old bTagBase (test program)
 *
 * The following functionality needs to be tested:
 *  -> Full implementation in the TopHiggs framework
 *     -> test random tagging
 *     -> test shape reweighting (not supported by file, yet!)
 *  -> Full integration in the TtZAnalysis framework
 *
 * Once tests are done, please report the status here (or just delete the lines)
 *
 */
class bTagSFBase : public bTagEfficiency  {

public:

	static bool debug;

	enum systematics {
		nominal,

		// shape uncertainties for norm. diff measurements
		heavyup, heavydown,
		lightup, lightdown,
		heavyuppt, heavydownpt,
		heavyupeta, heavydowneta,
		lightuppt, lightdownpt,
		lightupeta, lightdowneta,

		allup,alldown,

		// uncertainties on the discriminator shape reweighting still partially missing


		// leave this as last entry
		length_syst
	};

	bTagSFBase();
	~bTagSFBase();

private: //dont allow constructors
	bTagSFBase(const bTagSFBase&rhs):bTagEfficiency(),calib_(0),fullcalibration_(0),bccalibration_(0), udsgcalibration_(0),syst_(nominal),wpval_(0){}
	bTagSFBase& operator = (const bTagSFBase&rhs){return *this;}

public:

	void fillEff(const float& pt, const float&eta,
			const int& genpartflav, const float &disc,
			const float &puweight){
		bTagEfficiency::fillEff(pt,eta,genpartflav,isAboveWorkingPoint(disc),puweight);
	}

	/**
	 * loads full calibration with all SF for all flavours
	 * default choice!
	 */
	void loadSF(const std::string& filename, BTagEntry::OperatingPoint op,
			const std::string& tagger, const std::string measurementType="",
			const std::string systsourceup="", const std::string  systsourcedown="");

	/**
	 * loads flavour specific calibration. Overwrites default calibration only for specified flavour
	 */
	void loadBCSF(const std::string& filename, BTagEntry::OperatingPoint op,
			const std::string& tagger, const std::string measurementType="",
			const std::string systsourceup="", const std::string  systsourcedown="");


	void loadUDSGSF(const std::string& filename, BTagEntry::OperatingPoint op,
			const std::string& tagger, const std::string measurementType="",
			const std::string systsourceup="", const std::string  systsourcedown="");



	void setSystematics(systematics sys);

	bool isAboveWorkingPoint(const float & val)const {return val>wpval_;}

	//don't use this in heavy-load loops, just to get a value ONCE
	float getWPDiscrValue(const std::string& tagger, BTagEntry::OperatingPoint op)const;

	//for usage in loops
	const float& getWPDiscrValue()const{return wpval_;}

	//only for couts or something like that
	std::string getWorkingPointString()const;



	bool jetIsTagged(const float& pt, const float& eta, const int & genPartonFlavor,
			const float & tagValue, const unsigned int & seed) const;


	float getJetDiscrShapeWeight(const float & pt, const float& eta,
			const int & genPartonFlavor, const float& jetdiscr)const;
	/**
	 * resets everything
	 */
	void clear();

private:
	BTagCalibration * calib_;
	class calibrations{
	public:
		calibrations();
		calibrations(BTagCalibration *,BTagEntry::OperatingPoint,std::string,std::string , std::string );
		~calibrations();

		const BTagCalibrationReader* up()const{return up_;}
		const BTagCalibrationReader* down()const{return down_;}
		const BTagCalibrationReader* central()const{return central_;}
	private:
		BTagCalibrationReader *up_,*down_,*central_;
		calibrations(const calibrations&):up_(0),down_(0),central_(0){}
		calibrations& operator = (const calibrations&){return *this;}

	};
	calibrations * fullcalibration_, *bccalibration_, *udsgcalibration_;
	systematics syst_;

	std::string wpstring_;
	float wpval_;
	std::string filename_;

	//public: // debug
	//returns scale factor/shape weight depending on the systematics set
	float jetSF(const BTagEntry::JetFlavor& jetflav, const float& pt, const float& eta, const float& discr)const;

	void loadSFToCalib(calibrations** cal, const std::string& filename, BTagEntry::OperatingPoint op,
			const std::string& tagger, const std::string measurementType,
			const std::string systsourceup, const std::string  systsourcedown);

	//per sample efficiencies are only loaded


};

}

#endif /* BTAGSFBASE_H_ */
