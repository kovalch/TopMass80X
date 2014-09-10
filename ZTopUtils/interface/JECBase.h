#ifndef JECBASE_H
#define JECBASE_H

#include "../ext/interface/JetCorrectorParameters.h"
#include "../ext/interface/JetCorrectionUncertainty.h"
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <string>

namespace ztop {

/**
 *
 *
 *
 *
 *
 WHATEVER you add as functions, please don't use exit() in case an error occurs.
 replace it with either:
 - throw an exception (throw std::logic_error("sometext") or std::runtime_error("");)
 - return something (-1 or another int for dubugging)

 */

class JECBase {
public:

	JECBase() {
		is2012_ = true;
		totalunc_ = 0;
		noupdown_=0;
	}
	JECBase(const ztop::JECBase &);
	JECBase & operator =(const ztop::JECBase &);
	~JECBase();

	void setFile(std::string pathToFile, bool quiet = false);
	void setSystematics(std::string); //! up, down, no
	void setIs2012(bool is) {
		is2012_ = is;
		std::cout << "JEC mode changed; set File again!" << std::endl;
	}


	/**
	 * ADDS! a source with name to the sources to be varied
	 */
	void setSource(const std::string&);
	/**
	 * resets configuration as far as sources to be varied are concerned
	 */
	void clearSources(){sources_.clear();}


	/**
	 * returns a vector of available source names
	 */
	std::vector<std::string> getSourceNames()const;

	/**
	 * Applies uncertainties only if they match the flavour specified here
	 * This function ADDS and entry! to clear, use clearRestrictToFlavour()
	 */
	void restrictToFlavour(int flav);

	/**
	 * Clears all flavour restrictions added before
	 */
	void clearRestrictToFlavour();

	/**
	 * applies the uncertainties on jet quantities.
	 * jetFlavour is the absolute of the one given by PAT! If none specified (<9998), uncertainties will be
	 * applied to all in the same manner. (default for backward compatibility)
	 * The same is true if no restriction is specified in restrictToAbsFlavour();
	 */
	void applyJECUncertainties(float & pt, float& eta, float & phi, float& m, int jetFlavour=-9999);

protected:
	std::vector<unsigned int> & sources() {
		return sources_;
	}
	std::string pathToFile_;
	std::vector<ztop::JetCorrectionUncertainty*> vsrc_;
	ztop::JetCorrectionUncertainty* totalunc_;
	int noupdown_;
	std::vector<unsigned int> sources_;
	std::map<std::string,unsigned int> sourcenames_;
	bool is2012_;

	std::vector<int> restricttoflav_;

	void copyFrom(const ztop::JECBase &);

};

}
#endif
