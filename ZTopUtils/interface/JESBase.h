#ifndef JESBASE_H
#define JESBASE_H

#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <string>
#include "TString.h"

class JetCorrectionUncertainty;
class FactorizedJetCorrector;

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

class JESBase {
public:

    JESBase() {
        is2012_ = true;
        totalunc_ = 0;
        noupdown_=0;
    }
    JESBase(const ztop::JESBase &);
    JESBase & operator =(const ztop::JESBase &);
    ~JESBase();

    void setFile(std::string pathToFile, bool quiet = false);
    void setSystematics(std::string); //! up, down, no
    void setIs2012(bool is) {
        is2012_ = is;
        std::cout << "JES mode changed; set File again!" << std::endl;
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
    void applyUncertainties(float & pt, float& eta, float & phi, float& m, int jetFlavour=-9999);
    
    /// Get systematic uncertainty factor to be directly multiplied with jet P4
    double uncertaintyFactor(const double& pt, const double& eta, const int jetFlavour =-9999)const;
    
    /**
     * provides correction value due to JES, which should be used for pt, et correction of the particular uncorrected jet used
     * as the input for this function
     */
    double correctionForUncorrectedJet(const double& jetInitialArea, const double& jetInitialEta, const double& jetInitialPt, const double& rho);
    
    /**
     * sets files which will be used in getCorrectionValueForUncorrectedJet() for getting correction values due to JES,
     * boolean "isMC" triggers usage of "L2L3Residual" JES corrections for Data
     */
    void configureFactorizedJetCorrector(const TString* filePathL1, const TString* filePathL2, const TString* filePathL3, const TString* filePathL2L3, const bool& isMC);
    
    /// Configure factorised jet corrector from a vector of filenames, each representing one correction in the order of the vector
    void configureFactorizedJetCorrector(const std::vector<std::string>& v_filename);
    
protected:
    std::vector<unsigned int> & sources() {
        return sources_;
    }
    std::string pathToFile_;
    std::vector<JetCorrectionUncertainty*> vsrc_;
    JetCorrectionUncertainty* totalunc_;
    int noupdown_;
    std::vector<unsigned int> sources_;
    std::map<std::string,unsigned int> sourcenames_;
    bool is2012_;

    std::vector<int> restricttoflav_;

    void copyFrom(const ztop::JESBase &);
        
    FactorizedJetCorrector *jetCorrector_;

};

}
#endif
