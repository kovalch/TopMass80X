/*
 * muonRochCorr.h
 *
 *  Created on: Jan 24, 2014
 *      Author: kiesej
 */

#ifndef MUONROCHCORR_H_
#define MUONROCHCORR_H_

#include "../ext/src/rochcor2012jan22.h"

namespace ztop{
/**
 * this is only a wrapper class that does not depend on CMSSW
 * does not use z pt correction (for now)
 */
class muonRochCorr{
public:
    muonRochCorr();
    muonRochCorr(int seed);
    enum systematics{nominal,sys_up,sys_down};

    void setIsMC(bool ismc){isMC_=ismc;}

    void setSystematics(systematics sys){syst_=sys;}

    /**
     * this does not do anything right now.. maybe this changes
     */
    void setRunOpt(int opt){runopt_=opt;}

    void correctP4(TLorentzVector & , const int & charge, float& qterr);
    void correctP4(TLorentzVector & , const int & charge);

private:
    bool isMC_;
    systematics syst_;
    rochcor2012 corrector_;
    int runopt_;

};

}



#endif /* MUONROCHCORR_H_ */
