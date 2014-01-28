/*
 * muonRochCorr.cc
 *
 *  Created on: Jan 24, 2014
 *      Author: kiesej
 */




#include "../interface/muonRochCorr.h"

namespace ztop{

muonRochCorr::muonRochCorr():isMC_(true),syst_(nominal),runopt_(0){
    /* nothing */
}

muonRochCorr::muonRochCorr(int seed ):isMC_(true),syst_(nominal),corrector_(seed),runopt_(0){
    /* nothing */
}

void muonRochCorr::correctP4(TLorentzVector & vec, const int & charge, float& qterr){
    if(isMC_){
        corrector_.momcor_mc(vec,(float)charge,runopt_,qterr);
    }
    else{
        corrector_.momcor_data(vec,(float)charge,runopt_,qterr);
    }
}

void muonRochCorr::correctP4(TLorentzVector & vec, const int & charge){
    float tmp=1;
    correctP4(vec,charge, tmp);
}


}//namespace
