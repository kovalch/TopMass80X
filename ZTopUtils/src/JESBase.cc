#include "../interface/JESBase.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include <stdexcept>
#include "unistd.h"

namespace ztop {

void JESBase::copyFrom(const ztop::JESBase & old) {
    totalunc_ = 0;
    is2012_ = old.is2012_;
    if (!old.pathToFile_.empty())
        setFile(old.pathToFile_, true); // JESUncertainties... don't provide a copy contructor..
    noupdown_ = old.noupdown_;
    sources_ = old.sources_;
    sourcenames_=old.sourcenames_;
    restricttoflav_=old.restricttoflav_;
}

JESBase::JESBase(const ztop::JESBase & old) {
    copyFrom(old);
}

JESBase & JESBase::operator =(const ztop::JESBase & old) {
    copyFrom(old);
    return *this;
}

JESBase::~JESBase() {

    if (totalunc_)
        delete totalunc_;
    for (unsigned int i = 0; i < vsrc_.size(); i++) {
        if (vsrc_.at(i))
            delete vsrc_.at(i);
    }

}

void JESBase::setSource(const std::string& str){
    std::map< std::string,unsigned int>::const_iterator mit=sourcenames_.find(str);
    if(mit!=sourcenames_.end()){
        sources_.push_back(mit->second);
        return;
    }
    std::cout << "JESBase::setSource: available JES sources: " << std::endl;
    for(mit=sourcenames_.begin();mit!=sourcenames_.end();++mit)
        std::cout << mit->first << std::endl;
    throw std::runtime_error("JESBase::setSource: Source unknown, chose one of the above");
}

void JESBase::setFile(std::string pathToFile, bool quiet) {
    if (pathToFile.empty()) {
        std::cout << "JESBase::setFile: path empty!" << std::endl;
        throw std::runtime_error("JESBase::setFile: path empty!");
    }
    std::ifstream check_file(pathToFile.data());
    if (!check_file.good()) {
        std::cout << "JESBase::setFile: cannot open file!" << pathToFile
                << std::endl;
        throw std::runtime_error("JESBase::setFile: cannot open file!");
    }
    if (!quiet)
        std::cout << "setting JES uncertainties file to: " << pathToFile
        << std::endl;
    pathToFile_ = pathToFile;

    //delete old JES

    if (totalunc_)
        delete totalunc_;
    for (unsigned int i = 0; i < vsrc_.size(); i++) {
        if (vsrc_[i])
            delete vsrc_[i];
    }
    vsrc_.clear();

    const int nsrc = 16;
    const char* srcnames[nsrc] = { "Absolute", //0                         //uncorr
            "HighPtExtra",    //1
            "SinglePion",     //2
            "Flavor",         //3
            "Time",           //4
            "RelativeJEREC1", //5
            "RelativeJEREC2", //6
            "RelativeJERHF",  //7
            "RelativeStatEC2",  //8
            "RelativeStatHF", //9
            "RelativeFSR",    //10
            "PileUpDataMC",   //11                         //uncorr
            "PileUpOOT",      //12
            "PileUpPt",       //13
            "PileUpBias",     //14                         //uncorr
            "PileUpJetRate" }; //15

    const int nsrc12 = 41;
    const char* srcnames12[nsrc12] = {
            "AbsoluteStat", //0
            "AbsoluteScale", //1
            "AbsoluteFlavMap", //2
            "AbsoluteMPFBias", //3
            "HighPtExtra",     //4
            "SinglePionECAL",  //5
            "SinglePionHCAL",  //6
            "FlavorQCD",         //7
            "Time",                //8
            "RelativeJEREC1",     //9
            "RelativeJEREC2",     //10
            "RelativeJERHF",    //11
            "RelativePtBB",        //12
            "RelativePtEC1",     //13
            "RelativePtEC2",     //14
            "RelativePtHF",     //15
            "RelativeFSR",        //16
            "RelativeStatEC2",     //17
            "RelativeStatHF",    //18
            "PileUpDataMC",        //19
            "PileUpPtBB",         //20
            "PileUpPtEC",         //21
            "PileUpPtHF",        //22
            "PileUpBias",        //23
            "SubTotalPileUp",    //24
            "SubTotalRelative",    //25
            "SubTotalPt",        //26
            "SubTotalMC",        //27
            "Total",            //28
            "TotalNoFlavor",    //29
            "FlavorZJet",        //30
            "FlavorPhotonJet",    //31
            "FlavorPureGluon",    //32
            "FlavorPureQuark",    //33
            "FlavorPureCharm",    //34
            "FlavorPureBottom", //35
            "CorrelationGroupMPFInSitu", //36
            "CorrelationGroupIntercalibration", //37
            "CorrelationGroupbJES", //38
            "CorrelationGroupFlavor", //39
            "CorrelationGroupUncorrelated",//40
    };

    if (is2012_) {
        for (int isrc = 0; isrc < nsrc12; isrc++) {
            const char *name = srcnames12[isrc];
            sourcenames_[name] = isrc;
            bool got = true;
            JetCorrectionUncertainty *unc = 0;
            try {
                unc = new JetCorrectionUncertainty(
                        JetCorrectorParameters(pathToFile.data(), name));
            }
            catch(std::runtime_error &rte) {
                std::cout << "JESBase::setFile: Uncertainty for source " << name
                        << " not found! Skipping" << std::endl;
                got = false;
                sleep(2);
            }
            if (got)
                vsrc_.push_back(unc);
        }
    } else {
        for (int isrc = 0; isrc < nsrc; isrc++) {
            const char *name = srcnames[isrc];
            sourcenames_[name] = isrc;
            bool got = true;
            JetCorrectionUncertainty *unc = 0;
            try {
                unc = new JetCorrectionUncertainty(
                        JetCorrectorParameters(pathToFile.data(), name));
            }
            catch(std::runtime_error &rte) {
                std::cout << "JESBase::setFile: Uncertainty for source " << name
                        << " not found! Skipping" << std::endl;
                got = false;
                sleep(2);
            }
            if (got)
                vsrc_.push_back(unc);
        }
    }
    totalunc_ = new JetCorrectionUncertainty(
            JetCorrectorParameters(pathToFile.data(), "Total"));
}

void JESBase::setSystematics(std::string set) {
    if (set == "up") {
        noupdown_ = 1;
        std::cout << "JESBase::setSystematics: Systematics set to: " << set
                << std::endl;
    } else if (set == "down") {
        noupdown_ = -1;
        std::cout << "JESBase::setSystematics: Systematics set to: " << set
                << std::endl;
    } else if (set == "no") {
        noupdown_ = 0;
        std::cout << "JESBase::setSystematics: Systematics set to: " << set
                << std::endl;
    } else {
        std::cout << "JESBase::setSystematics: String " << set
                << " not allowed. available options: up, down, no" << std::endl;
    }
}


std::vector<std::string> JESBase::getSourceNames()const{
    std::vector<std::string>  out;
    for(std::map<std::string,unsigned int > :: const_iterator mit=sourcenames_.begin();mit!=sourcenames_.end();++mit){
        out.push_back(mit->first);
    }
    return out;
}



void JESBase::restrictToFlavour(int flav){
    restricttoflav_.push_back(flav);
}

void JESBase::clearRestrictToFlavour(){
    restricttoflav_.clear();
}


void JESBase::applyUncertainties(float & pt, float& eta, float & phi, float& m, int jetFlavour) {

    if (noupdown_ == 0) // no variation
        return;

    if (!(totalunc_)) { // nothing set. exit?!?
        std::cout << "JESBase::applyUncertainties: no inputfile set, exit"
                << std::endl;
        throw std::logic_error(
                "JESBase::applyUncertainties: no inputfile set, exit");
    }
    if(restricttoflav_.size()>0){
        if(std::find(restricttoflav_.begin(),restricttoflav_.end(),jetFlavour) == restricttoflav_.end())
            return;
    }

    bool up = false;
    if (noupdown_ > 0)
        up = true;

    double dunc = 0;

    if (sources_.size() < 1) { //total
        totalunc_->setJetPt(pt);
        totalunc_->setJetEta(eta);
        dunc = totalunc_->getUncertainty(up);
    } else if (sources_.size() <= vsrc_.size()) {
        for (unsigned int i = 0; i < sources_.size(); i++) {
            if (sources_[i] < vsrc_.size()) { //for a spec source
                JetCorrectionUncertainty *unc = vsrc_[sources_[i]];
                unc->setJetPt(pt);
                unc->setJetEta(eta);
                double uncert = unc->getUncertainty(up);
                dunc = sqrt(dunc * dunc + uncert * uncert);
            } else {
                std::cout << "JESBase::applyUncertainties: source "
                        << sources_[i] << " doesn't exist." << std::endl;
            }
        }
    } else {
        std::cout
        << "JESBase::applyUncertainties: too many sources; must be below "
        << vsrc_.size() - 1 << "." << std::endl;
    }
    if (up){
        pt*=  (1 + dunc);
        m*=  (1 + dunc);
    }
    else{
        pt*=  (1 - dunc);
        m*=  (1 - dunc);
    }
}



double JESBase::uncertaintyFactor(const double& pt, const double& eta, const int jetFlavour)const
{
    if (noupdown_ == 0) // no variation
        return 1.;

    if (!(totalunc_)) { // nothing set. exit?!?
        std::cout << "JESBase::uncertaintyFactor: no inputfile set, exit"
                << std::endl;
        throw std::logic_error(
                "JESBase::uncertaintyFactor: no inputfile set, exit");
    }
    
    if(restricttoflav_.size() > 0){
        if(std::find(restricttoflav_.begin(),restricttoflav_.end(),jetFlavour) == restricttoflav_.end())
            return 1.;
    }

    bool up = false;
    if (noupdown_ > 0)
        up = true;

    double dunc = 0.;

    if (sources_.size() < 1) { //total
        totalunc_->setJetPt(pt);
        totalunc_->setJetEta(eta);
        dunc = totalunc_->getUncertainty(up);
    } else if (sources_.size() <= vsrc_.size()) {
        for (unsigned int i = 0; i < sources_.size(); i++) {
            if (sources_[i] < vsrc_.size()) { //for a spec source
                JetCorrectionUncertainty *unc = vsrc_[sources_[i]];
                unc->setJetPt(pt);
                unc->setJetEta(eta);
                double uncert = unc->getUncertainty(up);
                dunc = sqrt(dunc * dunc + uncert * uncert);
            } else {
                std::cout << "JESBase::uncertaintyFactor: source "
                        << sources_[i] << " doesn't exist." << std::endl;
            }
        }
    } else {
        std::cout
        << "JESBase::uncertaintyFactor: too many sources; must be below "
        << vsrc_.size() - 1 << "." << std::endl;
    }
    
    if(up) return 1. + dunc;
    else return 1. - dunc;
}



double JESBase::correctionForUncorrectedJet(const double& jetInitialArea, const double& jetInitialEta, const double& jetInitialPt, const double& rho) {

        jetCorrector_->setJetA(jetInitialArea);
        jetCorrector_->setJetEta(jetInitialEta);
        jetCorrector_->setJetPt(jetInitialPt);
        jetCorrector_->setRho(rho);
        
        return jetCorrector_->getCorrection();
}

void JESBase::configureFactorizedJetCorrector(const TString* filePathL1, const TString* filePathL2, const TString* filePathL3, const TString* filePathL2L3, const bool& isMC) {
            
        std::vector<JetCorrectorParameters> vJetPar;
        JetCorrectorParameters *jetParL1, *jetParL2, *jetParL3, *jetParL2L3;
        jetParL1 = new JetCorrectorParameters(filePathL1->Data());
        jetParL2 = new JetCorrectorParameters(filePathL2->Data());
        jetParL3 = new JetCorrectorParameters(filePathL3->Data());
        if(!isMC) jetParL2L3 = new JetCorrectorParameters(filePathL2L3->Data());
        else jetParL2L3 = 0;
        
        //Order must correspond to one provided by JEC  group: "https://twiki.cern.ch/twiki/bin/view/CMS/IntroToJEC"
        vJetPar.push_back(*jetParL1);
        vJetPar.push_back(*jetParL2);
        vJetPar.push_back(*jetParL3);
        if(!isMC) vJetPar.push_back(*jetParL2L3);
        
        jetCorrector_ = new FactorizedJetCorrector(vJetPar);
}



void JESBase::configureFactorizedJetCorrector(const std::vector<std::string>& v_filename)
{
    std::vector<JetCorrectorParameters> v_jetPar;
    for(const std::string& filename : v_filename){
        v_jetPar.push_back(JetCorrectorParameters(filename));
    }
    jetCorrector_ = new FactorizedJetCorrector(v_jetPar);
}



}

