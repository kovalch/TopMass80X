#include <vector>
#include "../interface/JERBase.h"
#include "../interface/JECBase.h"
#include "../interface/bTagBase.h"
#include "../interface/RecoilCorrector.h"

#include "../ext/interface/JetCorrectorParameters.h"
#include "../ext/interface/JetCorrectionUncertainty.h"
//#include "../ext/interface/FactorizedJetCorrector.h"
#include "../ext/interface/SimpleJetCorrectionUncertainty.h"
//#include "../ext/interface/SimpleJetCorrector.h"
#include "../ext/interface/JetCorrectorParameters.h"

namespace {
struct dict {

    std::vector<double> klklkl;
    std::vector<TString> kkj;
    std::pair<TString, std::vector<double> > jj;
    std::vector<std::pair<TString, std::vector<double> > > j;
    std::vector<long int> kl;

    std::vector<TH2D> ijij;
    std::map<TString, std::vector<TH2D> > kokd;
    std::vector<TH1D> ifffjij;

    std::map<std::string, std::vector<TH2D> > sdfdsf;
    std::map<std::string, std::vector<float> > formedianmap;
    std::vector<std::vector<std::vector<TH1D> > > forshapRWhistos;

    ztop::RecoilCorrector corrector;
    ztop::bTagBase ked;
    ztop::JECBase kkss;
    ztop::JERBase ksdsk;
    std::vector<ztop::RecoilCorrector> correctorV;
    std::vector<ztop::bTagBase> kedV;
    std::vector<ztop::JECBase> kkssV;
    std::vector<ztop::JERBase> ksdskV;

    ///external

    ztop::JetCorrectorParameters corr;
    ztop::JetCorrectorParameters::Definitions def;
    ztop::JetCorrectorParameters::Record record;
    std::vector<ztop::JetCorrectorParameters> corrv;
    std::vector<ztop::JetCorrectorParameters::Record> recordv;
    ztop::JetCorrectorParametersCollection coll;
    ztop::JetCorrectorParametersCollection::pair_type pair_type;
    ztop::JetCorrectorParametersCollection::collection_type colltype;
    std::vector<ztop::JetCorrectorParametersCollection> collv;

};
}

