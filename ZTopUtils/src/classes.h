#include <vector>
#include "../interface/JERBase.h"
#include "../interface/JECBase.h"
#include "../interface/bTagBase.h"

#include "../ext/interface/JetCorrectorParameters.h"
#include "../ext/interface/JetCorrectionUncertainty.h"
#include "../ext/interface/FactorizedJetCorrector.h"
#include "../ext/interface/SimpleJetCorrectionUncertainty.h"
#include "../ext/interface/SimpleJetCorrector.h"
#include "../ext/interface/JetCorrectorParameters.h"


namespace
{
  struct dict {

    std::vector<double> klklkl;
    std::vector<TString> kkj;
    std::pair<TString,std::vector<double> > jj;
    std::vector<std::pair<TString,std::vector<double> > > j;
    std::vector<long int> kl;

    std::vector<TH2D> ijij;
    std::map<TString,std::vector<TH2D> > kokd;


    

    ztop::bTagBase ked;
    ztop::JECBase kkss;
    ztop::JERBase ksdsk;
    std::vector<ztop::bTagBase> kedV;
    std::vector<ztop::JECBase> kkssV;
    std::vector<ztop::JERBase> ksdskV;


///external


JetCorrectorParameters corr;
JetCorrectorParameters::Definitions def;
JetCorrectorParameters::Record record;
std::vector<JetCorrectorParameters> corrv;
std::vector<JetCorrectorParameters::Record> recordv;
JetCorrectorParametersCollection coll;
JetCorrectorParametersCollection::pair_type pair_type;
JetCorrectorParametersCollection::collection_type colltype;
std::vector<JetCorrectorParametersCollection> collv;

  };
}

