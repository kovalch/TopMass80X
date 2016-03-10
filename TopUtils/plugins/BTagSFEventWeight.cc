#include "FWCore/Utilities/interface/EDMException.h"
#include "TopAnalysis/TopUtils/plugins/BTagSFEventWeight.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


BTagSFEventWeight::BTagSFEventWeight(const edm::ParameterSet& cfg):
  jets_                   ( cfg.getParameter<edm::InputTag>    ( "jets"   ) ),
  bTagAlgo_               ( cfg.getParameter<std::string>      ("bTagAlgo") ),
  version_                ( cfg.getParameter<std::string>      ("version"  ) ),
  newRecipe_              ( cfg.getParameter<bool>             ("newRecipe") ),
  maxJets_                ( cfg.getParameter<int>              ("maxJets" ) ),
  sysVar_                 ( cfg.getParameter<std::string>      ("sysVar"  ) ),
  shapeVarPtThreshold_    ( cfg.getParameter<double>           ("shapeVarPtThreshold"  ) ),
  shapeVarEtaThreshold_   ( cfg.getParameter<double>           ("shapeVarEtaThreshold"  ) ),
  uncertaintySFb_         ( cfg.getParameter<double>           ("uncertaintySFb"  ) ),
  shapeDistortionFactor_  ( cfg.getParameter<double>           ("shapeDistortionFactor"  ) ),
  verbose_                ( cfg.getParameter<int>              ("verbose" ) ),
  filename_               ( cfg.getParameter<edm::FileInPath>  ("filename") ),
  noHistograms_           ( cfg.getParameter<bool>             ("noHistograms"))
{
  consumes<edm::View< pat::Jet > >(jets_);
  produces<double>();
  
  // set the edges of the last histo bin
  maxPtDB_     = 240.;
  maxPtMisTag_ = 520.;
  maxEta_      = 2.4;
  maxPt2012_   = 800.; // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript

  if(!noHistograms_){
    // load TFile Service
    edm::Service<TFileService> fs;
    if( !fs ){
      throw edm::Exception( edm::errors::Configuration,
			    "TFile Service is not registered in cfg file" );
    }
    /// booking of histogram for b tag eff SF
    hists_["effBTagEventSF"]     = fs->make<TH1F>( "effBTagEventSF", "effBTagEventSF", 100, 0.5, 1.5 );
    hists_["effBTagEventSFMean"] = fs->make<TH1F>( "effBTagEventSFMean", "effBTagEventSFMean", 1, 0, 1 );
  }
  
  /// getting efficiency histos from input files
  if(filename_.location()){
    file_ = new TFile((TString)filename_.fullPath());
    if(!(file_->IsZombie())){
      if(verbose_>=1) std::cout<<filename_.fullPath()<<" opened"<<std::endl;
      
	effHists_["NumBJetsPtEta"]       = (TH2F*) file_->Get("bTagEff/NumBJetsPtEta")->Clone();
	effHists_["NumBJetsTaggedPtEta"] = (TH2F*) file_->Get("bTagEff/NumBJetsTaggedPtEta")->Clone();
	effHists_["EffBJetsTaggedPtEta"] = (TH2F*) file_->Get("bTagEff/EffBJetsTaggedPtEta")->Clone();
	effHists_["NumCJetsPtEta"]       = (TH2F*) file_->Get("bTagEff/NumCJetsPtEta")->Clone();
	effHists_["NumCJetsTaggedPtEta"] = (TH2F*) file_->Get("bTagEff/NumCJetsTaggedPtEta")->Clone();
	effHists_["EffCJetsTaggedPtEta"] = (TH2F*) file_->Get("bTagEff/EffCJetsTaggedPtEta")->Clone();
	effHists_["NumLJetsPtEta"]       = (TH2F*) file_->Get("bTagEff/NumLJetsPtEta")->Clone();
	effHists_["NumLJetsTaggedPtEta"] = (TH2F*) file_->Get("bTagEff/NumLJetsTaggedPtEta")->Clone();
	effHists_["EffLJetsTaggedPtEta"] = (TH2F*) file_->Get("bTagEff/EffLJetsTaggedPtEta")->Clone();
	
	/// re-calculation of b tag efficiencies as input might be corrupted due to hadd
	if(effHists_.count("NumBJetsPtEta") && effHists_.count("NumBJetsTaggedPtEta") && effHists_.count("EffBJetsTaggedPtEta") &&
	  effHists_.count("NumCJetsPtEta") && effHists_.count("NumCJetsTaggedPtEta") && effHists_.count("EffCJetsTaggedPtEta") &&
	  effHists_.count("NumBJetsPtEta") && effHists_.count("NumBJetsTaggedPtEta") && effHists_.count("EffBJetsTaggedPtEta")) {
	  
	  effHists_.find("EffBJetsTaggedPtEta")->second->Reset();
	  effHists_.find("EffCJetsTaggedPtEta")->second->Reset();
	  effHists_.find("EffLJetsTaggedPtEta")->second->Reset();
	
	  effHists_.find("EffBJetsTaggedPtEta")->second->Divide(effHists_.find("NumBJetsTaggedPtEta")->second, 
	  effHists_.find("NumBJetsPtEta")->second,1,1,"B");
	  effHists_.find("EffCJetsTaggedPtEta")->second->Divide(effHists_.find("NumCJetsTaggedPtEta")->second, 
	  effHists_.find("NumCJetsPtEta")->second,1,1,"B");
	  effHists_.find("EffLJetsTaggedPtEta")->second->Divide(effHists_.find("NumLJetsTaggedPtEta")->second, 
	  effHists_.find("NumLJetsPtEta")->second,1,1,"B");
	 }
	 else{
	   std::cout<<"Eff.Histos not found!!!!! Efficiencies cannot be taken from this file!!! Default taken!"<<std::endl;
	   filename_ = edm::FileInPath();
	 }
    }
    else{
      std::cout<<filename_.fullPath()<<" not found!!!!! Efficiencies cannot be taken from this file!!! Default taken!"<<std::endl;
      filename_ = edm::FileInPath();
    }
  }
  
  /// load map from database
  measureMap_["BTAGBEFFCORR"]=PerformanceResult::BTAGBEFFCORR;
  measureMap_["BTAGBERRCORR"]=PerformanceResult::BTAGBERRCORR;
  measureMap_["BTAGLEFFCORR"]=PerformanceResult::BTAGLEFFCORR;
  measureMap_["BTAGLERRCORR"]=PerformanceResult::BTAGLERRCORR;
  
  mayConsume<edm::View< pat::Jet >>(jets_);
}

BTagSFEventWeight::~BTagSFEventWeight()
{
  if(filename_.location()) {if(!(file_->IsZombie())) file_->Close();}
}

void
BTagSFEventWeight::produce(edm::Event& evt, const edm::EventSetup& setup)
{
    //Setup measurement from database
//   setup.get<BTagPerformanceRecord>().get( "BTAG"+bTagAlgo_, perfHBTag);
//   setup.get<BTagPerformanceRecord>().get( "MISTAG"+bTagAlgo_, perfHMisTag);
  
  edm::Handle<edm::View< pat::Jet > > jets;
  evt.getByLabel(jets_, jets);
//   edm::Handle<std::vector<pat::Jet> > jets;
//   evt.getByToken(jets_, jets);
  


  double pt, eta;
  std::vector<double> oneMinusBEffies(0) , oneMinusBEffies_scaled(0);
  std::vector<double> oneMinusBMistags(0), oneMinusBMistags_scaled(0);
  
  double effBTagEvent_unscaled = 1.;
  double effBTagEvent_scaled   = 1.;

  
  if (!newRecipe_) {
    for(edm::View<pat::Jet>::const_iterator jet = jets->begin();jet != jets->end(); ++jet) {
      pt  = jet->pt();
      eta = std::abs(jet->eta());
      if(jet->partonFlavour() == 5 || jet->partonFlavour() == -5){
	      oneMinusBEffies               .push_back(1.- effBTag(pt, eta));
	      oneMinusBEffies_scaled        .push_back(1.- (effBTag(pt, eta) * effBTagSF(pt, eta, false)));
      }

      else if(jet->partonFlavour() == 4 || jet->partonFlavour() == -4){
	      oneMinusBMistags               .push_back(1.- effBTagCjet(pt, eta));
	      oneMinusBMistags_scaled        .push_back(1.-(effBTagCjet(pt, eta) * effBTagSF(pt, eta, true))); // ATTENTION: btag SF used for c-jets to with 2x the uncertainty!
      }
  
      else{
        oneMinusBMistags               .push_back(1.- effMisTag(pt, eta));
        oneMinusBMistags_scaled        .push_back(1.-(effMisTag(pt, eta) * effMisTagSF(pt, eta)));
      }
    }
     
  effBTagEvent_unscaled = effBTagEvent( oneMinusBEffies, oneMinusBMistags );
  effBTagEvent_scaled   = effBTagEvent( oneMinusBEffies_scaled, oneMinusBMistags_scaled );
  }
  
  // Use new recipe: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods#1a_Event_reweighting_using_scale
  else {
    if (bTagAlgo_ != "CSVM") std::cout<< "WARNING!!! New recipe only implemented for CSVM!!! CHECK!!!"<<std::endl;
    for(edm::View<pat::Jet>::const_iterator jet = jets->begin();jet != jets->end(); ++jet) {
      if (jet - jets->begin() == maxJets_) break;
      pt  = jet->pt();
      eta = std::abs(jet->eta());
      
      // tagged jets
      if (jet->bDiscriminator("combinedSecondaryVertexBJetTags")>0.679) {    //FIXME hardgecoded?? up-to-date? (dont think so...)
        if (abs(jet->partonFlavour()) == 5) {
          effBTagEvent_unscaled *= effBTag(pt, eta);
          effBTagEvent_scaled   *= effBTag(pt, eta) * effBTagSF(pt, eta, false);
        }
        else if (abs(jet->partonFlavour()) == 4) {
          effBTagEvent_unscaled *= effBTagCjet(pt, eta);
          effBTagEvent_scaled   *= effBTagCjet(pt, eta) * effBTagSF(pt, eta, true);
        }
        else {
          effBTagEvent_unscaled *= effMisTag(pt, eta);
          effBTagEvent_scaled   *= effMisTag(pt, eta) * effMisTagSF(pt, eta);
        }
      }
      
      // untagged jets
      else {
        if (abs(jet->partonFlavour()) == 5) {
          effBTagEvent_unscaled *= 1.- effBTag(pt, eta);
          effBTagEvent_scaled   *= 1.-(effBTag(pt, eta) * effBTagSF(pt, eta, false));
        }
        else if (abs(jet->partonFlavour()) == 4) {
          effBTagEvent_unscaled *= 1.- effBTagCjet(pt, eta);
          effBTagEvent_scaled   *= 1.-(effBTagCjet(pt, eta) * effBTagSF(pt, eta, true));
        }
        else {
          effBTagEvent_unscaled *= 1.- effMisTag(pt, eta);
          effBTagEvent_scaled   *= 1.-(effMisTag(pt, eta) * effMisTagSF(pt, eta));
        }
      }
    }
  }


  double effBTagEventSF = effBTagEvent_scaled / effBTagEvent_unscaled;
  // Catch inf and nan
  if (effBTagEventSF>9999. || effBTagEventSF!=effBTagEventSF) effBTagEventSF = 0;

  if(verbose_>=1) std::cout<<"effBTagEvent_unscaled= "<<effBTagEvent_unscaled
                    <<" effBTagEvent_scaled = " <<effBTagEvent_scaled
                    <<" effBTagEventSF ="       <<effBTagEventSF << std::endl;

  if(!noHistograms_) hists_.find("effBTagEventSF" )->second->Fill( effBTagEventSF );

//nataliia is working on this, atm set weight on 1... FIXME
effBTagEventSF=1.;

  std::auto_ptr<double> bTagSFEventWeight(new double);
  *bTagSFEventWeight = effBTagEventSF;    
  evt.put(bTagSFEventWeight);  
}

//--------------------------------------------------------------------------

/// Default SF values taken from database wrt. PAS BTV-11-001 (pTrel method),
/// or from PAS BTV-11-004
/// or from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript for 2012 ABCD
/// Values for eff. from user-defined histo as a function of pt and eta.

double BTagSFEventWeight::effBTagSF11004(double x)
{
  // function from PAS 11-004; x = jetPt
  if     (bTagAlgo_=="SSVHEM") return 0.896462*((1.+(0.00957275*x))/(1.+(0.00837582*x)));
  else if(bTagAlgo_=="CSVM")   return 0.6981*((1.+(0.414063*x))/(1.+(0.300155*x)));
  else if(bTagAlgo_=="JPM")    return 0.90806*((1.+(0.000236997*x))/(1.+(5.49455e-05*x)));
  else { 
    std::cout<< "WARNING!!! b tag SF for "<< bTagAlgo_ <<" not in code!!! CHECK!!!"<<std::endl;
    return 1.; 
  }
}

double BTagSFEventWeight::effBTagSF2012(double x)
{
  // function from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript; x = jetPt
  if     (bTagAlgo_=="CSVM") return (0.938887+(0.00017124*x))+(-2.76366e-07*(x*x));
  else if(bTagAlgo_=="CSVT") return (0.927563+(1.55479e-05*x))+(-1.90666e-07*(x*x));
  else { 
    std::cout<< "WARNING!!! b tag SF for "<< bTagAlgo_ << " for effBTagSF2012 not in code!!! CHECK!!!" << std::endl;
    return 1.; 
  }
}

double BTagSFEventWeight::effBTagSFerr11004(double x)
{
  // function from PAS 11-004; x = jetPt
  // pt binning
  double pt[] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 670};
  // corresponding SFb uncertainties
  double SFb_errorSSVHEM[] = {
    0.0316234,
    0.0310149,
    0.02381,
    0.0223228,
    0.023461,
    0.0202517,
    0.0156249,
    0.0214799,
    0.0399369,
    0.0416666,
    0.0431031,
    0.0663209,
    0.0687731,
    0.0793305 };
  double SFb_errorCSVM[] = {
    0.0295675,
    0.0295095,
    0.0210867,
    0.0219349,
    0.0227033,
    0.0204062,
    0.0185857,
    0.0256242,
    0.0383341,
    0.0409675,
    0.0420284,
    0.0541299,
    0.0578761,
    0.0655432 };
  double SFb_errorJPM[] = {
    0.0352594,
    0.0353008,
    0.0299008,
    0.0276606,
    0.0292312,
    0.0336607,
    0.0284701,
    0.029544,
    0.0358872,
    0.0367869,
    0.0375048,
    0.0597367,
    0.0653152,
    0.074242 };
  /// look for index corresponding to pt
  int iBin = -1;
  // if pt<30 use 12%
  if (x<pt[0]) return 0.12;
  for(int i=0; i<14; i++) {
    if (x>pt[i] && x<pt[i+1]) {
      iBin =i;
      break;
    }
  }
  double factor = 1.;
  if(iBin<0){
    // if pt>670 use SFb(670) and twice its error
    iBin=13;
    factor=2;
  }
  // FIXME for 8TeV: errors enlarged by factor 1.5
  factor*=1.5;
  if(bTagAlgo_=="SSVHEM") return factor * SFb_errorSSVHEM[iBin];
  if(bTagAlgo_=="CSVM")   return factor * SFb_errorCSVM[iBin];
  if(bTagAlgo_=="JPM")    return factor * SFb_errorJPM[iBin];
  else { 
    std::cout<< "WARNING!!! b tag SF for "<< bTagAlgo_ <<" not in code!!! CHECK!!!"<<std::endl;
    return 1.; 
  }
}

double BTagSFEventWeight::effBTagSFerr2012(double x)
{
  // function from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript (dataset=ABCD); x = jetPt
  // pt binning
  double pt[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
  // corresponding SFb uncertainties
  double SFb_errorCSVM[] = {
    0.0415707,
    0.0204209,
    0.0223227,
    0.0206655,
    0.0199325,
    0.0174121,
    0.0202332,
    0.0182446,
    0.0159777,
    0.0218531,
    0.0204688,
    0.0265191,
    0.0313175,
    0.0415417,
    0.0740446,
    0.0596716 };

  double SFb_errorCSVT[] = {
    0.0515703,
    0.0264008,
    0.0272757,
    0.0275565,
    0.0248745,
    0.0218456,
    0.0253845,
    0.0239588,
    0.0271791,
    0.0273912,
    0.0379822,
    0.0411624,
    0.0786307,
    0.0866832,
    0.0942053,
    0.102403 };

  /// look for index corresponding to pt
  int iBin = -1;
  for(int i=0; i<15; i++) {
    if (x>pt[i] && x<pt[i+1]) {
      iBin =i;
      break;
    }
  }
  double factor = 1.;
  if(iBin<0){
    // outside the quoted pt range: use twice the error
    factor=2;
    // if pt>800: use SFb(800)
    if(x>800) iBin=14;
    // if pt<20: use SFb(20)
    if(x<20 ) iBin=0;
  } 
  if     (bTagAlgo_=="CSVM")   return factor * SFb_errorCSVM[iBin];
  else if(bTagAlgo_=="CSVT")   return factor * SFb_errorCSVT[iBin];
  else { 
    std::cout<< "WARNING!!! b tag SF for "<< bTagAlgo_ << " in effBTagSFerr2012 not in code!!! CHECK!!!"<<std::endl;
    return 1.; 
  }
}
double BTagSFEventWeight::effBTagSF76X(double x)
{
  // function from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X; x = jetPt
  // pt binning
  double pt[] = {30, 60, 120, 320};
  // corresponding SFb uncertainties
  // corresponding SFb uncertainties
    double SFb_CSVL[] = {
    0.995,
    0.995,
    0.964};
  double SFb_CSVM[] = {
    0.966,
    0.971,
    0.927};

  double SFb_CSVT[] = {
    0.94,
    0.959,
    0.917};

  /// look for index corresponding to pt
  int iBin = -1;
  for(int i=0; i<3; i++) {
    if (x>pt[i] && x<pt[i+1]) {
      iBin =i;
      break;
    }
  }
//   double factor = 1.;
//   if(iBin<0){
//     // outside the quoted pt range: use twice the error
//     factor=2;
//     // if pt>800: use SFb(800)
//     if(x>800) iBin=14;
//     // if pt<20: use SFb(20)
//     if(x<20 ) iBin=0;
//   } 
  if     (bTagAlgo_=="CSVM")   return SFb_CSVM[iBin];
  else if(bTagAlgo_=="CSVT")   return SFb_CSVT[iBin];
  else if(bTagAlgo_=="CSVL")   return SFb_CSVL[iBin];
  
  else { 
    std::cout<< "WARNING!!! b tag SF for "<< bTagAlgo_ << " in effBTagSF76X not in code!!! CHECK!!!"<<std::endl;
    return 1.; 
  }
}
double BTagSFEventWeight::effBTagSFerr76X(double x)
{
  // function from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X; x = jetPt
  // pt binning
  double pt[] = {30, 60, 120, 320};
  // corresponding SFb uncertainties
    double SFb_errorCSVL[] = {
    0.054,
    0.033,
    0.044};
  double SFb_errorCSVM[] = {
    0.032,
    0.017,
    0.028};

  double SFb_errorCSVT[] = {
    0.037,
    0.020,
    0.032};

  /// look for index corresponding to pt
  int iBin = -1;
  for(int i=0; i<3; i++) {
    if (x>pt[i] && x<pt[i+1]) {
      iBin =i;
      break;
    }
  }
  double factor = 1.;
//   if(iBin<0){
//     // outside the quoted pt range: use twice the error
//     factor=2;
//     // if pt>800: use SFb(800)
//     if(x>800) iBin=14;
//     // if pt<20: use SFb(20)
//     if(x<20 ) iBin=0;
//   } 
  if     (bTagAlgo_=="CSVM")   return factor * SFb_errorCSVM[iBin];
  else if(bTagAlgo_=="CSVT")   return factor * SFb_errorCSVT[iBin];
  else if(bTagAlgo_=="CSVL")   return factor * SFb_errorCSVL[iBin];
  
  else { 
    std::cout<< "WARNING!!! b tag SF for "<< bTagAlgo_ << " in effBTagSFerr76X not in code!!! CHECK!!!"<<std::endl;
    return 1.; 
  }
}
double BTagSFEventWeight::effMisTagSF11004(double x, double jetEta, TString meanminmax)
{
  // function from PAS 11-004; x = jetPt
  // meanminmax = "mean" -> central value; = "min" -> down variation; = "max" -> up variation
  if(bTagAlgo_=="SSVHEM"){
    if(jetEta>=0. && jetEta <=0.8 && x< 670.){
      if( meanminmax == "mean" ) return  ((0.86318+(0.000801639*x))+(-1.64119e-06*(x*x)))+(2.59121e-10*(x*(x*x)));
      if( meanminmax == "min" )  return  ((0.790364+(0.000463086*x))+(-4.35934e-07*(x*x)))+(-9.08296e-10*(x*(x*x)));
      if( meanminmax == "max" )  return  ((0.935969+(0.0011402*x))+(-2.84645e-06*(x*x)))+(1.42654e-09*(x*(x*x)));
    }
    else if(jetEta>0.8 && jetEta <=1.6 && x< 670.){
      if( meanminmax == "mean" ) return  ((0.958973+(-0.000269555*x))+(1.381e-06*(x*x)))+(-1.87744e-09*(x*(x*x)));
      if( meanminmax == "min" )  return  ((0.865771+(-0.000279908*x))+(1.34144e-06*(x*x)))+(-1.75588e-09*(x*(x*x)));
      if( meanminmax == "max" )  return  ((1.0522+(-0.000259296*x))+(1.42056e-06*(x*x)))+(-1.999e-09*(x*(x*x)));
    }
    else if(jetEta>1.6 && jetEta <=2.4 && x< 670.){
      if( meanminmax == "mean" ) return  ((0.923033+(-0.000898227*x))+(4.74565e-06*(x*x)))+(-6.11053e-09*(x*(x*x)));
      if( meanminmax == "min" )  return  ((0.828021+(-0.000731926*x))+(4.19613e-06*(x*x)))+(-5.81379e-09*(x*(x*x)));
      if( meanminmax == "max" )  return  ((1.01812+(-0.00106483*x))+(5.29518e-06*(x*x)))+(-6.40728e-09*(x*(x*x)));
    }
    else if(jetEta>=0. && jetEta <=2.4 && x> 670.){
      x=670.;
      if( meanminmax == "mean" ) return  ((0.890254+(0.000553319*x))+(-1.29993e-06*(x*x)))+(4.19294e-10*(x*(x*x)));
      if( meanminmax == "min" )  return  ((0.817099+(0.000421567*x))+(-9.46432e-07*(x*x)))+(1.62339e-10*(x*(x*x)));
      if( meanminmax == "max" )  return  ((0.963387+(0.000685092*x))+(-1.65343e-06*(x*x)))+(6.76249e-10*(x*(x*x)));
    }
  }
  else if(bTagAlgo_=="CSVM"){
    double val=0;
    if(jetEta>=0. && jetEta <=0.8 && x< 670.){
      if( meanminmax == "mean")  val=  ((1.06182 +(0.000617034*x))+(-1.5732e-06* (x*x)))+( 3.02909e-10*(x*(x*x)));
      if( meanminmax == "min" )  val=  ((0.972455+(7.51396e-06*x))+( 4.91857e-07*(x*x)))+(-1.47661e-09*(x*(x*x)));
      if( meanminmax == "max" )  val=  ((1.15116 +(0.00122657 *x))+(-3.63826e-06*(x*x)))+( 2.08242e-09*(x*(x*x)));
    }
    else if(jetEta>0.8 && jetEta <=1.6 && x< 670.){
      if( meanminmax == "mean")  val=  ((1.111  +(-9.64191e-06*x))+( 1.80811e-07*(x*x)))+(-5.44868e-10*(x*(x*x)));
      if( meanminmax == "min" )  val=  ((1.02055+(-0.000378856*x))+( 1.49029e-06*(x*x)))+(-1.74966e-09*(x*(x*x)));
      if( meanminmax == "max" )  val=  ((1.20146+(0.000359543 *x))+(-1.12866e-06*(x*x)))+( 6.59918e-10*(x*(x*x)));
    }
    else if(jetEta>1.6 && jetEta <=2.4 && x< 670.){
      if( meanminmax == "mean")  val=  ((1.08498 +(-0.000701422*x))+(3.43612e-06*(x*x)))+(-4.11794e-09*(x*(x*x)));
      if( meanminmax == "min" )  val=  ((0.983476+(-0.000607242*x))+(3.17997e-06*(x*x)))+(-4.01242e-09*(x*(x*x)));
      if( meanminmax == "max" )  val=  ((1.18654 +(-0.000795808*x))+(3.69226e-06*(x*x)))+(-4.22347e-09*(x*(x*x)));
    }
    else if(jetEta>=0. && jetEta <=2.4 && x> 670.){
      x=670.;
      if( meanminmax == "mean")  val=  ((1.04318 +(0.000848162*x))+(-2.5795e-06* (x*x)))+(1.64156e-09*(x*(x*x)));
      if( meanminmax == "min" )  val=  ((0.962627+(0.000448344*x))+(-1.25579e-06*(x*x)))+(4.82283e-10*(x*(x*x)));
      if( meanminmax == "max" )  val=  ((1.12368 +(0.00124806 *x))+(-3.9032e-06* (x*x)))+(2.80083e-09*(x*(x*x)));
    }
    // FIXME: 8TeV correction factor: central value scaled, rel. uncertainty stays the same
    val*=(1.10422 + -0.000523856*x + 1.14251e-06*x*x);
    if(val!=0) return val;

  }
  else if(bTagAlgo_=="JPM"){
    if(jetEta>=0. && jetEta <=0.8 && x< 670.){
      if( meanminmax == "mean" ) return  ((0.970028+(0.00118179*x))+(-4.23119e-06*(x*x)))+(3.61065e-09*(x*(x*x)));
      if( meanminmax == "min" )  return  ((0.840326+(0.000626372*x))+(-2.08293e-06*(x*x)))+(1.57604e-09*(x*(x*x)));
      if( meanminmax == "max" )  return  ((1.09966+(0.00173739*x))+(-6.37946e-06*(x*x)))+(5.64527e-09*(x*(x*x)));
    }
    else if(jetEta>0.8 && jetEta <=1.6 && x< 670.){
      if( meanminmax == "mean" ) return  ((0.918387+(0.000898595*x))+(-2.00643e-06*(x*x)))+(1.26486e-09*(x*(x*x)));
      if( meanminmax == "min" )  return  ((0.790843+(0.000548016*x))+(-6.70941e-07*(x*x)))+(1.90355e-11*(x*(x*x)));
      if( meanminmax == "max" )  return  ((1.0459+(0.00124924*x))+(-3.34192e-06*(x*x)))+(2.51068e-09*(x*(x*x)));
    }
    else if(jetEta>1.6 && jetEta <=2.4 && x< 670.){
      if( meanminmax == "mean" ) return  ((0.790103+(0.00117865*x))+(-2.07334e-06*(x*x)))+(1.42608e-09*(x*(x*x)));
      if( meanminmax == "min" )  return  ((0.667144+(0.00105593*x))+(-1.43608e-06*(x*x)))+(5.24039e-10*(x*(x*x)));
      if( meanminmax == "max" )  return  ((0.913027+(0.00130143*x))+(-2.71061e-06*(x*x)))+(2.32812e-09*(x*(x*x)));
    }
    else if(jetEta>=0. && jetEta <=2.4 && x> 670.){
      x=670.;
      if( meanminmax == "mean" ) return  ((0.871294+(0.00215201*x))+(-6.77675e-06*(x*x)))+(5.79197e-09*(x*(x*x)));
      if( meanminmax == "min" )  return  ((0.7654+(0.00149792*x))+(-4.47192e-06*(x*x)))+(3.67664e-09*(x*(x*x)));
      if( meanminmax == "max" )  return  ((0.977076+(0.00280638*x))+(-9.08158e-06*(x*x)))+(7.9073e-09*(x*(x*x)));
    }
  }
  else { 
    std::cout<< "WARNING!!! b tag SF for "<< bTagAlgo_ <<" not in code!!! CHECK!!!"<<std::endl;
    return 1.; 
  }
  return -1111.;
}

double BTagSFEventWeight::effMisTagSF76X(double x, double jetEta, TString meanminmax)
{
  // function from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X; cMVAv2.csv
  // x = jetPt; meanminmax = "mean" -> central value; = "min" -> down variation; = "max" -> up variation
  if(bTagAlgo_=="CSVM"){
    double val =0, valMean=0, valMin=0, valMax=0;
    bool outOfRange = false;
    if(x<20) {x=20; outOfRange = true;}
    if(jetEta>=0. && jetEta <=0.8){
      if(x>1000.) {x=1000.; outOfRange = true;}
      valMean= (((1.10175+(-0.000724208*x))+(1.46283e-06*(x*x)))+(-8.31007e-10*(x*(x*x))));
      valMin= (((1.02873+(-0.000988431*x))+(2.01385e-06*(x*x)))+(-1.15193e-09*(x*(x*x))));
      valMax= (((1.17475+(-0.000458922*x))+(9.09541e-07*(x*x)))+(-5.09248e-10*(x*(x*x))));
    }
    else if(jetEta>0.8 && jetEta <=1.6){
      if(x>1000.) {x=1000.; outOfRange = true;}
      valMean= (((1.10358+(-0.000664145*x))+(1.05713e-06*(x*x)))+(-6.61479e-10*(x*(x*x))));
      valMin= (((1.034+(-0.000884165*x))+(1.571e-06*(x*x)))+(-9.71575e-10*(x*(x*x))));
      valMax= (((1.17314+(-0.000443333*x))+(5.41428e-07*(x*x)))+(-3.50681e-10*(x*(x*x))));
    }
    else if(jetEta>1.6 && jetEta <=2.4){
      if(x>1000.) {x=1000.; outOfRange = true;}
      valMean= (((1.01978+(-0.000154877*x))+(5.94704e-08*(x*x)))+(-1.45985e-10*(x*(x*x))));
      valMin= (((0.961861+(-0.000325026*x))+(4.83892e-07*(x*x)))+(-4.31829e-10*(x*(x*x))));
      valMax= (((1.07769+(1.53589e-05*x))+(-3.6544e-07*(x*x)))+(1.40114e-10*(x*(x*x))));
    }
    if( meanminmax == "mean") val= valMean;
    if( meanminmax == "min" ) val= outOfRange? valMin-(valMean-valMin) : valMin; // if outOfRange -> 2x the uncertainty
    if( meanminmax == "max" ) val= outOfRange? valMax+(valMax-valMean) : valMax; // if outOfRange -> 2x the uncertainty
    if(val!=0) return val;
  }
  else if(bTagAlgo_=="CSVT"){
    double val =0, valMean=0, valMin=0, valMax=0;
    bool outOfRange = false;
    if(x<20) {x=20; outOfRange = true;}
    if(jetEta>=0. && jetEta <=2.4){
      if(x>1000.) {x=1000.; outOfRange = true;}
      valMean= (((1.07557+(-0.000176198*x))+(-5.11449e-07*(x*x)))+(6.80495e-10*(x*(x*x))));
      valMin = (((0.922907+(-0.000413914*x))+(1.48537e-07*(x*x)))+(2.15697e-10*(x*(x*x))));
      valMax = (((1.22825+(6.09115e-05*x))+(-1.16904e-06*(x*x)))+(1.14407e-09*(x*(x*x))) );
    }
    if( meanminmax == "mean") val= valMean;
    if( meanminmax == "min" ) val= outOfRange? valMin-(valMean-valMin) : valMin; // if outOfRange -> 2x the uncertainty
    if( meanminmax == "max" ) val= outOfRange? valMax+(valMax-valMean) : valMax; // if outOfRange -> 2x the uncertainty
    if(val!=0) return val;
  }
  else if(bTagAlgo_=="CSVL"){
    double val =0, valMean=0, valMin=0, valMax=0;
    bool outOfRange = false;
    if(x<20) {x=20; outOfRange = true;}
    if(jetEta>=0. && jetEta <=0.3){
      if(x>1000.) {x=1000.; outOfRange = true;}
      valMean= (((0.874921+(0.000885525*x))+(-6.68503e-07*(x*x)))+(7.17689/x));
      valMin= (((0.842791+(0.000809607*x))+(-5.88337e-07*(x*x)))+(7.10941/x));
      valMax= (((0.903676+(0.000971152*x))+(-7.56441e-07*(x*x)))+(7.50067/x));
    }
    else if(jetEta>0.3 && jetEta <=0.6){
      if(x>1000.) {x=1000.; outOfRange = true;}
      valMean= (((0.859206+(0.000613401*x))+(-4.20717e-07*(x*x)))+(6.4676/x));
      valMin= (((0.829917+(0.0005496*x))+(-3.52954e-07*(x*x)))+(6.39268/x));
      valMax= (((0.885499+(0.000685967*x))+(-4.95495e-07*(x*x)))+(6.77182/x));
    }
    else if(jetEta>0.6 && jetEta <=0.9){
      if(x>1000.) {x=1000.; outOfRange = true;}
      valMean= (((0.85091+(0.00044326*x))+(-2.39676e-07*(x*x)))+(6.44336/x));
      valMin= (((0.822458+(0.000390438*x))+(-1.82749e-07*(x*x)))+(6.37768/x));
      valMax= (((0.876324+(0.000505268*x))+(-3.04015e-07*(x*x)))+(6.73902/x));
    }
    else if(jetEta>0.9 && jetEta <=1.2){
      if(x>1000.) {x=1000.; outOfRange = true;}
      valMean= (((0.884116+(0.000452952*x))+(-2.89454e-07*(x*x)))+(5.58513/x));
      valMin= (((0.854151+(0.000410262*x))+(-2.40923e-07*(x*x)))+(5.57/x));
      valMax= (((0.911327+(0.000504077*x))+(-3.44914e-07*(x*x)))+(5.80246/x));
    }
    else if(jetEta>1.2 && jetEta <=1.5){
      if(x>1000.) {x=1000.; outOfRange = true;}
      valMean= (((1.05409+(0.000190909*x))+(-2.65215e-07*(x*x)))+(-1.33829/x));
      valMin= (((1.01354+(0.000167333*x))+(-2.16554e-07*(x*x)))+(-0.922633/x));
      valMax= (((1.09533+(0.000212092*x))+(-3.11972e-07*(x*x)))+(-1.80293/x));
    }
    else if(jetEta>1.5 && jetEta <=1.8){
      if(x>1000.) {x=1000.; outOfRange = true;}
      valMean= (((1.02708+(0.000193904*x))+(-3.316e-07*(x*x)))+(5.82174e-11*(x*(x*x))));
      valMin= (((1.00578+(1.62714e-05*x))+(1.17942e-07*(x*x)))+(-2.35896e-10*(x*(x*x))));
      valMax= (((1.04838+(0.000371425*x))+(-7.81332e-07*(x*x)))+(3.52593e-10*(x*(x*x))));
    }
    else if(jetEta>1.8 && jetEta <=2.4){
      if(x>1000.) {x=1000.; outOfRange = true;}
      valMean= (((0.9986+(0.000458216*x))+(-8.50529e-07*(x*x)))+(3.64084e-10*(x*(x*x))));
      valMin= (((0.978148+(0.000286299*x))+(-3.80785e-07*(x*x)))+(2.90521e-11*(x*(x*x))));
      valMax= (((1.01905+(0.0006298*x))+(-1.32015e-06*(x*x)))+(6.99514e-10*(x*(x*x))));
    }
    if( meanminmax == "mean") val= valMean;
    if( meanminmax == "min" ) val= outOfRange? valMin-(valMean-valMin) : valMin; // if outOfRange -> 2x the uncertainty
    if( meanminmax == "max" ) val= outOfRange? valMax+(valMax-valMean) : valMax; // if outOfRange -> 2x the uncertainty
    if(val!=0) return val;
  }
  else { 
    std::cout<< "WARNING!!! b tag SF for "<< bTagAlgo_ << " in effMisTagSF76X not in code!!! CHECK!!!"<<std::endl;
    return 1.; 
  }
  return -1111.;
}

double BTagSFEventWeight::effMisTagSF2012(double x, double jetEta, TString meanminmax)
{
  // function from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript; dataset=ABCD 
  // x = jetPt; meanminmax = "mean" -> central value; = "min" -> down variation; = "max" -> up variation
  if(bTagAlgo_=="CSVM"){
    double val =0, valMean=0, valMin=0, valMax=0;
    bool outOfRange = false;
    if(x<20) {x=20; outOfRange = true;}
    if(jetEta>=0. && jetEta <=0.8){
      if(x>1000.) {x=1000.; outOfRange = true;}
      valMean= ((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x)));
      valMin= ((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x)));
      valMax= ((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x)));
    }
    else if(jetEta>0.8 && jetEta <=1.6){
      if(x>1000.) {x=1000.; outOfRange = true;}
      valMean= ((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)));
      valMin= ((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x)));
      valMax= ((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x)));
    }
    else if(jetEta>1.6 && jetEta <=2.4){
      if(x>850.) {x=850.; outOfRange = true;}
      valMean= ((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)));
      valMin= ((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x)));
      valMax= ((1.15575+(0.000693344*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x)));
    }
    if( meanminmax == "mean") val= valMean;
    if( meanminmax == "min" ) val= outOfRange? valMin-(valMean-valMin) : valMin; // if outOfRange -> 2x the uncertainty
    if( meanminmax == "max" ) val= outOfRange? valMax+(valMax-valMean) : valMax; // if outOfRange -> 2x the uncertainty
    if(val!=0) return val;
  }
  else if(bTagAlgo_=="CSVT"){
    double val =0, valMean=0, valMin=0, valMax=0;
    bool outOfRange = false;
    if(x<20) {x=20; outOfRange = true;}
    if(jetEta>=0. && jetEta <=2.4){
      if(x>1000.) {x=1000.; outOfRange = true;}
      valMean= ((1.00462+(0.00325971*x))+(-7.79184e-06*(x*x)))+(5.22506e-09*(x*(x*x)));
      valMin = ((0.845757+(0.00186422*x))+(-4.6133e-06*(x*x)))+(3.21723e-09*(x*(x*x)));
      valMax = ((1.16361+(0.00464695*x))+(-1.09467e-05*(x*x)))+(7.21896e-09*(x*(x*x)));
    }
    if( meanminmax == "mean") val= valMean;
    if( meanminmax == "min" ) val= outOfRange? valMin-(valMean-valMin) : valMin; // if outOfRange -> 2x the uncertainty
    if( meanminmax == "max" ) val= outOfRange? valMax+(valMax-valMean) : valMax; // if outOfRange -> 2x the uncertainty
    if(val!=0) return val;
  }
  else { 
    std::cout<< "WARNING!!! b tag SF for "<< bTagAlgo_ << " in effMisTagSF2012 not in code!!! CHECK!!!"<<std::endl;
    return 1.; 
  }
  return -1111.;
}

// b tag eff. from MC as a function of jet pt, eta
double BTagSFEventWeight::effBTag(double jetPt, double jetEta)
{
  double result = -1111.;
  // if histo file exists, take value from there; else return a default value
  if(filename_.location()) {
    TH2F* his = effHists_.find("EffBJetsTaggedPtEta")->second;
    // ensure that pt is in accepted range of BTV DB
    if(jetPt  >= maxPt2012_) jetPt  = maxPt2012_-1.;
    if(jetEta >= maxEta_   ) jetEta = maxEta_-0.1;
    result = his->GetBinContent( his->FindBin(jetPt, jetEta) );
  }
  else {result = 0.7; std::cout<< "WARNING!!! b tag eff. is ALWAYS 0.7!!! CHECK!!!"<<std::endl; }
  if(verbose_>=2) std::cout<< "effBTag= "<<result<<std::endl;
  return result;
}

// b tag eff. SF as a function of jet pt, eta
double BTagSFEventWeight::effBTagSF(double jetPt, double jetEta, bool isCjet)
{
  double result = -1111., error = -1111.;
  const BtagPerformance & perf = *(perfHBTag.product()); //FIXME
  BinningPointByMap measurePoint;
  if(version_=="DB11-001"){
    /// either take SF from BTV database...
      // ensure that pt is in accepted range
      if(jetPt >= maxPtDB_) jetPt = maxPtDB_-1.;
      if(jetEta >= maxEta_) jetEta = maxEta_-0.1;
      measurePoint.insert(BinningVariables::JetEt, jetPt);
      measurePoint.insert(BinningVariables::JetAbsEta, jetEta);      
      if(perf.isResultOk( measureMap_[ "BTAGBEFFCORR" ], measurePoint))
	result = perf.getResult( measureMap_[ "BTAGBEFFCORR" ], measurePoint);
      else {
	std::cout << "ERROR! B-tag SF could not be taken from DB! b-tag SF is taken as 1!" << std::endl;
	result = 1.;
      }
  }
  else if(version_=="11-004"){
    /// ...or by hand from 11-004 (Moriond recommendation)
    result = effBTagSF11004(jetPt);
  }
  else if(version_=="2012"){
    /// ...or by hand from 2012 (EPS 2013 recommendation)
    result = effBTagSF2012(jetPt);
  }
  else if(version_=="76X"){
    ///by hand from 76X https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X/cMVAv2.csv 	
    result = effBTagSF76X(jetPt);
  }
  if(uncertaintySFb_<0.){
    if(version_=="DB11-001"){
      /// either take SF from BTV database...
      if(perf.isResultOk( measureMap_[ "BTAGBERRCORR" ], measurePoint))
	error = perf.getResult( measureMap_[ "BTAGBERRCORR" ], measurePoint);
      else {
	std::cout << "ERROR! B-tag SF err could not be taken from DB! b-tag SF err is taken as 0.1!" << std::endl;
	error = 0.1;
      }
    }
    else if(version_=="11-004"){
      /// ...or by hand from 11-004 (Moriond recommendation)
      error = effBTagSFerr11004(jetPt);
    }
    else if(version_=="2012"){
      /// ...or by hand from 2012 (EPS 2013 recommendation)
      error = effBTagSFerr2012(jetPt);
    }
    else if(version_=="76X"){
      ///by hand from 76X https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X/cMVAv2.csv 
      error = effBTagSFerr76X(jetPt);
    }
  }
  else     error = uncertaintySFb_;
  // twice the error for c-jets
  if(isCjet) error*=2;

  /// different versions of sys. variations
  if(sysVar_ == "bTagSFUp")   result += error;
  else if(sysVar_ == "bTagSFDown") result -= error;
  else if(sysVar_ == "bTagSFShapeUpPt"){
    if(jetPt<shapeVarPtThreshold_) result += (shapeDistortionFactor_*error);
    else                           result -= (shapeDistortionFactor_*error);
  }
  else if(sysVar_ == "bTagSFShapeDownPt"){
    if(jetPt<shapeVarPtThreshold_) result -= (shapeDistortionFactor_*error);
    else                           result += (shapeDistortionFactor_*error);
  }
  else if(sysVar_ == "bTagSFShapeUpEta"){
    if(fabs(jetEta)<shapeVarEtaThreshold_) result += (shapeDistortionFactor_*error);
    else                                   result -= (shapeDistortionFactor_*error);
  }
  else if(sysVar_ == "bTagSFShapeDownEta"){
    if(fabs(jetEta)<shapeVarEtaThreshold_) result -= (shapeDistortionFactor_*error);
    else                                   result += (shapeDistortionFactor_*error);
  }
  if(verbose_>=2) std::cout<< "effBTagSF= "<<result<<" +/- "<<error<< "------ shapeDistortionFactor_=" << shapeDistortionFactor_ << "------ shapeDistortionFactor_*error=" << shapeDistortionFactor_*error <<std::endl;
  return result;
}

// b tag eff. from MC for c jets as a function of jet pt, eta;
// as first step: take average of b and mis eff.
double BTagSFEventWeight::effBTagCjet(double jetPt, double jetEta)
{
  double result = -1111.;
  // if histo file exists, take value from there; else return a default value
  if(filename_.location()) {
    TH2F* his = effHists_.find("EffCJetsTaggedPtEta")->second;
    // ensure that pt is in accepted range
    if(jetPt  >= maxPt2012_) jetPt  = maxPt2012_-1.;
    if(jetEta >= maxEta_   ) jetEta = maxEta_-0.1;
    result = his->GetBinContent( his->FindBin(jetPt, jetEta) );
  }
  else {result = 0.35; std::cout<< "WARNING!!! c tag eff. is ALWAYS 0.35!!! CHECK!!!"<<std::endl; }
  if(verbose_>=2) std::cout<< "effBTagCjet= "<<result<<std::endl;
  return result;
}

// mistag eff. from MC as a function of jet pt, eta
double BTagSFEventWeight::effMisTag(double jetPt, double jetEta)
{
  double result = -1111.;
  // if histo file exists, take value from there; else return a default value
  if(filename_.location()) {
    TH2F* his = effHists_.find("EffLJetsTaggedPtEta")->second;
    // ensure that pt is in accepted range
    if(jetPt >= maxPtMisTag_) jetPt = maxPtMisTag_-1.;
    if(jetEta >= maxEta_) jetEta = maxEta_-0.1;
    result = his->GetBinContent( his->FindBin(jetPt, jetEta) );
  }
  else {result = 0.01; std::cout<< "WARNING!!! b tag eff. is ALWAYS 0.01!!! CHECK!!!"<<std::endl; }
  if(verbose_>=2) std::cout<< "effMisTag= "<<result<<std::endl;
  return result;
}

// mistag eff. SF as a function of jet pt, eta
double BTagSFEventWeight::effMisTagSF(double jetPt, double jetEta)
{
  double result = -1111., error = -1111.;
  const BtagPerformance & perf = *(perfHMisTag.product()); //FIXME
  BinningPointByMap measurePoint;
  if(version_=="DB11-001"){
    /// either take SF from BTV database...
    if(jetPt >= maxPtMisTag_) jetPt = maxPtMisTag_-1.;
    if(jetEta >= maxEta_) jetEta = maxEta_-0.1;
    measurePoint.insert(BinningVariables::JetEt, jetPt);
    measurePoint.insert(BinningVariables::JetAbsEta, jetEta);
    if(perf.isResultOk( measureMap_[ "BTAGLEFFCORR" ], measurePoint))
	result = perf.getResult( measureMap_[ "BTAGLEFFCORR" ], measurePoint);
    else result = 1.;
    if(perf.isResultOk( measureMap_[ "BTAGLERRCORR" ], measurePoint))
	error = perf.getResult( measureMap_[ "BTAGLERRCORR" ], measurePoint);
    else error = 0.1;
    if(sysVar_ == "misTagSFUp")   result += error;
    if(sysVar_ == "misTagSFDown") result -= error;
    if(verbose_>=2) std::cout<< "effMisTagSF= "<<result<<" +/- "<<error<<std::endl;
  }
  else if(version_=="11-004"){
    /// ...or by hand from 11-004 (Moriond recommendation)
    if(jetEta >= maxEta_) jetEta = maxEta_-0.1;
    TString                       meanminmax = "mean";
    if(sysVar_ == "misTagSFUp"  ) meanminmax = "max";
    if(sysVar_ == "misTagSFDown") meanminmax = "min";
    result = effMisTagSF11004(jetPt, jetEta, meanminmax);
    if(verbose_>=2) std::cout<< "effMisTagSF= "<<effMisTagSF11004(jetPt, jetEta, "mean")<<" + "<<effMisTagSF11004(jetPt, jetEta, "max")
	  <<" - "<<effMisTagSF11004(jetPt, jetEta, "min")<<std::endl;
  }
  else if(version_=="2012"){
    /// ...or by hand from 11-004 (Moriond recommendation)
    if(jetEta >= maxEta_) jetEta = maxEta_-0.1;
    TString                       meanminmax = "mean";
    if(sysVar_ == "misTagSFUp"  ) meanminmax = "max";
    if(sysVar_ == "misTagSFDown") meanminmax = "min";
    result = effMisTagSF2012(jetPt, jetEta, meanminmax);
    if(verbose_>=2) std::cout<< "effMisTagSF= "<<effMisTagSF2012(jetPt, jetEta, "mean")<<" + "<<effMisTagSF2012(jetPt, jetEta, "max")
	  <<" - "<<effMisTagSF2012(jetPt, jetEta, "min")<<std::endl;
  }
  else if(version_=="76X"){
    /// ...or by hand from 11-004 (Moriond recommendation)
    if(jetEta >= maxEta_) jetEta = maxEta_-0.1;
    TString                       meanminmax = "mean";
    if(sysVar_ == "misTagSFUp"  ) meanminmax = "max";
    if(sysVar_ == "misTagSFDown") meanminmax = "min";
    result = effMisTagSF76X(jetPt, jetEta, meanminmax);
    if(verbose_>=2) std::cout<< "effMisTagSF= "<<effMisTagSF76X(jetPt, jetEta, "mean")<<" + "<<effMisTagSF76X(jetPt, jetEta, "max")
	  <<" - "<<effMisTagSF76X(jetPt, jetEta, "min")<<std::endl;
  }
  return result;
}

//--------------------------------------------------------------------------

// calculate event b tag efficiency for >=2 b tags
double BTagSFEventWeight::effBTagEvent(std::vector<double> &oneMinusBEffies,
				       std::vector<double> &oneMinusBMistags)
{
  double bTaggingEfficiency = 1.;
  double tmp = 1.;

  if(verbose_) std::cout << oneMinusBEffies.size() << ": " << std::flush;

  for(std::vector<double>::const_iterator eff = oneMinusBEffies.begin(); eff != oneMinusBEffies.end(); ++eff){
    tmp *= (*eff);
    if(verbose_) std::cout << 1.-(*eff) << ", ";
  }
  if(verbose_) std::cout << oneMinusBMistags.size() << ": " << std::flush;
  for(std::vector<double>::const_iterator mis = oneMinusBMistags.begin(); mis != oneMinusBMistags.end(); ++mis){
    tmp *= (*mis);
    if(verbose_) std::cout << 1.-(*mis) << ", ";
  }
  bTaggingEfficiency -= tmp;
  for(std::vector<double>::const_iterator eff = oneMinusBEffies.begin(); eff != oneMinusBEffies.end(); ++eff){
    tmp = 1.-(*eff);
    for(std::vector<double>::const_iterator eff2 = oneMinusBEffies.begin(); eff2 != oneMinusBEffies.end(); ++eff2){
      if(eff != eff2) tmp *= (*eff2);
    }
    for(std::vector<double>::const_iterator mis = oneMinusBMistags.begin(); mis != oneMinusBMistags.end(); ++mis){
      tmp *= (*mis);
    }
    bTaggingEfficiency -= tmp;
  }
  for(std::vector<double>::const_iterator mis = oneMinusBMistags.begin(); mis != oneMinusBMistags.end(); ++mis){
    tmp = 1.-(*mis);
    for(std::vector<double>::const_iterator eff = oneMinusBEffies.begin(); eff != oneMinusBEffies.end(); ++eff){
      tmp *= (*eff);
    }
    for(std::vector<double>::const_iterator mis2 = oneMinusBMistags.begin(); mis2 != oneMinusBMistags.end(); ++mis2){
      if(mis != mis2) tmp *= (*mis2);
    }
    bTaggingEfficiency -= tmp;
  }
  if(verbose_) std::cout << " -> " << bTaggingEfficiency << std::endl;
  return bTaggingEfficiency;
}

// executed at the end after looping over all events
void BTagSFEventWeight::endJob() 
{
  if(!noHistograms_) {
    double effBTagEventSFMean = hists_.find("effBTagEventSF" )->second->GetMean();
    hists_.find("effBTagEventSFMean" )->second->Fill(0.5, effBTagEventSFMean );
    if(verbose_>=1) std::cout<<"Mean effBTagEventSF = "<<effBTagEventSFMean<<std::endl;
  }
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( BTagSFEventWeight );

