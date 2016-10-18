#include "FWCore/Utilities/interface/EDMException.h"
#include "TopAnalysis/TopUtils/plugins/BTagSFEventWeight.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"

//custom headers
#include "BTagCalibrationStandalone.h"

#include <stdexcept>
// std::string cmsswBase = (getenv ("CMSSW_BASE"));
// setup calibration reader
// // BTagCalibration calib("CSVv2",cmsswBase+"/src/TopAnalysis/TopUtils/plugins/CSVv2.csv");
// BTagCalibrationReader reader(&calib,BTagEntry::OP_MEDIUM,"mujets","central");// "mujets" for b and c jets 
// BTagCalibrationReader reader_up(&calib,BTagEntry::OP_MEDIUM,"mujets","up");// "mujets" for b and c jets 
// BTagCalibrationReader reader_down(&calib,BTagEntry::OP_MEDIUM,"mujets","down");// "mujets" for b and c jets 
// 
// BTagCalibrationReader reader_incl(&calib,BTagEntry::OP_MEDIUM,"incl","central");// "incl" for light jets
// BTagCalibrationReader reader_incl_up(&calib,BTagEntry::OP_MEDIUM,"incl","up");// "incl" for light jets
// BTagCalibrationReader reader_incl_down(&calib,BTagEntry::OP_MEDIUM,"incl","down");// "incl" for light jets


BTagSFEventWeight::BTagSFEventWeight(const edm::ParameterSet& cfg):
  jets_                   ( cfg.getParameter<edm::InputTag>    ( "jets"   ) ),
  version_                ( cfg.getParameter<std::string>      ("version"  ) ),
  bTagAlgo_               ( cfg.getParameter<std::string>      ("bTagAlgo") ),
  discrCut_               ( cfg.getParameter<double>           ("discrCut" ) ),
  newRecipe_              ( cfg.getParameter<bool>             ("newRecipe") ),
  maxJets_                ( cfg.getParameter<int>              ("maxJets" ) ),
  sysVar_                 ( cfg.getParameter<std::string>      ("sysVar"  ) ),
  verbose_                ( cfg.getParameter<int>              ("verbose" ) ),
  filename_               ( cfg.getParameter<edm::FileInPath>  ("filename") ),
  csv_filename_           ( cfg.getParameter<edm::FileInPath>  ("csv_filename") ),
  noHistograms_           ( cfg.getParameter<bool>             ("noHistograms"))
{
  // setup calibration readers 
  BTagCalibration calib(version_ ,csv_filename_.fullPath());
  if ( (version_ == "CSVv2") || (version_ == "JP")){
    reader = BTagCalibrationReader(&calib,BTagEntry::OP_MEDIUM,"comb","central");// "mujets" for b and c jets 
    reader_up = BTagCalibrationReader(&calib,BTagEntry::OP_MEDIUM,"comb","up");// "mujets" for b and c jets 
    reader_down = BTagCalibrationReader(&calib,BTagEntry::OP_MEDIUM,"comb","down");// "mujets" for b and c jets 
  } else if (version_ == "cMVAv2") {
    reader = BTagCalibrationReader(&calib,BTagEntry::OP_MEDIUM,"ttbar","central");// "ttbar" for b and c jets 
    reader_up = BTagCalibrationReader(&calib,BTagEntry::OP_MEDIUM,"ttbar","up");// "ttbar" for b and c jets 
    reader_down = BTagCalibrationReader(&calib,BTagEntry::OP_MEDIUM,"ttbar","down");// "ttbar" for b and c jets 
  }   
  reader_incl = BTagCalibrationReader(&calib,BTagEntry::OP_MEDIUM,"incl","central");// "incl" for light jets
  reader_incl_up = BTagCalibrationReader(&calib,BTagEntry::OP_MEDIUM,"incl","up");// "incl" for light jets
  reader_incl_down = BTagCalibrationReader(&calib,BTagEntry::OP_MEDIUM,"incl","down");// "incl" for light jets

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
  std::cout<<csv_filename_.fullPath()<<" opened"<<std::endl;
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
	  effHists_.count("NumLJetsPtEta") && effHists_.count("NumLJetsTaggedPtEta") && effHists_.count("EffLJetsTaggedPtEta")) {
	  
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
  edm::Handle<edm::View< pat::Jet > > jets;
  evt.getByLabel(jets_, jets);
                             
  double pt, eta;                    
  std::vector<double> oneMinusBEffies(0) , oneMinusBEffies_scaled(0);
  std::vector<double> oneMinusBMistags(0), oneMinusBMistags_scaled(0);
  
  double effBTagEvent_unscaled = 1.;                        
  double effBTagEvent_scaled   = 1.;

  for(edm::View<pat::Jet>::const_iterator jet = jets->begin();jet != jets->end(); ++jet) {
    if (jet - jets->begin() == maxJets_) break; //!????????????????
    pt  = jet->pt();
    eta = std::abs(jet->eta());
    int flavor = std::abs(jet->partonFlavour());
    
    // tagged jets
    if (jet->bDiscriminator(bTagAlgo_)>discrCut_) {
      if (flavor == 5) {//b
	effBTagEvent_unscaled *= effBTag(pt, eta);
	effBTagEvent_scaled   *= effBTag(pt, eta) * effBTagSF(pt,eta);
  //       std::cout << " SF = "<<   reader.eval(BTagEntry::FLAV_B,eta,pt)<<std::endl;   
      }
      else if (flavor == 4) {//c
	effBTagEvent_unscaled *= effBTagCjet(pt, eta);
	effBTagEvent_scaled   *= effBTagCjet(pt, eta) * effBTagCjetSF(pt,eta);
      }
      else {
	effBTagEvent_unscaled *= effMisTag(pt, eta);
	effBTagEvent_scaled   *= effMisTag(pt, eta) * effMisTagSF(pt,eta);
      }
    }
    // untagged jets
    else {
      if (flavor == 5) {
	effBTagEvent_unscaled *= 1.- effBTag(pt, eta);
	effBTagEvent_scaled   *= 1.-(effBTag(pt, eta) * effBTagSF(pt,eta));
      }
      else if (flavor == 4) {
	effBTagEvent_unscaled *= 1.- effBTagCjet(pt, eta);
	effBTagEvent_scaled   *= 1.-(effBTagCjet(pt, eta) *  effBTagCjetSF(pt,eta));
      }
      else {
	effBTagEvent_unscaled *= 1.- effMisTag(pt, eta);
	effBTagEvent_scaled   *= 1.-(effMisTag(pt, eta) * effMisTagSF(pt,eta));
      }
    }
  }
  
  double effBTagEventSF = effBTagEvent_scaled / effBTagEvent_unscaled;
//   std::cout<<"effBTagEvent_unscaled= "<<effBTagEvent_unscaled
//                     <<" effBTagEvent_scaled = " <<effBTagEvent_scaled
//                     <<" effBTagEventSF ="       <<effBTagEventSF << std::endl;

  // Catch inf and nan
//   std::cout<<"effBTagEventSF  "<<effBTagEventSF<<std::endl;
  if (effBTagEventSF>9999. || effBTagEventSF!=effBTagEventSF) effBTagEventSF = 0;

  if(verbose_>=1) std::cout<<"effBTagEvent_unscaled= "<<effBTagEvent_unscaled
                    <<" effBTagEvent_scaled = " <<effBTagEvent_scaled
                    <<" effBTagEventSF ="       <<effBTagEventSF << std::endl;

  if(!noHistograms_) hists_.find("effBTagEventSF" )->second->Fill( effBTagEventSF );

  std::auto_ptr<double> bTagSFEventWeight(new double);
  *bTagSFEventWeight = effBTagEventSF;    

  evt.put(bTagSFEventWeight);  
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
double BTagSFEventWeight::effBTagSF(double jetPt, double jetEta)
{
  double result = -1111.;
  // ensure that pt is in accepted range of BTV DB
  if(jetPt  >= maxPt2012_) jetPt  = maxPt2012_-1.;
  if(jetEta >= maxEta_   ) jetEta = maxEta_-0.1;
  result = reader.eval(BTagEntry::FLAV_B,jetEta,jetPt);
  
  /// different versions of sys. variations
  if(sysVar_ == "bTagSFUp")   result = reader_up.eval(BTagEntry::FLAV_B,jetEta,jetPt);
  else if(sysVar_ == "bTagSFDown") result = reader_down.eval(BTagEntry::FLAV_B,jetEta,jetPt);
 
  if(verbose_>=2) std::cout<< "effBTagSF= "<<result <<std::endl;
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

double BTagSFEventWeight::effBTagCjetSF(double jetPt, double jetEta)
{
  double result = -1111.;
  // ensure that pt is in accepted range of BTV DB
  if(jetPt  >= maxPt2012_) jetPt  = maxPt2012_-1.;
  if(jetEta >= maxEta_   ) jetEta = maxEta_-0.1;
  result = reader.eval(BTagEntry::FLAV_C,jetEta,jetPt);
  
  /// different versions of sys. variations
  if(sysVar_ == "bTagCjetSFUp")   result = reader_up.eval(BTagEntry::FLAV_C,jetEta,jetPt);
  else if(sysVar_ == "bTagCjetSFDown") result = reader_down.eval(BTagEntry::FLAV_C,jetEta,jetPt);
 
  if(verbose_>=2) std::cout<< "effBTagCjetSF= "<<result <<std::endl;
  
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
  double result = -1111.;
  // ensure that pt is in accepted range of BTV DB
  if(jetPt  >= maxPt2012_) jetPt  = maxPt2012_-1.;
  if(jetEta >= maxEta_   ) jetEta = maxEta_-0.1;
  result = reader_incl.eval(BTagEntry::FLAV_UDSG,jetEta,jetPt);
  
  /// different versions of sys. variations
  if(sysVar_ == "misTagSFUp")   result = reader_incl_up.eval(BTagEntry::FLAV_UDSG,jetEta,jetPt);
  else if(sysVar_ == "misTagSFDown") result = reader_incl_down.eval(BTagEntry::FLAV_UDSG,jetEta,jetPt);
 
  if(verbose_>=2) std::cout<< "effBTagCjetSF= "<<result <<std::endl;
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

