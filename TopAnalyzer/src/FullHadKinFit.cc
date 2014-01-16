#include "TFile.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "TopAnalysis/TopAnalyzer/interface/FullHadKinFit.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TFitConstraintMGaus.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"

/// default constructor
FullHadKinFit::FullHadKinFit(const edm::ParameterSet& cfg)
{

  gRandom->SetSeed(0);
  rnd.SetSeed(0);

  src_          = cfg.getParameter<edm::InputTag>("src");
  smearOnly_    = cfg.getParameter<int>("smearOnly");
  smear_        = cfg.getParameter<double>("smear");
  smearType_    = cfg.getParameter<std::string>("smearType");
  resol_        = cfg.getParameter<double>("resol");
  resolType_    = cfg.getParameter<std::string>("resolType");
  bTags_        = cfg.getParameter<unsigned int>("bTags");
  onlyGenMatch_ = cfg.getParameter<bool>("onlyGenMatch");
  saveAllCombo_ = cfg.getParameter<bool>("saveAllCombo");
  exchangePart_ = cfg.getParameter<bool>("exchangePart");
  fitAsFitInpt_ = cfg.getParameter<bool>("fitAsFitInpt");
  noRadiations_ = cfg.getParameter<bool>("noRadiations");
  recoPtCut_    = cfg.getParameter<double>("recoPtCut");
  genPtCut_     = cfg.getParameter<double>("genPtCut");

  if(resolType_ != "gaussian") {
    if(cfg.exists("udscResolutions") && cfg.exists("bResolutions")){
      std::vector<edm::ParameterSet> functionSets_ = cfg.getParameter <std::vector<edm::ParameterSet> >("udscResolutions");
      for(std::vector<edm::ParameterSet>::const_iterator iSet = functionSets_.begin(); iSet != functionSets_.end(); ++iSet){
	if(iSet->exists("bin")) binsUdsc_.push_back(iSet->getParameter<std::string>("bin"));
	else if(functionSets_.size()==1) binsUdsc_.push_back("");
	else throw cms::Exception("WrongConfig") << "Parameter 'bin' is needed if more than one PSet is specified!\n";
      
	funcEtUdsc_.push_back(iSet->getParameter<std::string>("et"));
	funcEtaUdsc_.push_back(iSet->getParameter<std::string>("eta"));
	funcPhiUdsc_.push_back(iSet->getParameter<std::string>("phi"));
      }
      functionSets_.clear();
      functionSets_ = cfg.getParameter <std::vector<edm::ParameterSet> >("bResolutions");
      for(std::vector<edm::ParameterSet>::const_iterator iSet = functionSets_.begin(); iSet != functionSets_.end(); ++iSet){
	if(iSet->exists("bin")) binsB_.push_back(iSet->getParameter<std::string>("bin"));
	else if(functionSets_.size()==1) binsB_.push_back("");
	else throw cms::Exception("WrongConfig") << "Parameter 'bin' is needed if more than one PSet is specified!\n";
      
	funcEtB_.push_back(iSet->getParameter<std::string>("et"));
	funcEtaB_.push_back(iSet->getParameter<std::string>("eta"));
	funcPhiB_.push_back(iSet->getParameter<std::string>("phi"));
      }
    }
    else if(cfg.exists("udscResolutions") || cfg.exists("bResolutions")){
      if(cfg.exists("udscResolutions")) throw cms::Exception("WrongConfig") << "Parameter 'bResolutions' is needed if parameter 'udscResolutions' is defined!\n";
      else                              throw cms::Exception("WrongConfig") << "Parameter 'udscResolutions' is needed if parameter 'bResolutions' is defined!\n";
    }
    else {
      edm::LogWarning("Config") << "Parameters 'udscResolutions' and 'bResolutions' are not defined, taking standard resolutions\n";
      binsUdsc_.push_back("1.4<=abs(eta)");
      funcEtUdsc_.push_back("4.8^2 + (0.89^2)*et + (0.043^2)*(et)^2");
      funcEtaUdsc_.push_back("1.773/((et)^2) + 0.034/et + 0.000356");
      funcPhiUdsc_.push_back("2.908/((et)^2) + 0.021/et + 0.000259");
      binsUdsc_.push_back("1.4>abs(eta)");
      funcEtUdsc_.push_back("5.6^2 + (1.25^2)*et + (0.033^2)*(et)^2");
      funcEtaUdsc_.push_back("1.215/((et)^2) + 0.037/et + 0.0007941");
      funcPhiUdsc_.push_back("6.65/((et)^2) + 0.04/et + 0.0000849");
      binsB_.push_back("1.4<=abs(eta)");
      funcEtB_.push_back("4.8^2 + (0.89^2)*et + (0.043^2)*(et)^2");
      funcEtaB_.push_back("1.773/((et)^2) + 0.034/et + 0.000356");
      funcPhiB_.push_back("2.908/((et)^2) + 0.021/et + 0.000259");
      binsB_.push_back("1.4>abs(eta)");
      funcEtB_.push_back("5.6^2 + (1.25^2)*et + (0.033^2)*(et)^2");
      funcEtaB_.push_back("1.215/((et)^2) + 0.037/et + 0.0007941");
      funcPhiB_.push_back("6.65/((et)^2) + 0.04/et + 0.0000849");
    }
  }
}

/// do generator matching
void
FullHadKinFit::genMatcher(const edm::View<reco::LeafCandidate>& src)
{
  std::vector<const reco::GenParticle*> parts = dynCastVector<reco::GenParticle>( src );
  //std::vector<const reco::GenParticle*> parts(src.size());
  //std::transform(src.begin(),src.end(),parts.begin(),F<reco::GenParticle>());
  //parts.erase(std::remove_if(parts.begin(), parts.end(),std::bind2nd(std::equal_to<const reco::GenParticle*>(), 0)),parts.end());

  if(parts.size()){

    genParticles.resize(6);
    
    for(std::vector<const reco::GenParticle*>::const_iterator obj = parts.begin(); obj != parts.end(); ++obj){
      
      if((*obj)->status() == 3){
	
	if((*obj)->mother()){
	  if((*obj)->mother()->mother()){
	    if((*obj)->pdgId() ==   5 && (*obj)->mother()->mother()->pdgId() ==  6) genParticles[B]    = *obj;
	    if((*obj)->pdgId() ==  -5 && (*obj)->mother()->mother()->pdgId() == -6) genParticles[BBar] = *obj;
	  }
	  if((*obj)->pdgId() < 5 && (*obj)->pdgId() > -5 && (*obj)->pdgId() != 0) {
	    if((*obj)->pdgId() > 0 && (*obj)->mother()->pdgId() ==  24) genParticles[LQ]    = *obj;
	    if((*obj)->pdgId() < 0 && (*obj)->mother()->pdgId() ==  24) genParticles[LQBar] = *obj;
	    if((*obj)->pdgId() > 0 && (*obj)->mother()->pdgId() == -24) genParticles[LP]    = *obj;
	    if((*obj)->pdgId() < 0 && (*obj)->mother()->pdgId() == -24) genParticles[LPBar] = *obj;
	  }
	}
      }
    }
  }
}

/// find ME level partons not belonging to top decay directly
int
FullHadKinFit::foundRadiation(const edm::View<reco::LeafCandidate>& src, std::vector<const reco::GenParticle*>& genParts)
{
  std::vector<const reco::GenParticle*> parts = dynCastVector<reco::GenParticle>( src );

  int radiations = 0;

  for(std::vector<const reco::GenParticle*>::const_iterator obj = parts.begin(); obj != parts.end(); ++obj){
    if((*obj)->status() == 3){

      bool found = false;
      std::vector<const reco::GenParticle*>::const_iterator lowestPtParton = genParts.begin();
      //std::cout  << "PDG-ID: " << (*obj)->pdgId() << "; pt: " << (*obj)->pt() << "; eta: " << (*obj)->eta() << "; phi: " << (*obj)->phi() << std::endl;
      for(std::vector<const reco::GenParticle*>::const_iterator genPart = genParts.begin(); genPart != genParts.end(); ++genPart){
	//std::cout << genPart - genParts.begin() << std::endl;
	//std::cout << "PDG-ID: " << (*genPart)->pdgId() << "; pt: " << (*genPart)->pt() << "; eta: " << (*genPart)->eta() "; phi: " << (*genPart)->phi() << std::endl;
	if((*genPart)->pt() < (*lowestPtParton)->pt()) lowestPtParton = genPart;
	if((*obj)->pt() == (*genPart)->pt() && (*obj)->eta() == (*genPart)->eta() && (*obj)->phi() == (*genPart)->phi()) found = true;
      }
      if(!found && (*obj)->pdgId() != 6 && (*obj)->pdgId() != -6 && (*obj)->pdgId() != 24 && (*obj)->pdgId() != -24 && (*obj)->pdgId() != 2212 && (*obj)->pt() != 0){
	//std::cout  << "PDG-ID: " << (*obj)->pdgId() << "; pt: " << (*obj)->pt() << "; eta: " << (*obj)->eta() << "; phi: " << (*obj)->phi();
	//if((*obj)->mother()){
	//  //std::cout << "; mother: " << (*obj)->mother()->pdgId();
	//  //if((*obj)->mother()->mother()){
	//  //  //std::cout << "; grandmother: " << (*obj)->mother()->mother()->pdgId();
	//  //}
	//}
	//std::cout << std::endl;

	++radiations;

	if(exchangePart_) genParts[lowestPtParton - genParts.begin()] = (*obj);
      }
    }
  }
  return radiations;
}

/// smear input vectors for kinematic fit
void
FullHadKinFit::smear(std::vector< TLorentzVector >& vecs){
  Double_t phi, pt, eta;//, e;
  for(std::vector< TLorentzVector >::iterator vec = vecs.begin(); vec != vecs.end(); ++vec){
    if( smearOnly_ == -1 || smearOnly_ == vec - vecs.begin() ){
      double pt_  = vec->Pt();
      double eta_ = vec->Eta();
      double phi_ = vec->Phi();
      if(smearType_ == "gaussian"){
	pt = -1.0;
	while( pt < 0 ){//|| fabs(pt_-pt) < 3.0 * smear_ * 0.3 * sqrt( pt_ ) ){
	  pt = rnd.Gaus( pt_, smear_ * 0.3 * sqrt( pt_ ) );
	}			     		 
	eta = rnd.Gaus( eta_, smear_ * 0.5 / sqrt( pt_ ) );
	phi = rnd.Gaus( phi_, smear_ * 0.5 / sqrt( pt_ ) );
      }
      else{
	TH1D* tempHist = dynamic_cast<TH2F*>(hists_.find("bResolution_pt")->second)->ProjectionY("",0,-1,"o");
	Int_t iBin = tempHist->FindBin(pt_);
	delete tempHist;
	if(vec - vecs.begin() == 0 || vec - vecs.begin() == 3){
	  pt = -1.0;
	  while( pt < 0 ){
	    TH1D* tempHist = dynamic_cast<TH2F*>(hists_.find("bResolution_pt")->second)->ProjectionX("",iBin,iBin,"o");
	    (tempHist->GetEntries() != 0) ? pt = tempHist->GetRandom() * pt_ * smear_  + pt_ : pt = pt_;
	    //std::cout << " -> pt(B): " << pt << std::endl;
	    delete tempHist;
	  }
	  TH1D* tempHist = dynamic_cast<TH2F*>(hists_.find("bResolution_eta")->second)->ProjectionX("",iBin,iBin,"o");
	  (tempHist->GetEntries() != 0) ? eta = tempHist->GetRandom() * smear_  + eta_ : eta = eta_;
	  delete tempHist;
	  tempHist = dynamic_cast<TH2F*>(hists_.find("bResolution_phi")->second)->ProjectionX("",iBin,iBin,"o");
	  (tempHist->GetEntries() != 0) ? phi = tempHist->GetRandom() * smear_  + phi_ : phi = phi_;
	  delete tempHist;
	}
	else{
	  pt = -1.0;
	  while( pt < 0 ){
	    TH1D* tempHist = dynamic_cast<TH2F*>(hists_.find("udscResolution_pt")->second)->ProjectionX("",iBin,iBin,"o");
	    (tempHist->GetEntries() != 0) ? pt = tempHist->GetRandom() * pt_ * smear_  + pt_ : pt = pt_;
	    //std::cout << " -> pt(L): " << pt << std::endl;
	    delete tempHist;
	  }
	  TH1D* tempHist = dynamic_cast<TH2F*>(hists_.find("udscResolution_eta")->second)->ProjectionX("",iBin,iBin,"o");
	  (tempHist->GetEntries() != 0) ? eta = tempHist->GetRandom() * smear_  + eta_ : eta = eta_;
	  delete tempHist;
	  tempHist = dynamic_cast<TH2F*>(hists_.find("udscResolution_phi")->second)->ProjectionX("",iBin,iBin,"o");
	  (tempHist->GetEntries() != 0) ? phi = tempHist->GetRandom() * smear_  + phi_ : phi = phi_;
	  delete tempHist;
	}
      }
      if( phi >= TMath::Pi() ) phi -= 2* TMath::Pi();
      if( phi < -TMath::Pi() ) phi += 2* TMath::Pi();
      //e = vec->E() * ( pt / pt_ );
      //std::cout << "-------------------" << std::endl;
      //std::cout << "pt:  " << vec->Pt()  << " +- " << 0.5 * sqrt( vec->Pt() )                          << " -> " << pt  << std::endl;
      //std::cout << "eta: " << vec->Eta() << " +- " << 0.2 * vec->Eta() * sqrt( vec->Pt() ) / vec->Pt() << " -> " << eta << std::endl;
      //std::cout << "phi: " << vec->Phi() << " +- " << 0.1                                              << " -> " << phi << std::endl;
      //std::cout << "e:   " << vec->E()   << " +- " << 0.5 * sqrt( vec->E() )                           << " -> " << e   << std::endl;
      //if(pt > e) std::cout << "--==!! ERROR !!==--" << std::endl << "pt > e : " << pt << " > " << e << std::endl;
      if(vec - vecs.begin() == 0 || vec - vecs.begin() == 3){
	dynamic_cast<TH2F*>(hists_.find("pt_b" )->second)->Fill( ( pt - pt_ )/pt_, pt_ );
	dynamic_cast<TH2F*>(hists_.find("eta_b")->second)->Fill( eta - eta_, pt_ );
	dynamic_cast<TH2F*>(hists_.find("phi_b")->second)->Fill( deltaPhi(phi, phi_), pt_ );
      }
      else{
	dynamic_cast<TH2F*>(hists_.find("pt_udsc" )->second)->Fill( ( pt - pt_ )/pt_, pt_ );
	dynamic_cast<TH2F*>(hists_.find("eta_udsc")->second)->Fill( eta - eta_, pt_ );
	dynamic_cast<TH2F*>(hists_.find("phi_udsc")->second)->Fill( deltaPhi(phi, phi_), pt_ );
      }
      dynamic_cast<TH2F*>(hists_.find("pt" )->second)->Fill( ( pt - pt_ )/pt_, pt_ );
      dynamic_cast<TH2F*>(hists_.find("eta")->second)->Fill( eta - eta_, pt_ );
      dynamic_cast<TH2F*>(hists_.find("phi")->second)->Fill( deltaPhi(phi, phi_), pt_ );
      //vec->SetPtEtaPhiE( pt, eta, phi, e );
      vec->SetPtEtaPhiM( pt, eta, phi, 0 );
    }
  }
  //std::cout << "------------------------------------------" << std::endl;
}

/// function to do PseudoBTagging with no / one / two BTags
bool
FullHadKinFit::doPseudoBTagging(std::vector< unsigned int > combi)
{
  if( bTags_ == 0 ){
    return true;
  }
  if( bTags_ == 1 ){
    if(combi[0] == 0 || combi[3] == 3 || combi[0] == 3 || combi[3] == 0) return true;
  }
  if( bTags_ == 2 ){
    if( (combi[0] == 0 && combi[3] == 3) || (combi[0] == 3 || combi[3] == 0) ) return true;
  }
  if( bTags_ > 2 ){
    throw cms::Exception("Configuration") << "Wrong number of bTags (" << bTags_ << " bTags not supported)!\n";
    return true;
  }
  return false;
}

/// function to find types of jet-combinations in KinFits (1 right, 2 branches right, but inner-branche particles mixup, 3 inter-branch mixup, 4 missing jet)
int
FullHadKinFit::comboType(std::vector< unsigned int > combi)
{
  int bFit         = combi[B];
  int lightQFit    = combi[LQ];
  int lightQBarFit = combi[LQBar];
  int bBarFit      = combi[BBar];
  int lightPFit    = combi[LP];
  int lightPBarFit = combi[LPBar];

  int bGen         = 0;
  int lightQGen    = 1;
  int lightQBarGen = 2;
  int bBarGen      = 3;
  int lightPGen    = 4;
  int lightPBarGen = 5;

  std::list<int> fit;
  std::list<int> fitBar;
  std::list<int> gen;
  std::list<int> genBar;

  fit.push_back(combi[B]);
  fit.push_back(combi[LQ]);
  fit.push_back(combi[LQBar]);

  fitBar.push_back(combi[BBar]);
  fitBar.push_back(combi[LP]);
  fitBar.push_back(combi[LPBar]);

  gen.push_back(0);
  gen.push_back(1);
  gen.push_back(2);

  genBar.push_back(3);
  genBar.push_back(4);
  genBar.push_back(5);

  if( ( bFit == bGen && bBarFit == bBarGen ) || ( bFit == bBarGen && bBarFit == bGen ) ) {

    if( bFit == bGen &&
	( ( lightQFit == lightQGen    && lightQBarFit == lightQBarGen ) ||
	  ( lightQFit == lightQBarGen && lightQBarFit == lightQGen    ) ) &&
	( ( lightPFit == lightPGen    && lightPBarFit == lightPBarGen ) ||
	  ( lightPFit == lightPBarGen && lightPBarFit == lightPGen    ) ) )
      return 1;
    
    if( bBarFit == bGen &&
	( ( lightPFit == lightQGen    && lightPBarFit == lightQBarGen ) ||
	  ( lightPFit == lightQBarGen && lightPBarFit == lightQGen    ) ) &&
	( ( lightQFit == lightPGen    && lightQBarFit == lightPBarGen ) ||
	  ( lightQFit == lightPBarGen && lightQBarFit == lightPGen    ) ) )
      return 1;
  }

  fit.sort();
  fitBar.sort();
  gen.sort();
  genBar.sort();

  if( ( fit == gen && fitBar == genBar ) || ( fit == genBar && fitBar == gen ) ) return 2;

  fit.splice(fit.end(),fitBar);
  gen.splice(gen.end(),genBar);

  fit.sort();
  fitBar.sort();

  if( fit == gen ) return 3;

  return 4;
}

/// definition of resolutions to be used in kinematic fit
std::vector<FullHadKinFit::KinFitResolutions>
FullHadKinFit::getResolutions(const std::vector< TLorentzVector >* vecs, const std::string whichResolution)
{
  std::vector<KinFitResolutions> result;
  for(std::vector< TLorentzVector >::const_iterator vec = vecs->begin(); vec != vecs->end(); ++vec){
    KinFitResolutions resolu;
    if(resolType_ == "gaussian"){
      double et = vec->Et();
      resolu.Et  = pow(resol_ * 0.3 * sqrt( et ), 2.);
      resolu.Eta = pow(resol_ * 0.5 / sqrt( et ), 2.);
      resolu.Phi = pow(resol_ * 0.5 / sqrt( et ), 2.);
      result.push_back(resolu);
    }
    else{
      int selectedBin=-1;
      reco::LeafCandidate candidate;
      if(whichResolution == "udsc"){
	for(unsigned int i=0; i<binsUdsc_.size(); ++i){
	  StringCutObjectSelector<reco::LeafCandidate> select_(binsUdsc_[i]);
	  candidate = reco::LeafCandidate( 0, reco::LeafCandidate::LorentzVector(vec->Px(), vec->Py(), vec->Pz(), vec->E()) );
	  if(select_(candidate)){
	    selectedBin = i;
	    break;
	  }
	}
	if(selectedBin>=0){
	  resolu.Et  = pow(resol_ * StringObjectFunction<reco::LeafCandidate>(funcEtUdsc_ [selectedBin]).operator()(candidate), 2.);
	  resolu.Eta = pow(resol_ * StringObjectFunction<reco::LeafCandidate>(funcEtaUdsc_[selectedBin]).operator()(candidate), 2.);
	  resolu.Phi = pow(resol_ * StringObjectFunction<reco::LeafCandidate>(funcPhiUdsc_[selectedBin]).operator()(candidate), 2.);
	}
	else{
	  resolu.Et  = 0.00000001;
	  resolu.Eta = 0.00000001;
	  resolu.Phi = 0.00000001;
	}
	
	//std::cout << "LQ: Et: " << vec->Et() << "+/-" << resolu.Et << ", Eta: " << vec->Eta() << "+/-" << resolu.Eta << ", Phi: " << vec->Phi() << "+/-" << resolu.Phi << std::endl;
	
      }
      else if(whichResolution == "b"){
	for(unsigned int i=0; i<binsB_.size(); ++i){
	  StringCutObjectSelector<reco::LeafCandidate> select_(binsB_[i]);
	  candidate = reco::LeafCandidate( 0, reco::LeafCandidate::LorentzVector(vec->Px(), vec->Py(), vec->Pz(), vec->E()) );
	  if(select_(candidate)){
	    selectedBin = i;
	    break;
	  }
	}
	if(selectedBin>=0){
	  resolu.Et  = pow(resol_ * StringObjectFunction<reco::LeafCandidate>(funcEtB_ [selectedBin]).operator()(candidate), 2.);
	  resolu.Eta = pow(resol_ * StringObjectFunction<reco::LeafCandidate>(funcEtaB_[selectedBin]).operator()(candidate), 2.);
	  resolu.Phi = pow(resol_ * StringObjectFunction<reco::LeafCandidate>(funcPhiB_[selectedBin]).operator()(candidate), 2.);
	}
	else{
	  resolu.Et  = 0.00000001;
	  resolu.Eta = 0.00000001;
	  resolu.Phi = 0.00000001;
	}

	//std::cout << "LQ: Et: " << vec->Et() << "+/-" << resolu.Et << ", Eta: " << vec->Eta() << "+/-" << resolu.Eta << ", Phi: " << vec->Phi() << "+/-" << resolu.Phi << std::endl;

      }
      else throw cms::Exception("WrongConfig") << "Only 'udsc' and 'b' are supported types in function FullHadKinFit::getResolutions()!\n";
      result.push_back(resolu);
    }
  }
  return result;
}

/// function to do lower pt cut
bool
FullHadKinFit::cutPt(std::vector< TLorentzVector >& vecs, double ptCut)
{
  for(std::vector< TLorentzVector >::const_iterator vec = vecs.begin(); vec != vecs.end(); ++vec){
    if(vec->Pt() < ptCut) return true;
  }
  return false;
}

/// initiate histograms
void
FullHadKinFit::beginJob()
{
  edm::Service<TFileService> fs;
  if( !fs ) throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

  eventCounter = 0;
  fullHadCounter = 0;
  fitCounter = 0;

  hists_["chi2"         ] = fs->make<TH1F>( "chi2"          , "Chi2;#chi^{2};" , 100 , 0.  , 1.  );
  hists_["prob"         ] = fs->make<TH1F>( "prob"          , "Prob;Prob;"     , 100 , 0.  , 1.  );
  hists_["combo"        ] = fs->make<TH1F>( "combo"         , "Combo;Combo;"   ,   4 , 0.5 , 4.5 );
  hists_["ptGenFit"     ] = fs->make<TH2F>( "ptGenFit"      , "ptGenFit;#frac{pt_{fit}}{pt_{gen}};pt_{gen}", 500 , 0. , 10. , 1000 , 0. , 1000. );
  hists_["ptGenFitOther"] = fs->make<TH2F>( "ptGenFitOther" , "ptGenFit;#frac{pt_{fit}}{pt_{gen}};pt_{gen}", 500 , 0. , 10. , 1000 , 0. , 1000. );
  hists_["ptGenRec"     ] = fs->make<TH2F>( "ptGenRec"      , "ptGenRec;#frac{pt_{rec}}{pt_{gen}};pt_{gen}", 500 , 0. , 10. , 1000 , 0. , 1000. );
  hists_["ptGenRecOther"] = fs->make<TH2F>( "ptGenRecOther" , "ptGenRec;#frac{pt_{rec}}{pt_{gen}};pt_{gen}", 500 , 0. , 10. , 1000 , 0. , 1000. );
  hists_["ptFitRec"     ] = fs->make<TH2F>( "ptFitRec"      , "ptFitRec;#frac{pt_{rec}}{pt_{fit}};pt_{fit}", 500 , 0. , 10. , 1000 , 0. , 1000. );
  hists_["ptFitRecOther"] = fs->make<TH2F>( "ptFitRecOther" , "ptFitRec;#frac{pt_{rec}}{pt_{fit}};pt_{fit}", 500 , 0. , 10. , 1000 , 0. , 1000. );
  hists_["mTop"         ] = fs->make<TH1F>( "mTop"          , "mTop;m_{top};"   , 500 , 0. , 400. );
  hists_["mTopRight"    ] = fs->make<TH1F>( "mTopRight"     , "mTop;m_{top};"   , 500 , 0. , 400. );
  hists_["mTopWrong"    ] = fs->make<TH1F>( "mTopWrong"     , "mTop;m_{top};"   , 500 , 0. , 400. );
  hists_["mW"           ] = fs->make<TH1F>( "mW"            , "mW;m_{W};"       , 500 , 0. , 200. );

  hists_["pullEtaGenFit"     ] = fs->make<TH2F>( "pullEtaGenFit"      , "pullEtaGenFit;#frac{#eta_{fit}}{#eta_{gen}};pt_{gen}", 600 , -2. , 4. , 1000 , 0. , 1000. );
  hists_["pullEtaGenFitOther"] = fs->make<TH2F>( "pullEtaGenFitOther" , "pullEtaGenFit;#frac{#eta_{fit}}{#eta_{gen}};pt_{gen}", 600 , -2. , 4. , 1000 , 0. , 1000. );
  hists_["pullEtaGenRec"     ] = fs->make<TH2F>( "pullEtaGenRec"      , "pullEtaGenRec;#frac{#eta_{rec}}{#eta_{gen}};pt_{gen}", 600 , -2. , 4. , 1000 , 0. , 1000. );
  hists_["pullEtaGenRecOther"] = fs->make<TH2F>( "pullEtaGenRecOther" , "pullEtaGenRec;#frac{#eta_{rec}}{#eta_{gen}};pt_{gen}", 600 , -2. , 4. , 1000 , 0. , 1000. );
  hists_["pullEtaFitRec"     ] = fs->make<TH2F>( "pullEtaFitRec"      , "pullEtaFitRec;#frac{#eta_{rec}}{#eta_{fit}};pt_{fit}", 600 , -2. , 4. , 1000 , 0. , 1000. );
  hists_["pullEtaFitRecOther"] = fs->make<TH2F>( "pullEtaFitRecOther" , "pullEtaFitRec;#frac{#eta_{rec}}{#eta_{fit}};pt_{fit}", 600 , -2. , 4. , 1000 , 0. , 1000. );

  hists_["pullPhiGenFit"     ] = fs->make<TH2F>( "pullPhiGenFit"      , "pullPhiGenFit;#frac{#phi_{fit}}{#phi_{gen}};pt_{gen}", 600 , -2. , 4. , 1000 , 0. , 1000. );
  hists_["pullPhiGenFitOther"] = fs->make<TH2F>( "pullPhiGenFitOther" , "pullPhiGenFit;#frac{#phi_{fit}}{#phi_{gen}};pt_{gen}", 600 , -2. , 4. , 1000 , 0. , 1000. );
  hists_["pullPhiGenRec"     ] = fs->make<TH2F>( "pullPhiGenRec"      , "pullPhiGenRec;#frac{#phi_{rec}}{#phi_{gen}};pt_{gen}", 600 , -2. , 4. , 1000 , 0. , 1000. );
  hists_["pullPhiGenRecOther"] = fs->make<TH2F>( "pullPhiGenRecOther" , "pullPhiGenRec;#frac{#phi_{rec}}{#phi_{gen}};pt_{gen}", 600 , -2. , 4. , 1000 , 0. , 1000. );
  hists_["pullPhiFitRec"     ] = fs->make<TH2F>( "pullPhiFitRec"      , "pullPhiFitRec;#frac{#phi_{rec}}{#phi_{fit}};pt_{fit}", 600 , -2. , 4. , 1000 , 0. , 1000. );
  hists_["pullPhiFitRecOther"] = fs->make<TH2F>( "pullPhiFitRecOther" , "pullPhiFitRec;#frac{#phi_{rec}}{#phi_{fit}};pt_{fit}", 600 , -2. , 4. , 1000 , 0. , 1000. );

  hists_["mWGenPart"  ] = fs->make<TH1F>( "mWGenPart"  , "mW;m_{W};"    , 500 , 0. , 200. );
  hists_["mTopGenPart"] = fs->make<TH1F>( "mTopGenPart", "mTop;m_{top};", 500 , 0. , 400. );

  hists_["pullMTopGenFit"] = fs->make<TH1F>( "pullMTopGenFit", "mTop;m_{top,gen}-m_{top,fit};", 500 , -100. , 100. );
  hists_["pullMTopGenRec"] = fs->make<TH1F>( "pullMTopGenRec", "mTop;m_{top,gen}-m_{top,rec};", 500 , -100. , 100. );
  hists_["pullMTopFitRec"] = fs->make<TH1F>( "pullMTopFitRec", "mTop;m_{top,fit}-m_{top,rec};", 500 , -100. , 100. );
  hists_["pullMWGenFit"]   = fs->make<TH1F>( "pullMWGenFit"  , "mW;m_{W,gen}-m_{W,fit};"      , 500 ,  -50. ,  50. );
  hists_["pullMWGenRec"]   = fs->make<TH1F>( "pullMWGenRec"  , "mW;m_{W,gen}-m_{W,rec};"      , 500 ,  -50. ,  50. );
  hists_["pullMWFitRec"]   = fs->make<TH1F>( "pullMWFitRec"  , "mW;m_{W,fit}-m_{W,rec};"      , 500 ,  -50. ,  50. );

  hists_["pullPtWGenFit"  ] = fs->make<TH1F>( "pullPtWGenFit"  , "pt-pull W  ;#frac{pt_{gen}-pt_{fit}}{pt_{gen}};", 400 , -2. , 2. );
  hists_["pullPtTopGenFit"] = fs->make<TH1F>( "pullPtTopGenFit", "pt-pull Top;#frac{pt_{gen}-pt_{fit}}{pt_{gen}};", 400 , -2. , 2. );
  hists_["pullPtWGenRec"  ] = fs->make<TH1F>( "pullPtWGenRec"  , "pt-pull W  ;#frac{pt_{gen}-pt_{rec}}{pt_{gen}};", 400 , -2. , 2. );
  hists_["pullPtTopGenRec"] = fs->make<TH1F>( "pullPtTopGenRec", "pt-pull Top;#frac{pt_{gen}-pt_{rec}}{pt_{gen}};", 400 , -2. , 2. );
  hists_["pullPtWFitRec"  ] = fs->make<TH1F>( "pullPtWFitRec"  , "pt-pull W  ;#frac{pt_{Fit}-pt_{rec}}{pt_{fit}};", 400 , -2. , 2. );
  hists_["pullPtTopFitRec"] = fs->make<TH1F>( "pullPtTopFitRec", "pt-pull Top;#frac{pt_{Fit}-pt_{rec}}{pt_{fit}};", 400 , -2. , 2. );

  hists_["pullEtaWGenFit"  ] = fs->make<TH1F>( "pullEtaWGenFit"  , "#eta-pull W  ;#eta_{gen}-#eta_{fit};", 400 , -2. , 2. );
  hists_["pullEtaTopGenFit"] = fs->make<TH1F>( "pullEtaTopGenFit", "#eta-pull Top;#eta_{gen}-#eta_{fit};", 400 , -2. , 2. );
  hists_["pullEtaWGenRec"  ] = fs->make<TH1F>( "pullEtaWGenRec"  , "#eta-pull W  ;#eta_{gen}-#eta_{rec};", 400 , -2. , 2. );
  hists_["pullEtaTopGenRec"] = fs->make<TH1F>( "pullEtaTopGenRec", "#eta-pull Top;#eta_{gen}-#eta_{rec};", 400 , -2. , 2. );
  hists_["pullEtaWFitRec"  ] = fs->make<TH1F>( "pullEtaWFitRec"  , "#eta-pull W  ;#eta_{Fit}-#eta_{rec};", 400 , -2. , 2. );
  hists_["pullEtaTopFitRec"] = fs->make<TH1F>( "pullEtaTopFitRec", "#eta-pull Top;#eta_{Fit}-#eta_{rec};", 400 , -2. , 2. );

  hists_["pullPhiWGenFit"  ] = fs->make<TH1F>( "pullPhiWGenFit"  , "#phi-pull W  ;#phi_{gen}-#phi_{fit};", 400 , -2. , 2. );
  hists_["pullPhiTopGenFit"] = fs->make<TH1F>( "pullPhiTopGenFit", "#phi-pull Top;#phi_{gen}-#phi_{fit};", 400 , -2. , 2. );
  hists_["pullPhiWGenRec"  ] = fs->make<TH1F>( "pullPhiWGenRec"  , "#phi-pull W  ;#phi_{gen}-#phi_{rec};", 400 , -2. , 2. );
  hists_["pullPhiTopGenRec"] = fs->make<TH1F>( "pullPhiTopGenRec", "#phi-pull Top;#phi_{gen}-#phi_{rec};", 400 , -2. , 2. );
  hists_["pullPhiWFitRec"  ] = fs->make<TH1F>( "pullPhiWFitRec"  , "#phi-pull W  ;#phi_{Fit}-#phi_{rec};", 400 , -2. , 2. );
  hists_["pullPhiTopFitRec"] = fs->make<TH1F>( "pullPhiTopFitRec", "#phi-pull Top;#phi_{Fit}-#phi_{rec};", 400 , -2. , 2. );

  hists_["prob_0tail"] = fs->make<TH1F>( "prob_0tail", "0Tail;Prob;"   , 100 , 0. , 1. );
  hists_["prob_1tail"] = fs->make<TH1F>( "prob_1tail", "1Tail;Prob;"   , 100 , 0. , 1. );
  hists_["prob_2tail"] = fs->make<TH1F>( "prob_2tail", "2Tail;Prob;"   , 100 , 0. , 1. );
  hists_["prob_3tail"] = fs->make<TH1F>( "prob_3tail", "3Tail;Prob;"   , 100 , 0. , 1. );
  hists_["prob_4tail"] = fs->make<TH1F>( "prob_4tail", "4Tail;Prob;"   , 100 , 0. , 1. );
  hists_["prob_gaus" ] = fs->make<TH1F>( "prob_gaus" , "Gaussian;Prob;", 100 , 0. , 1. );

  //TFile * file = new TFile("./analyzePatObjectResolution_Jets_new.root");
  TFile * file = new TFile("./JetResolution.root");
  
  if(smearType_ == "calo"){
    hists_["bResolution_pt"    ] = (TH2F*)((file->Get("analyzePatJetResolutions/pt_b"))->Clone("bResolution_pt"));
    hists_["bResolution_eta"   ] = (TH2F*)((file->Get("analyzePatJetResolutions/eta_b"))->Clone("bResolution_eta"));
    hists_["bResolution_phi"   ] = (TH2F*)((file->Get("analyzePatJetResolutions/phi_b"))->Clone("bResolution_phi"));
    hists_["udscResolution_pt" ] = (TH2F*)((file->Get("analyzePatJetResolutions/pt_udsc"))->Clone("udscResolution_pt"));
    hists_["udscResolution_eta"] = (TH2F*)((file->Get("analyzePatJetResolutions/eta_udsc"))->Clone("udscResolution_eta"));
    hists_["udscResolution_phi"] = (TH2F*)((file->Get("analyzePatJetResolutions/phi_udsc"))->Clone("udscResolution_phi"));
  }
  else if(smearType_ == "pf"){
    hists_["bResolution_pt"    ] = (TH2F*)((file->Get("analyzePatJetResolutionsPF/pt_b"))->Clone("bResolution_pt"));
    hists_["bResolution_eta"   ] = (TH2F*)((file->Get("analyzePatJetResolutionsPF/eta_b"))->Clone("bResolution_eta"));
    hists_["bResolution_phi"   ] = (TH2F*)((file->Get("analyzePatJetResolutionsPF/phi_b"))->Clone("bResolution_phi"));
    hists_["udscResolution_pt" ] = (TH2F*)((file->Get("analyzePatJetResolutionsPF/pt_udsc"))->Clone("udscResolution_pt"));
    hists_["udscResolution_eta"] = (TH2F*)((file->Get("analyzePatJetResolutionsPF/eta_udsc"))->Clone("udscResolution_eta"));
    hists_["udscResolution_phi"] = (TH2F*)((file->Get("analyzePatJetResolutionsPF/phi_udsc"))->Clone("udscResolution_phi"));
  }
  else if(smearType_ == "gaussian"){
  }
  else{
    throw cms::Exception("WrongConfig") << "Parameter 'smearType' not defined properly, needs to be *calo*, *pf* or *gaussian*!\n";
  }

  delete file;

  // for validation of resolutions

  hists_["pt"]  = fs->make<TH2F>( "pt"  , "pt"  , 2000 , -1. , 10.  , 1000, 0., 1000.);
  hists_["eta"] = fs->make<TH2F>( "eta" , "eta" ,  200 , -0.5 , 0.5 , 1000, 0., 1000.);
  hists_["phi"] = fs->make<TH2F>( "phi" , "phi" ,  200 , -0.5 , 0.5 , 1000, 0., 1000.);

  hists_["pt_b"]  = fs->make<TH2F>( "pt_b"  , "pt_b"  , 2000 , -1. , 10.  , 1000, 0., 1000.);
  hists_["eta_b"] = fs->make<TH2F>( "eta_b" , "eta_b" ,  200 , -0.5 , 0.5 , 1000, 0., 1000.);
  hists_["phi_b"] = fs->make<TH2F>( "phi_b" , "phi_b" ,  200 , -0.5 , 0.5 , 1000, 0., 1000.);

  hists_["pt_udsc"]  = fs->make<TH2F>( "pt_udsc"  , "pt_udsc"  , 2000 , -1. , 10.  , 1000, 0., 1000.);
  hists_["eta_udsc"] = fs->make<TH2F>( "eta_udsc" , "eta_udsc" ,  200 , -0.5 , 0.5 , 1000, 0., 1000.);
  hists_["phi_udsc"] = fs->make<TH2F>( "phi_udsc" , "phi_udsc" ,  200 , -0.5 , 0.5 , 1000, 0., 1000.);
}

/// do the fitting
std::list< FullHadKinFit::KinFitResult >
FullHadKinFit::fit(std::vector< TLorentzVector >& vecs)
{
  // store results of the kinematic fit
  std::list< KinFitResult > fitResults;

  std::vector< TLorentzVector > genVecs = vecs;

  if(smear_) smear(vecs);

  std::vector< TLorentzVector > smearVecs = vecs;

  if(recoPtCut_ != 0. && cutPt(vecs,recoPtCut_)){
    std::vector< unsigned int > combi(6);
    combi[B]     = 0;
    combi[LQ]    = 1;
    combi[LQBar] = 2;
    combi[BBar]  = 3;
    combi[LP]    = 4;
    combi[LPBar] = 5;

    KinFitResult result;
    result.GenVecs = genVecs;
    result.IniVecs = smearVecs;
    result.FitVecs = smearVecs;
    result.Chi2 = -1.;
    result.NDF  = -1.;
    result.Prob = -1.;
    result.JetCombi = combi;
    fitResults.push_back( result );
    return fitResults;
  }

  TMatrixD m0(3,3);
  TMatrixD m1(3,3);
  TMatrixD m2(3,3);
  TMatrixD m3(3,3);
  TMatrixD m4(3,3);
  TMatrixD m5(3,3);
  m0.Zero();
  m1.Zero();
  m2.Zero();
  m3.Zero();
  m4.Zero();
  m5.Zero();

  std::vector<KinFitResolutions> resolutionsUdsc = getResolutions(&vecs, "udsc");

  //In this example the covariant matrix depends on the transverse energy and eta of the jets
  m0(0,0) = resolutionsUdsc[0].Et; // et
  m0(1,1) = resolutionsUdsc[0].Eta; // eta
  m0(2,2) = resolutionsUdsc[0].Phi; // phi

  m1(0,0) = resolutionsUdsc[1].Et; // et
  m1(1,1) = resolutionsUdsc[1].Eta; // eta
  m1(2,2) = resolutionsUdsc[1].Phi; // phi

  m2(0,0) = resolutionsUdsc[2].Et; // et
  m2(1,1) = resolutionsUdsc[2].Eta; // eta
  m2(2,2) = resolutionsUdsc[2].Phi; // phi

  m3(0,0) = resolutionsUdsc[3].Et; // et
  m3(1,1) = resolutionsUdsc[3].Eta; // eta
  m3(2,2) = resolutionsUdsc[3].Phi; // phi

  m4(0,0) = resolutionsUdsc[4].Et; // et
  m4(1,1) = resolutionsUdsc[4].Eta; // eta
  m4(2,2) = resolutionsUdsc[4].Phi; // phi

  m5(0,0) = resolutionsUdsc[5].Et; // et
  m5(1,1) = resolutionsUdsc[5].Eta; // eta
  m5(2,2) = resolutionsUdsc[5].Phi; // phi

  std::vector< TMatrixD > matsUdsc;
  matsUdsc.push_back( m0 );
  matsUdsc.push_back( m1 );
  matsUdsc.push_back( m2 );
  matsUdsc.push_back( m3 );
  matsUdsc.push_back( m4 );
  matsUdsc.push_back( m5 );

  m0.Zero();
  m1.Zero();
  m2.Zero();
  m3.Zero();
  m4.Zero();
  m5.Zero();

  std::vector<KinFitResolutions> resolutionsB = getResolutions(&vecs, "b");

  //In this example the covariant matrix depends on the transverse energy and eta of the jets
  m0(0,0) = resolutionsB[0].Et; // et
  m0(1,1) = resolutionsB[0].Eta; // eta
  m0(2,2) = resolutionsB[0].Phi; // phi

  m1(0,0) = resolutionsB[1].Et; // et
  m1(1,1) = resolutionsB[1].Eta; // eta
  m1(2,2) = resolutionsB[1].Phi; // phi

  m2(0,0) = resolutionsB[2].Et; // et
  m2(1,1) = resolutionsB[2].Eta; // eta
  m2(2,2) = resolutionsB[2].Phi; // phi

  m3(0,0) = resolutionsB[3].Et; // et
  m3(1,1) = resolutionsB[3].Eta; // eta
  m3(2,2) = resolutionsB[3].Phi; // phi

  m4(0,0) = resolutionsB[4].Et; // et
  m4(1,1) = resolutionsB[4].Eta; // eta
  m4(2,2) = resolutionsB[4].Phi; // phi

  m5(0,0) = resolutionsB[5].Et; // et
  m5(1,1) = resolutionsB[5].Eta; // eta
  m5(2,2) = resolutionsB[5].Phi; // phi

  std::vector< TMatrixD > matsB;
  matsB.push_back( m0 );
  matsB.push_back( m1 );
  matsB.push_back( m2 );
  matsB.push_back( m3 );
  matsB.push_back( m4 );
  matsB.push_back( m5 );

  std::vector< unsigned int > combi(6);
  combi[B]     = 0;
  combi[LQ]    = 1;
  combi[LQBar] = 2;
  combi[BBar]  = 3;
  combi[LP]    = 4;
  combi[LPBar] = 5;

  for(unsigned int cnt=0; cnt<TMath::Factorial(combi.size()); ++cnt){
    if( ( combi[LQ] < combi[LQBar]  ||
	  combi[LP] < combi[LPBar]  ||
	  combi[B]  < combi[BBar] ) &&
	doPseudoBTagging(combi) ){

      TFitParticleEtEtaPhi *cand_b     = new TFitParticleEtEtaPhi( "cand_b"    , "cand_b"    , &vecs[combi[B]]    , &matsB[combi[B]]     );
      TFitParticleEtEtaPhi *cand_bBar  = new TFitParticleEtEtaPhi( "cand_bBar" , "cand_bBar" , &vecs[combi[BBar]] , &matsB[combi[BBar]]  );
      TFitParticleEtEtaPhi *cand_lQ    = new TFitParticleEtEtaPhi( "cand_lQ"   , "cand_lQ"   , &vecs[combi[LQ]]   , &matsUdsc[combi[LQ]]    );
      TFitParticleEtEtaPhi *cand_lQBar = new TFitParticleEtEtaPhi( "cand_lQBar", "cand_lQBar", &vecs[combi[LQBar]], &matsUdsc[combi[LQBar]] );
      TFitParticleEtEtaPhi *cand_lP    = new TFitParticleEtEtaPhi( "cand_lP"   , "cand_lP"   , &vecs[combi[LP]]   , &matsUdsc[combi[LP]]    );
      TFitParticleEtEtaPhi *cand_lPBar = new TFitParticleEtEtaPhi( "cand_lPBar", "cand_lPBar", &vecs[combi[LPBar]], &matsUdsc[combi[LPBar]] );

      TFitConstraintM *mCons1 = new TFitConstraintM( "WMassConstraint1", "WMass-Constraint1", 0, 0 , 80.4 );
      //TFitConstraintMGaus *mCons1 = new TFitConstraintMGaus( "WMassConstraint1", "WMass-Constraint1", 0, 0 , 80.4 , 2.141 );
      mCons1->addParticles1( cand_lQ, cand_lQBar );
      TFitConstraintM *mCons2 = new TFitConstraintM( "WMassConstraint2", "WMass-Constraint2", 0, 0 , 80.4 );
      //TFitConstraintMGaus *mCons2 = new TFitConstraintMGaus( "WMassConstraint2", "WMass-Constraint2", 0, 0 , 80.4 , 2.141 );
      mCons2->addParticles1( cand_lP, cand_lPBar );
      TFitConstraintM *mCons3 = new TFitConstraintM( "EqualTopMasses", "EqualTopMasses" , 0 , 0 , 0 );
      mCons3->addParticles1( cand_b   , cand_lQ, cand_lQBar );
      mCons3->addParticles2( cand_bBar, cand_lP, cand_lPBar );

      //TFitConstraintM *mCons4 = new TFitConstraintM( "TopMassConstraint1", "TopMass-Constraint1", 0, 0 , 173.0 );
      //mCons4->addParticles1( cand_b   , cand_lQ, cand_lQBar );
      //TFitConstraintM *mCons5 = new TFitConstraintM( "TopMassConstraint2", "TopMass-Constraint2", 0, 0 , 173.0 );
      //mCons5->addParticles1( cand_bBar, cand_lP, cand_lPBar );


      //Definition of the fitter
      TKinFitter* fitter = new TKinFitter("fitter", "fitter");
      //Add measured particles(jets)
      fitter->addMeasParticle( cand_b     );
      fitter->addMeasParticle( cand_bBar  );
      fitter->addMeasParticle( cand_lQ    );
      fitter->addMeasParticle( cand_lQBar );
      fitter->addMeasParticle( cand_lP    );
      fitter->addMeasParticle( cand_lPBar );
      //Add constraints
      fitter->addConstraint( mCons1 );
      fitter->addConstraint( mCons2 );
      fitter->addConstraint( mCons3 );
      //fitter->addConstraint( mCons4 );
      //fitter->addConstraint( mCons5 );

      //Set convergence criteria
      fitter->setMaxNbIter( 300 );
      fitter->setMaxDeltaS( 1e-2 );
      fitter->setMaxF( 1e-1 );
      fitter->setVerbosity(1);

      //Perform the fit
      fitter->fit();

      if(fitter->getStatus() == 0){
	    
	std::vector< TLorentzVector > helperVec;
	helperVec.resize(6);
	helperVec[B]     = *cand_b    ->getCurr4Vec();
	helperVec[LQ]    = *cand_lQ   ->getCurr4Vec();
	helperVec[LQBar] = *cand_lQBar->getCurr4Vec();
	helperVec[BBar]  = *cand_bBar ->getCurr4Vec();
	helperVec[LP]    = *cand_lP   ->getCurr4Vec();
	helperVec[LPBar] = *cand_lPBar->getCurr4Vec();

	KinFitResult result;
	result.GenVecs = genVecs;
	result.IniVecs = smearVecs;
	result.FitVecs = helperVec;
	result.Chi2 = fitter->getS();
	result.NDF = fitter->getNDF();
	result.Prob = TMath::Prob(result.Chi2, result.NDF);
	result.JetCombi = combi;
	fitResults.push_back( result );
      }

      delete cand_b;
      delete cand_bBar;
      delete cand_lQ;
      delete cand_lQBar;
      delete cand_lP;
      delete cand_lPBar;
      delete mCons1;
      delete mCons2;
      delete mCons3;
      //delete mCons4;
      //delete mCons5;
      delete fitter;
    }
    if(onlyGenMatch_) break;
    std::next_permutation( combi.begin(), combi.end() );
  }
  return fitResults;
}

// analyze kinFit and fill histograms
void
FullHadKinFit::analyze(const edm::Event& event, const edm::EventSetup&)
{
  ++eventCounter;

  edm::Handle<edm::View< reco::LeafCandidate > > src;
  event.getByLabel(src_, src);

  genMatcher(*src);

  if(!genParticles.size()) { genParticles.clear(); return; }
  if(!genParticles[B] || !genParticles[BBar] || !genParticles[LQ] || !genParticles[LQBar] || !genParticles[LP] || !genParticles[LPBar]) { genParticles.clear(); return; }
  int nmbOfRad = foundRadiation(*src, genParticles);
  if( !(!exchangePart_ &&  noRadiations_ && nmbOfRad == 0) &&
      !(!exchangePart_ && !noRadiations_                 ) &&
      !( exchangePart_ && !noRadiations_ && nmbOfRad  > 0) ) { genParticles.clear(); return; }

  ++fullHadCounter;

  std::vector< TLorentzVector > vecs;
  vecs.resize(6);
  vecs[B]     = TLorentzVector( genParticles[B]    ->px(), genParticles[B]    ->py(), genParticles[B]    ->pz(), genParticles[B]    ->energy() );
  vecs[LQ]    = TLorentzVector( genParticles[LQ]   ->px(), genParticles[LQ]   ->py(), genParticles[LQ]   ->pz(), genParticles[LQ]   ->energy() );
  vecs[LQBar] = TLorentzVector( genParticles[LQBar]->px(), genParticles[LQBar]->py(), genParticles[LQBar]->pz(), genParticles[LQBar]->energy() );
  vecs[BBar]  = TLorentzVector( genParticles[BBar] ->px(), genParticles[BBar] ->py(), genParticles[BBar] ->pz(), genParticles[BBar] ->energy() );
  vecs[LP]    = TLorentzVector( genParticles[LP]   ->px(), genParticles[LP]   ->py(), genParticles[LP]   ->pz(), genParticles[LP]   ->energy() );
  vecs[LPBar] = TLorentzVector( genParticles[LPBar]->px(), genParticles[LPBar]->py(), genParticles[LPBar]->pz(), genParticles[LPBar]->energy() );

  if(genPtCut_ != 0. && cutPt(vecs,genPtCut_)) { genParticles.clear(); return; }

  std::list< KinFitResult > fitResults;
  if(!fitAsFitInpt_) fitResults = fit(vecs);
  else{
    std::list< KinFitResult > preFitResults = fit(vecs);
    if(preFitResults.size()){
      preFitResults.sort();
      fitResults = fit(preFitResults.front().FitVecs);
    }
    else{
      std::cout << "First fit did not converge, not doing a second one!" << std::endl;
      genParticles.clear();
      return;
    }
  }

  fitResults.remove_if(negativeChi2());

  if(!fitResults.size()) { genParticles.clear(); return; }
  ++fitCounter;

  fitResults.sort();

  for(std::list< KinFitResult >::const_iterator fitResult = fitResults.begin(); fitResult != fitResults.end(); ++fitResult){

    TLorentzVector wGen1 = fitResult->GenVecs[LQ] + fitResult->GenVecs[LQBar];
    TLorentzVector wGen2 = fitResult->GenVecs[LP] + fitResult->GenVecs[LPBar];
    TLorentzVector wRec1 = fitResult->IniVecs[LQ] + fitResult->IniVecs[LQBar];
    TLorentzVector wRec2 = fitResult->IniVecs[LP] + fitResult->IniVecs[LPBar];
    TLorentzVector wFit1 = fitResult->FitVecs[LQ] + fitResult->FitVecs[LQBar];
    TLorentzVector wFit2 = fitResult->FitVecs[LP] + fitResult->FitVecs[LPBar];
    
    TLorentzVector topGen1 = fitResult->GenVecs[LQ] + fitResult->GenVecs[LQBar] + fitResult->GenVecs[B]   ;
    TLorentzVector topGen2 = fitResult->GenVecs[LP] + fitResult->GenVecs[LPBar] + fitResult->GenVecs[BBar];
    TLorentzVector topRec1 = fitResult->IniVecs[LQ] + fitResult->IniVecs[LQBar] + fitResult->IniVecs[B]   ;
    TLorentzVector topRec2 = fitResult->IniVecs[LP] + fitResult->IniVecs[LPBar] + fitResult->IniVecs[BBar];
    TLorentzVector topFit1 = fitResult->FitVecs[LQ] + fitResult->FitVecs[LQBar] + fitResult->FitVecs[B]   ;
    TLorentzVector topFit2 = fitResult->FitVecs[LP] + fitResult->FitVecs[LPBar] + fitResult->FitVecs[BBar];

    double wMassGen1  = (wGen1).M(); double wPtGen1    = (wGen1).Pt(); double wEtaGen1   = (wGen1).Eta(); double wPhiGen1   = (wGen1).Phi();
    double wMassGen2  = (wGen2).M(); double wPtGen2    = (wGen2).Pt(); double wEtaGen2   = (wGen2).Eta(); double wPhiGen2   = (wGen2).Phi();
    double wMassRec1  = (wRec1).M(); double wPtRec1    = (wRec1).Pt(); double wEtaRec1   = (wRec1).Eta(); double wPhiRec1   = (wRec1).Phi();
    double wMassRec2  = (wRec2).M(); double wPtRec2    = (wRec2).Pt(); double wEtaRec2   = (wRec2).Eta(); double wPhiRec2   = (wRec2).Phi();
    double wMass1     = (wFit1).M(); double wPt1       = (wFit1).Pt(); double wEta1      = (wFit1).Eta(); double wPhi1      = (wFit1).Phi();
    double wMass2     = (wFit2).M(); double wPt2       = (wFit2).Pt(); double wEta2      = (wFit2).Eta(); double wPhi2      = (wFit2).Phi();
   
    double topMassGen1  = (topGen1).M(); double topPtGen1    = (topGen1).Pt(); double topEtaGen1   = (topGen1).Eta(); double topPhiGen1   = (topGen1).Phi();
    double topMassGen2  = (topGen2).M(); double topPtGen2    = (topGen2).Pt(); double topEtaGen2   = (topGen2).Eta(); double topPhiGen2   = (topGen2).Phi();
    double topMassRec1  = (topRec1).M(); double topPtRec1    = (topRec1).Pt(); double topEtaRec1   = (topRec1).Eta(); double topPhiRec1   = (topRec1).Phi();
    double topMassRec2  = (topRec2).M(); double topPtRec2    = (topRec2).Pt(); double topEtaRec2   = (topRec2).Eta(); double topPhiRec2   = (topRec2).Phi();
    double topMass1     = (topFit1).M(); double topPt1       = (topFit1).Pt(); double topEta1      = (topFit1).Eta(); double topPhi1      = (topFit1).Phi();
    double topMass2     = (topFit2).M(); double topPt2       = (topFit2).Pt(); double topEta2      = (topFit2).Eta(); double topPhi2      = (topFit2).Phi();

    hists_.find("mWGenPart")  ->second->Fill( wMassGen1 );
    hists_.find("mWGenPart")  ->second->Fill( wMassGen2 );
    hists_.find("mTopGenPart")->second->Fill( topMassGen1 );
    hists_.find("mTopGenPart")->second->Fill( topMassGen2 );

    hists_.find("pullMTopGenFit")->second->Fill( topMassGen1 - topMass1    ); hists_.find("pullMTopGenFit")->second->Fill( topMassGen2 - topMass2    );
    hists_.find("pullMTopGenRec")->second->Fill( topMassGen1 - topMassRec1 ); hists_.find("pullMTopGenRec")->second->Fill( topMassGen2 - topMassRec2 );
    hists_.find("pullMTopFitRec")->second->Fill( topMass1    - topMassRec1 ); hists_.find("pullMTopFitRec")->second->Fill( topMass2    - topMassRec2 );
    hists_.find("pullMWGenFit"  )->second->Fill( wMassGen1   - wMass1      ); hists_.find("pullMWGenFit"  )->second->Fill( wMassGen2   - wMass2      );
    hists_.find("pullMWGenRec"  )->second->Fill( wMassGen1   - wMassRec1   ); hists_.find("pullMWGenRec"  )->second->Fill( wMassGen2   - wMassRec2   );
    hists_.find("pullMWFitRec"  )->second->Fill( wMass1      - wMassRec1   ); hists_.find("pullMWFitRec"  )->second->Fill( wMass2      - wMassRec2   );
				   
    hists_.find("pullPtWGenFit"  )->second->Fill( (wPtGen1   - wPt1     )/wPtGen1   ); hists_.find("pullPtWGenFit"  )->second->Fill( (wPtGen2   - wPt2     )/wPtGen2   );
    hists_.find("pullPtTopGenFit")->second->Fill( (topPtGen1 - topPt1   )/topPtGen1 ); hists_.find("pullPtTopGenFit")->second->Fill( (topPtGen2 - topPt2   )/topPtGen2 );
    hists_.find("pullPtWGenRec"  )->second->Fill( (wPtGen1   - wPtRec1  )/wPtGen1   ); hists_.find("pullPtWGenRec"  )->second->Fill( (wPtGen2   - wPtRec2  )/wPtGen2   );
    hists_.find("pullPtTopGenRec")->second->Fill( (topPtGen1 - topPtRec1)/topPtGen1 ); hists_.find("pullPtTopGenRec")->second->Fill( (topPtGen2 - topPtRec2)/topPtGen2 );
    hists_.find("pullPtWFitRec"  )->second->Fill( (wPt1      - wPtRec1  )/wPt1      ); hists_.find("pullPtWFitRec"  )->second->Fill( (wPt2      - wPtRec2  )/wPt2      );
    hists_.find("pullPtTopFitRec")->second->Fill( (topPt1    - topPtRec1)/topPt1    ); hists_.find("pullPtTopFitRec")->second->Fill( (topPt2    - topPtRec2)/topPt2    );
				   
    hists_.find("pullEtaWGenFit"  )->second->Fill( wEtaGen1   - wEta1      ); hists_.find("pullEtaWGenFit"  )->second->Fill( wEtaGen2   - wEta2      );
    hists_.find("pullEtaTopGenFit")->second->Fill( topEtaGen1 - topEta1    ); hists_.find("pullEtaTopGenFit")->second->Fill( topEtaGen2 - topEta2    );
    hists_.find("pullEtaWGenRec"  )->second->Fill( wEtaGen1   - wEtaRec1   ); hists_.find("pullEtaWGenRec"  )->second->Fill( wEtaGen2   - wEtaRec2   );
    hists_.find("pullEtaTopGenRec")->second->Fill( topEtaGen1 - topEtaRec1 ); hists_.find("pullEtaTopGenRec")->second->Fill( topEtaGen2 - topEtaRec2 );
    hists_.find("pullEtaWFitRec"  )->second->Fill( wEta1      - wEtaRec1   ); hists_.find("pullEtaWFitRec"  )->second->Fill( wEta2      - wEtaRec2   );
    hists_.find("pullEtaTopFitRec")->second->Fill( topEta1    - topEtaRec1 ); hists_.find("pullEtaTopFitRec")->second->Fill( topEta2    - topEtaRec2 );
				   
    hists_.find("pullPhiWGenFit"  )->second->Fill( wPhiGen1   - wPhi1      ); hists_.find("pullPhiWGenFit"  )->second->Fill( wPhiGen2   - wPhi2      );
    hists_.find("pullPhiTopGenFit")->second->Fill( topPhiGen1 - topPhi1    ); hists_.find("pullPhiTopGenFit")->second->Fill( topPhiGen2 - topPhi2    );
    hists_.find("pullPhiWGenRec"  )->second->Fill( wPhiGen1   - wPhiRec1   ); hists_.find("pullPhiWGenRec"  )->second->Fill( wPhiGen2   - wPhiRec2   );
    hists_.find("pullPhiTopGenRec")->second->Fill( topPhiGen1 - topPhiRec1 ); hists_.find("pullPhiTopGenRec")->second->Fill( topPhiGen2 - topPhiRec2 );
    hists_.find("pullPhiWFitRec"  )->second->Fill( wPhi1      - wPhiRec1   ); hists_.find("pullPhiWFitRec"  )->second->Fill( wPhi2      - wPhiRec2   );
    hists_.find("pullPhiTopFitRec")->second->Fill( topPhi1    - topPhiRec1 ); hists_.find("pullPhiTopFitRec")->second->Fill( topPhi2    - topPhiRec2 );
 
    hists_.find("chi2" )->second->Fill(fitResult->Chi2);
    hists_.find("prob" )->second->Fill(fitResult->Prob);
    int combinationType = comboType(fitResult->JetCombi);
    hists_.find("combo")->second->Fill(combinationType);
    hists_.find("mTop" )->second->Fill( topMass1 );
    hists_.find("mW")->second->Fill( wMass1 );
    hists_.find("mW")->second->Fill( wMass2 );
    if(combinationType == 1 || combinationType == 2) hists_.find("mTopRight")->second->Fill( topMass1 );
    else                                             hists_.find("mTopWrong")->second->Fill( topMass1 );
	    
    unsigned int tailCounter = 0;
    if(fabs(wMassGen1 - 80.4) > 3.5)     ++tailCounter;
    if(fabs(wMassGen2 - 80.4) > 3.5)     ++tailCounter;
    if(fabs(topMassGen1 - 172.0) > 12.0) ++tailCounter;
    if(fabs(topMassGen2 - 172.0) > 12.0) ++tailCounter;

    if     (tailCounter == 0) hists_.find("prob_0tail")->second->Fill(fitResult->Prob);
    else if(tailCounter == 1) hists_.find("prob_1tail")->second->Fill(fitResult->Prob);
    else if(tailCounter == 2) hists_.find("prob_2tail")->second->Fill(fitResult->Prob);
    else if(tailCounter == 3) hists_.find("prob_3tail")->second->Fill(fitResult->Prob);
    else if(tailCounter == 4) hists_.find("prob_4tail")->second->Fill(fitResult->Prob);
    else throw cms::Exception("CalculationError") << "Counted higher than possible! 'tailCounter' = " << tailCounter << " , but has to be less than 5!\n";

    bool nonGaussianSmearing = false;

    std::vector< TLorentzVector >::const_iterator vec_gen = fitResult->GenVecs.begin();
    std::vector< TLorentzVector >::const_iterator vec_rec = fitResult->IniVecs.begin();
    for(std::vector< TLorentzVector >::const_iterator vec = fitResult->FitVecs.begin() ; vec != fitResult->FitVecs.end() ; ++vec, ++vec_gen, ++vec_rec){
      if(vec_gen == fitResult->GenVecs.end() || vec_rec == fitResult->IniVecs.end()) break;
      
      if(tailCounter > 0 || (-vec_gen->Pt()/266 + 2) < ((vec_rec->Pt() - vec_gen->Pt())/vec_gen->Pt())) nonGaussianSmearing = true;

      if(smearOnly_ == -1 || smearOnly_ == vec - fitResult->FitVecs.begin()){
	hists_.find("ptGenRec")->second->Fill(vec_rec->Pt()/vec_gen->Pt(), vec_gen->Pt());
	hists_.find("ptGenFit")->second->Fill(vec->Pt()    /vec_gen->Pt(), vec_gen->Pt());
	hists_.find("ptFitRec")->second->Fill(vec_rec->Pt()/vec    ->Pt(), vec    ->Pt());
	
	hists_.find("pullEtaGenRec")->second->Fill(vec_rec->Eta()/vec_gen->Eta(), vec_gen->Pt());
	hists_.find("pullEtaGenFit")->second->Fill(vec->Eta()    /vec_gen->Eta(), vec_gen->Pt());
	hists_.find("pullEtaFitRec")->second->Fill(vec_rec->Eta()/vec    ->Eta(), vec    ->Pt());
	
	hists_.find("pullPhiGenRec")->second->Fill(vec_rec->Phi()/vec_gen->Phi(), vec_gen->Pt());
	hists_.find("pullPhiGenFit")->second->Fill(vec->Phi()    /vec_gen->Phi(), vec_gen->Pt());
	hists_.find("pullPhiFitRec")->second->Fill(vec_rec->Phi()/vec    ->Phi(), vec    ->Pt());
      }
      else{
	hists_.find("ptGenRecOther")->second->Fill(vec_rec->Pt()/vec_gen->Pt(), vec_gen->Pt());
	hists_.find("ptGenFitOther")->second->Fill(vec->Pt()    /vec_gen->Pt(), vec_gen->Pt());
	hists_.find("ptFitRecOther")->second->Fill(vec_rec->Pt()/vec    ->Pt(), vec    ->Pt());
	
	hists_.find("pullEtaGenRecOther")->second->Fill(vec_rec->Eta()/vec_gen->Eta(), vec_gen->Pt());
	hists_.find("pullEtaGenFitOther")->second->Fill(vec->Eta()    /vec_gen->Eta(), vec_gen->Pt());
	hists_.find("pullEtaFitRecOther")->second->Fill(vec_rec->Eta()/vec    ->Eta(), vec    ->Pt());
	
	hists_.find("pullPhiGenRecOther")->second->Fill(vec_rec->Phi()/vec_gen->Phi(), vec_gen->Pt());
	hists_.find("pullPhiGenFitOther")->second->Fill(vec->Phi()    /vec_gen->Phi(), vec_gen->Pt());
	hists_.find("pullPhiFitRecOther")->second->Fill(vec_rec->Phi()/vec    ->Phi(), vec    ->Pt());
      }
    }
    if(!nonGaussianSmearing) hists_.find("prob_gaus")->second->Fill(fitResult->Prob);
    if(!saveAllCombo_) break;
  }
  genParticles.clear();
}

void
FullHadKinFit::endJob()
{
  std::cout << " Total Events   : " << eventCounter   << std::endl;
  std::cout << " Fully hadronic : " << fullHadCounter << std::endl;
  std::cout << " Fitted events  : " << fitCounter     << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( FullHadKinFit );
