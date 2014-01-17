#include "TopAnalysis/TopAnalyzer/interface/SingleObjectJetTrigger.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TMath.h"

//#include "DataFormats/PatCandidates/interface/TriggerFilter.h"

/// default constructor
SingleObjectJetTrigger::SingleObjectJetTrigger(const edm::ParameterSet& cfg)
{
  patTriggerEvent_   = cfg.getParameter<edm::InputTag>("patTriggerEvent");
  patTrigger_        = cfg.getParameter<edm::InputTag>("patTrigger");
  jets_              = cfg.getParameter<edm::InputTag>("jets");
  triggerMatchedJets_= cfg.getParameter<edm::InputTag>("triggerMatchedJets");
  width_             = cfg.getParameter<double>("width");
}

/// calculate the trigger efficiency for QuadJet (or DiJet) triggers based on single object trigger efficiency
double 
SingleObjectJetTrigger::triggerEfficiency(const edm::Handle<pat::JetCollection> jets, const bool diJet)
{
  double effy_ = 1.; //   X0,       X1,       X2,       X3,       X4,       X5,       X6,       X7,       X8,       X9,
  double effies_[100] = { 0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        // 0X
			  0,        0,        0,        0,        0,        0,        0,        0,        0,        0,	      // 1X
			  0.103389, 0.13367,  0.158888, 0.176171, 0.181272, 0.19171,  0.194168, 0.206794, 0.22368,  0.269992, // 2X
			  0.609654, 0.768382, 0.837768, 0.903173, 0.957051, 0.980878, 0.985011, 0.988416, 0.989534, 0.991462, // 3X
			  0.992501, 0.994188, 0.994625, 0.995198, 0.995896, 0.99635,  0.997574, 0.997853, 0.99841,  0.998711, // 4X
			  0.99913,  0.999176, 0.999175, 0.999411, 0.999525, 0.999681, 0.99987,  0.999767, 0.999799, 0.999829, // 5X
			  0.999792, 0.99993,  0.999892, 0.999964, 1,        0.999924, 1,        0.99996,  1,        1,	      // 6X
			  1,        0.999957, 1,        0.999955, 0.999955, 1,        1,        0.999951, 1,        0.999948, // 7X
			  1,        0.999946, 1,        1,        1,        1,        1,        0.999875, 1,        1,	      // 8X
			  0.999863, 1,        1,        1,        1,        1,        1,        1,        1,        1 };      // 9X
  int ptBin = 100;

  double result_i = 0.;
  bool result_i_set = false;
  
  if(diJet){

    for(unsigned int i = 0; i < jets->size(); ++i){
    
      ptBin = floor(jets->at(i).pt());
      if(!result_i_set) { result_i_set = true; result_i = 1.; }
      ptBin < 100 ? result_i *= ( 1. - effies_[ptBin] ) : result_i *= 0.;
    }
    effy_ -= result_i;

    for(unsigned int i = 0; i < jets->size(); ++i){
    
      double p_i = 1.;
      ptBin = floor(jets->at(i).pt());
      if(ptBin < 100) p_i = effies_[ptBin];

      double result_k = 0.;
      bool result_k_set = false;
      for(unsigned int k = 0; k < jets->size(); ++k){
    
	if(k != i){
	  ptBin = floor(jets->at(k).pt());
	  if(!result_k_set) { result_k_set = true; result_k = 1.; }
	  ptBin < 100 ? result_k *= ( 1. - effies_[ptBin] ) : result_k *= 0.;
	}
      }
      effy_ -= p_i * result_k;
    }
    
    return effy_;
  }

  for(unsigned int i = 0; i < jets->size(); ++i){
    
    ptBin = floor(jets->at(i).pt());
    if(!result_i_set) { result_i_set = true; result_i = 1.; }
    ptBin < 100 ? result_i *= ( 1. - effies_[ptBin] ) : result_i *= 0.;
  }
  effy_ -= result_i;

  for(unsigned int i = 0; i < jets->size(); ++i){
    
    double p_i = 1.;
    ptBin = floor(jets->at(i).pt());
    if(ptBin < 100) p_i = effies_[ptBin];

    double result_k = 0.;
    bool result_k_set = false;
    for(unsigned int k = 0; k < jets->size(); ++k){
    
      if(k != i){
	ptBin = floor(jets->at(k).pt());
	if(!result_k_set) { result_k_set = true; result_k = 1.; }
	ptBin < 100 ? result_k *= ( 1. - effies_[ptBin] ) : result_k *= 0.;
      }
    }
    effy_ -= p_i * result_k;
  }

  for(int i = 0; i < (int)jets->size() - 1; ++i){
    
    double p_i = 1.;
    ptBin = floor(jets->at(i).pt());
    if(ptBin < 100) p_i = effies_[ptBin];

    for(int k = i + 1; k < (int)jets->size(); ++k){
    
      double p_k = 1.;
      ptBin = floor(jets->at(k).pt());
      if(ptBin < 100) p_k = effies_[ptBin];
      
      double result_l = 0.;
      bool result_l_set = false;
      for(int l = 0; l < (int)jets->size(); ++l){
	
	if(l != i && l != k){
	  ptBin = floor(jets->at(l).pt());
	  if(!result_l_set) { result_l_set = true; result_l = 1.; }
	  ptBin < 100 ? result_l *= ( 1. - effies_[ptBin] ) : result_l *= 0.;
	}
      }
      effy_ -= p_i * p_k * result_l;
    }
  }

  for(int i = 0; i < (int)jets->size() - 2; ++i){
    
    double p_i = 1.;
    ptBin = floor(jets->at(i).pt());
    if(ptBin < 100) p_i = effies_[ptBin];

    for(int k = i + 1; k < (int)jets->size() - 1; ++k){
    
      double p_k = 1.;
      ptBin = floor(jets->at(k).pt());
      if(ptBin < 100) p_k = effies_[ptBin];
      
      for(int l = k + 1; l < (int)jets->size(); ++l){
	
	double p_l = 1.;
	ptBin = floor(jets->at(l).pt());
	if(ptBin < 100) p_l = effies_[ptBin];
	
	double result_m = 1.;
	bool result_m_set = false;
	for(int m = 0; m < (int)jets->size(); ++m){
	  
	  if(m != i && m != k && m != l){
	    ptBin = floor(jets->at(m).pt());
	    if(!result_m_set) { result_m_set = true; result_m = 1.; }
	    ptBin < 100 ? result_m *= ( 1. - effies_[ptBin] ) : result_m *= 0.;
	  }
	}
	effy_ -= p_i * p_k * p_l * result_m;
      }
    }
  }

  /*
  for(unsigned int idx = 0; idx <= 4; ++idx){

    if(!fiveJet && idx == 4) break;
    if(jets->size() > idx) ptBin = floor(jets->at(idx).pt()); // dirty little hack, probably needs a fix at some time
    ptBin > 99 ? effy_ *= 1. : effy_ *= effies_[ptBin];
  }
  */
  /*
  effy_ *= 0.5*TMath::Erf((jets->at(0).pt()-29.9592)/2.64442961111111)+0.5;
  effy_ *= 0.5*TMath::Erf((jets->at(1).pt()-29.9592)/2.64442961111111)+0.5;
  effy_ *= 0.5*TMath::Erf((jets->at(2).pt()-29.9592)/2.64442961111111)+0.5;
  effy_ *= 0.5*TMath::Erf((jets->at(3).pt()-29.9592)/2.64442961111111)+0.5;
  if (fiveJet) effy_ *= 0.5*TMath::Erf((jets->at(4).pt()-29.9592)/2.64442961111111)+0.5;
  */
  //std::cout << effy_ << std::endl;
  return effy_;
}

/// calculate the trigger efficiency for QuadJet (or DiJet) triggers based on single object trigger efficiency with smeared trigger jet pt's
double
SingleObjectJetTrigger::triggerEfficiencySmear(const std::vector<double> jets, const bool diJet)
{
  double effy_ = 1.; //   X0,       X1,       X2,       X3,       X4,       X5,       X6,       X7,       X8,       X9,
  double effies_[100] = { 0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        // 0X
			  0,        0,        0,        0,        0,        0,        0,        0,        0,        0,	      // 1X
			  0.103389, 0.13367,  0.158888, 0.176171, 0.181272, 0.19171,  0.194168, 0.206794, 0.22368,  0.269992, // 2X
			  0.609654, 0.768382, 0.837768, 0.903173, 0.957051, 0.980878, 0.985011, 0.988416, 0.989534, 0.991462, // 3X
			  0.992501, 0.994188, 0.994625, 0.995198, 0.995896, 0.99635,  0.997574, 0.997853, 0.99841,  0.998711, // 4X
			  0.99913,  0.999176, 0.999175, 0.999411, 0.999525, 0.999681, 0.99987,  0.999767, 0.999799, 0.999829, // 5X
			  0.999792, 0.99993,  0.999892, 0.999964, 1,        0.999924, 1,        0.99996,  1,        1,	      // 6X
			  1,        0.999957, 1,        0.999955, 0.999955, 1,        1,        0.999951, 1,        0.999948, // 7X
			  1,        0.999946, 1,        1,        1,        1,        1,        0.999875, 1,        1,	      // 8X
			  0.999863, 1,        1,        1,        1,        1,        1,        1,        1,        1 };      // 9X
  int ptBin = 100;

  double result_i = 0.;
  bool result_i_set = false;

  if(diJet){

    for(unsigned int i = 0; i < jets.size(); ++i){
    
      ptBin = floor(jets.at(i));
      if(!result_i_set) { result_i_set = true; result_i = 1.; }
      ptBin < 100 ? result_i *= ( 1. - effies_[ptBin] ) : result_i *= 0.;
    }
    effy_ -= result_i;

    for(unsigned int i = 0; i < jets.size(); ++i){
    
      double p_i = 1.;
      ptBin = floor(jets.at(i));
      if(ptBin < 100) p_i = effies_[ptBin];

      double result_k = 0.;
      bool result_k_set = false;
      for(unsigned int k = 0; k < jets.size(); ++k){
    
	if(k != i){
	  ptBin = floor(jets.at(k));
	  if(!result_k_set) { result_k_set = true; result_k = 1.; }
	  ptBin < 100 ? result_k *= ( 1. - effies_[ptBin] ) : result_k *= 0.;
	}
      }
      effy_ -= p_i * result_k;
    }
    return effy_;
  }

  for(unsigned int i = 0; i < jets.size(); ++i){
    
    ptBin = floor(jets.at(i));
    if(!result_i_set) { result_i_set = true; result_i = 1.; }
    ptBin < 100 ? result_i *= ( 1. - effies_[ptBin] ) : result_i *= 0.;
  }
  effy_ -= result_i;

  for(unsigned int i = 0; i < jets.size(); ++i){
    
    double p_i = 1.;
    ptBin = floor(jets.at(i));
    if(ptBin < 100) p_i = effies_[ptBin];

    double result_k = 0.;
    bool result_k_set = false;
    for(unsigned int k = 0; k < jets.size(); ++k){
    
      if(k != i){
	ptBin = floor(jets.at(k));
	if(!result_k_set) { result_k_set = true; result_k = 1.; }
	ptBin < 100 ? result_k *= ( 1. - effies_[ptBin] ) : result_k *= 0.;
      }
    }
    effy_ -= p_i * result_k;
  }

  for(int i = 0; i < (int)jets.size() - 1; ++i){
    
    double p_i = 1.;
    ptBin = floor(jets.at(i));
    if(ptBin < 100) p_i = effies_[ptBin];

    for(int k = i + 1; k < (int)jets.size(); ++k){
    
      double p_k = 1.;
      ptBin = floor(jets.at(k));
      if(ptBin < 100) p_k = effies_[ptBin];
      
      double result_l = 0.;
      bool result_l_set = false;
      for(int l = 0; l < (int)jets.size(); ++l){
	
	if(l != i && l != k){
	  ptBin = floor(jets.at(l));
	  if(!result_l_set) { result_l_set = true; result_l = 1.; }
	  ptBin < 100 ? result_l *= ( 1. - effies_[ptBin] ) : result_l *= 0.;
	}
      }
      effy_ -= p_i * p_k * result_l;
    }
  }

  for(int i = 0; i < (int)jets.size() - 2; ++i){
    
    double p_i = 1.;
    ptBin = floor(jets.at(i));
    if(ptBin < 100) p_i = effies_[ptBin];

    for(int k = i + 1; k < (int)jets.size() - 1; ++k){
    
      double p_k = 1.;
      ptBin = floor(jets.at(k));
      if(ptBin < 100) p_k = effies_[ptBin];
      
      for(int l = k + 1; l < (int)jets.size(); ++l){
	
	double p_l = 1.;
	ptBin = floor(jets.at(l));
	if(ptBin < 100) p_l = effies_[ptBin];
	
	double result_m = 1.;
	bool result_m_set = false;
	for(int m = 0; m < (int)jets.size(); ++m){
	  
	  if(m != i && m != k && m != l){
	    ptBin = floor(jets.at(m));
	    if(!result_m_set) { result_m_set = true; result_m = 1.; }
	    ptBin < 100 ? result_m *= ( 1. - effies_[ptBin] ) : result_m *= 0.;
	  }
	}
	effy_ -= p_i * p_k * p_l * result_m;
      }
    }
  }

  /*
  for(unsigned int idx = 0; idx <= 4; ++idx){

    if(!fiveJet && idx == 4) break;
    if(jets.size() > idx) ptBin = floor(jets.at(idx)->pt()); // dirty little hack, probably needs a fix at some time
    ptBin > 99 ? effy_ *= 1. : effy_ *= effies_[ptBin];
  }
  */
  /*
  effy_ *= 0.5*TMath::Erf((jets.at(0)->pt()-29.9592)/2.64442961111111)+0.5;
  effy_ *= 0.5*TMath::Erf((jets.at(1)->pt()-29.9592)/2.64442961111111)+0.5;
  effy_ *= 0.5*TMath::Erf((jets.at(2)->pt()-29.9592)/2.64442961111111)+0.5;
  effy_ *= 0.5*TMath::Erf((jets.at(3)->pt()-29.9592)/2.64442961111111)+0.5;
  if (fiveJet) effy_ *= 0.5*TMath::Erf((jets.at(4)->pt()-29.9592)/2.64442961111111)+0.5;
  */
  //std::cout << effy_ << std::endl;
  return effy_;
}

/// calculate distance to nearest jet
double
SingleObjectJetTrigger::minDR(const edm::Handle<pat::JetCollection> jets, const unsigned int iJet)
{
  double jet_phi = (*jets)[iJet].phi();
  double jet_eta = (*jets)[iJet].eta();
  double minDR = 9999.;
  for(unsigned int idx = 0; idx < jets->size(); ++idx){
    if ( idx == iJet ) continue;
    if ( deltaR( (*jets)[idx].eta() , (*jets)[idx].phi() , jet_eta, jet_phi ) < minDR ) minDR = deltaR( (*jets)[idx].eta() , (*jets)[idx].phi() , jet_eta, jet_phi );
  }
  return minDR;
}

/// initiate histograms
void
SingleObjectJetTrigger::beginJob()
{
  edm::Service<TFileService> fs;
  if( !fs ) throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

  // pt jets
  hists_["pt" ] = fs->make<TH1F>( "pt"  , "pt"  , 100,  0. , 100. );
  // pt 1. jet
  hists_["pt1"] = fs->make<TH1F>( "pt1" , "pt1" , 100,  0. , 100. );
  // pt 2. jet
  hists_["pt2"] = fs->make<TH1F>( "pt2" , "pt2" , 100,  0. , 100. );
  // pt 3. jet
  hists_["pt3"] = fs->make<TH1F>( "pt3" , "pt3" , 100,  0. , 100. );
  // pt 4. jet
  hists_["pt4"] = fs->make<TH1F>( "pt4" , "pt4" , 100,  0. , 100. );
  // pt 5. jet
  hists_["pt5"] = fs->make<TH1F>( "pt5" , "pt5" , 100,  0. , 100. );
  // pt 6. jet
  hists_["pt6"] = fs->make<TH1F>( "pt6" , "pt6" , 100,  0. , 100. );
  // jet multiplicity
  hists_["n"] = fs->make<TH1F>( "n" , "n" , 10,  -0.5 , 9.5 );
  // pt jets
  hists_["pt_match" ] = fs->make<TH1F>( "pt_match"  , "pt_match"  , 100,  0. , 100. );
  // pt 1. jet
  hists_["pt1_match"] = fs->make<TH1F>( "pt1_match" , "pt1_match" , 100,  0. , 100. );
  // pt 2. jet
  hists_["pt2_match"] = fs->make<TH1F>( "pt2_match" , "pt2_match" , 100,  0. , 100. );
  // pt 3. jet
  hists_["pt3_match"] = fs->make<TH1F>( "pt3_match" , "pt3_match" , 100,  0. , 100. );
  // pt 4. jet
  hists_["pt4_match"] = fs->make<TH1F>( "pt4_match" , "pt4_match" , 100,  0. , 100. );
  // pt 5. jet
  hists_["pt5_match"] = fs->make<TH1F>( "pt5_match" , "pt5_match" , 100,  0. , 100. );
  // pt 6. jet
  hists_["pt6_match"] = fs->make<TH1F>( "pt6_match" , "pt6_match" , 100,  0. , 100. );
  // pt jets
  hists_["pt_trigger_match" ] = fs->make<TH1F>( "pt_trigger_match"  , "pt_trigger_match"  , 100,  0. , 100. );
  // pt 1. jet
  hists_["pt1_trigger_match"] = fs->make<TH1F>( "pt1_trigger_match" , "pt1_trigger_match" , 100,  0. , 100. );
  // pt 2. jet
  hists_["pt2_trigger_match"] = fs->make<TH1F>( "pt2_trigger_match" , "pt2_trigger_match" , 100,  0. , 100. );
  // pt 3. jet
  hists_["pt3_trigger_match"] = fs->make<TH1F>( "pt3_trigger_match" , "pt3_trigger_match" , 100,  0. , 100. );
  // pt 4. jet
  hists_["pt4_trigger_match"] = fs->make<TH1F>( "pt4_trigger_match" , "pt4_trigger_match" , 100,  0. , 100. );
  // pt 5. jet
  hists_["pt5_trigger_match"] = fs->make<TH1F>( "pt5_trigger_match" , "pt5_trigger_match" , 100,  0. , 100. );
  // pt 6. jet
  hists_["pt6_trigger_match"] = fs->make<TH1F>( "pt6_trigger_match" , "pt6_trigger_match" , 100,  0. , 100. );
  // pt trigger jets
  hists_["pt_trigger" ] = fs->make<TH1F>( "pt_trigger"  , "pt_trigger"  , 100,  0. , 100. );
  // pt 1. trigger jet
  hists_["pt1_trigger"] = fs->make<TH1F>( "pt1_trigger" , "pt1_trigger" , 100,  0. , 100. );
  // pt 2. trigger jet
  hists_["pt2_trigger"] = fs->make<TH1F>( "pt2_trigger" , "pt2_trigger" , 100,  0. , 100. );
  // pt 3. trigger jet
  hists_["pt3_trigger"] = fs->make<TH1F>( "pt3_trigger" , "pt3_trigger" , 100,  0. , 100. );
  // pt 4. trigger jet
  hists_["pt4_trigger"] = fs->make<TH1F>( "pt4_trigger" , "pt4_trigger" , 100,  0. , 100. );
  // pt 5. trigger jet
  hists_["pt5_trigger"] = fs->make<TH1F>( "pt5_trigger" , "pt5_trigger" , 100,  0. , 100. );
  // pt 6. trigger jet
  hists_["pt6_trigger"] = fs->make<TH1F>( "pt6_trigger" , "pt6_trigger" , 100,  0. , 100. );
  // trigger jet multiplicity
  hists_["n_trigger"] = fs->make<TH1F>( "n_trigger" , "n_trigger" , 10,  -0.5 , 9.5 );
  // pt jets
  hists_["pt_effy" ] = fs->make<TH1F>( "pt_effy"  , "pt_effy"  , 100,  0. , 100. );
  // pt 1. jet
  hists_["pt1_effy"] = fs->make<TH1F>( "pt1_effy" , "pt1_effy" , 100,  0. , 100. );
  // pt 2. jet
  hists_["pt2_effy"] = fs->make<TH1F>( "pt2_effy" , "pt2_effy" , 100,  0. , 100. );
  // pt 3. jet
  hists_["pt3_effy"] = fs->make<TH1F>( "pt3_effy" , "pt3_effy" , 100,  0. , 100. );
  // pt 4. jet
  hists_["pt4_effy"] = fs->make<TH1F>( "pt4_effy" , "pt4_effy" , 100,  0. , 100. );
  // pt 5. jet
  hists_["pt5_effy"] = fs->make<TH1F>( "pt5_effy" , "pt5_effy" , 100,  0. , 100. );
  // pt 6. jet
  hists_["pt6_effy"] = fs->make<TH1F>( "pt6_effy" , "pt6_effy" , 100,  0. , 100. );
  // Trigger vs. offline p_{T} (GeV)
  hists2D_["ptTrigOff" ] = fs->make<TH2F>( "ptTrigOff" , "Trigger vs. offline p_{T} (GeV)", 60,  0., 300., 60,  0., 300. );
  // Trigger vs. offline eta
  hists2D_["etaTrigOff"] = fs->make<TH2F>( "etaTrigOff", "Trigger vs. offline #eta"       , 50, -2.5, 2.5, 50, -2.5, 2.5 );
  // Trigger vs. offline phi
  hists2D_["phiTrigOff"] = fs->make<TH2F>( "phiTrigOff", "Trigger vs. offline #phi"       , 64, -3.2, 3.2, 64, -3.2, 3.2 );
  // eta of 1. jet pt < 40 GeV, triggered pt > 40 GeV
  hists_["eta_1"] = fs->make<TH1F>( "eta_1" , "eta_1" , 100,  -5. , 5. );
  // eta of 2. jet pt < 40 GeV, triggered pt > 40 GeV
  hists_["eta_2"] = fs->make<TH1F>( "eta_2" , "eta_2" , 100,  -5. , 5. );
  // eta of 3. jet pt < 40 GeV, triggered pt > 40 GeV
  hists_["eta_3"] = fs->make<TH1F>( "eta_3" , "eta_3" , 100,  -5. , 5. );
  // eta of 4. jet pt < 40 GeV, triggered pt > 40 GeV
  hists_["eta_4"] = fs->make<TH1F>( "eta_4" , "eta_4" , 100,  -5. , 5. );
  // eta of jets pt < 40 GeV, triggered pt > 40 GeV
  hists_["eta"] = fs->make<TH1F>( "eta" , "eta" , 100,  -5. , 5. );
  // eta of jets pt < 30 GeV, triggered pt > 30 GeV
  hists_["eta_l"] = fs->make<TH1F>( "eta_l" , "eta_l" , 100,  -5. , 5. );
  // eta of jets pt < 30 GeV, triggered pt > 30 GeV
  hists_["eta_h"] = fs->make<TH1F>( "eta_h" , "eta_h" , 100,  -5. , 5. );
  // phi of jets pt < 30 GeV, triggered pt > 30 GeV
  hists_["phi_l"] = fs->make<TH1F>( "phi_l" , "phi_l" , 63,  -3.15 , 3.15 );
  // phi of jets pt < 30 GeV, triggered pt > 30 GeV
  hists_["phi_h"] = fs->make<TH1F>( "phi_h" , "phi_h" , 63,  -3.15 , 3.15 );
  // emf of 1. jet pt < 40 GeV, triggered pt > 40 GeV
  hists_["emf_1"] = fs->make<TH1F>( "emf_1" , "emf_1" , 100,  0. , 1. );
  // emf of 2. jet pt < 40 GeV, triggered pt > 40 GeV
  hists_["emf_2"] = fs->make<TH1F>( "emf_2" , "emf_2" , 100,  0. , 1. );
  // emf of 3. jet pt < 40 GeV, triggered pt > 40 GeV
  hists_["emf_3"] = fs->make<TH1F>( "emf_3" , "emf_3" , 100,  0. , 1. );
  // emf of 4. jet pt < 40 GeV, triggered pt > 40 GeV
  hists_["emf_4"] = fs->make<TH1F>( "emf_4" , "emf_4" , 100,  0. , 1. );
  // emf of jets pt < 40 GeV, triggered pt > 40 GeV
  hists_["emf"] = fs->make<TH1F>( "emf" , "emf" , 100,  0. , 1. );
  // emf of jets pt < 30 GeV, triggered pt > 30 GeV
  hists_["emf_l"] = fs->make<TH1F>( "emf_l" , "emf_l" , 100,  0. , 1. );
  // emf of jets pt < 30 GeV, triggered pt > 30 GeV
  hists_["emf_h"] = fs->make<TH1F>( "emf_h" , "emf_h" , 100,  0. , 1. );
  // minDR of 1. jet pt < 40 GeV, triggered pt > 40 GeV
  hists_["minDR_1"] = fs->make<TH1F>( "minDR_1" , "minDR_1" , 150,  0. , 3. );
  // minDR of 2. jet pt < 40 GeV, triggered pt > 40 GeV
  hists_["minDR_2"] = fs->make<TH1F>( "minDR_2" , "minDR_2" , 150,  0. , 3. );
  // minDR of 3. jet pt < 40 GeV, triggered pt > 40 GeV
  hists_["minDR_3"] = fs->make<TH1F>( "minDR_3" , "minDR_3" , 150,  0. , 3. );
  // minDR of 4. jet pt < 40 GeV, triggered pt > 40 GeV
  hists_["minDR_4"] = fs->make<TH1F>( "minDR_4" , "minDR_4" , 150,  0. , 3. );
  // minDR of jets pt < 40 GeV, triggered pt > 40 GeV
  hists_["minDR"] = fs->make<TH1F>( "minDR" , "minDR" , 150,  0. , 3. );
  // minDR of jets pt < 30 GeV, triggered pt > 30 GeV
  hists_["minDR_l"] = fs->make<TH1F>( "minDR_l" , "minDR_l" , 150,  0. , 3. );
  // minDR of jets pt > 30 GeV, triggered pt > 30 GeV
  hists_["minDR_h"] = fs->make<TH1F>( "minDR_h" , "minDR_h" , 150,  0. , 3. );
  // dR of 1. jet pt < 30 GeV to triggered jet pt > 30 GeV
  hists_["dR2Trig_1"] = fs->make<TH1F>( "dR2Trig_1" , "dR2Trig_1" , 50,  0. , 0.5 );
  // dR of 2. jet pt < 30 GeV to triggered jet pt > 30 GeV
  hists_["dR2Trig_2"] = fs->make<TH1F>( "dR2Trig_2" , "dR2Trig_2" , 50,  0. , 0.5 );
  // dR of 3. jet pt < 30 GeV to triggered jet pt > 30 GeV
  hists_["dR2Trig_3"] = fs->make<TH1F>( "dR2Trig_3" , "dR2Trig_3" , 50,  0. , 0.5 );
  // dR of 4. jet pt < 30 GeV to triggered jet pt > 30 GeV
  hists_["dR2Trig_4"] = fs->make<TH1F>( "dR2Trig_4" , "dR2Trig_4" , 50,  0. , 0.5 );
  // dR of jets pt < 30 GeV to triggered jet pt > 30 GeV
  hists_["dR2Trig_l"] = fs->make<TH1F>( "dR2Trig_l" , "dR2Trig_l" , 50,  0. , 0.5 );
  // dR of jets pt > 30 GeV to triggered jet pt > 30 GeV
  hists_["dR2Trig_h"] = fs->make<TH1F>( "dR2Trig_h" , "dR2Trig_h" , 50,  0. , 0.5 );
  // number of constituents of jets pt < 30 GeV, triggered jet pt > 30 GeV
  hists_["nConst_l"] = fs->make<TH1F>( "nConst_l" , "nConst_l" , 50,  0. , 50. );
  // number of constituents of jets 30 GeV < pt < 40 GeV, triggered jet pt > 30 GeV
  hists_["nConst_m1"] = fs->make<TH1F>( "nConst_m1" , "nConst_m1" , 50,  0. , 50. );
  // number of constituents of jets 40 GeV < pt < 50 GeV, triggered jet pt > 30 GeV
  hists_["nConst_m2"] = fs->make<TH1F>( "nConst_m2" , "nConst_m2" , 50,  0. , 50. );
  // number of constituents of jets pt > 50 GeV, triggered jet pt > 30 GeV
  hists_["nConst_h"] = fs->make<TH1F>( "nConst_h" , "nConst_h" , 50,  0. , 50. );


  /**
     ansatz of depending probabilities
  **/
 
  // trigger matched 2. jets if trigger matched 1. jet above 40 GeV
  hists_["depPt2"] = fs->make<TH1F>( "depPt2" , "depPt2" , 100 , 0. , 100. );
  // trigger matched 3. jets if trigger matched 2. jet above 40 GeV
  hists_["depPt3"] = fs->make<TH1F>( "depPt3" , "depPt3" , 100 , 0. , 100. );
  // trigger matched 4. jets if trigger matched 3. jet above 40 GeV
  hists_["depPt4"] = fs->make<TH1F>( "depPt4" , "depPt4" , 100 , 0. , 100. );


  /**
     pt distribution for triggered / untriggered events which pass the SOT
  **/

  // QuadJet30 pt distribution of 1. jet for triggered / untriggered events
  hists_["QJ30_1_trig"] = fs->make<TH1F>( "QJ30_1_trig" , "QuadJet30_1_triggered" , 100,  0. , 100. );
  hists_["QJ30_1_untrig"] = fs->make<TH1F>( "QJ30_1_untrig" , "QuadJet30_1_untriggered" , 100,  0. , 100. );
  // QuadJet30 pt distribution of 2. jet for triggered / untriggered events
  hists_["QJ30_2_trig"] = fs->make<TH1F>( "QJ30_2_trig" , "QuadJet30_2_triggered" , 100,  0. , 100. );
  hists_["QJ30_2_untrig"] = fs->make<TH1F>( "QJ30_2_untrig" , "QuadJet30_2_untriggered" , 100,  0. , 100. );
  // QuadJet30 pt distribution of 3. jet for triggered / untriggered events
  hists_["QJ30_3_trig"] = fs->make<TH1F>( "QJ30_3_trig" , "QuadJet30_3_triggered" , 100,  0. , 100. );
  hists_["QJ30_3_untrig"] = fs->make<TH1F>( "QJ30_3_untrig" , "QuadJet30_3_untriggered" , 100,  0. , 100. );
  // QuadJet30 pt distribution of 4. jet for triggered / untriggered events
  hists_["QJ30_4_trig"] = fs->make<TH1F>( "QJ30_4_trig" , "QuadJet30_5_triggered" , 100,  0. , 100. );
  hists_["QJ30_4_untrig"] = fs->make<TH1F>( "QJ30_4_untrig" , "QuadJet30_4_untriggered" , 100,  0. , 100. );
  // QuadJet40 pt distribution of 1. jet for triggered / untriggered events
  hists_["QJ40_1_trig"] = fs->make<TH1F>( "QJ40_1_trig" , "QuadJet40_1_triggered" , 100,  0. , 100. );
  hists_["QJ40_1_untrig"] = fs->make<TH1F>( "QJ40_1_untrig" , "QuadJet40_1_untriggered" , 100,  0. , 100. );
  // QuadJet40 pt distribution of 2. jet for triggered / untriggered events
  hists_["QJ40_2_trig"] = fs->make<TH1F>( "QJ40_2_trig" , "QuadJet40_2_triggered" , 100,  0. , 100. );
  hists_["QJ40_2_untrig"] = fs->make<TH1F>( "QJ40_2_untrig" , "QuadJet40_2_untriggered" , 100,  0. , 100. );
  // QuadJet40 pt distribution of 3. jet for triggered / untriggered events
  hists_["QJ40_3_trig"] = fs->make<TH1F>( "QJ40_3_trig" , "QuadJet40_3_triggered" , 100,  0. , 100. );
  hists_["QJ40_3_untrig"] = fs->make<TH1F>( "QJ40_3_untrig" , "QuadJet40_3_untriggered" , 100,  0. , 100. );
  // QuadJet40 pt distribution of 4. jet for triggered / untriggered events
  hists_["QJ40_4_trig"] = fs->make<TH1F>( "QJ40_4_trig" , "QuadJet40_5_triggered" , 100,  0. , 100. );
  hists_["QJ40_4_untrig"] = fs->make<TH1F>( "QJ40_4_untrig" , "QuadJet40_4_untriggered" , 100,  0. , 100. );

  /**
     Efficiencies vs. reco pt
  **/

  // QuadJet30 efficiency for 1. jet and normalization histogram
  hists_["QJ30_1"]   = fs->make<TH1F>( "QJ30_1"   , "QuadJet30_1"      , 100,  0. , 100. );
  hists_["QJ30_1_n"] = fs->make<TH1F>( "QJ30_1_n" , "QuadJet30_1_norm" , 100,  0. , 100. );
  hists_["QJ30_1_e"] = fs->make<TH1F>( "QJ30_1_e" , "QuadJet30_1_effy" , 100,  0. , 100. );
  // QuadJet30 efficiency for 2. jet and normalization histogram
  hists_["QJ30_2"]   = fs->make<TH1F>( "QJ30_2"   , "QuadJet30_2"      , 100,  0. , 100. );
  hists_["QJ30_2_n"] = fs->make<TH1F>( "QJ30_2_n" , "QuadJet30_2_norm" , 100,  0. , 100. );
  hists_["QJ30_2_e"] = fs->make<TH1F>( "QJ30_2_e" , "QuadJet30_2_effy" , 100,  0. , 100. );
  // QuadJet30 efficiency for 3. jet and normalization histogram
  hists_["QJ30_3"]   = fs->make<TH1F>( "QJ30_3"   , "QuadJet30_3"      , 100,  0. , 100. );
  hists_["QJ30_3_n"] = fs->make<TH1F>( "QJ30_3_n" , "QuadJet30_3_norm" , 100,  0. , 100. );
  hists_["QJ30_3_e"] = fs->make<TH1F>( "QJ30_3_e" , "QuadJet30_3_effy" , 100,  0. , 100. );
  // QuadJet30 efficiency for 4. jet and normalization histogram
  hists_["QJ30_4"]   = fs->make<TH1F>( "QJ30_4"   , "QuadJet30_4"      , 100,  0. , 100. );
  hists_["QJ30_4_n"] = fs->make<TH1F>( "QJ30_4_n" , "QuadJet30_5_norm" , 100,  0. , 100. );
  hists_["QJ30_4_e"] = fs->make<TH1F>( "QJ30_4_e" , "QuadJet30_4_effy" , 100,  0. , 100. );
  // QuadJet40 efficiency for 1. jet and normalization histogram
  hists_["QJ40_1"]   = fs->make<TH1F>( "QJ40_1"   , "QuadJet40_1"      , 100,  0. , 100. );
  hists_["QJ40_1_n"] = fs->make<TH1F>( "QJ40_1_n" , "QuadJet40_1_norm" , 100,  0. , 100. );
  hists_["QJ40_1_e"] = fs->make<TH1F>( "QJ40_1_e" , "QuadJet40_1_effy" , 100,  0. , 100. );
  // QuadJet40 efficiency for 2. jet and normalization histogram
  hists_["QJ40_2"]   = fs->make<TH1F>( "QJ40_2"   , "QuadJet40_2"      , 100,  0. , 100. );
  hists_["QJ40_2_n"] = fs->make<TH1F>( "QJ40_2_n" , "QuadJet40_2_norm" , 100,  0. , 100. );
  hists_["QJ40_2_e"] = fs->make<TH1F>( "QJ40_2_e" , "QuadJet40_2_effy" , 100,  0. , 100. );
  // QuadJet40 efficiency for 3. jet and normalization histogram
  hists_["QJ40_3"]   = fs->make<TH1F>( "QJ40_3"   , "QuadJet40_3"      , 100,  0. , 100. );
  hists_["QJ40_3_n"] = fs->make<TH1F>( "QJ40_3_n" , "QuadJet40_3_norm" , 100,  0. , 100. );
  hists_["QJ40_3_e"] = fs->make<TH1F>( "QJ40_3_e" , "QuadJet40_3_effy" , 100,  0. , 100. );
  // QuadJet40 efficiency for 4. jet and normalization histogram
  hists_["QJ40_4"]   = fs->make<TH1F>( "QJ40_4"   , "QuadJet40_4"      , 100,  0. , 100. );
  hists_["QJ40_4_n"] = fs->make<TH1F>( "QJ40_4_n" , "QuadJet40_5_norm" , 100,  0. , 100. );
  hists_["QJ40_4_e"] = fs->make<TH1F>( "QJ40_4_e" , "QuadJet40_4_effy" , 100,  0. , 100. );

  /**
     imaginary new self defined DiJet40 trigger for testing purpose only
  **/

  // DiJet40 efficiency for 1. jet and normalization histogram
  hists_["DJ40_1"]   = fs->make<TH1F>( "DJ40_1"   , "DiJet40_1"      , 100,  0. , 100. );
  hists_["DJ40_1_n"] = fs->make<TH1F>( "DJ40_1_n" , "DiJet40_1_norm" , 100,  0. , 100. );
  hists_["DJ40_1_e"] = fs->make<TH1F>( "DJ40_1_e" , "DiJet40_1_effy" , 100,  0. , 100. );
  // DiJet40 efficiency for 2. jet and normalization histogram
  hists_["DJ40_2"]   = fs->make<TH1F>( "DJ40_2"   , "DiJet40_2"      , 100,  0. , 100. );
  hists_["DJ40_2_n"] = fs->make<TH1F>( "DJ40_2_n" , "DiJet40_2_norm" , 100,  0. , 100. );
  hists_["DJ40_2_e"] = fs->make<TH1F>( "DJ40_2_e" , "DiJet40_2_effy" , 100,  0. , 100. );

  /*
  // QuadJet80303030 efficiency for 1. jet and normalization histogram
  hists_["QJ8333_1"]   = fs->make<TH1F>( "QJ8333_1"   , "QuadJet80303030_1"      , 100,  0. , 100. );
  hists_["QJ8333_1_n"] = fs->make<TH1F>( "QJ8333_1_n" , "QuadJet80303030_1_norm" , 100,  0. , 100. );
  hists_["QJ8333_1_e"] = fs->make<TH1F>( "QJ8333_1_e" , "QuadJet80303030_1_effy" , 100,  0. , 100. );
  // QuadJet80303030 efficiency for 2. jet and normalization histogram
  hists_["QJ8333_2"]   = fs->make<TH1F>( "QJ8333_2"   , "QuadJet80303030_2"      , 100,  0. , 100. );
  hists_["QJ8333_2_n"] = fs->make<TH1F>( "QJ8333_2_n" , "QuadJet80303030_2_norm" , 100,  0. , 100. );
  hists_["QJ8333_2_e"] = fs->make<TH1F>( "QJ8333_2_e" , "QuadJet80303030_2_effy" , 100,  0. , 100. );
  // QuadJet80303030 efficiency for 3. jet and normalization histogram
  hists_["QJ8333_3"]   = fs->make<TH1F>( "QJ8333_3"   , "QuadJet80303030_3"      , 100,  0. , 100. );
  hists_["QJ8333_3_n"] = fs->make<TH1F>( "QJ8333_3_n" , "QuadJet80303030_3_norm" , 100,  0. , 100. );
  hists_["QJ8333_3_e"] = fs->make<TH1F>( "QJ8333_3_e" , "QuadJet80303030_3_effy" , 100,  0. , 100. );
  // QuadJet80303030 efficiency for 4. jet and normalization histogram
  hists_["QJ8333_4"]   = fs->make<TH1F>( "QJ8333_4"   , "QuadJet80303030_4"      , 100,  0. , 100. );
  hists_["QJ8333_4_n"] = fs->make<TH1F>( "QJ8333_4_n" , "QuadJet80303030_5_norm" , 100,  0. , 100. );
  hists_["QJ8333_4_e"] = fs->make<TH1F>( "QJ8333_4_e" , "QuadJet80303030_4_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 1. jet and normalization histogram
  hists_["FJ30_1"]   = fs->make<TH1F>( "FJ30_1"   , "FiveJet30_1"      , 100,  0. , 100. );
  hists_["FJ30_1_n"] = fs->make<TH1F>( "FJ30_1_n" , "FiveJet30_1_norm" , 100,  0. , 100. );
  hists_["FJ30_1_e"] = fs->make<TH1F>( "FJ30_1_e" , "FiveJet30_1_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 2. jet and normalization histogram
  hists_["FJ30_2"]   = fs->make<TH1F>( "FJ30_2"   , "FiveJet30_2"      , 100,  0. , 100. );
  hists_["FJ30_2_n"] = fs->make<TH1F>( "FJ30_2_n" , "FiveJet30_2_norm" , 100,  0. , 100. );
  hists_["FJ30_2_e"] = fs->make<TH1F>( "FJ30_2_e" , "FiveJet30_2_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 3. jet and normalization histogram
  hists_["FJ30_3"]   = fs->make<TH1F>( "FJ30_3"   , "FiveJet30_3"      , 100,  0. , 100. );
  hists_["FJ30_3_n"] = fs->make<TH1F>( "FJ30_3_n" , "FiveJet30_3_norm" , 100,  0. , 100. );
  hists_["FJ30_3_e"] = fs->make<TH1F>( "FJ30_3_e" , "FiveJet30_3_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 4. jet and normalization histogram
  hists_["FJ30_4"]   = fs->make<TH1F>( "FJ30_4"   , "FiveJet30_4"      , 100,  0. , 100. );
  hists_["FJ30_4_n"] = fs->make<TH1F>( "FJ30_4_n" , "FiveJet30_4_norm" , 100,  0. , 100. );
  hists_["FJ30_4_e"] = fs->make<TH1F>( "FJ30_4_e" , "FiveJet30_4_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 5. jet and normalization histogram
  hists_["FJ30_5"]   = fs->make<TH1F>( "FJ30_5"   , "FiveJet30_5"      , 100,  0. , 100. );
  hists_["FJ30_5_n"] = fs->make<TH1F>( "FJ30_5_n" , "FiveJet30_5_norm" , 100,  0. , 100. );
  hists_["FJ30_5_e"] = fs->make<TH1F>( "FJ30_5_e" , "FiveJet30_5_effy" , 100,  0. , 100. );
  */
  /**
     Efficiencies vs. trigger pt
  **/

  // QuadJet30 efficiency for 1. jet and normalization histogram
  hists_["QJ30_trigger_1"]   = fs->make<TH1F>( "QJ30_trigger_1"   , "QuadJet30_trigger_1"      , 100,  0. , 100. );
  hists_["QJ30_trigger_1_n"] = fs->make<TH1F>( "QJ30_trigger_1_n" , "QuadJet30_trigger_1_norm" , 100,  0. , 100. );
  hists_["QJ30_trigger_1_e"] = fs->make<TH1F>( "QJ30_trigger_1_e" , "QuadJet30_trigger_1_effy" , 100,  0. , 100. );
  // QuadJet30 efficiency for 2. jet and normalization histogram
  hists_["QJ30_trigger_2"]   = fs->make<TH1F>( "QJ30_trigger_2"   , "QuadJet30_trigger_2"      , 100,  0. , 100. );
  hists_["QJ30_trigger_2_n"] = fs->make<TH1F>( "QJ30_trigger_2_n" , "QuadJet30_trigger_2_norm" , 100,  0. , 100. );
  hists_["QJ30_trigger_2_e"] = fs->make<TH1F>( "QJ30_trigger_2_e" , "QuadJet30_trigger_2_effy" , 100,  0. , 100. );
  // QuadJet30 efficiency for 3. jet and normalization histogram
  hists_["QJ30_trigger_3"]   = fs->make<TH1F>( "QJ30_trigger_3"   , "QuadJet30_trigger_3"      , 100,  0. , 100. );
  hists_["QJ30_trigger_3_n"] = fs->make<TH1F>( "QJ30_trigger_3_n" , "QuadJet30_trigger_3_norm" , 100,  0. , 100. );
  hists_["QJ30_trigger_3_e"] = fs->make<TH1F>( "QJ30_trigger_3_e" , "QuadJet30_trigger_3_effy" , 100,  0. , 100. );
  // QuadJet30 efficiency for 4. jet and normalization histogram
  hists_["QJ30_trigger_4"]   = fs->make<TH1F>( "QJ30_trigger_4"   , "QuadJet30_trigger_4"      , 100,  0. , 100. );
  hists_["QJ30_trigger_4_n"] = fs->make<TH1F>( "QJ30_trigger_4_n" , "QuadJet30_trigger_5_norm" , 100,  0. , 100. );
  hists_["QJ30_trigger_4_e"] = fs->make<TH1F>( "QJ30_trigger_4_e" , "QuadJet30_trigger_4_effy" , 100,  0. , 100. );
  // QuadJet40 efficiency for 1. jet and normalization histogram
  hists_["QJ40_trigger_1"]   = fs->make<TH1F>( "QJ40_trigger_1"   , "QuadJet40_trigger_1"      , 100,  0. , 100. );
  hists_["QJ40_trigger_1_n"] = fs->make<TH1F>( "QJ40_trigger_1_n" , "QuadJet40_trigger_1_norm" , 100,  0. , 100. );
  hists_["QJ40_trigger_1_e"] = fs->make<TH1F>( "QJ40_trigger_1_e" , "QuadJet40_trigger_1_effy" , 100,  0. , 100. );
  // QuadJet40 efficiency for 2. jet and normalization histogram
  hists_["QJ40_trigger_2"]   = fs->make<TH1F>( "QJ40_trigger_2"   , "QuadJet40_trigger_2"      , 100,  0. , 100. );
  hists_["QJ40_trigger_2_n"] = fs->make<TH1F>( "QJ40_trigger_2_n" , "QuadJet40_trigger_2_norm" , 100,  0. , 100. );
  hists_["QJ40_trigger_2_e"] = fs->make<TH1F>( "QJ40_trigger_2_e" , "QuadJet40_trigger_2_effy" , 100,  0. , 100. );
  // QuadJet40 efficiency for 3. jet and normalization histogram
  hists_["QJ40_trigger_3"]   = fs->make<TH1F>( "QJ40_trigger_3"   , "QuadJet40_trigger_3"      , 100,  0. , 100. );
  hists_["QJ40_trigger_3_n"] = fs->make<TH1F>( "QJ40_trigger_3_n" , "QuadJet40_trigger_3_norm" , 100,  0. , 100. );
  hists_["QJ40_trigger_3_e"] = fs->make<TH1F>( "QJ40_trigger_3_e" , "QuadJet40_trigger_3_effy" , 100,  0. , 100. );
  // QuadJet40 efficiency for 4. jet and normalization histogram
  hists_["QJ40_trigger_4"]   = fs->make<TH1F>( "QJ40_trigger_4"   , "QuadJet40_trigger_4"      , 100,  0. , 100. );
  hists_["QJ40_trigger_4_n"] = fs->make<TH1F>( "QJ40_trigger_4_n" , "QuadJet40_trigger_5_norm" , 100,  0. , 100. );
  hists_["QJ40_trigger_4_e"] = fs->make<TH1F>( "QJ40_trigger_4_e" , "QuadJet40_trigger_4_effy" , 100,  0. , 100. );

  /**
     imaginary new self defined DiJet40 trigger for testing purpose only
  **/

  // DiJet40 efficiency for 1. jet and normalization histogram
  hists_["DJ40_trigger_1"]   = fs->make<TH1F>( "DJ40_trigger_1"   , "DiJet40_trigger_1"      , 100,  0. , 100. );
  hists_["DJ40_trigger_1_n"] = fs->make<TH1F>( "DJ40_trigger_1_n" , "DiJet40_trigger_1_norm" , 100,  0. , 100. );
  hists_["DJ40_trigger_1_e"] = fs->make<TH1F>( "DJ40_trigger_1_e" , "DiJet40_trigger_1_effy" , 100,  0. , 100. );
  // DiJet40 efficiency for 2. jet and normalization histogram
  hists_["DJ40_trigger_2"]   = fs->make<TH1F>( "DJ40_trigger_2"   , "DiJet40_trigger_2"      , 100,  0. , 100. );
  hists_["DJ40_trigger_2_n"] = fs->make<TH1F>( "DJ40_trigger_2_n" , "DiJet40_trigger_2_norm" , 100,  0. , 100. );
  hists_["DJ40_trigger_2_e"] = fs->make<TH1F>( "DJ40_trigger_2_e" , "DiJet40_trigger_2_effy" , 100,  0. , 100. );

  /*
  // QuadJet80303030 efficiency for 1. jet and normalization histogram
  hists_["QJ8333_trigger_1"]   = fs->make<TH1F>( "QJ8333_trigger_1"   , "QuadJet80303030_trigger_1"      , 100,  0. , 100. );
  hists_["QJ8333_trigger_1_n"] = fs->make<TH1F>( "QJ8333_trigger_1_n" , "QuadJet80303030_trigger_1_norm" , 100,  0. , 100. );
  hists_["QJ8333_trigger_1_e"] = fs->make<TH1F>( "QJ8333_trigger_1_e" , "QuadJet80303030_trigger_1_effy" , 100,  0. , 100. );
  // QuadJet80303030 efficiency for 2. jet and normalization histogram
  hists_["QJ8333_trigger_2"]   = fs->make<TH1F>( "QJ8333_trigger_2"   , "QuadJet80303030_trigger_2"      , 100,  0. , 100. );
  hists_["QJ8333_trigger_2_n"] = fs->make<TH1F>( "QJ8333_trigger_2_n" , "QuadJet80303030_trigger_2_norm" , 100,  0. , 100. );
  hists_["QJ8333_trigger_2_e"] = fs->make<TH1F>( "QJ8333_trigger_2_e" , "QuadJet80303030_trigger_2_effy" , 100,  0. , 100. );
  // QuadJet80303030 efficiency for 3. jet and normalization histogram
  hists_["QJ8333_trigger_3"]   = fs->make<TH1F>( "QJ8333_trigger_3"   , "QuadJet80303030_trigger_3"      , 100,  0. , 100. );
  hists_["QJ8333_trigger_3_n"] = fs->make<TH1F>( "QJ8333_trigger_3_n" , "QuadJet80303030_trigger_3_norm" , 100,  0. , 100. );
  hists_["QJ8333_trigger_3_e"] = fs->make<TH1F>( "QJ8333_trigger_3_e" , "QuadJet80303030_trigger_3_effy" , 100,  0. , 100. );
  // QuadJet80303030 efficiency for 4. jet and normalization histogram
  hists_["QJ8333_trigger_4"]   = fs->make<TH1F>( "QJ8333_trigger_4"   , "QuadJet80303030_trigger_4"      , 100,  0. , 100. );
  hists_["QJ8333_trigger_4_n"] = fs->make<TH1F>( "QJ8333_trigger_4_n" , "QuadJet80303030_trigger_5_norm" , 100,  0. , 100. );
  hists_["QJ8333_trigger_4_e"] = fs->make<TH1F>( "QJ8333_trigger_4_e" , "QuadJet80303030_trigger_4_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 1. jet and normalization histogram
  hists_["FJ30_trigger_1"]   = fs->make<TH1F>( "FJ30_trigger_1"   , "FiveJet30_trigger_1"      , 100,  0. , 100. );
  hists_["FJ30_trigger_1_n"] = fs->make<TH1F>( "FJ30_trigger_1_n" , "FiveJet30_trigger_1_norm" , 100,  0. , 100. );
  hists_["FJ30_trigger_1_e"] = fs->make<TH1F>( "FJ30_trigger_1_e" , "FiveJet30_trigger_1_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 2. jet and normalization histogram
  hists_["FJ30_trigger_2"]   = fs->make<TH1F>( "FJ30_trigger_2"   , "FiveJet30_trigger_2"      , 100,  0. , 100. );
  hists_["FJ30_trigger_2_n"] = fs->make<TH1F>( "FJ30_trigger_2_n" , "FiveJet30_trigger_2_norm" , 100,  0. , 100. );
  hists_["FJ30_trigger_2_e"] = fs->make<TH1F>( "FJ30_trigger_2_e" , "FiveJet30_trigger_2_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 3. jet and normalization histogram
  hists_["FJ30_trigger_3"]   = fs->make<TH1F>( "FJ30_trigger_3"   , "FiveJet30_trigger_3"      , 100,  0. , 100. );
  hists_["FJ30_trigger_3_n"] = fs->make<TH1F>( "FJ30_trigger_3_n" , "FiveJet30_trigger_3_norm" , 100,  0. , 100. );
  hists_["FJ30_trigger_3_e"] = fs->make<TH1F>( "FJ30_trigger_3_e" , "FiveJet30_trigger_3_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 4. jet and normalization histogram
  hists_["FJ30_trigger_4"]   = fs->make<TH1F>( "FJ30_trigger_4"   , "FiveJet30_trigger_4"      , 100,  0. , 100. );
  hists_["FJ30_trigger_4_n"] = fs->make<TH1F>( "FJ30_trigger_4_n" , "FiveJet30_trigger_4_norm" , 100,  0. , 100. );
  hists_["FJ30_trigger_4_e"] = fs->make<TH1F>( "FJ30_trigger_4_e" , "FiveJet30_trigger_4_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 5. jet and normalization histogram
  hists_["FJ30_trigger_5"]   = fs->make<TH1F>( "FJ30_trigger_5"   , "FiveJet30_trigger_5"      , 100,  0. , 100. );
  hists_["FJ30_trigger_5_n"] = fs->make<TH1F>( "FJ30_trigger_5_n" , "FiveJet30_trigger_5_norm" , 100,  0. , 100. );
  hists_["FJ30_trigger_5_e"] = fs->make<TH1F>( "FJ30_trigger_5_e" , "FiveJet30_trigger_5_effy" , 100,  0. , 100. );
  */
  rnd.SetSeed(0);
}

/// analyze triggers and fill histograms
void
SingleObjectJetTrigger::analyze(const edm::Event& event, const edm::EventSetup&)
{
  //edm::Handle<std::vector<pat::TriggerFilter> > triggerFilters;
  //event.getByLabel(patTrigger_, triggerFilters);
  //edm::Handle<std::vector<pat::TriggerPath> > triggerPaths;
  //event.getByLabel(patTrigger_, triggerPaths);
  //edm::Handle<std::vector<pat::TriggerObject> > triggerObjects;
  //event.getByLabel(patTrigger_, triggerObjects);
  edm::Handle<pat::TriggerEvent> triggerEvent;
  event.getByLabel(patTriggerEvent_, triggerEvent);
  //edm::Handle<edm::View<pat::Jet> > jets;
  edm::Handle<pat::JetCollection> jets;
  event.getByLabel(jets_, jets);
  //edm::Handle<edm::View<pat::Jet> > triggerMatchedJets;
  //event.getByLabel(triggerMatchedJets_, triggerMatchedJets);

  if( !triggerEvent.failedToGet() ) {

    //bool L1Triggered = false;
    bool HLTriggered = false;
    bool singleJetTriggered = false;
    /*
    std::cout << "-------------" << std::endl;
    for(pat::TriggerObjectCollection::const_iterator obj = triggerEvent->objects()->begin(); obj != triggerEvent->objects()->end(); ++obj){
      if( obj->hasFilterId(85) ) std::cout << "HLTriggerJet found with pt = " << obj->pt() << std::endl;
      //std::cout << obj->pdgId() << " :: " << obj->pt()  << " :: " << obj->isJet() << std::endl;
    }
    */
    pat::TriggerFilterRefVector filters = triggerEvent->acceptedFilters();
    for(pat::TriggerFilterRefVector::const_iterator filter = filters.begin(); filter != filters.end(); ++filter){
      if((*filter)->label() == "hlt1jet30") singleJetTriggered = true;
      if((*filter)->label() == "hlt4jet30") HLTriggered = true;
      //if((*filter)->label() == "hltL1sQuadJet30") L1Triggered = true;
   }
    
    //if(L1Triggered){

      // PAT trigger helper for trigger matching information
      const pat::helper::TriggerMatchHelper matchHelper;

      const pat::TriggerObjectMatch * HLTriggerMatch( triggerEvent->triggerObjectMatchResult( "jetTriggerMatchHLTJets" ) );
      unsigned int jetCounter = 1;
      hists_.find("n")->second->Fill( jets->size() );

      double effy = 1.;
      double jet1_pt = 0;
      double jet2_pt = 0;
      double jet3_pt = 0;
      double jet4_pt = 0;
      //double jet5_pt = 0;
      double HLTriggerMatchedJet1_pt = 0;
      double HLTriggerMatchedJet2_pt = 0;
      double HLTriggerMatchedJet3_pt = 0;
      //double HLTriggerMatchedJet4_pt = 0;
      //double HLTriggerMatchedJet5_pt = 0;
      double HLTriggerJet1_pt = 0;
      double HLTriggerJet2_pt = 0;
      double HLTriggerJet3_pt = 0;
      double HLTriggerJet4_pt = 0;
      //double HLTriggerJet5_pt = 0;
    
      for (size_t iJet = 0; iJet < jets->size(); ++iJet, ++jetCounter) {
	double jet_pt = (*jets)[iJet].pt();
	hists_.find("pt")->second->Fill( jet_pt );
	if     (jetCounter == 1) { hists_.find("pt1")->second->Fill( jet_pt ); jet1_pt = jet_pt; }
	else if(jetCounter == 2) { hists_.find("pt2")->second->Fill( jet_pt ); jet2_pt = jet_pt; }
	else if(jetCounter == 3) { hists_.find("pt3")->second->Fill( jet_pt ); jet3_pt = jet_pt; }
	else if(jetCounter == 4) { hists_.find("pt4")->second->Fill( jet_pt ); jet4_pt = jet_pt; }
	else if(jetCounter == 5) { hists_.find("pt5")->second->Fill( jet_pt ); } //jet5_pt = jet_pt; }
	else if(jetCounter == 6)   hists_.find("pt6")->second->Fill( jet_pt );

	const reco::CandidateBaseRef candBaseRef( pat::JetRef( jets , iJet) );
	const pat::TriggerObjectRef HLTrigRef( matchHelper.triggerMatchObject( candBaseRef, HLTriggerMatch, event, *triggerEvent ) );

	if ( HLTrigRef.isAvailable() ) {
	  double HLTriggerMatchedJet_pt = candBaseRef->pt();
	  double HLTriggerJet_pt = HLTrigRef->pt();
	  double dR2Trig = deltaR( candBaseRef->eta() , candBaseRef->phi() , HLTrigRef->eta(), HLTrigRef->phi() );
	  double minDR_ = minDR(jets, iJet);
	  double eta_ = (*jets)[iJet].eta();
	  double phi_ = (*jets)[iJet].phi();
	  double emf_ = (*jets)[iJet].emEnergyFraction();
	  hists_.find("pt_match")->second->Fill( HLTriggerMatchedJet_pt );
	  hists_.find("pt_trigger_match")->second->Fill( HLTriggerJet_pt );
	  hists2D_.find("ptTrigOff")->second->Fill( HLTriggerMatchedJet_pt, HLTriggerJet_pt );
	  hists2D_.find("etaTrigOff")->second->Fill( candBaseRef->eta(), HLTrigRef->eta() );
	  hists2D_.find("phiTrigOff")->second->Fill( candBaseRef->phi(), HLTrigRef->phi() );
	  if ( HLTriggerMatchedJet_pt <  30 ) {
	    hists_.find("dR2Trig_l")->second->Fill( dR2Trig );
	    hists_.find("nConst_l")->second->Fill( (*jets)[iJet].nConstituents() );
	    hists_.find("minDR_l")->second->Fill( minDR_ );
	    hists_.find("eta_l")->second->Fill( eta_ );
	    hists_.find("phi_l")->second->Fill( phi_ );
	    hists_.find("emf_l")->second->Fill( emf_ );
	  }
	  if ( HLTriggerMatchedJet_pt >= 30 &&  HLTriggerMatchedJet_pt < 40 ) hists_.find("nConst_m1")->second->Fill( (*jets)[iJet].nConstituents() );
	  if ( HLTriggerMatchedJet_pt >= 40 &&  HLTriggerMatchedJet_pt < 50 ) hists_.find("nConst_m2")->second->Fill( (*jets)[iJet].nConstituents() );

	  if ( HLTriggerMatchedJet_pt >= 30 ) {
	    hists_.find("dR2Trig_h")->second->Fill( dR2Trig );
	    if ( HLTriggerMatchedJet_pt >=  50 ) hists_.find("nConst_h")->second->Fill( (*jets)[iJet].nConstituents() );
	    hists_.find("minDR_h")->second->Fill( minDR_ );
	    hists_.find("eta_h")->second->Fill( eta_ );
	    hists_.find("phi_h")->second->Fill( phi_ );
	    hists_.find("emf_h")->second->Fill( emf_ );
	  }
	  if     (jetCounter == 1) {
	    hists_.find("pt1_match")->second->Fill( HLTriggerMatchedJet_pt );
	    hists_.find("pt1_trigger_match")->second->Fill( HLTriggerJet_pt );
	    HLTriggerMatchedJet1_pt = HLTriggerMatchedJet_pt; 
	    if ( HLTriggerMatchedJet_pt < 30 ) hists_.find("dR2Trig_1")->second->Fill( dR2Trig );
	  }
	  else if(jetCounter == 2) { 
	    hists_.find("pt2_match")->second->Fill( HLTriggerMatchedJet_pt ); 
	    hists_.find("pt2_trigger_match")->second->Fill( HLTriggerJet_pt ); 
	    HLTriggerMatchedJet2_pt = HLTriggerMatchedJet_pt; 
	    if ( HLTriggerMatchedJet1_pt > 40 ) hists_.find("depPt2")->second->Fill( HLTriggerMatchedJet_pt );
	    if ( HLTriggerMatchedJet_pt < 30 ) hists_.find("dR2Trig_2")->second->Fill( dR2Trig );
	  }
	  else if(jetCounter == 3) { 
	    hists_.find("pt3_match")->second->Fill( HLTriggerMatchedJet_pt ); 
	    hists_.find("pt3_trigger_match")->second->Fill( HLTriggerJet_pt ); 
	    HLTriggerMatchedJet3_pt = HLTriggerMatchedJet_pt; 
	    if ( HLTriggerMatchedJet2_pt > 40 ) hists_.find("depPt3")->second->Fill( HLTriggerMatchedJet_pt );
	    if ( HLTriggerMatchedJet_pt < 30 ) hists_.find("dR2Trig_3")->second->Fill( dR2Trig );
	  }
	  else if(jetCounter == 4) { 
	    hists_.find("pt4_match")->second->Fill( HLTriggerMatchedJet_pt ); 
	    hists_.find("pt4_trigger_match")->second->Fill( HLTriggerJet_pt ); 
	    if ( HLTriggerMatchedJet3_pt > 40 ) hists_.find("depPt4")->second->Fill( HLTriggerMatchedJet_pt );
	    if ( HLTriggerMatchedJet_pt < 30 ) hists_.find("dR2Trig_4")->second->Fill( dR2Trig );
	  }
	  else if(jetCounter == 5) { 
	    hists_.find("pt5_match")->second->Fill( HLTriggerMatchedJet_pt ); 
	    hists_.find("pt5_trigger_match")->second->Fill( HLTriggerJet_pt ); 
	  }
	  else if(jetCounter == 6) { 
	    hists_.find("pt6_match")->second->Fill( HLTriggerMatchedJet_pt ); 
	    hists_.find("pt6_trigger_match")->second->Fill( HLTriggerJet_pt ); 
	  }
	}
      }
      
      if(singleJetTriggered){
	
	pat::TriggerObjectRefVector objects = triggerEvent->filterObjects("hlt1jet30");
	std::vector<double> objects_pt;
	
	for(unsigned int i = 0; i < objects.size(); ++i){
	  objects_pt.push_back(rnd.Gaus(objects[i]->pt(), width_));
	}

	if(objects.size() >= 2){
	  if(objects[0]->pt() >= 40. && objects[1]->pt() >= 40.) {
	    
	    //effy = triggerEfficiency(jets, true);
	    effy = triggerEfficiencySmear(objects_pt, true);
	    
	    if(jets->size() >= 1) {
	      hists_.find("DJ40_1")->second->Fill(   objects_pt[0], effy );
	      hists_.find("DJ40_1_n")->second->Fill( objects_pt[0] );
	    }
	    if(jets->size() >= 2) {
	      hists_.find("DJ40_2")->second->Fill(   objects_pt[1], effy );
	      hists_.find("DJ40_2_n")->second->Fill( objects_pt[1] );
	    }
	    
	    hists_.find("DJ40_trigger_1")->second->Fill( objects[0]->pt(), effy ); hists_.find("DJ40_trigger_1_n")->second->Fill( objects[0]->pt() );
	    hists_.find("DJ40_trigger_2")->second->Fill( objects[1]->pt(), effy ); hists_.find("DJ40_trigger_2_n")->second->Fill( objects[1]->pt() );
	    
	  }
	}
      }
      
      bool QJ40Triggered = false;

      if(HLTriggered){
	
	pat::TriggerObjectRefVector objects = triggerEvent->filterObjects("hlt4jet30");
	std::vector<double> objects_pt;
	
	for(unsigned int i = 0; i < objects.size(); ++i){
	  objects_pt.push_back(rnd.Gaus(objects[i]->pt(), width_));
	}

	hists_.find("n_trigger")->second->Fill( objects_pt.size() );
	unsigned int HLTriggerJetCounter = 1;
	for (size_t iTrigJet = 0; iTrigJet < objects_pt.size(); ++iTrigJet, ++HLTriggerJetCounter) {
	  double HLTriggerJet_pt = objects_pt[iTrigJet];
	  hists_.find("pt_trigger")->second->Fill( HLTriggerJet_pt );
	  if(HLTriggerJetCounter == 1) { hists_.find("pt1_trigger")->second->Fill( HLTriggerJet_pt ); HLTriggerJet1_pt = HLTriggerJet_pt; }
	  if(HLTriggerJetCounter == 2) { hists_.find("pt2_trigger")->second->Fill( HLTriggerJet_pt ); HLTriggerJet2_pt = HLTriggerJet_pt; }
	  if(HLTriggerJetCounter == 3) { hists_.find("pt3_trigger")->second->Fill( HLTriggerJet_pt ); HLTriggerJet3_pt = HLTriggerJet_pt; }
	  if(HLTriggerJetCounter == 4) { hists_.find("pt4_trigger")->second->Fill( HLTriggerJet_pt ); HLTriggerJet4_pt = HLTriggerJet_pt; }
	  if(HLTriggerJetCounter == 5) { hists_.find("pt5_trigger")->second->Fill( HLTriggerJet_pt ); } //HLTriggerJet5_pt = HLTriggerJet_pt; }
	  if(HLTriggerJetCounter == 6)   hists_.find("pt6_trigger")->second->Fill( HLTriggerJet_pt );  
	}

	if(objects[0]->pt() >= 30. && objects[1]->pt() >= 30. && objects[2]->pt() >= 30. && objects[3]->pt() >= 30.) {
	  
	  //effy = triggerEfficiency(jets, false);
	  effy = triggerEfficiencySmear(objects_pt, true);
	  
	  if(jets->size() >= 1) {
	    hists_.find("QJ30_1")->second->Fill( objects_pt[0], effy );
	    hists_.find("QJ30_1_n")->second->Fill( objects_pt[0] );
	    hists_.find("QJ30_1_trig")->second->Fill( objects_pt[0] );
	  }
	  if(jets->size() >= 2) {
	    hists_.find("QJ30_2")->second->Fill( objects_pt[1], effy );
	    hists_.find("QJ30_2_n")->second->Fill( objects_pt[1] );
	    hists_.find("QJ30_2_trig")->second->Fill( objects_pt[1] );
	  }
	  if(jets->size() >= 3) {
	    hists_.find("QJ30_3")->second->Fill( objects_pt[2], effy );
	    hists_.find("QJ30_3_n")->second->Fill( objects_pt[2] );
	    hists_.find("QJ30_3_trig")->second->Fill( objects_pt[2] );
	  }
	  if(jets->size() >= 4) {
	    hists_.find("QJ30_4")->second->Fill( objects_pt[3], effy );
	    hists_.find("QJ30_4_n")->second->Fill( objects_pt[3] );
	    hists_.find("QJ30_4_trig")->second->Fill( objects_pt[3] );
	  }

	  hists_.find("QJ30_trigger_1")->second->Fill( HLTriggerJet1_pt, effy ); hists_.find("QJ30_trigger_1_n")->second->Fill( HLTriggerJet1_pt );
	  hists_.find("QJ30_trigger_2")->second->Fill( HLTriggerJet2_pt, effy ); hists_.find("QJ30_trigger_2_n")->second->Fill( HLTriggerJet2_pt );
	  hists_.find("QJ30_trigger_3")->second->Fill( HLTriggerJet3_pt, effy ); hists_.find("QJ30_trigger_3_n")->second->Fill( HLTriggerJet3_pt );
	  hists_.find("QJ30_trigger_4")->second->Fill( HLTriggerJet4_pt, effy ); hists_.find("QJ30_trigger_4_n")->second->Fill( HLTriggerJet4_pt );
	  /*
	    if(objects[0]->pt() >= 80.) {
	    if(jets->size() >= 1) { hists_.find("QJ8333_1")->second->Fill( jet1_pt, effy ); hists_.find("QJ8333_1_n")->second->Fill( jet1_pt ); }
	    if(jets->size() >= 2) { hists_.find("QJ8333_2")->second->Fill( jet2_pt, effy ); hists_.find("QJ8333_2_n")->second->Fill( jet2_pt ); }
	    if(jets->size() >= 3) { hists_.find("QJ8333_3")->second->Fill( jet3_pt, effy ); hists_.find("QJ8333_3_n")->second->Fill( jet3_pt ); }
	    if(jets->size() >= 4) { hists_.find("QJ8333_4")->second->Fill( jet4_pt, effy ); hists_.find("QJ8333_4_n")->second->Fill( jet4_pt ); }
	    
	    hists_.find("QJ8333_trigger_1")->second->Fill( HLTriggerJet1_pt, effy ); hists_.find("QJ8333_trigger_1_n")->second->Fill( HLTriggerJet1_pt );
	    hists_.find("QJ8333_trigger_2")->second->Fill( HLTriggerJet2_pt, effy ); hists_.find("QJ8333_trigger_2_n")->second->Fill( HLTriggerJet2_pt );
	    hists_.find("QJ8333_trigger_3")->second->Fill( HLTriggerJet3_pt, effy ); hists_.find("QJ8333_trigger_3_n")->second->Fill( HLTriggerJet3_pt );
	    hists_.find("QJ8333_trigger_4")->second->Fill( HLTriggerJet4_pt, effy ); hists_.find("QJ8333_trigger_4_n")->second->Fill( HLTriggerJet4_pt );
	    }
	  */
	  if(objects[0]->pt() >= 40. && objects[1]->pt() >= 40. && objects[2]->pt() >= 40. && objects[3]->pt() >= 40.) {

	    QJ40Triggered = true;
	    if(jets->size() >= 1) {
	      hists_.find("QJ40_1")->second->Fill( objects_pt[0], effy );
	      hists_.find("QJ40_1_n")->second->Fill( objects_pt[0] );
	      hists_.find("QJ40_1_trig")->second->Fill( objects_pt[0] );
	    }
	    if(jets->size() >= 2) {
	      hists_.find("QJ40_2")->second->Fill( objects_pt[1], effy );
	      hists_.find("QJ40_2_n")->second->Fill( objects_pt[1] );
	      hists_.find("QJ40_2_trig")->second->Fill( objects_pt[1] );
	    }
	    if(jets->size() >= 3) {
	      hists_.find("QJ40_3")->second->Fill( objects_pt[2], effy );
	      hists_.find("QJ40_3_n")->second->Fill( objects_pt[2] );
	      hists_.find("QJ40_3_trig")->second->Fill( objects_pt[2] );
	    }
	    if(jets->size() >= 4) {
	      hists_.find("QJ40_4")->second->Fill( objects_pt[3], effy );
	      hists_.find("QJ40_4_n")->second->Fill( objects_pt[3] );
	      hists_.find("QJ40_4_trig")->second->Fill( objects_pt[3] );
	    }

	    hists_.find("QJ40_trigger_1")->second->Fill( HLTriggerJet1_pt, effy ); hists_.find("QJ40_trigger_1_n")->second->Fill( HLTriggerJet1_pt );
	    hists_.find("QJ40_trigger_2")->second->Fill( HLTriggerJet2_pt, effy ); hists_.find("QJ40_trigger_2_n")->second->Fill( HLTriggerJet2_pt );
	    hists_.find("QJ40_trigger_3")->second->Fill( HLTriggerJet3_pt, effy ); hists_.find("QJ40_trigger_3_n")->second->Fill( HLTriggerJet3_pt );
	    hists_.find("QJ40_trigger_4")->second->Fill( HLTriggerJet4_pt, effy ); hists_.find("QJ40_trigger_4_n")->second->Fill( HLTriggerJet4_pt );

	    double eta, emf, minDR_;

	    if(jets->size() >= 1 && jet1_pt < 40) { 
	      eta = (*jets)[0].eta();
	      emf = (*jets)[0].emEnergyFraction();
	      minDR_ = minDR(jets, 0);
	      hists_.find("eta_1")->second->Fill(eta);
	      hists_.find("emf_1")->second->Fill(emf); 
	      hists_.find("minDR_1")->second->Fill(minDR_); 
	      hists_.find("eta")->second->Fill(eta);
	      hists_.find("emf")->second->Fill(emf); 
	      hists_.find("minDR")->second->Fill(minDR_); 
	    }
	    if(jets->size() >= 2 && jet2_pt < 40) { 
	      eta = (*jets)[1].eta();
	      emf = (*jets)[1].emEnergyFraction();
	      minDR_ = minDR(jets, 1);
	      hists_.find("eta_2")->second->Fill(eta); 
	      hists_.find("emf_2")->second->Fill(emf); 
	      hists_.find("minDR_2")->second->Fill(minDR_); 
	      hists_.find("eta")->second->Fill(eta); 
	      hists_.find("emf")->second->Fill(emf); 
	      hists_.find("minDR")->second->Fill(minDR_); 
	    }
	    if(jets->size() >= 3 && jet3_pt < 40) { 
	      eta = (*jets)[2].eta();
	      emf = (*jets)[2].emEnergyFraction();
	      minDR_ = minDR(jets, 2);
	      hists_.find("eta_3")->second->Fill(eta); 
	      hists_.find("emf_3")->second->Fill(emf); 
	      hists_.find("minDR_3")->second->Fill(minDR_); 
	      hists_.find("eta")->second->Fill(eta); 
	      hists_.find("emf")->second->Fill(emf); 
	      hists_.find("minDR")->second->Fill(minDR_); 
	    }
	    if(jets->size() >= 4 && jet4_pt < 40) { 
	      eta = (*jets)[3].eta();
	      emf = (*jets)[3].emEnergyFraction();
	      minDR_ = minDR(jets, 3);
	      hists_.find("eta_4")->second->Fill(eta); 
	      hists_.find("emf_4")->second->Fill(emf); 
	      hists_.find("minDR_4")->second->Fill(minDR_); 
	      hists_.find("eta")->second->Fill(eta); 
	      hists_.find("emf")->second->Fill(emf); 
	      hists_.find("minDR")->second->Fill(minDR_); 
	    }

	  }
	  /*
	    if(objects.size() >= 5){
	    if(objects[4]->pt() >= 30.){

	    effy = triggerEfficiency(jets, true);

	    if(jets->size() >= 1) { hists_.find("FJ30_1")->second->Fill( jet1_pt, effy ); hists_.find("FJ30_1_n")->second->Fill( jet1_pt ); }
	    if(jets->size() >= 2) { hists_.find("FJ30_2")->second->Fill( jet2_pt, effy ); hists_.find("FJ30_2_n")->second->Fill( jet2_pt ); }
	    if(jets->size() >= 3) { hists_.find("FJ30_3")->second->Fill( jet3_pt, effy ); hists_.find("FJ30_3_n")->second->Fill( jet3_pt ); }
	    if(jets->size() >= 4) { hists_.find("FJ30_4")->second->Fill( jet4_pt, effy ); hists_.find("FJ30_4_n")->second->Fill( jet4_pt ); }
	    if(jets->size() >= 5) { hists_.find("FJ30_5")->second->Fill( jet5_pt, effy ); hists_.find("FJ30_5_n")->second->Fill( jet5_pt ); }

	    hists_.find("FJ30_trigger_1")->second->Fill( HLTriggerJet1_pt, effy ); hists_.find("FJ30_trigger_1_n")->second->Fill( HLTriggerJet1_pt );
	    hists_.find("FJ30_trigger_2")->second->Fill( HLTriggerJet2_pt, effy ); hists_.find("FJ30_trigger_2_n")->second->Fill( HLTriggerJet2_pt );
	    hists_.find("FJ30_trigger_3")->second->Fill( HLTriggerJet3_pt, effy ); hists_.find("FJ30_trigger_3_n")->second->Fill( HLTriggerJet3_pt );
	    hists_.find("FJ30_trigger_4")->second->Fill( HLTriggerJet4_pt, effy ); hists_.find("FJ30_trigger_4_n")->second->Fill( HLTriggerJet4_pt );
	    hists_.find("FJ30_trigger_5")->second->Fill( HLTriggerJet5_pt, effy ); hists_.find("FJ30_trigger_5_n")->second->Fill( HLTriggerJet5_pt );
	    }
	    }
	  */
	}
      }
      else if(!HLTriggered){

	if(jets->size() >= 1) {
	  hists_.find("QJ30_1_untrig")->second->Fill( jet1_pt );
	}
	if(jets->size() >= 2) {
	  hists_.find("QJ30_2_untrig")->second->Fill( jet2_pt );
	}
	if(jets->size() >= 3) {
	  hists_.find("QJ30_3_untrig")->second->Fill( jet3_pt );
	}
	if(jets->size() >= 4) {
	  hists_.find("QJ30_4_untrig")->second->Fill( jet4_pt );
	}
	if(!QJ40Triggered){
	  if(jets->size() >= 1) {
	    hists_.find("QJ40_1_untrig")->second->Fill( jet1_pt );
	  }
	  if(jets->size() >= 2) {
	    hists_.find("QJ40_2_untrig")->second->Fill( jet2_pt );
	  }
	  if(jets->size() >= 3) {
	    hists_.find("QJ40_3_untrig")->second->Fill( jet3_pt );
	  }
	  if(jets->size() >= 4) {
	    hists_.find("QJ40_4_untrig")->second->Fill( jet4_pt );
	  }
	}
      }
  }
}
//}

void
SingleObjectJetTrigger::endJob()
{
  std::vector<std::string> hists2Norm;
  hists2Norm.push_back("pt");
  hists2Norm.push_back("pt1");
  hists2Norm.push_back("pt2");
  hists2Norm.push_back("pt3");
  hists2Norm.push_back("pt4");
  hists2Norm.push_back("pt5");
  hists2Norm.push_back("pt6");
  for( std::vector<std::string>::const_iterator hist = hists2Norm.begin(); hist != hists2Norm.end(); ++hist ){
    hists_.find(*hist + "_effy")->second->Divide(hists_.find(*hist + "_match")->second, hists_.find(*hist)->second);
  }

  hists2Norm.clear();
  hists2Norm.push_back("QJ30_1");
  hists2Norm.push_back("QJ40_1");
  hists2Norm.push_back("QJ30_2");
  hists2Norm.push_back("QJ40_2");
  hists2Norm.push_back("QJ30_3");
  hists2Norm.push_back("QJ40_3");
  hists2Norm.push_back("QJ30_4");
  hists2Norm.push_back("QJ40_4");
  hists2Norm.push_back("QJ30_trigger_1");
  hists2Norm.push_back("QJ40_trigger_1");
  hists2Norm.push_back("QJ30_trigger_2");
  hists2Norm.push_back("QJ40_trigger_2");
  hists2Norm.push_back("QJ30_trigger_3");
  hists2Norm.push_back("QJ40_trigger_3");
  hists2Norm.push_back("QJ30_trigger_4");
  hists2Norm.push_back("QJ40_trigger_4");
  hists2Norm.push_back("DJ40_1");
  hists2Norm.push_back("DJ40_2");
  hists2Norm.push_back("DJ40_trigger_1");
  hists2Norm.push_back("DJ40_trigger_2");

  //hists2Norm.push_back("QJ8333_1");
  //hists2Norm.push_back("QJ8333_2");
  //hists2Norm.push_back("QJ8333_3");
  //hists2Norm.push_back("QJ8333_4");
  //hists2Norm.push_back("QJ8333_trigger_1");
  //hists2Norm.push_back("QJ8333_trigger_2");
  //hists2Norm.push_back("QJ8333_trigger_3");
  //hists2Norm.push_back("QJ8333_trigger_4");
  //hists2Norm.push_back("FJ30_1");
  //hists2Norm.push_back("FJ30_2");
  //hists2Norm.push_back("FJ30_3");
  //hists2Norm.push_back("FJ30_4");
  //hists2Norm.push_back("FJ30_5");
  //hists2Norm.push_back("FJ30_trigger_1");
  //hists2Norm.push_back("FJ30_trigger_2");
  //hists2Norm.push_back("FJ30_trigger_3");
  //hists2Norm.push_back("FJ30_trigger_4");
  //hists2Norm.push_back("FJ30_trigger_5");
  for( std::vector<std::string>::const_iterator hist = hists2Norm.begin(); hist != hists2Norm.end(); ++hist ){
    hists_.find(*hist + "_e")->second->Divide(hists_.find(*hist)->second, hists_.find(*hist + "_n")->second);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( SingleObjectJetTrigger );
