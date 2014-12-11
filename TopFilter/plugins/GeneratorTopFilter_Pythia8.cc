#include <memory>  

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "TopQuarkAnalysis/TopSkimming/interface/TtDecayChannelSelector.h"

/**
  \class   GeneratorTopFilter GeneratorTopFilter.cc "TopAnalysis/TopFilter/plugins/GeneratorTopFilter.h"

  \brief   EDFilter to select ttbar or single top decays on generator level
 
  filter for single top or ttbar decay channels on generator level.
  decay channels are internally coded by integer numbers:
  1 stands for a hadronically decaying top
  2 is top to electron
  3 muon
  4 tau further decaying to hadrons
  5 tau to electron
  6 tau to muon.
  in two digit numbers the first digit is the top and the second the antitop.
  so e.g. dimuon channel is 33 or single electron is 21 and 12
  
  The scheme listed above is valid in case that both tops decay to b quarks. For the rare cases where at least one of them decays to another quark,
  the numbering scheme is the same, but a third digit is added in front, in order to indicate these cases:
  1xx means that one of the two tops does not decay to b quarks
  2xx means that both tops do not decay to b quarks
*/



using namespace std;



class GeneratorTopFilter_Pythia8 : public edm::EDFilter {
 public:
  /// default constructor
  explicit GeneratorTopFilter_Pythia8(const edm::ParameterSet&);
  /// default destructor
  ~GeneratorTopFilter_Pythia8();
  
 private:
  /// read in which channels are going to be selected and print out
  virtual void beginJob();
  /// the real filtering procedure is implemendted in analyze(), wDecay() and tauDecay()
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  /// gives a summary how many events have been selected
  virtual void endJob();
   
  /// finds out if event passes or is blocked. Remark: if invert_selection is set as true in 
  /// config then the meaning of blocked and passed swapes
  bool analyze(edm::Event&, const edm::EventSetup&);
  /// gives back an index corresponding to the decay of the W  
  /// -1 in case of an error
  ///  1 for hadron
  ///  2 for electron
  ///  3 for muon
  ///  4 - 7 for tau (depending on output of tauDecay() )
  int wDecay(); 
  /// gives back an index corresponding to the decay of the tau 
  ///  4 for hadronic tau decay
  ///  5 for electron
  ///  6 for muon
  ///  7 for unknown
  int tauDecay();
  
  /// as source the gen particle collection shoud be used
  edm::InputTag src_;
  /// how many events do you want to be selected?
  int nEvents_;
  /// Should top decays where a top does not decay to a bottom quark being rejected?
  bool rejectNonBottomDecaysOfTops_;
  /// string to save channel short cut from config
  std::vector<std::string> channels_;
  /// counter for number of events
  int nEvts; 
  /// counter for filtered events
  int nHits;

  /// bool from config: select this channel?
  bool SingleTopHadronic_;
  /// bool from config: select this channel?
  bool SingleTopElectron_;
  /// bool from config: select this channel?
  bool SingleTopMuon_; 
  /// bool from config: select this channel?
  bool SingleTopTau_;
  /// bool from config: select this channel?
  bool TtbarHadronic_;
  /// bool from config: select this channel?
  bool TtbarSemiElec_;
  /// bool from config: select this channel?
  bool TtbarSemiMuon_; 
  /// bool from config: select this channel?
  bool TtbarSemiTau_;
  /// bool from config: select this channel?
  bool TtbarElecElec_;
  /// bool from config: select this channel?
  bool TtbarMuonMuon_;
  /// bool from config: select this channel?
  bool TtbarTauTau_;
  /// bool from config: select this channel?
  bool TtbarElecMuon_;
  /// bool from config: select this channel?
  bool TtbarElecTau_;    
  /// bool from config: select this channel?
  bool TtbarMuonTau_; 
  /// bool from config: select this tau decay?
  bool Tau2Had_;
  /// bool from config: select this tau decay?
  bool Tau2Elec_; 
  /// bool from config: select this tau decay?
  bool Tau2Muon_;
  /// bool from config: invert the selection?
  bool invert_;   
  
  /// vector to store the index numbers of the channels to be selected  
  vector<int> selectedChannels; 
  /// length of vector<int> selectedChannels;
  unsigned int selChSize; 
  
  /// pointer on bquark from top decay (due to CKM, almost all decay to this)
  const reco::Candidate* bquark;
  /// pointer on non-bquark from top decay (rare cases where top decays to d or s)
  const reco::Candidate* otherQuark;
  /// pointer on W boson
  const reco::Candidate* wbos;
  /// pointer on tau
  const reco::Candidate* tau;
  /// pointer on first daughter from W
  const reco::Candidate* wdaugh1;
  /// pointer on second daughter from W
  const reco::Candidate* wdaugh2; 
      
};



GeneratorTopFilter_Pythia8::GeneratorTopFilter_Pythia8(const edm::ParameterSet& cfg) :
  src_              (cfg.getParameter<edm::InputTag>("src"       )),
  nEvents_          (cfg.getParameter<int> ("n_events"           )),
  rejectNonBottomDecaysOfTops_(cfg.getParameter<bool>("rejectNonBottomDecaysOfTops")),
  channels_         (cfg.getParameter<std::vector<std::string> >("channels")),    
  SingleTopHadronic_(cfg.getParameter<bool>("SingleTop_Hadronic" )),  
  SingleTopElectron_(cfg.getParameter<bool>("SingleTop_Electron" )),  
  SingleTopMuon_    (cfg.getParameter<bool>("SingleTop_Muon"     )),  
  SingleTopTau_     (cfg.getParameter<bool>("SingleTop_Tau"      )),  
  TtbarHadronic_    (cfg.getParameter<bool>("Ttbar_Hadronic"     )), 
  TtbarSemiElec_    (cfg.getParameter<bool>("Ttbar_SemiElec"     )),   
  TtbarSemiMuon_    (cfg.getParameter<bool>("Ttbar_SemiMuon"     )),   
  TtbarSemiTau_     (cfg.getParameter<bool>("Ttbar_SemiTau"      )),   
  TtbarElecElec_    (cfg.getParameter<bool>("Ttbar_ElecElec"     )),   
  TtbarMuonMuon_    (cfg.getParameter<bool>("Ttbar_MuonMuon"     )),   
  TtbarTauTau_      (cfg.getParameter<bool>("Ttbar_TauTau"       )), 
  TtbarElecMuon_    (cfg.getParameter<bool>("Ttbar_ElecMuon"     )),   
  TtbarElecTau_     (cfg.getParameter<bool>("Ttbar_ElecTau"      )),   
  TtbarMuonTau_     (cfg.getParameter<bool>("Ttbar_MuonTau"      )),
  Tau2Had_          (cfg.getParameter<bool>("Tau_To_Hadron"      )),   
  Tau2Elec_         (cfg.getParameter<bool>("Tau_To_Electron"    )),
  Tau2Muon_         (cfg.getParameter<bool>("Tau_To_Muon"        )), 
  invert_           (cfg.getParameter<bool>("invert_selection"   ))              
{
    nEvts = 0;
    nHits = 0;
    produces<int>("decayMode");
}



GeneratorTopFilter_Pythia8::~GeneratorTopFilter_Pythia8() { 
}



bool GeneratorTopFilter_Pythia8::filter(edm::Event& evt, const edm::EventSetup& es) {
    nEvts++;
    bool result = analyze( evt, es);
    if(result) nHits++;
    return result;
}



bool GeneratorTopFilter_Pythia8::analyze(edm::Event& evt, const edm::EventSetup& es) {  
  if(nEvents_!=-1){
    if(nHits >= nEvents_){
      return false;      
    }
  }
  
  bool passed = true;
  bool blocked = false;
  
  if(invert_){
    blocked = true;
    passed  = false;    
  }
        
  int decayVl = 0;
  
  bquark = 0;
  otherQuark = 0;
  wbos = 0;
  tau = 0;
  wdaugh1 = 0;
  wdaugh2 = 0;


  edm::Handle<edm::View<reco::GenParticle>> genParticles;
  evt.getByLabel(src_, genParticles);
  for(edm::View<reco::GenParticle>::const_iterator cand = genParticles->begin(); cand!=genParticles->end(); ++cand) {

    if(abs(cand->pdgId())!=6) continue; 

    //    std::cout << "Particle ID: " << cand->pdgId()  << std::endl;
    //    std::cout << "Status     : " << cand->status() << std::endl;

if (abs(cand->pdgId())==6 && (cand->status()>69 || cand->status()<59)) continue;

    
    for(size_t i=0; i<cand->numberOfDaughters(); ++i){

      //      std::cout << "     Daughter PID  : " << cand->daughter(i)->pdgId()  << std::endl;
      //      std::cout << "     Dauther Status: " << cand->daughter(i)->status() << std::endl;

      const int absPdgId = std::abs(cand->daughter(i)->pdgId());
      /////FIXME: This should be implemented Recursively///////////////
      if(absPdgId==24){
         const reco::Candidate* wbos_temp = cand->daughter(i);
                for(size_t j=0; j<wbos_temp->numberOfDaughters(); ++j){
                //                std::cout << "          W-Daughter PID  : " << wbos_temp->daughter(j)->pdgId()  << std::endl;
                //                std::cout << "          W-Dauther Status: " << wbos_temp->daughter(j)->status() << std::endl;
                        if (std::abs(wbos_temp->daughter(j)->pdgId())==24 && wbos_temp->daughter(j)->status()== 52){
                                if (wbos_temp->daughter(j)->numberOfDaughters()>=2){
                                wbos = wbos_temp->daughter(j);
                                }
                                else wbos=wbos_temp->daughter(j)->daughter(0);
                        }
                        else if(wbos_temp->daughter(j)->status()==23 || wbos_temp->daughter(j)->status()==1){
                                wbos = wbos_temp;
                        }

                }
      }


      else if(absPdgId==5){
	bquark = cand->daughter(i);
      }
      else if(absPdgId==1 || absPdgId==3){
        otherQuark = cand->daughter(i);
      }
    }
    // In few events the top decays to s or d, but these are not incorrect events, thus they should get a top decay number
    // CAVEAT: if in the MC the W is not in the decay chain but decays hadronically, there will be a quark of type d or s:
    // these cases cannot be treated properly, since it is not clear which quark comes from the top decay and not from the W,
    // thus it will not get a top decay number
    if(bquark==0){ 
      if(wbos!=0 && otherQuark!=0){
        decayVl += 100;
      }
      else{
        edm::LogInfo("TopDecaySelector") << "Top does not decay to b, and also not to d or s quark.";
        return blocked;
      }
    }
    // in some MC the w is not in the decay chain but directly its decay products
    if(wbos==0){
      // search for top daughters which are neither the b nor a radiated gluon
      for(size_t i=0; i<cand->numberOfDaughters(); ++i){        
        if(abs(cand->daughter(i)->pdgId())!=5 && cand->daughter(i)->pdgId()!=21){
	  wdaugh1 = cand->daughter(i);
	  break;
	}
      }
      for(size_t i=0; i<cand->numberOfDaughters(); ++i){   	  
        if(abs(cand->daughter(i)->pdgId())!=5 && cand->daughter(i)->pdgId()!=21 && 
	 cand->daughter(i)!=wdaugh1){	     
	   wdaugh2 = cand->daughter(i);	  	  
	}	  	
      }			        
    }
    else{
      // if w is in decay chain search for its daughters
      if(wbos->numberOfDaughters()>=2){
        wdaugh1 = wbos->daughter(0);
	wdaugh2 = wbos->daughter(1);
      }   	      
    }
    if(wdaugh1==0 || wdaugh2==0){
      edm::LogWarning("TopDecaySelector") << "Top has no W daughter. Probable Bug. Should be investigated   "<<"W Status:  "<<wbos->status()<<"W Boson Daughter:  "<< wbos->daughter(0)->pdgId();
      return blocked;	
    }
    // calculate the decay value from the w decays
    int wNumber = wDecay();
    
    if(wNumber<0){ 
      edm::LogInfo("TopDecaySelector") << "Could not identify W-decay.";
      return blocked;
    }
    
    if(wNumber==7){
      // if it happens that the wNumber is 7 that means that no tau daughter
      // is in the genParticle collection. This then is an error of the MC.
      // If at least one tau decay mode is blocked it is assumed that only
      // special tau decays shoul be selected. So in this case "stable taus"
      // are blocked.
      edm::LogInfo("TopDecaySelector") << "Tau from top decay does not decay.";
      if(!Tau2Had_ || !Tau2Elec_ || ! Tau2Muon_) return blocked;
    }
    
    if(cand->pdgId()==6){ 
      decayVl = decayVl + 10*wNumber;
    }
    else{
      decayVl = decayVl + wNumber;	  
    }
  }

  // Store the decay mode
  std::auto_ptr<int> decay(new int);
  *decay = decayVl;
  evt.put(decay, "decayMode");
  
  // Do filtering based on selected channels, and on whether to allow non-bottom decays of top
  if(!rejectNonBottomDecaysOfTops_) decayVl = decayVl%100;
  for(unsigned int i=0; i<selChSize; i++){
    if(decayVl==selectedChannels.at(i)){
        return passed;
    }
  }
  return blocked;
}   



int GeneratorTopFilter_Pythia8::wDecay(){
  if((abs(wdaugh1->pdgId())<6)||(abs(wdaugh2->pdgId())<6)){	  
    return 1; 
  }
  if((abs(wdaugh1->pdgId())==11)||(abs(wdaugh2->pdgId())==11)){	  
    return 2;	    
  }
  if((abs(wdaugh1->pdgId())==13)||(abs(wdaugh2->pdgId())==13)){	  
    return 3;
  }
  if((abs(wdaugh1->pdgId())==15)||(abs(wdaugh2->pdgId())==15)){
    tau = wdaugh2;
    if(abs(wdaugh1->pdgId())==15){    	  
      tau = wdaugh1;
    }
    return tauDecay();          
  }

  // if w does not decay into quarks or leptons
  return -1;
}



int GeneratorTopFilter_Pythia8::tauDecay(){  
  for(unsigned int i=0; i<tau->numberOfDaughters();++i){
    if(abs(tau->daughter(i)->pdgId())>99 || abs(tau->daughter(i)->pdgId())<10){ 
      return 4;
    }
    if(abs(tau->daughter(i)->pdgId())==11){ 
      return 5;
    }
    if(abs(tau->daughter(i)->pdgId())==13){ 
      return 6;
    }
    if(abs(tau->daughter(i)->pdgId())==15||abs(tau->daughter(i)->pdgId())==24){
      tau = tau->daughter(i);
      return tauDecay();
    }  
  }
  // if tau is stable
  return 7;   
}



void GeneratorTopFilter_Pythia8::beginJob() {
  
  // check if channel short cuts are given in config
  if(channels_.size()>0){
    for(size_t i=0;i<channels_.size();++i){
      std::string chn = channels_[i];
      if(chn == "ElectronElectron"){ 
	selectedChannels.push_back(22);
      }	 
      else if(chn == "ElectronMuon"){ 
	selectedChannels.push_back(23); 
	selectedChannels.push_back(32); 	  
      }
      else if(chn == "MuonElectron"){ // this is not forseen but it is clear what is meant 
	  selectedChannels.push_back(23); 
	  selectedChannels.push_back(32); 	  
      }	  
      else if(chn == "MuonMuon"){ 
	  selectedChannels.push_back(33); 
      }
      else if(chn == "ElectronElectronViaTau"){  
	  selectedChannels.push_back(25); 
	  selectedChannels.push_back(52); 
	  selectedChannels.push_back(55);	  
      }	  	  
      else if(chn == "ElectronMuonViaTau"){ 
	  selectedChannels.push_back(26); 
	  selectedChannels.push_back(62); 
	  selectedChannels.push_back(35); 
	  selectedChannels.push_back(53); 
	  selectedChannels.push_back(56); 
	  selectedChannels.push_back(65); 	  
      }	       
      else if(chn == "MuonElectronViaTau"){ // this is not forseen but it is clear what is meant  
	  selectedChannels.push_back(26); 
	  selectedChannels.push_back(62); 
	  selectedChannels.push_back(35); 
	  selectedChannels.push_back(53); 
	  selectedChannels.push_back(56); 
	  selectedChannels.push_back(65); 	  
      }  
      else if(chn == "MuonMuonViaTau"){ 
	  selectedChannels.push_back(36); 
	  selectedChannels.push_back(63); 
	  selectedChannels.push_back(66); 	  
      }
      else if(chn == "SingleElectron"){ 
	  selectedChannels.push_back(12); 
	  selectedChannels.push_back(21); 	  
      }	  	  
      else if(chn == "SingleMuon"){ 
	  selectedChannels.push_back(13); 
	  selectedChannels.push_back(31); 	  
      }
      else if(chn == "FullHadronic"){ 
	  selectedChannels.push_back(11); 	  
      }
      else if(chn == "none"){
          //empty
      }
      else{
	edm::LogError("GeneratorTopFilter") << "Unknown channel short cut in configuration: " << chn;
	throw cms::Exception("Configuration Error");	  	      
      }
    }
  }

  else{
    // build vector of selected decay channel numbers
    if(SingleTopHadronic_){
      selectedChannels.push_back(1);  
      selectedChannels.push_back(10);
    }

    if(SingleTopElectron_){
      selectedChannels.push_back(2);  
      selectedChannels.push_back(20);
    }

    if(SingleTopMuon_){
      selectedChannels.push_back(3);  
      selectedChannels.push_back(30);
    } 

    if(SingleTopTau_){
      if(Tau2Had_){
	selectedChannels.push_back( 4);  
	selectedChannels.push_back(40);
      }
      if(Tau2Elec_){      
	selectedChannels.push_back( 5);  
	selectedChannels.push_back(50); 
      }
      if(Tau2Muon_){                 
	selectedChannels.push_back( 6);  
	selectedChannels.push_back(60);
      }       
    }

    if(TtbarHadronic_){
      selectedChannels.push_back(11);
    }

    if(TtbarSemiElec_){
      selectedChannels.push_back(12); 
      selectedChannels.push_back(21);
    }

    if(TtbarSemiMuon_){
      selectedChannels.push_back(13); 
      selectedChannels.push_back(31);
    } 

    if(TtbarSemiTau_){
      if(Tau2Had_){  
	selectedChannels.push_back(14); 
	selectedChannels.push_back(41);
      }
      if(Tau2Elec_){     
	selectedChannels.push_back(15); 
	selectedChannels.push_back(51);  
      }
      if(Tau2Muon_){        
	selectedChannels.push_back(16); 
	selectedChannels.push_back(61); 
      }     
    }

    if(TtbarElecElec_){
      selectedChannels.push_back(22);
    }

    if(TtbarMuonMuon_){
      selectedChannels.push_back(33);
    }

    if(TtbarTauTau_){
      if(Tau2Had_){   
	selectedChannels.push_back(44);
      }
      if(Tau2Elec_){     
	selectedChannels.push_back(55); 
      }
      if(Tau2Muon_){         
	selectedChannels.push_back(66);
      }
      if(Tau2Had_ && Tau2Elec_){           
	selectedChannels.push_back(45);
	selectedChannels.push_back(54);
      } 
      if(Tau2Had_ && Tau2Muon_){   
	selectedChannels.push_back(46);    
	selectedChannels.push_back(64);
      } 
      if(Tau2Elec_ && Tau2Muon_){          
	selectedChannels.push_back(56);    
	selectedChannels.push_back(65);
      }          
    }

    if(TtbarElecMuon_){
      selectedChannels.push_back(23); 
      selectedChannels.push_back(32);
    }

    if(TtbarElecTau_){
      if(Tau2Had_){  
	selectedChannels.push_back(24); 
	selectedChannels.push_back(42);
      }
      if(Tau2Elec_){     
	selectedChannels.push_back(25); 
	selectedChannels.push_back(52);  
      }
      if(Tau2Muon_){        
	selectedChannels.push_back(26); 
	selectedChannels.push_back(62); 
      }    
    }  

    if(TtbarMuonTau_){
      if(Tau2Had_){  
	selectedChannels.push_back(34); 
	selectedChannels.push_back(43);
      }
      if(Tau2Elec_){     
	selectedChannels.push_back(35); 
	selectedChannels.push_back(53);  
      }
      if(Tau2Muon_){        
	selectedChannels.push_back(36); 
	selectedChannels.push_back(63); 
      }     
    } 
  }
  
  selChSize = selectedChannels.size();
  
  edm::LogInfo log("GeneratorTopFilter");
  
  bool saidTau     = false;
  bool saidTauTau  = false;
  bool saidTauMuon = false;
  bool saidTauElec = false;
  bool saidTauHadr = false;
  
  log << "==============================================================\n";
  if(invert_==false) log << "List of selected decay channels:\n";
  else log << "All decay channels selected except of:\n";
  for(unsigned int i=0; i<selChSize; i++){
    int channel = selectedChannels[i];
    int first = channel/10;
    int second = channel%10;
    
    if(first == 0 || second==0){
      if(first==1) log << "single top to hadrons\n";
      else if(first==2) log << "single top to electron\n";
      else if(first==3) log << "single top to muon\n";
      else if(first>3 && !saidTau){ 
        log << "single top to tau\n";
        saidTau = true;
      }
    }
  
    else if(channel==11) log << "fullhadronic\n";
  
    else if((first == second) || (first >= 4 && second >= 4)){
      if(channel==22) log << "dileptonic e e\n";
      else if(channel==33) log << "dileptonic mu mu\n";
      else if(!saidTauTau){ 
        log << "dileptonic tau tau\n";
	saidTauTau = true;
      }
    }
    
    else if(first > second){
      if(second>1){
	if(first==3) log << "dilepton mu e\n";
	else if(first>=4 && second==3 && !saidTauMuon){
	  log << "dilepton tau mu\n";
	  saidTauMuon = true;
	}
	else if(second==2 && !saidTauElec){
	  log << "tau e\n";
	  saidTauElec = true; 
	}
      }
      else{
        log << "semileptonic with ";
	if(first==2) log << "electron\n";
	else if(first==3) log << "muon\n";
	else if(!saidTauHadr){  
	  log << "tau\n";
	  saidTauHadr = true;
	}
      }
    }      
  } 
  if(!Tau2Had_ || !Tau2Elec_ || !Tau2Muon_ ){
    log << "\nTaus decaying to ";
    if(!Tau2Had_ ) log << "hadrons ";
    if(!Tau2Elec_) log << "electrons ";
    if(!Tau2Muon_)  log << "muons ";
    log << "are excluded form selection.\n";
  }     
  log << "==============================================================\n";   
}



void GeneratorTopFilter_Pythia8::endJob() {  
  edm::LogInfo log("GeneratorTopFilter");
  log << "==============================================================\n";
  log << "       Selected: " << nHits  << " of " << nEvts << "\n";
  log << "==============================================================\n";   
}



#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GeneratorTopFilter_Pythia8);
