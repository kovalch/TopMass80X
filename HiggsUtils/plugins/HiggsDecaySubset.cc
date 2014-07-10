// -*- C++ -*-
//
// Package:    HiggsDecaySubset
// Class:      HiggsDecaySubset
// 
/**\class HiggsDecaySubset HiggsDecaySubset.cc TopAnalysis/HiggsUtils/plugins/HiggsDecaySubset.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Johannes Hauk,,,DESY
//         Created:  Tue Jan 15 14:35:52 CET 2013
// $Id: HiggsDecaySubset.cc,v 1.5 2013/02/19 17:44:01 hauk Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TopQuarkAnalysis/TopEventProducers/interface/TopDecaySubset.h"

//
// class declaration
//

class HiggsDecaySubset : public edm::EDProducer {
   public:
      explicit HiggsDecaySubset(const edm::ParameterSet&);
      ~HiggsDecaySubset();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      
      /// find Higgs in list of input particles
      std::vector<const reco::GenParticle*> findHiggs(const reco::GenParticleCollection&);
      
      /// check the decay chain for the exploited shower model
      TopDecaySubset::ShowerModel checkShowerModel(const std::vector<const reco::GenParticle*>&)const;
      
      /// clear references
      void clearReferences();
      
      /// fill output vector for full decay chain 
      void fillListing(const std::vector<const reco::GenParticle*>&, reco::GenParticleCollection&);
      
      /// calculate lorentz vector from input
      reco::Particle::LorentzVector p4(const reco::GenParticle::const_iterator, int);
      
      /// fill vector including all radiations from quarks
      void addRadiation(int&, const reco::GenParticle::const_iterator, reco::GenParticleCollection&);
      
      /// fill references for output vector
      void fillReferences(const reco::GenParticleRefProd&, reco::GenParticleCollection&);
      
      /// recursively fill vector for all further decay particles of a given particle
      void addDaughters(int& idx, const reco::GenParticle::const_iterator part, reco::GenParticleCollection& target, bool recursive=true);
      
      // ----------member data ---------------------------
      
      /// input tag for the genParticle source
      edm::InputTag genParticleSource_;
      
      /// add radiation or not?
      bool addRadiation_;
      
      /// parton shower mode (filled in checkShowerModel)
      TopDecaySubset::ShowerModel showerModel_;
      
      /// management of daughter indices for fillRefs
      std::map<int,std::vector<int> > refs_;
      /// index in new evt listing of parts with daughters; 
      /// has to be set to -1 in produce to deliver consistent 
      /// results!
      int motherPartIdx_;                    
      
      /// mode of decaySubset creation 
      TopDecaySubset::FillMode fillMode_;
      
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
HiggsDecaySubset::HiggsDecaySubset(const edm::ParameterSet& iConfig): 
genParticleSource_(iConfig.getParameter<edm::InputTag>("src")),
addRadiation_(iConfig.getParameter<bool>("addRadiation")),
showerModel_(TopDecaySubset::kStart)
{
    // Mapping of the corresponding fillMode;
    // see FillMode enumerator of TopDecaySubset for available modes
    const std::string mode = iConfig.getParameter<std::string>("fillMode");
    if(mode=="kME")
        fillMode_ = TopDecaySubset::kME;
    else if(mode=="kStable")
        fillMode_ = TopDecaySubset::kStable;
    else
        throw cms::Exception("Configuration")<<mode<<" is not a supported FillMode!\n";
    
    // Produces a set of GenParticles following
    // the decay branch of Higgs
    produces<reco::GenParticleCollection>(); 
}


HiggsDecaySubset::~HiggsDecaySubset()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
HiggsDecaySubset::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // Create target vector
    std::auto_ptr<reco::GenParticleCollection> target(new reco::GenParticleCollection);
    
    // Get source collection
    edm::Handle<reco::GenParticleCollection> src;
    iEvent.getByLabel(genParticleSource_, src);
    
    // Find Higgs in list of input particles
    const std::vector<const reco::GenParticle*> v_higgs = this->findHiggs(*src);
    
    // Determine shower model (only in first event)
    if(showerModel_==TopDecaySubset::kStart) showerModel_ = this->checkShowerModel(v_higgs);
    
    // Find decay chain
    if(showerModel_!=TopDecaySubset::kNone){
        // get ref product from the event
        const reco::GenParticleRefProd ref = iEvent.getRefBeforePut<reco::GenParticleCollection>();
        // clear existing refs
        this->clearReferences();
        // fill output
        this->fillListing(v_higgs, *target);
        // fill references
        this->fillReferences(ref, *target);
    }
    
    // Write vectors to the event
    iEvent.put(target);
}



// Find Higgs in list of input particles
std::vector<const reco::GenParticle*>
HiggsDecaySubset::findHiggs(const reco::GenParticleCollection& v_particle)
{
    std::vector<const reco::GenParticle*> v_higgs;
    for(reco::GenParticleCollection::const_iterator i_particle=v_particle.begin(); i_particle!=v_particle.end(); ++i_particle){
        if(i_particle->pdgId()==25 && i_particle->status()==3){
            //std::cout<<"\t"<<i_particle->pdgId()<<" "<<i_particle->status()<<" "<<i_particle->numberOfMothers()<<" "<<i_particle->numberOfDaughters()<<std::endl;
            v_higgs.push_back( &(*i_particle) );
        }
    }
    if(v_higgs.size() != 1) edm::LogWarning("decayChain")<<"Not exactly 1 Higgs found in event, but: "<<v_higgs.size();
    return v_higgs;
}



// check the decay chain for the exploited shower model
TopDecaySubset::ShowerModel
HiggsDecaySubset::checkShowerModel(const std::vector<const reco::GenParticle*>& v_higgs)const
{
    for(std::vector<const reco::GenParticle*>::const_iterator i_higgs = v_higgs.begin(); i_higgs != v_higgs.end(); ++i_higgs){
        const reco::GenParticle* higgs = *i_higgs;
        // check for kHerwig type showers: here the status 3 Higgs will 
        // have a single status 2 Higgs as daughter, which has again 3 
        // or more status 2 daughters
        // WARNING: I did not explicitely check for Higgs in HERWIG, only in PYTHIA, so HERWIG criterion can be wrong
        if(higgs->numberOfDaughters() == 1){
            if(higgs->begin()->pdgId()==higgs->pdgId() && higgs->begin()->status()==2 && higgs->begin()->numberOfDaughters()>1)
                return TopDecaySubset::kHerwig;
        }
        // Check for kPythia type showers: here the status 3 Higgs will 
        // have all decay products and a status 2 Higgs as daughters
        // the status 2 Higgs will be w/o further daughters
        if(higgs->numberOfDaughters() > 1){
            bool daughterHiggsHasNoDaughters(false);
            for(reco::GenParticle::const_iterator i_higgsDaughter = higgs->begin(); i_higgsDaughter != higgs->end(); ++i_higgsDaughter){
                //std::cout<<"\t"<<i_higgsDaughter->pdgId()<<" "<<i_higgsDaughter->status()<<" "<<i_higgsDaughter->numberOfDaughters()<<std::endl;
                if(i_higgsDaughter->pdgId()==25 && i_higgsDaughter->status()==2 && i_higgsDaughter->numberOfDaughters()==0) daughterHiggsHasNoDaughters = true;
            }
            if(daughterHiggsHasNoDaughters) return TopDecaySubset::kPythia;
        }
    }
    // If neither Herwig nor Pythia like
    if(v_higgs.size() == 0)
        edm::LogWarning("decayChain")
            <<" Failed to find Higgs in decay chain. Will assume that this a \n"
            <<" non-Higgs sample and produce an empty decaySubset.\n";
    else
        throw edm::Exception(edm::errors::LogicError,
            " Can not find back any of the supported hadronization models. Models \n"
            " which are supported are:                                            \n"
            " Pythia  LO(+PS): Higgs(status 3) --> particle(status 3) antiparticle(status 3)\n"
            " Herwig NLO(+PS): Higgs(status 2) --> Higgs(status 3) --> Higgs(status 2)  \n");
    return TopDecaySubset::kNone;
}



// Clear references
void
HiggsDecaySubset::clearReferences()
{
    // Clear vector of references 
    refs_.clear();  
    // Set idx for mother particles to a start value
    // of -1 (the first entry will raise it to 0)
    motherPartIdx_ = -1;
}



void 
HiggsDecaySubset::fillListing(const std::vector<const reco::GenParticle*>& v_higgs, reco::GenParticleCollection& target)
{
    // Determine status flag of the new 
    // particle depending on the FillMode
    unsigned int statusFlag;
    fillMode_ == TopDecaySubset::kME ? statusFlag=3 : statusFlag=2;
    
    for(std::vector<const reco::GenParticle*>::const_iterator i_higgs = v_higgs.begin(); i_higgs != v_higgs.end(); ++i_higgs){
        const reco::GenParticle* higgs = *i_higgs;
        
        // Do we need a case distiction for additional radiation off the Higgs? I do not think so, so p4(i_higgs, statusFlag) is not implemented, and no further particles are kept in case of addRadiation_
        //std::auto_ptr<reco::GenParticle> higgsPtr(new reco::GenParticle(higgs->threeCharge(), this->p4(i_higgs, statusFlag), higgs->vertex(), higgs->pdgId(), statusFlag, false));
        std::auto_ptr<reco::GenParticle> higgsPtr(new reco::GenParticle(higgs->threeCharge(), (*i_higgs)->p4(), higgs->vertex(), higgs->pdgId(), statusFlag, false));
        target.push_back(*higgsPtr);
        ++motherPartIdx_;
        
        // Keep the higgs index for the map to manage the daughter refs
        const int iHiggs = motherPartIdx_; 
        std::vector<int> higgsDaughters;
        
        // Sanity check
        if(showerModel_!=TopDecaySubset::kPythia && higgs->begin()==higgs->end())
            throw edm::Exception(edm::errors::LogicError,
                "showerModel_!=kPythia && higgs->begin()==higgs->end()\n");
        
        // Iterate over Higgs daughters
        for(reco::GenParticle::const_iterator i_higgsDaughter = ((showerModel_==TopDecaySubset::kPythia)?higgs->begin():higgs->begin()->begin()); i_higgsDaughter != ((showerModel_==TopDecaySubset::kPythia)?higgs->end():higgs->begin()->end()); ++i_higgsDaughter){
            //std::cout<<"\t"<<i_higgsDaughter->pdgId()<<" "<<i_higgsDaughter->status()<<" "<<i_higgsDaughter->numberOfDaughters()<<std::endl;
            
            // If particle is beauty or other lighter quark
            if(i_higgsDaughter->status()==3 && std::abs(i_higgsDaughter->pdgId())<=5){ 
                std::auto_ptr<reco::GenParticle> bPtr(new reco::GenParticle(i_higgsDaughter->threeCharge(), this->p4(i_higgsDaughter, statusFlag), i_higgsDaughter->vertex(), i_higgsDaughter->pdgId(), statusFlag, false));
                target.push_back(*bPtr);
                // increment & push index of the higgs daughter
                higgsDaughters.push_back(++motherPartIdx_); 
                if(addRadiation_){
                    this->addRadiation(motherPartIdx_, i_higgsDaughter, target); 
                }
            }
            
            // Sanity check
            if(showerModel_!=TopDecaySubset::kPythia && i_higgsDaughter->begin()==i_higgsDaughter->end())
                throw edm::Exception(edm::errors::LogicError,
                    "showerModel_!=kPythia && i_higgsDaughter->begin()==i_higgsDaughter->end()\n");
            
            reco::GenParticle::const_iterator buffer = showerModel_==TopDecaySubset::kPythia ? i_higgsDaughter : i_higgsDaughter->begin();
            const int status(buffer->status());
            const int absPdgId(std::abs(buffer->pdgId()));
            
            // If particle is a muon
            if(status==3 && absPdgId==13){
                std::auto_ptr<reco::GenParticle> muonPtr(new reco::GenParticle(buffer->threeCharge(), p4(buffer, statusFlag), buffer->vertex(), buffer->pdgId(), statusFlag, true));
                target.push_back(*muonPtr);
                // Increment & push index of the top daughter
                higgsDaughters.push_back(++motherPartIdx_);
                if(addRadiation_){
                    this->addRadiation(motherPartIdx_, buffer, target);
                }
            }
            
            // If particle is a tau lepton
            else if(status==3 && absPdgId==15){
                std::auto_ptr<reco::GenParticle> tauPtr(new reco::GenParticle(buffer->threeCharge(), p4(buffer, statusFlag), buffer->vertex(), buffer->pdgId(), statusFlag, true));
                target.push_back(*tauPtr);
                // Increment & push index of the top daughter
                higgsDaughters.push_back(++motherPartIdx_);
                if(addRadiation_){
                    this->addRadiation(motherPartIdx_, buffer, target);
                }
                // If further decay chain of tau becomes interesting, add it here (would need tau lepton index to manage daughter refs)
            }
            
            // If particle is a gluon
            else if(status==3 && absPdgId==21){
                std::auto_ptr<reco::GenParticle> gluonPtr(new reco::GenParticle(buffer->threeCharge(), p4(buffer, statusFlag), buffer->vertex(), buffer->pdgId(), statusFlag, true));
                target.push_back(*gluonPtr);
                // Increment & push index of the top daughter
                higgsDaughters.push_back(++motherPartIdx_);
                if(addRadiation_){
                    this->addRadiation(motherPartIdx_, buffer, target);
                }
            }
            
            // If particle is a photon
            else if(status==3 && absPdgId==22){
                std::auto_ptr<reco::GenParticle> photonPtr(new reco::GenParticle(buffer->threeCharge(), p4(buffer, statusFlag), buffer->vertex(), buffer->pdgId(), statusFlag, true));
                target.push_back(*photonPtr);
                // Increment & push index of the top daughter
                higgsDaughters.push_back(++motherPartIdx_);
                if(addRadiation_){
                    this->addRadiation(motherPartIdx_, buffer, target);
                }
            }
            
            // If particle is a Z boson
            else if(status==3 && absPdgId==23){
                std::auto_ptr<reco::GenParticle> zPtr(new reco::GenParticle(buffer->threeCharge(), p4(buffer, statusFlag), buffer->vertex(), buffer->pdgId(), statusFlag, true));
                target.push_back(*zPtr);
                // Increment & push index of the top daughter
                higgsDaughters.push_back(++motherPartIdx_);
                if(addRadiation_){
                    addRadiation(motherPartIdx_, buffer, target); 
                }
                // If further decay chain of Z becomes interesting, add it here (would need Z boson index to manage daughter refs)
            }
            
            // If particle is a W boson
            else if(status==3 && absPdgId==24){
                std::auto_ptr<reco::GenParticle> wPtr(new reco::GenParticle(buffer->threeCharge(), p4(buffer, statusFlag), buffer->vertex(), buffer->pdgId(), statusFlag, true));
                target.push_back(*wPtr);
                // Increment & push index of the top daughter
                higgsDaughters.push_back(++motherPartIdx_);
                if(addRadiation_){
                    addRadiation(motherPartIdx_, buffer, target); 
                }
                // If further decay chain of W becomes interesting, add it here (would need W boson index to manage daughter refs)
            }
            
            // Here would be the space to implement other Higgs decay particles and their decay chains (e.g. top, invisible, ...)
        }
        
        // Add potential sisters of the Higgs
        //std::cout<<"Number of mothers: "<<higgs->numberOfMothers()<<std::endl;
        //std::cout<<"\t"<<higgs->mother(0)->pdgId()<<" "<<higgs->mother(0)->status()<<" "<<higgs->mother(0)->numberOfDaughters()<<std::endl;
        //std::cout<<"\t"<<higgs->mother(1)->pdgId()<<" "<<higgs->mother(1)->status()<<" "<<higgs->mother(1)->numberOfDaughters()<<"\n\n";
        if(higgs->numberOfMothers()>0 && higgs->pdgId()==25){
            // Loop over all daughters of the Higgs mother i.e.
            // the Higgs and its potential sisters
            for(reco::GenParticle::const_iterator i_higgsSister = higgs->mother()->begin(); i_higgsSister!=higgs->mother()->end(); ++i_higgsSister){
                //std::cout<<"\t"<<i_higgsSister->pdgId()<<" "<<i_higgsSister->status()<<" "<<i_higgsSister->numberOfDaughters()<<std::endl;
                // Add all further particles but the Higgs and potential 
                // cases where the mother of the Higgs has itself as daughter
                if(i_higgsSister->pdgId()!=higgs->pdgId() && i_higgsSister->pdgId()!=higgs->mother()->pdgId()){
                    reco::GenParticle* cand = new reco::GenParticle(i_higgsSister->threeCharge(), i_higgsSister->p4(), i_higgsSister->vertex(), i_higgsSister->pdgId(), i_higgsSister->status(), false);
                    std::auto_ptr<reco::GenParticle> sisterPtr(cand);
                    target.push_back(*sisterPtr);
                    if(i_higgsSister->begin() != i_higgsSister->end()){
                        // In case the sister has daughters increment
                        // and add the first generation of daughters
                        this->addDaughters(++motherPartIdx_, i_higgsSister->begin(), target, false);
                    }
                }
            }
        }
        
        // Fill the map for the administration of daughter indices
        refs_[iHiggs] = higgsDaughters;
    }
}



// Calculate lorentz vector from input
reco::Particle::LorentzVector 
HiggsDecaySubset::p4(const reco::GenParticle::const_iterator part, int statusFlag)
{
    // Calculate the four vector for all higgs daughters from 
    // their daughters including additional radiation 
    if(statusFlag == 3){
        // Return 4 momentum as it is
        return part->p4();
    }
    reco::Particle::LorentzVector vec;
    for(reco::GenParticle::const_iterator p = part->begin(); p != part->end(); ++p){
        if(p->status()<=2 && std::abs(p->pdgId())==24){
            vec = p->p4();
        }
        else{
            if(p->status() <= 2){
                // Sum up the p4 of all stable particles 
                // (of status 1 or 2)
                vec += p->p4();
            }
            else{
                if(p->status() == 3){
                    // If the particle is unfragmented (i.e.
                    // status 3) descend by one level
                    vec += p4(p, statusFlag);   
                }
            }
        }
    }
    return vec;
}



// Fill vector including all radiations from quarks
void 
HiggsDecaySubset::addRadiation(int& idx, const reco::GenParticle::const_iterator part, reco::GenParticleCollection& target)
{
    std::vector<int> daughters;
    int idxBuffer = idx;
    for(reco::GenParticle::const_iterator daughter = part->begin(); daughter != part->end(); ++daughter){
        if(daughter->status()<=2 && daughter->pdgId()!=part->pdgId()){
            // Skip comment lines and make sure that
            // the particle is not double counted as 
            // daughter of itself
            std::auto_ptr<reco::GenParticle> ptr(new reco::GenParticle(daughter->threeCharge(), daughter->p4(), daughter->vertex(), daughter->pdgId(), daughter->status(), false));
            target.push_back(*ptr);
            daughters.push_back(++idx); // Push index of daughter
        }
    }
    if(daughters.size()){
        refs_[idxBuffer] = daughters;
    }
}



// Fill references for output vector
void 
HiggsDecaySubset::fillReferences(const reco::GenParticleRefProd& ref, reco::GenParticleCollection& sel)
{ 
    int idx = 0;
    for(reco::GenParticleCollection::iterator p = sel.begin(); p != sel.end(); ++p, ++idx){
        //find daughter reference vectors in refs_ and add daughters
        std::map<int, std::vector<int> >::const_iterator daughters=refs_.find(idx);
        if(daughters != refs_.end()){
            for(std::vector<int>::const_iterator daughter = daughters->second.begin(); daughter != daughters->second.end(); ++daughter){
                reco::GenParticle* part = dynamic_cast<reco::GenParticle*>(&(*p));
                if(part == 0)
                    throw edm::Exception( edm::errors::InvalidReference, "Not a GenParticle");
                part->addDaughter(reco::GenParticleRef(ref, *daughter));
                sel[*daughter].addMother(reco::GenParticleRef(ref, idx));
            }
        }
    }
}



// Recursively fill vector for all further decay particles of a given particle
void 
HiggsDecaySubset::addDaughters(int& idx, const reco::GenParticle::const_iterator part, reco::GenParticleCollection& target, const bool recursive)
{
    std::vector<int> daughters;
    int idxBuffer = idx;
    for(reco::GenParticle::const_iterator daughter = part->begin(); daughter != part->end(); ++daughter){
        std::auto_ptr<reco::GenParticle> ptr(new reco::GenParticle(daughter->threeCharge(), daughter->p4(), daughter->vertex(), daughter->pdgId(), daughter->status(), false));
        target.push_back(*ptr);
        // Increment & push index of daughter
        daughters.push_back(++idx);
        // Continue recursively if desired
        if(recursive){
            this->addDaughters(idx, daughter, target);  
        }
    }
    if(daughters.size()){
        refs_[idxBuffer] = daughters;
    }
}




// ------------ method called once each job just before starting event loop  ------------
void 
HiggsDecaySubset::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HiggsDecaySubset::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
HiggsDecaySubset::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
HiggsDecaySubset::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HiggsDecaySubset::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HiggsDecaySubset::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HiggsDecaySubset::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // prefer to set the parameters without a default value:
  // an exception is thrown when the parameter is not defined in the config files, instead of silently using the default given here
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<std::string>("fillMode");
  desc.add<bool>("addRadiation");
  descriptions.add("decaySubsetHiggs", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiggsDecaySubset);
