// -*- C++ -*-
//
// Package:    MadgraphWDecayProducer
// Class:      MadgraphWDecayProducer
// 
/**\class MadgraphWDecayProducer MadgraphWDecayProducer.cc TopAnalysis/TopUtils/plugins/MadgraphWDecayProducer.cc

 Description: Stores for MadGraph samples the decay mode of each matrix-element-level W boson,
              in order to correct branching ratios (which are set all equal in MadGraph to 1/9)

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Johannes Hauk,01a/O2.115,3139,
//         Created:  Wed Jul 16 14:22:23 CEST 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TopQuarkAnalysis/TopEventProducers/interface/TopDecaySubset.h"





//
// class declaration
//

class MadgraphWDecayProducer : public edm::EDProducer{
public:
    explicit MadgraphWDecayProducer(const edm::ParameterSet&);
    ~MadgraphWDecayProducer();
    
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
private:
    virtual void beginJob();
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob();
    
    virtual void beginRun(edm::Run&, edm::EventSetup const&);
    virtual void endRun(edm::Run&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
    
    /// find Higgs in list of input particles
    std::vector<const reco::GenParticle*> findW(const reco::GenParticleCollection& v_particle)const;
    
    /// Check decay chain for the exploited shower model
    TopDecaySubset::ShowerModel checkShowerModel(const std::vector<const reco::GenParticle*>& v_particle)const;
    
    // ----------member data ---------------------------
    
    /// Input tag for the genParticle source
    const edm::InputTag genParticleSource_;
    
    /// Parton shower mode (filled in checkShowerModel)
    TopDecaySubset::ShowerModel showerModel_;
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
MadgraphWDecayProducer::MadgraphWDecayProducer(const edm::ParameterSet& iConfig):
genParticleSource_(iConfig.getParameter<edm::InputTag>("src")),
showerModel_(TopDecaySubset::kStart)
{
    produces<std::vector<int> >("madgraphWDecay");
    //produces<std::vector<int> >();
}



MadgraphWDecayProducer::~MadgraphWDecayProducer()
{}



//
// member functions
//


// ------------ method called to produce the data  ------------
void
MadgraphWDecayProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // Create target vector
    std::auto_ptr<std::vector<int> > target(new std::vector<int>);
    
    // Get source collection
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel(genParticleSource_, genParticles);
    
    // Find W's in list of input particles
    const std::vector<const reco::GenParticle*> v_w = this->findW(*genParticles);
    
    // Determine shower model (only in first event containing a W)
    if(v_w.size() && showerModel_==TopDecaySubset::kStart) showerModel_ = this->checkShowerModel(v_w);
    
    // Find W decay modes
    for(std::vector<const reco::GenParticle*>::const_iterator i_particle = v_w.begin(); i_particle != v_w.end(); ++i_particle){
        const reco::GenParticle* particle = *i_particle;
        
        // Sanity check
        if(showerModel_!=TopDecaySubset::kPythia && particle->begin()==particle->end())
            throw edm::Exception(edm::errors::LogicError,
                "showerModel_!=kPythia && particle->begin()==particle->end()\n");
        
        // Iterate over W daughters
        for(reco::GenParticle::const_iterator i_daughter = ((showerModel_==TopDecaySubset::kPythia)?particle->begin():particle->begin()->begin()); i_daughter != ((showerModel_==TopDecaySubset::kPythia)?particle->end():particle->begin()->end()); ++i_daughter){
            //std::cout<<"\t"<<i_daughter->pdgId()<<" "<<i_daughter->status()<<" "<<i_daughter->numberOfDaughters()<<std::endl;
            const int absPdgId = std::abs(i_daughter->pdgId());
            if(absPdgId == 24) continue;
            else if(absPdgId>=1 && absPdgId<=6){ // hadronic decay
                target->push_back(1);
                break;
            }
            else if(absPdgId>=11 && absPdgId<=12){ // electronic decay
                target->push_back(2);
                break;
            }
            else if(absPdgId>=13 && absPdgId<=14){ // muonic decay
                target->push_back(3);
                break;
            }
            else if(absPdgId>=15 && absPdgId<=16){ // tauonic decay
                target->push_back(4);
                break;
            }
            else{ // decay not identified
                edm::LogInfo("MadgraphWDecayProducer")<<"Decay of W could not be classified\n"
                    <<"PDG ID of daughter is: "<<i_daughter->pdgId();
                target->push_back(0);
                break;
            }
        }
    }
    
    iEvent.put(target, "madgraphWDecay");
}



std::vector<const reco::GenParticle*> MadgraphWDecayProducer::findW(const reco::GenParticleCollection& v_particle)const
{
    std::vector<const reco::GenParticle*> v_w;
    
    for(reco::GenParticleCollection::const_iterator i_particle = v_particle.begin(); i_particle != v_particle.end(); ++i_particle){
        if(std::abs(i_particle->pdgId())==24 && i_particle->status()==3){
            //std::cout<<"\t"<<i_particle->pdgId()<<" "<<i_particle->status()<<" "<<i_particle->numberOfMothers()<<" "<<i_particle->numberOfDaughters()<<std::endl;
            v_w.push_back( &(*i_particle) );
        }
    }
    
    return v_w;
}



TopDecaySubset::ShowerModel MadgraphWDecayProducer::checkShowerModel(const std::vector<const reco::GenParticle*>& v_particle)const
{
    for(std::vector<const reco::GenParticle*>::const_iterator i_particle = v_particle.begin(); i_particle != v_particle.end(); ++i_particle){
        const reco::GenParticle* particle = *i_particle;
        // check for kHerwig type showers: here the status 3 W will 
        // have a single status 2 W as daughter, which has again 3 
        // or more status 2 daughters
        // WARNING: I did not explicitely check for W in HERWIG, only in PYTHIA, so HERWIG criterion can be wrong
        if(particle->numberOfDaughters() == 1){
            if(particle->begin()->pdgId()==particle->pdgId() && particle->begin()->status()==2 && particle->begin()->numberOfDaughters()>1)
                return TopDecaySubset::kHerwig;
        }
        // Check for kPythia type showers: here the status 3 W will 
        // have all decay products and a status 2 W as daughters
        // the status 2 W will be w/o further daughters
        // (or at maximum one as observed in specific ttbar MadSpin samples -- it has itself once more as status 2 daughter)
        if(particle->numberOfDaughters() > 1){
            bool daughterSameHasNoDaughters(false);
            for(reco::GenParticle::const_iterator i_daughter = particle->begin(); i_daughter != particle->end(); ++i_daughter){
                //std::cout<<"\t"<<i_daughter->pdgId()<<" "<<i_daughter->status()<<" "<<i_daughter->numberOfDaughters()<<std::endl;
                //if(std::abs(i_daughter->pdgId())==24 && i_daughter->status()==2 && i_daughter->numberOfDaughters()==0) daughterSameHasNoDaughters = true;
                if(std::abs(i_daughter->pdgId())==24 && i_daughter->status()==2 && i_daughter->numberOfDaughters()<=1) daughterSameHasNoDaughters = true;
            }
            if(daughterSameHasNoDaughters) return TopDecaySubset::kPythia;
        }
    }
    
    // If neither Herwig nor Pythia like
    if(v_particle.size() == 0)
        // This should never be reached, since function is only called in case W's were found
        edm::LogWarning("decayChain")
            <<" Failed to find W in decay chain. Will assume that this a \n"
            <<" non-W sample and produce an empty decaySubset.\n";
    else
        throw edm::Exception(edm::errors::LogicError,
            " Can not find back any of the supported hadronization models. Models \n"
            " which are supported are:                                            \n"
            " Pythia  LO(+PS): W(status 3) --> particle(status 3) antiparticle(status 3)\n"
            " Herwig NLO(+PS): W(status 2) --> W(status 3) --> W(status 2)  \n");
    return TopDecaySubset::kNone;
}



// ------------ method called once each job just before starting event loop  ------------
void 
MadgraphWDecayProducer::beginJob()
{}



// ------------ method called once each job just after ending the event loop  ------------
void 
MadgraphWDecayProducer::endJob()
{}



// ------------ method called when starting to processes a run  ------------
void 
MadgraphWDecayProducer::beginRun(edm::Run&, edm::EventSetup const&)
{}



// ------------ method called when ending the processing of a run  ------------
void 
MadgraphWDecayProducer::endRun(edm::Run&, edm::EventSetup const&)
{}



// ------------ method called when starting to processes a luminosity block  ------------
void 
MadgraphWDecayProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{}



// ------------ method called when ending the processing of a luminosity block  ------------
void 
MadgraphWDecayProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MadgraphWDecayProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src");
    descriptions.add("MadgraphWDecayProducer", desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(MadgraphWDecayProducer);
