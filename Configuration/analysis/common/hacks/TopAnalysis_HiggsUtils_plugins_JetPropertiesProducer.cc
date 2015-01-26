// -*- C++ -*-
//
// Package:    JetPropertiesProducer
// Class:      JetPropertiesProducer
// 
/**\class JetPropertiesProducer JetPropertiesProducer.cc TopAnalysis/HiggsUtils/plugins/JetPropertiesProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Johannes Hauk,,,DESY
//         Created:  Tue Mar 12 14:29:43 CET 2013
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
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetupRecordKey.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/BaseTagInfo.h"

#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputerRecord.h"
#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputer.h"

#include "TopAnalysis/HiggsUtils/interface/JetProperties.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/Point3D.h"



//
// class declaration
//

class JetPropertiesProducer : public edm::EDProducer
{
public:
    explicit JetPropertiesProducer(const edm::ParameterSet&);
    ~JetPropertiesProducer();
    
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
private:
    virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    
    virtual void beginRun(edm::Run&, edm::EventSetup const&);
    virtual void endRun(edm::Run&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
    
    // ----------member data ---------------------------
    
    const edm::ParameterSet parameterSet_;
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
JetPropertiesProducer::JetPropertiesProducer(const edm::ParameterSet& iConfig):
parameterSet_(iConfig)
{
    produces<std::vector<JetProperties> >();
}



JetPropertiesProducer::~JetPropertiesProducer()
{}



//
// member functions
//

// ------------ method called to produce the data  ------------
void
JetPropertiesProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::auto_ptr<std::vector<JetProperties> > v_jetProperties(new std::vector<JetProperties>);
    
    
    edm::InputTag jets(parameterSet_.getParameter<edm::InputTag>("src"));
    
    edm::Handle<pat::JetCollection> jetHandle;
    //edm::Handle<edm::View< pat::Jet > > jetHandle;
    iEvent.getByLabel(jets, jetHandle);
    
    edm::InputTag pfs(parameterSet_.getParameter<edm::InputTag>("pfCands"));
    // This does not work anymore in CMSSW_7_2_X, since method GenericMVAJetTagComputer::passEventSetup() does not exist anymore
    // Get a JetTagComputer in order to retrieve the pT-corrected secondary vertex mass used by the CSV algorithm
    // Check that the ES record really exists
/*    const GenericMVAJetTagComputer* computer(0);
    const edm::eventsetup::EventSetupRecordKey eventSetupRecordKey = JetTagComputerRecord::keyForClass();
    if(iSetup.find(eventSetupRecordKey)){
        edm::ESHandle<JetTagComputer> computerHandle;
        edm::InputTag secondaryVertexComputer(parameterSet_.getParameter<edm::InputTag>("svComputer"));
        iSetup.get<JetTagComputerRecord>().get(secondaryVertexComputer.label(), computerHandle);
        computer = dynamic_cast<const GenericMVAJetTagComputer*>(computerHandle.product());
        computer->passEventSetup(iSetup);
    }*/
    
    // Get the interaction vertices
    edm::InputTag vertices(parameterSet_.getParameter<edm::InputTag>("primaryVertexInputTag"));
    edm::Handle<reco::VertexCollection> vertexHandle;
    
   iEvent.getByLabel(vertices, vertexHandle);
   const reco::VertexCollection *vertexCollection = vertexHandle.product();
    
    // The sum of the multiplicities of the jetSelectedTracks for all the jets before the jet that is currently analysed
    // This variable is used so that we can give the proper value to the jetSecondaryVertexTrackMatchToSelectedTrackIndex
    int jetSelectedTrackMultiplicity_total = 0; 
    
    // The sum of the multiplicities of the jetSecondaryVertex for all the jets before the jet that is currently analysed
    // This variable is used so that we can give the proper value to the jetSecondaryVertexTrackVertexIndex
    int jetSecondaryVertexMultiplicity_total = 0;
    
    // The sum of the multiplicities of the jetPfCandidateTrack for all the jets before the jet that is currently analysed
    // This variable is used so that we can give the proper value to the jetPfCandidatePtrToRecoTrackIndex
    int jetPfCandidateMultiplicity_total = 0;
    
    for(std::vector<pat::Jet>::const_iterator i_jet = jetHandle->begin(); i_jet != jetHandle->end(); ++i_jet){
        // PfCandidateTrack related variables definition
        std::vector<math::PtEtaPhiMLorentzVectorD> jetPfCandidateTrack;
        std::vector<int> jetPfCandidateTrackCharge;
        // The trackId represents the type of particle: 1->charged hadron, 2->electron, 3->muon
        std::vector<int> jetPfCandidateTrackId;
        std::vector<const reco::Track*> jetPfCandidatePtrToRecoTrack;
        std::vector<int>  jetPfCandidatePtrToRecoTrackIndex;
        
        // Jet charge as given by PAT (weighted by global pt)
        const double jetChargeGlobalPtWeighted(i_jet->jetCharge());
        
        // Jet charge as weighted sum
        // weight = projection of charged pflow object momentum on jet axis 
        double jetPx = i_jet->px();
        double jetPy = i_jet->py();
        double jetPz = i_jet->pz();
        
        const std::vector<reco::PFCandidatePtr> &pfConstituents = i_jet->getPFConstituents();
        
        double sumMomentum = 0.;
        double sumMomentumQ = 0.; 
        
        // Access trackIPTagInfo->impact parameter
        std::vector<double> jetSelectedTrackIPValue;
        std::vector<double> jetSelectedTrackIPSignificance;
        std::vector<int> jetSelectedTrackCharge;
        std::vector<math::PtEtaPhiMLorentzVectorD> jetSelectedTrack;
        
        // Access secondaryVertexTagInfo->secondaryVertex
        std::vector<math::PtEtaPhiMLorentzVectorD> jetSecondaryVertex;
        std::vector<int> jetSecondaryVertexTrackVertexIndex;
        std::vector<double> jetSecondaryVertexFlightDistanceValue;
        std::vector<double> jetSecondaryVertexFlightDistanceSignificance;
        
        // Access vertex information: some pfCandidates are associated to certain primary vertices in event reconstruction through a weight (>0). If no weight>0 found, a check via closest z distance is performed to match the pfCandidate to a vertex
        // The jetPfCandidatePrimaryVertexId variable is filled as follows: 
        // -1 -> more than one vertex associated to the same pfCandidate, 
        // 0 -> vertex zero associated through weight value, 
        // 1 -> vertex zero associated through min z distance, 
        // 2 -> vertex different from zero associated through min z distance, 
        // 3 -> vertex different from zero associated through weight value.
        std::vector<int> jetPfCandidatePrimaryVertexId;
        
        // pT-corrected secondary vertex mass used by the CSV algorithm (if there's no secondary vertex in the jet then it is set to -8888.)
        double jetSecondaryVertexPtCorrectedMass = -8888.;
        
        // Find the Index of the jetSelectedTracks that are matched to the jetSecondaryVertexTracks
        std::vector<int> jetSecondaryVertexTrackMatchToSelectedTrackIndex;
        
        for(std::vector<reco::PFCandidatePtr>::const_iterator i_candidate = pfConstituents.begin(); i_candidate != pfConstituents.end(); ++i_candidate){
            const int charge = (*i_candidate)->charge();
            
            if(charge == 0) continue;
            
            jetPfCandidateTrack.push_back( (*i_candidate)->polarP4());
            jetPfCandidateTrackCharge.push_back(charge);
            jetPfCandidateTrackId.push_back((*i_candidate)->particleId());
            
            const double constituentPx = (*i_candidate)->px();
            const double constituentPy = (*i_candidate)->py();
            const double constituentPz = (*i_candidate)->pz();
            const double product = constituentPx*jetPx + constituentPy*jetPy + constituentPz*jetPz;
            
            sumMomentum += product;
            sumMomentumQ += static_cast<double>(charge)*product;
            
            // Access reco::Track collection references
            reco::TrackRef refToPfTrack; 
            refToPfTrack = (*i_candidate)->trackRef();
            const reco::Track* pfCandPtr = refToPfTrack.get();
            
            // Store pfCandidate pointers for further matching to selectedTracks
            jetPfCandidatePtrToRecoTrack.push_back(&*pfCandPtr);
            jetPfCandidatePtrToRecoTrackIndex.push_back(i_candidate-pfConstituents.begin());
            
            // Access information about the vertex for the pfCandidate
            int vertexIndex = -1;
            unsigned int nFoundVertex = 0;
            float bestweight=0;
            double dzmin = 10000;
            
            for (reco::VertexCollection::const_iterator vtx = vertexCollection->begin();vtx != vertexCollection->end(); ++vtx)
            {
                // Look for the vertex with the maximum weight
                float w = vtx->trackWeight(refToPfTrack);
                if (w > bestweight)
                {
                    bestweight=w;
                    vertexIndex = vtx-vertexCollection->begin();
                    nFoundVertex++;
                }
            } 
            if (nFoundVertex==1 && vertexIndex==0) jetPfCandidatePrimaryVertexId.push_back(0);
                
            else if (nFoundVertex==1&& vertexIndex!=0) jetPfCandidatePrimaryVertexId.push_back(3);
                
            else if (nFoundVertex>1) jetPfCandidatePrimaryVertexId.push_back(-1);
               
            // If no maximum weight is found, look for the vertex closest in z
            else if (nFoundVertex<1)
            {
                double ztrack = (*i_candidate)->vertex().z();
                for(reco::VertexCollection::const_iterator vtx2 = vertexCollection->begin();vtx2 != vertexCollection->end(); ++vtx2)
                {
                    if(std::abs(ztrack - vtx2->z())<dzmin)
                    {
                        dzmin = std::abs(ztrack - vtx2->z());
                        vertexIndex = vtx2 - vertexCollection->begin();
                    }
                }
                if (vertexIndex!=0) jetPfCandidatePrimaryVertexId.push_back(2);
                else if (vertexIndex==0) jetPfCandidatePrimaryVertexId.push_back(1);
            }
        }
        
        const double jetChargeRelativePtWeighted(sumMomentum>0 ? sumMomentumQ/sumMomentum : 0);
        
        // Vector of pfCandidate indices of same size as selectedTracks filled when both pfCandidate and selectedTrack point to the same reco::Track
        // -1 means that there is no matching
        std::vector<int> jetSelectedTrackMatchToPfCandidateIndex; 
        
        if (i_jet->hasTagInfo("impactParameter")){
            
            const reco::TrackIPTagInfo* trackIPTagInfo = i_jet->tagInfoTrackIP("impactParameter");
            if (trackIPTagInfo != NULL){
                
                const reco::TrackRefVector jetSelectedTracks = trackIPTagInfo->selectedTracks();
                
                int iSelectedTrack=0;
                for(reco::TrackRefVector::const_iterator i_selectedTrack = jetSelectedTracks.begin(); i_selectedTrack != jetSelectedTracks.end(); ++i_selectedTrack){
                    
                    // Access the reco::track collection reference for the i_selectedTrack:
                    const reco::Track* selectedTrackPtr = i_selectedTrack->get();
                    
                    // Find among the pfCandidate pointers if there's any matching the studied selectedTrack
                    std::vector< const reco::Track* >::const_iterator foundTrack = std::find(jetPfCandidatePtrToRecoTrack.begin(), jetPfCandidatePtrToRecoTrack.end(), selectedTrackPtr);
                    int pfCandidateIndex = (foundTrack == jetPfCandidatePtrToRecoTrack.end()) ? -1 : (foundTrack - jetPfCandidatePtrToRecoTrack.begin())+jetPfCandidateMultiplicity_total;
                    
                    jetSelectedTrackMatchToPfCandidateIndex.push_back(pfCandidateIndex);
                    
                    const double ipValue = trackIPTagInfo->impactParameterData().at(iSelectedTrack).ip3d.value();
                    const double ipSignificance = trackIPTagInfo->impactParameterData().at(iSelectedTrack).ip3d.significance();
                    
                    jetSelectedTrackIPValue.push_back(ipValue);
                    jetSelectedTrackIPSignificance.push_back(ipSignificance);
                    jetSelectedTrackCharge.push_back((*i_selectedTrack)->charge());
                    
                    // Create the track LV
                    const double trackPt = (*i_selectedTrack)->pt();
                    const double trackEta = (*i_selectedTrack)->eta();
                    const double trackPhi = (*i_selectedTrack)->phi();
                    const double trackM = 0.13957018; //mass of the track is agreed to be the pion mass by default
                    
                    const math::PtEtaPhiMLorentzVectorD jetSelectedTrackLV (trackPt, trackEta, trackPhi, trackM);
                    jetSelectedTrack.push_back(jetSelectedTrackLV);
                    
                    if (i_jet->hasTagInfo("secondaryVertex")){
                        
                        const reco::SecondaryVertexTagInfo* secondaryVertexTagInfo = i_jet->tagInfoSecondaryVertex("secondaryVertex");
                        if (secondaryVertexTagInfo != NULL){
                            
                            unsigned int nVertex = secondaryVertexTagInfo->nVertices();
                           // This does not work anymore in CMSSW_7_2_X, since method GenericMVAJetTagComputer::passEventSetup() does not exist anymore 
                          /*  if(computer){
                                // Pass the TagInfos to the JetTagComputer
                                std::vector<const reco::BaseTagInfo*>  baseTagInfos;
                                const JetTagComputer::TagInfoHelper helper(baseTagInfos);
                                baseTagInfos.push_back(trackIPTagInfo);
                                baseTagInfos.push_back(secondaryVertexTagInfo);
                                const reco::TaggingVariableList vars = computer->taggingVariables(helper);
                                
                                // Retrieve the pT-corrected secondary vertex mass used by the CSV algorithm
                                if(nVertex > 0)
                                    jetSecondaryVertexPtCorrectedMass = vars.checkTag(reco::btau::vertexMass) ? vars.get(reco::btau::vertexMass) : -9999.;
                            }
                            else jetSecondaryVertexPtCorrectedMass = -9990.;
                            */
                            for(size_t iVertex=0; iVertex<nVertex; ++iVertex){   
                                
                                const reco::Vertex& secVertex = secondaryVertexTagInfo->secondaryVertex(iVertex);
                                
                                // Create the Secondary Vertex LV
                                double secondaryVertexP4Pt = secVertex.p4().pt();
                                double secondaryVertexP4Eta = secVertex.p4().eta();
                                double secondaryVertexP4Phi = secVertex.p4().phi();
                                double secondaryVertexP4Mass = secVertex.p4().mass();
                                
                                const math::PtEtaPhiMLorentzVectorD jetSecondaryVertexLV (secondaryVertexP4Pt, secondaryVertexP4Eta, secondaryVertexP4Phi, secondaryVertexP4Mass);
                                
                                const double flightDistanceValue = secondaryVertexTagInfo->flightDistance(iVertex).value();
                                const double flightDistanceSignificance = secondaryVertexTagInfo->flightDistance(iVertex).significance();
                                
                                if (iSelectedTrack==0){
                                    jetSecondaryVertex.push_back(jetSecondaryVertexLV);
                                    jetSecondaryVertexFlightDistanceValue.push_back(flightDistanceValue);
                                    jetSecondaryVertexFlightDistanceSignificance.push_back(flightDistanceSignificance);
                                }
                                
                                const reco::TrackRefVector jetSecondaryVertexTracks = secondaryVertexTagInfo->vertexTracks(iVertex);
                                
                                for(reco::TrackRefVector::const_iterator i_secondaryVertexTrack = jetSecondaryVertexTracks.begin(); i_secondaryVertexTrack != jetSecondaryVertexTracks.end(); ++i_secondaryVertexTrack){
                                    
                                    // Access the reco::track collection reference for the i_secondaryVertexTrack:
                                    const reco::Track* secondaryVertexTrackPtr = i_secondaryVertexTrack->get();
                                    
                                    // Keep only the selectedTracks that are matched to the i_secondaryVertexTrack
                                    if(secondaryVertexTrackPtr != selectedTrackPtr) continue;
                                    
                                    // Assign the proper index to the jetSecondaryVertexTracks in order to find the corresponding index of the matched jetSelectedTracks
                                    jetSecondaryVertexTrackMatchToSelectedTrackIndex.push_back(iSelectedTrack+jetSelectedTrackMultiplicity_total);
                                    
                                    // Assign the proper index for the jetSecondaryVertexTracks in order to find the index of the Secondary Vertex they belong to
                                    jetSecondaryVertexTrackVertexIndex.push_back(iVertex+jetSecondaryVertexMultiplicity_total);
                                    
                                    // Since the matching between the jetSecondaryVertexTracks and the jetSelectedTracks is unique, as soon as we find one matching pair,
                                    // there's no need to search for an other track of the jetSecondaryVertexTracks matched to the same track of the jetSelectedTracks
                                    break;
                                }  
                            }
                            
                        }
                        
                    }
                    ++iSelectedTrack;
                }
                
            }
            
        }
        
        // Find the total #SelectedTracks per event for all the jets
        jetSelectedTrackMultiplicity_total = jetSelectedTrackMultiplicity_total + jetSelectedTrack.size();
        
        // Find the total #SecondaryVertices per event for all the jets
        jetSecondaryVertexMultiplicity_total = jetSecondaryVertexMultiplicity_total + jetSecondaryVertex.size();
        
        // Find the total #pfCandidates per event for all the jets
        jetPfCandidateMultiplicity_total = jetPfCandidateMultiplicity_total + jetPfCandidateTrack.size();
        
        // Access Lorentz vector and PDG ID of parton associated to jet by PAT
        // If it does not exist, this can be identified by PDG ID =0
        int jetAssociatedPartonPdgId(0);
        math::PtEtaPhiMLorentzVectorD jetAssociatedParton(0.,0.,0.,0.);
        const reco::GenParticle * genParton = i_jet->genParton();
        if(genParton){
            jetAssociatedPartonPdgId = genParton->pdgId();
            jetAssociatedParton = genParton->polarP4();
        }
        
        JetProperties jetProperties(jetChargeGlobalPtWeighted, jetChargeRelativePtWeighted, jetAssociatedPartonPdgId, jetAssociatedParton, jetPfCandidateTrack, jetPfCandidateTrackCharge,jetPfCandidateTrackId, jetPfCandidatePrimaryVertexId, jetSelectedTrackMatchToPfCandidateIndex, jetSelectedTrack, jetSelectedTrackIPValue, jetSelectedTrackIPSignificance, jetSelectedTrackCharge, jetSecondaryVertexTrackMatchToSelectedTrackIndex, jetSecondaryVertexTrackVertexIndex, jetSecondaryVertex, jetSecondaryVertexFlightDistanceValue, jetSecondaryVertexFlightDistanceSignificance, jetSecondaryVertexPtCorrectedMass);
        v_jetProperties->push_back(jetProperties);
        
        edm::LogVerbatim log("JetPropertiesProducer");
        log<<"   ---   Jet Properties   ---   \n";
        log<<"Jet charge global pt weighted:   "<<jetProperties.jetChargeGlobalPtWeighted()<<"\n";
        log<<"Jet charge relative pt weighted: "<<jetProperties.jetChargeRelativePtWeighted()<<"\n";
        log<<"Jet associated parton PDG ID:    "<<jetProperties.jetAssociatedPartonPdgId()<<"\n";
        log<<"Jet associated parton pt: "<<jetProperties.jetAssociatedParton().pt()<<"\n";
        log<<"Jet associated parton eta: "<<jetProperties.jetAssociatedParton().eta()<<"\n";
        log<<"Jet associated parton phi: "<<jetProperties.jetAssociatedParton().phi()<<"\n";
        log<<"Jet associated parton mass: "<<jetProperties.jetAssociatedParton().M()<<"\n";
        log<<"   --------------------------   \n\n";    
        
        jetPfCandidateTrack.clear();
        jetPfCandidateTrackCharge.clear();
        jetPfCandidateTrackId.clear();
        jetSelectedTrackMatchToPfCandidateIndex.clear();
        jetSelectedTrack.clear();
        jetSelectedTrackIPValue.clear();
        jetSelectedTrackIPSignificance.clear();
        jetSelectedTrackCharge.clear();
        jetSecondaryVertex.clear();
        jetSecondaryVertexFlightDistanceValue.clear();
        jetSecondaryVertexFlightDistanceSignificance.clear();
        jetSecondaryVertexTrackVertexIndex.clear();
        jetSecondaryVertexTrackMatchToSelectedTrackIndex.clear();
        jetPfCandidatePrimaryVertexId.clear();
    }
    
    iEvent.put(v_jetProperties);
    
}



// ------------ Method called once each job just before starting event loop  ------------
void 
JetPropertiesProducer::beginJob()
{}



// ------------ Method called once each job just after ending the event loop  ------------
void 
JetPropertiesProducer::endJob()
{}



// ------------ Method called when starting to processes a run  ------------
void 
JetPropertiesProducer::beginRun(edm::Run&, edm::EventSetup const&)
{}



// ------------ Method called when ending the processing of a run  ------------
void 
JetPropertiesProducer::endRun(edm::Run&, edm::EventSetup const&)
{}



// ------------ Method called when starting to processes a luminosity block  ------------
void 
JetPropertiesProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{}



// ------------ Method called when ending the processing of a luminosity block  ------------
void 
JetPropertiesProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{}



// ------------ Method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetPropertiesProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    // Prefer to set the parameters without a default value:
    // An exception is thrown when the parameter is not defined in the config files, instead of silently using the default given here
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src");
    desc.add<edm::InputTag>("pfCands");
    desc.add<edm::InputTag>("svComputer");
    desc.add<edm::InputTag>("primaryVertexInputTag");
    descriptions.add("jetProperties", desc);
}



// Define this as a plug-in
DEFINE_FWK_MODULE(JetPropertiesProducer);
