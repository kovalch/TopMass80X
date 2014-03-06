/**\class GenHFHadronAnalyzer GenHFHadronAnalyzer.cc
* @brief Finds the origin of each heavy flavour hadron and associated jets to it
*
* Starting from each consituent of each jet, tracks back in chain to find heavy flavour hadrons.
* From each hadron traces back until finds the b quark and its mother.
* For each hadron identifies the jet to which it was injected as a ghost hadron.
*
* The description of the run-time parameters can be found at fillDescriptions()
*
* The description of the products can be found at GenHFHadronAnalyzer()
*/

// Original Author:  Nazar Bartosik,DESY

// system include files
#include <memory>
#include <algorithm>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// added by me
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"

//
// class declaration
//

class GenHFHadronAnalyzer : public edm::EDAnalyzer
{
public:
    explicit GenHFHadronAnalyzer ( const edm::ParameterSet& );
    ~GenHFHadronAnalyzer();

//    static void fillDescriptions ( edm::ConfigurationDescriptions& descriptions );

private:
    virtual void beginJob() ;
    virtual void analyze( const edm::Event&, const edm::EventSetup& );
    virtual void endJob() ;

    virtual void beginRun ( edm::Run&, edm::EventSetup const& );
    virtual void endRun ( edm::Run&, edm::EventSetup const& );
    virtual void beginLuminosityBlock ( edm::LuminosityBlock&, edm::EventSetup const& );
    virtual void endLuminosityBlock ( edm::LuminosityBlock&, edm::EventSetup const& );

    std::vector<int> findHadronJets ( const reco::GenJetCollection& genJets,  std::vector<int> &hadIndex, std::vector<reco::GenParticle> &hadMothersGenPart, 
                                      std::vector<std::vector<int> > &hadMothersIndices, std::vector<int> &hadLeptonIndex, 
                                      std::vector<int> &hadLeptonHadIndex, std::vector<int> &hadFlavour, 
                                      std::vector<int> &hadFromTopWeakDecay, std::vector<int> &hadFromBHadron );
    typedef const reco::Candidate* pCRC;
    int analyzeMothers ( const reco::Candidate* thisParticle, pCRC *hadron, pCRC *lepton, int& topDaughterQId, int& topBarDaughterQId, 
                         std::vector<const reco::Candidate*> &hadMothers, std::vector<std::vector<int> > &hadMothersIndices, 
                         std::set<const reco::Candidate*> *analyzedParticles, const int prevPartIndex );
    bool putMotherIndex ( std::vector<std::vector<int> > &hadMothersIndices, int partIndex, int mothIndex );
    bool isHadron ( const int flavour, const reco::Candidate* thisParticle );
    bool isHadronPdgId ( const int flavour, const int pdgId );
    bool hasHadronDaughter ( const int flavour, const reco::Candidate* thisParticle );
    int isInList ( std::vector<const reco::Candidate*> particleList, const reco::Candidate* particle );
    int isInList ( std::vector<int> list, const int value );
    bool findInMothers (int idx, std::vector<int> &mothChains, std::vector<std::vector<int> > &hadMothersIndices, 
                        std::vector<reco::GenParticle> &hadMothers, int status, int pdgId, bool pdgAbs, int stopId, int firstLast, bool verbose );
    bool isNeutralPdg ( int pdgId );

    bool fixExtraSameFlavours(
        const unsigned int hadId, const std::vector<int> &hadIndices, const std::vector<reco::GenParticle> &hadMothers, 
        const std::vector<std::vector<int> > &hadMothersIndices, const std::vector<int> &isFromTopWeakDecay, 
        const std::vector<std::vector<int> > &LastQuarkIds, const std::vector<std::vector<int> > &LastQuarkMotherIds, 
        std::vector<int> &lastQuarkIndices, std::vector<int> &hadronFlavour, std::set<int> &checkedHadronIds, const int lastQuarkIndex);

// ----------member data ---------------------------
    edm::InputTag ttGenEvent_;
    edm::InputTag genJets_;
    int flavour_;
    bool noBBbarResonances_;
    bool onlyJetClusteredHadrons_;
    bool doValidationPlotsForImprovedHadronMatching_;
    bool toPlot;

    std::string flavourStr_;  // Name of the flavour specified in config file

    edm::ESHandle<ParticleDataTable> pdt_;


// -- testing histograms
    TH1I *h_nHad, *h_hadJetIndex;
    TH2D *h2_Q2MothPdg, *h2_Q2Moth_2Pdg;
    TH1D *h_jetdPt, *h_jetDeltaPt;
    TH2D *h2_QQMothSamePdg;
    TH1D *h_QQMothEFrac, *h_QHQPtFrac, *h_QQMothdR, *h_QQMothdR_min, *h_QHQdR;
    TH1D *h_QQMothEFrac_ambig, *h_QHQPtFrac_ambig, *h_QQMothdR_ambig, *h_QQMothdR_min_ambig, *h_QHQdR_ambig;
    TH1D *h_jetHad_dR, *h_jetHadOld_dR, *h_jetHadOldNoMatch_dR, *h_jetHad_dR_noFlav, *h_jetHadOldNoMatch_dR_noFlav;

    TH1I *h_nFirstQ,*h_firstQstatus, *h_firstQPdg;
    TH2D *h2_firstQFlav_HadFlav, *h2_firstQFlav_QMotherFlav;
    TH1I *h_nLastQ,*h_lastQstatus, *h_lastQPdg, *h_nHadTop, *h_nHadHiggs, *h_nHadGluon, *h_nHadZ, *h_nHadNo;
    TH2D *h2_lastQFlav_HadFlav, *h2_lastQFlav_QMotherFlav, *h2_lastQPdg_QMotherPdg;
    TH1D *h_lastQ_lastQ_dRmin, *h_hadJet_dRmin, *h_lastQ12_had_dRratio;
    TH1D *h_dEHadJet_T, *h_dEHadJet_H, *h_dEHadJet_Z, *h_dEHadJet_G, *h_dEHadJet_P, *h_dEHadJet_Q, *h_dEHadJet_Other, *h_dEHadJet[7];
    TH1I *h_nHadInJet, *h_hadFlavour;
    TH1I *h_nLeptonsInHad, *h_nHadNotFromTDecay, *h_nHadNotFromTHDecay, *h_hadIsFromBHadron;
    TH1I *h_nJetNotFromTDecay;
    TH1I *h_topDaughterQuarkFlavour, *h_hadFromTDecayPdg, *h_hadNotFromTDecayPdg;
    TH1D *h_hadNotClusteredInJet_Pt, *h_hadNoDecayClusteredInJet_Pt, *h_hadNoDecayButHadClusteredInJet_Pt, *h_hadAll_Pt;


};



//
// constructors and destructor
//
/**
* @brief constructor initialising producer products and config parameters
*
* Output generated by this producer:
* <TABLE>
* <TR><TH> name                </TH><TH> type                           </TH><TH> description </TH> </TR>
* <TR><TD> BHadJetIndex        </TD><TD> std::vector<int>               </TD><TD> position of the jet identified as b-hadron jet in the input GenJetCollection </TD></TR>
* <TR><TD> BHadrons            </TD><TD> std::vector<reco::GenParticle> </TD><TD> vector of the identified b-hadrons (definition of b-hadron and anti-b-hadron below) </TD></TR>
* <TR><TD> BHadronFromTopB     </TD><TD> std::vector<bool>              </TD><TD> true if the corresponding b-hadron originates from the ttGenEvent b-quark </TD></TR>
* <TR><TD> BHadronVsJet        </TD><TD> std::vector<int>               </TD><TD> matrix of which b-hadron appears in which GenJet, access by [iJet * BHadrons.size() + iBhadron] </TD></TR>
*
* <TR><TD> AntiBHadJetIndex    </TD><TD> std::vector<int>               </TD><TD> position of the jet identified as anti-b-hadron jet in the input GenJetCollection </TD></TR>
* <TR><TD> AntiBHadrons        </TD><TD> std::vector<reco::GenParticle> </TD><TD> vector of the identified anti-b-hadrons (definition of b-hadron and anti-b-hadron below) </TD></TR>
* <TR><TD> AntiBHadronFromTopB </TD><TD> std::vector<bool>              </TD><TD> true if the corresponding anti-b-hadron originates from the ttGenEvent anti-b-quark </TD></TR>
* <TR><TD> AntiBHadronVsJet    </TD><TD> std::vector<int>               </TD><TD> matrix of which anti-b-hadron appears in which GenJet, access by [iJet * AntiBHadrons.size() + iBhadron] </TD></TR>
*
* </TABLE>
*
* @warning Definition of b-hadron and anti-b-hadron: The term b-hadron and anti-b-hadron is in reference to the quark content and <b>not</b> to distinguish particles from anti-particles.
* Here a b-hadron contains a b-quark and an anti-b-hadron contains an anti-b-quark.
* For mesons this means an inversion with respect to the PDG definition, as mesons actually contain anti-b-quarks and anti-mesons contain b-quarks.
*
*/
GenHFHadronAnalyzer::GenHFHadronAnalyzer ( const edm::ParameterSet& cfg )
{

    ttGenEvent_        = cfg.getParameter<edm::InputTag> ( "ttGenEvent" );
    genJets_           = cfg.getParameter<edm::InputTag> ( "genJets" );
    flavour_           = cfg.getParameter<int> ( "flavour" );
    onlyJetClusteredHadrons_ = cfg.getParameter<bool> ( "onlyJetClusteredHadrons" );
    doValidationPlotsForImprovedHadronMatching_ = cfg.getParameter<bool> ( "doValidationPlotsForImprovedHadronMatching" );

    toPlot = doValidationPlotsForImprovedHadronMatching_;


    flavour_ = abs ( flavour_ ); // Make flavour independent of sign given in configuration
    if ( flavour_==5 ) {
        flavourStr_="B";
    } else if ( flavour_==4 ) {
        flavourStr_="C";
    } else {
        edm::LogError ( "GenHFHadronAnalyzer" ) << "Flavour option must be 4 (c-jet) or 5 (b-jet), but is: " << flavour_ << ". Correct this!";
    }


    if ( !doValidationPlotsForImprovedHadronMatching_ ) {
        return;
    }

// Initializing ROOT file for storing histograms
    edm::Service<TFileService> fs;
    if( !fs ) throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

    h_nHad = fs->make<TH1I> ( "nHad", "Nr. of hadrons;N", 10, 0, 10 );
    h_hadJetIndex = fs->make<TH1I> ( "hadJetIndex", "Index of hadron jet;index;hadrons", 102, -2, 100 );
    h_jetdPt = fs->make<TH1D> ( "jetdPt", "#DeltaPt(Clust/dR);#DeltaPt(Clust/dR)", 50, 0, 3 );
    h_jetDeltaPt = fs->make<TH1D> ( "jetDeltaPt", "#DeltaPt(Clust-dR);#DeltaPt [GeV]", 50, -100, 100 );

    h2_QQMothSamePdg = fs->make<TH2D> ( "QQMothSamePdg", "PdgId of quark and mother quark with same abs(pdgId);Daughter PdgId;Mother PdgId", 14, -7, 7,14,-7,7 );
    h_QQMothEFrac = fs->make<TH1D> ( "QQMothEFrac", "E fraction of daughter;#frac{E}{E^{mother}}", 50, 0., 10.0 );
    h_QHQPtFrac = fs->make<TH1D> ( "QHQPtFrac", "Pt fraction of highest Pt daughter;#frac{Pt^{highest}-Pt}{Pt^{highest}}", 50, 0., 1.01 );
    h_QQMothdR = fs->make<TH1D> ( "QQMothdR", "dR between mother and same pdgId daughters;dR", 100, 0.0, 5.0 );
    h_QQMothdR_min = fs->make<TH1D> ( "QQMothdR_min", "Min. dR between mother and same pdgId daughters;dR", 100, 0.0, 5.0 );
    h_QHQdR = fs->make<TH1D> ( "QHQdR", "dR between mother and highest Pt daughter;dR", 100, 0.0, 5.0 );

    h_QQMothEFrac_ambig = fs->make<TH1D> ( "QQMothEFrac_ambig", "E fraction of daughter;#frac{E}{E^{mother}}", 50, 0., 10.0 );
    h_QHQPtFrac_ambig = fs->make<TH1D> ( "QHQPtFrac_ambig", "Pt fraction of highest Pt daughter;#frac{Pt^{highest}-Pt}{Pt^{highest}}", 50, 0., 1.01 );
    h_QQMothdR_ambig = fs->make<TH1D> ( "QQMothdR_ambig", "dR between mother and same pdgId daughters;dR", 100, 0.0, 5.0 );
    h_QQMothdR_min_ambig = fs->make<TH1D> ( "QQMothdR_min_ambig", "Min. dR between mother and same pdgId daughters;dR", 100, 0.0, 5.0 );
    h_QHQdR_ambig = fs->make<TH1D> ( "QHQdR_ambig", "dR between mother and highest Pt daughter;dR", 100, 0.0, 5.0 );

    h_jetHad_dR = fs->make<TH1D> ( "jetHad_dR", "dR between hadron and jet matched by jet clustering algo;dR_{hadron}^{jet}", 50, 0, 4 );
    h_jetHadOld_dR = fs->make<TH1D> ( "jetHadOld_dR", "dR between hadron and jet matched by decay products;dR_{hadron}^{jet}", 50, 0, 4 );
    h_jetHadOldNoMatch_dR = fs->make<TH1D> ( "jetHadOldNoMatch_dR", "dR between hadron and jet matched by decay products;dR_{hadron}^{jet}", 50, 0, 4 );
    h_jetHad_dR_noFlav = fs->make<TH1D> ( "jetHad_dR_noFlav", "dR between hadron and jet matched by jet clustering algo;dR_{hadron}^{jet}", 50, 0, 4 );
    h_jetHadOldNoMatch_dR_noFlav = fs->make<TH1D> ( "jetHadOldNoMatch_dR_noFlav", "dR between hadron and jet matched by decay products;dR_{hadron}^{jet}", 50, 0, 4 );

    h_hadJet_dRmin = fs->make<TH1D> ( "hadJet_dRmin", "Min dR between hadron and any jet;dR_{hadron}^{jet}", 50, 0, 4 );

    h_nFirstQ = fs->make<TH1I> ( "nFirstQ", "Number of the 1-st quarks;N", 10, 0, 10 );
    h_firstQstatus = fs->make<TH1I> ( "firstQstatus", "Status of the 1-st quark;status", 5, 0, 5 );
    h_firstQPdg = fs->make<TH1I> ( "firstQPdg", "PdgId of the 1-st quark;PdgId", 60, -30, 30 );
    h2_firstQFlav_HadFlav = fs->make<TH2D> ( "firstQFlav_HadFlav", "Charge of the 1-st quark and hadron;Hadron charge;Quark charge", 5, -2, 3, 5,-2,3 );
    h2_firstQFlav_QMotherFlav = fs->make<TH2D> ( "firstQMotherFlav_QMotherFlav", "Charge of the 1-st quark and its mother;Quark charge;Quark's mother charge", 5, -2, 3, 5, -2, 3 );
    h_nLastQ = fs->make<TH1I> ( "nLastQ", "Number of the last quarks;N", 10, 0, 10 );
    h_lastQstatus = fs->make<TH1I> ( "lastQstatus", "Status of the last quark;status", 5, 0, 5 );
    h_lastQPdg = fs->make<TH1I> ( "lastQPdg", "PdgId of the last quark;PdgId", 60, -30, 30 );
    h2_lastQFlav_HadFlav = fs->make<TH2D> ( "lastQFlav_HadFlav", "Charge of the last quark and hadron;Hadron charge;Quark charge", 5, -2, 2,5,-2,2 );
    h2_lastQFlav_QMotherFlav = fs->make<TH2D> ( "lastQFlav_QMotherFlav", "Charge of the last quark and its mother;Quark charge;Quark's mother charge", 5, -2, 3,5,-2,3 );
    h2_lastQPdg_QMotherPdg = fs->make<TH2D> ( "lastQPdg_QMotherPdg", "PdgId of the last quark and its mother;Quark pdgId;Quark's mother pdgId", 13, -6, 7,35,-8,27 );
    h_lastQ_lastQ_dRmin = fs->make<TH1D> ( "lastQ_lastQ_dRmin", "dR_{min} between hadron and closest last quark;dR_{hadron}^{Quark}", 50, 0, 4 );
    h_lastQ12_had_dRratio = fs->make<TH1D> ( "lastQ12_had_dRratio", "(dR2-dR1)/dR2 between hadron and 2 closest last quarks;(dR2-dR1)/dR2;Hadrons", 40, 0, 1 );

    h_nHadTop = fs->make<TH1I> ( "nHadTop", "Nr. of hadrons from Top;N", 10, 0, 10 );
    h_nHadHiggs = fs->make<TH1I> ( "nHadHiggs", "Nr. of hadrons from Higgs;N", 10, 0, 10 );
    h_nHadGluon = fs->make<TH1I> ( "nHadGluon", "Nr. of hadrons from Gluon;N", 10, 0, 10 );
    h_nHadZ = fs->make<TH1I> ( "nHadZ", "Nr. of hadrons from Z;N", 10, 0, 10 );
    h_nHadNo = fs->make<TH1I> ( "nHadNo", "Nr. of hadrons without flavour;N", 10, 0, 10 );

    h_nHadInJet = fs->make<TH1I> ( "nHadInJet", "N hadrons associated to a jet;N hadrons;Jets", 10, 0, 10 );
    h_hadFlavour = fs->make<TH1I> ( "hadFlavour", "Flavour of the hadron origin;Flavour;Hadrons", 60, -30, 30 );
    
    h_nLeptonsInHad = fs->make<TH1I> ( "nLeptonInHad", "Number of leptons in hadron;N;Hadrons", 10, 0, 10 );
    h_nHadNotFromTDecay = fs->make<TH1I> ( "nHadNotFromTDecay", "Number of hadrons not from Top decay;N;Events", 10, 0, 10 );
    h_nJetNotFromTDecay = fs->make<TH1I> ( "nJetNotFromTDecay", "Number of jets not from Top decay;N;Events", 10, 0, 10 );
    h_nHadNotFromTHDecay = fs->make<TH1I> ( "nHadNotFromTHDecay", "Number of hadrons not from Top/Higgs decay;N;Events", 10, 0, 10 );
    h_hadNotFromTDecayPdg = fs->make<TH1I> ( "hadNotFromTDecayPdg", "PdgId of hadrons not from Top decay;PdgId;Hadrons", 52, -26, 26 );
    h_hadFromTDecayPdg = fs->make<TH1I> ( "hadFromTDecayPdg", "PdgId of hadrons from Top decay;PdgId;Hadrons", 52, -26, 26 );
    h_hadIsFromBHadron = fs->make<TH1I> ( "hadIsFromBHadron", "Whether hadron comes from a b-hadron;0-No 1-Yes;Hadrons", 3, 0, 3 );
    
    h_topDaughterQuarkFlavour = fs->make<TH1I> ( "topDaughterQuarkFlavour", "PdgId of the top daughter quark;PdgId;Events", 20, -10, 10 );
    
    h_hadAll_Pt = fs->make<TH1D> ( "hadAll_Pt", "Pt of all hadrons;Pt;Hadrons", 50, 0, 200 );
    h_hadNotClusteredInJet_Pt = fs->make<TH1D> ( "hadNotClusteredInJet_Pt", "Pt of hadrons not clustered to any jet;Pt;Hadrons", 50, 0, 200 );
    h_hadNoDecayClusteredInJet_Pt= fs->make<TH1D> ( "hadNoDecayClusteredInJet_Pt", "Pt of hadrons without decay products clustered to any jet;Pt;Hadrons", 50, 0, 200 );
    h_hadNoDecayButHadClusteredInJet_Pt= fs->make<TH1D> ( "hadNoDecayButHadClusteredInJet_Pt", "Pt of hadrons clustered to jets but decay products not clustered;Pt;Hadrons", 50, 0, 200 );
}

GenHFHadronAnalyzer::~GenHFHadronAnalyzer()
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
/**
* @brief description of the run-time parameters
*
* <TABLE>
* <TR><TH> name                </TH><TH> description </TH> </TR>
* <TR><TD> ttGenEvent          </TD><TD> input collection of TtGenEvent, used to identify the b-quark from top </TD></TR>
* <TR><TD> genJets             </TD><TD> input GenJet collection </TD></TR>
* <TR><TD> noBBbarResonances   </TD><TD> exclude resonances to be identified as hadrons </TD></TR>
* <TR><TD> onlyJetClusteredHadrons   </TD><TD> Whether only hadrons, injected to jets, shold be analyzed. Runs x1000 faster in Sherpa. Hadrons not clustered to jets will not be identified.
* <TR><TD> resolveParticleName </TD><TD> print particle name during warning and debug output instead of PDG ID </TD></TR>
* <TR><TD> flavour      </TD><TD> flavour of weakly decaying hadron that the jets should be matched to (5-b, 4-c) </TD></TR>
* </TABLE>
*
*/
// void GenHFHadronAnalyzer::fillDescriptions ( edm::ConfigurationDescriptions& descriptions )
// {

//     edm::ParameterSetDescription desc;
//     desc.add<edm::InputTag> ( "ttGenEvent",edm::InputTag ( "genEvt" ) )->setComment ( "Input collection of TtGenEvent, used to identify the quark from top" );
//     desc.add<edm::InputTag> ( "genJets",edm::InputTag ( "ak5GenJets","","HLT" ) )->setComment ( "Input GenJet collection" );
//     desc.add<bool> ( "noBBbarResonances",true )->setComment ( "Exclude resonances to be identified as hadrons" );
//     desc.add<bool> ( "onlyJetClusteredHadrons",false )->setComment ( "Whether only hadrons, injected to jets, shold be analyzed. Runs x1000 faster in Sherpa. Hadrons not clustered to jets will not be identified." );
//     desc.add<int> ( "flavour",5 )->setComment ( "Flavour of weakly decaying hadron that should be added to the jet for futher tagging. (4-c, 5-b)" );
//     descriptions.add ( "matchGenHFHadron",desc );
// }



//
// member functions
//

// ------------ method called to produce the data  ------------
void GenHFHadronAnalyzer::analyze ( const edm::Event& evt, const edm::EventSetup& setup )
{

    setup.getData ( pdt_ );

    using namespace edm;

//     printf("Run: %d\tLumi: %d\tEvent: %d\n",(int)evt.eventAuxiliary().run(), (int)evt.eventAuxiliary().luminosityBlock(), (int)evt.eventAuxiliary().event());

    edm::Handle<TtGenEvent> genEvt;
    evt.getByLabel ( ttGenEvent_, genEvt );

    edm::Handle<reco::GenJetCollection> genJets;
    evt.getByLabel ( genJets_, genJets );

// Hadron matching variables
    std::auto_ptr<std::vector<reco::GenParticle> > hadMothers ( new std::vector<reco::GenParticle> );
    std::auto_ptr<std::vector<std::vector<int> > > hadMothersIndices ( new std::vector<std::vector<int> > );
    std::auto_ptr<std::vector<int> > hadIndex ( new std::vector<int> );
    std::auto_ptr<std::vector<int> > hadFlavour ( new std::vector<int> );
    std::auto_ptr<std::vector<int> > hadJetIndex ( new std::vector<int> );
    std::auto_ptr<std::vector<int> > hadLeptonIndex ( new std::vector<int> );
    std::auto_ptr<std::vector<int> > hadLeptonHadIndex ( new std::vector<int> );
    std::auto_ptr<std::vector<int> > hadFromTopWeakDecay ( new std::vector<int> );
    std::auto_ptr<std::vector<int> > hadFromBHadron ( new std::vector<int> );

    LogDebug ( flavourStr_+"Jet (new)" ) << "searching for "<< flavourStr_ <<"-jets in " << genJets_;
    *hadJetIndex = findHadronJets ( *genJets, *hadIndex, *hadMothers, *hadMothersIndices, *hadLeptonIndex, *hadLeptonHadIndex, *hadFlavour, *hadFromTopWeakDecay, *hadFromBHadron );

}

// ------------ method called once each job just before starting event loop  ------------
void GenHFHadronAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void GenHFHadronAnalyzer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void GenHFHadronAnalyzer::beginRun ( edm::Run&, edm::EventSetup const& )
{
}

// ------------ method called when ending the processing of a run  ------------
void
GenHFHadronAnalyzer::endRun ( edm::Run&, edm::EventSetup const& )
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void GenHFHadronAnalyzer::beginLuminosityBlock ( edm::LuminosityBlock&, edm::EventSetup const& )
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void GenHFHadronAnalyzer::endLuminosityBlock ( edm::LuminosityBlock&, edm::EventSetup const& )
{
}

// ------------ helper functions -------------


/**
* @brief identify the jets that contain b-hadrons
*
* All jets originating from a b-hadron with the right b (c) content (b or anti-b) are identified in the GenJetCollection.
* The b (c) jet is identified by searching for a hadron of corresponding flavour in the jet. Hadron are put in jets
* by "TopAnalysis.TopUtils.GenJetParticles" plugin.
* For each hadron all mothers from all levels and chains are analyzed to find the quark or gluon from which the hadron has originated.
* This is done by searching through the generator particle decay tree starting from the hadron, performing checks for flavour and kinematic consistency.
* The hadrons that are not last in the decay chain (don't decay weakly) are skipped.
*
* b-bbar (c-cbar) resonances can either be considered as hadrons or not depending on the configuration.
*
*
* To allow performance studies following informations are tracked and returned:
* \arg the GenParticles identified as b-hadrons
* \arg if the b-hadron originates from the top-b-quark
* \arg a matrix showing which jets contain particles from which b-hadron (remnant from original approach in TopAnalysis)
*
* @param[in] genJets the GenJetCollection to be searched
* @param[out] hadIndex vector of indices of found hadrons in hadMothers
* @param[out] hadMothers vector of all mothers at all levels of each found hadron
* @param[out] hadMothersIndices connection between each particle from hadMothers and its mothers
* @param[out] hadLeptonIndex vector of indices representing leptons in hadMothers
* @param[out] hadLeptonHadIndex index of hadron associated to each lepton
* @param[out] hadFlavour flavour of each found hadron
* @param[out] hadFromTopWeakDecay flag showing whether the hadron comes from the products of top quark weak decay (works only with B-Hadrons)
* @returns vector of jets being matched to each hadron by jet clustering algorithm
*/
std::vector<int> GenHFHadronAnalyzer::findHadronJets ( const reco::GenJetCollection& genJets, std::vector<int> &hadIndex, 
                                                      std::vector<reco::GenParticle> &hadMothers, std::vector<std::vector<int> > &hadMothersIndices, 
                                                      std::vector<int> &hadLeptonIndex, std::vector<int> &hadLeptonHadIndex, std::vector<int> &hadFlavour, 
                                                      std::vector<int> &hadFromTopWeakDecay, std::vector<int> &hadFromBHadron )
{

    std::vector<int> result;
    std::vector<const reco::Candidate*> hadMothersCand;
    std::vector<int> hadFromJetWithClusteredHadrons;

    const unsigned int nJets = genJets.size();
    bool hadVsJet[100][400]= {{false}};				// Whether jet contains decay products of the hadron (max 100 hadrons, 250 jets foreseen)


    int topDaughterQId = -1;
    int topBarDaughterQId= -1;
    for ( size_t iJet = 0; iJet < nJets; ++iJet )  {
        const reco::GenJet* thisJet = & ( genJets[iJet] );
        std::vector<const reco::GenParticle*> particles = thisJet->getGenConstituents();
        
        bool hasClusteredHadron = false;
        // Skipping jets that don't have clustered hadrons
        for(const reco::GenParticle* particle : particles) {
            if(!isHadronPdgId(flavour_, particle->pdgId())) continue;
            hasClusteredHadron = true; 
            break;
        }
        //############################ WILL WORK X500 FASTER IN SHERPA
        //########## HADRONS NOT CLUSTERED TO ANY JET ARE LOST (~1-2%)
        if(onlyJetClusteredHadrons_) {
            if(!hasClusteredHadron) continue;
        }   // If jets without clustered hadrons should be skipped

    // printf(" genJet: %d pt: %.4f\teta: %.4f\tphi: %.4f\n",(int)iJet,thisJet->pt(),thisJet->eta(),thisJet->phi());
        for ( unsigned int iParticle = 0; iParticle < particles.size(); ++iParticle ) {
            const reco::GenParticle* thisParticle = particles[iParticle];
            const reco::Candidate* hadron = 0;
            const reco::Candidate* lepton = 0;
            
            if ( thisParticle->status() > 1 ) continue;    // Skipping non-final state particles (e.g. bHadrons)

//             printf("   particle: %d\tpdgId: %d\n", iParticle, thisParticle->pdgId());
            int hadronIndex = analyzeMothers ( thisParticle, &hadron, &lepton, topDaughterQId, topBarDaughterQId, hadMothersCand, hadMothersIndices, 0, -1 );           
//             printf("nHadMothers: %d  daughter1: %d  daughter2: %d\n", (int)hadMothersCand.size(), topDaughterQId, topBarDaughterQId);
            if ( hadron ) { // Putting hadron index to the list if it is not yet
                // Storing the index of the hadron and lepton
                int hadListIndex = isInList(hadIndex, hadronIndex);
                if ( hadListIndex<0 ) {
                    hadIndex.push_back ( hadronIndex );
                    hadListIndex = hadIndex.size()-1;
                    if(hasClusteredHadron) hadFromJetWithClusteredHadrons.push_back(1); else hadFromJetWithClusteredHadrons.push_back(0);
                }
                // Identifying the lepton in the hadron jet
                if ( lepton ) {
                    int leptonId = isInList(hadMothersCand, lepton);
                    if(isInList(hadLeptonIndex, leptonId)<0) {
                        hadLeptonIndex.push_back( leptonId );
                        hadLeptonHadIndex.push_back( hadListIndex );
                    }
                }
            }   // If hadron has been found in the chain
        }   // End of loop over jet consituents
    }   // End of loop over jets
    

    for ( int i=0; i< ( int ) hadMothersCand.size(); i++ ) {
        hadMothers.push_back ( ( *dynamic_cast<const reco::GenParticle*> ( hadMothersCand.at(i) ) ) );
    }

// Checking mothers of hadrons in order to assign flags (where the hadron comes from)
    unsigned int nHad = hadIndex.size();

// Checking number of different flavours of hadrons in the event
    int nFlavHadrons[5] = {0};		// Numbers of hadrons (0-Top, 1-Higgs, 2-Gluon, 3-Z)

    std::vector<std::vector<int> > LastQuarkMotherIds;
    std::vector<std::vector<int> > LastQuarkIds;
    std::vector<int> lastQuarkIndices(nHad, -1);
    

// Looping over all hadrons
    for ( unsigned int hadNum=0; hadNum<nHad; hadNum++ ) {
        int hadIdx = hadIndex.at(hadNum);   // Index of hadron in the hadMothers
        const reco::GenParticle* hadron = &hadMothers.at(hadIdx);
        int jetIndex = -1;
        
        // Checking daughters of the hadron
//         printf("Checking hadron %d out of %d\n", hadNum, nHad);
        for(int mothId=0; mothId<(int)hadMothersIndices.size(); ++mothId) {
//             printf(" Checking particle %d out of %d: nMoth: %d\n", mothId, (int)hadMothersIndices.size(), (int)hadMothersIndices.at(mothId).size());
//             printf("  Particle (%d)\n", hadMothers.at(mothId).pdgId());
            int daughterId = hadMothersIndices.at(mothId).at(0);
            if(daughterId != hadIdx) continue;
//             printf("  Daughter (%d) of hadron (%d)\n", hadMothers.at(mothId).pdgId(), hadMothers.at(hadIdx).pdgId());
        }

// Looping over all jets to match them to current hadron
        for ( unsigned int jetNum = 0; jetNum < nJets; jetNum++ ) {
// Checking whether jet contains this hadron in it (that was put in jet by clustering algorithm)
            std::vector<const reco::GenParticle*> particles = genJets.at(jetNum).getGenConstituents();
            for ( unsigned int partNum=0; partNum<particles.size(); partNum++ ) {
                const reco::GenParticle* particle = particles.at(partNum);
                if ( particle->status() <2 ) {
                    continue;    // Skipping final state particles
                }
                if ( !isHadron ( flavour_,particle ) ) {
                    continue;
                }
// Checking whether hadron and particle in jet are identical
                if ( hadron->pdgId() !=particle->pdgId() || fabs ( hadron->eta()-particle->eta() ) >0.00001 || fabs ( hadron->phi()-particle->phi() ) >0.00001 ) {
                    continue;
                }
                jetIndex=jetNum;
                break;
            }   // End of loop over jet constituents
            if ( jetIndex>=0 ) {
                break;
            }
        }   // End of loop over jets

        result.push_back ( jetIndex );  // Putting jet index to the result list

        std::vector <int> FirstQuarkId;
        std::vector <int> LastQuarkId;
        std::vector <int> LastQuarkMotherId;

        int hadFlav = hadMothers.at(hadIdx).pdgId() <0?-1:1; // Charge of the hadron (-1,1)
        if ( abs ( hadMothers.at(hadIdx).pdgId() ) /1000 < 1 ) {
            hadFlav*=-1;    // Inverting flavour of hadron if it is a meson
        }

        // Creating a list of pdgIds of neutral hadrons
        int bHadNeutralPdg_[] = {511, 513, 515, 531, 533, 535, 551, 553, 555, 5122, 5142, 5212,
          5214, 5232, 5322, 5324, 5342, 5412, 5414, 5432, 5434, 5522, 5524, 5542, 5544 };
        std::vector<int> bHadNeutralPdg (bHadNeutralPdg_, bHadNeutralPdg_ + sizeof(bHadNeutralPdg_) / sizeof(int) );
        bool verbose = false;

        // Checking whether this is a neutral hadron
        // if(std::find(bHadNeutralPdg.begin(), bHadNeutralPdg.end(), std::abs(hadMothers[hadIdx].pdgId()) ) != bHadNeutralPdg.end() ) {
        //   verbose = false;
        // }

        if(verbose) printf("\n\n\nChecking a hadron: %d ##########################################\n", hadMothers[hadIdx].pdgId());

// Searching only first quark in the chain with the same flavour as hadron
        findInMothers ( hadIdx, FirstQuarkId, hadMothersIndices, hadMothers, 0, hadFlav*flavour_, false, -1, 1, verbose );
        if(verbose) {
          printf("Found %d b quarks Hadron: %d ###########################################\n", (int)FirstQuarkId.size(), hadMothers[hadIdx].pdgId());
          for(int i=0; i<(int)FirstQuarkId.size(); i++) {
            printf("   %d. pdgId: %d\tstatus: %d", i, hadMothers[FirstQuarkId[i]].pdgId(), hadMothers[FirstQuarkId[i]].status());
            // if(hadFlav*hadMothers[FirstQuarkId[i]].pdgId() < 0) printf(" WRONG CHARGE");
            printf("\n");
          }
        }
//      if((int)FirstQuarkId.size()<1) findInMothers ( hadIdx, FirstQuarkId, hadMothersIndices, hadMothers, 0, hadFlav*flavour_, true, 1, true );
        // if((int)FirstQuarkId.size()<1) printf(" NO QUARK FOUND!!!\n");
        // return result;


        // Filling histograms for each 1-st quark
        for ( unsigned int qId=0; qId<FirstQuarkId.size(); qId++ ) {
            int qFlav = ( hadMothers[FirstQuarkId[qId]].pdgId() <0 ) ?-1:1;
            if ( toPlot ) {
                h_firstQstatus->Fill ( hadMothers[FirstQuarkId[qId]].status() );
                h_firstQPdg->Fill ( hadMothers[FirstQuarkId[qId]].pdgId() );
                h2_firstQFlav_HadFlav->Fill ( hadFlav,qFlav );
            }

            // Getting mothers of the first quark
            std::vector<int> FirstQMotherId = hadMothersIndices.at ( FirstQuarkId[qId] );
            for ( unsigned int qmId=0; qmId<FirstQMotherId.size(); qmId++ ) {
                int qmFlav = ( hadMothers.at(FirstQMotherId.at(qmId)).pdgId() <0 ) ?-1:1;
                if ( toPlot ) {
                    h2_firstQFlav_QMotherFlav->Fill ( qFlav,qmFlav );
                }
            }		// End of loop over all mothers of the first quark of the hadron

            // Finding last quark of the hadron starting from the first quark
//             printf("   Looking for last b quark --------------------------------\n");
            findInMothers ( FirstQuarkId.at(qId), LastQuarkId, hadMothersIndices, hadMothers, 0, qFlav*flavour_, false, -1, 2, false );

        }       // End of loop over all first quarks of the hadron
        // printf("   Found %d last b quarks +++++++++++++++++++++++++++++++\n", (int)LastQuarkId.size());
        // for(int i=0; i<(int)LastQuarkId.size(); i++) printf("      %d. qIdx: %d\tpdg: %d\tpt: %.3f\n", i, LastQuarkId.at(i), hadMothers.at(LastQuarkId.at(i)).pdgId(), hadMothers.at(LastQuarkId.at(i)).pt());
        
        
//         printf("First quarks: %d\n", (int)FirstQuarkId.size());
//         for(int qId : FirstQuarkId) printf("   %d.\tPdg: %d\tPt: %.3f\n", qId, hadMothers.at(qId).pdgId(), hadMothers.at(qId).pt());
//         printf("Last quarks: %d\n", (int)LastQuarkId.size());
//         for(int qId : LastQuarkId) {
//             printf("   %d.\tPdg: %d\tPt: %.3f\n", qId, hadMothers.at(qId).pdgId(), hadMothers.at(qId).pt());
//             int mother1Id = hadMothersIndices.at(qId).at(0);
//             if(mother1Id>=0) {
//                 const reco::GenParticle mother1 = hadMothers.at(mother1Id);
//                 printf("Mother 1: Pdg: %d\tPt: %.3f\n", mother1.pdgId(), mother1.pt());
//                 int mother2Id = hadMothersIndices.at(mother1Id).at(0);
//                 if(mother2Id>=0) {
//                     const reco::GenParticle mother2 = hadMothers.at(mother2Id);
//                     printf("Mother 2: Pdg: %d\tPt: %.3f\n", mother2.pdgId(), mother2.pt());
//                 } else printf("NO MOTHER 2");
//             } else printf("NO MOTHER 1");
//         }
//         if(hadFromTopWeakDecay.at(hadNum)) printf("From TOP\n"); else printf("NOT From TOP\n");
//         printf("\n\n");

        // Setting initial flavour of the hadron
        int hadronFlavour = 0;

        std::vector<std::pair<double, int> > lastQuark_dR_id_pairs;

        // Finding the closest quark in dR
        for ( unsigned int qId=0; qId<LastQuarkId.size(); qId++ ) {
            int qIdx = LastQuarkId.at(qId);
            int qFlav = ( hadMothers.at ( qIdx ).pdgId() <0 ) ?-1:1;

            // Checking the dR between hadron and quark
            float dR = deltaR ( hadMothers[hadIdx].eta(),hadMothers[hadIdx].phi(),hadMothers[qIdx].eta(),hadMothers[qIdx].phi() );
            // printf("  qIdx: %d dR: %.3f\n", qIdx, dR);
            std::pair<double, int> dR_hadId_pair(dR,qIdx);
            lastQuark_dR_id_pairs.push_back(dR_hadId_pair);

            if ( toPlot ) {
                h_lastQstatus->Fill ( hadMothers.at ( qIdx ).status() );
                h_lastQPdg->Fill ( hadMothers.at ( qIdx ).pdgId() );
                h2_lastQFlav_HadFlav->Fill ( hadFlav,qFlav );
            }
        }       // End of loop over all last quarks of the hadron

        std::sort(lastQuark_dR_id_pairs.begin(), lastQuark_dR_id_pairs.end());

        if(lastQuark_dR_id_pairs.size()>1) {
            double dRratio = (lastQuark_dR_id_pairs.at(1).first - lastQuark_dR_id_pairs.at(0).first)/lastQuark_dR_id_pairs.at(1).first;
            h_lastQ12_had_dRratio->Fill(dRratio);
            int qIdx_closest = lastQuark_dR_id_pairs.at(0).second;
            LastQuarkId.clear();
            if(dRratio>0.7) LastQuarkId.push_back(qIdx_closest); 
            else for(std::pair<double, int> qIdDrPair : lastQuark_dR_id_pairs) LastQuarkId.push_back(qIdDrPair.second);
        }
//         if(LastQuarkId.size()>0) {
        for(int qIdx : LastQuarkId) {
            int qmIdx = hadMothersIndices.at ( qIdx ).at(0);
            LastQuarkMotherId.push_back( qmIdx );
            // Filling the plots
            int qFlav = ( hadMothers.at ( qIdx ).pdgId() <0 ) ?-1:1;
            int qmFlav = ( hadMothers.at ( qmIdx ).pdgId() <0 ) ?-1:1;
            if ( toPlot ) {
                h2_lastQFlav_QMotherFlav->Fill ( qFlav,qmFlav );
                h2_lastQPdg_QMotherPdg->Fill ( hadMothers.at ( qIdx ).pdgId(),hadMothers.at ( qmIdx ).pdgId() );
            }
        }


        if((int)LastQuarkId.size()>0) lastQuarkIndices.at(hadNum) = 0;     // Setting the first quark in array as a candidate if it exists

        LastQuarkIds.push_back( LastQuarkId );

        LastQuarkMotherIds.push_back ( LastQuarkMotherId );

        // if ( toPlot ) {
        //     h_lastQ_lastQ_dRmin->Fill ( dRmin );
        // }

//         printf("   Ordered last b quarks, hadron: %d +++++++++++++++++++\n", hadNum);
//         for(int i=0; i<(int)LastQuarkId.size(); i++) printf("      %d. qIdx: %d\tpdg: %d\tpt: %.3f\tdR: %.3f\n", i, LastQuarkId.at(i), hadMothers.at(LastQuarkId.at(i)).pdgId(), hadMothers.at(LastQuarkId.at(i)).pt(), lastQuark_dR_id_pairs.at(i).first);
       
        if(LastQuarkMotherId.size()<1) {
            hadronFlavour = 0;
        } else {
            int qIdx = LastQuarkId.at( lastQuarkIndices.at(hadNum) );
            int qFlav = ( hadMothers.at(qIdx).pdgId() < 0 ) ? -1 : 1;
            hadronFlavour = qFlav*std::abs( hadMothers.at( LastQuarkMotherId.at( lastQuarkIndices.at(hadNum) ) ).pdgId() );
        }
        hadFlavour.push_back(hadronFlavour);    // Adding hadron flavour to the list of flavours
// 		printf(" 1. hadronFlavour: %d (%d)\n",hadFlavour.at(hadNum),hadFlav);
        
//         printf("nLastQuarks: %d  nHadMothers: %d  tQ: %d  tbarQId: %d\n", (int) LastQuarkId.size(), (int)hadMothers.size(), topDaughterQId, topBarDaughterQId);
        
        // Checking whether hadron comes from the Top weak decay
        int isFromTopWeakDecay = 1;
        std::vector <int> checkedParticles;
        if(hadFlavour.at(hadNum)!=0) {
            int lastQIndex = LastQuarkId.at(lastQuarkIndices.at(hadNum));
            bool fromTB = topDaughterQId>=0?findInMothers( lastQIndex, checkedParticles, hadMothersIndices, hadMothers, -1, 0, false, topDaughterQId, 2, false ) : false;
            checkedParticles.clear();
            bool fromTbarB = topBarDaughterQId>=0?findInMothers( lastQIndex, checkedParticles, hadMothersIndices, hadMothers, -1, 0, false, topBarDaughterQId, 2, false):false;
            if(!fromTB && !fromTbarB) {
                isFromTopWeakDecay = 0;
            }
        } else isFromTopWeakDecay = 2;
        hadFromTopWeakDecay.push_back(isFromTopWeakDecay);
        int isFromBHadron = findInMothers( hadIdx, checkedParticles, hadMothersIndices, hadMothers, 0, 555555, true, -1, 1, false )?1:0;
        h_hadIsFromBHadron->Fill((int)isFromBHadron);
        

        if(LastQuarkMotherId.size()>0) {
            std::set<int> checkedHadronIds;
//             int qIdx = LastQuarkId.at( lastQuarkIndices.at(hadNum) );
//             int qFlav = ( hadMothers.at(qIdx).pdgId() < 0 ) ? -1 : 1;
            fixExtraSameFlavours(hadNum, hadIndex, hadMothers, hadMothersIndices, hadFromTopWeakDecay, LastQuarkIds, LastQuarkMotherIds, lastQuarkIndices, hadFlavour, checkedHadronIds, 0);
        }
//         printf(" 2. hadronFlavour: %d (%d)\n",hadFlavour.at(hadNum),hadFlav);
        
        if(isFromTopWeakDecay==0) h_hadNotFromTDecayPdg->Fill(hadFlavour.at(hadNum)); else if(isFromTopWeakDecay==1) h_hadFromTDecayPdg->Fill(hadFlavour.at(hadNum));
//         if(std::abs(hadFlavour.at(hadNum))==6 && isFromTopWeakDecay==0) printf("ATTENTION: Top not from Top\n\n");


        int jetIdx1 = result[hadNum];   // Index of jet from matching by jet clustering algorithm

        int jetIdx2 = -1;
        float maxPt = 0.0;
        float dR_min=999.9;
// Looping over jets containing decay products of hadron to find the highest Pt one
        for ( unsigned int jetNum=0; jetNum<nJets; jetNum++ ) {
            if ( toPlot ) {
                float dR = deltaR ( genJets[jetNum].eta(),genJets[jetNum].phi(),hadron->eta(),hadron->phi() );
                if ( dR<dR_min ) {
                    dR_min=dR;
                }
            }
            if ( !hadVsJet[hadNum][jetNum] ) {
                continue;
            }
            if ( genJets[jetNum].pt() <=maxPt ) {
                continue;
            }
            maxPt=genJets[jetNum].pt();
            jetIdx2=jetNum;
        }
        

        if ( toPlot ) {
            h_hadJetIndex->Fill ( jetIdx1 );
            h_hadJet_dRmin->Fill ( dR_min );
            if ( jetIdx1>=0 ) {
                h_jetHad_dR->Fill ( deltaR ( genJets[jetIdx1].eta(),genJets[jetIdx1].phi(),hadron->eta(),hadron->phi() ) );
            }
            if ( jetIdx2>=0 && jetIdx1>=0 ) {
                h_jetdPt->Fill ( genJets[jetIdx1].pt() /genJets[jetIdx2].pt() );
                h_jetDeltaPt->Fill ( genJets[jetIdx1].pt()-genJets[jetIdx2].pt() );
            }
            if ( jetIdx2>=0 ) {
                h_jetHadOld_dR->Fill ( deltaR ( genJets[jetIdx2].eta(),genJets[jetIdx2].phi(),hadron->eta(),hadron->phi() ) );
            }
// Setting index to highest Pt jet containing decay products of the hadron if no jet was associated by jet clustering algorithm
            if ( jetIdx2>=0 && jetIdx1<0 ) {
                h_jetHadOldNoMatch_dR->Fill ( deltaR ( genJets[jetIdx2].eta(),genJets[jetIdx2].phi(),hadron->eta(),hadron->phi() ) );
            }
            if ( jetIdx2>=0 && jetIdx1<0 && abs ( hadronFlavour ) ==0 ) {
                h_jetHadOldNoMatch_dR_noFlav->Fill ( deltaR ( genJets[jetIdx2].eta(),genJets[jetIdx2].phi(),hadron->eta(),hadron->phi() ) );
            }
            if ( abs ( hadronFlavour ) ==0 && jetIdx1>=0 ) {
                h_jetHad_dR_noFlav->Fill ( deltaR ( genJets[jetIdx1].eta(),genJets[jetIdx1].phi(),hadron->eta(),hadron->phi() ) );
            }

            h_nFirstQ->Fill ( FirstQuarkId.size() );
            h_nLastQ->Fill ( LastQuarkId.size() );
        }
    }		// End of loop over all hadrons

    int nHadNotFromTDecay = 0;
    int nHadNotFromTHDecay = 0;
    std::vector<int> hadNotFromTJetId;
        

// Checking number leptons and different flavours of hadrons in the event
    for ( unsigned int i=0; i<hadFlavour.size(); i++ ) {
        // Checking N leptons
        int nLeps = 0;
        for ( unsigned int lepHadId : hadLeptonHadIndex ) {
            if(i==lepHadId) nLeps++;
        }
        h_nLeptonsInHad->Fill( nLeps );
        
        h_hadAll_Pt->Fill(hadMothers.at(i).pt());
        if(result.at(i)<0 ) {
            h_hadNotClusteredInJet_Pt->Fill(hadMothers.at(i).pt());
            if(hadFromJetWithClusteredHadrons.at(i)==0) h_hadNoDecayClusteredInJet_Pt->Fill(hadMothers.at(i).pt());
        } else if(result.at(i)>=0 && hadFromJetWithClusteredHadrons.at(i)==0) h_hadNoDecayButHadClusteredInJet_Pt->Fill(hadMothers.at(i).pt());
        
        // Checking whether the hadron is from Top/Higgs
        if(hadFromTopWeakDecay.at(i)==0 && std::abs(hadFlavour.at(i))!=6) {
            nHadNotFromTDecay++;
            if(std::abs(hadFlavour.at( i )) != 25) nHadNotFromTHDecay++;
            if(std::find(hadNotFromTJetId.begin(), hadNotFromTJetId.end(), result.at(i)) == hadNotFromTJetId.end()) hadNotFromTJetId.push_back(result.at(i));
        }
        
        // Checking flavours
        h_hadFlavour->Fill( hadFlavour.at(i) );
        // printf("HADRON %d : FLAVOUR %d\n", i, hadFlavour.at(i));
        switch ( abs ( hadFlavour.at ( i ) ) ) {
        case 6:
            nFlavHadrons[0]++;	// Top
            break;
        case 25:
            nFlavHadrons[1]++;	// Higgs
            break;
        case 21:
            nFlavHadrons[2]++;	// Gluon
            break;
        case 23:
            nFlavHadrons[3]++;	// Z
            break;
        default:
            nFlavHadrons[4]++;	// No
            break;
        }
    }
    
    
    h_nJetNotFromTDecay->Fill((int)hadNotFromTJetId.size());

    // if( nFlavHadrons[0] > 3 ) printf("nHadTop > 3\n\n");
    // if( nFlavHadrons[1] < 2 ) printf("nHadHiggs < 2**********************\n\n");

    if ( toPlot ) {
        h_nHad->Fill ( nHad );
        h_nHadTop->Fill ( nFlavHadrons[0] );
        h_nHadHiggs->Fill ( nFlavHadrons[1] );
        h_nHadGluon->Fill ( nFlavHadrons[2] );
        h_nHadZ->Fill ( nFlavHadrons[3] );
        h_nHadNo->Fill ( nFlavHadrons[4] );
        h_nHadNotFromTDecay->Fill( nHadNotFromTDecay );
        h_nHadNotFromTHDecay->Fill( nHadNotFromTHDecay );
        
        if(topDaughterQId>=0) h_topDaughterQuarkFlavour->Fill(hadMothers.at(topDaughterQId).pdgId()); else h_topDaughterQuarkFlavour->Fill(0);
        if(topBarDaughterQId>=0) h_topDaughterQuarkFlavour->Fill(hadMothers.at(topBarDaughterQId).pdgId()); else h_topDaughterQuarkFlavour->Fill(0);

// Filling number of hadrons associated to a single jet
        for ( unsigned int iJet=0; iJet<nJets; iJet++ ) {
            int N = 0;
            for ( unsigned int i=0; i<result.size(); i++ ) if ( result[i]== ( int ) iJet ) {
                    N++;
                }
            h_nHadInJet->Fill ( N );
        }

    }		// If validation histograms should be filled

    return result;
}


/**
* @brief Check if the cpecified particle is already in the list of particles
*
* @param[in] particleList - list of particles to be checked
* @param[in] particle - particle that should be checked
*
* returns true if a particle is already in the list
*/
int GenHFHadronAnalyzer::isInList ( std::vector<const reco::Candidate*> particleList, const reco::Candidate* particle )
{
    for ( unsigned int i = 0; i<particleList.size(); i++ )
        if ( particleList[i]==particle ) {
            return i;
        }

    return -1;
}

int GenHFHadronAnalyzer::isInList ( std::vector<int> list, const int value )
{
    for ( unsigned int i = 0; i<list.size(); i++ )
        if ( list.at(i)==value ) {
            return i;
        }

    return -1;
}


/**
* @brief Check the pdgId of a given particle if it is a hadron
*
* @param[in] flavour - flavour of a hadron that is being searched (5-B, 4-C)
* @param[in] thisParticle - a particle that is to be analysed
*
* @returns if the particle is a hadron of specified flavour
*/
bool GenHFHadronAnalyzer::isHadron ( const int flavour, const reco::Candidate* thisParticle )
{
    return isHadronPdgId(flavour, thisParticle->pdgId());
}

/**
* @brief Check the pdgId if it represents a hadron of particular flavour
*
* @param[in] flavour - flavour of a hadron that is being searched (5-B, 4-C)
* @param[in] pdgId - pdgId to be checked
*
* @returns if the pdgId represents a hadron of specified flavour
*/
bool GenHFHadronAnalyzer::isHadronPdgId ( const int flavour, const int pdgId )
{
    int flavour_abs = std::abs(flavour);
    if(flavour_abs > 5 || flavour_abs < 1) return false;
    int pdgId_abs = std::abs(pdgId);

    if ( pdgId_abs / 1000 == flavour_abs // baryons
            || ( pdgId_abs / 100 % 10 == flavour_abs // mesons
                 && ! ( noBBbarResonances_ && pdgId_abs / 10 % 100 == 11*flavour_abs ) // but not a resonance
               )
       ) {
        return true;
    } else {
        return false;
    }
}

/**
* @brief Check if the particle has bHadron among daughters
*
* @param[in] flavour - flavour of a hadron that is being searched (5-B, 4-C)
* @param[in] thisParticle - a particle that is to be analysed
*
* @returns if the particle has a hadron among its daughters
*/
bool GenHFHadronAnalyzer::hasHadronDaughter ( const int flavour, const reco::Candidate* thisParticle )
{
// Looping through daughters of the particle
    bool hasDaughter = false;
    for ( int k=0; k< ( int ) thisParticle->numberOfDaughters(); k++ ) {
        if ( !isHadron ( flavour, thisParticle->daughter ( k ) ) ) {
            continue;
        }
        hasDaughter = true;
        break;
    }
    return hasDaughter;
}

/**
* @brief do a recursive search for the mother particles until the b-quark is found or the absolute mother is found
*
* the treatment of b-bar resonances depends on the global parameter noBBbarResonances_
*
* @param[in] thisParticle current particle from which starts the search of the hadron and all its mothers up to proton
* @param[out] hadron the last hadron in the decay chain (that decays weekly)
* @param[out] lepton lepton found in the current decay chain
* @param[out] top top quark found in the decay chain
* @param[out] hadMothers list of all particles starting with hadron and ending with proton
* @param[out] hadMothersIndices list of i-vectors containing j-indices representing particles that are mothers of each i-particle from hadMothers
* @param[out] analyzedParticles list of particles analysed in the chain (loop detection)
* @param[out] prevPartIndex index of the previous particle in the current chain (loop detection)
*
* @returns index of the hadron in the hadMothers list. -1 if no hadron found
*/
int GenHFHadronAnalyzer::analyzeMothers ( const reco::Candidate* thisParticle, pCRC *hadron, pCRC *lepton, int& topDaughterQId, int& topBarDaughterQId,
                                          std::vector<const reco::Candidate*> &hadMothers, std::vector<std::vector<int> > &hadMothersIndices, 
                                          std::set<const reco::Candidate*> *analyzedParticles, const int prevPartIndex )
{

    int hadronIndex=-1;	// Index of the hadron that is returned by this function
// Storing the first hadron has been found in the chain when going up from the final particle of the jet
    if ( *hadron == 0 // find only the first b-hadron on the way (the one that decays weekly)
            && isHadron ( flavour_, thisParticle ) // is a hadron
            && !hasHadronDaughter ( flavour_, thisParticle ) // has no hadron daughter (decays weekly)
       ) {
        *hadron = thisParticle;

        int index = isInList ( hadMothers, thisParticle );

        if ( index<0 ) { // If hadron is not in the list of mothers yet
            hadMothers.push_back ( thisParticle );
            hadronIndex=hadMothers.size()-1;
        } else {	    // If hadron is in the list of mothers already
            hadronIndex=index;
        }
    }
    // Checking if the particle is a lepton
    if(!*hadron){
        int absPdg = std::abs(thisParticle->pdgId());
        if(absPdg==11 || absPdg==13) {
            const reco::Candidate* mother1 = 0;
            if(thisParticle->numberOfMothers()>0) mother1 = thisParticle->mother(0);
            if(mother1 && isHadron(flavour_, mother1) ) {
                *lepton = thisParticle;
                if(isInList(hadMothers, thisParticle)<0) hadMothers.push_back(thisParticle);
            }   // If the lepton's mother is a hadron
            else if(mother1 && std::abs(mother1->pdgId()) == 15) {
                const reco::Candidate* mother2 = 0;
                if(mother1->numberOfMothers()>0) mother2 = mother1->mother(0);
                if(mother2 && isHadron(flavour_, mother2)) {
                    *lepton = thisParticle;
                    if(isInList(hadMothers, thisParticle)<0) hadMothers.push_back(thisParticle);
                }   // If the tau's mother is a hadron
            }   // If the lepton's mother is a tau lepton
        }   // If this is a lepton
    }   // If no hadron found yet


    //################################################### FOR SHERPA
    // if(hadronIndex > -1) return hadronIndex;
    //################################################### FOR SHERPA

    int partIndex = -1;   // Index of particle being checked in the list of mothers
    partIndex = isInList ( hadMothers, thisParticle );
    // printf("      hadronIndex: %d\tpartIndex: %d|%d (%d)\tpdgId: %d\n", hadronIndex, partIndex, prevPartIndex, (int)hadMothers.size(), thisParticle->pdgId());

// Checking whether this particle is already in the chain of analyzed particles in order to identify a loop
    bool isLoop = false;
    if ( !analyzedParticles ) {
        analyzedParticles = new std::set<const reco::Candidate*>;
    }
    for ( unsigned int i=0; i<analyzedParticles->size(); i++ ) {
        if ( analyzedParticles->count ( thisParticle ) <=0 ) {
            continue;
        }
        isLoop = true;
        break;
    }
    

// If a loop has been detected
    if ( isLoop ) {
        if ( prevPartIndex>=0 ) {
            putMotherIndex ( hadMothersIndices, prevPartIndex, -1 );    // Setting mother index of previous particle to -1
        }
        return hadronIndex;		// Stopping further processing of the current chain
    }
    analyzedParticles->insert ( thisParticle );


    

// Putting the mothers to the list of mothers
    for ( size_t iMother = 0; iMother < thisParticle->numberOfMothers(); ++iMother ) {
        const reco::Candidate* mother = thisParticle->mother ( iMother );
        int mothIndex = isInList ( hadMothers, mother );
        if ( mothIndex == partIndex && partIndex>=0 ) continue;		// Skipping the mother that is its own daughter

        int nSameDaugh=0;

// Checking whether this mother has other daughters that have the same pdgId [Only if plotting option is enabled]
        if ( toPlot && abs ( thisParticle->pdgId() ) ==flavour_ && abs ( thisParticle->pdgId() ) == abs ( mother->pdgId() ) ) {
            double highPt = -1.0;
            int highPtId = -1;
// Loop over daughters of the mother to look how many daughters with the same pdgId it has
            float dR_min=999.9;
            for ( unsigned int iDaughter = 0; iDaughter<mother->numberOfDaughters(); iDaughter++ ) {
                const reco::Candidate* daughter = mother->daughter ( iDaughter );
                if ( abs ( daughter->pdgId() ) !=abs ( thisParticle->pdgId() ) ) {
                    continue;
                }
                nSameDaugh++;
            }
            for ( unsigned int iDaughter = 0; iDaughter<mother->numberOfDaughters(); iDaughter++ ) {
                const reco::Candidate* daughter = mother->daughter ( iDaughter );
                if ( abs ( daughter->pdgId() ) !=abs ( thisParticle->pdgId() ) ) {
                    continue;
                }
                if ( daughter->pdgId() *mother->pdgId() < 0 ) {
                    continue;
                }
                if ( daughter->pt() >highPt ) {
                    highPt=daughter->pt();
                    highPtId=iDaughter;
                }
                float dR=deltaR ( daughter->eta(),daughter->phi(),mother->eta(),mother->phi() );

                h_QQMothEFrac->Fill ( daughter->energy() /mother->energy() );
                h_QQMothdR->Fill ( dR );
                if ( nSameDaugh>1 ) {
                    h_QQMothEFrac_ambig->Fill ( daughter->energy() /mother->energy() );
                    h_QQMothdR_ambig->Fill ( dR );
                }

                if ( dR<dR_min ) {
                    dR_min=dR;
                }
            }

            if ( highPt<0.0 ) {
                continue;
            }
            for ( unsigned int iDaughter = 0; iDaughter<mother->numberOfDaughters(); iDaughter++ ) {
                const reco::Candidate* daughter = mother->daughter ( iDaughter );
                if ( daughter->pdgId() *mother->pdgId() < 0 ) {
                    continue;    // Skipping particles with wrong charge
                }
                if ( ( int ) iDaughter==highPtId ) {
                    if ( toPlot ) {
                        h_QHQdR->Fill ( deltaR ( daughter->eta(),daughter->phi(),mother->eta(),mother->phi() ) );
                        if ( nSameDaugh>1 ) {
                            h_QHQdR_ambig->Fill ( deltaR ( daughter->eta(),daughter->phi(),mother->eta(),mother->phi() ) );
                        }
                    }
                    continue;		// Skipping the highest Pt particle
                }
                h_QHQPtFrac->Fill ( ( highPt-daughter->pt() ) /highPt );
                if ( nSameDaugh>1 ) {
                    h_QHQPtFrac_ambig->Fill ( ( highPt-daughter->pt() ) /highPt );
                }
            }
            h_QQMothdR_min->Fill ( dR_min );
            if ( nSameDaugh>1 ) {
                h_QQMothdR_min_ambig->Fill ( dR_min );
            }
            h2_QQMothSamePdg->Fill ( thisParticle->pdgId(),mother->pdgId() );
        }	// If both particle and mother are quarks and have same abs(pdgId)

// If this mother isn't yet in the list and hadron or lepton is in the list
        if ( mothIndex<0 && ((*hadron) !=0 || (*lepton)) ) {
            hadMothers.push_back ( mother );
            mothIndex=hadMothers.size()-1;
        }
// If hadron has already been found in current chain and the mother isn't a duplicate of the particle being checked
        if ( (*hadron || *lepton) && mothIndex!=partIndex && partIndex>=0 ) {
//             printf(" Setting mother (%d) for particle (%d)\n", mother->pdgId(), thisParticle->pdgId());
            putMotherIndex ( hadMothersIndices, partIndex, mothIndex );			// Putting the index of mother for current particle
        }
        int index = analyzeMothers ( mother, hadron, lepton, topDaughterQId, topBarDaughterQId, hadMothers, hadMothersIndices, analyzedParticles, partIndex );
        hadronIndex = index<0?hadronIndex:index;
        
        // Setting the id of the particle which is a quark from the top decay
        if(*hadron && std::abs(mother->pdgId())==6) {
//             printf("pdg: %d  pt: %.3f\n", thisParticle->pdgId(), thisParticle->pt());
            int& bId = mother->pdgId() < 0 ? topBarDaughterQId : topDaughterQId;
            int thisFlav = std::abs(thisParticle->pdgId());
            if( bId<0){
                if(thisFlav <= 5) bId = partIndex; 
            } else {
                int bIdFlav = std::abs(hadMothers.at(bId)->pdgId());
                if( bIdFlav != 5 && thisFlav == 5) bId = partIndex;
                else if( thisFlav == 5 && thisParticle->pt() > hadMothers.at(bId)->pt() ) bId = partIndex;
            }
        }
    }			// End of loop over mothers

    analyzedParticles->erase ( thisParticle );		// Removing current particle from the current chain that is being analyzed

    if ( partIndex<0 ) {
        return hadronIndex;    // Safety check
    }

// Adding -1 to the list of mother indices for current particle if it has no mothers (for consistency between numbering of indices and mothers)
    if ( ( int ) thisParticle->numberOfMothers() <=0 && ( *hadron ) !=0 ) {
        putMotherIndex ( hadMothersIndices, partIndex, -1 );
    }


    return hadronIndex;

}

/**
* @brief puts mother index to the list of mothers of particle, if it isn't there already
*
* @param[in] hadMothersIndices vector of indices of mothers for each particle
* @param[in] partIndex index of the particle for which the mother index should be stored
* @param[in] mothIndex index of mother that should be stored for current particle
*
*/
bool GenHFHadronAnalyzer::putMotherIndex ( std::vector<std::vector<int> > &hadMothersIndices, int partIndex, int mothIndex )
{
    // Putting vector of mothers indices for the given particle
    bool inList=false;
    if ( partIndex<0 ) {
        return false;
    }

    while ( ( int ) hadMothersIndices.size() <=partIndex ) { // If there is no list of mothers for current particle yet
        std::vector<int> mothersIndices;
        hadMothersIndices.push_back ( mothersIndices );
    }

    std::vector<int> *hadMotherIndices=&hadMothersIndices.at ( partIndex );
// Removing other mothers if particle must have no mothers
    if ( mothIndex==-1 ) {
        hadMotherIndices->clear();
    } else {
// Checking if current mother is already in the list of theParticle's mothers
        for ( int k=0; k< ( int ) hadMotherIndices->size(); k++ ) {
            if ( hadMotherIndices->at ( k ) !=mothIndex && hadMotherIndices->at ( k ) !=-1 ) {
                continue;
            }
            inList=true;
            break;
        }
    }
// Adding current mother to the list of mothers of this particle
    if ( !inList ) {
        hadMotherIndices->push_back ( mothIndex );
    }

    return inList;
}


/**
* @brief helper function to find indices of particles with particular pdgId and status from the list of mothers of a given particle
*
* @param[in] idx index of particle, mothers of which should be searched
* @param[in] mothChains vector of indices where the found mothers are put
* @param[in] hadMothersIndices list of indices pointing to mothers of each particle from list of mothers
* @param[in] hadMothers vector of all hadron mother particles of all levels
* @param[in] status status of mother that is being looked for
* @param[in] pdgId pdgId of mother that is being looked for
* @param[in] pdgAbs option, whether the sign of pdgId matters
* @param[in] stopId id of the particle after which the checking should stop
* @param[in] firstLast should all(0), first(1) or last(2) occurances of the searched particle be stored
* @param[in] verbose option to print every particle that is processed during the search
*
*/

bool GenHFHadronAnalyzer::findInMothers ( int idx, std::vector<int> &mothChains, std::vector<std::vector<int> > &hadMothersIndices, std::vector<reco::GenParticle> &hadMothers, int status, int pdgId, bool pdgAbs=false, int stopId=-1, int firstLast=0, bool verbose=false	)
{
    bool foundStopId = false;
    int pdg_1 = hadMothers.at ( idx ).pdgId();
    int partCharge = ( hadMothers.at ( idx ).pdgId() >0 ) ?1:-1;
// Inverting charge if mother is a b(c) meson
    if ( abs ( hadMothers.at ( idx ).pdgId() ) /1000 < 1 && ( abs ( hadMothers.at ( idx ).pdgId() ) /100%10 == 4 || abs ( hadMothers.at ( idx ).pdgId() ) /100%10 == 5 ) ) {
        partCharge*=-1;
    }

    if ( ( int ) hadMothersIndices.size() <=idx ) {
        if ( verbose ) {
            printf ( " Stopping checking particle %d. No mothers are stored.\n",idx );
        }
        return false;	  // Skipping if no mothers are stored for this particle
    }

    std::vector<int> mothers = hadMothersIndices.at ( idx );
    unsigned int nMothers = mothers.size();
    bool isCorrect=false;		// Whether current particle is what is being searched
    if ( verbose ) {
        if ( abs ( hadMothers.at ( idx ).pdgId() ) ==2212 ) {
            printf ( "Chk:  %d\tpdg: %d\tstatus: %d",idx, hadMothers.at ( idx ).pdgId(), hadMothers.at ( idx ).status() );
        } else {
            printf ( " Chk:  %d(%d mothers)\tpdg: %d\tstatus: %d\tPt: %.3f\tEta: %.3f",idx, nMothers, hadMothers.at ( idx ).pdgId(), hadMothers.at ( idx ).status(), hadMothers.at ( idx ).pt(),hadMothers.at ( idx ).eta() );
        }
    }
    bool hasCorrectMothers = true;
    if(nMothers<1) hasCorrectMothers=false; else if(mothers.at(0)<0) hasCorrectMothers=false;
    // for ( unsigned int i=0; i<nMothers; i++ ) {
    //     int motherId = mothers.at ( i );
    //     if ( motherId<0 ) {
    //         continue;
    //     }
    //     if ( !isNeutralPdg ( hadMothers.at ( motherId ).pdgId() ) && hadMothers.at ( motherId ).pdgId() * hadMothers.at ( idx ).pdgId() < 0 && !pdgAbs && !isHadronPdgId(pdgId, hadMothers.at( motherId ).pdgId() ) ) {
    //         continue;    // Skipping if mother of the particle has wrong charge
    //     }
    //     hasCorrectMothers=true;
    //     break;
    // }
    if(!hasCorrectMothers) {
        if(verbose) printf("    NO CORRECT MOTHER\n");
        return false;
    }
    // Stopping if the particular particle has been found
    if(stopId>=0 && idx == stopId) return true;
    
    // Stopping if the hadron of particular flavour has been found
    if(pdgId%111111==0 && pdgId!=0) {
        if(isHadronPdgId(pdgId/111111, hadMothers.at(idx).pdgId())) {
            return true;
        }
    }
    
// Checking whether current mother satisfies selection criteria
    if ( ( ( hadMothers.at ( idx ).pdgId() == pdgId && pdgAbs==false )
            || ( abs ( hadMothers.at ( idx ).pdgId() ) == abs ( pdgId ) && pdgAbs==true ) )
            && ( hadMothers.at ( idx ).status() == status || status==0 )
            && hasCorrectMothers ) {
        isCorrect=true;
        bool inList=false;
        for ( unsigned int k=0; k<mothChains.size(); k++ ) if ( mothChains[k]==idx ) {
                inList=true;    // Checking whether isn't already in the list
                break;
            }
        if ( !inList && mothers.at ( 0 ) >=0 && ( hadMothers.at ( idx ).pdgId() *pdgId>0 || !pdgAbs ) ) {		// If not in list and mother of this quark has correct charge
            if ( firstLast==0 || firstLast==1 ) {
                mothChains.push_back ( idx );
            }
            if ( verbose ) {
                printf ( "   *" );
            }
        }
        if ( verbose ) {
            printf ( "   +++" );
        }
    }
    if ( verbose ) {
        printf ( "\n" );
    }

    if ( isCorrect && firstLast==1 ) {
        return false;   // Stopping if only the first particle in the chain is looked for
    }

// Checking next level mothers
    for ( unsigned int i=0; i<nMothers; i++ ) {
        int idx2 = mothers[i];
        if ( idx2<0 ) {
            continue;    // Skipping if mother's id is -1 (no mother), that means current particle is a proton
        }
        if ( idx2==idx ) {
            continue;    // Skipping if particle is stored as its own mother
        }
        int pdg_2 = hadMothers[idx2].pdgId();
        // Inverting the flavour if bb oscillation identified
        if ( isHadronPdgId(pdgId, pdg_1) && isHadronPdgId(pdgId, pdg_2) &&  pdg_1*pdg_2 < 0 ) {
            pdgId*=-1;
            if(verbose) printf("######### Inverting flavour of the hadron\n");
        }
        if ( firstLast==2 && isCorrect && (
                    ( abs ( hadMothers[idx2].pdgId() ) != abs ( pdgId ) && pdgAbs==true ) ||
                    ( hadMothers[idx2].pdgId() != pdgId && pdgAbs==false ) ) ) {		// If only last occurance must be stored and mother has different flavour
            if ( verbose ) {
                printf ( "Checking mother %d out of %d mothers once more to store it as the last quark\n",i,nMothers );
            }
            foundStopId = findInMothers ( idx, mothChains, hadMothersIndices, hadMothers, 0, pdgId, pdgAbs, stopId, 1, verbose );
        }

// Checking next level mother
        if ( verbose ) {
            printf ( "Checking mother %d out of %d mothers with pdgId: %d\n",i,nMothers,pdgId );
        }
        if(firstLast==2 && pdg_1 != pdg_2) continue;    // Requiring the same flavour when finding the last quark
        foundStopId = findInMothers ( idx2, mothChains, hadMothersIndices, hadMothers, status, pdgId, pdgAbs, stopId, firstLast, verbose );
    }
    
    return foundStopId;
}


/**
* @brief Check whether a given pdgId represents neutral particle
*
* @param[in] pdgId
* @param[in] thisParticle - a particle that is to be analysed
*
* @returns if the particle has a hadron among its daughters
*/
bool GenHFHadronAnalyzer::isNeutralPdg ( int pdgId )
{
    const int max = 5;
    int neutralPdgs[max]= {9,21,22,23,25};
    for ( int i=0; i<max; i++ ) if ( abs ( pdgId ) ==neutralPdgs[i] ) {
            return true;
        }
    return false;
}


/**
 * Finds hadrons that have the same flavour and origin and resolve this ambiguity
 * @method fixExtraSameFlavours
 * @param  hadId                Index of the hadron being checked
 * @param  hadIndex             Vector of indices of each hadron pointing to the hadMothers
 * @param  hadMothers           Vector of gen particles (contain hadrons and all its ancestors)
 * @param  hadMothersIndices    Vector of indices for each particle from hadMothers
 * @param  LastQuarkIds         Vector of indices of last quarks for each hadron
 * @param  LastQuarkMotherIds   Vector of indices of mothers for each last quark from LastQuarkIds
 * @param  lastQuakIndices      Vector of indices pointing to particular index from the LastQuarkIds and LastQuarkMotherIds to be used for each hadron
 * @param  lastQuarkIndex       Index from the LastQuarkIds and LastQuarkMotherIds for this particular hadron with index hadId
 * @return Whether other mother with unique association has been found for the hadrons
 */
bool GenHFHadronAnalyzer::fixExtraSameFlavours(
    const unsigned int hadId, const std::vector<int> &hadIndices, const std::vector<reco::GenParticle> &hadMothers, 
    const std::vector<std::vector<int> > &hadMothersIndices, const std::vector<int> &isFromTopWeakDecay, 
    const std::vector<std::vector<int> > &LastQuarkIds, const std::vector<std::vector<int> > &LastQuarkMotherIds, 
    std::vector<int> &lastQuarkIndices, std::vector<int> &hadronFlavour, std::set<int> &checkedHadronIds, const int lastQuarkIndex)
{
    if(checkedHadronIds.count(hadId) != 0) return false;      // Hadron already checked previously and should be skipped
    checkedHadronIds.insert(hadId);                           // Putting hadron to the list of checked ones in this run

    if(lastQuarkIndex<0) return false;
    if((int)LastQuarkIds.at(hadId).size()<lastQuarkIndex+1) return false;
    int LastQuarkId = LastQuarkIds.at(hadId).at(lastQuarkIndex);
    int LastQuarkMotherId = LastQuarkMotherIds.at( hadId ).at( lastQuarkIndex );
    int qmFlav = hadMothers.at(LastQuarkId).pdgId() < 0 ? -1 : 1;
    int hadFlavour = qmFlav*std::abs( hadMothers.at( LastQuarkMotherId ).pdgId() );
//     int hadFlavour = hadronFlavour.at(hadId);
    bool ambiguityResolved = true;
    // If last quark has inconsistent flavour with its mother, setting the hadron flavour to gluon
    if( (hadMothers.at(LastQuarkId).pdgId()*hadMothers.at(LastQuarkMotherId).pdgId() < 0 && !isNeutralPdg(hadMothers.at(LastQuarkMotherId).pdgId())) || 
        // If particle not coming from the Top weak decay has Top flavour
        (std::abs(hadronFlavour.at(hadId))==6 && isFromTopWeakDecay.at(hadId)==0) ) {
        if((int)LastQuarkIds.at(hadId).size()>lastQuarkIndex+1) fixExtraSameFlavours(hadId, hadIndices, hadMothers, hadMothersIndices, isFromTopWeakDecay, LastQuarkIds, LastQuarkMotherIds, lastQuarkIndices, hadronFlavour, checkedHadronIds, lastQuarkIndex+1); 
        else hadronFlavour.at(hadId) = qmFlav*21;
//         printf("Setting flavour of hadron %d to %d (wrong mother flavour)\n", hadId, hadronFlavour.at(hadId));
        return true;
    }

//     printf("    1. hadId: %d  LastQuarkId: %d hadFlavour: %d\n", hadId, LastQuarkId, hadFlavour);
    int nSameFlavourHadrons = 0;
    // Looping over all previous hadrons
    for(unsigned int iHad = 0; iHad<hadronFlavour.size(); iHad++) {
        if(iHad==hadId) continue;
        int theLastQuarkIndex = lastQuarkIndices.at(iHad);
        if(theLastQuarkIndex<0) continue;
        if((int)LastQuarkMotherIds.at( iHad ).size() <= theLastQuarkIndex) continue;
        int theLastQuarkMotherId = LastQuarkMotherIds.at( iHad ).at( theLastQuarkIndex );
        int theHadFlavour = hadronFlavour.at(iHad);
        // Skipping hadrons with different flavour
        if(theHadFlavour==0 || std::abs(theHadFlavour)==21) continue;
//         printf("    2. hadId: %d flavour: %d|%d motherId: %d(%.3f)|%d(%.3f)\n", iHad, hadFlavour, theHadFlavour, LastQuarkMotherId, hadMothers.at(LastQuarkMotherId).pt(), theLastQuarkMotherId, hadMothers.at(theLastQuarkMotherId).pt());
        if(theHadFlavour != hadFlavour || theLastQuarkMotherId != LastQuarkMotherId) continue;
        ambiguityResolved = false;
        nSameFlavourHadrons++;
        
        // Checking other b-quark mother candidates of this hadron
        if((int)LastQuarkIds.at(hadId).size() > lastQuarkIndex+1) {
            if(fixExtraSameFlavours(hadId, hadIndices, hadMothers, hadMothersIndices, isFromTopWeakDecay, LastQuarkIds, LastQuarkMotherIds, lastQuarkIndices, hadronFlavour, checkedHadronIds, lastQuarkIndex+1) ) {
                ambiguityResolved = true;
                break;
            }
        } else
        // Checking other b-quark mother candidates of the particular previous hadron
        if((int)LastQuarkIds.at(iHad).size() > theLastQuarkIndex+1) {
            if(fixExtraSameFlavours(iHad, hadIndices, hadMothers, hadMothersIndices, isFromTopWeakDecay, LastQuarkIds, LastQuarkMotherIds, lastQuarkIndices, hadronFlavour, checkedHadronIds, theLastQuarkIndex+1) ) {
                ambiguityResolved = true;
                break;
            };
        } 

    }       // End of loop over all previous hadrons
    checkedHadronIds.erase(hadId);      // Removing the hadron from the list of checked hadrons
    if(nSameFlavourHadrons>0 && !ambiguityResolved) {
        hadronFlavour.at(hadId) = qmFlav*21;
//         printf("    3. Setting flavour of hadron %d to %d\n", hadId, hadronFlavour.at(hadId));

        return true;
    }
    lastQuarkIndices.at(hadId) = lastQuarkIndex;
    hadronFlavour.at(hadId) = hadFlavour;
//     printf("   Setting flavour of hadron %d to %d\n", hadId, hadronFlavour.at(hadId));
    return true;
}


//define this as a plug-in
DEFINE_FWK_MODULE ( GenHFHadronAnalyzer );
