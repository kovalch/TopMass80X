/*
 * BRegJetEventAnalyzer.cc
 *
 *  Created on: May 14, 2013
 *      Author: kirschen
 */

//#include <memory>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TopMass/TopEventTree/interface/TreeRegistryService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "TopMass/TopEventTree/plugins/BRegJetEventAnalyzer.h"

BRegJetEventAnalyzer::BRegJetEventAnalyzer(const edm::ParameterSet& cfg):
jets_        (cfg.getParameter<edm::InputTag>("jets")),
//allJets_     (cfg.getParameter<edm::InputTag>("allJets")),
//noPtEtaJets_ (cfg.getParameter<edm::InputTag>("noPtEtaJets")),
rho_tag_(cfg.getParameter<edm::InputTag>("rho_tag")),
rho25_tag_(cfg.getParameter<edm::InputTag>("rho25_tag")),
gluonTagName_  (cfg.getParameter<edm::InputTag>("gluonTagSrc").encode()),
kJetMAX_(cfg.getParameter<int>("maxNJets")),
BRegJet(0),
mva_name_(cfg.getParameter<std::string>("mva_name")),
mva_path_(cfg.getParameter<std::string>("mva_path")),
writeOutVariables_(cfg.getParameter<bool>("writeOutVariables")),
checkedIsPFJet(false), checkedJERSF(false), checkedJESSF(false), checkedTotalSF(false), checkedQGTag(false), checkedJESTotUnc(false),
       isPFJet(false),     hasJERSF(false),     hasJESSF(false),     hasTotalSF(false),     hasQGTag(false),     hasJESTotUnc(false) 
{

  //  outputJets_ = jets_.label();



//  //Initialize reader (adjusted with weight-files during first create)
//  if(mvaPath!=""){
//  reader_ = new TMVA::Reader( "!Color:Silent" );    
//  reader_->AddVariable( "Tlj_JetPt", &jetPt_ );
//  reader_->AddVariable( "Tlj_JetEta", &jetEta_ );
//  reader_->AddVariable( "Tlj_JetArea", &jetArea_);
//  reader_->AddVariable( "Tlj_Rho", &rho_ );
//  reader_->AddSpectator( "JetMETL1L2L3:=Tlj_L1L2L3JetPt/Tlj_JetPt",  &jetPt_); //not too reasonable
//
//  reader_->BookMVA(name,mvaPath);
//  //  }
//  }
//  else reader_=0;


}

void
BRegJetEventAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{





  //////////////////////////////////////////////////////////////////////////
  //  Get the computer for the CSV
  ////////////////////////////////////////////////////////////////////////

//  ESHandle<JetTagComputer> handle;
//  iSetup.get<JetTagComputerRecord>().get("combinedSecondaryVertex", handle);
//  const GenericMVAJetTagComputer *computer = dynamic_cast<const GenericMVAJetTagComputer*>(handle.product());



  //instantiate a tagging variable computer for unification of some calculations like vertex mass corrections
  edm::ESHandle<JetTagComputer> computerHandle;;
  setup.get<JetTagComputerRecord>().get( "combinedSecondaryVertex", computerHandle );
  const GenericMVAJetTagComputer *computer = dynamic_cast<const GenericMVAJetTagComputer*>( computerHandle.product() );
  if (!computer){
    edm::LogError("DataLost")<<"computer missing !!!"<<std::endl;
    exit(1);
  }
  computer->passEventSetup(setup);



  //////////////////////////////////////////////////////////////////////////////
  // INIT JetEvent
  ////////////////////////////////////////////////////////////////////////////

  BRegJet->init();

  //////////////////////////////////////////////////////////////////////////
  // JETS
  ////////////////////////////////////////////////////////////////////////

  edm::Handle<std::vector<pat::Jet> > jets;
  evt.getByLabel(jets_, jets);

  std::vector< std::string > JECSets   = jets->begin()->availableJECSets();
  std::vector< std::string > JECLevels = jets->begin()->availableJECLevels ();
  //  for(unsigned int i =0; i< JECSets.size();i++)std::cout << "JECSets: " << JECSets.at(i) << std::endl;
  //  for(unsigned int i =0; i< JECLevels.size();i++)std::cout << "JECLevels: " << JECLevels.at(i) << std::endl;


  unsigned short jetIndex = 0;
  for(std::vector< pat::Jet >::const_iterator ijet = jets->begin(); ijet != jets->end(); ++ijet, ++jetIndex) {
    // write only kJetMAX_ jets into the event
    if(jetIndex == kJetMAX_) break;

    // check only once per module run if the needed collections are available
    if(!checkedIsPFJet) { checkedIsPFJet = true; isPFJet = ijet->isPFJet(); }
    if(!checkedJESTotUnc  ) { checkedJESTotUnc   = true; if(ijet->hasUserFloat("jesTotUnc"      )) hasJESTotUnc   = true; }

    if(isPFJet){
      BRegJet->fChargedHadron .push_back(ijet->chargedHadronEnergyFraction());
      BRegJet->nChargedPFConstituents.push_back(ijet->chargedHadronMultiplicity () + ijet->electronMultiplicity () +ijet->muonMultiplicity () );

      const std::vector < reco::PFCandidatePtr > PFConstituents = ijet->getPFConstituents ();
      //    std::cout << "starting to list PF constituents:" <<std::endl; 
      BRegJet->leadingChargedConstPt.push_back(-1);
      for (std::vector<reco::PFCandidatePtr>::const_iterator constituent = PFConstituents.begin(); constituent!=PFConstituents.end(); ++constituent){
	//      std::cout << (*constituent)->pt() << " " << (*constituent)->charge() <<std::endl; 
	if(TMath::Abs((*constituent)->charge())>0){
	  BRegJet->leadingChargedConstPt.back() =  (*constituent)->pt();
	  break;
	}
      }
    }//end isPFJet
  

    //Jet area
    BRegJet->jetArea.push_back(ijet->jetArea());
    if( ijet->jetArea() < 0 ) {
      edm::LogError("BadArea") << "Area negative!";
    }
    if( ijet->jetArea() != ijet->jetArea() ) {
      edm::LogError("BadArea") << "Area is NaN!";
    }

    //secondary vertex information
    //    std::cout << "trying to retrieve secondaryVertexTagInfos, first check if it exists" << std::endl;
    //    std::cout << "hasTagInfo('')" << ijet->hasTagInfo("") << std::endl;
    //    std::cout << "hasTagInfo('secondaryVertex')" << ijet->hasTagInfo("secondaryVertex") << std::endl;
    //    std::cout << "hasTagInfo('secondaryVertexTagInfos')" << ijet->hasTagInfo("secondaryVertexTagInfos") << std::endl;
    //    std::cout << "hasTagInfo('secondaryVertexTagInfosAOD')" << ijet->hasTagInfo("secondaryVertexTagInfosAOD") << std::endl;
    //    std::cout << "hasTagInfo('recoBase')" << ijet->hasTagInfo("recoBase") << std::endl;
    const reco::SecondaryVertexTagInfo &svTagInfo = *ijet->tagInfoSecondaryVertex();
    //    std::cout << "test" << svTagInfo.nVertices() << std::endl;
    //    const reco::SecondaryVertexTagInfo &svTagInfo = *ijet->tagInfoSecondaryVertex();
    //    //    const reco::SecondaryVertexTagInfo &svTagInfo = *ijet->tagInfoSecondaryVertex("secondaryVertex");
    //    //    if(*ijet->tagInfoSecondaryVertex("secondaryVertex")==NULL)std::cout << "something went wrong: no secondaryVertexTagInfos found" << std::endl;
    BRegJet->nSV.push_back(svTagInfo.nVertices());
    if(svTagInfo.nVertices()>0){
      BRegJet->SVChi2.push_back(svTagInfo.secondaryVertex(0).chi2());         
      BRegJet->SV3DLength.push_back(svTagInfo.flightDistance(0).value());
      BRegJet->SV3DLengthError.push_back(svTagInfo.flightDistance(0).error());       
      
      //create more advanced variables using Computer
      std::vector<const reco::BaseTagInfo*>  baseTagInfos;
      baseTagInfos.push_back( ijet->tagInfoTrackIP ("impactParameter"  /*&(*ipTagInfo[thisJetRef])  */)); 
      baseTagInfos.push_back( ijet->tagInfoSecondaryVertex("secondaryVertex" /* &(*svTagInfo[thisJetRef]) */));
      JetTagComputer::TagInfoHelper helper(baseTagInfos);
      reco::TaggingVariableList vars = computer->taggingVariables(helper);
      
      if(vars.checkTag(reco::btau::vertexMass)) {
	//	std::cout << "angeblich hat es geklappt" << std::endl;
	BRegJet->SVMass.push_back( vars.get(reco::btau::vertexMass));

      }
      else  BRegJet->SVMass.push_back( -9999 );
      
      const reco::Vertex &vertex = svTagInfo.secondaryVertex(0);
      BRegJet->SVPt.push_back( vertex.p4().pt());
    }
    else{
      BRegJet->SVChi2.push_back          ( -1);
      BRegJet->SV3DLength.push_back      ( -1);
      BRegJet->SV3DLengthError.push_back ( -1);
      BRegJet->SVMass.push_back          ( -1);
      BRegJet->SVPt.push_back            ( -1); 
    }

    //check if daughters are present 
    if(ijet->numberOfDaughters() &&
       ijet->daughterPtr(0).isAvailable()) { 
      BRegJet->EtWeightedSigmaPhi.push_back( ijet->phiphiMoment() > 0 ? sqrt(ijet->phiphiMoment()) : 0);
      BRegJet->EtWeightedSigmaEta.push_back( ijet->etaetaMoment() > 0 ? sqrt(ijet->etaetaMoment()) : 0);
    } else {
      BRegJet->EtWeightedSigmaPhi.push_back( 0);
      BRegJet->EtWeightedSigmaEta.push_back( 0);
    }
    
    //JESUncertainty
    if(hasJESTotUnc  ) BRegJet->jesTotUnc   .push_back(ijet->userFloat("jesTotUnc"      ));
  

    //rho and rho25
    edm::Handle<double> pRho;
    evt.getByLabel(rho_tag_,pRho);
    edm::Handle<double> pRho25;
    evt.getByLabel(rho25_tag_,pRho25);

    if(*pRho < 0) edm::LogError("BadRho") << "Rho negative!";
    if(*pRho != *pRho) edm::LogError("BadRho") << "Rho is NaN!";
    if(*pRho25 < 0) edm::LogError("BadRho25") << "Rho25 negative!";
    if(*pRho25 != *pRho25) edm::LogError("BadRho25") << "Rho25 is NaN!";

    BRegJet->Rho   .push_back(*pRho);
    BRegJet->Rho25 .push_back(*pRho25);


//    //add breg-info to patjet //wont work as this is no EDProducer
//    ijet->addUserFloat("BRegResult", 0.1)
  }

  if(writeOutVariables_)trs->Fill();
}

void
BRegJetEventAnalyzer::beginJob()
{
  if( !trs ) throw edm::Exception( edm::errors::Configuration, "TreeRegistryService is not registered in cfg file!" );

  BRegJet = new BRegJetEvent();
  trs->Branch("BRegJet", BRegJet);
}

void
BRegJetEventAnalyzer::endJob()
{
}

BRegJetEventAnalyzer::~BRegJetEventAnalyzer()
{
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BRegJetEventAnalyzer);
