/*
 * BRegJetEventAnalyzer.cc
 *
 *  Created on: May 14, 2013
 *      Author: kirschen
 */

//#include <memory>
//#include "TMVA/Factory.h"
//#include "TMVA/Tools.h"
//#include "TMVA/Config.h"
//#include "TMVA/Reader.h"
//#include "TMVA/MethodCuts.h"
#include "TMVA/Reader.h"
// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TopMass/TopEventTree/interface/TreeRegistryService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

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
isPFJet(false),     hasJERSF(false),     hasJESSF(false),     hasTotalSF(false),     hasQGTag(false),     hasJESTotUnc(false),
tempJetPtCorr_(0),tempJetMt_(0),tempJetEta_(0),
tempfChargedHadrons_(0),tempfElectrons_(0), tempfMuons_(0),
tempBTagCSV_(0), tempnChargedHadrons_(0),tempnPFConstituents_(0),
readerJetPtRaw_(0), readerJetArea_(0), readerJetEtWEightedSigmaPhi_(0),
readerJesUncert_(0), readerSVtx3dLength_(0), readerSVtx3dLengthError_(0),
readerSVtxMass_(0), readerSVtxPt_(0), readerlChTrackPt_(0), readerRho25_(0)
{

	//Initialize reader (adjusted with weight-files during first create)
	if(mva_path_!=""){
		reader_ = new TMVA::Reader( "!Color:Silent" );

		reader_->AddVariable( "Tlj_jetpt",                       &readerJetPtRaw_                    );
		reader_->AddVariable( "Tlj_l1l2l3jetpt",                 &tempJetPtCorr_                     );
		reader_->AddVariable( "Tlj_jetmt",                       &tempJetMt_                         );
		reader_->AddVariable( "Tlj_jeteta",                      &tempJetEta_                        );
		reader_->AddVariable( "Tlj_jetarea",                     &readerJetArea_                     );
		reader_->AddVariable( "Tlj_jetetweightedsigmaphi",       &readerJetEtWEightedSigmaPhi_       );
		reader_->AddVariable( "Tlj_jesUncert",                   &readerJesUncert_                   );
		reader_->AddVariable( "Tlj_SVtx3dLength",                &readerSVtx3dLength_                );
		reader_->AddVariable( "Tlj_SVtx3dLengthError",           &readerSVtx3dLengthError_           );
		reader_->AddVariable( "Tlj_SVtxMass",                    &readerSVtxMass_                    );
		reader_->AddVariable( "Tlj_SVtxPt",                      &readerSVtxPt_                      );
		reader_->AddVariable( "Tlj_lChTrackPt",                  &readerlChTrackPt_                  );
		reader_->AddVariable( "Tlj_CHF",                         &tempfChargedHadrons_               );
		reader_->AddVariable( "Tlj_FElectrons",                  &tempfElectrons_                    );
		reader_->AddVariable( "Tlj_FMuons",                      &tempfMuons_                        );
		reader_->AddVariable( "Tlj_btag",                        &tempBTagCSV_                       );
		reader_->AddVariable( "Tlj_NTot",                        &tempnPFConstituents_               );
		reader_->AddVariable( "Tlj_rho25",                       &readerRho25_                       );
		reader_->AddVariable( "Tlj_NChargedHadrons",             &tempnChargedHadrons_               );
		reader_->AddVariable( "Tlj_NChargedPFConst",             &tempnPFConstituents_               );

		reader_->BookMVA(mva_name_,mva_path_);
	}
	else reader_=0;


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

//	std::vector< std::string > JECSets   = jets->begin()->availableJECSets();
//	std::vector< std::string > JECLevels = jets->begin()->availableJECLevels ();
//	for(unsigned int i =0; i< JECSets.size();i++)std::cout << "JECSets: " << JECSets.at(i) << std::endl;
//	for(unsigned int i =0; i< JECLevels.size();i++)std::cout << "JECLevels: " << JECLevels.at(i) << std::endl;
//	std::cout << jets->begin()->currentJECLevel() << std::endl;

	unsigned short jetIndex = 0;
	for(std::vector< pat::Jet >::const_iterator ijet = jets->begin(); ijet != jets->end(); ++ijet, ++jetIndex) {
		// write only kJetMAX_ jets into the event
		if(jetIndex == kJetMAX_) break;

		// check only once per module run if the needed collections are available
		if(!checkedIsPFJet) { checkedIsPFJet = true; isPFJet = ijet->isPFJet(); }
		if(!checkedJESTotUnc  ) { checkedJESTotUnc   = true; if(ijet->hasUserFloat("jesTotUnc"      )) hasJESTotUnc   = true; }

		if(isPFJet){
//			BRegJet->fChargedHadron .push_back(ijet->chargedHadronEnergyFraction());
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
		else{
//			BRegJet->fChargedHadron          .push_back(-1);
			BRegJet->nChargedPFConstituents  .push_back(-1);
			BRegJet->leadingChargedConstPt   .push_back(-1);
		}

		BRegJet->jetPtRaw.push_back(ijet->correctedP4("Uncorrected").Pt());


		//Jet area
		BRegJet->jetArea.push_back(ijet->jetArea());
		if( ijet->jetArea() < 0 ) {
			edm::LogError("BadArea") << "Area negative!";
		}
		if( ijet->jetArea() != ijet->jetArea() ) {
			edm::LogError("BadArea") << "Area is NaN!";
		}

		const reco::SecondaryVertexTagInfo &svTagInfo = *ijet->tagInfoSecondaryVertex();
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
		else BRegJet->jesTotUnc   .push_back(-1);

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

		const reco::SoftLeptonTagInfo &softMuonTagInfo = *ijet->tagInfoSoftLepton("softMuon");
		const reco::SoftLeptonTagInfo &softElectronTagInfo = *ijet->tagInfoSoftLepton("softElectron");

		if(softMuonTagInfo.leptons()>0){
//			std::cout << softMuonTagInfo.leptons() << std::endl;
			BRegJet->nSoftMuons.push_back(softMuonTagInfo.leptons());
			edm::Handle<reco::MuonCollection> allmuons;
			evt.getByLabel("muons",allmuons);
			// loop over all muons in tag info and match them to muon collection
			// create pt sorted map of muons
			std::map<Double_t, const reco::Muon*> muMap;
			std::map<const reco::Muon*, UInt_t> tagInfoMap;

			for(UInt_t iMuon = 0; iMuon < softMuonTagInfo.leptons(); iMuon++){
				//loop on muon collection
				//count muons
				UInt_t iMuCounter = 0;
				for(size_t i=0; i < allmuons->size(); ++i){
					const reco::Muon & mu = (*allmuons)[i];
					reco::TrackRef globTrack = mu.globalTrack();
					reco::TrackRef softLepTrackRef =  (softMuonTagInfo.lepton(iMuon)).castTo<reco::TrackRef>();
					if( globTrack == softLepTrackRef ){   // found a matched muon
						iMuCounter++;
						muMap[mu.globalTrack()->pt()] = &mu ;
						tagInfoMap[&mu] = iMuon;
					}

				}
				if(iMuCounter != 1){
					std::cout<<"ERROR: iMuCounter != 1 this should never happen" << std::endl;
					exit(1);
				}
			}

			if( muMap.size() != softMuonTagInfo.leptons()){
				std::cout<<"ERROR: softMuonTagInfo.leptons():  this should never happen"<< std::endl;
				exit(1);
			}

			//consider only the Muon with highest pt from soft lepton tag info
			std::map<Double_t, const reco::Muon*>::reverse_iterator it= muMap.rbegin();
			const reco::Muon * mu =  it->second;

			//filling additional variables for leading muon
			BRegJet->SoftMuonPt        .push_back(softMuonTagInfo.lepton(tagInfoMap[mu])->pt());
			BRegJet->SoftMuonPtRel     .push_back(softMuonTagInfo.properties(tagInfoMap[mu]).ptRel);
			BRegJet->SoftMuonRatioRel  .push_back(softMuonTagInfo.properties(tagInfoMap[mu]).ratioRel);
			BRegJet->SoftMuonDeltaR    .push_back(softMuonTagInfo.properties(tagInfoMap[mu]).deltaR);

			if(softMuonTagInfo.leptons()>1){
				it++;
				const reco::Muon * secondMu =  it->second;
//				std::cout << "Warning: More than one muon in muontaginfo!" << std::endl;
				assert(softMuonTagInfo.lepton(tagInfoMap[mu])->pt()>softMuonTagInfo.lepton(tagInfoMap[secondMu])->pt());
//				std::cout << "At least lepton0 pt > lepton1 pt" << std::endl;
			}
		}
		else{
			BRegJet->nSoftMuons.push_back( 0);
			BRegJet->SoftMuonPt        .push_back( 0);
			BRegJet->SoftMuonPtRel     .push_back( 0);
			BRegJet->SoftMuonRatioRel  .push_back( 0);
			BRegJet->SoftMuonDeltaR    .push_back( 0);
		}


		if(softElectronTagInfo.leptons()>0){
			BRegJet->nSoftElectrons.push_back(softElectronTagInfo.leptons());

			edm::Handle<reco::GsfElectronCollection> allelectrons;
			evt.getByLabel("gsfElectrons",allelectrons);
			// loop over all electrons in tag info and match them to gsfelectron collection
			// create pt sorted map of electrons
			std::map<Double_t, const reco::GsfElectron*> elMap;
			std::map<const reco::GsfElectron*, UInt_t> tagInfoMap;

			for(UInt_t iElectron = 0; iElectron < softElectronTagInfo.leptons(); iElectron++)
			{
				//loop on electron collection
				//count electrons
				UInt_t iElCounter = 0;
				for(size_t i=0; i < allelectrons->size(); ++i){
					const reco::GsfElectron & el = (*allelectrons)[i];
					reco::GsfTrackRef gsftrack = el.gsfTrack();
					if(gsftrack.isNull() && !gsftrack.isAvailable()){
						std::cout<<"track is null"<<std::endl;
						continue;
					}
					reco::GsfTrackRef softLepTrack =  (softElectronTagInfo.lepton(iElectron)).castTo<reco::GsfTrackRef>();
					if( gsftrack == softLepTrack){
						iElCounter++;
						elMap[gsftrack->pt()] = &el ;
						tagInfoMap[&el] = iElectron;
					}
				}

				if(iElCounter != 1){
					std::cout<<"ERROR: iElCounter= "<<iElCounter<<"this should never happen" << std::endl;
					exit(1);
				}
			}


			if( elMap.size() != softElectronTagInfo.leptons()){
				std::cout<<"ERROR: softElectronTagInfo.leptons():  this should never happen"<< std::endl;
				exit(1);
			}

			std::map<Double_t, const reco::GsfElectron*>::reverse_iterator it= elMap.rbegin();
		    const reco::GsfElectron * el =  it->second;

			//filling additional variables for leading electron
			BRegJet->SoftElectronPt        .push_back(softElectronTagInfo.lepton(tagInfoMap[el])->pt());
			BRegJet->SoftElectronPtRel     .push_back(softElectronTagInfo.properties(tagInfoMap[el]).ptRel);
			BRegJet->SoftElectronRatioRel  .push_back(softElectronTagInfo.properties(tagInfoMap[el]).ratioRel);
			BRegJet->SoftElectronDeltaR    .push_back(softElectronTagInfo.properties(tagInfoMap[el]).deltaR);


			if(softElectronTagInfo.leptons()>1){
				it++;
				const reco::GsfElectron * secondEl =  it->second;
//				std::cout << "Warning: More than one electron in electrontaginfo!" << std::endl;
				assert(softElectronTagInfo.lepton(tagInfoMap[el])->pt()>softElectronTagInfo.lepton(tagInfoMap[secondEl])->pt());
//				std::cout << "At least lepton0 pt > lepton1 pt" << std::endl;
			}

		}
		else{
			BRegJet->nSoftElectrons.push_back( 0);
			BRegJet->SoftElectronPt        .push_back( 0);
			BRegJet->SoftElectronPtRel     .push_back( 0);
			BRegJet->SoftElectronRatioRel  .push_back( 0);
			BRegJet->SoftElectronDeltaR    .push_back( 0);
		}

		//temporary variables needed to calculate b-regression result
		tempJetPtCorr_ = ijet->pt();
		tempJetMt_     = ijet->mt();
		tempJetEta_    = ijet->eta();
		if(isPFJet){
			tempfChargedHadrons_ = ijet->chargedHadronEnergyFraction();
			tempfElectrons_      = ijet->electronEnergyFraction();
			tempfMuons_          = ijet->muonEnergyFraction();
			tempnChargedHadrons_ = ijet->chargedHadronMultiplicity ();
			tempnPFConstituents_ = (ijet->chargedHadronMultiplicity () +ijet->neutralHadronMultiplicity ()
					+ijet->photonMultiplicity ()  +ijet->electronMultiplicity () +ijet->muonMultiplicity ()
					+ijet->HFHadronMultiplicity ()  +ijet->HFEMMultiplicity () );//NTot
		}
		tempBTagCSV_ = ijet->bDiscriminator("combinedSecondaryVertexBJetTags");

		//		   factory->AddVariable( "Tlj_NChargedHadrons", "NCh", "units", 'F' );


		Double_t mvaValue = (reader_->EvaluateRegression(mva_name_))[0];
		BRegJet->BRegResult.push_back(mvaValue);

	}

	if(writeOutVariables_)trs->Fill();
	else{
		  std::vector<double> tempBRegResult = BRegJet->BRegResult;
		  BRegJet->init();
		  BRegJet->BRegResult = tempBRegResult;
		  trs->Fill();
	}
}

void
BRegJetEventAnalyzer::beginJob()
{
	if( !trs ) throw edm::Exception( edm::errors::Configuration, "TreeRegistryService is not registered in cfg file!" );

	BRegJet = new BRegJetEvent();
	trs->Branch("BRegJet.", BRegJet);
}

void
BRegJetEventAnalyzer::endJob()
{

}

BRegJetEventAnalyzer::~BRegJetEventAnalyzer()
{
	delete reader_;
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BRegJetEventAnalyzer);
