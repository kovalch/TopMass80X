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
GBRmva_path_(cfg.getParameter<std::string>("GBRmva_path")),
writeOutVariables_(cfg.getParameter<bool>("writeOutVariables")),
checkedIsPFJet(false), checkedJERSF(false), checkedJESSF(false), checkedTotalSF(false), checkedQGTag(false), checkedJESTotUnc(false),
isPFJet(false),     hasJERSF(false),     hasJESSF(false),     hasTotalSF(false),     hasQGTag(false),     hasJESTotUnc(false),
tempJetPtCorr_(0),tempJetMt_(0),tempJetEta_(0),
tempfChargedHadrons_(0),tempfElectrons_(0), tempfMuons_(0),
tempBTagCSV_(0), tempnChargedHadrons_(0),tempnChargedPFConstituents_(0),tempnPFConstituents_(0),
readerJetPtRaw_(0), readerJetArea_(0), readerJetEtWEightedSigmaPhi_(0),
readerJesUncert_(0), readerSVtx3dLength_(0), readerSVtx3dLengthError_(0),
readerSVtxMass_(0), readerSVtxPt_(0), readerlChTrackPt_(0), readerRho25_(0),
readerSoftMuonPt_(0), readerSoftMuonRatioRel_(0), readerSoftMuonDeltaR_(0),
readerSoftElectronPt_(0), readerSoftElectronRatioRel_(0),
readerSoftElectronDeltaR_(0)

{

	TFile *fgbropt = new TFile(GBRmva_path_.c_str(),"READ");;
	gbropt_ = (GBRForest*)fgbropt->Get("gbrtrain");
	varlist_ = (std::vector<std::string>*)fgbropt->Get("varlist");


	//Initialize reader for TMVA (adjusted with weight-files during first create)
	if(mva_path_!=""){
		reader_ = new TMVA::Reader( "!Color:Silent" );
		//		reader_ = new TMVA::Reader( "!Color" );


		//collect pointers to variables holding
		variablePointer_.push_back(  &tempJetPtCorr_                     );//ok
		variablePointer_.push_back(  &tempJetEta_                        );//ok
		variablePointer_.push_back(  &readerRho25_                       );//ok, was wrong
		variablePointer_.push_back(  &readerJetArea_                     );//ok, was wrong
		variablePointer_.push_back(  &readerJetEtWEightedSigmaPhi_       );//ok
		//5 width and PU
		variablePointer_.push_back(  &readerlChTrackPt_                  );//ok
		variablePointer_.push_back(  &tempfChargedHadrons_               );//ok
		variablePointer_.push_back(  &tempfElectrons_                    );//ok
		variablePointer_.push_back(  &tempfMuons_                        );//ok
		//9 charged constituents
		variablePointer_.push_back(  &readerSVtx3dLength_                );//ok
		variablePointer_.push_back(  &readerSVtx3dLengthError_           );//ok
		variablePointer_.push_back(  &readerSVtxMass_                    );//ok
		variablePointer_.push_back(  &readerSVtxPt_                      );//ok
		variablePointer_.push_back(  &tempBTagCSV_                       );//ok
		//14 SVtx and btag
		variablePointer_.push_back(  &readerSoftMuonPt_             );  //added in recently
		variablePointer_.push_back(  &readerSoftMuonRatioRel_       );  //added in recently
		variablePointer_.push_back(  &readerSoftMuonDeltaR_         );  //added in recently
		variablePointer_.push_back(  &readerSoftElectronPt_         );  //added in recently
		variablePointer_.push_back(  &readerSoftElectronRatioRel_   );  //added in recently
		variablePointer_.push_back(  &readerSoftElectronDeltaR_     );  //added in recently
		//20 soft lepton variables
		variablePointer_.push_back(  &readerJetPtRaw_                    );//ok
		variablePointer_.push_back(  &tempJetMt_                         );//ok
		variablePointer_.push_back(  &readerJesUncert_                   );//ok
		variablePointer_.push_back(  &tempnPFConstituents_               );//ok
		variablePointer_.push_back(  &tempnChargedHadrons_               );//ok
		variablePointer_.push_back(  &tempnChargedPFConstituents_        );//ok

		UInt_t nvars = varlist_->size();

		assert(varlist_->size()==variablePointer_.size());


		vals_ = new Float_t[nvars];


		for (UInt_t ivar=0; ivar<varlist_->size(); ++ivar) {
//reader does not accept double			reader_->AddVariable( varlist_->at(ivar),                variablePointer_.at(ivar));//ok
			reader_->AddVariable( varlist_->at(ivar),                &vals_[ivar]);//ok
		}
//reader does not accept double		reader_->AddSpectator("BRegJet.genJetPt",/*dummy*/variablePointer_.at(0));
//reader does not accept double		reader_->AddSpectator("BRegJet.genPartonPt",/*dummy*/variablePointer_.at(0));
		reader_->AddSpectator("BRegJet.genJetPt",/*dummy*/&vals_[0]);
		reader_->AddSpectator("BRegJet.genPartonPt",/*dummy*/&vals_[0]);

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

	unsigned short jetIndex = 0;
	for(std::vector< pat::Jet >::const_iterator ijet = jets->begin(); ijet != jets->end(); ++ijet, ++jetIndex) {
		// write only kJetMAX_ jets into the event
		if(jetIndex == kJetMAX_) break;

		// check only once per module run if the needed collections are available
		if(!checkedIsPFJet) { checkedIsPFJet = true; isPFJet = ijet->isPFJet(); }
		if(!checkedJESTotUnc  ) { checkedJESTotUnc   = true; if(ijet->hasUserFloat("jesTotUnc"      )) hasJESTotUnc   = true; }

		if(isPFJet){
//			BRegJet->fChargedHadron .push_back(ijet->chargedHadronEnergyFraction());
		  tempnChargedPFConstituents_ = ijet->chargedHadronMultiplicity () + ijet->electronMultiplicity () +ijet->muonMultiplicity ();

			const std::vector < reco::PFCandidatePtr > PFConstituents = ijet->getPFConstituents ();
			//    std::cout << "starting to list PF constituents:" <<std::endl;
			readerlChTrackPt_ = -1;
			for (std::vector<reco::PFCandidatePtr>::const_iterator constituent = PFConstituents.begin(); constituent!=PFConstituents.end(); ++constituent){
				//      std::cout << (*constituent)->pt() << " " << (*constituent)->charge() <<std::endl;
				if(TMath::Abs((*constituent)->charge())>0){
				  readerlChTrackPt_ = (*constituent)->pt();
				  break;
				}
			}
		}//end isPFJet
		else{
//			BRegJet->fChargedHadron          .push_back(-1);
    		        tempnChargedPFConstituents_ = -999;
			readerlChTrackPt_ =  -999;
			
		}
		BRegJet->nChargedPFConstituents  .push_back(tempnChargedPFConstituents_);
		BRegJet->leadingChargedConstPt   .push_back(readerlChTrackPt_);


		readerJetPtRaw_=ijet->correctedP4("Uncorrected").Pt();
		BRegJet->jetPtRaw.push_back(readerJetPtRaw_);


		//Jet area
		readerJetArea_ = ijet->jetArea();
		BRegJet->jetArea.push_back(readerJetArea_);
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
			readerSVtx3dLength_ = svTagInfo.flightDistance(0).value();
			readerSVtx3dLengthError_ = svTagInfo.flightDistance(0).error();

			//create more advanced variables using Computer
			std::vector<const reco::BaseTagInfo*>  baseTagInfos;
			baseTagInfos.push_back( ijet->tagInfoTrackIP ("impactParameter"  /*&(*ipTagInfo[thisJetRef])  */));
			baseTagInfos.push_back( ijet->tagInfoSecondaryVertex("secondaryVertex" /* &(*svTagInfo[thisJetRef]) */));
			JetTagComputer::TagInfoHelper helper(baseTagInfos);
			reco::TaggingVariableList vars = computer->taggingVariables(helper);

			if(vars.checkTag(reco::btau::vertexMass)) {
				//	std::cout << "angeblich hat es geklappt" << std::endl;
			  readerSVtxMass_= vars.get(reco::btau::vertexMass);

			}
			else  readerSVtxMass_= -999;

			const reco::Vertex &vertex = svTagInfo.secondaryVertex(0);
			readerSVtxPt_ =  vertex.p4().pt();
		}
		else{
			BRegJet->SVChi2.push_back          ( -999);
			readerSVtx3dLength_ = -999;
			readerSVtx3dLengthError_ = -999;
			readerSVtxMass_ = -999;
			readerSVtxPt_ = -999;
			
		}

		BRegJet->SV3DLength.push_back      (readerSVtx3dLength_);
		BRegJet->SV3DLengthError.push_back (readerSVtx3dLengthError_ );
		BRegJet->SVMass.push_back(readerSVtxMass_  );
		BRegJet->SVPt.push_back            ( readerSVtxPt_);

		//check if daughters are present
		if(ijet->numberOfDaughters() &&	ijet->daughterPtr(0).isAvailable()) {
		  readerJetEtWEightedSigmaPhi_= ijet->phiphiMoment() > 0 ? sqrt(ijet->phiphiMoment()) : 0;
		  BRegJet->EtWeightedSigmaEta.push_back( ijet->etaetaMoment() > 0 ? sqrt(ijet->etaetaMoment()) : 0);
		} else {
		  readerJetEtWEightedSigmaPhi_= -999;
		  BRegJet->EtWeightedSigmaEta.push_back( -999);
		}
		BRegJet->EtWeightedSigmaPhi.push_back( readerJetEtWEightedSigmaPhi_);

		//JESUncertainty
		if(hasJESTotUnc  ) readerJesUncert_=ijet->userFloat("jesTotUnc"      );
		else readerJesUncert_ = -999;
		BRegJet->jesTotUnc   .push_back(readerJesUncert_);

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

		readerRho25_ = *pRho25;
		BRegJet->Rho25 .push_back(readerRho25_);


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
			//reader variables set consistently later

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
			BRegJet->nSoftMuons.push_back( -999);
			BRegJet->SoftMuonPt        .push_back( -999);
			BRegJet->SoftMuonPtRel     .push_back( -999);
			BRegJet->SoftMuonRatioRel  .push_back( -999);
			BRegJet->SoftMuonDeltaR    .push_back( -999);
		}

		readerSoftMuonPt_ = BRegJet->SoftMuonPt.back();
		readerSoftMuonRatioRel_ = BRegJet->SoftMuonRatioRel.back();
		readerSoftMuonDeltaR_ = BRegJet->SoftMuonDeltaR.back();



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
		    // reader variables set consistently after if/else

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
			BRegJet->nSoftElectrons.push_back( -999);
			BRegJet->SoftElectronPt        .push_back( -999);
			BRegJet->SoftElectronPtRel     .push_back( -999);
			BRegJet->SoftElectronRatioRel  .push_back( -999);
			BRegJet->SoftElectronDeltaR    .push_back( -999);
		}

	    readerSoftElectronPt_ = BRegJet->SoftElectronPt.back();
	    readerSoftElectronRatioRel_ = BRegJet->SoftElectronRatioRel.back();
	    readerSoftElectronDeltaR_ = BRegJet->SoftElectronDeltaR.back();



		if(ijet->genJet())BRegJet->genJetPt.push_back(ijet->genJet()->pt());
		else BRegJet->genJetPt.push_back(-999);

		if(ijet->genParton())BRegJet->genPartonPt.push_back(ijet->genParton()->pt());
		else BRegJet->genPartonPt.push_back(-999);


		//temporary variables needed to calculate b-regression result
		tempJetPtCorr_ = ijet->pt();
		tempJetMt_     = ijet->mt();
		tempJetEta_    = ijet->eta();
		BRegJet->jetPtCorr.push_back(tempJetPtCorr_);
		BRegJet->jetMt.push_back(tempJetMt_);
		BRegJet->jetEta.push_back(tempJetEta_);

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


		for (UInt_t ivar=0; ivar<variablePointer_.size(); ++ivar) {
			vals_[ivar]=(float) *variablePointer_.at(ivar);
//			std::cout << varlist_->at(ivar) <<  ": " << (float) *variablePointer_.at(ivar) << std::endl;

		}
		Double_t mvaValue = (reader_->EvaluateRegression(mva_name_))[0];
		BRegJet->BRegResult.push_back(mvaValue);

		BRegJet->BRegGBRTrainResult.push_back(gbropt_->GetResponse(vals_));
//std::cout << "mvaValue: " << mvaValue << " BRegJet->BRegGBRTrainResult.back(): " << BRegJet->BRegGBRTrainResult.back() << std::endl;

	}

	if(writeOutVariables_)trs->Fill();
	else{
		  std::vector<double> tempBRegResult = BRegJet->BRegResult;
		  std::vector<double> tempBRegGBRTrainResult = BRegJet->BRegGBRTrainResult;
		  BRegJet->init();
		  BRegJet->BRegResult = tempBRegResult;
		  BRegJet->BRegGBRTrainResult = tempBRegGBRTrainResult;
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
	delete	gbropt_;
	delete vals_;
	delete varlist_;

}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BRegJetEventAnalyzer);
