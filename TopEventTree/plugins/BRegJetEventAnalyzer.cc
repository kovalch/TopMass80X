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

#include "FWCore/MessageLogger/interface/MessageLogger.h"



#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"



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
JECUncSrcFile_(cfg.getParameter<edm::FileInPath>("JECUncSrcFile") ),
checkedIsPFJet(false), checkedJERSF(false), checkedJESSF(false), checkedTotalSF(false), checkedQGTag(false), checkedJESTotUnc(false), checkedBRegResult(false),
isPFJet(false),     hasJERSF(false),     hasJESSF(false),     hasTotalSF(false),     hasQGTag(false),     hasJESTotUnc(false), hasBRegResult(false),
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
//		variablePointer_.push_back(  &readerlChTrackPt_                  );//ok
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
//		variablePointer_.push_back(  &readerJesUncert_                   );//ok
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

std::vector<double> BRegJetEventAnalyzer::fillBRegJetAndReturnGBRResults(const edm::Event& evt, const edm::EventSetup& setup){
	fillBRegJet(evt, setup);
	return BRegJet->BRegGBRTrainResult;
}


std::vector<double> BRegJetEventAnalyzer::fillBRegJetAndReturnTMVAResults(const edm::Event& evt, const edm::EventSetup& setup){
	fillBRegJet(evt, setup);
	return BRegJet->BRegResult;
}

void BRegJetEventAnalyzer::fillBRegJet(const edm::Event& evt, const edm::EventSetup& setup){

	//////////////////////////////////////////////////////////////////////////
	//  PUJetId-variables
	////////////////////////////////////////////////////////////////////////

	// input variables
	edm::Handle<edm::ValueMap<StoredPileupJetIdentifier> > vmap;
	evt.getByLabel("puJetIdChs", vmap);





	//////////////////////////////////////////////////////////////////////////
	//  QGTaggerInfo
	////////////////////////////////////////////////////////////////////////


	edm::Handle<edm::ValueMap<float> >  QGTagsHandleAxis1MLP ;
	edm::Handle<edm::ValueMap<float> >  QGTagsHandleAxis2MLP ;
	edm::Handle<edm::ValueMap<float> >  QGTagsHandlePtDMLP ;
	edm::Handle<edm::ValueMap<float> >  QGTagsHandleMultMLP ;
	edm::Handle<edm::ValueMap<float> >  QGTagsHandleQGMLP ;

	evt.getByLabel("QGTagger","axis1MLP", QGTagsHandleAxis1MLP);
	evt.getByLabel("QGTagger","axis2MLP", QGTagsHandleAxis2MLP);
	evt.getByLabel("QGTagger","ptDMLP", QGTagsHandlePtDMLP);
	evt.getByLabel("QGTagger","multMLP", QGTagsHandleMultMLP);
//	std::cout << "trying to getbylabel" << std::endl;
	  evt.getByLabel("QGTagger","qgMLP", QGTagsHandleQGMLP);

	  edm::Handle<edm::ValueMap<float> >  QGTagsHandleMLP;
	  edm::Handle<edm::ValueMap<float> >  QGTagsHandleLikelihood;
	  evt.getByLabel("QGTagger","qgMLP", QGTagsHandleMLP);
	  evt.getByLabel("QGTagger","qgLikelihood", QGTagsHandleLikelihood);

//	++floatedmValueMap "QGTagger" "axis1MLP" "topMass" (productId = 5:100)
//	++floatedmValueMap "QGTagger" "axis2Likelihood" "topMass" (productId = 5:101)
//	++floatedmValueMap "QGTagger" "axis2MLP" "topMass" (productId = 5:102)
//	++floatedmValueMap "QGTagger" "multLikelihood" "topMass" (productId = 5:103)
//	++floatedmValueMap "QGTagger" "multMLP" "topMass" (productId = 5:104)
//	++floatedmValueMap "QGTagger" "ptDLikelihood" "topMass" (productId = 5:105)
//	++floatedmValueMap "QGTagger" "ptDMLP" "topMass" (productId = 5:106)
//	++floatedmValueMap "QGTagger" "qgLikelihood" "topMass" (productId = 5:107)
//	++floatedmValueMap "QGTagger" "qgMLP" "topMass" (productId = 5:108)


//    1 occurrences of key floatedmValueMap + "QGTagger" + "axis1MLP" "topMass"
//    1 occurrences of key floatedmValueMap + "QGTagger" + "axis2Likelihood" "topMass"
//    1 occurrences of key floatedmValueMap + "QGTagger" + "axis2MLP" "topMass"
//    1 occurrences of key floatedmValueMap + "QGTagger" + "multLikelihood" "topMass"
//    1 occurrences of key floatedmValueMap + "QGTagger" + "multMLP" "topMass"
//    1 occurrences of key floatedmValueMap + "QGTagger" + "ptDLikelihood" "topMass"
//    1 occurrences of key floatedmValueMap + "QGTagger" + "ptDMLP" "topMass"
//    1 occurrences of key floatedmValueMap + "QGTagger" + "qgLikelihood" "topMass"
//    1 occurrences of key floatedmValueMap + "QGTagger" + "qgMLP" "topMass"



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


	edm::Handle<edm::View<pat::Jet> > jets;
//	edm::Handle<std::vector<pat::Jet> > jets;
	evt.getByLabel(jets_, jets);

	double readerRlbReco=-999;
	std::vector<double> highestPtBJetsPt;
	std::vector<double> highestPtOtherJetsPt;

	unsigned short jetIndex = 0;
	unsigned short nB = 0;
	unsigned short nOther = 0;
	for(edm::View<pat::Jet>::const_iterator ijet=jets->begin(); ijet!=jets->end(); ++ijet,++jetIndex){
//	for(std::vector< pat::Jet >::const_iterator ijet = jets->begin(); ijet != jets->end(); ++ijet, ++jetIndex) {
		// write only kJetMAX_ jets into the event
		if(jetIndex == kJetMAX_) break;

		bool isB=false;

//		if(std::abs(ijet->partonFlavour())==5)isB=true;
		if(ijet->bDiscriminator("combinedSecondaryVertexBJetTags")>0.679)isB=true;//at least b-tagged //WARNING:HARD_CODED... SHOULD RETRIEVE THIS FROM SOME CONFIG!!!!

		if(isB){
			nB++;
			highestPtBJetsPt.push_back(ijet->pt());
		}
		else{
			nOther++;
			highestPtOtherJetsPt.push_back(ijet->pt());
		}
		if(nB>=2&&nOther>=2) break;
	}

	if(nB>=2&&nOther>=2){
		readerRlbReco = (highestPtBJetsPt.at(0)+highestPtBJetsPt.at(1))/(highestPtOtherJetsPt.at(0)+highestPtOtherJetsPt.at(1));
	}



	jetIndex = 0;
	nB = 0;
	for(edm::View<pat::Jet>::const_iterator ijet=jets->begin(); ijet!=jets->end(); ++ijet,++jetIndex){

//	for(std::vector< pat::Jet >::const_iterator ijet = jets->begin(); ijet != jets->end(); ++ijet, ++jetIndex) {




		// write only kJetMAX_ jets into the event
		if(jetIndex == kJetMAX_) break;

		// check only once per module run if the needed collections are available
		if(!checkedIsPFJet) { checkedIsPFJet = true; isPFJet = ijet->isPFJet(); }
		if(!checkedJESTotUnc  ) { checkedJESTotUnc   = true; if(ijet->hasUserFloat("jesTotUnc"      )) hasJESTotUnc   = true; }
		if(!checkedBRegResult ) { checkedBRegResult  = true; if(ijet->hasUserFloat("BRegResult"     )) hasBRegResult  = true; }


//		   UsersFavoriteData const & usersData = myMuon.userData<UsersFavoriteData>("myData");


//	     edm::RefToBase<pat::Jet> jetRef = jets->refAt(jetIndex);
////	    		 (edm::Ref<PFJetCollection>(jets,ijet));
//
//	     if (QGTagsHandleMLP.isValid()) std::cout << "MLP: " << (*QGTagsHandleMLP)[jetRef] << std::endl;
//	     if (QGTagsHandleLikelihood.isValid()) std::cout << "Likelihood: " << (*QGTagsHandleLikelihood)[jetRef] << std::endl;
//
//	     if (QGTagsHandleMLP.isValid()) std::cout << "francescos handle valid" << std::endl;
//	     if (QGTagsHandleQGMLP.isValid()) std::cout << "Hennings handle valid" << std::endl;


//        edm::RefToBase<pat::Jet> refToJetWithValueMaps = jetRef;
		edm::RefToBase<pat::Jet> refToJetWithValueMaps;

		if(hasBRegResult){
//			std::cout << "using stored reference" << std::endl;
			refToJetWithValueMaps =  * ijet->userData<edm::RefToBase<pat::Jet> >("refToJetWithValueMaps") ;
		}
		else {
//			std::cout << "using reference from stored collection" <<t std::endl;
			refToJetWithValueMaps =  jets->refAt(jetIndex);
		}
//		if (QGTagsHandleAxis1MLP.isValid()) std::cout << "Axis1: " << (*QGTagsHandleAxis1MLP)[refToJetWithValueMaps] << std::endl;
//		if (QGTagsHandleAxis2MLP.isValid()) std::cout << "Axis2: " << (*QGTagsHandleAxis2MLP)[refToJetWithValueMaps] << std::endl;
//		if (QGTagsHandlePtDMLP.isValid()) std::cout << "PtD: " << (*QGTagsHandlePtDMLP)[refToJetWithValueMaps] << std::endl;
//		if (QGTagsHandleMultMLP.isValid()) std::cout << "Mult: " << (*QGTagsHandleMultMLP)[refToJetWithValueMaps] << std::endl;
//		if (QGTagsHandleQGMLP.isValid()) std::cout << "QGMLP: " << (*QGTagsHandleQGMLP)[refToJetWithValueMaps] << std::endl;

//		std::cout << "einfach so QGMLP: " << (*QGTagsHandleQGMLP)[refToJetWithValueMaps] << std::endl;
//		std::cout << "einfach so Axis1: " << (*QGTagsHandleAxis1MLP)[refToJetWithValueMaps] << std::endl;
//		BRegJet->jetPtRaw
		  BRegJet->QGaxis1.push_back((*QGTagsHandleAxis1MLP)[refToJetWithValueMaps]);
		  BRegJet->QGaxis2.push_back((*QGTagsHandleAxis2MLP)[refToJetWithValueMaps]);
		  BRegJet->QGMult .push_back((*QGTagsHandleMultMLP)[refToJetWithValueMaps]);
		  BRegJet->QGPtD  .push_back((*QGTagsHandlePtDMLP)[refToJetWithValueMaps]);
		  BRegJet->QGMLP  .push_back((*QGTagsHandleQGMLP)[refToJetWithValueMaps]);
//		std::cout << "doing bregjeteventanalyser..." << std::endl;


		PileupJetIdentifier puIdentifier;
		// Or read it from the value map
		puIdentifier = (*vmap)[refToJetWithValueMaps];
//		std::cout << "puIdentifier.ptD() " << puIdentifier.ptD() << std::endl;

		BRegJet->PUIddZ       .push_back(puIdentifier.dZ());
		BRegJet->PUIddRMean   .push_back(puIdentifier.dRMean());
		BRegJet->PUIddr2Mean  .push_back(puIdentifier.dR2Mean());
		BRegJet->PUIdfrac01   .push_back(puIdentifier.frac01());
		BRegJet->PUIdfrac02   .push_back(puIdentifier.frac02());
		BRegJet->PUIdfrac03   .push_back(puIdentifier.frac03());
		BRegJet->PUIdfrac04   .push_back(puIdentifier.frac04());
		BRegJet->PUIdfrac05   .push_back(puIdentifier.frac05());
		BRegJet->PUIdbeta     .push_back(puIdentifier.beta());
		BRegJet->PUIdbetaStar .push_back(puIdentifier.betaStar());
		BRegJet->PUIdptD      .push_back(puIdentifier.ptD());

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


	    // GetTotal JES uncertainty
	    // get the uncertainty parameters from file, see
	    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECUncertaintySources
	    JetCorrectorParameters* TotalParam = new JetCorrectorParameters(JECUncSrcFile_.fullPath(), "Total");
	    // instantiate the jec uncertainty object
	    JetCorrectionUncertainty* TotalDeltaJEC = new JetCorrectionUncertainty(*TotalParam);
	    TotalDeltaJEC->setJetEta(ijet->eta()); TotalDeltaJEC->setJetPt(ijet->pt());
	    //set total JESUncertainty
		readerJesUncert_ = TotalDeltaJEC->getUncertainty(true);
		BRegJet->jesTotUnc   .push_back(readerJesUncert_);
		delete TotalParam; 
		delete TotalDeltaJEC;

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

		bool isOneOfLeading2B = false;
		if(std::abs(ijet->partonFlavour())==5)isOneOfLeading2B=true;
		if(isOneOfLeading2B){
			if(nB>2)isOneOfLeading2B=false;
			nB++;
		}

		BRegJet->OneOfLeading2B.push_back(isOneOfLeading2B);
		BRegJet->RlbReco.push_back(readerRlbReco);


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

//		if(hasBRegResult)std::cout << " BRegJet->BRegGBRTrainResult.back(): " << BRegJet->BRegGBRTrainResult.back() << " ijet->userFloat(BRegResult); " << ijet->userFloat("BRegResult"      ) << std::endl;

		if(hasBRegResult&&(float)BRegJet->BRegGBRTrainResult.back()!=ijet->userFloat("BRegResult")){
//TMVA		if(hasBRegResult&&(float)BRegJet->BRegResult.back()!=ijet->userFloat("BRegResult")){
		    edm::LogError msg("BRegression");
		    msg << "B regression result stored in jet and recalculated do not match. Please check your configuration accordingly \n";
		    throw cms::Exception("Configuration Error");
//
//			assert(BRegJet->BRegGBRTrainResult.back()==ijet->userFloat(BRegResult));
		}
//		std::cout << " ->BRegGBRTrainResult: " << BRegJet->BRegGBRTrainResult.back() << " >userFloat(BRegResult); " << ijet->userFloat("BRegResult") << ">userFloat(jesTotUnc)" << ijet->userFloat("jesTotUnc"      ) << " jetPtCorr " << BRegJet->jetPtCorr.back()  << " jetEta " << BRegJet->jetEta.back() << std::endl;
	}


}


void
BRegJetEventAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{

//	 std::cout << "==============================running analyze, before fill====================================" << std::endl;
	 fillBRegJet(evt,  setup);
//	 std::cout << "==============================running analyze, after fill ====================================" << std::endl;


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

// new BRegJetEvent is instantiated (should only be done once per job)
void
BRegJetEventAnalyzer::iniBRegEvent()
{
	BRegJet = new BRegJetEvent();
}

void
BRegJetEventAnalyzer::beginJob()
{
	if( !trs ) throw edm::Exception( edm::errors::Configuration, "TreeRegistryService is not registered in cfg file!" );

	iniBRegEvent();
	//BRegJet = new BRegJetEvent();
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
