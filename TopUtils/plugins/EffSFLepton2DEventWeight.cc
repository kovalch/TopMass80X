#include "FWCore/Utilities/interface/EDMException.h"
#include "TopAnalysis/TopUtils/plugins/EffSFLepton2DEventWeight.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
 
// void
// EffSFLepton2DEventWeight::beginJob()
// {
//   if( histoname_==" " ) throw edm::Exception( edm::errors::Configuration, "Not correct histname!" );

// }

EffSFLepton2DEventWeight::EffSFLepton2DEventWeight(const edm::ParameterSet& cfg):
  particles_            ( cfg.getParameter<edm::InputTag>    ("particles") ),
  sysVar_               ( cfg.getParameter<std::string>      ("sysVar"   ) ),
  verbose_              ( cfg.getParameter<int>              ("verbose"  ) ),
  filename_             ( cfg.getParameter<edm::FileInPath>  ("filename" ) ),
  histoname_            ( cfg.getParameter<std::string>      ("histoname" ) ),
  additionalSystErr_    ( cfg.getParameter<double>      ("additionalSystErr") ),
  etaMaxVal_ 		( cfg.getParameter<double>      ("etaMaxVal") )
{
  produces<double>();
  
  mayConsume<edm::View<reco::Candidate>>(particles_);

  // debugging control outpu
  if(verbose_>=2){
    std::cout << "executing EffSFLepton2DEventWeight with the following options"<< std::endl;   
    std::cout << "particles_            : " << particles_             << std::endl;
    std::cout << "sysVar_               : " << sysVar_                << std::endl;
    std::cout << "verbose_              : " << verbose_               << std::endl;
    std::cout << "filename_             : " << filename_              << std::endl;
    std::cout << "histoname_            : " << histoname_              << std::endl;
    std::cout << "additionalSystErr_     : " << additionalSystErr_      << std::endl;
    
  }
  
  // load TFile Service
  edm::Service<TFileService> fs;
  if( !fs ){
    throw edm::Exception( edm::errors::Configuration,
			  "TFile Service is not registered in cfg file" );
  }
  /// booking of histogram for muon ID+trigger eff SF
  hists_["lepEffSF"] = fs->make<TH1F>( "lepEffSF", "lepEffSF", 40, 0., 2.);
  
  /// get SF 2D histo from specified input file and input histoname
  if(filename_.location()){
    file_ = new TFile((TString)filename_.fullPath());
    if(!(file_->IsZombie())){
      if(verbose_>=1) std::cout<<filename_.fullPath()<<" opened"<<std::endl;
      if(!( histoname_==" ")){
	hists2D = (TH2F*) file_->Get(histoname_.c_str());      
      	if(verbose_>=1) std::cout<< histoname_ <<" histo used"<<std::endl;
	
        /// get pt and eta dependend 2D histo
	/// considered range for 2D SFhisto
	//check if the histo axis are not swapped; x-axis: pt, y-axis:  eta 

	if (hists2D->GetXaxis()->GetXmax() > etaMaxVal_){
	  ptmin =hists2D->GetXaxis()->GetXmin();
	  ptmax =hists2D->GetXaxis()->GetXmax();
	  etamin=hists2D->GetYaxis()->GetXmin();
	  etamax=hists2D->GetYaxis()->GetXmax();
	} else {
	  ptmin =hists2D->GetYaxis()->GetXmin();
	  ptmax =hists2D->GetYaxis()->GetXmax();
	  etamin=hists2D->GetXaxis()->GetXmin();
	  etamax=hists2D->GetXaxis()->GetXmax(); 
	}

	if(verbose_>=1){
	  std::cout << "x range: " << ptmin  << ".." << ptmax  << std::endl;
	  std::cout << "y range: " << etamin << ".." << etamax << std::endl;
	}
      }
    }
    else{
      std::cout<<filename_.fullPath()<<" not found!!!!! Efficiencies cannot be taken from this file!!! Default taken"<<std::endl;
      filename_ = edm::FileInPath();
    }
  }
}

EffSFLepton2DEventWeight::~EffSFLepton2DEventWeight()
{
  if(filename_.location()) {if(!(file_->IsZombie())) file_->Close();}
}

void
EffSFLepton2DEventWeight::produce(edm::Event& evt, const edm::EventSetup& setup)
{
  // get object collection
  // (designed for lepton or electron collection)
  edm::Handle<edm::View<reco::Candidate> > particles;
  evt.getByLabel(particles_, particles);
  
  // variables used
  double pt, eta;
  double result     =  1.0;
  double errorUp    = -1.0;
  double errorDown  = -1.0;
  double errorCorr    = 0.;
  int numPart       =  0;
  if(numPart>1){
    std::cout << "WARNING in EffSFLepton2DEventWeight: more than one";
    std::cout << " object in input collection " << particles_;
    std::cout << " found! the module was designed for lepton eff with";
    std::cout << " only 1 input object. SF for the last entry in the";
    std::cout << " collection will be used!" << std::endl;
      }
  // loop object collection
  for(edm::View<reco::Candidate>::const_iterator part=particles->begin(); part!=particles->end(); ++part){

    // get kinematics
    pt  = part->pt();
    eta = part->eta();
    // count number of particles
    numPart++;
    // debug output
    if(verbose_>=2){
      std::cout << "particle " << numPart << " in particles_ collection: ";
      std::cout << "pt=" << pt << ", eta=" << eta << std::endl;
    }
    /// pt and eta dependent eff. SF
    if(filename_.location() && sysVar_!="FlatEffSF") {
      /// search for corresponding pt/eta bin of this particle
      int binPt =-1;
      int binEta=-1;
      if(     pt<ptmin  ) pt=ptmin;
      else if(pt>=ptmax  ) pt=ptmax-0.001;
      if(  eta<-1*etamax) eta=-1*etamax;
      else if(eta>etamax) eta=etamax;


      if (hists2D->GetXaxis()->GetXmax() > etaMaxVal_){
	binPt = hists2D->GetXaxis()->FindBin(pt );
	if (etamin == 0){
	  binEta= hists2D->GetYaxis()->FindBin(std::abs(eta));
	} else if (etamin < 0){
	  binEta= hists2D->GetYaxis()->FindBin(eta);
	}
      } else {
	binPt = hists2D->GetYaxis()->FindBin(pt );
	if (etamin == 0){
	  binEta= hists2D->GetXaxis()->FindBin(std::abs(eta));
	} else if (etamin < 0){
	  binEta= hists2D->GetXaxis()->FindBin(eta);
	}
      }


       
      // get SF and errors for this 
      if(binPt==-1||binEta==-1){
	std::cout << "ERROR in EffSFLepton2DEventWeight: can not identify bin in 2D effSF histo for ";
	std::cout << "pt=" << pt << "& eta=" << eta << std::endl;
	std::cout << "will use SF 1.+/-0. " << std::endl;
	result    =1.;
	errorUp   =0.;
	errorDown =0.;
      }
      else{
	if (hists2D->GetXaxis()->GetXmax() > etaMaxVal_){
	  result    = hists2D->GetBinContent(binPt, binEta);
	  errorUp   = hists2D->GetBinError  (binPt, binEta);
	} else {
	  result    = hists2D->GetBinContent(binEta, binPt);
	  errorUp   = hists2D->GetBinError  (binEta, binPt);  
	}
	errorDown = errorUp; // asymmetric errors not possible with 2D histos
      }
      if(verbose_>1) std::cout << "bin(pt,eta)=("<< binPt << "," << binEta << ")" << std::endl;
    }
    else{// flat SF
      result = 1.0;
    }
    // debug output
    if(verbose_>=1){
      std::cout << "loaded values from 2D histo" << std::endl;
      std::cout << "result   : " << result    << std::endl;
      std::cout << "errorUp  : " << errorUp   << std::endl;
      std::cout << "errorDown: " << errorDown << std::endl;
    }
    /// systematic variations for trigger eff. SF (normalisation and shape uncertainties)
    
    errorCorr = additionalSystErr_; 
//     std::cout<<numPart<< ": pt=" <<pt<< "; eta=" <<eta<< "; SF= "<<result<<"; errorUp="<< errorUp <<"; errorDown="<< errorDown<< "; errorCorr= "<< errorCorr << std::endl;
//     if(result==0 || result==INFINITY) LogDebug ( "EffSFLepton2DEventWeight" ) << " SF= "<<result<<"; errorUp="<< errorUp <<"; errorDown="<< errorDown;
    
    if     (sysVar_ == "combinedEffSFNormUpStat")   result += errorUp + errorCorr;
    else if(sysVar_ == "combinedEffSFNormDownStat") result -= errorDown + errorCorr;
    
    if(verbose_>=1) std::cout<<numPart<< ": pt=" <<pt<< "; eta=" <<eta<< "; SF= "<<result<<"; errorUp="<< errorUp <<"; errorDown="<< errorDown<<std::endl;
    hists_.find("lepEffSF" )->second->Fill(result);
  // break in order to have only one event weight (the one of the leading part.) in case of more part. in the event
  break;
  }
  std::auto_ptr<double> SFEventWeight(new double);
  *SFEventWeight = result;    
  evt.put(SFEventWeight);
}

// executed at the end after looping over all events
void 
    EffSFLepton2DEventWeight::endJob() 
{
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( EffSFLepton2DEventWeight );

