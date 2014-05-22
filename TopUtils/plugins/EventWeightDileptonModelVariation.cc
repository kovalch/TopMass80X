#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TMath.h"
#include <TH2.h>
#include <TF1.h>
#include <Math/Boost.h>

class EventWeightDileptonModelVariation : public edm::EDProducer {

public:
    explicit EventWeightDileptonModelVariation ( const edm::ParameterSet& );

private:
    virtual void produce ( edm::Event&, const edm::EventSetup& );
    virtual void beginJob();
    double getTopPtWeight();
    double getTopPtDataWeight();
    double getTopPtSoftenWeight();
    double getTopPtHardenWeight();
    double getTopPtTTsysWeight();
    double getTopEtaWeight();
    double getTopMassWeight();
    double getTopMassWeightLandau();
    double getTopRapidityWeight();

    edm::InputTag ttGenEvent_;
    edm::Handle<TtGenEvent> genEvt;
    std::string weightVariable_;
    double slope_;
    double weight1x_;
    double minWeight_;
    double maxWeight_;
    double landauMPV_;
    //double landauConstant_;
    double landauSigma_;
    double landauMoveX_;
    TH1 *modelWeightSum_;
    TH1 *pt_, *pttt_;
    TH1 *mass_;
    TH1 *rap_, *raptt_;
    TH1 *weightVsMass_, *weightVsPt_;
    TH1 *pttopWeightedEvwgt_, *ytopWeightedEvwgt_, *ptttbarWeightedEvwgt_, *yttbarWeightedEvwgt_, *mttbarWeightedEvwgt_;
    TH1 *evwght_;
    TF1* function;

};

EventWeightDileptonModelVariation::EventWeightDileptonModelVariation(const edm::ParameterSet& cfg):
    ttGenEvent_(cfg.getParameter<edm::InputTag>("ttGenEvent")),
    weightVariable_(cfg.getParameter<std::string>("weightVariable")),
    slope_(cfg.getParameter<double>("slope")),
    weight1x_(cfg.getParameter<double>("weight1x")),
    minWeight_(cfg.getParameter<double>("minWeight")),
    maxWeight_(cfg.getParameter<double>("maxWeight")),
    landauMPV_(cfg.getParameter<double>("landauMPV")),
    //landauConstant_(cfg.getParameter<double>("landauConstant")),
    landauSigma_(cfg.getParameter<double>("landauSigma")),
    landauMoveX_(cfg.getParameter<double>("landauMoveX"))
{
    produces<double>("");
    function=0;
}

void EventWeightDileptonModelVariation::beginJob()
{
    edm::Service<TFileService> fs;

    if( !fs ){
        throw edm::Exception( edm::errors::Configuration,
                            "TFile Service is not registered in cfg file" );
    }
    modelWeightSum_ = fs->make<TH1F>("modelWeightSum", "model weights", 2, 0., 2.);
    modelWeightSum_->GetXaxis()->SetBinLabel(1, "NEvents");
    modelWeightSum_->GetXaxis()->SetBinLabel(2, "Weighted");
    evwght_= fs->make<TH1F>("eventweight", "event weight", 1000, 0., 10);
    // quantities without weight
    pt_    = fs->make<TH1F>("pttop"        , "pttop"        , 80 , 0   , 800 );
    mass_  = fs->make<TH1F>("ttbarmass"    , "ttbarmass"    , 100, 0   , 2000);
    rap_   = fs->make<TH1F>("rapidity"     , "rapidity"     , 50 , -2.5, 2.5 );
    pttt_  = fs->make<TH1F>("ptttbar"      , "ptttbar"      , 120, 0   , 1200);
    raptt_ = fs->make<TH1F>("rapidityttbar", "rapidityttbar", 50 , -2.5, 2.5 );    
    // weight as function of the relevant quantity
    weightVsMass_ = fs->make<TH1F>("SFttbarmass", "SF(m^{t#bar{t}})", 1500, 300., 1800.);
    weightVsPt_   = fs->make<TH1F>("SFpttop"    , "SF(p_{T}^{top})" , 800 , 0.  , 800. );
    // quantities including weight
    pttopWeightedEvwgt_  = fs->make<TH1F>("pttopweighted"  , "p_{T}^{top} (event weight)"     , 80 , 0   , 800 );
    ytopWeightedEvwgt_   = fs->make<TH1F>("ytopweighted"   , "y_{top} (event weight)"         , 50 , -2.5, 2.5 );
    ptttbarWeightedEvwgt_= fs->make<TH1F>("ptttbarweighted", "p_{T}^{t#bar{t}} (event weight)", 120, 0   , 1200);
    yttbarWeightedEvwgt_ = fs->make<TH1F>("yttbarweighted" , "y^{t#bar{t}} (event weight)"    , 50 , -2.5, 2.5 );
    mttbarWeightedEvwgt_ = fs->make<TH1F>("mttbarweighted" , "m^{t#bar{t}} (event weight)"    , 100, 0   , 2000);
}


double EventWeightDileptonModelVariation::getTopPtWeight()
{
    // get weights for both tops
    double weight1 = (1+(genEvt->top()   ->pt()-weight1x_)*slope_);
    double weight2 = (1+(genEvt->topBar()->pt()-weight1x_)*slope_);
    // draw reweighting function
    if(!function){
      function=new TF1("function", "1+(x-[0])*[1]", 0., 800.);
      function->SetParameter(0, weight1x_);
      function->SetParameter(1, slope_  );
      weightVsPt_->Add((TF1*)function->Clone());
    }
    return weight1 * weight2;
}

double EventWeightDileptonModelVariation::getTopPtDataWeight()
{
  // get values
  double ptTop =genEvt->top()   ->pt();
  double ptATop=genEvt->topBar()->pt();
  // get weights for both tops
  double weightTop =TMath::Exp(slope_*(weight1x_-ptTop ));
  double weightATop=TMath::Exp(slope_*(weight1x_-ptATop));
  // draw reweighting function
  if(!function) {
    function=new TF1("function", "TMath::Exp([1]*([0]-x))", 0., 800.);
    function->SetParameter(0, weight1x_);
    function->SetParameter(1, slope_   );
    weightVsPt_->Add((TF1*)function->Clone());
  }
  return TMath::Sqrt(weightTop * weightATop);
}

double EventWeightDileptonModelVariation::getTopPtSoftenWeight()
{
  // get values
  double ptTop =genEvt->top()   ->pt();
  double ptATop=genEvt->topBar()->pt();
  // get weights for both tops
  double weightTop =TMath::Exp(weight1x_+slope_*ptTop );
  double weightATop=TMath::Exp(weight1x_+slope_*ptATop);
  // draw reweighting function
  if(!function) {
    function=new TF1("function", "TMath::Exp([0]+x*[1])", 0., 800.);
    function->SetParameter(0, weight1x_);
    function->SetParameter(1, slope_  );
    weightVsPt_->Add((TF1*)function->Clone());
  }
  return TMath::Sqrt(weightTop * weightATop);
}

double EventWeightDileptonModelVariation::getTopPtHardenWeight()
{
  // get values
  double ptTop =genEvt->top()   ->pt();
  double ptATop=genEvt->topBar()->pt();
  // get weights for both tops
  double weightTop =-1.*TMath::Exp(weight1x_+slope_*ptTop )+2.;
  double weightATop=-1.*TMath::Exp(weight1x_+slope_*ptATop)+2.;
  // draw reweighting function
  if(!function) {
    function=new TF1("function", "-1.*TMath::Exp([0]+x*[1])+2", 0., 800.);
    function->SetParameter(0, weight1x_);
    function->SetParameter(1, slope_  );
    weightVsPt_->Add((TF1*)function->Clone());
  }
  return TMath::Sqrt(weightTop * weightATop);
}

double EventWeightDileptonModelVariation::getTopPtTTsysWeight()
{
  // get values
  ROOT::Math::Boost CoMBoostGenTtbar((genEvt->top()->p4() + genEvt->topBar()->p4()).BoostToCM());
  reco::Particle::LorentzVector TopBoosted = genEvt->top()->p4();
  TopBoosted = CoMBoostGenTtbar(TopBoosted);
  double ptTop=TopBoosted.pt();
  // get weights for both tops
  double weightTop =TMath::Exp(slope_*(weight1x_-ptTop ));
  // draw reweighting function
  if(!function) {
    function=new TF1("function", "TMath::Exp([1]*([0]-x))", 0., 800.);
    function->SetParameter(0, weight1x_);
    function->SetParameter(1, slope_   );
    weightVsPt_->Add((TF1*)function->Clone());
  }
  return weightTop;
}

double EventWeightDileptonModelVariation::getTopEtaWeight()
{
    double weight = (1+(std::abs(genEvt->top()->eta())-weight1x_)*slope_) 
                  * (1+(std::abs(genEvt->topBar()->eta())-weight1x_)*slope_);
    return weight;
}

double EventWeightDileptonModelVariation::getTopRapidityWeight()
{
    double weight1 = (1+(std::abs(genEvt->top()   ->rapidity())-weight1x_)*slope_);
    double weight2 = (1+(std::abs(genEvt->topBar()->rapidity())-weight1x_)*slope_);
    return weight1 * weight2;
}

double EventWeightDileptonModelVariation::getTopMassWeight()
{
    double mass = (genEvt->top()->p4() + genEvt->topBar()->p4()).M();
    double weight = (1+(mass-weight1x_)*slope_);
    // draw reweighting function
    if(!function) {
      function=new TF1("function", "1+(x-[0])*[1]", 300., 1800.);
      function->SetParameter(0, weight1x_);
      function->SetParameter(1, slope_  );
      weightVsMass_->Add((TF1*)function->Clone());
    }
    return weight;
}

double EventWeightDileptonModelVariation::getTopMassWeightLandau()
{
    double mass = (genEvt->top()->p4() + genEvt->topBar()->p4()).M();
    double weight =
        (TMath::Landau(mass, landauMPV_ + landauMoveX_, landauSigma_) / 
         TMath::Landau(mass, landauMPV_, landauSigma_));
    return weight;
    
}

void EventWeightDileptonModelVariation::produce(edm::Event& evt, const edm::EventSetup& setup)
{
    std::auto_ptr<double> eventWeight(new double);
    evt.getByLabel(ttGenEvent_, genEvt );
    // get event weight
    if (genEvt.failedToGet() || weightVariable_.empty()) *eventWeight = 1;
    else if (!weightVariable_.compare("toppt")) *eventWeight = getTopPtWeight();
    else if (!weightVariable_.compare("topeta")) *eventWeight = getTopEtaWeight();
    else if (!weightVariable_.compare("toprapidity")) *eventWeight = getTopRapidityWeight();
    else if (!weightVariable_.compare("ttbarmass")) *eventWeight = getTopMassWeight();
    else if (!weightVariable_.compare("ttbarmasslandau")) *eventWeight = getTopMassWeightLandau();
    else if (!weightVariable_.compare("data")) *eventWeight = getTopPtDataWeight();
    else if (!weightVariable_.compare("topptTTsys")) *eventWeight = getTopPtTTsysWeight();
    else if (!weightVariable_.compare("topptHard")) *eventWeight = getTopPtSoftenWeight();
    else if (!weightVariable_.compare("topptSoft")) *eventWeight = getTopPtHardenWeight();

    else edm::LogError("EventWeightDilepton") << "don't know which variation to take, " 
        << weightVariable_ << " is invalid!" << std::endl;
    // consider minimum and maximum weights
    *eventWeight = std::min(std::max(*eventWeight, minWeight_), maxWeight_);
    // fill unweighted plots
    mass_ ->Fill((genEvt->top()->p4()+genEvt->topBar()->p4()).M()       );
    raptt_->Fill((genEvt->top()->p4()+genEvt->topBar()->p4()).Rapidity());
    pttt_ ->Fill((genEvt->top()->p4()+genEvt->topBar()->p4()).Pt()      );
    rap_ ->Fill(genEvt->top()   ->rapidity());
    rap_ ->Fill(genEvt->topBar()->rapidity());
    pt_  ->Fill(genEvt->top()   ->pt()      );
    pt_  ->Fill(genEvt->topBar()->pt()      );
    // fill weighted plots
    mttbarWeightedEvwgt_ ->Fill((genEvt->top()->p4()+genEvt->topBar()->p4()).M()       , *eventWeight);
    yttbarWeightedEvwgt_ ->Fill((genEvt->top()->p4()+genEvt->topBar()->p4()).Rapidity(), *eventWeight);
    ptttbarWeightedEvwgt_->Fill((genEvt->top()->p4()+genEvt->topBar()->p4()).Pt()      , *eventWeight);
    ytopWeightedEvwgt_ ->Fill(genEvt->top()   ->rapidity(), *eventWeight);
    ytopWeightedEvwgt_ ->Fill(genEvt->topBar()->rapidity(), *eventWeight);
    pttopWeightedEvwgt_->Fill(genEvt->top()   ->pt()      , *eventWeight);
    pttopWeightedEvwgt_->Fill(genEvt->topBar()->pt()      , *eventWeight);
    // fill weight plots to get average weight from
    modelWeightSum_->Fill(0);
    modelWeightSum_->Fill(1, *eventWeight);
    evwght_->Fill(*eventWeight);
//     std::cout << "putting " << *eventWeight << "\n";
    evt.put(eventWeight);
}

DEFINE_FWK_MODULE( EventWeightDileptonModelVariation );
