#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <string>
#include <TH1D.h>

class CountEventAnalyzer : public edm::EDAnalyzer {

public:

    CountEventAnalyzer(const edm::ParameterSet&);
    ~CountEventAnalyzer();

protected:
    void beginJob();
    void analyze(const edm::Event&, const edm::EventSetup&);
    void endJob();

private:
    TH1* unweightedEvents_;
    TH1* weightedEvents_;
    TH1* pdfEvents_;
    bool includepdfweights_;
    edm::InputTag pdfweights_;
};

CountEventAnalyzer::CountEventAnalyzer( const edm::ParameterSet& iConfig ) :
    includepdfweights_(iConfig.getParameter<bool>("includePDFWeights")),
    pdfweights_(iConfig.getParameter<edm::InputTag>("pdfWeights"))
{
}

CountEventAnalyzer::~CountEventAnalyzer() {
}

void CountEventAnalyzer::beginJob() {
    edm::Service<TFileService> fs;
    if( !fs ){
        throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );
    }
    unweightedEvents_ = fs->make<TH1I>( "unweightedEvents", "N of Events", 1,0.5,1.5);
    weightedEvents_ = fs->make<TH1D>( "weightedEvents", "N of weighted events", 1,0.5,1.5);
    pdfEvents_ = fs->make<TH1D>("pdfEventWeights", "PDF event weights", 200, -0.5, 199.5);
}

void CountEventAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& context) {
    unweightedEvents_->Fill(1);
    if (iEvent.isRealData()) return; //we are done for data
    
    edm::InputTag genEventInfoTag("generator");
    edm::Handle<GenEventInfoProduct> evt_info;
    iEvent.getByLabel(genEventInfoTag, evt_info);
    weightedEvents_->Fill(1, evt_info->weight());
        
    if (includepdfweights_) {
        edm::Handle<std::vector<double> > weightHandle;
        iEvent.getByLabel(pdfweights_, weightHandle);
        for (size_t i = 0; i < (*weightHandle).size(); ++i) {
            pdfEvents_->Fill(i, weightHandle->at(i));
        }
    }
}

void CountEventAnalyzer::endJob() {
}

DEFINE_FWK_MODULE( CountEventAnalyzer );
