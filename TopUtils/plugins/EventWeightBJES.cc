#include "FWCore/Utilities/interface/EDMException.h"
#include "TopAnalysis/TopUtils/plugins/EventWeightBJES.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Common/interface/View.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#define IS_BHADRON_PDGID(id) ( ((abs(id)/100)%10 == 5) || (abs(id) >= 5000 && abs(id) <= 5999) )
#define IS_NEUTRINO_PDGID(id) ( (abs(id) == 12) || (abs(id) == 14) || (abs(id) == 16) )

EventWeightBJES::EventWeightBJES(const edm::ParameterSet& cfg):
  genJets_                (cfg.getParameter<edm::InputTag>("genJets")),
  genParticles_           (cfg.getParameter<edm::InputTag>("genParticles")),
  nuDecayFractionSource_  (cfg.getParameter<double>("nuDecayFractionSource")),
  nuDecayFractionTarget_  (cfg.getParameter<double>("nuDecayFractionTarget")),
  fragSourceFile_         (cfg.getParameter<edm::FileInPath>("fragSourceFile")),
  fragTargetFile_         (cfg.getParameter<edm::FileInPath>("fragTargetFile"))
{
  consumes<std::vector< reco::GenJet > >(genJets_);
  consumes<reco::GenParticleCollection>(genParticles_);
  produces<double>();

  // load TFile Service
  edm::Service<TFileService> fs;
  if( !fs ){
    throw edm::Exception( edm::errors::Configuration,
		    "TFile Service is not registered in cfg file");
  }
  /// booking of histograms
  hists["eventWeight"]           = fs->make<TH1F>("eventWeight", "eventWeight", 1000, 0, 10 );
  hists["genBHadronNuDecay"]     = fs->make<TH1F>("genBHadronNuDecay", "genBHadronNuDecay", 2, 0, 2);
  hists["genBHadronPtFraction"]  = fs->make<TH1F>("genBHadronPtFraction", "genBHadronPtFraction", 100, 0, 2);

  sourceFile = new TFile((TString)fragSourceFile_.fullPath());
  targetFile = new TFile((TString)fragTargetFile_.fullPath());

  sourceHist = (TH1F*) sourceFile->Get("EventWeightBJES/genBHadronPtFraction")->Clone();
  targetHist = (TH1F*) targetFile->Get("EventWeightBJES/genBHadronPtFraction")->Clone();

  if (sourceHist->GetNbinsX() != targetHist->GetNbinsX()) std::cout << "Incompatible b-fragmentation histograms: Number of bins not equal" << std::endl;

  sourceHist->Scale(1./sourceHist->Integral());
  targetHist->Scale(1./targetHist->Integral());

  hists["weightHist"] = (TH1F*) targetHist->Clone();
  hists.find("weightHist")->second->Divide(sourceHist);

  std::cout << "Weights for b-fragmentation" << std::endl;
  for (int i = 0; i < sourceHist->GetNbinsX(); ++i) {
    std::cout << std::setiosflags(std::ios::left)
              << "[" << i << "] "
              << std::setw(10) << hists.find("weightHist")->second->GetBinContent(i);
  }
  std::cout << std::endl;
}

EventWeightBJES::~EventWeightBJES()
{
}

void
EventWeightBJES::produce(edm::Event& evt, const edm::EventSetup& setup)
{
  edm::Handle<std::vector< reco::GenJet > > genJets;
  evt.getByLabel(genJets_, genJets);

  double eventWeight = 1.;

  //////////////////////////////////////////////////////////////////////////
  // GENPARTICLES
  ////////////////////////////////////////////////////////////////////////

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByLabel(genParticles_, genParticles);
  for(size_t i = 0; i < genParticles->size(); ++ i) {
    const reco::GenParticle & p = (*genParticles)[i];
    if (p.pt() == 0) continue;

    int id = p.pdgId();
    if (!IS_BHADRON_PDGID(id)) continue;

    int n = p.numberOfDaughters();

    bool hasBDaughter = false;
    bool hasNuDaughter = false;
    double genBHadronNuDecay = 0.;

    for(int j = 0; j < n; ++j) {
      const reco::Candidate * d = p.daughter( j );
      int dauId = d->pdgId();
      if (IS_BHADRON_PDGID(dauId)) {
        hasBDaughter = true;
        break;
      }
      if (IS_NEUTRINO_PDGID(dauId)) hasNuDaughter = true;
    }

    // Weakly decaying B hadron
    if (!hasBDaughter) {

      // Neutrino fraction weight
      if (hasNuDaughter) {
        genBHadronNuDecay = 1.;
        eventWeight *= nuDecayFractionTarget_/nuDecayFractionSource_;
      }
      else {
        eventWeight *= (1.-nuDecayFractionTarget_)/(1.-nuDecayFractionSource_);
      }
      hists.find("genBHadronNuDecay")->second->Fill( genBHadronNuDecay );

      // Fragmentation weight
      for (std::vector< reco::GenJet >::const_iterator ijet = genJets->begin(); ijet != genJets->end(); ++ijet) {
        if (p.pt() == 0 || ijet->pt() == 0) continue;
        double deta = p.eta() - ijet->eta();
        double dphi = TVector2::Phi_mpi_pi(p.phi() - ijet->phi());
        double dr   = sqrt( deta*deta + dphi*dphi );

        // Simple dR match of hadron and GenJet
        if (dr < 0.5) {
          double xb = p.pt()/ijet->pt();
          hists.find("genBHadronPtFraction")->second->Fill(xb);
          if (xb < 2.) {
            eventWeight *= hists.find("weightHist")->second->GetBinContent(
              hists.find("weightHist")->second->GetXaxis()->FindFixBin(xb)
            );
          }
          break;
        }
      }
    }
  }

  std::auto_ptr<double> EventWeightBJES(new double);
  *EventWeightBJES = eventWeight;
  evt.put(EventWeightBJES);

  hists.find("eventWeight" )->second->Fill( eventWeight );
}

// executed at the end after looping over all events
void EventWeightBJES::endJob()
{
  double genBHadronNuDecayMean = hists.find("genBHadronNuDecay" )->second->GetMean();
  std::cout << "B hadron semi-leptonic BR = " << genBHadronNuDecayMean << std::endl;

  double genBHadronPtFractionMean = hists.find("genBHadronPtFraction" )->second->GetMean();
  std::cout << "B hadron pT fraction = " << genBHadronPtFractionMean << std::endl;
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( EventWeightBJES );
