#ifndef SherpaGenEventAnalyzer_h
#define SherpaGenEventAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "TopMass/TopEventTree/interface/SherpaGenEvent.h"

//#include <utility>

class SherpaGenEventAnalyzer : public edm::EDAnalyzer {

  public:

  explicit SherpaGenEventAnalyzer(const edm::ParameterSet&);
  ~SherpaGenEventAnalyzer();

 private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  edm::Service<TreeRegistryService> trs;

  edm::InputTag genJets_;

  // THE SherpaGenEvent to store the information
  SherpaGenEvent* sherpa;

};

#endif
