#ifndef JetEventAnalyzer_h
#define JetEventAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "TopMass/TopEventTree/interface/JetEvent.h"


class JetEventAnalyzer : public edm::EDAnalyzer {

 public:

  explicit JetEventAnalyzer(const edm::ParameterSet&);
  ~JetEventAnalyzer();

 private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  edm::Service<TreeRegistryService> trs;

  edm::InputTag jets_;
  //edm::InputTag allJets_;
  //edm::InputTag noPtEtaJets_;
  edm::InputTag gluonTagSrc_;

  // max possible number of jets in events
  const int kJetMAX_;

  // THE JetEvent to store the information
  JetEvent* jet;

  // check only once per module run if the needed collections are available
  bool checkedJERSF, checkedJESSF, checkedTotalSF, checkedQGTag;
  bool     hasJERSF,     hasJESSF,     hasTotalSF,     hasQGTag;

};

#endif
