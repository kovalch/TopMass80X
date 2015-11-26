#ifndef EventHypothesisAnalyzer_h
#define EventHypothesisAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiLeptonicEvent.h"
#include "AnalysisDataFormats/TopObjects/interface/TtFullHadronicEvent.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "TopMass/TopEventTree/interface/TopEvent.h"


class DumpEvent : public edm::EDAnalyzer {
 public:
  explicit DumpEvent(const edm::ParameterSet&) {}
  ~DumpEvent() {}
 private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
};

#endif
