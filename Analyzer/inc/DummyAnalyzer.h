#include "MassAnalyzer.h"

class DummyAnalyzer : public MassAnalyzer {
  public:
    DummyAnalyzer(TChain* chain) : MassAnalyzer(chain) {};
    void Analyze(TString cuts);
};
