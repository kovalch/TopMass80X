#ifndef EventTreeCreator_h
#define EventTreeCreator_h

#include "TTree.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"

class EventTreeCreator : public edm::EDAnalyzer {

 public:

  explicit EventTreeCreator(const edm::ParameterSet&);
  ~EventTreeCreator();
  
 private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
  TTree* eventTree;
};

#endif
