#ifndef GeneratorWNJetsFilter_h
#define GeneratorWNJetsFilter_h

#include <TString.h>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class GeneratorWNJetsFilter : public edm::EDFilter {

 public:

  explicit GeneratorWNJetsFilter(const edm::ParameterSet&);
  ~GeneratorWNJetsFilter(){};
  
 private:

  virtual void beginJob();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

 private:
  
  int NJet_;
};

#endif
