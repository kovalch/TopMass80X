#ifndef SignalProcessIDAnalyzer_h
#define SignalProcessIDAnalyzer_h

#include <TString.h>

#include "TH1.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

/**
   \class   SignalProcessIDAnalyzer SignalProcessIDAnalyzer.h "TopAnalysis/TopAnalyzer/plugins/SignalProcessIDAnalyzer.h"

   \brief   EDAnalyzer to to plot ID of process

   EDAnalyzer to to plot properties of primary vertex
*/

using namespace std;
using namespace edm;

class SignalProcessIDAnalyzer : public EDAnalyzer {

  public:
    /// default constructor
    explicit SignalProcessIDAnalyzer(const ParameterSet&);
    /// default destructor
    ~SignalProcessIDAnalyzer();
    
  private:
    /// initiate n_TrigPaths
    virtual void beginJob();
    /// look which triggers have fired and compare to given set of triggers
    virtual void analyze(const Event&, const EventSetup&);
    /// empty
    virtual void endJob();
        
    TH1D* id_;    
};

#endif
