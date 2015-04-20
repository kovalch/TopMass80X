#ifndef TriggerFilter_h
#define TriggerFilter_h

#include <string>
#include <vector>
#include <iostream>

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDMException.h"

/**
   \class   TriggerFilter TriggerFilter.h "TopAnalysis/TopFilter/plugins/TriggerFilter.h"

   \brief   EDFilter for trigger

   This filter filters for triggers given in the config file.
*/

class TriggerFilter : public edm::EDFilter {

public:
    /// default constructor
    explicit TriggerFilter(const edm::ParameterSet&);
    /// default destructor
    ~TriggerFilter();

private:
    /// initiate n_TrigPaths
    virtual void beginJob();
    /// look which triggers have fired and compare to given set of triggers
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    /// empty
    virtual void endJob();

    /// triger result input collection
    edm::InputTag trigResults_;
    /// triggers to be selected given in config
    std::vector<std::string> hltPaths_;
    /// triggers which may not have fired
    std::vector<std::string> vetoPaths_;
    /// number of trigger paths given in config
    int n_TrigPaths;
    /// number of forbidden trigger paths given in config
    int n_VetoPaths;
    /// print list of available triggers including fired/not fired info
    bool printTriggers_;

};

#endif
