//#include <memory>

// user include files
//#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TopMass/TopEventTree/interface/TreeRegistryService.h"

#include "AnalysisDataFormats/TopObjects/interface/TtSemiLepEvtPartons.h"
#include "AnalysisDataFormats/TopObjects/interface/TtFullHadEvtPartons.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"

#include "TopMass/TopEventTree/plugins/DumpEvent.h"


void DumpEvent::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
	std::vector<edm::Provenance const*> plist;
	event.getAllProvenance(plist);

	for (std::vector<edm::Provenance const*>::const_iterator prov = plist.begin(); prov < plist.end(); ++prov)
		std::cout << (*prov)->className() << " " << (*prov)->branchName() << " " << (*prov)->moduleLabel() << " " << (*prov)->productInstanceName() << " " << (*prov)->processName() << std::endl;
}

void DumpEvent::endRun(edm::Run const& run, edm::EventSetup const&)
{
	std::vector<edm::Provenance const*> plist;
	run.getAllProvenance(plist);

	for (std::vector<edm::Provenance const*>::const_iterator prov = plist.begin(); prov < plist.end(); ++prov)
		std::cout << (*prov)->className() << " " << (*prov)->branchName() << " " << (*prov)->moduleLabel() << " " << (*prov)->productInstanceName() << " " << (*prov)->processName() << std::endl;
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DumpEvent);
