import FWCore.ParameterSet.Config as cms

#
# service to handle the event tree for mutiple modules
#

TreeRegistryService = cms.Service("TreeRegistryService",
    treeName  = cms.string("eventTree"),
    treeTitle = cms.string("Tree for UHH top-quark analysis\nParticles are in order {TTBar, HadTop, LepTop, HadW, LepW, HadB, LightQ, LightQBar, LepB, Lepton, Neutrino}")
)