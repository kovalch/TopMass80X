import FWCore.ParameterSet.Config as cms

#
# module to add extra information for b-regression studies
#
analyzeBRegJets = cms.EDAnalyzer("BRegJetEventAnalyzer",
    jets         = cms.InputTag("goodJetsPF30"),
    #allJets      = cms.InputTag("patJets"),
    #noPtEtaJets  = cms.InputTag("noPtEtaJetsPF"),
    rho_tag      = cms.InputTag('kt6PFJetsUsePATWithPFlow','rho'),
    rho25_tag    = cms.InputTag('rho25kt6PFJets','rho'),
    mva_name     = cms.string("BDT::BDTG"),
    mva_path     = cms.string("/afs/naf.desy.de/user/k/kirschen/scratch/2013_05_TopMass/CMSSW_5_3_9/src/TopMass/TopEventTree/plugins/TMVARegression_BDTG.weights.xml"),
    writeOutVariables = cms.bool(True),
    gluonTagSrc  = cms.InputTag(""),
    maxNJets = cms.int32(20)
)
