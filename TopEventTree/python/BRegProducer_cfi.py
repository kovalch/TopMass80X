import FWCore.ParameterSet.Config as cms

addBRegProducer = cms.EDProducer("BRegProducer",
    jets            = cms.InputTag("selectedPatJets"),
    rho_tag      = cms.InputTag('kt6PFJets','rho'),
    rho25_tag    = cms.InputTag('rho25kt6PFJets','rho'),
    mva_name     = cms.string("BDT::BDTG"),
#    mva_path     = cms.string("/afs/naf.desy.de/user/k/kirschen/public/TMVARegTEST_26BDTG_BRegJet_genPartonPt_BRegJet_jetPtCorr_BDTG.weights.xml"),
#    GBRmva_path     = cms.string("/afs/naf.desy.de/user/k/kirschen/public/fgbrtraintestBRegJet_genPartonPt_BRegJet_jetPtCorr.root"),
#    mva_path     = cms.string("/afs/naf.desy.de/user/k/kirschen/public/TMVARegTEST_26BDTG_BRegJet_genJetPt_BRegJet_jetPtCorr_BDTG.weights.xml"),
#    GBRmva_path     = cms.string("/afs/naf.desy.de/user/k/kirschen/public/fgbrtraintestBRegJet_genJetPt_BRegJet_jetPtCorr.root"),
    mva_path     = cms.string("/afs/naf.desy.de/user/k/kirschen/public/TMVARegTEST_26BDTG_BRegJet_genPartonPt_BRegJet_jetPtCorr_BDTG.weights.xml"),
    GBRmva_path     = cms.string("/afs/naf.desy.de/user/k/kirschen/public/fgbrtraintestBRegJet_genPartonPt_BRegJet_jetPtCorr.root"),
    writeOutVariables = cms.bool(True),
    gluonTagSrc  = cms.InputTag(""),
    maxNJets = cms.int32(20)
)
