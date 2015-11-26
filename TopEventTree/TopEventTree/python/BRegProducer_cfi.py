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
#    mva_path     = cms.string("/afs/naf.desy.de/user/k/kirschen/public/TMVARegTEST_26BDTG_BRegJet_genPartonPt_BRegJet_jetPtCorr_BDTG.weights.xml"),
#    GBRmva_path     = cms.string("/afs/naf.desy.de/user/k/kirschen/public/fgbrtraintestBRegJet_genPartonPt_BRegJet_jetPtCorr.root"),
#    mva_path     = cms.string("/afs/naf.desy.de/user/k/kirschen/public/TMVARegTEST5311SemiLept_26BDTG_BRegJet_genJetPt_BRegJet_jetPtCorr_BDTG.weights.xml"),
#    GBRmva_path     = cms.string("/afs/naf.desy.de/user/k/kirschen/public/fgbrtraintestSemiLeptBRegJet_genJetPt_BRegJet_jetPtCorr_26.root"),
    mva_path     = cms.string("/afs/naf.desy.de/user/k/kirschen/public/_scratch_hh_dust_naf_cms_user_kirschen_BRegression_GC_TMVATreeOutput_TMVARegTEST5311QGPUIdSmallTrainingSemiLept_24BDTG_BRegJet_genJetPt_BRegJet_jetPtCorr_BDTG.weights.xml"),
#    GBRmva_path     = cms.string("/afs/naf.desy.de/user/k/kirschen/public/fgbrtraintest100SemiLeptPtFlatWeightPt300CutBRegJet_genJetPt_BRegJet_jetPtCorr_24.root"),
    GBRmva_path     = cms.string("/afs/naf.desy.de/user/k/kirschen/public/fgbrtraintest1000TresMin00250MinCutSig3TransitQ07SemiLeptAllMassAllLargeAllNewBRegSamplePtFlatWeightPt300CutBRegJet_genJetPt_BRegJet_jetPtCorr_24.root"),
    writeOutVariables = cms.bool(True),
    JECUncSrcFile     = cms.FileInPath("TopAnalysis/TopUtils/data/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt"),
    gluonTagSrc  = cms.InputTag(""),
    maxNJets = cms.int32(20)
)
