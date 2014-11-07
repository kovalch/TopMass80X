import FWCore.ParameterSet.Config as cms

writeNTuple = cms.EDAnalyzer('NTupleWriter',
    
    # Specification of analysed sample
    sampleName = cms.string("please give sampleName"),
    channelName = cms.string("please give channelName"),
    systematicsName = cms.string("please give systematicsName"),
    isMC = cms.bool(True),
    isTtbarSample = cms.bool(True),
    isHiggsSample = cms.bool(False),
    isZSample = cms.bool(False),
    isMadgraphSample = cms.bool(False),
    
    # Bools for inclusion/exclusion of specific info
    includeTrigger = cms.bool(False),
    includePdfWeights = cms.bool(False),
    saveHadronMothers = cms.bool(False),
    saveCHadronParticles = cms.bool(False),
    
    # Input sources for reco-level info
    electrons = cms.InputTag("electrons"),
    muons = cms.InputTag("muons"),
    jets = cms.InputTag("jets"),
    jetsForMet = cms.InputTag("jetsForMet"),
    jetsForMetUncorrected = cms.InputTag("jetsForMetUncorrected"),
    jetProperties = cms.InputTag("jetProperties"),
    met = cms.InputTag("met"),
    mvaMet  = cms.InputTag("mvaMet"),
    vertices = cms.InputTag('goodOfflinePrimaryVertices'),
    triggerResults = cms.InputTag('TriggerResults', '', 'HLT'),
    
    # Input sources for gen-level info
    pileupInfo = cms.InputTag("addPileupInfo"),
    genParticles = cms.InputTag("genParticles"),
    genJets = cms.InputTag("ak5GenJets", "", "SIM"),
    pdfWeights = cms.InputTag('pdfWeights:cteq66'),
    genEventTtbar = cms.InputTag("genEvt"),
    genEventHiggs = cms.InputTag("genEvtHiggs"),
    genZDecay = cms.InputTag("genZDecay"),
    BHadJetIndex = cms.InputTag("produceGenLevelBJets", "BHadJetIndex"),
    AntiBHadJetIndex = cms.InputTag("produceGenLevelBJets", "AntiBHadJetIndex"),
    BHadrons = cms.InputTag("produceGenLevelBJets", "BHadrons"),
    AntiBHadrons = cms.InputTag("produceGenLevelBJets", "AntiBHadrons"),
    BHadronFromTopB = cms.InputTag("produceGenLevelBJets", "BHadronFromTopB"),
    AntiBHadronFromTopB = cms.InputTag("produceGenLevelBJets", "AntiBHadronFromTopB"),
    BHadronVsJet = cms.InputTag("produceGenLevelBJets", "BHadronVsJet"),
    AntiBHadronVsJet = cms.InputTag("produceGenLevelBJets", "AntiBHadronVsJet"),
    genBHadPlusMothers = cms.InputTag("matchGenHFHadronJets", "genBHadPlusMothers"),
    genBHadPlusMothersIndices = cms.InputTag("matchGenHFHadronJets", "genBHadPlusMothersIndices"),
    genBHadIndex = cms.InputTag("matchGenHFHadronJets", "genBHadIndex"),
    genBHadFlavour = cms.InputTag("matchGenHFHadronJets", "genBHadFlavour"),
    genBHadJetIndex = cms.InputTag("matchGenHFHadronJets", "genBHadJetIndex"),
    genBHadLeptonIndex = cms.InputTag("matchGenHFHadronJets", "genBHadLeptonIndex"),
    genBHadLeptonHadronIndex = cms.InputTag("matchGenHFHadronJets", "genBHadLeptonHadronIndex"),
    genBHadLeptonViaTau = cms.InputTag("matchGenHFHadronJets", "genBHadLeptonViaTau"),
    genBHadFromTopWeakDecay = cms.InputTag("matchGenHFHadronJets", "genBHadFromTopWeakDecay"),
    genCHadPlusMothers = cms.InputTag("matchGenHFHadronJets", "genCHadPlusMothers"),
    genCHadPlusMothersIndices = cms.InputTag("matchGenHFHadronJets", "genCHadPlusMothersIndices"),
    genCHadIndex = cms.InputTag("matchGenHFHadronJets", "genCHadIndex"),
    genCHadJetIndex = cms.InputTag("matchGenHFHadronJets", "genCHadJetIndex"),
    genCHadLeptonIndex = cms.InputTag("matchGenHFHadronJets", "genCHadLeptonIndex"),
    genCHadLeptonHadronIndex = cms.InputTag("matchGenHFHadronJets", "genCHadLeptonHadronIndex"),
    genCHadLeptonViaTau = cms.InputTag("matchGenHFHadronJets", "genCHadLeptonViaTau"),
    genCHadFromTopWeakDecay = cms.InputTag("matchGenHFHadronJets", "genCHadFromTopWeakDecay"),
    genCHadBHadronId = cms.InputTag("matchGenHFCHadronJets", "genCHadBHadronId"),
    ttbarDecayMode = cms.InputTag("generatorTopFilter", "decayMode"),
    higgsDecayMode = cms.InputTag("generatorHiggsFilter", "higgsDecayMode"),
    madgraphWDecay = cms.InputTag("madgraphWDecayProducer", "madgraphWDecay"),
)

