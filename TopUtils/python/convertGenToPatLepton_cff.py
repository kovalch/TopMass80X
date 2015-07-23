import FWCore.ParameterSet.Config as cms

recoMET = cms.EDProducer("MyRecoMETProducer",
                           metSrc = cms.InputTag("genMetTrue"),
                           )

patMETs = cms.EDProducer("PATMETProducer",
    # input
    metSource  = cms.InputTag("recoMET"),

    # add user data
    userData = cms.PSet(
      # add custom classes here
      userClasses = cms.PSet(
        src = cms.VInputTag('')
      ),
      # add doubles here
      userFloats = cms.PSet(
        src = cms.VInputTag('')
      ),
      # add ints here
      userInts = cms.PSet(
        src = cms.VInputTag('')
      ),
      # add candidate ptrs here
      userCands = cms.PSet(
        src = cms.VInputTag('')
      ),
      # add "inline" functions here
      userFunctions = cms.vstring(),
      userFunctionLabels = cms.vstring()
    ),

    # muon correction
    addMuonCorrections = cms.bool(False),
    muonSource         = cms.InputTag(""),

    # mc matching configurables
    addGenMET    = cms.bool(False),
    genMETSource = cms.InputTag(""),

    # efficiencies
    addEfficiencies = cms.bool(False),
    efficiencies    = cms.PSet(),

    # resolution
    addResolutions  = cms.bool(False),
    resolutions     = cms.PSet(),
)

scaledMET = cms.EDFilter("PATMETSelector",
    src = cms.InputTag("patMETs"),
    cut = cms.string("")
)

genMuons = cms.EDFilter("CandViewSelector",
                       src = cms.InputTag("genParticles"),
                       cut = cms.string("status = 1 & abs(pdgId) = 13 & abs(eta) < 2.5 & pt > 20.")
                       )

recoMuons = cms.EDProducer("MyRecoMuonProducer",
                           muonSrc = cms.InputTag("genMuons"),
                           )

patMuons = cms.EDProducer("PATMuonProducer",
    # input
    muonSource      = cms.InputTag("recoMuons"),

    # use particle flow instead of std reco
    useParticleFlow =  cms.bool( False ),
    pfMuonSource    = cms.InputTag(""),

    # add user data
    userData = cms.PSet(
      # add custom classes here
      userClasses = cms.PSet(
        src = cms.VInputTag('')
      ),
      # add doubles here
      userFloats = cms.PSet(
        src = cms.VInputTag('')
      ),
      # add ints here
      userInts = cms.PSet(
        src = cms.VInputTag('')
      ),
      # add candidate ptrs here
      userCands = cms.PSet(
        src = cms.VInputTag('')
      ),
      # add "inline" functions here
      userFunctions = cms.vstring(),
      userFunctionLabels = cms.vstring()
    ),

    # embedding objects
    embedMuonBestTrack      = cms.bool(False),  ## embed in AOD externally stored muon best track from global pflow
    embedTunePMuonBestTrack = cms.bool(False),  ## embed in AOD externally stored muon best track from muon only
    forceBestTrackEmbedding = cms.bool(False), ## force embedding separately the best tracks even if they're already embedded e.g. as tracker or global tracks
    embedTrack          = cms.bool(False),  ## embed in AOD externally stored tracker track
    embedCombinedMuon   = cms.bool(False),  ## embed in AOD externally stored combined muon track
    embedStandAloneMuon = cms.bool(False),  ## embed in AOD externally stored standalone muon track
    embedPickyMuon      = cms.bool(False),  ## embed in AOD externally stored TeV-refit picky muon track
    embedTpfmsMuon      = cms.bool(False),  ## embed in AOD externally stored TeV-refit TPFMS muon track
    embedDytMuon        = cms.bool(False),  ## embed in AOD externally stored TeV-refit DYT muon track
    embedPFCandidate    = cms.bool(False),  ## embed in AOD externally stored particle flow candidate

    # embedding of muon MET corrections for caloMET
    embedCaloMETMuonCorrs = cms.bool(False),
    caloMETMuonCorrs = cms.InputTag(""),
    # embedding of muon MET corrections for tcMET
    embedTcMETMuonCorrs   = cms.bool(False), # removed from RECO/AOD!
    tcMETMuonCorrs   = cms.InputTag(""),

    # embed IsoDeposits
    isoDeposits = cms.PSet(
    ),

    userIsolation = cms.PSet(
    ),

    # mc matching
    addGenMatch   = cms.bool(False),
    embedGenMatch = cms.bool(False),
    genParticleMatch = cms.InputTag(""), ## particles source to be used for the matching

    # efficiencies
    addEfficiencies = cms.bool(False),
    efficiencies    = cms.PSet(),

    # resolution configurables
    addResolutions  = cms.bool(False),
    resolutions      = cms.PSet(),

    # high level selections
    embedHighLevelSelection = cms.bool(False),
    usePV                   = cms.bool(False),
    beamLineSrc             = cms.InputTag(""),
    pvSrc                   = cms.InputTag("")
)

looseMuons= cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string('')
)
tightMuons= cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("pt > 33 & abs(eta) < 2.1")
)

genElectrons = cms.EDFilter("CandViewSelector",
                       src = cms.InputTag("genParticles"),
                       cut = cms.string("status = 1 & abs(pdgId) = 11 & pt > 20 & abs(eta) < 2.5")
                       )

recoElectrons = cms.EDProducer("MyRecoMuonProducer",
                           muonSrc = cms.InputTag("genElectrons"),
                           )

patElectrons = cms.EDProducer("PATMuonProducer",
    # input
    muonSource      = cms.InputTag("recoElectrons"),

    # use particle flow instead of std reco
    useParticleFlow =  cms.bool( False ),
    pfMuonSource    = cms.InputTag(""),

    # add user data
    userData = cms.PSet(
      # add custom classes here
      userClasses = cms.PSet(
        src = cms.VInputTag('')
      ),
      # add doubles here
      userFloats = cms.PSet(
        src = cms.VInputTag('')
      ),
      # add ints here
      userInts = cms.PSet(
        src = cms.VInputTag('')
      ),
      # add candidate ptrs here
      userCands = cms.PSet(
        src = cms.VInputTag('')
      ),
      # add "inline" functions here
      userFunctions = cms.vstring(),
      userFunctionLabels = cms.vstring()
    ),

    # embedding objects
    embedMuonBestTrack      = cms.bool(False),  ## embed in AOD externally stored muon best track from global pflow
    embedTunePMuonBestTrack = cms.bool(False),  ## embed in AOD externally stored muon best track from muon only
    forceBestTrackEmbedding = cms.bool(False), ## force embedding separately the best tracks even if they're already embedded e.g. as tracker or global tracks
    embedTrack          = cms.bool(False),  ## embed in AOD externally stored tracker track
    embedCombinedMuon   = cms.bool(False),  ## embed in AOD externally stored combined muon track
    embedStandAloneMuon = cms.bool(False),  ## embed in AOD externally stored standalone muon track
    embedPickyMuon      = cms.bool(False),  ## embed in AOD externally stored TeV-refit picky muon track
    embedTpfmsMuon      = cms.bool(False),  ## embed in AOD externally stored TeV-refit TPFMS muon track
    embedDytMuon        = cms.bool(False),  ## embed in AOD externally stored TeV-refit DYT muon track
    embedPFCandidate    = cms.bool(False),  ## embed in AOD externally stored particle flow candidate

    # embedding of muon MET corrections for caloMET
    embedCaloMETMuonCorrs = cms.bool(False),
    caloMETMuonCorrs = cms.InputTag(""),
    # embedding of muon MET corrections for tcMET
    embedTcMETMuonCorrs   = cms.bool(False), # removed from RECO/AOD!
    tcMETMuonCorrs   = cms.InputTag(""),

    # embed IsoDeposits
    isoDeposits = cms.PSet(
    ),

    userIsolation = cms.PSet(
    ),

    # mc matching
    addGenMatch   = cms.bool(False),
    embedGenMatch = cms.bool(False),
    genParticleMatch = cms.InputTag(""), ## particles source to be used for the matching

    # efficiencies
    addEfficiencies = cms.bool(False),
    efficiencies    = cms.PSet(),

    # resolution configurables
    addResolutions  = cms.bool(False),
    resolutions      = cms.PSet(),

    # high level selections
    embedHighLevelSelection = cms.bool(False),
    usePV                   = cms.bool(False),
    beamLineSrc             = cms.InputTag(""),
    pvSrc                   = cms.InputTag("")
      )

looseElectrons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patElectrons"),
    cut = cms.string('')
)
tightElectronsEJ = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patElectrons"),
    cut = cms.string("pt > 33 & abs(eta) < 2.1")
)

convertGenToPatLeptons = cms.Sequence(recoMET*patMETs*genMuons*recoMuons*patMuons*looseMuons*tightMuons*genElectrons*recoElectrons*patElectrons*looseElectrons*tightElectronsEJ)

