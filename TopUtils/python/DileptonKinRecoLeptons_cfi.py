import FWCore.ParameterSet.Config as cms

dileptonKinRecoLeptons = cms.EDProducer('DileptonKinRecoLeptons',
    
    ## sources
    electrons = cms.InputTag("selectedPatElectrons"),
    muons = cms.InputTag("selectedPatMuons"),
    
    ## Require a specific channel
    ## options are "ee", "emu", "mumu"
    ## No selection for empty string ""
    filterChannel = cms.string(""),
    
    ## Requirements on dilepton mass
    ## Needs an even number of arguments, interpreted as intervals
    ## Empty vector means no selection
    excludeMasses = cms.vdouble(),
)
