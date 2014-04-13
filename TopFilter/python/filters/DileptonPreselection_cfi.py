import FWCore.ParameterSet.Config as cms

## Preselect events which have at least one dilepton pair fulfilling the criteria specified in the filter
dileptonPreselection = cms.EDFilter("DileptonPreselection",

    ## sources
    electrons = cms.InputTag("selectedPatElectrons"),
    muons = cms.InputTag("selectedPatMuons"),

    ## Requirements on charge of lepton pair
    ## -1 for unlike, 1 for like sign
    ## 0: no selection
    filterCharge = cms.int32(0),
    
    ## Require a specific channel
    ## options are "ee", "emu", "mumu"
    ## No selection for empty string ""
    filterChannel = cms.string(""),
    
    ## Requirements on dilepton mass
    ## Needs an even number of arguments, interpreted as intervals
    ## Empty vector means no selection
    excludeMasses = cms.vdouble(),
)


