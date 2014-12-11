import FWCore.ParameterSet.Config as cms

generatorTopFilter = cms.EDFilter("GeneratorTopFilter_Pythia8",

    src = cms.InputTag("genEvt"), 

    #how many event do you want to select?
    n_events = cms.int32(-1),
    
    rejectNonBottomDecaysOfTops = cms.bool(True),
    
    ## SHORT CUTS: 
    # if not empty all the other boolean parameters
    # will be ignored (except invert selection)
    channels = cms.vstring(),
    # supported: MuonMuon
    #            ElectronMuon
    #            ElectronElectron
    #            MuonMuonViaTau
    #            ElectronMuonViaTau
    #            ElectronElectronViaTau     
    #            SingleMuon
    #            SingleElectron  
    #            FullHadronic  
    
    ## MANUAL CONFIGURATION
    SingleTop_Hadronic = cms.bool(False),
    SingleTop_Electron = cms.bool(False),
    SingleTop_Muon     = cms.bool(False),
    SingleTop_Tau      = cms.bool(False),
    
    Ttbar_Hadronic = cms.bool(False),
    
    Ttbar_SemiElec = cms.bool(False),
    Ttbar_SemiMuon = cms.bool(False),
    Ttbar_SemiTau  = cms.bool(False),
    
    Ttbar_ElecElec = cms.bool(False),
    Ttbar_MuonMuon = cms.bool(False),
    Ttbar_TauTau   = cms.bool(False),
        
    Ttbar_ElecMuon = cms.bool(False),    
    Ttbar_ElecTau  = cms.bool(False),         
    Ttbar_MuonTau  = cms.bool(False),

    # if a channel is selected which contains
    # a tau you can choose a specific further
    # tau decay
    Tau_To_Hadron   = cms.bool(False), 
    Tau_To_Electron = cms.bool(False),       
    Tau_To_Muon     = cms.bool(False),
    
    # selection can be inverted
    invert_selection = cms.bool(False)   
)
