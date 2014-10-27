# analog to TopAnalysis/TopFilters/python/filters/GeneratorTopFilter_cfi.py

import FWCore.ParameterSet.Config as cms

generatorHiggsFilter = cms.EDFilter('GeneratorHiggsFilter',
    # input subset containing relevant Higgs information
    src = cms.InputTag("genEvtHiggs"),
    
    ## SHORT CUTS: 
    # if not empty all the other boolean parameters
    # will be ignored (except invert selection)
    channels = cms.vstring(),
    # supported: bb, tautau, WW, ZZ, none
    
    # selection can be inverted
    invert_selection = cms.bool(False),
)
