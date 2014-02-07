import FWCore.ParameterSet.Config as cms

#
# This small module can combine bools from the event content
# It returns true if all (if any) sources from "mustBeTrue" are true
# and all sources (if any) from "mustBeFalse" are false
# Else it returns false
# 
# This is NOT a filter itself so it does not stop events
# from being processed
#
combinedbools = cms.EDProducer("Boolcombiner",
                               mustBeFalse = cms.VInputTag(),
                               mustBeTrue = cms.VInputTag(),
                               debug=cms.bool(False)
)
