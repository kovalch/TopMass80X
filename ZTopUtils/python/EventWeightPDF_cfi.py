import FWCore.ParameterSet.Config as cms

EventWeightPDF = cms.EDProducer("EventWeightPDF",
                                genEventInfoTag=cms.InputTag("generator"),
                                lheEventInfoTag=cms.InputTag("externalLHEProducer"),
                                beginWeightID=cms.string("2001"),
                                PDFSetNames=cms.string("CT10"),
                                printLHE=cms.bool(False)
                                )


