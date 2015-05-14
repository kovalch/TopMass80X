import FWCore.ParameterSet.Config as cms

EventWeightPDF = cms.EDProducer("EventWeightPDF",
                                lheEventInfoTag=cms.InputTag("externalLHEProducer"),
                                PDFName=cms.string("CT10"),
                                nominalWeightID=cms.string("1001"),
                                beginWeightID=cms.string("2001"),
                                endWeightID=cms.string("2052"),
                                printLHE=cms.bool(False)
                                )


