import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/188E610D-9871-E411-BABD-002481E15008.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/76014B9F-B372-E411-B731-0025B3E0654E.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9C990A71-9771-E411-8F63-002590200A88.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BC22410E-9871-E411-A161-001E673971C5.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C4B6F654-9B71-E411-8219-001E673972C4.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/EC29E084-5D73-E411-B97B-002590200AE0.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F060F7F3-B171-E411-AE7A-002590200B38.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/208787DC-8D71-E411-BD05-001E6739672F.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/5CA9BF41-AA72-E411-A46B-002481E150EE.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/8A1F89D8-8D71-E411-BFDD-0025B3E06448.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/FCF68348-8E71-E411-BC8F-002590A37128.root' ] );


secFiles.extend( [
               ] )

