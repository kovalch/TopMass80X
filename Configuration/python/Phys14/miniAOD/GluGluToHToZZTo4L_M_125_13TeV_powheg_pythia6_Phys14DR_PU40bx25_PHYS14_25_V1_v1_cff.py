import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU40bx25_PHYS14_25_V1-v1/00000/06462EA1-1872-E411-9440-7845C4FC3A13.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU40bx25_PHYS14_25_V1-v1/00000/0E9EF24D-E272-E411-AD62-7845C4FC37B5.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU40bx25_PHYS14_25_V1-v1/00000/26D910E5-B073-E411-9AD0-00266CFAE24C.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU40bx25_PHYS14_25_V1-v1/00000/4C224F8D-1672-E411-9947-001D09FDD91E.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU40bx25_PHYS14_25_V1-v1/10000/DA1E9747-7072-E411-A9DE-848F69FD457A.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU40bx25_PHYS14_25_V1-v1/10000/DE9F98E1-9E73-E411-8302-008CFA0020D4.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU40bx25_PHYS14_25_V1-v1/10000/F61F3822-A671-E411-AD4B-3417EBE2F4CC.root' ] );


secFiles.extend( [
               ] )

