import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU40bx25_PHYS14_25_V1-v1/00000/1034F24B-1A78-E411-9D64-001E67398C1E.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU40bx25_PHYS14_25_V1-v1/00000/6856B40F-0C77-E411-893D-D8D385FF7678.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU40bx25_PHYS14_25_V1-v1/00000/A8D5A514-DD77-E411-95AB-002481E14D64.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU40bx25_PHYS14_25_V1-v1/00000/C20B68E7-0277-E411-85E5-001E67396A22.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU40bx25_PHYS14_25_V1-v1/00000/E601D324-FA76-E411-B2A3-0025B3E066A4.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU40bx25_PHYS14_25_V1-v1/00000/EC51B40A-0F77-E411-AB65-002590A831AA.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU40bx25_PHYS14_25_V1-v1/00000/ECF2804E-F476-E411-A097-001E67398110.root' ] );


secFiles.extend( [
               ] )

