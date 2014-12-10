import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/1AD11194-BA6E-E411-A8A9-7845C4FC39E3.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/340F5C3C-0F70-E411-847B-008CFA001F1C.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/4C3EEEC5-1570-E411-B513-008CFA002028.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/6C5145A2-CE6E-E411-8063-7845C4FC35F3.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/84C4A46F-0B70-E411-9AB3-848F69FD4DCC.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/90873375-1170-E411-BCE1-7845C4FC3AFD.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/A27E3450-D86E-E411-B134-7845C4FC3A40.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/A4906D4B-E06E-E411-86E3-00266CF9B9F0.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/B02C0AA1-EF6F-E411-8AD4-848F69FD29BB.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/B8F2ED64-E96E-E411-A110-00266CF25590.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/BA8637CD-D36E-E411-B737-3417EBE2F493.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/C8B35E41-1570-E411-8104-848F69FD46BB.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/CE5F2F08-C86E-E411-9AC9-3417EBE34BFD.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/D22A9A3D-8A6F-E411-89B0-7845C4F932D8.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/F8A0BA48-EA6F-E411-9538-00266CF250C4.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/201C0212-FF6F-E411-86C3-7845C4FC35E1.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/74CCB7F1-FB6F-E411-A833-7845C4FC361A.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/9E963B45-0370-E411-AAE6-00266CF9B254.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/BA80E262-1A70-E411-AAEE-848F69FD29BE.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/D459516D-2670-E411-8D2E-3417EBE2EEC6.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/DE96020F-FB6F-E411-8568-00266CF9C0F0.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/EEBA8311-F86F-E411-8047-7845C4FC361A.root' ] );


secFiles.extend( [
               ] )

