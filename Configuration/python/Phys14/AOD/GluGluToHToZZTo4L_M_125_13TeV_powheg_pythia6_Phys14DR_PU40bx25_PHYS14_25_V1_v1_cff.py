import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/0A6F432E-E571-E411-9E94-F04DA275C013.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/20A54270-F371-E411-A134-3417EBE2F478.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/26435C1D-F771-E411-9600-848F69FD45B6.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/2ADB53E8-0172-E411-9724-3417EBE33927.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/2E94EE8A-1672-E411-8B75-3417EBE2EEC6.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/42C934EA-0472-E411-A19C-00266CFAE890.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/46A581AB-1172-E411-93BF-00266CFAE69C.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/480B9FE7-0972-E411-BC14-7845C4FC3770.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/54D641CA-C171-E411-B542-180373FF8456.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/62E256D3-1A72-E411-BFA3-00266CF897A0.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/689EC8A3-FD71-E411-931F-008CFA002830.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/76C1C90B-FA71-E411-83A0-008CFA002830.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/76D18EC3-EF71-E411-BE09-00A0D1EEE3B0.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/784AE1ED-DE71-E411-8DF1-3417EBE2EFCE.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/92379103-E071-E411-A11E-848F69FD47A5.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/B4AD2E28-EA71-E411-B33D-7845C4FC3770.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/B87254A7-E472-E411-A637-00266CF9B318.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/BAF857F5-CA71-E411-A3D7-3417EBE3379E.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/E49221F6-1372-E411-8C25-00A0D1EE8E60.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/E6BB8257-D171-E411-AB8B-7845C4F914D4.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/F65746B6-0772-E411-AF64-7845C4FC3C83.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/FCD5FC29-0E72-E411-AB72-00A0D1EE8E60.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/10000/0A32899F-0172-E411-BA7D-848F69FD29D0.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/10000/1CDAF6F2-D072-E411-8592-848F69FD477B.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/10000/4A807A5E-8271-E411-94F0-7845C4FC3665.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/10000/501CD42F-8371-E411-BDCC-7845C4FC35EA.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/10000/52475553-8471-E411-BAA2-24B6FDFF2B44.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/10000/62C2C3CB-9471-E411-B108-7845C4FC370A.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/10000/6C04578C-9571-E411-80D6-00A0D1EE8A14.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/10000/6C426A77-FE71-E411-ABDB-848F69FD29D0.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/10000/6CB80924-8371-E411-A418-848F69FD2925.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/10000/74071614-9471-E411-8F6D-180373FF8DE0.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/10000/BA00997F-A271-E411-8902-848F69FD294C.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/10000/BEC7D61F-8671-E411-AD12-008CFA001334.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/10000/D845065D-B371-E411-A2AA-00266CF9B59C.root',
       '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/AODSIM/PU40bx25_PHYS14_25_V1-v1/10000/FE8785D1-F771-E411-AE47-3417EBE2F0B2.root' ] );


secFiles.extend( [
               ] )

