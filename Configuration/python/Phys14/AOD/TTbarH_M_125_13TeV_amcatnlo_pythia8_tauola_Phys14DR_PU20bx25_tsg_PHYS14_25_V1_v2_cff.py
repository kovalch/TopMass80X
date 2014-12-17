import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/1488A7F0-427F-E411-A025-001E67398408.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/18FF72AD-4E7F-E411-840E-001E67398156.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/20329E54-327F-E411-B47B-001E673972E7.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/38124FEF-347F-E411-ADC8-001E67398110.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/827B7EFC-3F7F-E411-97E2-002590200B3C.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/90086DA5-4E7F-E411-B3F3-001E67397698.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/AAF80C0C-397F-E411-BBA1-001E67397698.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/BE4AC439-437F-E411-8974-001E67397698.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/E6F49F81-2F7F-E411-8399-001E67397B11.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/F4D6C649-3C7F-E411-A8EF-001E67397698.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/FA27D654-2F7F-E411-8A3B-001E67397B11.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/10000/2ADE1A36-5F7F-E411-954F-0025905C38AA.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/10000/72E0F679-367F-E411-B42F-002590A8881E.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/10000/9E124FCD-457F-E411-8FC8-001E673969FA.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/10000/A869AFD0-587F-E411-89DC-002590A8881E.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/10000/B83B0538-287F-E411-BFAC-001E67397391.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/10000/FA6336EB-247F-E411-800C-001E673969FA.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/20000/A8CA477D-C77E-E411-ABD7-001E67396A22.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/30000/40E821B0-F57E-E411-858B-001E67397B11.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/30000/6E0C1CF6-EA7E-E411-9825-001E67397698.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/30000/C2362E1D-F77E-E411-A157-001E67397B11.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/40000/021E63B2-D87E-E411-8178-002590A370DC.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/40000/0E0C297C-BD7E-E411-98A4-002590A37106.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/40000/60F99499-D67E-E411-9727-002590A370DC.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/40000/7A22587F-D17E-E411-9B01-002590A4FFB8.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/40000/86DC2114-C67E-E411-A46A-002590A4FFB8.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/40000/9CDF211D-AE7E-E411-9D8B-001E6739722E.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/40000/AA80FDED-B67E-E411-B82C-002590A37106.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/40000/AC0AAB62-B37E-E411-892B-002590A370DC.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/40000/AE2100E8-C97E-E411-863A-002590A4FFB8.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/40000/EE456A27-C17E-E411-AB72-001E6739722E.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/40000/F69B157A-CE7E-E411-ABB5-002590A4FFB8.root' ] );


secFiles.extend( [
               ] )

