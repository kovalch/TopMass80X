import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/0456F76A-EA77-E411-81F8-0025B31E3D3C.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/10C19613-E176-E411-92F7-F04DA23BBCCA.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/18F3D0C9-0177-E411-81D0-002590A37114.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/1CA79CF9-F876-E411-BD55-001E67396E05.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/28D3938A-F776-E411-B51A-00259020083C.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/28D5EAFC-F976-E411-A45C-002590A887EE.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/32AE7681-FC76-E411-A940-002481E14D64.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/36B52DAC-0077-E411-AA51-001E67396793.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/3CD5D12A-0877-E411-8C96-001E67397D64.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/4057BA6C-F876-E411-B260-0025902008C8.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/4059ABE0-E676-E411-81E5-002590A4C6BE.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/422588BA-F376-E411-82B1-001E67397AE4.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/46324C48-F076-E411-8143-001E67396892.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/507A8AF3-D276-E411-87E5-002590A370B2.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/68F5DFAD-FD76-E411-AAE0-002590200AE0.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/6C4F596C-E576-E411-B749-002590A4C69A.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/6E125622-FF76-E411-A4CC-001E67397747.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/70A88054-0F77-E411-A2DC-002590A8312A.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/70D0E1C9-DE76-E411-B619-00259020083C.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/76324C40-D776-E411-8F75-001E67397215.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/766E4FF0-FB76-E411-8E3D-002590200A58.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/7C93250E-F676-E411-8D70-002590A3711C.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/981D0C32-1778-E411-85B6-001E67397431.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/9A8376E4-E876-E411-B51D-001E67396905.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/9E1F21FE-ED76-E411-AC2A-001E67396973.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/A453C147-F776-E411-ABB1-001E67397A76.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/B020C664-E376-E411-920A-001E67396581.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/B85B358C-F576-E411-81B3-001E6739801B.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/BE2C75D8-EC76-E411-9507-002590A371AE.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/C0D697DA-EE76-E411-A7E0-002590200ABC.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/C22944BE-FA76-E411-81DF-002590A80E1E.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/C2E4E0B5-F176-E411-8525-D8D385FF4A72.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/C42ABF54-EB76-E411-9397-002590A4C69A.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/CEFB1297-F276-E411-9F2A-002590200900.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/D4BF93D5-F076-E411-87DA-002590A831CC.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/DC13DCE5-EB76-E411-BF43-001E67398223.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/DC9DBCA5-E976-E411-AD32-0025B3E05BBE.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/F6C064BF-3A77-E411-A171-001E6739677A.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/F8EC1CC0-0377-E411-ADE9-001E673987D2.root',
       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/FA84864E-F476-E411-B233-0025B3E063FA.root' ] );


secFiles.extend( [
               ] )

