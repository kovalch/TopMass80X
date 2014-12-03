import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0A10D0E9-D470-E411-BD51-0025B3E06612.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1420928C-F670-E411-B1F0-001E673968BA.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1A2F4E4D-AB70-E411-95AD-001E67398B29.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1CBCA10C-A070-E411-A581-002590A371D4.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2A3B43BC-E470-E411-81F5-002590A80DD4.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2EC77996-EE70-E411-A4FB-001E6739689C.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/32353E9A-D270-E411-BB0F-001E67398412.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/364FD09A-4971-E411-BABE-002590A887FE.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3A5200DD-6671-E411-B174-001E67396A1D.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4638726A-A672-E411-AD24-001E67398110.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/562B1EE1-C370-E411-918E-001E67396707.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/58575C8E-AF70-E411-8E04-002590200A58.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5E98FDBC-C971-E411-90D3-002481E7632A.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5ECADA7A-5573-E411-8F0A-002590A81DAC.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/640B8AD2-DC70-E411-A9B2-001E67396BAD.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/642E99FC-BC70-E411-8034-001E673972CE.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7818F383-C970-E411-8EAB-002590200A28.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/82D8AB64-C070-E411-9ADB-001E673982E6.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/84F9A825-C670-E411-AD8C-001E67398DE5.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B26BEC47-E370-E411-B430-001E67396DBA.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B88D6A78-CB70-E411-B79C-001E673968F1.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BCDBD2DC-DE70-E411-BDE6-001E673986B5.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C6258943-C870-E411-9989-001E67396DBA.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CE8380A9-A570-E411-8045-002590200970.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F2E6C4F2-FE70-E411-8F7C-001E6739732D.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FA706D32-DA70-E411-A344-001E67396D6F.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FEB8198E-EA70-E411-A041-001E67398E62.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/4EE4C93D-E671-E411-B620-001E67398223.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/5C5512C8-E070-E411-A34E-002590147CA2.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/60BF9CF0-D870-E411-87C9-002590A370DC.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/869F5A0F-A170-E411-8897-0025B3E05CF2.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/92BF8E1D-AB70-E411-9A8B-001E67398156.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/92DA4D8B-7072-E411-BD2F-002481E154CE.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/9E59A5CF-C870-E411-86B0-001E67397215.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/9E6FDA9F-C170-E411-833B-002590A80DFA.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/AA24DDA9-D070-E411-A8CE-001E67397616.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/BA18CE03-B970-E411-9613-001E67396FD1.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/BE8FB7C4-B170-E411-898A-F04DA23BBCCA.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/C83F3A3A-0971-E411-94CA-0025902008A8.root',
       '/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/CA7AAB13-9170-E411-A34D-002590200AE4.root' ] );


secFiles.extend( [
               ] )