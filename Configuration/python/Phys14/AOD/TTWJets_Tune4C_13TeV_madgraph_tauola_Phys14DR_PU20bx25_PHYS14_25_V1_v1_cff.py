import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0476A930-2271-E411-A028-00266CF32EB8.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0651F1DB-1F71-E411-AE79-0025901D4B06.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/06D50EBB-4E71-E411-B245-002590AC4C8C.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0E609611-6071-E411-909F-0025901D4D20.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2ACBAB41-5771-E411-9E04-0025907DC9F8.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/36E3A819-D171-E411-95CD-002590AC4B54.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/427CD2A7-2471-E411-81F1-0025901D4936.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/444407F0-1E71-E411-96B9-002590A2CCE6.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/44CFC3AD-F670-E411-B673-00266CF2E2C8.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/60341755-5771-E411-95D6-002590DB9166.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/603A3850-2171-E411-AA1D-003048F0E82A.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/60D488EA-1E71-E411-B5BE-002590AC4C72.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/62287254-5771-E411-8359-002481E14F2A.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/628063D5-5A71-E411-9B6F-002590AC4C26.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7229B517-FF70-E411-938E-00266CFFA25C.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8A75F553-7D71-E411-AAAF-00266CF32F90.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/98D51B5A-7E71-E411-80BB-002481E0D3C0.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B239B7C2-0171-E411-8631-002590AC5482.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C6F23644-6071-E411-A8A7-00266CFFA7A8.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C6F5F1B5-4E71-E411-BC01-003048D437DA.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DA1A09FB-6871-E411-8D26-0025908325DC.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DAD9FF2B-1E71-E411-8603-0025901D4B06.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DC36B830-5771-E411-A248-003048D43982.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/EE5EBFE6-0A71-E411-835F-002590AC4CCC.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F0EC29CE-6971-E411-A610-00266CFFB390.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F273124E-7D71-E411-811B-002590DB05F4.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F2E37A37-5771-E411-AC7B-0025901D4090.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F8DBA31D-D171-E411-BF84-002590AC4CC8.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FA6E443A-5771-E411-B55B-00266CFFB7B8.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FCEDCA50-5771-E411-AE9C-003048F30422.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/3ED93C69-2271-E411-8037-0025907DC9BC.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/625057F3-7371-E411-9E24-00266CFFA754.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/6CAEE42C-BD71-E411-AC2A-0025901D47AA.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/84CF0732-0C71-E411-A92D-00266CF9B420.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/8E49B164-5B71-E411-A160-00266CF2E2C8.root',
       
'/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/9E57CF4D-BD71-E411-8751-0025907DCA0C.root' 
] );


secFiles.extend( [
               ] )

